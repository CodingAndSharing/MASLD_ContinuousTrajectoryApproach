"""Real, from-scratch regeneration of mergedcounts.csv from raw FASTQ.

mergedcounts_generation.py's --redownload only fetches small metadata/
supplementary files and prints ENA run accessions -- it stops well short of
downloading raw reads or running an alignment pipeline (see its docstring).
This script is the real thing: it downloads the actual raw paired-end FASTQ
for both cohorts and runs them through the paper's own pipeline --

    Raw FASTQ -> FastQC v0.11.9 -> HISAT2 v2.1.0 (GRCh38) -> HTSeq v0.11.1
    -> per-cohort gene-count matrices -> inner-joined -> low-count filter
    -> a fresh mergedcounts_from_raw_fastq.csv

-- entirely independently of the committed data/mergedcounts.csv.

Real, verified data volumes (queried live from ENA/GEO; see the two `fetch_
*_manifest()` functions below for how to reproduce these numbers yourself)
----------------------------------------------------------------------------
    Cohort            Samples   Raw FASTQ    Total reads
    UCAM (E-MTAB-9815)     58     21.6 GB          562M
    VCU/Sanyal (SRP197353) 78    316.6 GB        4.31B
    Combined              136   ~338 GB         ~4.87B

VCU/Sanyal is ~15x bigger than UCAM and was never sized in mergedcounts_
generation.py's docstring -- its raw reads live under SRA project SRP197353
(BioProject PRJNA542148), linked from GSE130970's family SOFT file, not in
the 2 Salmon-quantified supplementary files that script downloads.

Why this can only ever be a *best-effort* reproduction
----------------------------------------------------------------------------
Two pipeline parameters the paper's Methods never state exactly are guessed
here and clearly flagged where used:
  - GTF_RELEASE: which Ensembl GRCh38 annotation release HTSeq counted
    against (affects the exact gene set/IDs).
  - HTSEQ_STRANDEDNESS: library strandedness passed to `htseq-count -s`
    (wrong stranded-ness silently produces plausible-looking but wrong
    counts for a large fraction of genes).
Both are constants below -- override them if you know the real values.
Byte-identical reproduction of data/mergedcounts.csv is NOT expected; this
produces an independent matrix to validate/compare against it.

Resource cost -- read before running --full
----------------------------------------------------------------------------
~338 GB of downloads, a ~4.3 GB HISAT2 GRCh38 index, and (per informal
estimate, NOT benchmarked on any specific machine) roughly 4-12 days of
alignment on a machine with only ~8 GB RAM (HISAT2's GRCh38 index alone
needs ~4.3-4.9 GB resident, so realistically one alignment at a time).
--parallel only helps on a machine with enough RAM for N concurrent HISAT2
processes (~4-5 GB marginal each). Run --benchmark first to get a real
measured throughput number for your machine before committing to --full.

Where things are written
----------------------------------------------------------------------------
Everything downloaded/generated goes under DATA_DIR/temp/ (DATA_DIR follows
the same MASLD_DATA_DIR env var convention as mergedcounts_generation.py /
metadata_generation.py -- point it at sandbox/data_dev/ to keep multi-
hundred-GB downloads out of data/, which is not gitignored and IS tracked
by git). The final matrix is written to DATA_DIR/mergedcounts_from_raw_fastq.csv
-- a new file, so it never silently overwrites the committed mergedcounts.csv.

Usage
-----
    # Cheap, network-only sanity check (no downloads, no tools needed):
    python regenerate_from_raw_fastq.py --dry-run

    # Time ONE sample end-to-end and extrapolate a full-run estimate:
    MASLD_DATA_DIR=../../sandbox/data_dev python regenerate_from_raw_fastq.py --benchmark

    # The real, multi-day, ~338 GB run (do this on a bigger machine):
    MASLD_DATA_DIR=../../sandbox/data_dev python regenerate_from_raw_fastq.py --full --parallel 4
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
import os
import shutil
import subprocess
import sys
import tarfile
import time
import urllib.error
import urllib.request
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, asdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mergedcounts_generation import (  # noqa: E402
    DATA_DIR,
    EXPECTED_DATASET_COUNTS,
    METADATA_PATH,
    OUTLIER_SAMPLE_NAME,
    filter_low_expressed_genes,
    warn_and_confirm,
)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
TEMP_DIR = DATA_DIR / "temp"
RAW_FASTQ_DIR = TEMP_DIR / "raw_fastq"
REFERENCE_DIR = TEMP_DIR / "reference"
COUNTS_DIR = TEMP_DIR / "htseq_counts"
MANIFEST_DIR = TEMP_DIR / "manifests"
OUTPUT_PATH = DATA_DIR / "mergedcounts_from_raw_fastq.csv"

# ---------------------------------------------------------------------------
# Sample manifest sources (see module docstring for verified totals)
# ---------------------------------------------------------------------------
ARRAYEXPRESS_ACCESSION = "E-MTAB-9815"  # UCAM cohort, 58 samples
ARRAYEXPRESS_SDRF_URL = f"https://www.ebi.ac.uk/biostudies/files/{ARRAYEXPRESS_ACCESSION}/{ARRAYEXPRESS_ACCESSION}.sdrf.txt"

GEO_ACCESSION = "GSE130970"  # VCU/Sanyal cohort, 78 samples
GEO_SERIES_DIR = GEO_ACCESSION[:-3] + "nnn"
GEO_FAMILY_SOFT_URL = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{GEO_SERIES_DIR}/{GEO_ACCESSION}/soft/{GEO_ACCESSION}_family.soft.gz"
SRA_PROJECT_ACCESSION = "SRP197353"  # VCU/Sanyal raw reads (BioProject PRJNA542148)

ENA_FILEREPORT_URL = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport"
    "?accession={acc}&result=read_run&fields={fields}&format=tsv"
)

# ---------------------------------------------------------------------------
# Reference genome + annotation
# ---------------------------------------------------------------------------
# Official prebuilt HISAT2 GRCh38 genome index (no known-splice-sites/SNPs --
# matches a plain "HISAT2 v2.1.0 alignment to GRCh38", the level of detail
# the paper's Methods gives). ~4.3 GB compressed.
HISAT2_INDEX_URL = "https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz"
HISAT2_INDEX_PREFIX = REFERENCE_DIR / "grch38" / "genome"

# ASSUMPTION (paper does not state the exact annotation release used for
# HTSeq gene-level counting) -- override if you know the real one.
GTF_RELEASE = 110
GTF_URL = f"https://ftp.ensembl.org/pub/release-{GTF_RELEASE}/gtf/homo_sapiens/Homo_sapiens.GRCh38.{GTF_RELEASE}.gtf.gz"
GTF_PATH = REFERENCE_DIR / f"Homo_sapiens.GRCh38.{GTF_RELEASE}.gtf"

# ASSUMPTION (paper does not state library strandedness) -- "reverse" is the
# common case for Illumina TruSeq Stranded mRNA-seq. Override with
# --strandedness if you know the real protocol.
HTSEQ_STRANDEDNESS = "reverse"  # htseq-count -s : yes | no | reverse

REQUIRED_TOOLS = ["fastqc", "hisat2", "samtools", "htseq-count"]


@dataclass
class SampleRef:
    sample_name: str  # matches metadata.csv's "Sample name" / "Title"
    dataset: str  # "UCAM" or "SANYAL"
    run_accession: str  # ENA/SRA run accession (ERR.../SRR...)
    fastq_urls: list[str]  # 2 URLs (R1, R2), no scheme (ENA convention)
    fastq_bytes: int
    read_count: int


# ---------------------------------------------------------------------------
# Manifest fetching (network reads only -- no downloads of the actual data)
# ---------------------------------------------------------------------------
def _ena_filereport(accession: str, fields: str) -> list[dict[str, str]]:
    url = ENA_FILEREPORT_URL.format(acc=accession, fields=fields)
    with urllib.request.urlopen(url, timeout=60) as response:
        text = response.read().decode("utf-8", errors="replace")
    lines = text.strip().splitlines()
    if len(lines) < 2:
        return []
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:]]


def fetch_ucam_manifest() -> list[SampleRef]:
    """UCAM (E-MTAB-9815): 58 samples, one ENA run each, listed directly in
    the ArrayExpress SDRF. 'Source Name' matches metadata.csv's UCAM
    'Sample name' values (see metadata_generation.py's validate_ucam)."""
    with urllib.request.urlopen(ARRAYEXPRESS_SDRF_URL, timeout=30) as response:
        text = response.read().decode("utf-8", errors="replace")
    lines = text.splitlines()
    header = lines[0].split("\t")
    source_col = header.index("Source Name")
    run_col = header.index("Comment[ENA_RUN]")

    sample_by_run: dict[str, str] = {}
    for line in lines[1:]:
        cols = line.split("\t")
        if len(cols) > run_col and cols[run_col]:
            sample_by_run[cols[run_col]] = cols[source_col]

    # ENA's filereport endpoint silently returns 0 rows for a batch
    # comma-separated accession list (verified empirically) -- query one run
    # at a time instead; 58 requests is cheap.
    samples = []
    for run in sorted(sample_by_run):
        rows = _ena_filereport(run, "run_accession,fastq_ftp,fastq_bytes,read_count")
        if not rows:
            raise RuntimeError(f"ENA filereport returned nothing for run {run}")
        row = rows[0]
        fastq_urls = row["fastq_ftp"].split(";")
        fastq_bytes = sum(int(b) for b in row["fastq_bytes"].split(";") if b)
        samples.append(
            SampleRef(
                sample_name=sample_by_run[run],
                dataset="UCAM",
                run_accession=run,
                fastq_urls=fastq_urls,
                fastq_bytes=fastq_bytes,
                read_count=int(row["read_count"]),
            )
        )
    return samples


def fetch_vcu_manifest() -> list[SampleRef]:
    """VCU/Sanyal (GSE130970): 78 samples. GEO's family SOFT file gives each
    sample's title + linked SRA experiment (SRX); SRP197353's ENA filereport
    gives each run's (SRR) fastq_ftp/bytes/read_count. Joined on SRX."""
    with urllib.request.urlopen(GEO_FAMILY_SOFT_URL, timeout=30) as response:
        text = gzip.decompress(response.read()).decode("utf-8", errors="replace")

    title_by_srx: dict[str, str] = {}
    current_title = None
    for line in text.splitlines():
        if line.startswith("^SAMPLE"):
            current_title = None
        elif line.startswith("!Sample_title"):
            current_title = line.split("=", 1)[1].strip()
        elif line.startswith("!Sample_relation") and "term=SRX" in line:
            srx = line.rsplit("term=", 1)[1].strip()
            if current_title is not None:
                title_by_srx[srx] = current_title

    rows = _ena_filereport(
        SRA_PROJECT_ACCESSION, "run_accession,experiment_accession,fastq_ftp,fastq_bytes,read_count"
    )
    samples = []
    for row in rows:
        srx = row["experiment_accession"]
        if srx not in title_by_srx:
            continue
        fastq_urls = row["fastq_ftp"].split(";")
        fastq_bytes = sum(int(b) for b in row["fastq_bytes"].split(";") if b)
        samples.append(
            SampleRef(
                sample_name=title_by_srx[srx],
                dataset="SANYAL",
                run_accession=row["run_accession"],
                fastq_urls=fastq_urls,
                fastq_bytes=fastq_bytes,
                read_count=int(row["read_count"]),
            )
        )
    if len(samples) != EXPECTED_DATASET_COUNTS["SANYAL"]:
        print(
            f"[WARN] Expected {EXPECTED_DATASET_COUNTS['SANYAL']} VCU/Sanyal samples, "
            f"matched {len(samples)} SRX->title pairs."
        )
    return samples


def load_or_fetch_manifest(refresh: bool = False) -> list[SampleRef]:
    MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
    cache_path = MANIFEST_DIR / "samples.json"
    if cache_path.exists() and not refresh:
        raw = json.loads(cache_path.read_text())
        return [SampleRef(**s) for s in raw]

    print("[manifest] Fetching UCAM sample manifest from ArrayExpress/ENA ...")
    ucam = fetch_ucam_manifest()
    print(f"[manifest] UCAM: {len(ucam)} samples.")

    print("[manifest] Fetching VCU/Sanyal sample manifest from GEO/ENA ...")
    vcu = fetch_vcu_manifest()
    print(f"[manifest] VCU/Sanyal: {len(vcu)} samples.")

    samples = ucam + vcu
    cache_path.write_text(json.dumps([asdict(s) for s in samples], indent=2))
    return samples


# ---------------------------------------------------------------------------
# Reference setup
# ---------------------------------------------------------------------------
def ensure_hisat2_index() -> None:
    if list(REFERENCE_DIR.glob("grch38/genome.*.ht2")):
        return
    REFERENCE_DIR.mkdir(parents=True, exist_ok=True)
    tar_path = REFERENCE_DIR / "grch38_genome.tar.gz"
    print(f"[reference] Downloading HISAT2 GRCh38 index ({HISAT2_INDEX_URL}) -- ~4.3 GB ...")
    urllib.request.urlretrieve(HISAT2_INDEX_URL, tar_path)
    with tarfile.open(tar_path) as tar:
        tar.extractall(REFERENCE_DIR / "grch38")
    tar_path.unlink()


def ensure_gtf() -> None:
    if GTF_PATH.exists():
        return
    REFERENCE_DIR.mkdir(parents=True, exist_ok=True)
    gz_path = GTF_PATH.with_suffix(GTF_PATH.suffix + ".gz")
    print(f"[reference] Downloading Ensembl release {GTF_RELEASE} GTF ({GTF_URL}) ...")
    urllib.request.urlretrieve(GTF_URL, gz_path)
    with gzip.open(gz_path, "rb") as src, open(GTF_PATH, "wb") as dst:
        shutil.copyfileobj(src, dst)
    gz_path.unlink()


def check_required_tools() -> None:
    missing = [t for t in REQUIRED_TOOLS if shutil.which(t) is None]
    if missing:
        raise RuntimeError(
            f"Missing required tool(s) on PATH: {', '.join(missing)}. Install via, e.g.:\n"
            f"    pixi add -c bioconda -c conda-forge fastqc hisat2 samtools htseq"
        )


# ---------------------------------------------------------------------------
# Per-sample pipeline
# ---------------------------------------------------------------------------
@dataclass
class SampleTimings:
    sample_name: str
    download_sec: float
    fastqc_sec: float
    align_sec: float
    count_sec: float
    read_count: int
    fastq_bytes: int


def download_fastq(sample: SampleRef, dest_dir: Path) -> tuple[Path, Path]:
    dest_dir.mkdir(parents=True, exist_ok=True)
    if len(sample.fastq_urls) != 2:
        raise ValueError(f"{sample.sample_name}: expected 2 fastq URLs (paired-end), got {sample.fastq_urls}")
    r1_path = dest_dir / "R1.fastq.gz"
    r2_path = dest_dir / "R2.fastq.gz"
    for url, dest in [(sample.fastq_urls[0], r1_path), (sample.fastq_urls[1], r2_path)]:
        if dest.exists():
            continue
        full_url = url if url.startswith("http") else f"https://{url}"
        urllib.request.urlretrieve(full_url, dest)
    return r1_path, r2_path


def run_fastqc(r1: Path, r2: Path, dest_dir: Path) -> None:
    dest_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run(["fastqc", "-o", str(dest_dir), str(r1), str(r2)], check=True, capture_output=True)


def run_hisat2_align(r1: Path, r2: Path, dest_bam: Path, threads: int) -> None:
    dest_bam.parent.mkdir(parents=True, exist_ok=True)
    name_sorted = dest_bam.with_name(dest_bam.stem + ".namesorted.bam")
    hisat2 = subprocess.Popen(
        ["hisat2", "-x", str(HISAT2_INDEX_PREFIX), "-1", str(r1), "-2", str(r2), "-p", str(threads)],
        stdout=subprocess.PIPE,
    )
    sort = subprocess.run(
        ["samtools", "sort", "-n", "-@", str(threads), "-o", str(name_sorted), "-"],
        stdin=hisat2.stdout,
        check=True,
    )
    hisat2.stdout.close()
    hisat2.wait()
    if hisat2.returncode != 0:
        raise RuntimeError(f"hisat2 exited {hisat2.returncode}")
    name_sorted.rename(dest_bam)


def run_htseq_count(bam: Path, dest_tsv: Path, strandedness: str) -> None:
    dest_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(dest_tsv, "w") as out:
        subprocess.run(
            [
                "htseq-count", "-f", "bam", "-r", "name", "-s", strandedness,
                "-i", "gene_id", "-t", "exon", "-m", "union",
                str(bam), str(GTF_PATH),
            ],
            stdout=out, check=True,
        )


def process_sample(
    sample: SampleRef, threads: int, strandedness: str, keep_intermediates: bool, skip_fastqc: bool
) -> tuple[Path, SampleTimings]:
    sample_dir = RAW_FASTQ_DIR / sample.sample_name.replace("/", "_")
    bam_path = TEMP_DIR / "bam" / f"{sample.sample_name.replace('/', '_')}.bam"
    counts_path = COUNTS_DIR / f"{sample.sample_name.replace('/', '_')}.tsv"

    t0 = time.monotonic()
    r1, r2 = download_fastq(sample, sample_dir)
    t1 = time.monotonic()

    if not skip_fastqc:
        run_fastqc(r1, r2, TEMP_DIR / "fastqc" / sample.sample_name.replace("/", "_"))
    t2 = time.monotonic()

    run_hisat2_align(r1, r2, bam_path, threads)
    t3 = time.monotonic()

    run_htseq_count(bam_path, counts_path, strandedness)
    t4 = time.monotonic()

    if not keep_intermediates:
        shutil.rmtree(sample_dir, ignore_errors=True)
        bam_path.unlink(missing_ok=True)

    timings = SampleTimings(
        sample_name=sample.sample_name,
        download_sec=t1 - t0,
        fastqc_sec=t2 - t1,
        align_sec=t3 - t2,
        count_sec=t4 - t3,
        read_count=sample.read_count,
        fastq_bytes=sample.fastq_bytes,
    )
    return counts_path, timings


# ---------------------------------------------------------------------------
# Assembling the final matrix
# ---------------------------------------------------------------------------
def load_htseq_counts(counts_path: Path) -> pd.Series:
    """htseq-count output: 2-col TSV, gene_id \\t count, plus trailing
    __no_feature/__ambiguous/__too_low_aQual/__not_aligned/__alignment_not_unique
    summary rows (dropped -- these aren't genes)."""
    series = pd.read_csv(counts_path, sep="\t", header=None, index_col=0).iloc[:, 0]
    return series[~series.index.str.startswith("__")]


def assemble_cohort_matrix(sample_counts: dict[str, Path]) -> pd.DataFrame:
    columns = {name: load_htseq_counts(path) for name, path in sample_counts.items()}
    return pd.DataFrame(columns)


def merge_and_write(sample_counts: dict[str, Path], samples: list[SampleRef]) -> None:
    by_dataset: dict[str, dict[str, Path]] = {"UCAM": {}, "SANYAL": {}}
    for s in samples:
        if s.sample_name in sample_counts:
            by_dataset[s.dataset][s.sample_name] = sample_counts[s.sample_name]

    ucam_matrix = assemble_cohort_matrix(by_dataset["UCAM"])
    vcu_matrix = assemble_cohort_matrix(by_dataset["SANYAL"])
    print(f"[merge] UCAM matrix: {ucam_matrix.shape}, VCU/Sanyal matrix: {vcu_matrix.shape}")

    merged = ucam_matrix.join(vcu_matrix, how="inner")
    filtered = filter_low_expressed_genes(merged)

    metadata = pd.read_csv(METADATA_PATH, index_col=0)
    n_counts, n_meta = filtered.shape[1], len(sample_counts)
    print(f"[OK] Assembled {n_counts}/{len(samples)} requested samples, {filtered.shape[0]} genes.")
    if OUTLIER_SAMPLE_NAME not in filtered.columns:
        print(f"[note] Known outlier {OUTLIER_SAMPLE_NAME!r} not in this run's sample set.")

    filtered.to_csv(OUTPUT_PATH)
    print(f"[write] {OUTPUT_PATH} written ({filtered.shape[0]} genes x {filtered.shape[1]} samples).")
    print(
        f"[note] This used Ensembl release {GTF_RELEASE} annotation and "
        f"htseq-count -s {HTSEQ_STRANDEDNESS} (both assumptions -- see module docstring). "
        "Compare against data/mergedcounts.csv; differences may reflect these assumptions "
        "rather than a bug."
    )


# ---------------------------------------------------------------------------
# Modes
# ---------------------------------------------------------------------------
def run_dry_run(samples: list[SampleRef]) -> None:
    total_bytes = sum(s.fastq_bytes for s in samples)
    total_reads = sum(s.read_count for s in samples)
    by_dataset: dict[str, list[SampleRef]] = {}
    for s in samples:
        by_dataset.setdefault(s.dataset, []).append(s)
    print(f"\n[dry-run] {len(samples)} samples total, {total_bytes/1e9:.1f} GB, {total_reads:,} reads.")
    for dataset, group in by_dataset.items():
        b = sum(s.fastq_bytes for s in group)
        r = sum(s.read_count for s in group)
        print(f"  {dataset}: {len(group)} samples, {b/1e9:.1f} GB, {r:,} reads")
    print("No downloads performed, no tools required for this check.")


def run_benchmark(samples: list[SampleRef], threads: int, strandedness: str, skip_fastqc: bool) -> None:
    check_required_tools()
    ensure_hisat2_index()
    ensure_gtf()

    sample = min(samples, key=lambda s: s.fastq_bytes)
    print(f"\n[benchmark] Running the full pipeline on 1 sample: {sample.sample_name} "
          f"({sample.dataset}, {sample.fastq_bytes/1e9:.2f} GB, {sample.read_count:,} reads).")

    _, timings = process_sample(sample, threads=threads, strandedness=strandedness,
                                 keep_intermediates=False, skip_fastqc=skip_fastqc)

    total_sec = timings.download_sec + timings.fastqc_sec + timings.align_sec + timings.count_sec
    print(
        f"\n[benchmark] Sample {timings.sample_name}: download {timings.download_sec:.0f}s, "
        f"fastqc {timings.fastqc_sec:.0f}s, align {timings.align_sec:.0f}s, "
        f"count {timings.count_sec:.0f}s (total {total_sec:.0f}s)."
    )

    total_reads_all = sum(s.read_count for s in samples)
    total_bytes_all = sum(s.fastq_bytes for s in samples)
    align_reads_per_sec = timings.read_count / max(timings.align_sec, 1e-6)
    count_reads_per_sec = timings.read_count / max(timings.count_sec, 1e-6)
    download_bytes_per_sec = timings.fastq_bytes / max(timings.download_sec, 1e-6)

    est_download_hr = (total_bytes_all / download_bytes_per_sec) / 3600
    est_align_hr = (total_reads_all / align_reads_per_sec) / 3600
    est_count_hr = (total_reads_all / count_reads_per_sec) / 3600
    est_total_hr = est_download_hr + est_align_hr + est_count_hr

    print(
        f"\n[benchmark] Extrapolated to all {len(samples)} samples "
        f"({total_bytes_all/1e9:.0f} GB, {total_reads_all:,} reads), assuming this machine's "
        f"1-sample throughput holds and downloads/alignment run serially (no --parallel):\n"
        f"  download ~{est_download_hr:.1f} h, align ~{est_align_hr:.1f} h, "
        f"count ~{est_count_hr:.1f} h -> total ~{est_total_hr:.1f} h (~{est_total_hr/24:.1f} days).\n"
        f"  With --parallel N (N concurrent samples, RAM permitting), align/count scale down "
        f"roughly by N; download is bandwidth-bound and scales less cleanly."
    )


def run_full(
    samples: list[SampleRef], threads: int, parallel: int, strandedness: str,
    keep_intermediates: bool, skip_fastqc: bool, limit: int | None,
) -> None:
    check_required_tools()
    ensure_hisat2_index()
    ensure_gtf()

    if limit is not None:
        samples = samples[:limit]

    print(f"\n[full] Processing {len(samples)} samples with --parallel {parallel} ...")
    sample_counts: dict[str, Path] = {}
    completed = 0

    if parallel <= 1:
        for sample in samples:
            counts_path, timings = process_sample(sample, threads, strandedness, keep_intermediates, skip_fastqc)
            sample_counts[sample.sample_name] = counts_path
            completed += 1
            print(f"[full] ({completed}/{len(samples)}) {sample.sample_name} done "
                  f"(align {timings.align_sec:.0f}s, count {timings.count_sec:.0f}s).")
    else:
        with ProcessPoolExecutor(max_workers=parallel) as pool:
            futures = {
                pool.submit(process_sample, s, threads, strandedness, keep_intermediates, skip_fastqc): s
                for s in samples
            }
            for future in as_completed(futures):
                sample = futures[future]
                counts_path, timings = future.result()
                sample_counts[sample.sample_name] = counts_path
                completed += 1
                print(f"[full] ({completed}/{len(samples)}) {sample.sample_name} done "
                      f"(align {timings.align_sec:.0f}s, count {timings.count_sec:.0f}s).")

    merge_and_write(sample_counts, samples)


def main() -> None:
    args = _parse_args()

    if not args.dry_run and DATA_DIR == (Path(__file__).resolve().parents[2] / "data"):
        if not warn_and_confirm(
            f"Writing ~338 GB of raw FASTQ + reference data into {TEMP_DIR} -- this is inside "
            "data/, which IS tracked by git (unlike sandbox/data_dev/). Set MASLD_DATA_DIR to "
            "point at sandbox/data_dev/ instead unless you really mean to use data/.",
            "~338 GB + ~4.3 GB reference index",
            args.yes,
        ):
            print("[abort] Not proceeding without confirmation.")
            return

    samples = load_or_fetch_manifest(refresh=args.refresh_manifest)

    if args.dry_run:
        run_dry_run(samples)
    elif args.benchmark:
        run_benchmark(samples, threads=args.threads, strandedness=args.strandedness, skip_fastqc=args.skip_fastqc)
    elif args.full:
        run_full(
            samples, threads=args.threads, parallel=args.parallel, strandedness=args.strandedness,
            keep_intermediates=args.keep_intermediates, skip_fastqc=args.skip_fastqc, limit=args.limit,
        )
    else:
        print("Nothing to do -- pass --dry-run, --benchmark, or --full. See --help.")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--dry-run", action="store_true", help="Only fetch+print the sample manifest (no downloads).")
    mode.add_argument("--benchmark", action="store_true", help="Run the full pipeline on 1 sample and extrapolate.")
    mode.add_argument("--full", action="store_true", help="Run the full pipeline on all (or --limit) samples.")
    parser.add_argument("--limit", type=int, default=None, help="With --full, only process the first N samples.")
    parser.add_argument("--parallel", type=int, default=1, help="Concurrent samples for --full (RAM permitting).")
    parser.add_argument("--threads", type=int, default=8, help="HISAT2/samtools threads per sample.")
    parser.add_argument("--strandedness", choices=["yes", "no", "reverse"], default=HTSEQ_STRANDEDNESS,
                         help="htseq-count -s value (see module docstring: this is an assumption).")
    parser.add_argument("--keep-intermediates", action="store_true", help="Keep per-sample FASTQ/BAM after counting.")
    parser.add_argument("--skip-fastqc", action="store_true", help="Skip the FastQC QC step.")
    parser.add_argument("--refresh-manifest", action="store_true", help="Re-fetch the sample manifest instead of using the cached copy.")
    parser.add_argument("--yes", action="store_true", help="Skip the data/ vs sandbox/data_dev/ confirmation prompt.")
    return parser.parse_args()


if __name__ == "__main__":
    main()
