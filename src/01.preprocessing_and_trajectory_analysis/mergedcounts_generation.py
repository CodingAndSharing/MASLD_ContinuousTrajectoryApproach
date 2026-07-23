"""Documents and (re)generates data/mergedcounts.csv.

data/mergedcounts.csv is the raw, un-normalized UCAM + VCU/Sanyal gene-count
matrix read at the top of every script in this directory (e.g. "Batch
correction UCAM-VCU (for dataset + gender) after_removing_outlier.R",
"Normalization_and_PCA.R"). Until now no script in this repo produced it --
see sandbox/publication/NOTES_CLAUDE.md, item 4, and the "Provenance of the
UCAM/VCU upstream pipeline files" table in
sandbox/publication/methods_for_data.ipynb.

True provenance (paper's Methods, p.11 "RNA-seq analysis and data
processing", and Reporting Summary, p.30 -- see
sandbox/publication/MASDL_framework.pdf):

    Raw FASTQ
      - ArrayExpress E-MTAB-9815  (UCAM cohort,        58 patients)
      - GEO GSE130970             (VCU/Sanyal cohort,   4 controls + 74 patients)
    -> FastQC v0.11.9   quality control
    -> HISAT2 v2.1.0    alignment to GRCh38
    -> HTSeq v0.11.1    gene-level (Ensembl ID) counts, per sample, per cohort
    -> the two cohorts' count matrices concatenated into one 136-sample matrix
    -> low-count genes removed (mean count >= 1 across samples)  [recovered -- see below]
    -> data/mergedcounts.csv

The paper never names this file and never narrates the concatenation step
explicitly -- it sits implicitly between "gene-level counts were obtained
with HTSeq" (stated per dataset) and "Quantile normalization was applied
across datasets" (which presupposes one combined matrix). This script makes
that merge step explicit.

Low-count expression filter (recovered from the committed data)
---------------------------------------------------------------
The committed data/mergedcounts.csv is NOT raw HTSeq output: every one of its
17,090 genes has a mean count >= 1 across the 136 samples (minimum observed =
137/136 = 1.0074; zero all-zero rows; zero genes detected in <=1 sample). That
hard floor at mean = 1 is the fingerprint of a low-count filter -- keep a gene
iff rowMeans(counts) >= 1 (equivalently rowSums >= n_samples) -- applied during
the original, unstored HTSeq->merge step. No committed script implemented it,
even though Normalization_and_PCA.R's comment ("exclude low expressed counts")
refers to exactly this step. filter_low_expressed_genes() below makes it
explicit and reproducible: applied to the committed matrix it is a no-op (drops
0 of 17,090 genes), which both validates the recovered threshold and keeps the
output byte-identical; applied to a freshly merged raw matrix it does the real
filtering ahead of the rank-based quantile normalization that follows.

NOTE: this faithfully reproduces the ORIGINAL absolute-count threshold. Given
the 62-fold library-size range here (0.62M .. 38.5M reads), a depth-aware
filter (edgeR::filterByExpr / CPM) is preferable and would drop ~554 more
genes -- that upgrade belongs in the src_dev rewrite, not in this script.

What this script can and cannot reproduce from scratch
--------------------------------------------------------
Raw FASTQ files and the two cohorts' *separate* per-sample HTSeq count
matrices are not stored in this repository -- only the already-concatenated
136-sample matrix is (data/mergedcounts.csv, byte-identical to the copy at
data/ucam_sanyal/counts_matrix.csv, confirmed 2026-07-21). Producing the
separate matrices from scratch would require re-downloading FASTQ from
ArrayExpress/GEO and re-running FastQC/HISAT2/HTSeq, which is out of scope
for a lightweight script.

If the two per-cohort matrices ever become available locally, point
UCAM_RAW_COUNTS / VCU_RAW_COUNTS below at them and merge_from_raw_htseq_counts()
will perform the real concatenation. Until then, main() falls back to
validating the existing merged file against data/metadata.csv (sample count,
dataset composition, the excluded outlier) and can optionally rewrite it with
sample-name column headers instead of the generic placeholders it ships with.

Getting the per-cohort matrices: --redownload
------------------------------------------------
`--redownload` fetches the two cohorts' publicly deposited data into
data/raw/ so a later run can pick them up via UCAM_RAW_COUNTS/VCU_RAW_COUNTS:

  - GEO GSE130970 (VCU/Sanyal): 2 supplementary files, downloaded directly.
    Real sizes, verified via HTTP Content-Length on 2026-07-21:
      GSE130970_all_sample_salmon_tximport_TPM_entrez_gene_ID.csv.gz     6,027,512 B  (~5.75 MB)
      GSE130970_all_sample_salmon_tximport_counts_entrez_gene_ID.csv.gz 11,482,054 B (~10.95 MB)
      TOTAL                                                            17,509,566 B  (~16.7 MB)
    Note: these are the VCU/Sanyal study's own Salmon-based quantification,
    not this paper's HISAT2/HTSeq reprocessing -- see merge_from_raw_htseq_counts().

  - ArrayExpress E-MTAB-9815 (UCAM): ArrayExpress/BioStudies hosts only
    metadata for this record (IDF/SDRF files, a few hundred KB); the raw
    FASTQ themselves are deposited at ENA under 58 run accessions listed in
    the SDRF's ENA_RUN/FASTQ_URI columns (one run per UCAM sample). Real
    total size, verified via the ENA Portal filereport API on 2026-07-21,
    summed across all 58 runs (paired-end R1+R2):
      21,616,610,827 B (~21.6 GB) total raw FASTQ

Because these downloads can be large and this repo has no way to know the
user's bandwidth/disk constraints in advance, --redownload always prints a
warning with the (approximate) size and requires explicit confirmation --
either an interactive "yes", or --yes to skip the prompt for non-interactive
use -- before touching the network. Nothing is downloaded automatically
without that confirmation, and --redownload does nothing unless requested.

Usage
-----
    python mergedcounts_generation.py                    # validate only (no files changed)
    python mergedcounts_generation.py --write             # also rewrite mergedcounts.csv
                                                            # with sample-name columns
    python mergedcounts_generation.py --redownload         # fetch raw per-cohort data
                                                            # into data/raw/ (asks first)
    python mergedcounts_generation.py --redownload --yes   # same, no confirmation prompt
"""

import argparse
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[2]  # .../MASLD_ContinuousTrajectoryApproach
DATA_DIR = REPO_ROOT / "data"

METADATA_PATH = DATA_DIR / "metadata.csv"
MERGED_COUNTS_PATH = DATA_DIR / "mergedcounts.csv"  # the table this script produces
DUPLICATE_COPY_PATH = DATA_DIR / "ucam_sanyal" / "counts_matrix.csv"  # known duplicate, for cross-check

# Hypothetical per-cohort raw HTSeq outputs -- NOT present in this repo today.
# See module docstring: if you obtain these (e.g. by re-running the
# FastQC/HISAT2/HTSeq pipeline on the ArrayExpress/GEO FASTQ files), point
# these constants at them and merge_from_raw_htseq_counts() will do the real
# concatenation instead of falling back to the pre-merged file.
UCAM_RAW_COUNTS = DATA_DIR / "raw" / "ucam_htseq_counts.tsv"  # E-MTAB-9815, 58 samples
VCU_RAW_COUNTS = DATA_DIR / "raw" / "vcu_htseq_counts.tsv"  # GSE130970, 78 samples

EXPECTED_DATASET_COUNTS = {"UCAM": 58, "SANYAL": 78}  # Methods, p.11
OUTLIER_SAMPLE_NAME = "Sample 5"  # excluded downstream as an outlier (Supp. Fig. 1a,b)

# Low-count expression filter recovered from the committed matrix (see docstring):
# keep only genes whose MEAN count across all samples is >= this value. The
# committed mergedcounts.csv already satisfies this exactly (min mean = 137/136
# = 1.0074), so re-applying it is a no-op; on a freshly merged raw matrix it
# performs the real filtering.
MIN_MEAN_COUNT = 1.0

# ---------------------------------------------------------------------------
# Re-download sources
# ---------------------------------------------------------------------------
RAW_DIR = DATA_DIR / "raw"

GEO_ACCESSION = "GSE130970"  # VCU/Sanyal cohort
# GEO has no true single-archive bulk download for this series (it only
# bundles a *_RAW.tar when one was submitted; GSE130970 wasn't). Instead it
# lists 2 individual supplementary files -- fetch each directly by name.
GEO_SUPPLEMENTARY_FILES = [
    "GSE130970_all_sample_salmon_tximport_TPM_entrez_gene_ID.csv.gz",  # 6,027,512 B (~5.75 MB)
    "GSE130970_all_sample_salmon_tximport_counts_entrez_gene_ID.csv.gz",  # 11,482,054 B (~10.95 MB)
]
GEO_FILE_URL_TEMPLATE = "https://www.ncbi.nlm.nih.gov/geo/download/?acc={acc}&format=file&file={file}"
# Real, verified sizes (HTTP Content-Length, checked 2026-07-21): 6,027,512 +
# 11,482,054 = 17,509,566 bytes total, ~16.7 MB. Small -- these are the VCU/
# Sanyal study's own Salmon-based quantification (see file names), not this
# paper's HISAT2/HTSeq reprocessing.
GEO_ESTIMATED_SIZE = "~16.7 MB total (17,509,566 B; 2 files, verified 2026-07-21)"

ARRAYEXPRESS_ACCESSION = "E-MTAB-9815"  # UCAM cohort
ARRAYEXPRESS_SDRF_URL = f"https://www.ebi.ac.uk/biostudies/files/{ARRAYEXPRESS_ACCESSION}/{ARRAYEXPRESS_ACCESSION}.sdrf.txt"
ENA_FILEREPORT_URL_TEMPLATE = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport"
    "?accession={acc}&result=read_run&fields=run_accession,fastq_ftp,fastq_bytes&format=tsv"
)
# ArrayExpress/BioStudies hosts only metadata for this record (IDF/SDRF, a
# few hundred KB) -- the raw FASTQ are on ENA under 58 run accessions listed
# in the SDRF's ENA_RUN column (one run per UCAM sample), fetched per-run
# below via print_arrayexpress_instructions()/fetch_ena_run_urls().
# Real, verified total (ENA Portal filereport API, summed across all 58
# paired-end runs, checked 2026-07-21): 21,616,610,827 bytes, ~21.6 GB.
ARRAYEXPRESS_ESTIMATED_SIZE = "~21.6 GB total (21,616,610,827 B; 58 paired-end runs, verified 2026-07-21)"


def merge_from_raw_htseq_counts(ucam_path: Path, vcu_path: Path) -> pd.DataFrame:
    """Concatenate the two cohorts' HTSeq gene-count matrices into one table.

    Both inputs are expected to be HTSeq-style tables (rows = Ensembl gene
    ID, columns = per-sample raw counts), exactly as produced by
    hisat2 v2.1.0 (GRCh38) | htseq-count v0.11.1 (Methods p.11; Reporting
    Summary p.30).
    """
    ucam = pd.read_csv(ucam_path, sep=None, engine="python", index_col=0)
    vcu = pd.read_csv(vcu_path, sep=None, engine="python", index_col=0)

    # Genes are matched by Ensembl ID (the row index); an inner join keeps
    # only genes quantified in both cohorts, matching how a single combined
    # matrix would be built ahead of cross-dataset quantile normalization and
    # ComBat batch correction.
    return ucam.join(vcu, how="inner")


def filter_low_expressed_genes(counts: pd.DataFrame, min_mean: float = MIN_MEAN_COUNT) -> pd.DataFrame:
    """Remove low-expressed genes: keep a gene iff its mean count across all
    samples is >= `min_mean` (default 1.0), i.e. rowSums >= min_mean * n_samples.

    This reproduces the low-count filter baked into the committed
    data/mergedcounts.csv (see module docstring). The "exclude low expressed
    counts" comment in Normalization_and_PCA.R refers to this step, but no
    committed script performed it until now.

    Rationale: the next step (quantile normalization, in Normalization_and_PCA.R)
    is rank-based and therefore sensitive to a large mass of tied near-zero
    counts, which distort the mean-of-order-statistics reference distribution
    every sample is mapped onto. Dropping never-/barely-expressed genes first is
    what makes that normalization well-behaved.

    Boundary note: no gene in the committed matrix sits exactly at mean == 1
    (the minimum is 137/136), so `>=` vs `>` is indistinguishable from the data;
    the conventional `>=` form is used here.
    """
    means = counts.mean(axis=1)
    keep = means >= min_mean
    n_dropped = int((~keep).sum())
    print(
        f"[filter] Low-count filter (mean count >= {min_mean}): "
        f"kept {int(keep.sum())} / {len(keep)} genes (dropped {n_dropped})."
    )
    return counts.loc[keep]


def validate_against_metadata(counts: pd.DataFrame, metadata: pd.DataFrame) -> None:
    """Sanity-check the merged matrix against data/metadata.csv (raw/curated
    per-sample clinical annotation -- Supplementary Tables S1/S2)."""
    n_counts, n_meta = counts.shape[1], len(metadata)
    if n_counts != n_meta:
        raise ValueError(
            f"Sample count mismatch: {MERGED_COUNTS_PATH.name} has {n_counts} "
            f"columns, metadata.csv has {n_meta} rows."
        )

    dataset_counts = metadata["Dataset"].value_counts().to_dict()
    if dataset_counts != EXPECTED_DATASET_COUNTS:
        raise ValueError(
            f"Dataset composition mismatch: got {dataset_counts}, expected "
            f"{EXPECTED_DATASET_COUNTS} (Methods, p.11)."
        )

    if OUTLIER_SAMPLE_NAME not in set(metadata["Sample name"]):
        raise ValueError(
            f"Expected outlier {OUTLIER_SAMPLE_NAME!r} not found in "
            "metadata.csv (excluded downstream as an outlier; Supp. Fig. 1a,b)."
        )

    print(
        f"[OK] {n_counts} samples, {counts.shape[0]} Ensembl genes; "
        f"dataset composition {dataset_counts}; outlier "
        f"{OUTLIER_SAMPLE_NAME!r} present (excluded later, in "
        "'Batch correction UCAM-VCU (...).R', not here)."
    )


def label_columns_from_metadata(counts: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    """Reproduce the column-labelling every downstream script performs inline
    (e.g. `colnames(merged_counts) = template$Sample.name` in "Batch
    correction UCAM-VCU (for dataset + gender) after_removing_outlier.R",
    line 21), so the generated table is self-describing rather than relying
    on generic "Sample 1", "Sample 2", ... headers matched purely by column
    order."""
    labelled = counts.copy()
    labelled.columns = metadata["Sample name"].values
    return labelled


def warn_and_confirm(description: str, estimated_size: str, assume_yes: bool) -> bool:
    """Print a size/network warning and ask for explicit confirmation before
    any download starts. Returns True only if the user (or --yes) agreed."""
    print(
        "\n"
        "==================== DOWNLOAD WARNING ====================\n"
        f"About to fetch: {description}\n"
        f"Estimated size: {estimated_size}\n"
        "This requires network access and may take a long time and/or use\n"
        "significant disk space and bandwidth, depending on your connection.\n"
        "You can skip this safely -- nothing downstream requires it unless\n"
        "you want to regenerate mergedcounts.csv from the original raw data\n"
        "instead of using the copy already in data/mergedcounts.csv.\n"
        "============================================================\n"
    )
    if assume_yes:
        print("[--yes] Skipping confirmation prompt.")
        return True
    try:
        reply = input("Proceed with this download? [y/N]: ").strip().lower()
    except EOFError:
        # No interactive terminal available (e.g. piped/non-interactive run)
        # and --yes was not passed -- refuse rather than silently downloading.
        print("[abort] No interactive input available and --yes not passed; not downloading.")
        return False
    return reply in ("y", "yes")


def download_geo_supplementary(assume_yes: bool) -> None:
    """Download GEO GSE130970's 2 supplementary files (VCU/Sanyal cohort's
    own Salmon-based quantification) directly by name -- see
    GEO_SUPPLEMENTARY_FILES. This is the one download in this script that
    can actually be automated end-to-end -- GEO serves each file over plain
    HTTP with no authentication."""
    pending = [f for f in GEO_SUPPLEMENTARY_FILES if not (RAW_DIR / f).exists()]
    if not pending:
        print(f"[skip] All {len(GEO_SUPPLEMENTARY_FILES)} GEO supplementary files already present in {RAW_DIR.relative_to(REPO_ROOT)}/.")
        return

    if not warn_and_confirm(
        f"GEO {GEO_ACCESSION} supplementary files (VCU/Sanyal cohort): {', '.join(pending)}",
        GEO_ESTIMATED_SIZE,
        assume_yes,
    ):
        print("[skip] User declined; not downloading GEO supplementary files.")
        return

    RAW_DIR.mkdir(parents=True, exist_ok=True)
    for filename in pending:
        url = GEO_FILE_URL_TEMPLATE.format(acc=GEO_ACCESSION, file=filename)
        dest = RAW_DIR / filename
        print(f"[download] Fetching {url} -> {dest.relative_to(REPO_ROOT)} ...")
        try:
            urllib.request.urlretrieve(url, dest)
        except (urllib.error.URLError, OSError) as exc:
            print(f"[error] Download failed ({exc}). You can retry, or fetch it manually from {url}.")
            continue
    print(
        f"[OK] GEO supplementary files saved to {RAW_DIR.relative_to(REPO_ROOT)}/. These are "
        "processed matrices (Salmon-based), not HTSeq counts -- point VCU_RAW_COUNTS at an "
        "HTSeq-format equivalent if you re-run the paper's own alignment pipeline on the raw FASTQ."
    )


def fetch_ena_run_accessions() -> list[str]:
    """Parse ArrayExpress E-MTAB-9815's SDRF metadata file for the 58 ENA run
    accessions (column ENA_RUN) that hold the UCAM cohort's raw FASTQ --
    ArrayExpress/BioStudies itself only stores the SDRF/IDF metadata, not the
    reads themselves."""
    with urllib.request.urlopen(ARRAYEXPRESS_SDRF_URL, timeout=30) as response:
        text = response.read().decode("utf-8", errors="replace")
    header = text.splitlines()[0].split("\t")
    run_col = header.index("Comment[ENA_RUN]")
    accessions = []
    for line in text.splitlines()[1:]:
        cols = line.split("\t")
        if len(cols) > run_col and cols[run_col]:
            accessions.append(cols[run_col])
    return sorted(set(accessions))


def print_arrayexpress_instructions(assume_yes: bool) -> None:
    """ArrayExpress/BioStudies (E-MTAB-9815, UCAM cohort) hosts only
    metadata; the raw FASTQ live on ENA under the 58 run accessions listed in
    the SDRF. This looks those up and prints per-run ENA download URLs plus
    the real total size, rather than fetching ~21.6 GB of FASTQ
    automatically."""
    if not warn_and_confirm(
        f"ArrayExpress {ARRAYEXPRESS_ACCESSION} (UCAM cohort) raw FASTQ, hosted on ENA "
        "-- this step only looks up run accessions/sizes, it does not download the FASTQ",
        ARRAYEXPRESS_ESTIMATED_SIZE,
        assume_yes,
    ):
        print("[skip] User declined; not querying ArrayExpress/ENA.")
        return

    print(f"[info] Fetching SDRF metadata from {ARRAYEXPRESS_SDRF_URL} ...")
    try:
        run_accessions = fetch_ena_run_accessions()
    except (urllib.error.URLError, OSError) as exc:
        print(f"[warn] Could not reach BioStudies ({exc}).")
        return

    print(
        f"[OK] Found {len(run_accessions)} ENA run accessions for {ARRAYEXPRESS_ACCESSION} "
        f"(one per UCAM sample), e.g. {run_accessions[0]} .. {run_accessions[-1]}."
    )
    print(
        f"\n[manual step required] Raw FASTQ total ~21.6 GB across {len(run_accessions)} runs -- "
        "not downloaded automatically by this script. To fetch a given run's FASTQ, query:\n"
        f"    {ENA_FILEREPORT_URL_TEMPLATE.format(acc=run_accessions[0], file='')}\n"
        "which returns the fastq_ftp URL(s) and fastq_bytes size for that run. Save the "
        f"resulting per-sample HTSeq counts (after FastQC/HISAT2/HTSeq) under "
        f"{RAW_DIR.relative_to(REPO_ROOT)}/ and point UCAM_RAW_COUNTS at it."
    )


def redownload_raw_data(assume_yes: bool) -> None:
    """Entry point for --redownload: fetch (or point at) the two cohorts'
    publicly deposited data. Each source is gated by its own warning +
    confirmation in warn_and_confirm(); declining one does not block the
    other."""
    download_geo_supplementary(assume_yes)
    print_arrayexpress_instructions(assume_yes)


def main(write: bool = False, redownload: bool = False, assume_yes: bool = False) -> None:
    if redownload:
        redownload_raw_data(assume_yes)

    metadata = pd.read_csv(METADATA_PATH, index_col=0)

    if UCAM_RAW_COUNTS.exists() and VCU_RAW_COUNTS.exists():
        print(
            "[merge] Per-cohort raw HTSeq matrices found -- performing the "
            "real UCAM + VCU/Sanyal concatenation described in Methods p.11."
        )
        counts = merge_from_raw_htseq_counts(UCAM_RAW_COUNTS, VCU_RAW_COUNTS)
        # Real regeneration path: apply the low-count filter that the committed
        # matrix already reflects (see filter_low_expressed_genes / docstring),
        # so a freshly merged matrix matches the published one's gene set.
        counts = filter_low_expressed_genes(counts, MIN_MEAN_COUNT)
    else:
        print(
            "[fallback] Per-cohort raw HTSeq matrices are not present in "
            f"this repo ({UCAM_RAW_COUNTS.relative_to(REPO_ROOT)}, "
            f"{VCU_RAW_COUNTS.relative_to(REPO_ROOT)} not found). Falling "
            f"back to validating the existing, already-merged "
            f"{MERGED_COUNTS_PATH.name} instead of regenerating it from "
            "scratch."
        )
        counts = pd.read_csv(MERGED_COUNTS_PATH, index_col=0)

        # The committed matrix is already filtered; re-applying the recovered
        # threshold must therefore be a no-op. Assert it, so this script also
        # serves as a check that MIN_MEAN_COUNT still matches the shipped data.
        refiltered = filter_low_expressed_genes(counts, MIN_MEAN_COUNT)
        if refiltered.shape[0] != counts.shape[0]:
            raise ValueError(
                f"Low-count filter (mean >= {MIN_MEAN_COUNT}) dropped "
                f"{counts.shape[0] - refiltered.shape[0]} genes from the committed "
                f"{MERGED_COUNTS_PATH.name}, which should already be filtered. "
                "The recovered threshold no longer matches the shipped data."
            )
        print(f"[OK] Committed matrix already satisfies mean count >= {MIN_MEAN_COUNT} (filter is a no-op).")

        if DUPLICATE_COPY_PATH.exists():
            duplicate = pd.read_csv(DUPLICATE_COPY_PATH, index_col=0)
            if counts.equals(duplicate):
                print(f"[OK] Matches the known duplicate copy at {DUPLICATE_COPY_PATH.relative_to(REPO_ROOT)}.")
            else:
                print(
                    "[WARN] Differs from the duplicate copy at "
                    f"{DUPLICATE_COPY_PATH.relative_to(REPO_ROOT)} -- these "
                    "were byte-identical when last checked (2026-07-21)."
                )

    validate_against_metadata(counts, metadata)
    labelled = label_columns_from_metadata(counts, metadata)

    if write:
        labelled.to_csv(MERGED_COUNTS_PATH)
        print(f"[write] {MERGED_COUNTS_PATH.relative_to(REPO_ROOT)} regenerated with sample-name columns.")
    else:
        print(
            "[dry-run] Not writing to disk (pass --write to regenerate "
            f"{MERGED_COUNTS_PATH.relative_to(REPO_ROOT)} with labelled columns)."
        )


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--write", action="store_true", help="Regenerate mergedcounts.csv with sample-name columns."
    )
    parser.add_argument(
        "--redownload",
        action="store_true",
        help="Fetch the two cohorts' raw/supplementary data into data/raw/ (asks for confirmation first).",
    )
    parser.add_argument(
        "--yes", action="store_true", help="Skip the confirmation prompt for --redownload (non-interactive use)."
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    main(write=args.write, redownload=args.redownload, assume_yes=args.yes)
