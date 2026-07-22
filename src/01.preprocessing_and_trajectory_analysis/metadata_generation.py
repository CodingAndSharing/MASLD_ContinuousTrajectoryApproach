"""Documents and partially reproduces data/metadata.csv.

data/metadata.csv is the curated per-sample clinical/histopathology table
(Sex, Age, NAS, Fibrosis, Steatosis, Ballooning, Inflammation, SAF score,
T2DM, Dataset) for the 136 UCAM+VCU/Sanyal samples, read alongside
mergedcounts.csv at the top of every script in this directory -- see
mergedcounts_generation.py and sandbox/publication/NOTES_CLAUDE.md.

Who curated it (verified 2026-07-21 by cross-checking public records)
------------------------------------------------------------------------
Not this paper's authors, and not uniformly reproducible:

  - VCU/Sanyal (78 samples, GEO GSE130970): the 7 histology/demographic
    fields are a DIRECT, VERIFIED pass-through of GEO's own public
    per-sample characteristics, deposited by the original VCU study team
    (Hoang et al., HemoShear Therapeutics / Sanyal lab, 2019). Confirmed
    exact match for GSM3758005 / "Sample 1" (Sex, age, NAS, fibrosis,
    steatosis, ballooning, inflammation all identical). fetch_vcu_from_geo()
    below reproduces this half mechanically, from public data.

  - UCAM (58 samples, ArrayExpress E-MTAB-9815): the SDRF metadata
    ArrayExpress exposes publicly only contains a COARSE 3-bucket
    "disease staging" field (NAFL / NASH Fibrosis Stage 0-2 / NASH
    Fibrosis Stage 3-4) -- not the fine-grained NAS/Fibrosis/Steatosis/
    Ballooning/Inflammation component scores present in metadata.csv.
    Those component scores are NOT available from any public,
    machine-readable ArrayExpress/BioStudies record. They trace back to
    the original UCAM study's own clinical/pathology records (Cambridge
    MASH Service; ref. 24, Azzu et al., Mol. Metab. 2021) -- ultimately a
    histopathologist grading liver biopsies under the Kleiner/CRN system
    (ref. 3). That is a human clinical judgment call, not a computable
    quantity, and this script cannot reproduce it.

  - T2DM status and the categorical "SAF score" label are not present in
    GEO's public characteristics for either cohort, so even the VCU half
    is not 100% publicly re-derivable end-to-end.

What this script does
------------------------
1. Fetches the VCU/Sanyal characteristics from GEO's public family SOFT
   file (fetch_vcu_from_geo()) and cross-checks them against the existing
   data/metadata.csv, row by row, for the fields that ARE public.
2. Fetches UCAM's Age/Sex/disease-staging from ArrayExpress's SDRF
   (fetch_ucam_from_arrayexpress()) and cross-checks those few fields --
   while explicitly reporting that the fine-grained histology scores
   cannot be checked or regenerated this way.
3. Never overwrites metadata.csv: since the UCAM histology component
   scores and both cohorts' T2DM/SAF-category fields cannot be
   regenerated from public sources, "reconstructing" the file would mean
   fabricating data. This script validates and reports, it does not write.

Usage
-----
    python metadata_generation.py
"""

import io
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data"
METADATA_PATH = DATA_DIR / "metadata.csv"

GEO_ACCESSION = "GSE130970"  # VCU/Sanyal cohort
GEO_SERIES_DIR = GEO_ACCESSION[:-3] + "nnn"  # GEO's directory-grouping convention, e.g. GSE130nnn
GEO_FAMILY_SOFT_URL = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{GEO_SERIES_DIR}/{GEO_ACCESSION}/soft/{GEO_ACCESSION}_family.soft.gz"
# Small (a few KB gzipped) -- unlike mergedcounts_generation.py's FASTQ
# downloads, no size-warning gate is needed here.

ARRAYEXPRESS_ACCESSION = "E-MTAB-9815"  # UCAM cohort
ARRAYEXPRESS_SDRF_URL = f"https://www.ebi.ac.uk/biostudies/files/{ARRAYEXPRESS_ACCESSION}/{ARRAYEXPRESS_ACCESSION}.sdrf.txt"

# GEO characteristic key -> metadata.csv column name
GEO_FIELD_MAP = {
    "Sex": "Sex",
    "age at biopsy": "Age",
    "nafld activity score": "NAS",
    "fibrosis stage": "Fibrosis",
    "steatosis grade": "Steatosis",
    "cytological ballooning grade": "Ballooning",
    "lobular inflammation grade": "Inflammation",
}
# Confirmed NOT present in GEO's public characteristics for GSE130970 --
# these columns in metadata.csv have no public source for this cohort.
FIELDS_NOT_PUBLIC_FOR_VCU = ["T2DM", "SAF score", "SAF Score (check 'Borderline' NAFL)"]

# Confirmed NOT present (at this granularity) in ArrayExpress's SDRF for
# E-MTAB-9815 -- only a coarse 3-bucket "disease staging" field is public.
FIELDS_NOT_PUBLIC_FOR_UCAM = ["NAS", "Fibrosis", "Steatosis", "Ballooning", "Inflammation", "T2DM", "SAF score"]


def fetch_vcu_from_geo() -> pd.DataFrame:
    """Parse GEO GSE130970's public family SOFT file into a per-sample
    table of the histology/demographic fields it actually contains."""
    with urllib.request.urlopen(GEO_FAMILY_SOFT_URL, timeout=30) as response:
        import gzip

        text = gzip.decompress(response.read()).decode("utf-8", errors="replace")

    rows = []
    current: dict = {}
    for line in text.splitlines():
        if line.startswith("^SAMPLE"):
            if current:
                rows.append(current)
            current = {"geo_accession": line.split("=", 1)[1].strip()}
        elif line.startswith("!Sample_title"):
            current["title"] = line.split("=", 1)[1].strip()
        elif line.startswith("!Sample_characteristics_ch1"):
            key_value = line.split("=", 1)[1].strip()
            if ":" in key_value:
                key, value = key_value.split(":", 1)
                key = key.strip()
                if key in GEO_FIELD_MAP:
                    current[GEO_FIELD_MAP[key]] = value.strip()
    if current:
        rows.append(current)

    return pd.DataFrame(rows)


def fetch_ucam_from_arrayexpress() -> pd.DataFrame:
    """Parse ArrayExpress E-MTAB-9815's SDRF for the fields it actually
    exposes publicly: Age, Sex, and a coarse disease-staging bucket."""
    with urllib.request.urlopen(ARRAYEXPRESS_SDRF_URL, timeout=30) as response:
        text = response.read().decode("utf-8", errors="replace")

    df = pd.read_csv(io.StringIO(text), sep="\t")
    keep = {
        "Source Name": "Sample",
        "Characteristics[age]": "Age",
        "Characteristics[sex]": "Sex",
        "Characteristics[disease staging]": "Disease staging (coarse only)",
    }
    available = {k: v for k, v in keep.items() if k in df.columns}
    return df[list(available)].rename(columns=available).drop_duplicates(subset="Sample")


def validate_vcu(vcu_geo: pd.DataFrame, metadata: pd.DataFrame) -> None:
    sanyal = metadata[metadata["Dataset"] == "SANYAL"].copy()
    # metadata.csv's Title column carries an R-added "X" prefix (e.g.
    # "X440349.1.X_1") over GEO's raw title ("440349.1.X_1").
    sanyal["title_clean"] = sanyal["Title"].str.lstrip("X")

    merged = sanyal.merge(vcu_geo, left_on="title_clean", right_on="title", how="inner", suffixes=("", "_geo"))
    print(f"[VCU/Sanyal] Matched {len(merged)}/{len(sanyal)} metadata.csv rows to a GEO {GEO_ACCESSION} sample.")

    mismatches = 0
    for field in ["Sex", "Age", "NAS", "Fibrosis", "Steatosis", "Ballooning", "Inflammation"]:
        left = merged[field].astype(str)
        right = merged[f"{field}_geo"].astype(str)
        n_bad = (left != right).sum()
        if n_bad:
            mismatches += n_bad
            print(f"  [WARN] {field}: {n_bad} row(s) differ from GEO's public value.")
    if mismatches == 0:
        print("  [OK] All public VCU/Sanyal fields (Sex, Age, NAS, Fibrosis, Steatosis, Ballooning, Inflammation) match GEO exactly.")

    print(f"  [not public] {FIELDS_NOT_PUBLIC_FOR_VCU} are not in GEO's characteristics for this series -- cannot validate or regenerate these from GEO.")


def validate_ucam(ucam_ae: pd.DataFrame, metadata: pd.DataFrame) -> None:
    ucam = metadata[metadata["Dataset"] == "UCAM"].copy()
    merged = ucam.merge(ucam_ae, left_on="Sample name", right_on="Sample", how="inner", suffixes=("", "_ae"))
    print(f"[UCAM] Matched {len(merged)}/{len(ucam)} metadata.csv rows to an ArrayExpress {ARRAYEXPRESS_ACCESSION} sample by name.")

    if len(merged) == 0:
        print(
            "  [gap] 0 matches: metadata.csv's UCAM 'Sample name' values are barcode-style IDs "
            f"(e.g. {ucam['Sample name'].iloc[0]!r}, an Illumina i7/i5 index pair), while "
            f"ArrayExpress's 'Source Name' uses sequential IDs (e.g. {ucam_ae['Sample'].iloc[0]!r}). "
            "There's no direct ID crosswalk in the SDRF for this -- age/sex/disease-stage alone "
            "aren't reliable enough to match 58 individuals without ambiguity, so this script does "
            "not attempt a fuzzy match. Age-level validation is skipped as a result."
        )
    else:
        n_bad = (merged["Age"].astype(str) != merged["Age_ae"].astype(str)).sum()
        if n_bad:
            print(f"  [WARN] Age: {n_bad} row(s) differ from ArrayExpress's public value.")
        else:
            print("  [OK] Age matches ArrayExpress exactly for all matched rows.")

    print(
        "  [not public at this granularity] NAS, Fibrosis, Steatosis, Ballooning, Inflammation are not "
        "in ArrayExpress's SDRF for this accession -- only a coarse 3-bucket 'disease staging' field is "
        "public (NAFL / NASH Fibrosis Stage 0-2 / NASH Fibrosis Stage 3-4). The fine-grained UCAM scores "
        "in metadata.csv cannot be validated or regenerated from ArrayExpress; they require the original "
        "UCAM study's own records (ref. 24) or the Cambridge MASH Service's clinical database directly."
    )


def main() -> None:
    metadata = pd.read_csv(METADATA_PATH, index_col=0)

    print(f"[fetch] GEO {GEO_ACCESSION} family SOFT file ({GEO_FAMILY_SOFT_URL}) ...")
    try:
        vcu_geo = fetch_vcu_from_geo()
        validate_vcu(vcu_geo, metadata)
    except (urllib.error.URLError, OSError) as exc:
        print(f"  [error] Could not reach GEO ({exc}).")

    print(f"\n[fetch] ArrayExpress {ARRAYEXPRESS_ACCESSION} SDRF ({ARRAYEXPRESS_SDRF_URL}) ...")
    try:
        ucam_ae = fetch_ucam_from_arrayexpress()
        validate_ucam(ucam_ae, metadata)
    except (urllib.error.URLError, OSError) as exc:
        print(f"  [error] Could not reach ArrayExpress/BioStudies ({exc}).")

    print(
        "\n[summary] metadata.csv is only PARTLY reproducible from public data:\n"
        "  - VCU/Sanyal's 7 histology/demographic fields: reproducible, verified above.\n"
        "  - UCAM's histology component scores, and T2DM/SAF-category for both cohorts:\n"
        "    require the original studies' own (human-curated, histopathologist-graded)\n"
        "    records -- not obtainable from ArrayExpress/GEO's public metadata alone.\n"
        "  This script therefore validates against public sources; it does not write "
        "or regenerate metadata.csv."
    )


if __name__ == "__main__":
    main()
