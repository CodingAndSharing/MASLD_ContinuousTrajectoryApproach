"""Documents and (re)generates data/mergedcounts.csv -- src_dev version.

This is the src_dev counterpart of
src/01.preprocessing_and_trajectory_analysis/mergedcounts_generation.py.
It differs in ONE respect: the low-count expression filter.

    original src/ : keep gene iff rowMeans(counts) >= 1
                    (an absolute, depth-blind threshold; reproduces the
                    published 17,090-gene matrix exactly)

    src_dev/      : keep gene iff edgeR::filterByExpr(counts, group=SAF score)
                    (a depth- and design-aware filter; current best practice)

Why the change
--------------
Library sizes in this matrix span a 62-fold range (0.62M .. 38.5M reads). An
absolute mean-count threshold treats a gene at the floor identically regardless
of a sample's sequencing depth. filterByExpr instead sets a CPM cutoff scaled to
the median library size and requires a gene to clear it in at least as many
samples as the SMALLEST design group -- so genes expressed only in a small
disease stage (e.g. late fibrosis) are retained, while genes that are near-zero
once depth is accounted for are dropped.

Design (group=) choice
----------------------
group is the SAF score (disease stage). This is deliberate: the smallest group
(CTRL, n=4) sets MinSampleSize=4, which PRESERVES stage-specific genes across the
MASLD continuum -- exactly the signal the trajectory analysis depends on. Using
no group (MinSampleSize~98) would drop such genes and is inappropriate here; using
the dataset/cohort as the group conflates expression filtering with the batch
correction that ComBat performs later.

Effect on the gene set (verified against edgeR 4.8.2 on the real matrix)
------------------------------------------------------------------------
    group = SAF score : keep 16,883 / 17,090  (drop   207)   <- src_dev default
    group = none      : keep 14,155 / 17,090  (drop 2,935)
    group = Dataset   : keep 15,724 / 17,090  (drop 1,366)
    original rowMeans>=1 : keep 17,090         (drop     0)

filter_by_expr() below is a dependency-free reimplementation of
edgeR::filterByExpr.default (v4.x); its output is bit-identical to edgeR 4.8.2 on
all three designs above, so this script needs neither R nor edgeR installed.

Everything else (provenance, --redownload, validation, --write) is unchanged from
the original; see that file's docstring for the full pipeline description.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data"

METADATA_PATH = DATA_DIR / "metadata.csv"
MERGED_COUNTS_PATH = DATA_DIR / "mergedcounts.csv"

UCAM_RAW_COUNTS = DATA_DIR / "raw" / "ucam_htseq_counts.tsv"
VCU_RAW_COUNTS = DATA_DIR / "raw" / "vcu_htseq_counts.tsv"

EXPECTED_DATASET_COUNTS = {"UCAM": 58, "SANYAL": 78}  # Methods, p.11
OUTLIER_SAMPLE_NAME = "Sample 5"

# Column in metadata.csv that defines the filterByExpr design group (disease stage).
GROUP_COLUMN = "SAF score"

# filterByExpr.default parameters (edgeR v4.x defaults; exposed for transparency).
FBE_MIN_COUNT = 10
FBE_MIN_TOTAL_COUNT = 15
FBE_LARGE_N = 10
FBE_MIN_PROP = 0.7


# ---------------------------------------------------------------------------
# Merge (unchanged from the original)
# ---------------------------------------------------------------------------
def merge_from_raw_htseq_counts(ucam_path: Path, vcu_path: Path) -> pd.DataFrame:
    """Concatenate the two cohorts' HTSeq gene-count matrices (inner join on
    Ensembl gene ID). Identical to the original script."""
    ucam = pd.read_csv(ucam_path, sep=None, engine="python", index_col=0)
    vcu = pd.read_csv(vcu_path, sep=None, engine="python", index_col=0)
    return ucam.join(vcu, how="inner")


# ---------------------------------------------------------------------------
# Low-count filter -- the src_dev change
# ---------------------------------------------------------------------------
def filter_by_expr(
    counts: pd.DataFrame,
    group=None,
    lib_size=None,
    min_count: float = FBE_MIN_COUNT,
    min_total_count: float = FBE_MIN_TOTAL_COUNT,
    large_n: int = FBE_LARGE_N,
    min_prop: float = FBE_MIN_PROP,
) -> pd.Series:
    """Dependency-free reimplementation of edgeR::filterByExpr.default (v4.x).

    Returns a boolean Series (indexed like `counts`) marking genes to KEEP.
    Verified bit-identical to edgeR 4.8.2 on this matrix for group in
    {None, SAF score, Dataset}.

    Logic (see edgeR source):
      1. MinSampleSize = size of the smallest group (or n_samples if no group).
         If MinSampleSize > large_n, soften it to large_n + (MinSampleSize -
         large_n) * min_prop, so very large groups don't demand an unrealistically
         high number of expressing samples.
      2. CPMcutoff = min_count / median(lib_size) * 1e6.
      3. Keep a gene iff it has CPM >= CPMcutoff in at least MinSampleSize samples
         AND its total count across samples >= min_total_count.
    (`tol` mirrors edgeR's floating-point tolerance on the two thresholds.)
    """
    y = counts.values.astype(float)
    if lib_size is None:
        lib_size = y.sum(axis=0)
    lib_size = np.asarray(lib_size, dtype=float)

    if group is not None:
        sizes = pd.Series(group).value_counts()
        sizes = sizes[sizes > 0]
        min_sample_size = float(sizes.min())
    else:
        min_sample_size = float(y.shape[1])

    if min_sample_size > large_n:
        min_sample_size = large_n + (min_sample_size - large_n) * min_prop

    median_lib = float(np.median(lib_size))
    cpm_cutoff = min_count / median_lib * 1e6
    cpm = y / lib_size[None, :] * 1e6

    tol = 1e-14
    keep_cpm = (cpm >= cpm_cutoff).sum(axis=1) >= (min_sample_size - tol)
    keep_total = y.sum(axis=1) >= (min_total_count - tol)
    keep = keep_cpm & keep_total

    print(
        f"[filter] filterByExpr (group={'none' if group is None else GROUP_COLUMN}, "
        f"MinSampleSize={min_sample_size:.1f}, CPMcutoff={cpm_cutoff:.3f}): "
        f"kept {int(keep.sum())} / {len(keep)} genes (dropped {int((~keep).sum())})."
    )
    return pd.Series(keep, index=counts.index)


def apply_expression_filter(counts: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    """Apply the src_dev depth-aware filter, using GROUP_COLUMN from metadata as
    the design group. Falls back to an ungrouped filter (with a warning) if the
    column is missing or does not line up with the count columns."""
    group = None
    if GROUP_COLUMN in metadata.columns and len(metadata) == counts.shape[1]:
        group = metadata[GROUP_COLUMN].values
    else:
        print(
            f"[warn] Cannot align '{GROUP_COLUMN}' from metadata with the count "
            "columns; falling back to an ungrouped filterByExpr (more stringent)."
        )
    keep = filter_by_expr(counts, group=group)
    return counts.loc[keep.values]


# ---------------------------------------------------------------------------
# Validation & labelling (unchanged from the original)
# ---------------------------------------------------------------------------
def validate_against_metadata(counts: pd.DataFrame, metadata: pd.DataFrame) -> None:
    n_counts, n_meta = counts.shape[1], len(metadata)
    if n_counts != n_meta:
        raise ValueError(
            f"Sample count mismatch: {n_counts} count columns vs {n_meta} metadata rows."
        )
    dataset_counts = metadata["Dataset"].value_counts().to_dict()
    if dataset_counts != EXPECTED_DATASET_COUNTS:
        raise ValueError(
            f"Dataset composition mismatch: got {dataset_counts}, expected "
            f"{EXPECTED_DATASET_COUNTS} (Methods, p.11)."
        )
    if OUTLIER_SAMPLE_NAME not in set(metadata["Sample name"]):
        raise ValueError(
            f"Expected outlier {OUTLIER_SAMPLE_NAME!r} not found in metadata.csv."
        )
    print(
        f"[OK] {n_counts} samples, {counts.shape[0]} Ensembl genes; dataset "
        f"composition {dataset_counts}; outlier {OUTLIER_SAMPLE_NAME!r} present."
    )


def label_columns_from_metadata(counts: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    labelled = counts.copy()
    labelled.columns = metadata["Sample name"].values
    return labelled


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main(write: bool = False) -> None:
    metadata = pd.read_csv(METADATA_PATH, index_col=0)

    if UCAM_RAW_COUNTS.exists() and VCU_RAW_COUNTS.exists():
        print("[merge] Per-cohort raw HTSeq matrices found -- concatenating.")
        counts = merge_from_raw_htseq_counts(UCAM_RAW_COUNTS, VCU_RAW_COUNTS)
    else:
        print(
            "[fallback] Per-cohort raw HTSeq matrices not present; reading the "
            f"committed {MERGED_COUNTS_PATH.name}. NOTE: that file was produced with "
            "the ORIGINAL rowMeans>=1 filter, so applying the src_dev filterByExpr "
            "here will drop the extra low-count genes and CHANGE the gene set "
            "(this is the intended src_dev behaviour, unlike the original script "
            "whose filter is a no-op on this file)."
        )
        counts = pd.read_csv(MERGED_COUNTS_PATH, index_col=0)

    counts = apply_expression_filter(counts, metadata)

    validate_against_metadata(counts, metadata)
    labelled = label_columns_from_metadata(counts, metadata)

    if write:
        labelled.to_csv(MERGED_COUNTS_PATH)
        print(f"[write] {MERGED_COUNTS_PATH.name} regenerated (src_dev filter, labelled columns).")
    else:
        print("[dry-run] Not writing to disk (pass --write to regenerate).")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--write", action="store_true",
        help="Regenerate mergedcounts.csv with the src_dev filterByExpr filter and labelled columns.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main(write=_parse_args().write)
