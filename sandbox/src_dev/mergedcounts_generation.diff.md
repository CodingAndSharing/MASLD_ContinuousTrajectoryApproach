## `mergedcounts_generation.py` — low-count expression filter

**Pipeline stage:** preprocessing, step 0 (builds `data/mergedcounts.csv`, the raw
count matrix every script in `src/01.preprocessing_and_trajectory_analysis/` reads).

**Change class:** two-tier.
`src_dev` ships the filter in two forms — a *reproduce-exact* form that leaves the
published output byte-identical, and an *upgrade* form that changes the retained
gene set. They address different goals and are labelled separately below.

---

### 1. What the original `src/` did

The original script is a provenance/validation utility. Its only row-reducing
operation is the cohort merge:

```python
return ucam.join(vcu, how="inner")   # keep genes present in BOTH cohorts
```

That is a **gene-identity intersection**, not an **expression-level filter** — it
drops a gene only when one cohort never measured it, never because its counts are
low. No committed script performs low-count filtering, even though
`Normalization_and_PCA.R` carries the comment `# ... exclude low expressed counts`.
The filter it refers to was applied in the original, **unstored** HTSeq→merge step
and never entered version control.

### 2. Evidence the filter existed (recovered from the data)

The committed `data/mergedcounts.csv` is **not** raw HTSeq output. Reading the
matrix directly (17,090 genes × 136 samples):

| Property | Value | Interpretation |
|---|---|---|
| Minimum per-gene **mean** count | **1.0074** (= 137/136) | hard floor exactly at mean = 1 |
| Genes with mean count < 1.0 | **0** | nothing below the floor |
| All-zero genes | **0** | raw HTSeq output always has some |
| Genes detected in ≤ 1 sample | **0** | — |

A wall precisely at `mean = 1` with nothing to its left is the fingerprint of a
filter of the form **keep gene iff `rowMeans(counts) ≥ 1`** (equivalently
`rowSums ≥ n_samples = 136`). Raw counts ramp smoothly through zero; this does not.

![Filter signature: hard floor at mean count = 1]({{artifact:art_76f841a3-c903-4700-8d3c-b3f006cc5a60}})

### 3. Change A — reproduce-exact (`rowMeans ≥ 1`)

Makes the recovered filter **explicit and reproducible** without altering the
published output.

- Adds `MIN_MEAN_COUNT = 1.0` and `filter_low_expressed_genes(counts, min_mean)`.
- **Real-merge path:** the filter is applied immediately after concatenation, so a
  matrix regenerated from raw per-cohort counts matches the published gene set.
- **Fallback path** (the everyday run, which reads the committed matrix): the filter
  is re-applied as an **idempotency assertion** — if it ever drops a gene from the
  committed matrix, the script raises, because that would mean the recovered
  threshold no longer matches the shipped data.

**Effect on published output:** none. `kept 17090 / 17090 (dropped 0)`. The change
is documentary — the previously-invisible step is now in code and self-checking.

**Boundary note:** no gene sits exactly at `rowSum = 136`, so `≥` vs `>` is
unobservable from the data; the conventional `≥` form is used.

### 4. Change B — upgrade (`edgeR::filterByExpr`, group = disease stage)

Replaces the absolute-count threshold with a **depth- and design-aware** filter,
implemented in `src_dev/mergedcounts_generation_dev.py`.

**Why it is preferable.** Library sizes in this matrix span a **62-fold range
(0.62M – 38.5M reads)**. `rowMeans ≥ 1` is an unnormalized absolute-count rule: it
treats a gene at the floor identically whether a sample had 0.6M or 38M reads.
`filterByExpr` instead sets a CPM cutoff scaled to the median library size
(`10/median-lib-in-M ≈ 0.38 CPM`) and requires a gene to clear it in at least as
many samples as the **smallest design group** — current best practice ahead of
normalization and modelling.

**Design (`group=`) choice matters — and it is biological.** The number of genes
dropped depends entirely on the group argument (all figures below verified against
**edgeR 4.8.2** on the real 17,090-gene matrix):

| `group=` | MinSampleSize | Genes kept | Genes dropped |
|---|---|---|---|
| **SAF score (disease stage)** ← src_dev default | 4 (CTRL, smallest stage) | **16,883** | **207** |
| none (filterByExpr default) | 98.2 | 14,155 | 2,935 |
| Dataset (cohort) | 43.6 | 15,724 | 1,366 |

`src_dev` uses **`group = SAF score`** deliberately. The smallest stage group
(CTRL, n = 4) sets MinSampleSize = 4, which **preserves stage-specific genes** —
e.g. a gene switched on only in the F3/F4 fibrosis samples — that a trajectory
spanning all stages depends on. The ungrouped default (MinSampleSize ≈ 98) would
force such a gene to be expressed in ~98 of 136 samples and **drop it**, discarding
exactly the late-stage signal of interest; the cohort grouping conflates expression
filtering with the batch correction ComBat performs later.

**Effect on output — this changes the gene set:** 17,090 → **16,883** (−207
relative to the original filter). The dropped genes are near-noise once depth is
accounted for, and still contribute tied low values to the rank-based
quantile-normalization reference distribution.

**Implementation note.** `src_dev` ships a dependency-free reimplementation of
`edgeR::filterByExpr.default` (`filter_by_expr()`); its output is bit-identical to
edgeR 4.8.2 on all three designs above, so the script needs neither R nor edgeR
installed.

### 5. Why the filter matters mathematically

The next step (`Normalization_and_PCA.R`) is **quantile normalization**, which is
**rank-based**: every sample is mapped onto a common reference distribution built as
the mean of order statistics across samples. A large mass of tied near-zero counts
distorts that reference and injects sampling noise into the low-expression ranks.
Removing never-/barely-expressed genes *before* normalization is what keeps it
well-behaved — and, downstream, keeps those noisy genes out of the PC1 axis that is
interpreted as the disease trajectory. Change A makes this pre-condition explicit;
Change B strengthens it to a depth-aware criterion.

---

**Summary**

| | `src/` (original) | `src_dev/` Change A | `src_dev/` Change B |
|---|---|---|---|
| Low-count filter | implicit / not in code | `rowMeans ≥ 1`, explicit | `filterByExpr` / CPM |
| Depth-aware | — | no | yes |
| Genes retained | 17,090 | 17,090 | **16,883** |
| Published output | (baseline) | **unchanged** | **changed** (−207) |
| Purpose | — | reproduce + document | methodological upgrade |
