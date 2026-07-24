# MASLD_Continuous-trajectory-approach

By moving beyond conventional [stage-based MASLD classification](https://www.aasld.org/sites/default/files/2024-01/SLD-597%20AASLD%202023%20MASLD%20Decision%20Tree%20Digital%20%28sm%291.8.pdf), we uncover the sequence of **critical molecular events**, enhancing our understanding of MASLD/MASH pathophysiology. 
This approach enabled the **identification of novel trajectory-specific biomarkers**, offering a more refined and personalised strategy for managing MASLD patients.

![MASDL progression](https://www.mdpi.com/ijms/ijms-25-05640/article_deploy/html/images/ijms-25-05640-g001-550.jpg)

Refs: 
1.-[A data-driven framework reconstructs the molecular continuum of human MASLD progression](https://www.nature.com/articles/s42255-026-01543-7)
2.-[Progression of hepatic pathology](https://cdn.ncbi.nlm.nih.gov/pmc/blobs/57f4/12150358/41ee7ee629ab/DOM-27-3-g004.jpg)



Food for thought:

_How much population is affected by this?_
_When is needed clinical intervention?_
_Who will benefit of a therapy (stage)?_
_Are the Biomarkers/critical targets targeatable?_
_Which other potential diseases are we avoiding by targeting the MASDL?_


## Reproducibility

The original repository is archived on Zenodo for long-term reproducibility:
https://doi.org/10.5281/zenodo.19816413
This corresponds to GitHub release v1.0.


## Code structure

- The different **parts of the analysis** have been organised in specific files and/or sub-folders in the main 'src' folder. 
- The full analysis is organised **into sequential steps, numbered (0-12)** according to the order they should be executed (following the order of the results/figures/tables presented in the manuscript). 
- Each step corresponds to a specific part of the analysis **pipeline**, and can be run independently of the rest (using data generated in the previous steps).


0. **download_and_merge_data**: Documents and (re)generates `data/mergedcounts.csv` and `data/metadata.csv` -- the raw, un-normalised UCAM + VCU/Sanyal gene-count matrix and per-sample clinical/histopathology table every downstream script starts from. `mergedcounts_generation.py` validates the committed matrix against `metadata.csv` (and can re-merge it from separate per-cohort HTSeq matrices, if you have them); `metadata_generation.py` cross-checks the committed metadata against GEO/ArrayExpress' public records, since large parts of it (UCAM's histology component scores, T2DM/SAF-category for both cohorts) are human-curated and not otherwise reproducible from public sources. `regenerate_from_raw_fastq.py` goes further and does a real, from-scratch regeneration -- downloading the actual raw paired-end FASTQ for both cohorts and running FastQC -> HISAT2 -> HTSeq itself -- at the cost of ~338 GB of downloads and, per an unbenchmarked estimate, several days of alignment compute; see its module docstring before running it. `sandbox/notebooks/00_download_and_create_data.ipynb` walks through running all three against a disposable `sandbox/data_dev/` copy, without ever touching `data/` itself.

1. **preprocessing_and_trajectory_analysis**: File organisation, preprocessing (normalisation, PCA, batch correction) and main trajectory inference and analysis.
**Validation:** The following subfolders include the scripts used for validation of the preprocessing and **trajectory inference** part of the analysis. 
**a)Linking_disease_state_to_phenotypes**: Script producing the plots showing agreement of progressing disease characteristics with trajectory position, 
**b)Sex_Differences**: Script producing the correlation plots (disease and general population characteristics against the first 10 principal components) separately for males and females, 
**c)Error_bars_in_histology_per_SW**: Script producing the plots of the average histological scores per SW (including the error bars in each SW), 
**d)Correlation_analysis_normalised_vs_batchcorrected_counts**: Script producing all the supplementary plots validating successful batch effect correction, including correlations per disease stage and stratified correlation for all disease stages.


2. **finding_of_optimal_SWs_sequence**: Contains the scripts used to find the optimal set of Sliding Windows for patients' stratification in the pseudo-temporal space.

3. **wgcna_and_linear_modelling**: Contains the scripts used to identify the gene co-expression modules (running WGCNA) and to model the MASLD phenotypes using the calculated eigen-genes.

4. **deseq_and_msviper**: Contains the scripts that execute the differential expression and transcription factor analysis for both discrete and SW-based patients' stratification.
5. **enrichment_analysis**: Contains the functions used to perform TF- and pathway (Reactome) enrichment analysis in order to connect the gene co-expression modules with TFs and upstream signalling pathways.
6. **create_MASLD_network**: Contains all the steps used to create the MASLD network which constitute the main framework used for the downstream analysis. Firstly, a reference network is created based on available signalling and metabolic pathway databases. Then the individual upstream and downstream networks that are associated with the MASLD variables are constructed. Finally, these distinct networks are unified into the global MASLD network, whose edges are annotated with the functional similarities of interacting pairs.
7. **network_analysis**: Includes the functions used for network propagation on the MASLD network in order to extract the "up" and "down" signatures for each disease stage, as well as the interpretation of these signatures using the Reactome pathways database.
8. **cell_type_deconvolution**: Contains a simple script to aggregate the results from cell type deconvolution and calculate an average composition for each cell type.
9. **pseudo_bulk_analysis**: This is a worklfow of different tasks to generate pseudo bulk RNA-seq samples for each disease state, given the respective cell type proportions and then perform differential and network analysis to identify molecular changes along the disease trajectory. It is assumed that these changes have been derived only from the differentiations in cell type proportions. Finally the part of differentiated pathways that cannot be explained by this factor is defined as the 'deregulated' disease component and is associated with the different cell types to pinpoint disease-relevant cell type-specific processes.
10. **markers_analysis**: It contains scripts to retrieve data from ChEMBL and relate the 57 plasma markers (found in the MASLD network) with drug-target indications. Deregulated targets in the network at any stage of the disease are annotated with drugs so as to reverse the direction of their regulation. Data retrieval from ChEMBL should be performed using Python versions <= 3.8, due to the compatibility of the api service that is used.
11. **RandomForest_analysis**: Includes all Random Forest classification and regression analyses, applied to produce the biomarker sets described in the study. The training has been applied exclusively on the UCAM/VCU dataset. Additionally, we provide all the files that reproduce the analyses on the external validation transcriptomic datasets (EPoS, GUBRA, Fujiwara) as well as the validation on the plasma proteomic datasets.
12. **Genetic_evidence_57biomarkers**: Contains the script that produces all files and figures related to the GWAS catalog analysis including the gene/biomarker-trait category associations, and enrichment of MASLD biomarker genes on the trait categories (e.g., metabolic, liver, cardiovascular, inflammatory etc), as well as summary gwas files and associated genes/biomarkers with multiple gwas categories.

*Detailed comments are included throughout the scripts to guide users in reproducing all the results of our analysis.

## Data Availability

The analysis uses multiple human transcriptomic and proteomic datasets. All datasets used in our study are fully provided in the 'data/' directory of this repository.
All data are publicly available and have been included here to ensure full reproducibility of the analysis. No special permissions are required to access or use these datasets.
All the results and figures presented in the manuscript can be reproduced directly using the data and code provided in this repository. 

More specifically, the datasets are:
- UCAM/VCU dataset
- EPoS dataset
- GUBRA dataset
- Fujiwara dataset
- Plasma proteomics dataset



## Software requirements


#### Python libraries:

Python3: 3.12.3

- chembl_webresource_client: 0.10.9 (Python 3.8 environment)
- igraph: 0.11.4
- json: 2.0.9
- numpy: 1.26.4
- pandas: 2.1.4
- scipy: 1.11.4
- statmodels: 0.14.1



#### R libraries:

R: 4.4.1

- assertr: 3.0.1
- BiocParallel: 1.40.2
- biomaRt: 2.54.0
- broom: 1.0.10
- caret: 6.0-94
- cowplot: 1.1.1
- data.table: 1.17.8
- decoupleR: 2.12.0
- DESeq2: 1.46.0
- dplyr: 1.1.4
- flashClust: 1.1.2
- ggplot2: 4.0.1
- ggpubr: 0.6.0
- GO.db: 3.20.0
- GOSemSim: 2.29.2
- gridExtra: 2.3
- igraph: 2.2.1
- jsonlite: 2.0.0
- modelr: 0.1.11
- numbers: 0.8.5
- OmnipathR: 3.8.0
- PCSF: 0.99.1
- pheatmap: 1.0.12
- plotly: 4.11.0
- preprocessCore: 1.68.0
- pROC: 1.18.5
- purrr: 1.2.0
- randomForest: 4.6.14
- Rcpp: 1.1.0
- readr: 2.1.6
- readxl: 1.4.5
- reshape2: 1.4.4
- rlang: 1.1.6
- Seurat: 5.3.1
- SingleCellExperiment: 1.28.1
- slingshot: 3.22
- stats: 4.0.3
- stringr: 1.6.0
- sva: 3.38.0
- tidyr: 1.3.1
- tidyverse: 2.0.0
- viper: 1.40.0
- WGCNA: 1.73.0


## Sandbox pixi environment (Python + R, Jupyter)

`sandbox/pixi.toml` + `sandbox/pixi.lock` reproduce the stack above (conda-forge + bioconda, every dependency pinned with `==` to a version validated to solve together) -- no separate conda/mamba install, everything goes through the `pixi` CLI.

### Setup

```bash
cd sandbox
pixi install                  # creates/updates .pixi/envs/default from the committed lock file
pixi run python --version     # 3.12.3
pixi run Rscript --version    # R 4.4.1
```

Two R packages have no conda-forge/bioconda build and install from GitHub via `remotes` instead (see the `[tasks]` in `pixi.toml`):

```bash
pixi run install-pcsf         # IOR-Bioinformatics/PCSF -> 0.99.1
pixi run install-omnipathr    # saezlab/OmnipathR -> current (4.x; requested 3.8.0 is pre-rewrite
                               # and incompatible with current OmnipathR/Bioconductor APIs)
```

Also run this once, so the two Jupyter kernels below show up in the kernel picker with clear, unambiguous names instead of ipykernel's generic default (see below):

```bash
pixi run setup-jupyter-kernels
```

`sandbox/.pixi/config.toml` sets `run-post-link-scripts = "insecure"`, scoped to this project only -- required for `GO.db`, `org.Hs.eg.db`, and `bioconductor-genomeinfodbdata`, which ship as stub packages whose post-link script downloads the real annotation database at install time.

### Using this environment as the Jupyter kernel for `sandbox/notebooks/`

A single Jupyter kernel is inherently one language -- there is no kernel that natively runs both raw R and raw Python syntax in the same cell. This pixi env registers **two separate kernels** (both installed by `pixi install`, both named clearly by `pixi run setup-jupyter-kernels` above), and which one to pick depends on the notebook:

- **`00_download_and_create_data.ipynb` is pure Python** -- pick **"Python 3 (MASLD sandbox pixi)"**. Picking the R kernel here fails immediately (R can't run `import subprocess` etc.) -- this is the most likely cause if a notebook "doesn't run" right after a kernel switch.
- **`01_preprocessing_and_trajectory_analysis.ipynb` is R-heavy** (Seurat/slingshot/DESeq2) -- pick **"R (MASLD sandbox pixi)"**.

Before `setup-jupyter-kernels` was added, ipykernel's default display name was the generic **"Python 3 (ipykernel)"** -- indistinguishable from any other Python kernel on the machine, while the R one already showed up as "R (MASLD sandbox pixi)". That asymmetry is exactly what makes it easy to end up on the wrong kernel (or the R one by mistake) when looking for "the pixi environment's kernel" in the picker; re-run the task any time that ambiguity resurfaces (e.g. after wiping `.pixi/envs` and reinstalling).

```bash
cd sandbox
pixi run jupyter lab notebooks/
```

In the kernel picker choose the display name matching the notebook (see above). In VS Code: open the notebook, click the kernel selector (top right), choose "Select Another Kernel" -> "Jupyter Kernel..." and pick the same one -- VS Code discovers both because `pixi run jupyter lab` (or any `jupyter`/`python`/`Rscript` command run through `pixi run`/`pixi shell`) exposes the kernelspecs under `.pixi/envs/default/share/jupyter/kernels/`.

If VS Code's kernel picker does not show either Pixi-backed kernel at all:

1. Open the folder `sandbox/` in VS Code (not just the notebook file).
2. Press `Ctrl+Shift+P` (or `Cmd+Shift+P` on macOS) and run `Python: Select Interpreter`.
3. Choose the Pixi interpreter at `.pixi/envs/default/bin/python`.
4. Open your notebook and click the kernel selector in the top-right corner.
5. Choose the Pixi-backed kernel matching the notebook (see above).

If it still does not list it, run `pixi run jupyter lab notebooks/` once (as above) and reopen the notebook from that Jupyter session or reload VS Code -- this makes the kernelspecs under `.pixi/envs/default/share/jupyter/kernels/` visible to the editor.

If you ever wipe `.pixi/envs` and reinstall, the `ir` kernelspec may also need its R path hand-patched to the env's absolute `R` binary (`.pixi/envs/default/lib/R/bin/R`) so it still resolves if launched by a process without the pixi env active on `PATH` -- or just always launch Jupyter via `pixi run`, which works with a bare `R` too and is the common case.

**Diagnosing "wrong kernel" symptoms**: if a Python cell fails with something like

```
Error in parse(text = input): <text>:1:8: unexpected symbol
1: import importlib
```

that `Error in parse(text = input)` is R's own parser, not this environment -- the notebook is connected to the R kernel, not the Python one. Switch it via the kernel selector. Each `.ipynb`'s `metadata.kernelspec` records which kernel it should default to; `00_download_and_create_data.ipynb` has `kernelspec.name = "python3"` set explicitly for this reason. A notebook with no `kernelspec` at all (only a bare `language_info` hint) has nothing telling the editor which kernel to default to, and can end up silently reusing whatever was last manually selected -- if a notebook keeps reopening on the wrong kernel, check whether its metadata is missing this block.

**Mixing both languages in one notebook**: `%%R` cell magic (rpy2) -- keep the Python kernel, run individual cells as R, pass data back and forth. This is the actual way to get R and Python "in the same kernel session", since a kernel itself can't be bilingual:

```python
%load_ext rpy2.ipython
```
```python
%%R
library(Seurat)
print(R.version.string)
```

Pass Python -> R and back with `%%R -i my_python_df -o result_df`.

Both were tested directly against the Jupyter kernel protocol (not just `Rscript`) and against `rpy2.robjects` -- both execute correctly.
