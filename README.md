# Code and data repository

This repository contains the analysis code, gene sets and parameter files for:

**Single-cell and spatial transcriptomics reveal an immune-metabolic decoupled state carried by
plasma cells and fibroblasts in ovarian cancer**

## Layout

```
code_repository/
├── 00_gene_sets/            M1-M4, ICC45_full/shared28/unique17, random-set generators
├── 01_bulk/                 bulk cohorts, deconvolution, survival, meta-analysis
├── 02_single_cell/          scoring, carriers, CNV, trajectory, TF perturbation, drug sensitivity
├── 03_spatial/              Visium processing, NMF domains, NNLS, Moran, COMMOT
├── 04_revision/             new analyses added in revision (null models, UCell, depth adjustment,
│                            nested CV, pseudotime robustness, clinical meta-analysis)
├── 05_figures/              figure-generation scripts
├── environment/             Python and R package versions
└── README.md
```

## Order of execution

1. `00_gene_sets/build_gene_sets.R` – builds `mitoxy_sets.RDS` from the published mechanism genes.
2. `01_bulk/00_track0_tcga.R` – TCGA-OV module scores, state calls, survival.
3. `01_bulk/03_external_survival_meta.R` – external cohorts and meta-analysis.
4. `02_single_cell/20a_score_gse266577.R` – primary reference module scores.
5. `02_single_cell/44*`–`47*` – figures and carrier analyses.
6. `03_spatial/28*` – Visium, NMF, NNLS, Moran, COMMOT.
7. `04_revision/A1*`–`A7*` – null models, sensitivity, technical adjustment, UCell, nested CV,
   pseudotime robustness, meta-analysis (new in revision).

## Key parameters

| Parameter | Value |
|---|---|
| State cutoffs | M1 ≥ P75; M2 ≤ P25; M3 ≤ P25 (within source) |
| AddModuleScore | nbin = 24, control = 100 |
| UCell | maxRank = 1500, ties.method = "average" |
| Permutations | state 10,000; Moran 999; COMMOT 1,000; LR 1,000; TF 2,000; drug 1,000 |
| Random seeds | 42 (analysis), 20261105 (revision null models) |
| NNLS bootstrap | 1,000 |
| Pseudotime bins | 7 equal-frequency bins |
| Cox effect scale | per 1 SD of dec_res (SD = 0.956) |

## Data

All data are public (GEO accessions in Table 1). Processed module-score tables, state calls,
NNLS bootstraps, COMMOT permutations and meta-analysis inputs are provided in `data_processed/`
or available from the corresponding author.

## Software

Python 3.12 (scanpy, scikit-learn, scipy, anndata, scikit-survival, lifelines, statsmodels,
commot, matplotlib) and R 4.4.3 (Seurat 5.5.0, UCell 2.10.1, lme4, metafor, survival, GSVA,
data.table, Matrix). Exact versions are in `environment/`.
