# PBMC IFN-beta integration with Harmony

This Seurat workflow uses the IFNB-stimulated and control PBMC example distributed through SeuratData. It covers quality control, normalization, exploratory Harmony alignment, clustering, donor-level differential expression, and visualization.

## Data and environment

Create the Conda environment:

~~~bash
conda env create -f environment.yml
conda activate scRNA_integration
~~~

Install the SeuratData package and its IFNB dataset once:

~~~r
install.packages("remotes")
remotes::install_github("satijalab/seurat-data")
SeuratData::InstallData("ifnb")
~~~

The dataset includes a donor_id field. The DE step uses only donors represented in both CTRL and STIM for each cell type, aggregates raw counts by donor and condition, and fits the model ~ donor_id + stim. It skips a cell type if fewer than two complete donors remain.

## Run

Run from the repository root, in order:

~~~bash
Rscript scripts/01_load_and_setup.R
Rscript scripts/02_qc_and_filtering.R
Rscript scripts/03_integration_harmony.R
Rscript scripts/04_celltype_DE_IFNB.R
Rscript scripts/05_visualization.R
~~~

The pseudobulk results are written to results/de_celltype/IFNB_pseudobulk_DESeq2_by_celltype.csv and summarized in IFNB_pseudobulk_summary.csv. The visualization step writes the ISG dot plot and, when CD14 Mono results are available, a donor-level volcano plot.

## Interpretation

Harmony is run against stimulation status for exploratory alignment. Since this can remove treatment-associated signal, differential expression uses the uncorrected RNA count layer, not the Harmony embedding. The IFNB example has donor metadata, so cells are aggregated within donors before DE; cells are not treated as independent replicates.

Existing plots and tables are historical snapshots. Regenerate them with the current scripts before reporting their values. For the dataset and recommended sample-level DE approach, see the [Seurat IFNB integration vignette](https://satijalab.org/seurat/articles/integration_introduction) and [Seurat differential-expression vignette](https://satijalab.org/seurat/articles/de_vignette).

Code in this repository is licensed under the MIT License (see `LICENSE`). The SeuratData dataset and third-party software remain subject to their respective source terms.
