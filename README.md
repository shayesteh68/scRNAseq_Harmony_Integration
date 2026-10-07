# scRNA-seq Integration with Harmony - Human PBMC IFN-beta Dataset

**Condition-aware single-cell RNA-seq analysis in R (Seurat v5 + Harmony): integration of a two-condition PBMC dataset, graph-based clustering, UMAP, marker-based cell-type annotation, and cell-type-specific differential expression.**

![R](https://img.shields.io/badge/R-4.6.1-blue.svg)
![Seurat](https://img.shields.io/badge/Seurat-v5-green.svg)
![Harmony](https://img.shields.io/badge/Integration-Harmony-orange.svg)
![License](https://img.shields.io/badge/License-MIT-yellow.svg)

---

## Overview

This repository contains a reproducible single-cell RNA-seq (scRNA-seq) workflow implemented in **R** with **Seurat v5**, applied to the 10x Genomics **PBMC IFN-beta** dataset. The dataset contains two conditions - a **control (CTRL)** group and an **interferon-beta stimulated (STIM)** group - generated in separate batches.

When both conditions are analyzed together without correction, cells can separate in the embedding space for technical as well as biological reasons, and unsupervised clustering may follow the condition rather than the cell identity. This workflow therefore integrates the data with **Harmony**, which projects the PCA embedding into a shared space while correcting the specified condition/batch covariate, and only then performs clustering, visualization, and downstream statistics.

The practical goal of the project is to show a clean, methodologically defensible treatment of a central scRNA-seq question: **does the clustering survive integration, and does it still capture real cell types?**

---

## Pipeline Workflow

1. **Data loading, quality control and setup** - loads the IFN-beta PBMC dataset, applies QC filtering, and writes the filtered Seurat object to `data/processed/ifnb_filtered.rds`.
2. **Normalization and linear dimensionality reduction** - log-normalization, selection of highly variable genes, scaling, and PCA on the scaled data.
3. **Harmony integration and clustering** - `RunHarmony()` is applied to the PCA embedding using the **condition/batch variable** as the grouping covariate; the SNN graph, graph-based clusters, and UMAP are computed **on the Harmony embedding**; the integrated object is written to `data/processed/ifnb_harmony_integrated.rds`.
4. **Cell-type annotation** - clusters are annotated from canonical marker genes and visualized on the integrated UMAP.
5. **Cell-type-specific differential expression** - STIM vs. CTRL differential expression within each annotated cell type, summarized in `DE_summary_top_genes.csv` (top up- and down-regulated genes per cell type).

The analysis order is described above; the exact script files that implement it are listed in the **Scripts** section below, taken directly from the repository contents.

---

## Key Outputs

- `data/processed/ifnb_filtered.rds` - quality-controlled Seurat object
- `data/processed/ifnb_harmony_integrated.rds` - Harmony-integrated, clustered and annotated Seurat object
- `DE_summary_top_genes.csv` - top up- and down-regulated genes per cell type (STIM vs. CTRL). The top signature is dominated by interferon-stimulated genes (ISGs), the expected biological response to IFN-beta stimulation.

---

## Scripts
```text
scripts/01_load_and_setup.R
scripts/02_qc_and_filtering.R
scripts/03_integration_harmony.R
scripts/04_celltype_DE_IFNB.R
scripts/05_visualization.R

## Repository Structure

text
README.docx
README.md
data
data/processed
data/raw
environment.yml
figures
figures/de_celltype
figures/integration
figures/qc
fix_readme.sh
results
results/de_celltype
scripts
scripts/01_load_and_setup.R
scripts/02_qc_and_filtering.R
scripts/03_integration_harmony.R
scripts/04_celltype_DE_IFNB.R
scripts/05_visualization.R

## Requirements

- R 4.6.1
- Seurat v5, SeuratObject, harmony, dplyr, ggplot2
- SeuratData (source of the IFNB dataset)

---

## Author

**Narges Shayesteh**
- GitHub: [@shayesteh68](https://github.com/shayesteh68)
- LinkedIn: [linkedin.com/in/narges-shayesteh](https://linkedin.com/in/narges-shayesteh)
