# Single-Cell RNA-Seq Integration with Harmony

**PBMC control vs. IFN-β stimulation - Harmony batch integration, unsupervised clustering and per-cell-type differential expression**

![Language](https://img.shields.io/badge/language-R-blue.svg)
![Integration](https://img.shields.io/badge/integration-Harmony-orange.svg)
![Data](https://img.shields.io/badge/data-scRNA--seq-green.svg)

---

## Overview

This repository contains a reproducible single-cell RNA-sequencing (scRNA-seq) workflow written in **R**. Starting from raw PBMC count data profiled in two conditions - untreated **control** and **IFN-β stimulation** - the pipeline performs quality control, normalization, batch/condition integration with **Harmony**, graph-based clustering, and per-cell-type **differential expression (DE)** analysis between the two conditions.

The central question addressed here is not only "which cell types are present", but **"how does each individual cell population respond to IFN-β stimulation"**.

---

## Dataset and experimental design

- Two conditions: `control` (`ctrl`) and `IFN-β stimulated` (`stim`).
- Cells are annotated into **13 immune cell populations** (table below).
- Raw count data is **not** tracked in this repository. Place the input files under `data/raw/`; intermediate and processed objects produced by the pipeline are written to `data/processed/`.

---

## Cell populations profiled

| Cell population | Cells (total) | Control | IFN-β stimulated |
| :--- | ---: | ---: | ---: |
| CD14 Mono | 4362 | 2215 | 2147 |
| CD4 Naive T | 2504 | 978 | 1526 |
| CD4 Memory T | 1762 | 859 | 903 |
| CD16 Mono | 1044 | 507 | 537 |
| B | 977 | 407 | 570 |
| CD8 T | 813 | 352 | 461 |
| T activated | 633 | 300 | 333 |
| NK | 619 | 298 | 321 |
| DC | 472 | 258 | 214 |
| B Activated | 388 | 185 | 203 |
| Mk | 236 | 115 | 121 |
| pDC | 132 | 51 | 81 |
| Eryth | 55 | 23 | 32 |
| **TOTAL** | **13997** | **6548** | **7449** |

---

## Analysis pipeline

| Step | Script | Purpose |
| :---: | :--- | :--- |
| 1 | `scripts/01_load_and_setup.R` | Load the count data, build the Seurat object and set up the project |
| 2 | `scripts/02_qc_and_filtering.R` | Quality control and cell filtering; QC violin plots (`figures/qc/`) |
| 3 | `scripts/03_integration_harmony.R` | Integrate conditions/batches with Harmony; UMAP before and after integration |
| 4 | `scripts/04_celltype_DE_IFNB.R` | Per-cell-type DE analysis (control vs IFN-β), ISG signature dot plot, volcano plot, DE summary table |
| 5 | `scripts/05_visualization.R` | Final publication-ready figures |

---

## Repository structure
```text
scRNAseq_Harmony_Integration/
+-- data/
|   +-- raw/                 # input count data (place here)
|   +-- processed/           # intermediate / processed objects
+-- figures/
|   +-- qc/
|   |   +-- pre_filtering_violin.png
|   |   +-- post_filtering_violin.png
|   +-- integration/
|   |   +-- unintegrated_umap.png
|   |   +-- harmony_integrated_umap.png
|   |   +-- integration_comparison_portfolio.png
|   +-- de_celltype/
|       +-- ISG_signature_dotplot.png
|       +-- Monocyte_IFNB_Volcano.png
+-- results/
|   +-- de_celltype/
|       +-- DE_summary_top_genes.csv
+-- scripts/
|   +-- 01_load_and_setup.R
|   +-- 02_qc_and_filtering.R
|   +-- 03_integration_harmony.R
|   +-- 04_celltype_DE_IFNB.R
|   +-- 05_visualization.R
+-- .gitignore
+-- environment.yml
+-- README.md

---

## Requirements and environment

- **R** with the **Seurat** and **Harmony** packages, plus the analysis dependencies used by the pipeline.
- The complete, pinned environment is defined in `environment.yml` (Conda).

---

## Installation and usage

### 1. Clone the repository

bash
git clone https://github.com/shayesteh68/scRNAseq_Harmony_Integration.git
cd scRNAseq_Harmony_Integration

### 2. Create and activate the environment

bash
conda env create -f environment.yml
conda activate scRNA_integration

### 3. Run the pipeline in order

bash
Rscript scripts/01_load_and_setup.R
Rscript scripts/02_qc_and_filtering.R
Rscript scripts/03_integration_harmony.R
Rscript scripts/04_celltype_DE_IFNB.R
Rscript scripts/05_visualization.R

---

## Results

- **`results/de_celltype/DE_summary_top_genes.csv`** - one row per cell population, with the number of cells per condition and the top 10 up- and down-regulated genes (`top_up`, `top_dn`).
- **`figures/integration/`** - side-by-side UMAP comparison before and after Harmony integration.
- **`figures/de_celltype/`** - ISG signature dot plot across cell populations and the monocyte volcano plot.

---

## Author

**Narges Shayesteh**
- LinkedIn: [linkedin.com/in/narges-shayesteh](https://linkedin.com/in/narges-shayesteh)
- GitHub: [@shayesteh68](https://github.com/shayesteh68)