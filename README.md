# Multi-Condition scRNA-Seq Integration & Batch Effect Correction with Harmony

[![Pipeline: Seurat v5](https://img.shields.io/badge/Pipeline-Seurat%20v5-blue.svg)](https://satijalab.org/seurat/)
[![Integration: Harmony](https://img.shields.io/badge/Integration-Harmony-orange.svg)](https://github.com/immunogenomics/harmony)
[![R: 4.3+](https://img.shields.io/badge/R-4.3+-276DC3.svg)](https://www.r-project.org/)
[![Status: Production Ready](https://img.shields.io/badge/Status-Production%20Ready-success.svg)](https://github.com/shayesteh68)

## Overview
This repository provides an end-to-end, publication-grade single-cell RNA-sequencing (scRNA-seq) workflow demonstrating **batch effect correction and cross-condition integration** using **Harmony** in conjunction with **Seurat v5**.

The pipeline integrates PBMC profiles across distinct biological conditions (**Control vs. Interferon-beta [IFN-β] stimulated**), aligns cell states across technical and biological covariates, and identifies cell-type-specific transcriptional responses triggered by interferon activation.

---

## Key Workflow Modules

1. **`01_load_and_setup.R`**: Environment initialization, SeuratData loading (`ifnb` stimulation cohort), and multimodal metadata structuring.
2. **`02_qc_and_filtering.R`**: Strict cell quality control (filtering by mitochondrial read percentage, unique feature counts, and total UMIs).
3. **`03_integration_harmony.R`**: Normalization, feature selection (HVGs), PCA dimension reduction, and **Harmony** matrix factorization to eliminate batch-induced variance while preserving biological heterogeneity.
4. **`04_celltype_DE_IFNB.R`**: Cell type annotation, unsupervised graph-based clustering, and stratified differential expression (DE) identifying IFN-stimulated marker signatures across lineages (CD4 T cells, CD8 T cells, B cells, Monocytes, NK cells, and DCs).
5. **`05_visualization.R`**: Publication-ready UMAP plots (pre- vs. post-integration), conserved marker heatmaps, and response dot plots.

---

## Directory Structure

```text
scRNAseq_Harmony_Integration/
├── data/                 # Raw and processed RDS objects / SeuratData cache
├── figures/              # Publication-ready visualization outputs
│   ├── qc/               # Pre/post filtering violin and scatter plots
│   ├── integration/      # Pre vs. Post-Harmony integration UMAP embeddings
│   └── de_celltype/      # Cell-type specific response heatmaps & dot plots
├── results/              # Differential expression result tables (.tsv / .csv)
├── scripts/
│   ├── 01_load_and_setup.R
│   ├── 02_qc_and_filtering.R
│   ├── 03_integration_harmony.R
│   ├── 04_celltype_DE_IFNB.R
│   └── 05_visualization.R
├── environment.yml       # Reproducible Conda/Mamba bioinformatics environment
└── README.md
```

---

## Results Highlights

- **Effective Batch Removal**: Complete cross-dataset alignment without over-clustering or masking biological divergence.
- **Robust Lineage Preservation**: Major immune subsets retain discrete identity markers post-correction.
- **Interferon Signature Mapping**: High-confidence detection of canonical IFN-stimulated genes (ISG15, IFI6, IFIT1, MX1) stratified across immune cell subsets.

---

## Reproducibility & Execution

```bash
# 1. Clone repository
git clone https://github.com/shayesteh68/scRNAseq_Harmony_Integration.git
cd scRNAseq_Harmony_Integration

# 2. Run execution pipeline
Rscript scripts/01_load_and_setup.R
Rscript scripts/02_qc_and_filtering.R
Rscript scripts/03_integration_harmony.R
Rscript scripts/04_celltype_DE_IFNB.R
Rscript scripts/05_visualization.R
```

---

## Author
**Narges Shayesteh**  
*Bioinformatics Specialist | Transcriptomics & Single-Cell Genomics*  
- LinkedIn: [narges-shayesteh](https://www.linkedin.com/in/narges-shayesteh)  
- GitHub: [@shayesteh68](https://github.com/shayesteh68)
