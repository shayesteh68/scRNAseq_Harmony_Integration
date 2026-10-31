# PBMC IFN-beta integration with Harmony

This repository analyzes the public SeuratData IFNB PBMC example (13,999 cells and 14,053 features). It runs cell-level QC, exploratory Harmony alignment, clustering, donor-level pseudobulk differential expression, and plots. Differential expression uses raw RNA counts and paired donor models; it does not use Harmony-corrected values.

The IFNB example derives from the public PBMC stimulation study by Kang et al., “Multiplexed droplet single-cell RNA-sequencing using natural genetic variation,” *Nature Biotechnology* (2018), [doi:10.1038/nbt.4042](https://doi.org/10.1038/nbt.4042); GEO series [GSE96583](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583). This is an independent computational case study, not a commissioned analysis or a collaboration with the authors or their institutions.

## Environment and data

For the exact audited Linux x86_64 environment, create it from the lockfile:

```bash
conda create --name scRNA_integration --file environment-linux-64.lock.txt
conda activate scRNA_integration
```

`environment.yml` pins the direct R dependencies and can be used on other platforms with `conda env create -f environment.yml`.

Install the IFNB object once if it is not already available locally:

```r
install.packages("remotes")
remotes::install_github("satijalab/seurat-data")
SeuratData::InstallData("ifnb")
```

`scripts/01_load_and_setup.R` uses `data/raw/ifnb_raw.rds` when present; otherwise it loads the object from SeuratData. The IFNB object has stimulation and cell-type annotations but no donor IDs. The script therefore downloads the two high-confidence Demuxlet singlet tables referenced by the [Seurat differential-expression vignette](https://satijalab.org/seurat/articles/de_vignette), checks their MD5 digests, and matches barcodes within each condition. The inputs are `ye1.ctrl.8.10.sm.best` (MD5 `eaef34891cdc6ce3896de1113b4bdd85`) and `ye2.stim.8.10.sm.best` (MD5 `812d2a2142210a566737e9cc36a9876a`) from the [Demuxlet paper-code repository](https://github.com/yelabucsf/demuxlet_paper_code/tree/master/fig3). Donor matching relies on these external source files; 203 cells without a unique high-confidence assignment are excluded from pseudobulk DE, but remain in QC and visualization.

## Run

Run the scripts in order from the repository root:

```bash
Rscript scripts/01_load_and_setup.R
Rscript scripts/02_qc_and_filtering.R
Rscript scripts/03_integration_harmony.R
Rscript scripts/04_celltype_DE_IFNB.R
Rscript scripts/05_visualization.R
```

QC retains cells with 200-2,499 detected features. The dataset has no features named with the human `MT-` convention, so mitochondrial percentage cannot be computed and is not shown or used as a filter. This filter is not a dedicated doublet detector. In the audited run, 13,997 of 13,999 cells passed the feature-count filter.

The Harmony script deliberately aligns on `stim` only as an exploratory demonstration. Stimulation is the biological treatment, not a nuisance batch; Harmony can remove treatment-associated signal. Do not use the Harmony embedding to estimate IFN-beta effects. The DE script uses raw RNA counts, sums counts by donor and condition within each annotated cell type, requires at least two complete donors, and fits `~ donor_id + stim`. Before QC, 13,796 of 13,999 cells matched to eight donors; 203 lacked a unique high-confidence assignment. After QC, 13,997 cells remained, including 13,794 assigned cells. The paired donor-by-condition groups used in DE contain 13,790 cells across 13 cell types (6,547 CTRL and 7,243 STIM); four assigned cells were omitted because their donor lacked the other condition in that cell type. Eight paired donors were available in 12 cell types and six in Eryth. Benjamini-Hochberg adjustment is reported both within each cell type and across all tested cell-type/gene pairs (`padj_global_BH`): the run produced 7,231 discoveries using within-cell-type adjustment and 7,122 using the global adjustment. The monocyte volcano plot uses the global adjustment.

Outputs include:

- `results/de_celltype/IFNB_pseudobulk_DESeq2_by_celltype.csv`: gene-level estimates, within-cell-type adjusted p-values, global adjusted p-values, and donor counts.
- `results/de_celltype/IFNB_pseudobulk_summary.csv`: cell counts and significant-gene counts per cell type under both adjustment scopes.
- `results/de_celltype/DE_summary_top_genes.csv`: top significant genes by direction and cell type.
- `results/de_celltype/donor_assignment_counts.csv` and `donor_assignment_summary.txt`: donor matching coverage and source checksums.
- `figures/qc/`, `figures/integration/`, and `figures/de_celltype/`: QC, exploratory embedding, ISG, and monocyte DE plots.
- `run_logs/`: per-step console logs and the R session information from the audited execution.

## Interpretation and limits

This is a small public PBMC example, not a clinical validation cohort. Donor assignments depend on externally supplied Demuxlet tables and exact barcode matching; unmatched cells are excluded from donor-level tests. The cell-type labels are supplied with the source dataset and are not independently validated here. Per-cell-type and global adjusted p-values answer different multiple-testing questions; consult both before interpreting a gene. Harmony is intentionally run against treatment for illustration only, with the risk that it removes the biological signal of interest.

The integration and DE approach follows the [Seurat IFNB integration vignette](https://satijalab.org/seurat/articles/integration_introduction) and [differential-expression vignette](https://satijalab.org/seurat/articles/de_vignette). Code is licensed under MIT (see `LICENSE`); the dataset and third-party software retain their respective source terms.
