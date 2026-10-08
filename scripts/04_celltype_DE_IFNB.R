library(Seurat)
library(dplyr)
library(ggplot2)

dir.create("results/de_celltype", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/de_celltype", recursive = TRUE, showWarnings = FALSE)

# Load integrated object (Harmony already computed)
obj <- readRDS("data/processed/ifnb_harmony_integrated.rds")

# Safety checks
stopifnot("stim" %in% colnames(obj@meta.data))
stopifnot("seurat_annotations" %in% colnames(obj@meta.data))

# Make sure identities are cell types
Idents(obj) <- "seurat_annotations"

celltypes <- levels(Idents(obj))

all_summaries <- list()

# Run DE within each cell type: STIM vs CTRL
# Using RNA assay's normalized data; for pseudobulk later we can add a stronger analysis.
DefaultAssay(obj) <- "RNA"

for (ct in celltypes) {
  message("Processing cell type: ", ct)

  cells_ct <- WhichCells(obj, idents = ct)
  obj_ct <- subset(obj, cells = cells_ct)

  # Ensure both conditions exist in this cell type
  tab <- table(obj_ct$stim)
  if (!all(c("CTRL", "STIM") %in% names(tab))) {
    message("Skipping ", ct, " (missing CTRL or STIM). Counts: ", paste(names(tab), tab, collapse = ", "))
    next
  }

  # Set identities to condition within this cell type
  Idents(obj_ct) <- "stim"

  # DE: STIM vs CTRL
  markers <- FindMarkers(
    obj_ct,
    ident.1 = "STIM",
    ident.2 = "CTRL",
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.10
  )

  markers <- markers %>%
    tibble::rownames_to_column("gene") %>%
    arrange(p_val_adj, desc(avg_log2FC))

  out_csv <- file.path("results/de_celltype", paste0(gsub("[ /]", "_", ct), "_STIM_vs_CTRL.csv"))
  write.csv(markers, out_csv, row.names = FALSE)

  # Summary: top up/down genes
  top_up <- markers %>% filter(!is.na(p_val_adj)) %>% arrange(p_val_adj, desc(avg_log2FC)) %>% head(10)
  top_dn <- markers %>% filter(!is.na(p_val_adj)) %>% arrange(p_val_adj, avg_log2FC) %>% head(10)

  all_summaries[[ct]] <- list(
    celltype = ct,
    n_cells = ncol(obj_ct),
    n_ctrl = as.integer(tab["CTRL"]),
    n_stim = as.integer(tab["STIM"]),
    top_up = paste(top_up$gene, collapse = ";"),
    top_dn = paste(top_dn$gene, collapse = ";")
  )
}

summary_df <- dplyr::bind_rows(lapply(all_summaries, as.data.frame))
write.csv(summary_df, "results/de_celltype/DE_summary_top_genes.csv", row.names = FALSE)

message("Saved per-celltype DE tables to results/de_celltype/")
message("Saved summary to results/de_celltype/DE_summary_top_genes.csv")
