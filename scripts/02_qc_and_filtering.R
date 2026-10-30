library(Seurat)
library(ggplot2)
library(patchwork)

dir.create("figures/qc", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# Load raw updated object
ifnb <- readRDS("data/raw/ifnb_raw.rds")

mt_features <- grep("^MT-", rownames(ifnb), value = TRUE)
qc_features <- c("nFeature_RNA", "nCount_RNA")
if (length(mt_features)) {
  ifnb[["percent.mt"]] <- PercentageFeatureSet(ifnb, features = mt_features)
  qc_features <- c(qc_features, "percent.mt")
} else {
  message("No ^MT- features are present; mitochondrial percentage cannot be calculated and is not used for filtering.")
}

# Plot pre-filtering QC metrics
p_pre_vln <- VlnPlot(ifnb, features = qc_features, ncol = length(qc_features), group.by = "stim")
ggsave("figures/qc/pre_filtering_violin.png", plot = p_pre_vln, width = 10, height = 5, bg = "white")

# Filter by gene count and, when available, mitochondrial fraction. This is not a dedicated doublet caller.
if (length(mt_features)) {
  ifnb_filtered <- subset(ifnb, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
} else {
  ifnb_filtered <- subset(ifnb, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
}

# Plot post-filtering QC metrics
p_post_vln <- VlnPlot(ifnb_filtered, features = qc_features, ncol = length(qc_features), group.by = "stim")
ggsave("figures/qc/post_filtering_violin.png", plot = p_post_vln, width = 10, height = 5, bg = "white")

# Report cell counts
message("--- Cell Count Summary ---")
message("Original cell count: ", ncol(ifnb))
message("Filtered cell count: ", ncol(ifnb_filtered))

# Save processed dataset
saveRDS(ifnb_filtered, file = "data/processed/ifnb_filtered.rds")
message("Filtered dataset successfully saved to data/processed/ifnb_filtered.rds")
