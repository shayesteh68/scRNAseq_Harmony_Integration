library(Seurat)
library(ggplot2)
library(patchwork)

dir.create("figures/qc", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# Load raw updated object
ifnb <- readRDS("data/raw/ifnb_raw.rds")

# Calculate Mitochondrial Percentage
ifnb[["percent.mt"]] <- PercentageFeatureSet(ifnb, pattern = "^MT-")

# Plot pre-filtering QC metrics
p_pre_vln <- VlnPlot(ifnb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "stim")
ggsave("figures/qc/pre_filtering_violin.png", plot = p_pre_vln, width = 10, height = 5)

# Filter low-quality cells and doublets/debris
ifnb_filtered <- subset(ifnb, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Plot post-filtering QC metrics
p_post_vln <- VlnPlot(ifnb_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "stim")
ggsave("figures/qc/post_filtering_violin.png", plot = p_post_vln, width = 10, height = 5)

# Report cell counts
message("--- Cell Count Summary ---")
message("Original cell count: ", ncol(ifnb))
message("Filtered cell count: ", ncol(ifnb_filtered))

# Save processed dataset
saveRDS(ifnb_filtered, file = "data/processed/ifnb_filtered.rds")
message("Filtered dataset successfully saved to data/processed/ifnb_filtered.rds")
