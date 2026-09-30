library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)

dir.create("figures/integration", recursive = TRUE, showWarnings = FALSE)

# 1. Load Filtered Data
ifnb <- readRDS("data/processed/ifnb_filtered.rds")

# 2. Standard Preprocessing Workflow
ifnb <- NormalizeData(ifnb, verbose = FALSE)
ifnb <- FindVariableFeatures(ifnb, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
ifnb <- ScaleData(ifnb, verbose = FALSE)
ifnb <- RunPCA(ifnb, npcs = 30, verbose = FALSE)

# 3. Dimensionality Reduction Without Integration (Showing Batch Effect)
message("Generating Unintegrated UMAP...")
ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "pca", reduction.name = "umap_unintegrated", verbose = FALSE)

p1 <- DimPlot(ifnb, reduction = "umap_unintegrated", group.by = "stim") + 
  ggtitle("Unintegrated (Separated by Condition)")
p2 <- DimPlot(ifnb, reduction = "umap_unintegrated", group.by = "seurat_annotations", label = TRUE, repel = TRUE) + 
  ggtitle("Unintegrated (Cell Types)") + NoLegend()

ggsave("figures/integration/unintegrated_umap.png", plot = (p1 | p2), width = 14, height = 6)

# 4. Run Harmony Integration
message("Running Harmony Integration across stim condition...")
ifnb <- RunHarmony(ifnb, group.by.vars = "stim", reduction.use = "pca", dims.use = 1:30, verbose = FALSE)

# 5. UMAP & Clustering on Harmony Embeddings
message("Generating Integrated UMAP...")
ifnb <- RunUMAP(ifnb, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony", verbose = FALSE)
ifnb <- FindNeighbors(ifnb, reduction = "harmony", dims = 1:30, verbose = FALSE)
ifnb <- FindClusters(ifnb, resolution = 0.5, verbose = FALSE)

p3 <- DimPlot(ifnb, reduction = "umap_harmony", group.by = "stim") + 
  ggtitle("Harmony Integrated (Aligned by Condition)")
p4 <- DimPlot(ifnb, reduction = "umap_harmony", group.by = "seurat_annotations", label = TRUE, repel = TRUE) + 
  ggtitle("Harmony Integrated (Conserved Cell Types)") + NoLegend()

ggsave("figures/integration/harmony_integrated_umap.png", plot = (p3 | p4), width = 14, height = 6)

# Side-by-side comparison for portfolio figure
p_comparison <- (p1 | p3) / (p2 | p4)
ggsave("figures/integration/integration_comparison_portfolio.png", plot = p_comparison, width = 14, height = 12)

# Save the fully integrated object
saveRDS(ifnb, file = "data/processed/ifnb_harmony_integrated.rds")
message("Integrated dataset successfully saved to data/processed/ifnb_harmony_integrated.rds")
