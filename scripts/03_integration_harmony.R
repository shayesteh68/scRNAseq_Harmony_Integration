library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)

set.seed(42)
dir.create("figures/integration", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

set.seed(42)
ifnb <- readRDS("data/processed/ifnb_filtered.rds")
ifnb <- NormalizeData(ifnb, verbose = FALSE)
ifnb <- FindVariableFeatures(ifnb, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
ifnb <- ScaleData(ifnb, verbose = FALSE)
ifnb <- RunPCA(ifnb, npcs = 30, verbose = FALSE)

message("Generating unintegrated UMAP...")
ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "pca", reduction.name = "umap_unintegrated", seed.use = 42, verbose = FALSE)
p1 <- DimPlot(ifnb, reduction = "umap_unintegrated", group.by = "stim") +
  ggtitle("Unintegrated UMAP: control and stimulated")
p2 <- DimPlot(ifnb, reduction = "umap_unintegrated", group.by = "seurat_annotations") +
  ggtitle("Unintegrated UMAP: cell labels")

message("Running Harmony against stimulation status for exploratory cell alignment only.")
message("Alignment may remove IFN-associated biology; never use this embedding for treatment-effect inference.")
ifnb <- RunHarmony(ifnb, group.by.vars = "stim", reduction.use = "pca", dims.use = 1:30, verbose = FALSE)
ifnb <- RunUMAP(ifnb, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony", seed.use = 42, verbose = FALSE)
ifnb <- FindNeighbors(ifnb, reduction = "harmony", dims = 1:30, verbose = FALSE)
ifnb <- FindClusters(ifnb, resolution = 0.5, random.seed = 42, verbose = FALSE)

p3 <- DimPlot(ifnb, reduction = "umap_harmony", group.by = "stim") +
  ggtitle("Harmony space: exploratory alignment")
p4 <- DimPlot(ifnb, reduction = "umap_harmony", group.by = "seurat_annotations") +
  ggtitle("Harmony UMAP: cell labels")
ggsave("figures/integration/unintegrated_umap.png", p1 | p2, width = 16, height = 7, bg = "white")
ggsave("figures/integration/harmony_integrated_umap.png", p3 | p4, width = 16, height = 7, bg = "white")
ggsave("figures/integration/integration_comparison_portfolio.png", (p1 | p3) / (p2 | p4), width = 16, height = 12, bg = "white")

saveRDS(ifnb, "data/processed/ifnb_harmony_integrated.rds")
message("Exploratory alignment saved.")
