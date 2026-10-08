library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

dir.create("figures/de_celltype", recursive = TRUE, showWarnings = FALSE)

# 1. Load Integrated Object
obj <- readRDS("data/processed/ifnb_harmony_integrated.rds")
DefaultAssay(obj) <- "RNA"
Idents(obj) <- "seurat_annotations"

# 2. Canonical Type I Interferon Stimulated Genes (ISGs)
isg_genes <- c("ISG15", "IFIT1", "IFIT2", "IFIT3", "IFI6", "MX1", "OAS1", "STAT1", "CXCL10", "IRF7")

# Filter genes present in the dataset
isg_genes <- isg_genes[isg_genes %in% rownames(obj)]

# Create a combined identity: CellType_Condition
obj$celltype_stim <- paste(obj$seurat_annotations, obj$stim, sep = "_")
Idents(obj) <- "celltype_stim"

# 3. DotPlot: ISG Expression Across Cell Types and Conditions
p_dot <- DotPlot(obj, features = isg_genes, cols = c("lightgrey", "firebrick3")) + 
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  ) +
  labs(title = "Type I Interferon (IFN-β) Response Signature Across Cell Types")

ggsave("figures/de_celltype/ISG_signature_dotplot.png", plot = p_dot, width = 11, height = 9, dpi = 300)

# 4. Volcano Plot for CD14 Monocytes
mono_de_file <- "results/de_celltype/CD14_Mono_STIM_vs_CTRL.csv"
if (!file.exists(mono_de_file)) {
  # In case filename differs, search for Mono file
  mono_files <- list.files("results/de_celltype", pattern = "Mono.*STIM_vs_CTRL\\.csv", full.names = TRUE)
  if (length(mono_files) > 0) mono_de_file <- mono_files[1]
}

if (file.exists(mono_de_file)) {
  de_data <- read.csv(mono_de_file)
  de_data$neg_log10_padj <- -log10(de_data$p_val_adj + 1e-300)
  
  de_data$significance <- "Not Significant"
  de_data$significance[de_data$p_val_adj < 0.05 & de_data$avg_log2FC > 1] <- "Upregulated (STIM)"
  de_data$significance[de_data$p_val_adj < 0.05 & de_data$avg_log2FC < -1] <- "Downregulated (STIM)"
  
  top_labels <- de_data %>% 
    filter(significance != "Not Significant") %>% 
    arrange(p_val_adj, desc(abs(avg_log2FC))) %>% 
    head(10)

  p_volcano <- ggplot(de_data, aes(x = avg_log2FC, y = neg_log10_padj, color = significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Upregulated (STIM)" = "red2", "Downregulated (STIM)" = "dodgerblue2", "Not Significant" = "grey70")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_text(data = top_labels, aes(label = gene), vjust = -0.5, size = 3.5, color = "black", check_overlap = TRUE) +
    theme_bw(base_size = 12) +
    labs(
      title = "Volcano Plot: Differential Expression in Monocytes",
      subtitle = "IFN-β Stimulated vs Control",
      x = "Average log2 Fold Change",
      y = "-log10 Adjusted P-value"
    )

  ggsave("figures/de_celltype/Monocyte_IFNB_Volcano.png", plot = p_volcano, width = 8, height = 6, dpi = 300)
}

message("Visualization scripts finished successfully!")
