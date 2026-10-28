library(Seurat)
library(ggplot2)

dir.create("figures/de_celltype", recursive = TRUE, showWarnings = FALSE)
obj <- readRDS("data/processed/ifnb_harmony_integrated.rds")
DefaultAssay(obj) <- "RNA"

isg_genes <- c("ISG15", "IFIT1", "IFIT2", "IFIT3", "IFI6", "MX1", "OAS1", "STAT1", "CXCL10", "IRF7")
isg_genes <- isg_genes[isg_genes %in% rownames(obj)]
if (!length(isg_genes)) stop("None of the configured interferon-response genes are present.")

obj$celltype_stim <- paste(obj$seurat_annotations, obj$stim, sep = "_")
Idents(obj) <- "celltype_stim"
p_dot <- DotPlot(obj, features = isg_genes, cols = c("lightgrey", "firebrick3")) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.title = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)) +
  labs(title = "IFN-beta response genes by cell type and condition")
ggsave("figures/de_celltype/ISG_signature_dotplot.png", p_dot, width = 11, height = 9, dpi = 300)

de_file <- "results/de_celltype/IFNB_pseudobulk_DESeq2_by_celltype.csv"
if (file.exists(de_file)) {
  de <- read.csv(de_file, stringsAsFactors = FALSE)
  monocytes <- de[de$cell_type == "CD14 Mono" & !is.na(de$padj), , drop = FALSE]
  if (nrow(monocytes)) {
    monocytes$neg_log10_padj <- -log10(pmax(monocytes$padj, .Machine$double.xmin))
    monocytes$significant <- monocytes$padj < 0.05
    volcano <- ggplot(monocytes, aes(log2FoldChange, neg_log10_padj, color = significant)) +
      geom_point(alpha = 0.7, size = 1.5) +
      scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#d95f02")) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      labs(
        title = "CD14 Mono: paired donor pseudobulk DE",
        x = "DESeq2 log2 fold change (STIM versus CTRL)",
        y = "-log10 adjusted p-value",
        color = "padj < 0.05"
      ) +
      theme_bw()
    ggsave("figures/de_celltype/Monocyte_IFNB_Volcano.png", volcano, width = 8, height = 6, dpi = 300)
  } else {
    message("No CD14 Mono pseudobulk results available for a volcano plot.")
  }
} else {
  message("Run scripts/04_celltype_DE_IFNB.R first to create the pseudobulk results.")
}
