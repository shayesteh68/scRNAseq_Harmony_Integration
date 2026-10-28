#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(DESeq2)
})

input_rds <- "data/processed/ifnb_harmony_integrated.rds"
if (!file.exists(input_rds)) stop("Integrated object missing: ", input_rds)
dir.create("results/de_celltype", recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input_rds)
required <- c("stim", "donor_id", "seurat_annotations")
metadata <- obj[[]]
if (!all(required %in% colnames(metadata))) {
  stop("Object metadata must contain: ", paste(required, collapse = ", "))
}
if (anyNA(metadata$stim) || anyNA(metadata$donor_id) || anyNA(metadata$seurat_annotations)) {
  stop("Stimulation, donor, and cell-type metadata must not contain missing values.")
}
metadata$stim <- as.character(metadata$stim)
metadata$donor_id <- as.character(metadata$donor_id)
metadata$seurat_annotations <- as.character(metadata$seurat_annotations)
if (!all(metadata$stim %in% c("CTRL", "STIM"))) {
  stop("Expected stim values CTRL and STIM.")
}

DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj, assay = "RNA")
metadata <- obj[[]][colnames(obj), , drop = FALSE]
metadata$stim <- as.character(metadata$stim)
metadata$donor_id <- as.character(metadata$donor_id)
metadata$seurat_annotations <- as.character(metadata$seurat_annotations)
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
cell_types <- sort(unique(metadata$seurat_annotations))
result_rows <- list()
summary_rows <- list()

for (cell_type in cell_types) {
  cells <- rownames(metadata)[metadata$seurat_annotations == cell_type]
  cell_meta <- metadata[cells, , drop = FALSE]
  by_donor <- split(cell_meta$stim, cell_meta$donor_id)
  paired_donors <- names(by_donor)[
    vapply(by_donor, function(x) all(c("CTRL", "STIM") %in% x), logical(1))
  ]
  if (length(paired_donors) < 2L) {
    message("Skipping ", cell_type, ": fewer than two donors have both conditions.")
    next
  }

  cell_meta <- cell_meta[cell_meta$donor_id %in% paired_donors, , drop = FALSE]
  sample_key <- paste(cell_meta$donor_id, cell_meta$stim, sep = "__")
  sample_levels <- unique(sample_key)
  sample_meta <- data.frame(
    donor_id = cell_meta$donor_id[match(sample_levels, sample_key)],
    stim = cell_meta$stim[match(sample_levels, sample_key)],
    row.names = sample_levels,
    stringsAsFactors = FALSE
  )
  sample_meta$donor_id <- factor(sample_meta$donor_id)
  sample_meta$stim <- relevel(factor(sample_meta$stim, levels = c("CTRL", "STIM")), ref = "CTRL")

  design_matrix <- model.matrix(~ donor_id + stim, data = sample_meta)
  if (qr(design_matrix)$rank < ncol(design_matrix)) {
    message("Skipping ", cell_type, ": paired donor design is not full rank.")
    next
  }

  aggregation <- Matrix::sparse.model.matrix(
    ~ 0 + factor(sample_key, levels = sample_levels)
  )
  colnames(aggregation) <- sample_levels
  pseudobulk <- as.matrix(counts[, rownames(cell_meta), drop = FALSE] %*% aggregation)
  if (anyNA(pseudobulk) || any(!is.finite(pseudobulk)) || any(pseudobulk < 0) ||
      any(pseudobulk > .Machine$integer.max) || any(pseudobulk != floor(pseudobulk))) {
    stop("Pseudobulk counts must be finite non-negative integers within R's integer range.")
  }
  storage.mode(pseudobulk) <- "integer"
  if (!identical(colnames(pseudobulk), rownames(sample_meta))) {
    stop("Pseudobulk columns and sample metadata are out of order for ", cell_type, ".")
  }

  dds <- DESeqDataSetFromMatrix(
    countData = pseudobulk,
    colData = sample_meta,
    design = ~ donor_id + stim
  )
  keep_genes <- rowSums(counts(dds) >= 10) >= 2L
  dds <- dds[keep_genes, ]
  if (nrow(dds) < 2L) {
    message("Skipping ", cell_type, ": fewer than two genes pass the count filter.")
    next
  }
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("stim", "STIM", "CTRL"), alpha = 0.05)
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df$cell_type <- cell_type
  res_df$n_paired_donors <- length(paired_donors)
  res_df$significant_padj_lt_0_05 <- !is.na(res_df$padj) & res_df$padj < 0.05
  result_rows[[cell_type]] <- res_df

  summary_rows[[cell_type]] <- data.frame(
    cell_type = cell_type,
    n_paired_donors = length(paired_donors),
    n_cells_CTRL = sum(cell_meta$stim == "CTRL"),
    n_cells_STIM = sum(cell_meta$stim == "STIM"),
    n_genes_tested = nrow(dds),
    n_significant_padj_lt_0_05 = sum(res_df$significant_padj_lt_0_05),
    design = "~ donor_id + stim",
    contrast = "STIM versus CTRL",
    stringsAsFactors = FALSE
  )
  message(cell_type, ": paired pseudobulk DE completed for ", length(paired_donors), " donors.")
}

if (!length(result_rows)) stop("No cell type had enough paired donor samples for pseudobulk DE.")
all_results <- do.call(rbind, result_rows)
all_results <- all_results[order(all_results$cell_type, all_results$padj, na.last = TRUE), ]
write.csv(all_results, "results/de_celltype/IFNB_pseudobulk_DESeq2_by_celltype.csv", row.names = FALSE)
write.csv(do.call(rbind, summary_rows), "results/de_celltype/IFNB_pseudobulk_summary.csv", row.names = FALSE)
message("Saved donor-level pseudobulk DESeq2 results and summary.")
