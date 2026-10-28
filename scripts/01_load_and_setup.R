#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("SeuratData", quietly = TRUE)) {
  stop("Install SeuratData and its ifnb dataset as described in README.md.")
}
ifnb <- SeuratData::LoadData("ifnb")
ifnb <- UpdateSeuratObject(ifnb)
saveRDS(ifnb, file = "data/raw/ifnb_raw.rds")

message("Loaded the SeuratData ifnb object:")
print(ifnb)
