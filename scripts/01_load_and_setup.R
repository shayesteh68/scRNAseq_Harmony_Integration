#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/reference/ifnb_donor_metadata", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/de_celltype", recursive = TRUE, showWarnings = FALSE)
dir.create("run_logs", recursive = TRUE, showWarnings = FALSE)

raw_file <- "data/raw/ifnb_raw.rds"
if (file.exists(raw_file)) {
  message("Using the existing local IFNB object: ", raw_file)
  ifnb <- readRDS(raw_file)
} else {
  if (!requireNamespace("SeuratData", quietly = TRUE)) {
    stop("No local IFNB RDS found. Install SeuratData and its ifnb dataset as described in README.md.")
  }
  ifnb <- SeuratData::LoadData("ifnb")
}
ifnb <- UpdateSeuratObject(ifnb)

# SeuratData's IFNB object lacks donor IDs. Reuse the high-confidence Demuxlet
# assignments referenced by the official Seurat pseudobulk vignette.
donor_dir <- Sys.getenv("IFNB_DONOR_METADATA_DIR", "data/reference/ifnb_donor_metadata")
dir.create(donor_dir, recursive = TRUE, showWarnings = FALSE)
donor_sources <- c(
  "ye1.ctrl.8.10.sm.best" = "https://raw.githubusercontent.com/yelabucsf/demuxlet_paper_code/master/fig3/ye1.ctrl.8.10.sm.best",
  "ye2.stim.8.10.sm.best" = "https://raw.githubusercontent.com/yelabucsf/demuxlet_paper_code/master/fig3/ye2.stim.8.10.sm.best"
)
expected_md5 <- c(
  "ye1.ctrl.8.10.sm.best" = "eaef34891cdc6ce3896de1113b4bdd85",
  "ye2.stim.8.10.sm.best" = "812d2a2142210a566737e9cc36a9876a"
)
donor_paths <- file.path(donor_dir, names(donor_sources))
names(donor_paths) <- names(donor_sources)
for (name in names(donor_sources)) {
  path <- donor_paths[[name]]
  if (!file.exists(path)) {
    status <- tryCatch(
      download.file(donor_sources[[name]], path, mode = "wb", quiet = TRUE, method = "libcurl"),
      error = function(e) 1L
    )
    if (!identical(status, 0L) || !file.exists(path)) {
      stop("Failed to download donor metadata: ", donor_sources[[name]])
    }
  }
  observed_md5 <- unname(tools::md5sum(path))
  if (!identical(tolower(observed_md5), expected_md5[[name]])) {
    stop("Donor metadata checksum mismatch for ", name, "; preserve the file and inspect its source.")
  }
}

ctrl <- read.delim(donor_paths[["ye1.ctrl.8.10.sm.best"]], stringsAsFactors = FALSE, check.names = FALSE)
stim <- read.delim(donor_paths[["ye2.stim.8.10.sm.best"]], stringsAsFactors = FALSE, check.names = FALSE)
required_donor_columns <- c("BARCODE", "BEST")
if (!all(required_donor_columns %in% names(ctrl)) || !all(required_donor_columns %in% names(stim))) {
  stop("Demuxlet tables must contain BARCODE and BEST columns.")
}
ctrl_donors <- sort(unique(as.character(ctrl$BEST[grepl("^SNG-", ctrl$BEST)])))
stim_donors <- sort(unique(as.character(stim$BEST[grepl("^SNG-", stim$BEST)])))
if (!length(ctrl_donors) || !setequal(ctrl_donors, stim_donors)) {
  stop("CTRL and STIM Demuxlet files do not contain the same singlet donor IDs.")
}
ctrl$stim <- "CTRL"
stim$stim <- "STIM"
donor_info <- rbind(ctrl, stim)
donor_info$BARCODE <- gsub("-", ".", donor_info$BARCODE, fixed = TRUE)
donor_info <- donor_info[grepl("^SNG-", donor_info$BEST), , drop = FALSE]
donor_info <- donor_info[
  !duplicated(paste(donor_info$BARCODE, donor_info$stim)) &
    !duplicated(paste(donor_info$BARCODE, donor_info$stim), fromLast = TRUE), , drop = FALSE
]
donor_key <- paste(donor_info$BARCODE, donor_info$stim, sep = "__")
if (!nrow(donor_info) || anyDuplicated(donor_key)) {
  stop("No unique high-confidence Demuxlet barcode assignments were found.")
}
donor_map <- setNames(as.character(donor_info$BEST), donor_key)
cell_key <- paste(gsub("-", ".", colnames(ifnb), fixed = TRUE), ifnb$stim, sep = "__")
donor_id <- unname(donor_map[cell_key])
donor_id[is.na(donor_id)] <- "unknown"
ifnb <- AddMetaData(
  ifnb,
  metadata = data.frame(donor_id = donor_id, row.names = colnames(ifnb), stringsAsFactors = FALSE)
)
known <- ifnb$donor_id != "unknown"
paired <- tapply(ifnb$stim[known], ifnb$donor_id[known], function(x) all(c("CTRL", "STIM") %in% x))
message("High-confidence donor assignments: ", sum(known), "/", ncol(ifnb), " cells; donors represented: ", length(unique(ifnb$donor_id[known])))
message("Donors represented in both conditions: ", sum(paired))
if (!any(paired)) stop("No donor is represented in both CTRL and STIM after barcode matching.")
donor_counts <- as.data.frame(table(
  donor_id = factor(ifnb$donor_id[known], levels = sort(unique(ifnb$donor_id[known]))),
  stim = factor(ifnb$stim[known], levels = c("CTRL", "STIM"))
))
names(donor_counts)[names(donor_counts) == "Freq"] <- "n_cells"
write.csv(donor_counts, "results/de_celltype/donor_assignment_counts.csv", row.names = FALSE)
writeLines(c(
  "IFNB donor assignment summary",
  paste0("Cells: ", ncol(ifnb)),
  paste0("Unique high-confidence singlet assignments: ", sum(known)),
  paste0("Unassigned/ambiguous cells excluded from pseudobulk DE: ", sum(!known)),
  paste0("Demuxlet donors represented: ", length(unique(ifnb$donor_id[known]))),
  paste0("Donors represented in both CTRL and STIM: ", sum(paired)),
  "Barcode assignments were matched within stimulation condition.",
  paste0("Control source MD5: ", tolower(unname(tools::md5sum(donor_paths[["ye1.ctrl.8.10.sm.best"]])))),
  paste0("Stimulated source MD5: ", tolower(unname(tools::md5sum(donor_paths[["ye2.stim.8.10.sm.best"]]))))
), "results/de_celltype/donor_assignment_summary.txt")

saveRDS(ifnb, file = raw_file)
writeLines(capture.output(sessionInfo()), "run_logs/sessionInfo.txt")

message("Loaded the SeuratData ifnb object:")
print(ifnb)
