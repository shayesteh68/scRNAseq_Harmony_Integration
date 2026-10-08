library(Seurat)
library(SeuratData)
library(patchwork)

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)

# Load dataset from the installed package
data("ifnb")

# Update object structure to Seurat v5 format
message("Updating Seurat object structure to v5...")
ifnb <- UpdateSeuratObject(ifnb)

# Save standardized object to raw directory
saveRDS(ifnb, file = "data/raw/ifnb_raw.rds")

message("--- Dataset Summary ---")
print(ifnb)
