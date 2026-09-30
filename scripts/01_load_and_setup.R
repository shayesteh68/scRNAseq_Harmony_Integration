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

 your object`) یک رفتار کاملاً استاندارد در Seurat نسخه ۵ است؛ ساختار دیتاست `ifnb` متعلق به ساختار Seurat v3 است و باید با تابع رسمی `UpdateSeuratObject()` به ساختار ساختاری Seurat v5 ارتقا یابد.

---

### Step 1: Update Script `01_load_and_setup.R`

دستور زیر را در ترمینال بزن تا اسکریپت آپدیت شود:
```bash
cat <<'EOF' > scripts/01_load_and_setup.R
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

message("--- Experimental Conditions (stim) ---")
print(table(ifnb$stim))

message("--- Cell Type Annotations ---")
print(table(ifnb$seurat_annotations))

message("Raw dataset successfully updated and saved to data/raw/ifnb_raw.rds")
