###############################################################################
# Toil RNA-seq expected_count (UCSC Xena) → raw counts for CERVIX samples
# - Downloads the TcgaTargetGtex_gene_expected_count dataset (log2(x+1))
# - Keeps gene IDs as rownames and columns as sample IDs
# - Subsets to cervix samples based on current phenotype columns
# - Back-transforms to raw expected counts: raw = 2^log2 - 1
###############################################################################

## 0) Robust I/O for large files (vroom connection buffer ~128 MB)
Sys.setenv("VROOM_CONNECTION_SIZE" = as.character(131072 * 1024))

suppressPackageStartupMessages({
  library(UCSCXenaTools)
  library(dplyr)
  library(vroom)
})

## 1) Config
host            <- "https://toil.xenahubs.net"
dataset_chosen  <- "TcgaTargetGtex_gene_expected_count"
out_dir         <- "~/CESC_Network/1_Get_Data/"
dir.create(path.expand(out_dir), showWarnings = FALSE, recursive = TRUE)

#load(file = file.path(out_dir, "1_2_3_Image.RData"))

## 2) Build Xena hub and filter dataset
xh  <- XenaGenerate()
xh2 <- XenaFilter(xh, filterDatasets = dataset_chosen)
message("Filtered datasets in hub:"); print(datasets(xh2))
stopifnot(length(datasets(xh2)) == 1)

## 3) Create query (no per-sample subsetting here; we subset after reading)
xq <- XenaQuery(xh2)

## 4) Download with resume (curl). If file exists partially, resume.
dl <- XenaDownload(
  xq,
  destdir = out_dir,
  method  = "curl",    # use system curl for robust/resume
  extra   = "-C - -L", # resume and follow redirects
  force   = TRUE,
  max_try = 5
)

## 5) Prepare expression (matrix is log2(expected_count + 1))
expr_log_all <- XenaPrepare(dl)
# Expected shape: first column = gene id (named 'sample'), then many sample columns
stopifnot("sample" %in% colnames(expr_log_all))

## 6) Load phenotype and select current CERVIX columns
pheno <- vroom::vroom(
  paste0(host, "/download/TcgaTargetGTEX_phenotype.txt.gz"),
  show_col_types = FALSE
)

cervix_samples <- pheno %>%
  filter(grepl("CERVIX", detailed_category, ignore.case = TRUE) |
           grepl("CERVIX", `_primary_site`, ignore.case = TRUE) |
           grepl("CERVIX", `primary disease or tissue`, ignore.case = TRUE)) %>%
  pull(sample) %>%
  unique()
message("Cervix samples found: ", length(cervix_samples))

## 7) Preserve gene IDs as rownames and subset columns to cervix
# 7a) Extract gene IDs vector BEFORE converting to matrix
gene_ids <- expr_log_all$sample
stopifnot(length(gene_ids) == nrow(expr_log_all))

# 7b) Keep only numeric sample columns
expr_num <- dplyr::select(expr_log_all, -sample)

# 7c) Convert to a numeric matrix and assign rownames = gene IDs
expr_mat <- as.matrix(expr_num)
rownames(expr_mat) <- gene_ids
all(rownames(expr_mat) == gene_ids)  # TRUE

# 7d) Subset columns to the cervix sample intersection
keep_cols <- intersect(colnames(expr_mat), cervix_samples)
stopifnot(length(keep_cols) > 0)
expr_mat_sub <- expr_mat[, keep_cols, drop = FALSE]

## 8) Back-transform to raw expected counts from log2(x+1)
#    IMPORTANT: correct inverse is (2 ^ value) - 1
#    doi: 10.1038/s41525-020-00167-4. PMID: 33495453; PMCID: PMC7835362.
expr_raw <- 2^expr_mat_sub - 1
storage.mode(expr_raw) <- "double"
expr_raw <- data.frame(
  gene = rownames(expr_raw),
  expr_raw)

# --- Normalize cols of expr_raw to "TCGA-XX-XXXX-NN" ---
cn <- colnames(expr_raw)
cn[-1] <- gsub("\\.", "-", cn[-1])
cn[-1] <- sub("^(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}).*$", "\\1", cn[-1], perl = TRUE)
colnames(expr_raw) <- cn

dim(expr_raw)
#[1] 60498   319

## 9) Filter only metadata samples
Metadata <- vroom(file.path(out_dir,"1_1_1_Metadata.tsv"))
# --- Normalize specimenID of Metadata ---
Metadata$specimenID_portion <- sub("^(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}).*$", "\\1",
                                   Metadata$specimenID, perl = TRUE)
#Add GTEX samples to Metadata from colnames (expr_raw) ===
gtex_ids <- setdiff(grep("^GTEX-", colnames(expr_raw), value = TRUE),
                    Metadata$specimenID_portion)

if (length(gtex_ids) > 0) {
  meta_gtex <- tibble::tibble(
    specimenID         = gtex_ids,
    source             = "GTEX",
    case_id_12         = gtex_ids,
    sample_type        = "Solid Tissue Normal",
    HPV_type           = "Solid Tissue Normal",
    HPV_clade          = "Solid Tissue Normal",
    figo_stage         = "Solid Tissue Normal",
    primary_diagnosis  = "Solid Tissue Normal",
    race               = "Solid Tissue Normal",
    specimenID_portion = gtex_ids
  )
  Metadata <- dplyr::bind_rows(Metadata, meta_gtex)
}
saveRDS(Metadata, file = file.path(out_dir, "1_2_1_Metadata_full_TCGA_GTEX.rds"))
table(Metadata$source, Metadata$HPV_clade)
#       A7_clade A9_clade     Solid Tissue Normal
#GTEX        0        0                  10
#TCGA       66      202                   3

# --- Intersection and filtering ---
common_samples <- intersect(colnames(expr_raw)[-1], Metadata$specimenID_portion)
expr_raw_filt   <- expr_raw[, c("gene", common_samples), drop = FALSE]
dim(expr_raw_filt)
#[1] 60498   281

## 10) Save outputs
saveRDS(expr_raw_filt, file = file.path(out_dir, "1_2_2_Cervix_Toil_raw_counts.rds"))
#save.image(file = file.path(out_dir, "1_2_3_Image.RData"))

