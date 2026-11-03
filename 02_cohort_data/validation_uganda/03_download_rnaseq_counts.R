### Joel Ruiz
### Download the actual RNA-seq gene-count files (GDC STAR - Counts, open
### access) for the Uganda validation population (HPV clade A7 or A9, an
### open-access RNA-seq file available), and collapse them to one file per
### patient. Several patients have more than one file (independent aliquots/
### libraries of the same tumour, not lanes of the same run -- their mapping
### stats differ too much for that); the one with the most gene-assigned
### reads (deepest, most usable library) is kept per patient.
###
### Output: raw/counts/<file_id>.tsv (one per KEPT file) and
###         1_Output/3_1_raw_gene_counts_A7A9.rds -- a gene x patient matrix
###         of "unstranded" raw counts, with a gene_id/gene_name key, ready
###         for 03_preprocessing/validation_uganda/01_normalise_rnaseq.R.
suppressPackageStartupMessages({ library(dplyr); library(vroom); library(purrr) })

setwd("~/CESC_Network/02_cohort_data/validation_uganda")
dir.create("raw/counts", showWarnings = FALSE, recursive = TRUE)

meta  <- vroom("1_Output/1_3_HTMCP_clinical_HPV_HIV_with_rnaseq_flag.tsv", show_col_types = FALSE) %>%
  filter(has_rnaseq, HPV.clade %in% c("A7", "A9"))
files <- vroom("1_Output/1_2_GDC_RNAseq_file_manifest.tsv", show_col_types = FALSE) %>%
  filter(submitter_id %in% meta$Patient)
cat(sprintf("[1] %d patients (A7/A9, RNA-seq), %d candidate GDC files\n",
            n_distinct(files$submitter_id), nrow(files)))

## ---- download every candidate file (idempotent: skips what is already there)
dl_one <- function(file_id) {
  dest <- file.path("raw/counts", paste0(file_id, ".tsv"))
  if (!file.exists(dest) || file.size(dest) < 1e5)
    download.file(paste0("https://api.gdc.cancer.gov/data/", file_id), dest,
                  mode = "wb", quiet = TRUE)
  dest
}
files$path <- map_chr(seq_len(nrow(files)), function(i) {
  if (i %% 20 == 0) cat(sprintf("  ...%d/%d\n", i, nrow(files)))
  dl_one(files$file_id[i])
})

## ---- read one STAR-Counts file: skip the 4 N_* rows, keep unstranded counts
read_star <- function(path) {
  d <- vroom(path, skip = 1, show_col_types = FALSE) %>%
    filter(!startsWith(gene_id, "N_"))
  d %>% select(gene_id, gene_name, gene_type, unstranded)
}

## ---- for a patient with > 1 file, keep the one with the most gene-assigned
## reads (deepest library); genuinely different aliquots, not lanes to sum
depth <- map_dbl(files$path, ~ sum(read_star(.x)$unstranded))
files$depth <- depth
kept <- files %>% group_by(submitter_id) %>% slice_max(depth, n = 1, with_ties = FALSE) %>% ungroup()
cat(sprintf("[2] %d patients had >1 file; kept the deepest library for each\n",
            sum(table(files$submitter_id) > 1)))
vroom_write(kept %>% select(submitter_id, file_id, file_name, depth),
           "1_Output/3_0_selected_rnaseq_file_per_patient.tsv", delim = "\t")

## ---- build the combined gene x patient raw-count matrix --------------------
## Every GDC STAR-Counts file uses the same GENCODE v36 gene model in the same
## row order (checked below), so this is a straight cbind of the "unstranded"
## column across patients -- not a join -- which is what makes this fast for
## 109 files where a full_join per file was not.
d1 <- read_star(kept$path[1])
gene_id <- d1$gene_id; gene_name <- d1$gene_name
mat <- matrix(NA_real_, nrow = length(gene_id), ncol = nrow(kept),
             dimnames = list(gene_id, kept$submitter_id))
mat[, 1] <- d1$unstranded
for (i in 2:nrow(kept)) {
  d <- read_star(kept$path[i])
  stopifnot(identical(d$gene_id, gene_id))   # same gene model/order, or stop
  mat[, i] <- d$unstranded
  if (i %% 20 == 0) cat(sprintf("  ...%d/%d\n", i, nrow(kept)))
}
cat(sprintf("[3] Combined matrix: %d genes x %d patients\n", nrow(mat), ncol(mat)))
saveRDS(list(counts = mat, gene_id = gene_id, gene_name = gene_name),
        "1_Output/3_1_raw_gene_counts_A7A9.rds")
