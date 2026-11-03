### Joel Ruiz
### Download metadata for the HTMCP-CC (Ugandan, HIV+/- cervical cancer) validation
### cohort -- data only, no preprocessing yet.
###
### Sources:
###   - GDC project CGCI-HTMCP-CC (dbGaP phs000528), open-access case list +
###     RNA-seq (STAR - Counts) file manifest, via the GDC REST API.
###   - Per-patient HPV type / clade / HIV status: Supplementary Table 2 of
###     Gagliardi et al. 2020, Nat Genet 52:800-810 (doi 10.1038/s41588-020-0673-7),
###     fetched from the publisher (static-content.springer.com); GDC's own
###     structured clinical fields (follow_ups.molecular_tests.hpv_strain) exist
###     in the schema but are NOT populated for this project.
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(jsonlite); library(vroom); library(httr)
})

setwd("~/CESC_Network/02_cohort_data/validation_uganda")
raw_dir <- "raw"; out_dir <- "1_Output"
dir.create(raw_dir, showWarnings = FALSE); dir.create(out_dir, showWarnings = FALSE)

## ---- 1) Supplementary Table 2 (HPV type/clade, HIV, clinical) --------------
supp_url <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-020-0673-7/MediaObjects/41588_2020_673_MOESM3_ESM.xlsx"
supp_path <- file.path(raw_dir, "Gagliardi2020_SupplementaryTables.xlsx")
if (!file.exists(supp_path)) download.file(supp_url, supp_path, mode = "wb", quiet = TRUE)

t2 <- read_excel(supp_path, sheet = "Supplementary Table 2", skip = 3)
names(t2) <- make.names(names(t2), unique = TRUE)
t2 <- t2 %>% filter(!is.na(Patient), Patient != "")
message(sprintf("[1] Supplementary Table 2: %d patients, %d columns", nrow(t2), ncol(t2)))
vroom_write(t2, file.path(out_dir, "1_1_HTMCP_clinical_HPV_HIV_table.tsv"), delim = "\t")

## ---- 2) GDC case list + open-access RNA-seq (STAR - Counts) manifest -------
gdc_post <- function(endpoint, body) {
  r <- POST(paste0("https://api.gdc.cancer.gov/", endpoint),
            body = toJSON(body, auto_unbox = TRUE), encode = "raw",
            add_headers(`Content-Type` = "application/json", `Accept` = "application/json"))
  fromJSON(content(r, "text", encoding = "UTF-8"), simplifyVector = TRUE)
}

cases <- gdc_post("cases", list(
  filters = list(op = "in", content = list(field = "project.project_id", value = "CGCI-HTMCP-CC")),
  fields  = "submitter_id,case_id",
  size    = 300))$data$hits
message(sprintf("[2] GDC cases in CGCI-HTMCP-CC: %d", nrow(cases)))

files <- gdc_post("files", list(
  filters = list(op = "and", content = list(
    list(op = "in", content = list(field = "cases.project.project_id", value = "CGCI-HTMCP-CC")),
    list(op = "in", content = list(field = "data_type", value = "Gene Expression Quantification")),
    list(op = "in", content = list(field = "experimental_strategy", value = "RNA-Seq")),
    list(op = "in", content = list(field = "access", value = "open")))),
  fields  = "file_id,file_name,cases.submitter_id,cases.samples.sample_type,analysis.workflow_type",
  size    = 300))$data$hits
files_flat <- files %>%
  mutate(submitter_id = sapply(cases, function(x) x$submitter_id[1]),
         sample_type  = sapply(cases, function(x) x$samples[[1]]$sample_type[1])) %>%
  select(file_id, file_name, submitter_id, sample_type)
message(sprintf("[3] Open-access RNA-seq (STAR - Counts) gene-count files: %d (%d unique patients)",
                nrow(files_flat), n_distinct(files_flat$submitter_id)))
vroom_write(files_flat, file.path(out_dir, "1_2_GDC_RNAseq_file_manifest.tsv"), delim = "\t")

## ---- 3) Join: which patients have both clinical/HPV data AND an RNA-seq file
joined <- t2 %>%
  mutate(has_rnaseq = Patient %in% files_flat$submitter_id)
vroom_write(joined, file.path(out_dir, "1_3_HTMCP_clinical_HPV_HIV_with_rnaseq_flag.tsv"), delim = "\t")

cat("\n=== Cohort ===\n"); print(table(joined$Cohort, useNA = "ifany"))
cat("\n=== HPV clade ===\n"); print(table(joined$HPV.clade, useNA = "ifany"))
cat("\n=== HIV status ===\n"); print(table(joined$HIV.status, useNA = "ifany"))
cat("\n=== HPV clade x HIV status ===\n"); print(table(joined$HPV.clade, joined$HIV.status, useNA = "ifany"))
cat("\n=== Has RNA-seq file in GDC (open access) ===\n"); print(table(joined$has_rnaseq, useNA = "ifany"))
cat("\n=== HPV clade, restricted to patients WITH an RNA-seq file ===\n")
print(table(joined$HPV.clade[joined$has_rnaseq], joined$HIV.status[joined$has_rnaseq]))
