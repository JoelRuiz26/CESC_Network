###############################################################################
# Build TCGA-CESC metadata + HPV clade
# - Inputs:
#     * HPV_TCGA TSV: HPV_type of TCGA-CESC project were obtained from literature (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7459114/)
#                     supplementary PDF of doi: 10.1038/s41598-020-71300-7. PMID: 32868873; PMCID: PMC7459114.
#     * Output_query_TCGA: a data.frame/tibble with at least:
#           cases.submitter_id, sample_type, cases (specimen/sample barcode)
# - Outputs:
#     * metadata TSV
#     * A7 / A9 sample vectors as RDS in ~/CESC_Network/1_Get_Data/
###############################################################################
suppressPackageStartupMessages({
  library(dplyr)
  library(vroom)
  library(TCGAbiolinks)
  library(tidyr)
})

# --- Paths ---------------------------------------------------------------------
hpv_path   <- "~/CESC_Network/0_HPV_Distribution/1_1_HPV_Clade_Distrib.tsv"
out_meta   <- "~/CESC_Network/1_Get_Data/1_1_1_Metadata.tsv"
out_dir_id <- "~/CESC_Network/1_Get_Data/"
dir.create(path.expand(dirname(out_meta)), recursive = TRUE, showWarnings = FALSE)
dir.create(path.expand(out_dir_id),        recursive = TRUE, showWarnings = FALSE)

# --- 1) HPV clade (literature) -------------------------------------------------
# Expected columns: `TCGA Case ID` | tipos_HPV | Clado_filogenetico | frec | percent
HPV_TCGA <- vroom::vroom(hpv_path, show_col_types = FALSE)

# --- 2) TCGA clinical (per-case; includes FIGO if available) -------------------
clinical_tcga <- TCGAbiolinks::GDCquery_clinic("TCGA-CESC", "clinical") %>%
  dplyr::select(any_of(c("submitter_id", "figo_stage", "race", "primary_diagnosis")))

# --- 3) GDC sample barcodes for TCGA-CESC (per-sample) -------------------------
# We only retrieve metadata (no downloads).
# Columns we want: case submitter id (12-char), sample submitter id (barcode), sample_type
qry <- TCGAbiolinks::GDCquery(
  project       = "TCGA-CESC",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification"
)
res <- TCGAbiolinks::getResults(qry)

# Harmonize column names across possible GDC versions
# Try to extract: case_id_12 = cases.submitter_id ; specimenID = sample.submitter_id ; sample_type
case_col_candidates    <- c("cases.submitter_id", "submitter_id")
sample_col_candidates  <- c("sample.submitter_id", "sample_submitter_id", "cases")  # fallback "cases" in older tables
type_col_candidates    <- c("sample_type", "sample.type")

grab_first <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) return(NULL)
  df[[hit[1]]]
}

samples_tbl <- tibble(
  case_id_12 = grab_first(res, case_col_candidates),
  specimenID = grab_first(res, sample_col_candidates),
  sample_type = grab_first(res, type_col_candidates)
) %>%
  dplyr::filter(!is.na(case_id_12), !is.na(specimenID)) %>%
  distinct(case_id_12, specimenID, .keep_all = TRUE)

# --- 4) Join HPV (per-case) + clinical (per-case) onto sample-level barcodes ---
# Map per-case info to each sample via case_id_12.
meta_joined <- samples_tbl %>%
  left_join(clinical_tcga, by = c("case_id_12" = "submitter_id")) %>%
  left_join(
    HPV_TCGA %>%
      dplyr::select(`TCGA Case ID`, tipos_HPV, Clado_filogenetico, any_of(c("frec", "percent"))),
    by = c("case_id_12" = "TCGA Case ID")
  ) %>%
  # If a sample is "Solid Tissue Normal", mark staging/clade fields accordingly
  mutate(
    figo_stage          = ifelse(sample_type == "Solid Tissue Normal", "Solid Tissue Normal", figo_stage),
    tipos_HPV           = ifelse(sample_type == "Solid Tissue Normal", "Solid Tissue Normal", tipos_HPV),
    Clado_filogenetico  = ifelse(sample_type == "Solid Tissue Normal", "Solid Tissue Normal", Clado_filogenetico)
  ) %>%
  # Final tidy column names
  dplyr::rename(
    HPV_type  = tipos_HPV,
    HPV_clade = Clado_filogenetico
  ) %>%
  # Order columns for readability
  dplyr::select(
    specimenID, case_id_12, sample_type,
    HPV_type, HPV_clade, figo_stage,
    primary_diagnosis, race
  ) %>%
  # === Apply all filters and recoding here in one clean block ===
  filter(
    !grepl("metast", sample_type, ignore.case = TRUE),      # remove metastasis
    !(HPV_clade %in% "otro"),                              # drop "otro"
    !(HPV_type %in% "HPV70"),                              # drop HPV70
    !is.na(HPV_type)                                       # drop missing HPV type
  ) %>%
  mutate(
    HPV_clade = recode(HPV_clade, "A7" = "A7_clade", "A9" = "A9_clade"),
    source = if (exists("gtex_samples")) 
      if_else(specimenID %in% gtex_samples, "GTEX", "TCGA") 
    else "TCGA"
  ) %>%
  dplyr::select(specimenID, source, everything())

# --- 5) Quick clade counts (A7 / A9 / others) ---------------------------------
clade_counts <- meta_joined %>%
  group_by(HPV_clade) %>%
  summarise(n = n(), .groups = "drop")
print(clade_counts)
# A tibble: 5 × 2
#HPV_clade                n
#<chr>                  <int>
#1 A7                     68
#2 A9                    203
#3 Solid Tissue Normal     3
#4 otro                    5
#5 NA                     30


# --- 6) Save metadata ----------------------------------------------------------
vroom::vroom_write(meta_joined, out_meta, delim = "\t")

A7_ids <- meta_joined %>%
  filter(HPV_clade %in% "A7_clade") %>%
  pull(specimenID) %>%
  unique()

A9_ids <- meta_joined %>%
  filter(HPV_clade %in% "A9_clade") %>%
  pull(specimenID) %>%
  unique()

saveRDS(A7_ids, file = file.path(out_dir_id, "1_1_2_A7_Samples.rds"))
saveRDS(A9_ids, file = file.path(out_dir_id, "1_1_2_A9_Samples.rds"))
