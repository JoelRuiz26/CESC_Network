# ============================
# Batch OCTAD drug enrichment + n_present per set
# ============================
suppressPackageStartupMessages({
  library(octad)
  library(octad.db)   # for get_ExperimentHub_data()
  library(vroom)
  library(dplyr)
  library(purrr)
  library(readr)
})

#------------------------------
# Helper: load cmpd sets by target_type
# Returns a named list of vectors (drug sets), e.g. names = MeSH terms
#------------------------------
load_cmpd_sets <- function(target_type = c("chembl_targets","mesh","ChemCluster")){
  target_type <- match.arg(target_type)
  # EH ids match the source you pasted from octadDrugEnrichment()
  eh_map <- c(ChemCluster="EH7266", chembl_targets="EH7267", mesh="EH7268")
  x <- octad.db::get_ExperimentHub_data(eh_map[[target_type]])
  sets <- x$cmpd.sets
  names(sets) <- x$cmpd.set.names
  sets
}

#------------------------------
# Helper: compute n_present per set for ONE signature (one sRGES table)
# - srg: data.frame with at least pert_iname, sRGES
# - target_type: one of "chembl_targets","mesh","ChemCluster"
# Returns: data.frame(target, n_present, target_type, signature)
#------------------------------
n_present_per_set <- function(srg, target_type, signature){
  stopifnot(all(c("pert_iname","sRGES") %in% names(srg)))
  srg <- srg %>% mutate(pert_iname = toupper(trimws(pert_iname))) %>% distinct(pert_iname)
  sets <- load_cmpd_sets(target_type)
  sets <- lapply(sets, function(v) toupper(trimws(v)))
  tibble(
    target      = names(sets),
    n_present   = vapply(sets, function(v) sum(srg$pert_iname %in% v), integer(1)),
    target_type = target_type,
    signature   = signature
  )
}

#------------------------------
# Your original enrich_many(), unchanged in behavior
#------------------------------
enrich_many <- function(sRGES_list,
                        out_dir,
                        target_types = c("chembl_targets", "mesh", "ChemCluster"),
                        alpha = 0.01,
                        folder_prefix = "Enriquecimiento") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (is.data.frame(sRGES_list)) {
    sRGES_list <- list(item1 = sRGES_list)
  } else if (is.list(sRGES_list) && is.null(names(sRGES_list))) {
    names(sRGES_list) <- paste0("item", seq_along(sRGES_list))
  }
  
  csv_path <- function(base, sub, ttype) {
    file.path(base, sub, ttype, paste0("enriched_", ttype, ".csv"))
  }
  
  all_results <- list()
  sig_results <- list()
  
  for (nm in names(sRGES_list)) {
    srg <- sRGES_list[[nm]]
    stopifnot(all(c("pert_iname","sRGES") %in% colnames(srg)))
    enrich_sub <- paste0(folder_prefix, "_", nm)
    
    for (tt in target_types) {
      message(">> Enrichment: ", nm, " | target_type=", tt)
      octadDrugEnrichment(
        sRGES        = srg,
        target_type  = tt,
        enrichFolder = enrich_sub,
        outputFolder = out_dir,
        outputRank   = TRUE
      )
      
      fp <- csv_path(out_dir, enrich_sub, tt)
      if (!file.exists(fp)) {
        warning("   CSV not found (skipping): ", fp)
        next
      }
      
      df <- vroom::vroom(fp, show_col_types = FALSE) |>
        mutate(signature = nm, target_type = tt)
      
      all_results[[paste(nm, tt, sep = "_")]] <- df
      sig_results[[paste(nm, tt, sep = "_")]] <- df %>% filter(padj <= alpha)
    }
  }
  
  list(
    all = dplyr::bind_rows(all_results),
    sig = dplyr::bind_rows(sig_results)
  )
}

# --------------------
# Load your sRGES (A7, A9) and sort ascending (more negative first)
# --------------------

sRES_A7 <- readRDS("~/CESC_Network/6_OCTAD/6_3_RES/A7/RES_A7_common3_collapsed_FDA_Launched.rds") %>%
  arrange(sRGES)
sRES_A9 <- readRDS("~/CESC_Network/6_OCTAD/6_3_RES/A9/RES_A9_common3_collapsed_FDA_Launched.rds") %>%
  arrange(sRGES)
head(sRES_A7)

out_base <- "~/CESC_Network/6_OCTAD/6_4_EnrichedDrugs"

# Run enrichment 
res_enrich <- enrich_many(
  sRGES_list   = list(A7 = sRES_A7, A9 = sRES_A9),
  out_dir      = out_base,
  target_types = c("chembl_targets", "mesh", "ChemCluster"),
  alpha        = 0.01,
  folder_prefix= "Enriquecimiento"
)

# --------------------
# NEW: attach n_present per (signature, target_type, target)
# --------------------


# 1) Build an index of n_present for each signature and target_type
present_idx <- bind_rows(
  # A7
  n_present_per_set(sRES_A7, "chembl_targets", "A7"),
  n_present_per_set(sRES_A7, "mesh",            "A7"),
  n_present_per_set(sRES_A7, "ChemCluster",     "A7"),
  # A9
  n_present_per_set(sRES_A9, "chembl_targets", "A9"),
  n_present_per_set(sRES_A9, "mesh",            "A9"),
  n_present_per_set(sRES_A9, "ChemCluster",     "A9")
)

# 2) Join n_present to your enrichment tables
enrich_all_n <- res_enrich$all %>%
  left_join(present_idx, by = c("target","target_type","signature"))
# A tibble: 1 × 8
#n   min   q25   med   q75   max n_gt_0_5 n_gt_0_6
#<int> <dbl> <dbl> <dbl> <dbl> <dbl>    <int>    <int>
#  1   195 0.110 0.214 0.265 0.346 0.498        0        0

enrich_sig_n <- res_enrich$sig %>%
  left_join(present_idx, by = c("target","target_type","signature"))

# --------------------
# Inspect: same table as before but with n_present appended
# --------------------
dim(enrich_all_n); dim(enrich_sig_n)

###   0.3 and -0.20 RGES
#[1] 2356    7
#[1] 330   7

# Save if you want
saveRDS(enrich_all_n, file.path(out_base, "6_4_1_enrichment_all_with_n_present.rds"))
saveRDS(enrich_sig_n, file.path(out_base, "6_4_1_enrichment_sig_with_n_present.rds"))


