# =============================================================================
#  sRGES (RES) runner

# Input: -list of DGE (define the method used) and 
#       -patient sample IDs (included in TCGA/GTEX)
# Output: -List of drugs with RES significant scoring

# What it does
#   1) Aligns case IDs to OCTAD metadata (EH7274)
#   2) Computes cell-line similarity per signature using octad.whole (HDF5)
#   3) Runs octad::runsRGES for each (signature × method)
#   4) Savesthe filtered RDS per combo: RES_<signature>_<method>_cut*.rds
#

# Notes
#   - Requires the ~3GB HDF5 file for octad.whole.
# - You define "signatures" = sets of patient sample IDs (case cohorts).
# - For each signature, you specify one or more expression signatures (methods).
# =============================================================================

suppressPackageStartupMessages({
  library(octad)
  library(octad.db)
  library(dplyr)
  library(tibble)
})

# ------------------------- Parameters (edit if needed) ------------------------
medcor_thr   <- 0.35     # select cell lines with medcor > this threshold
sRGES_cutoff <- -0.20    # keep drugs with sRGES <= cutoff
permutations <- 10000    # permutations for runsRGES
choose_fda   <- TRUE     # restrict to FDA-approved compounds
set.seed(123)

# ------------------------------ Paths ----------------------------------------
base_dir <- "~/CESC_Network/6_OCTAD"
h5_dir   <- file.path(base_dir, "6_1_DGE_signature")
out_dir  <- file.path(base_dir, "6_3_RES")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# HDF5 for octad.whole
h5_file  <- file.path(h5_dir, "octad.counts.and.tpm.h5")
h5_url   <- "https://chenlab-data-public.s3-us-west-2.amazonaws.com/octad/octad.counts.and.tpm.h5"
if (!file.exists(h5_file)) {
  message(">> Downloading HDF5 (~3 GB) to: ", h5_file)
  dir.create(dirname(h5_file), recursive = TRUE, showWarnings = FALSE)
  download.file(h5_url, destfile = h5_file, mode = "wb", quiet = FALSE)
}
if (!file.exists(h5_file)) stop("HDF5 not found: ", h5_file)

# ------------------------------ USER CONFIG ----------------------------------
# Map each signature to one or more DE signatures (methods).
# You can have 1 method, 3 methods, or 20—this loop will handle it.
sig_cfg <- list(
  A7 = list(
    cases_rds = file.path(base_dir, "6_0_A7_Samples.rds"),
    expr = list(
      DESeq2 = file.path(base_dir, "6_1_DGE_signature/6_1_0_Output_rds/6_1_4_DE_A7_signific_DESeq2.rds"),
      EdgeR  = file.path(base_dir, "6_1_DGE_signature/6_1_0_Output_rds/6_1_4_DE_A7_signific_EdgeR.rds"),
      limma  = file.path(base_dir, "6_1_DGE_signature/6_1_0_Output_rds/6_1_4_DE_A7_signific_limma.rds")
    )
  ),
  A9 = list(
    cases_rds = file.path(base_dir, "6_0_A9_Samples.rds"),
    expr = list(
      DESeq2 = file.path(base_dir, "6_1_DGE_signature/6_1_0_Output_rds/6_1_4_DE_A9_signific_DESeq2.rds"),
      EdgeR  = file.path(base_dir, "6_1_DGE_signature/6_1_0_Output_rds/6_1_4_DE_A9_signific_EdgeR.rds"),
      limma  = file.path(base_dir, "6_1_DGE_signature/6_1_0_Output_rds/6_1_4_DE_A9_signific_limma.rds")
    )
  )
  # Add more signatures as needed...
)

# ---------------------- Align case IDs with OCTAD metadata -------------------
# EH7274 has the sample IDs OCTAD expects (TCGA/GTEx/etc.)
phenoDF <- octad.db::get_ExperimentHub_data("EH7274")

# =============================== FUNCTIONS ===================================

# (1) Compute cell-line similarity for one signature
#     Returns:
#       - sim_tbl  : data.frame with medcor and an added 'cell_id' column
#       - cells    : character vector of selected cell IDs (medcor > threshold)
compute_similarity_for_signature <- function(case_ids, medcor_thr, h5_file) {
  sim_tbl <- octad::computeCellLine(case_id = case_ids, source = "octad.whole", file = h5_file)
  sim_tbl <- sim_tbl %>%
    dplyr::mutate(cell_id = rownames(sim_tbl)) %>%
    dplyr::arrange(dplyr::desc(medcor))
  cells <- sim_tbl %>% dplyr::filter(medcor > medcor_thr) %>% dplyr::pull(cell_id)
  list(sim_tbl = sim_tbl, cells = cells)
}

# (2) Run runsRGES for one (signature × method) and save filtered results
run_res_for_signature <- function(signature, method, sig_path, sim_tbl, cells,
                                  permutations, choose_fda, sRGES_cutoff, out_dir) {
  if (!file.exists(sig_path)) stop("DE signature RDS not found: ", sig_path)
  
  # Load DE signature as-is (your schema already matches what runsRGES expects)
  dz_sig <- readRDS(sig_path) %>%
    dplyr::filter(!is.na(Symbol)) %>%
    dplyr::distinct(Symbol, .keep_all = TRUE)
  
  # Warm up EH datasets (prevents 'lincs_drug_prediction not found' in some setups)
  suppressMessages({
    invisible(octad.db::get_ExperimentHub_data("EH7270"))  # LINCS metadata
    invisible(octad.db::get_ExperimentHub_data("EH7271"))  # LINCS signatures
    if (isTRUE(choose_fda)) invisible(octad.db::get_ExperimentHub_data("EH7269"))  # FDA list
  })
  
  sig_dir  <- file.path(out_dir, signature, method)
  dir.create(sig_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(sig_dir, sprintf("RES_%s_%s_cut%.2f.rds", signature, method, sRGES_cutoff))
  
  if (length(cells) == 0) {
    message(sprintf(">> [%s - %s] 0 cells over medcor_thr=%.2f. Saving empty.", signature, method, medcor_thr))
    saveRDS(tibble::tibble(), out_file)
    return(invisible(tibble::tibble()))
  }
  
  message(sprintf(">> [%s - %s] runsRGES on %d cells...", signature, method, length(cells)))
  res <- octad::runsRGES(
    dz_signature      = dz_sig,
    cells             = cells,
    weight_cell_line  = sim_tbl,     # OK: rownames still present after mutate()
    permutations      = permutations,
    choose_fda_drugs  = choose_fda,
    output            = TRUE,        # also writes the helpful CSVs in sig_dir
    outputFolder      = sig_dir
  )
  
  kept <- tibble::as_tibble(res) %>%
    dplyr::filter(sRGES <= sRGES_cutoff) %>%
    dplyr::arrange(sRGES) %>%
    dplyr::mutate(signature = signature, method = method)
  
  saveRDS(kept, out_file)
  message(sprintf(">> [%s - %s] saved: %s (n=%d)", signature, method, out_file, nrow(kept)))
  invisible(kept)
}

# ============================== MAIN WORKFLOW ================================
for (signature in names(sig_cfg)) {
  
  # Load and align case IDs to OCTAD’s metadata
  cases_file <- sig_cfg[[signature]]$cases_rds
  if (!file.exists(cases_file)) {
    warning("Missing cases_rds for signature: ", signature, " -> ", cases_file)
    next
  }
  case_ids_raw <- readRDS(cases_file)
  case_ids     <- intersect(case_ids_raw, phenoDF$sample.id)
  message(sprintf("\n== Signature %s: %d raw cases | %d aligned in OCTAD",
                  signature, length(case_ids_raw), length(case_ids)))
  
  # (1) Similarity for this signature
  sim_out <- compute_similarity_for_signature(case_ids, medcor_thr, h5_file)
  message(sprintf("   Selected cell lines: %d (medcor > %.2f)", length(sim_out$cells), medcor_thr))
  
  # (2) RES for each method mapped to this signature
  for (method in names(sig_cfg[[signature]]$expr)) {
    sig_path <- sig_cfg[[signature]]$expr[[method]]
    run_res_for_signature(
      signature    = signature,
      method       = method,
      sig_path     = sig_path,
      sim_tbl      = sim_out$sim_tbl,
      cells        = sim_out$cells,
      permutations = permutations,
      choose_fda   = choose_fda,
      sRGES_cutoff = sRGES_cutoff,
      out_dir      = out_dir
    )
  }
}

message("\n>> DONE. Outputs under: ", out_dir)
