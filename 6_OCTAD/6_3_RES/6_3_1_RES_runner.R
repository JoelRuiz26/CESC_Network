#============================================================================= 
# sRGES (RES) runner
# Inputs:
#   - Disease expression signatures (DGE) per method (DESeq2/EdgeR/limma)
#   - Patient sample IDs (present in TCGA/GTEx)
# Outputs:
#   - Filtered drug rankings per (signature × method):  RES_<sig>_<method>_cut*.rds
#   - FULL drug rankings per (signature × method):      RES_<sig>_<method>_FULL.rds
#   - NEW: Similar cell lines used per signature:       SimilarCellLines_<sig>_medcor>*.rds
#
# Steps:
#   1) Align case IDs to OCTAD metadata (EH7274)
#   2) Compute cell-line similarity per signature (octad.whole HDF5)
#   3) Run octad::runsRGES per (signature × method)
#   4) Save filtered + FULL rankings; NEW: save similar cell lines (character)
#
# Notes:
#   - Requires ~3GB HDF5 (octad.whole)
#   - "signatures" are patient cohorts; each may run multiple DGE methods
#=============================================================================

suppressPackageStartupMessages({
  library(octad)
  library(octad.db)
  library(dplyr)
  library(tibble)
})

# ------------------------- Parameters -------------------------
medcor_thr   <- 0.31     # threshold on median correlation (medcor) to keep similar cell lines
sRGES_cutoff <- -0.20    # filtering cutoff for sRGES (more negative is better)
permutations <- 10000     # permutations for runsRGES
choose_fda   <- TRUE       # restrict to FDA drugs if TRUE
set.seed(123)

# ------------------------------ Paths ------------------------
base_dir <- "~/CESC_Network/6_OCTAD"
h5_dir   <- file.path(base_dir, "6_1_DGE_signature")
out_dir  <- file.path(base_dir, "6_3_RES")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# HDF5 (octad.whole)
h5_file  <- file.path(h5_dir, "octad.counts.and.tpm.h5")
h5_url   <- "https://chenlab-data-public.s3-us-west-2.amazonaws.com/octad/octad.counts.and.tpm.h5"
if (!file.exists(h5_file)) {
  message(">> Downloading HDF5 (~3 GB) to: ", h5_file)
  dir.create(dirname(h5_file), recursive = TRUE, showWarnings = FALSE)
  download.file(h5_url, destfile = h5_file, mode = "wb", quiet = FALSE)
}
if (!file.exists(h5_file)) stop("HDF5 not found: ", h5_file)

# ------------------------------ USER CONFIG -------------------
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
)

# ---------------------- Align case IDs with OCTAD --------------
# EH7274 contains OCTAD phenotype metadata; used to keep only aligned sample IDs
phenoDF <- octad.db::get_ExperimentHub_data("EH7274")

# =============================== FUNCTIONS ====================

# (1) Compute cell-line similarity for one signature
#     Returns:
#       - sim_tbl: data.frame with medcor per cell line (sorted desc)
#       - cells:   character vector of cell_ids with medcor > medcor_thr
compute_similarity_for_signature <- function(case_ids, medcor_thr, h5_file) {
  sim_tbl <- octad::computeCellLine(case_id = case_ids, source = "octad.small")#source = "octad.whole", file = h5_file)
  sim_tbl <- sim_tbl %>%
    dplyr::mutate(cell_id = rownames(sim_tbl)) %>%
    dplyr::arrange(dplyr::desc(medcor))
  cells <- sim_tbl %>% dplyr::filter(medcor > medcor_thr) %>% dplyr::pull(cell_id)
  list(sim_tbl = sim_tbl, cells = cells)
}

# (2) Run runsRGES for one (signature × method) and save FILTERED + FULL tables
run_res_for_signature <- function(signature, method, sig_path, sim_tbl, cells,
                                  permutations, choose_fda, sRGES_cutoff, out_dir) {
  if (!file.exists(sig_path)) stop("DE signature RDS not found: ", sig_path)
  
  # Load disease signature; keep unique gene Symbols only
  dz_sig <- readRDS(sig_path) %>%
    dplyr::filter(!is.na(Symbol)) %>%
    dplyr::distinct(Symbol, .keep_all = TRUE)
  
  # Ensure supporting data is available (as in OCTAD pipelines)
  suppressMessages({
    invisible(octad.db::get_ExperimentHub_data("EH7270"))
    invisible(octad.db::get_ExperimentHub_data("EH7271"))
    if (isTRUE(choose_fda)) invisible(octad.db::get_ExperimentHub_data("EH7269"))
  })
  
  # Output paths
  sig_dir   <- file.path(out_dir, signature, method)
  dir.create(sig_dir, recursive = TRUE, showWarnings = FALSE)
  out_kept  <- file.path(sig_dir, sprintf("RES_%s_%s_cut%.2f.rds",  signature, method, sRGES_cutoff))
  out_full  <- file.path(sig_dir, sprintf("RES_%s_%s_FULL.rds",     signature, method))   # NEW
  
  # If no similar cell lines pass the threshold, save empty placeholders
  if (length(cells) == 0) {
    message(sprintf(">> [%s - %s] 0 cells over medcor_thr=%.2f. Saving empty.", signature, method, medcor_thr))
    saveRDS(tibble::tibble(), out_kept)
    saveRDS(tibble::tibble(), out_full)   # NEW
    return(invisible(tibble::tibble()))
  }
  
  # Run sRGES
  message(sprintf(">> [%s - %s] runsRGES on %d cells...", signature, method, length(cells)))
  res <- octad::runsRGES(
    dz_signature      = dz_sig,
    cells             = cells,
    weight_cell_line  = sim_tbl,
    permutations      = permutations,
    choose_fda_drugs  = choose_fda,
    output            = TRUE,            # also writes CSVs to sig_dir
    outputFolder      = sig_dir
  )
  
  # Save FULL ranking (unfiltered)
  res_tbl <- tibble::as_tibble(res)
  saveRDS(res_tbl, out_full)
  
  # Save filtered table (sRGES <= cutoff), augmented with signature/method
  kept <- res_tbl %>%
    dplyr::filter(sRGES <= sRGES_cutoff) %>%
    dplyr::arrange(sRGES) %>%
    dplyr::mutate(signature = signature, method = method)
  
  saveRDS(kept, out_kept)
  message(sprintf(">> [%s - %s] saved: %s (n=%d) + FULL: %s (n=%d)",
                  signature, method, out_kept, nrow(kept), out_full, nrow(res_tbl)))
  invisible(kept)
}

# ============================== MAIN LOOP =====================
for (signature in names(sig_cfg)) {
  
  # 1) Load cohort case IDs and keep only those present in OCTAD
  cases_file <- sig_cfg[[signature]]$cases_rds
  if (!file.exists(cases_file)) {
    warning("Missing cases_rds for signature: ", signature, " -> ", cases_file)
    next
  }
  case_ids_raw <- readRDS(cases_file)
  case_ids     <- intersect(case_ids_raw, phenoDF$sample.id)
  message(sprintf("\n== Signature %s: %d raw cases | %d aligned in OCTAD",
                  signature, length(case_ids_raw), length(case_ids)))
  
  # 2) Compute cell-line similarity (returns both the table and the filtered vector)
  sim_out <- compute_similarity_for_signature(case_ids, medcor_thr, h5_file)
  message(sprintf("   Selected cell lines: %d (medcor > %.2f)", length(sim_out$cells), medcor_thr))
  
  # --- NEW: Save the "pull" of similar cell lines (character) per signature ---
  #     This does NOT affect downstream results; it only adds a small side artifact.
  sig_base_dir <- file.path(out_dir, signature)
  dir.create(sig_base_dir, recursive = TRUE, showWarnings = FALSE)
  similar_cells_chr <- as.character(sim_out$cells)  # already ordered by medcor desc
  out_cells_rds <- file.path(sig_base_dir, sprintf("SimilarCellLines_%s_medcor_gt%.2f.rds", signature, medcor_thr))
  saveRDS(similar_cells_chr, out_cells_rds)
  message(sprintf("   Saved similar cell lines (character): %s [n=%d]",
                  out_cells_rds, length(similar_cells_chr)))
  # (Optional) Save the full similarity table:
  # out_simtbl_rds <- file.path(sig_base_dir, sprintf("CellLineSimilarityTable_%s_FULL.rds", signature))
  # saveRDS(sim_out$sim_tbl, out_simtbl_rds)
  # ---------------------------------------------------------------------------
  
  # 3) Run sRGES per DGE method using the same similarity outputs
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

message("\n>> DONE. Filtered + FULL rankings under: ", out_dir)
