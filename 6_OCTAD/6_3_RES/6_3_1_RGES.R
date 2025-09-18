suppressPackageStartupMessages({
  library(octad); library(octad.db); library(dplyr); library(tibble)
})

# ---- Parameters ----
medcor_thr   <- 0.3       # minimum median correlation to keep a cell line
permutations <- 10000     # number of permutations for runsRGES
choose_fda   <- TRUE      # restrict to FDA drugs if TRUE (uses EH7269)
set.seed(123)

# ---- Paths ----
base_dir <- "~/CESC_Network/6_OCTAD"
out_dir  <- file.path(base_dir, "6_3_RES")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Per-signature config: cohort file + DGE result per method ----
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

# ---- OCTAD phenotype metadata (to align sample IDs) ----
phenoDF <- octad.db::get_ExperimentHub_data("EH7274")

# ---- Computes cell-line similarity for one signature ----
# Returns:
#   - sim_tbl: data.frame with medcor per cell line (sorted desc)
#   - cells: character vector of cell_ids with medcor > medcor_thr
compute_similarity_for_signature <- function(case_ids, medcor_thr) {
  # source = "octad.small" uses precomputed resources (no H5 file required)
  sim_tbl <- octad::computeCellLine(case_id = case_ids, source = "octad.small")
  sim_tbl <- sim_tbl %>% mutate(cell_id = rownames(sim_tbl)) %>% arrange(desc(medcor))
  cells <- sim_tbl %>% filter(medcor > medcor_thr) %>% pull(cell_id)
  list(sim_tbl = sim_tbl, cells = cells)
}

# ---- Run runsRGES for one (signature × method); save ONLY FULL.rds ----
run_res_for_signature <- function(signature, method, sig_path, sim_tbl, cells,
                                  permutations, choose_fda, out_dir) {
  if (!file.exists(sig_path)) stop("DE signature RDS not found: ", sig_path)
  
  # Load disease signature and keep unique gene Symbols
  dz_sig <- readRDS(sig_path) %>% filter(!is.na(Symbol)) %>% distinct(Symbol, .keep_all = TRUE)
  
  # Ensure supporting resources are available (LINCS metadata and signatures; FDA list if requested)
  suppressMessages({
    invisible(octad.db::get_ExperimentHub_data("EH7270")) # lincs_sig_info
    invisible(octad.db::get_ExperimentHub_data("EH7271")) # lincs_signatures
    if (isTRUE(choose_fda)) invisible(octad.db::get_ExperimentHub_data("EH7269")) # DRH/FDA
  })
  
  # Output location for this signature × method
  sig_dir  <- file.path(out_dir, signature, method)
  dir.create(sig_dir, recursive = TRUE, showWarnings = FALSE)
  out_full <- file.path(sig_dir, sprintf("RES_%s_%s_FULL.rds", signature, method))
  
  # If no similar cell lines pass the threshold, save an empty FULL and return
  if (length(cells) == 0) {
    message(sprintf(">> [%s - %s] 0 cells > medcor %.2f. Saving empty FULL.", signature, method, medcor_thr))
    saveRDS(tibble(), out_full); return(invisible(tibble()))
  }
  
  # Run sRGES (output=TRUE keeps OCTAD’s CSV artifacts such as sRGES_FDAapproveddrugs.csv)
  message(sprintf(">> [%s - %s] runsRGES on %d cells...", signature, method, length(cells)))
  res <- octad::runsRGES(
    dz_signature     = dz_sig,
    cells            = cells,
    weight_cell_line = sim_tbl,
    permutations     = permutations,
    choose_fda_drugs = choose_fda,
    output           = TRUE,
    outputFolder     = sig_dir
  )
  
  # Save ONLY the unfiltered ranking as FULL.rds
  res_tbl <- tibble::as_tibble(res)
  saveRDS(res_tbl, out_full)
  message(sprintf(">> [%s - %s] saved FULL: %s (n=%d)", signature, method, out_full, nrow(res_tbl)))
  invisible(res_tbl)
}

# ---- Main loop over signatures and methods ----
for (signature in names(sig_cfg)) {
  # Load cohort IDs and align them to OCTAD phenotype table
  cases_file <- sig_cfg[[signature]]$cases_rds
  if (!file.exists(cases_file)) { warning("Missing cases_rds: ", cases_file); next }
  case_ids_raw <- readRDS(cases_file)
  case_ids     <- intersect(case_ids_raw, phenoDF$sample.id)
  message(sprintf("\n== Signature %s: %d raw | %d aligned", signature, length(case_ids_raw), length(case_ids)))
  
  # Compute similarity and select cell lines
  sim_out <- compute_similarity_for_signature(case_ids, medcor_thr)
  message(sprintf("   Selected cell lines: %d (medcor > %.2f)", length(sim_out$cells), medcor_thr))
  
  # --- NEW: Save the selected similar cell lines as a simple character vector ---
  # This mirrors your old behavior; file name includes the exact medcor threshold.
  sig_base_dir <- file.path(out_dir, signature)
  dir.create(sig_base_dir, recursive = TRUE, showWarnings = FALSE)
  similar_cells_chr <- as.character(sim_out$cells)  # already ordered by medcor desc
  out_cells_rds <- file.path(sig_base_dir, sprintf("SimilarCellLines_%s_medcor_gt%.2f.rds", signature, medcor_thr))
  saveRDS(similar_cells_chr, out_cells_rds)
  message(sprintf("   Saved similar cell lines: %s [n=%d]", out_cells_rds, length(similar_cells_chr)))
  # (Optional) If you also want the full similarity table per signature, uncomment:
  # out_simtbl_rds <- file.path(sig_base_dir, sprintf("CellLineSimilarityTable_%s_FULL.rds", signature))
  # saveRDS(sim_out$sim_tbl, out_simtbl_rds)
  
  # Run sRGES for each DGE method using the same similarity results
  for (method in names(sig_cfg[[signature]]$expr)) {
    run_res_for_signature(signature, method, sig_cfg[[signature]]$expr[[method]],
                          sim_out$sim_tbl, sim_out$cells,
                          permutations, choose_fda, out_dir)
  }
}

message("\n>> DONE. FULL rankings (only .rds saved per method) under: ", out_dir)
