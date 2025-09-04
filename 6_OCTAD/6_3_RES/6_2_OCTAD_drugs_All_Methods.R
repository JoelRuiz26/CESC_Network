# =========================================
# OCTAD pipeline for 3 DE methods (edgeR / DESeq2 / limma)
# - Reuses your original flow:
#   * computeCellLine (per clade) -> runsRGES (LINCS) -> filter sRGES
#   * octadDrugEnrichment (CHEMBL targets)
#   * topLineEval (per similar cell line)
# - Reads the significant signatures saved previously by method & clade:
#   6_1_0_Signatures/6_1_4_DE_<CLADE>_signific_<METHOD>.rds
# =========================================

suppressPackageStartupMessages({
  library(octad)
  library(octad.db)
  library(dplyr)
  library(vroom)
  library(tibble)
  library(signatureSearch)
})

setwd("~/CESC_Network/6_OCTAD/")

# ---------- Parameters (tune if needed) ----------
medcor_thr   <- 0.31     # threshold for selecting similar cell lines
sRGES_cutoff <- -0.25    # keep drugs with sRGES <= this value
permutations <- 10000
choose_fda   <- TRUE

# The 3 DE methods to run (must match the RDS filenames you saved):
methods <- c("DESeq2", "EdgeR", "limma")

# ---------- Load clade sample IDs ----------
sample_clado_A7 <- readRDS("6_0_A7_Samples.rds")
sample_clado_A9 <- readRDS("6_0_A9_Samples.rds")

# ---------- Metadata (ExperimentHub: EH7274) ----------
phenoDF <- octad.db::get_ExperimentHub_data("EH7274")
case_A7 <- intersect(sample_clado_A7, phenoDF$sample.id)
case_A9 <- intersect(sample_clado_A9, phenoDF$sample.id)

# ---------- (Optional) HDF5 for 'octad.whole' ----------
h5_path <- file.path(getwd(), "octad.counts.and.tpm.h5")
if (!file.exists(h5_path)) {
  url <- "https://chenlab-data-public.s3-us-west-2.amazonaws.com/octad/octad.counts.and.tpm.h5"
  message(">> Downloading HDF5 (~3 GB) to: ", h5_path)
  download.file(url, destfile = h5_path, mode = "wb", quiet = FALSE)
}
stopifnot(file.exists(h5_path))

# =========================================
# 1) Compute cell-line similarity (once per clade)
#    Using 'octad.small' as in your example
# =========================================
compute_lines_for_clade <- function(case_ids) {
  sim <- computeCellLine(case_id = case_ids, source = "octad.small")  # or source="octad.whole", file=h5_path
  sim <- dplyr::mutate(sim, cell_id = rownames(sim))
  keep <- dplyr::filter(sim, medcor > medcor_thr) |> dplyr::pull(cell_id)
  list(sim = sim, keep = keep)
}

sim_A7 <- compute_lines_for_clade(case_A7)
sim_A9 <- compute_lines_for_clade(case_A9)

message("A7 similar lines (medcor > ", medcor_thr, "): ", length(sim_A7$keep))
message("A9 similar lines (medcor > ", medcor_thr, "): ", length(sim_A9$keep))

# =========================================
# 2) Helpers to load signatures & run sRGES per method/clade
# =========================================

# Read one significant signature by clade & method
read_sig <- function(clade, method_tag) {
  f <- file.path("6_1_0_Signatures",
                 sprintf("6_1_4_DE_%s_signific_%s.rds", clade, method_tag))
  if (!file.exists(f)) stop("Signature not found: ", f)
  sig <- readRDS(f)
  # Minimal safety checks for runsRGES()
  need <- c("Symbol", "log2FoldChange", "regulation")
  if (!all(need %in% names(sig))) {
    stop("Signature missing required cols (Symbol, log2FoldChange, regulation): ", f)
  }
  sig
}

# Run sRGES for one clade/method
run_srges <- function(sig_df, lines_keep, sim_weights) {
  runsRGES(
    dz_signature     = sig_df,
    cells            = lines_keep,
    weight_cell_line = sim_weights,
    permutations     = permutations,
    choose_fda_drugs = choose_fda,
    output           = FALSE
  )
}

# Wrapper per clade to iterate the 3 methods end-to-end
do_clade_all_methods <- function(clade_label, case_ids, sim_obj) {
  # Storage for outputs per method
  out <- list()
  
  for (meth in methods) {
    message("---- ", clade_label, " | ", meth, " ----")
    
    # 2.1 Load signature for this method
    sig <- read_sig(clade_label, meth)
    
    # 2.2 sRGES (LINCS)
    res <- run_srges(sig, sim_obj$keep, sim_obj$sim)
    
    # 2.3 Filter sRGES (<= cutoff) and sort
    res_f <- res %>% dplyr::filter(sRGES <= sRGES_cutoff) %>% dplyr::arrange(sRGES)
    
    # 2.4 Save results (method-specific filenames)
    saveRDS(res,  file = sprintf("6_2_1_OCTAD_%s_%s_results_full.rds",   clade_label, meth))
    saveRDS(res_f, file = sprintf("6_2_1_OCTAD_%s_%s_results_%0.2f.rds", clade_label, meth, abs(sRGES_cutoff)))
    
    message(sprintf("%s | %s  -> full: %d rows | filtered (<= %.2f): %d rows",
                    clade_label, meth, nrow(res), sRGES_cutoff, nrow(res_f)))
    
    # 2.5 Enrichment (CHEMBL targets), keep separate folders per method
    enrich_dir <- sprintf("6_2_2_Enriquecimiento_%s_%s", clade_label, meth)
    octadDrugEnrichment(
      sRGES        = res_f,
      target_type  = "chembl_targets",
      enrichFolder = enrich_dir,
      outputFolder = getwd()
    )
    
    # Load enriched table back for quick inspection & filtered version (padj<0.05)
    enriched_path <- file.path(getwd(), enrich_dir, "chembl_targets", "enriched_chembl_targets.csv")
    if (file.exists(enriched_path)) {
      enriched_tbl <- vroom::vroom(enriched_path, show_col_types = FALSE)
      enriched_flt <- dplyr::filter(enriched_tbl, padj < 0.05)
    } else {
      warning("Enrichment output not found: ", enriched_path)
      enriched_tbl <- tibble()
      enriched_flt <- tibble()
    }
    
    # 2.6 In-silico validation per similar cell line (topLineEval)
    #     Uses a custom function you mentioned.
    if (file.exists("topLineEval_custom.R")) {
      source("topLineEval_custom.R")
      base_out <- file.path(getwd(), sprintf("6_2_3_Validation_%s_%s", clade_label, meth))
      if (!dir.exists(base_out)) dir.create(base_out, recursive = TRUE)
      
      for (cline in sim_obj$keep) {
        out_dir_line <- file.path(base_out, cline)
        if (!dir.exists(out_dir_line)) dir.create(out_dir_line, recursive = TRUE)
        message(sprintf("Validating %s | %s | line: %s", clade_label, meth, cline))
        tryCatch({
          topLineEval(
            topline      = cline,
            mysRGES      = res,           # full sRGES table (as in your code)
            outputFolder = out_dir_line
          )
        }, error = function(e) {
          message("  ⚠️  Error with ", cline, " : ", e$message)
        })
      }
    } else {
      message("Note: 'topLineEval_custom.R' not found; skipping validation step for ", clade_label, " | ", meth)
    }
    
    # Store in list
    out[[meth]] <- list(
      signature_file = sprintf("6_1_0_Signatures/6_1_4_DE_%s_signific_%s.rds", clade_label, meth),
      sRGES_full     = res,
      sRGES_filtered = res_f,
      enrich_all     = enriched_tbl,
      enrich_flt     = enriched_flt
    )
  }
  
  out
}

# =========================================
# 3) Execute for both clades (A7, A9)
# =========================================
res_A7_all <- do_clade_all_methods("A7", case_A7, sim_A7)
res_A9_all <- do_clade_all_methods("A9", case_A9, sim_A9)

# Save workspace snapshot
save.image("6_2_4_OCTAD_all_methods.RData")
#load("~/CESC_Network/6_OCTAD/6_2_4_OCTAD_all_methods.RData")

# -------- Optional quick summaries --------
for (meth in methods) {
  if (!is.null(res_A7_all[[meth]])) {
    cat("A7 |", meth, "| kept drugs (<= ", sRGES_cutoff, "): ",
        nrow(res_A7_all[[meth]]$sRGES_filtered), "\n", sep = "")
  }
  if (!is.null(res_A9_all[[meth]])) {
    cat("A9 |", meth, "| kept drugs (<= ", sRGES_cutoff, "): ",
        nrow(res_A9_all[[meth]]$sRGES_filtered), "\n", sep = "")
  }
}
