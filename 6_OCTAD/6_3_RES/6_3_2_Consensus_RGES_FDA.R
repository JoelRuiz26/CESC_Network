# ============================================================
# Collapse sRGES (FULL and FDA-Launched) per signature (A7/A9)
# Requires per method:
#   base_res/<sig>/<method>/RES_<sig>_<method>_FULL.rds
#   base_res/<sig>/<method>/sRGES_FDAapproveddrugs.csv
# Outputs (per signature):
#   RES_<sig>_common_3methods_FULL.rds
#   RES_<sig>_common3_collapsed_FULL.rds
#   RES_<sig>_common_3methods_FDA_Launched.rds
#   RES_<sig>_common3_collapsed_FDA_Launched.rds
# Notes:
#   - Robust to missing methods/files.
#   - Uses mean(sRGES) to collapse; switch to median if preferred.
#   - Trims names, enforces required columns, and keeps only drugs seen in ≥ k methods.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
})

# ---------------- Config ----------------
base_res   <- "~/CESC_Network/6_OCTAD/6_3_RES"
signatures <- c("A7","A9")
methods    <- c("DESeq2","EdgeR","limma")
k_methods  <- 3  # keep drugs present in ≥3 methods
dir.create(base_res, showWarnings = FALSE, recursive = TRUE)

# -------------- Helpers ---------------

# read_full_rds(sig, mth)
# - Reads RES_<sig>_<mth>_FULL.rds.
# - Validates columns, trims names, tags with signature/method.
# - Returns tibble with (pert_iname, sRGES, signature, method) or NULL if missing/invalid.
read_full_rds <- function(sig, mth) {
  fp <- file.path(base_res, sig, mth, sprintf("RES_%s_%s_FULL.rds", sig, mth))
  if (!file.exists(fp)) return(NULL)
  x <- tryCatch(readRDS(fp), error = function(e) NULL)
  if (is.null(x)) return(NULL)
  x <- as_tibble(x)
  req <- c("pert_iname","sRGES")
  if (!all(req %in% names(x))) return(NULL)
  x %>%
    mutate(
      pert_iname = str_trim(as.character(pert_iname)),
      signature  = sig,
      method     = mth
    ) %>%
    dplyr::select(pert_iname, sRGES, signature, method)
}

# read_fda_launched_csv(sig, mth)
# - Reads sRGES_FDAapproveddrugs.csv.
# - Filters clinical_phase == "Launched".
# - Harmonizes to same 4 columns for downstream union.
read_fda_launched_csv <- function(sig, mth) {
  fp <- file.path(base_res, sig, mth, "sRGES_FDAapproveddrugs.csv")
  if (!file.exists(fp)) return(NULL)
  x <- tryCatch(readr::read_csv(fp, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(x)) return(NULL)
  req <- c("pert_iname","clinical_phase","sRGES")
  if (!all(req %in% names(x))) return(NULL)
  x %>%
    mutate(
      pert_iname     = str_trim(as.character(pert_iname)),
      clinical_phase = as.character(clinical_phase)
    ) %>%
    filter(clinical_phase == "Launched") %>%
    transmute(pert_iname, sRGES, signature = sig, method = mth)
}

# bind_methods_safe(sig, reader_fun)
# - Applies a reader across all methods and binds only non-null frames.
# - Keeps the pipeline resilient to missing outputs.
bind_methods_safe <- function(sig, reader_fun) {
  lst <- lapply(methods, function(m) reader_fun(sig, m))
  lst <- Filter(Negate(is.null), lst)
  if (length(lst) == 0) return(NULL)
  bind_rows(lst)
}

# keep_common_k(df, k)
# - Returns the set of pert_iname present in at least k distinct methods.
# - Avoids “one-off” hits and stabilizes consensus across methods.
keep_common_k <- function(df, k = 3) {
  df %>%
    group_by(pert_iname) %>%
    summarise(n_methods = n_distinct(method), .groups = "drop") %>%
    filter(n_methods >= k) %>%
    pull(pert_iname)
}

# collapse_mean(df_sig_common)
# - Collapses replicate rows per drug using mean(sRGES).
# - Use median if your distribution is heavy-tailed or you want outlier-robust estimates.
collapse_mean <- function(df_sig_common) {
  df_sig_common %>%
    group_by(pert_iname) %>%
    summarise(sRGES = mean(sRGES, na.rm = TRUE), .groups = "drop") %>%
    arrange(sRGES)
}

# save_sig_outputs(sig, tag, df_all)
# - Given a long table (all methods stacked), restricts to drugs in ≥ k_methods,
#   saves the “common” long table and the collapsed table.
# - Tag is either "FULL" or "FDA_Launched" to name outputs clearly.
save_sig_outputs <- function(sig, tag, df_all, out_dir = file.path(base_res, sig)) {
  if (is.null(df_all) || nrow(df_all) == 0) {
    message(sprintf(".. [%s - %s] no data, skipping save.", sig, tag)); return(invisible(NULL))
  }
  keep <- keep_common_k(df_all, k_methods)
  df_common <- df_all %>%
    filter(pert_iname %in% keep) %>%
    arrange(pert_iname, method, sRGES)
  df_collapsed <- collapse_mean(df_common)
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_common    <- file.path(out_dir, sprintf("RES_%s_common_%dmethods_%s.rds", sig, k_methods, tag))
  out_collapsed <- file.path(out_dir, sprintf("RES_%s_common3_collapsed_%s.rds", sig, tag))
  saveRDS(df_common,    out_common)
  saveRDS(df_collapsed, out_collapsed)
  
  message(sprintf(".. [%s - %s] saved: %s (n=%d) & %s (n=%d)",
                  sig, tag, basename(out_common), nrow(df_common),
                  basename(out_collapsed), nrow(df_collapsed)))
}

# -------------- Main -------------------
for (sig in signatures) {
  message(sprintf("\n== Signature %s ==", sig))
  
  # 1) FULL: read all available FULL.rds, then save common+collapsed
  full_all <- bind_methods_safe(sig, read_full_rds)
  if (!is.null(full_all)) {
    full_all <- arrange(full_all, signature, method, sRGES)
    save_sig_outputs(sig, "FULL", full_all)
  } else {
    message(sprintf(".. [%s - FULL] no FULL.rds found in any method.", sig))
  }
  
  # 2) FDA Launched: read per-method CSVs, filter Launched, then save common+collapsed
  fda_all <- bind_methods_safe(sig, read_fda_launched_csv)
  if (!is.null(fda_all)) {
    fda_all <- arrange(fda_all, signature, method, sRGES)
    save_sig_outputs(sig, "FDA_Launched", fda_all)
  } else {
    message(sprintf(".. [%s - FDA_Launched] no sRGES_FDAapproveddrugs.csv found in any method.", sig))
  }
}

message("\n>> DONE. Collapsed outputs per signature saved under: ", base_res)

