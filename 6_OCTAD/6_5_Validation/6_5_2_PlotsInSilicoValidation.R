# =====================
# In Silico Validation Plots
# =====================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
})

# ---------------------
# Base folders (new layout)
# ---------------------
base_dir      <- "~/CESC_Network/6_OCTAD/6_5_Validation"
base_dir_A7   <- file.path(base_dir, "Validation_A7")
base_dir_A9   <- file.path(base_dir, "Validation_A9")
plots_out_dir <- file.path(base_dir, "Plots")
dir.create(plots_out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------
# Robust reader for in-silico TSVs inside CellLineEval/
# ---------------------
read_insilico_file <- function(cell_line, base_clade_dir, type = c("ic50", "auc")) {
  type <- match.arg(type)
  cell_dir <- file.path(base_clade_dir, cell_line, "CellLineEval")
  if (!dir.exists(cell_dir)) return(NULL)
  
  # Expected filename
  expected <- file.path(cell_dir, sprintf("%s_%s_insilico_data.tsv", cell_line, type))
  candidate <- NULL
  if (file.exists(expected)) {
    candidate <- expected
  } else {
    # Fallback: find any file containing _ic50_ or _auc_
    patt <- sprintf("_%s_", type)
    hits <- list.files(cell_dir, pattern = patt, full.names = TRUE)
    if (length(hits) > 0) candidate <- hits[[1]]
  }
  if (is.null(candidate)) return(NULL)
  
  df <- tryCatch(read_tsv(candidate, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(df)) return(NULL)
  
  df$cell_line <- cell_line
  df
}

# ---------------------
# List cell line folders for a clade
# ---------------------
list_cell_lines <- function(base_clade_dir) {
  if (!dir.exists(base_clade_dir)) return(character())
  list.dirs(base_clade_dir, full.names = FALSE, recursive = FALSE)
}

# ---------------------
# Read and prepare clade-specific tables
# ---------------------
prep_clade_data <- function(base_clade_dir) {
  cell_lines <- list_cell_lines(base_clade_dir)
  
  ic50_tbl <- map_dfr(cell_lines, ~ read_insilico_file(.x, base_clade_dir, "ic50"))
  auc_tbl  <- map_dfr(cell_lines, ~ read_insilico_file(.x, base_clade_dir, "auc"))
  
  # Normalizar columnas
  if (!is.null(ic50_tbl) && nrow(ic50_tbl) > 0) {
    if (!"medIC50" %in% names(ic50_tbl)) {
      cand <- intersect(c("medianIC50","MedianIC50","medic50","MEDIC50"), names(ic50_tbl))
      if (length(cand) == 1) ic50_tbl <- ic50_tbl %>% rename(medIC50 = !!cand)
    }
    if (!"StronglyPredicted" %in% names(ic50_tbl)) {
      ic50_tbl <- ic50_tbl %>% mutate(StronglyPredicted = sRGES <= -0.20)
    }
    ic50_tbl <- ic50_tbl %>%
      filter(!is.na(medIC50), is.finite(medIC50), !is.na(sRGES))   
  }
  
  if (!is.null(auc_tbl) && nrow(auc_tbl) > 0) {
    if (!"medauc" %in% names(auc_tbl)) {
      cand <- intersect(c("medianAUC","MedianAUC","medAUC","MEDAUC","AUC","auc"), names(auc_tbl))
      if (length(cand) == 1) auc_tbl <- auc_tbl %>% rename(medauc = !!cand)
    }
    if (!"StronglyPredicted" %in% names(auc_tbl)) {
      auc_tbl <- auc_tbl %>% mutate(StronglyPredicted = sRGES <= -0.20)
    }
    auc_tbl <- auc_tbl %>%
      filter(!is.na(medauc), is.finite(medauc), !is.na(sRGES))   # <- ya NO filtra por sRGES < -0.20
  }
  
  list(ic50 = ic50_tbl, auc = auc_tbl)
}


data_A7 <- prep_clade_data(base_dir_A7)
data_A9 <- prep_clade_data(base_dir_A9)

# ---------------------
# Plotters
# ---------------------
# ========= Helpers para estadísticas =========
# ===== stats helpers =====
stats_ic50_by_panel <- function(tbl) {
  if (is.null(tbl) || nrow(tbl) == 0) return(NULL)
  tbl %>%
    dplyr::group_by(cell_line) %>%
    dplyr::summarise(
      n    = dplyr::n(),
      r    = suppressWarnings(cor(sRGES, log10(medIC50), method = "pearson")),
      rho  = suppressWarnings(cor(sRGES, log10(medIC50), method = "spearman")),
      beta = tryCatch(coef(lm(log10(medIC50) ~ sRGES))[2], error = function(e) NA_real_),
      R2   = tryCatch(summary(lm(log10(medIC50) ~ sRGES))$r.squared, error = function(e) NA_real_),
      .groups = "drop"
    ) %>% dplyr::mutate(label = sprintf("r=%.2f, ρ=%.2f, β=%.2f, R²=%.2f (n=%d)", r, rho, beta, R2, n))
}

stats_auc_by_panel <- function(tbl) {
  if (is.null(tbl) || nrow(tbl) == 0) return(NULL)
  tbl %>%
    dplyr::group_by(cell_line) %>%
    dplyr::summarise(
      n    = dplyr::n(),
      r    = suppressWarnings(cor(sRGES, medauc, method = "pearson")),
      rho  = suppressWarnings(cor(sRGES, medauc, method = "spearman")),
      beta = tryCatch(coef(lm(medauc ~ sRGES))[2], error = function(e) NA_real_),
      R2   = tryCatch(summary(lm(medauc ~ sRGES))$r.squared, error = function(e) NA_real_),
      .groups = "drop"
    ) %>% dplyr::mutate(label = sprintf("r=%.2f, ρ=%.2f, β=%.2f, R²=%.2f (n=%d)", r, rho, beta, R2, n))
}

stats_ic50_global <- function(tbl) {
  if (is.null(tbl) || nrow(tbl) == 0) return(NULL)
  fit <- lm(log10(medIC50) ~ sRGES, data = tbl)
  r   <- suppressWarnings(cor(tbl$sRGES, log10(tbl$medIC50), method = "pearson"))
  rho <- suppressWarnings(cor(tbl$sRGES, log10(tbl$medIC50), method = "spearman"))
  data.frame(label = sprintf("r=%.2f, ρ=%.2f, β=%.2f, R²=%.2f (n=%d)",
                             r, rho, coef(fit)[2], summary(fit)$r.squared, nrow(tbl)))
}

stats_auc_global <- function(tbl) {
  if (is.null(tbl) || nrow(tbl) == 0) return(NULL)
  fit <- lm(medauc ~ sRGES, data = tbl)
  r   <- suppressWarnings(cor(tbl$sRGES, tbl$medauc, method = "pearson"))
  rho <- suppressWarnings(cor(tbl$sRGES, tbl$medauc, method = "spearman"))
  data.frame(label = sprintf("r=%.2f, ρ=%.2f, β=%.2f, R²=%.2f (n=%d)",
                             r, rho, coef(fit)[2], summary(fit)$r.squared, nrow(tbl)))
}

# ===== en los globales, 1 sola recta =====
plot_global_ic50 <- function(tbl, title_prefix) {
  if (is.null(tbl) || nrow(tbl) == 0) return(NULL)
  ann <- stats_ic50_global(tbl)
  ggplot(tbl, aes(x = sRGES, y = log10(medIC50), color = StronglyPredicted)) +
    geom_point(alpha = 0.7) +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black", fill = "grey80") +
    annotate("text", x = -Inf, y = Inf, label = ann$label,
             hjust = -0.05, vjust = 1.1, size = 4) +
    labs(title = sprintf("%s: Global Correlation sRGES vs log10(IC50)", title_prefix),
         y = "log10(medIC50)", x = "sRGES") +
    theme_minimal(base_size = 14)
}

plot_global_auc <- function(tbl, title_prefix) {
  if (is.null(tbl) || nrow(tbl) == 0) return(NULL)
  ann <- stats_auc_global(tbl)
  ggplot(tbl, aes(x = sRGES, y = medauc, color = StronglyPredicted)) +
    geom_point(alpha = 0.7) +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black", fill = "grey80") +
    annotate("text", x = -Inf, y = Inf, label = ann$label,
             hjust = -0.05, vjust = 1.1, size = 4) +
    labs(title = sprintf("%s: Global Correlation sRGES vs AUC", title_prefix),
         y = "AUC", x = "sRGES") +
    theme_minimal(base_size = 14)
}


# ---------------------
# Build plots
# ---------------------
p_ic50_A7       <- plot_facets_ic50(data_A7$ic50, "A7")
p_auc_A7        <- plot_facets_auc (data_A7$auc,  "A7")
p_ic50_A9       <- plot_facets_ic50(data_A9$ic50, "A9")
p_auc_A9        <- plot_facets_auc (data_A9$auc,  "A9")

p_ic50_A7_glob  <- plot_global_ic50(data_A7$ic50, "A7")
p_auc_A7_glob   <- plot_global_auc (data_A7$auc,  "A7")
p_ic50_A9_glob  <- plot_global_ic50(data_A9$ic50, "A9")
p_auc_A9_glob   <- plot_global_auc (data_A9$auc,  "A9")

# ---------------------
# Save at maximum resolution (vector PDF + 600 dpi PNG)
# (No explicit width/height to respect your preference)
# ---------------------
save_plot_all <- function(plot_obj, stem) {
  if (is.null(plot_obj)) return(invisible(NULL))
  pdf_out  <- file.path(plots_out_dir, paste0(stem, ".pdf"))
  png_out  <- file.path(plots_out_dir, paste0(stem, ".png"))
  # Vector (best for journals)
  ggsave(filename = pdf_out, plot = plot_obj, device = cairo_pdf)
  # High-DPI raster
  ggsave(filename = png_out, plot = plot_obj, dpi = 1200)
  invisible(list(pdf = pdf_out, png = png_out))
}





save_plot_all(p_ic50_A7,      "A7_facet_ic50")
save_plot_all(p_auc_A7,       "A7_facet_auc")
save_plot_all(p_ic50_A9,      "A9_facet_ic50")
save_plot_all(p_auc_A9,       "A9_facet_auc")

save_plot_all(p_ic50_A7_glob, "A7_global_ic50")
save_plot_all(p_auc_A7_glob,  "A7_global_auc")
save_plot_all(p_ic50_A9_glob, "A9_global_ic50")
save_plot_all(p_auc_A9_glob,  "A9_global_auc")

# ---------------------
# Print to device (optional interactive visualization)
# ---------------------
if (!is.null(p_ic50_A7))      print(p_ic50_A7)
if (!is.null(p_auc_A7))       print(p_auc_A7)
if (!is.null(p_ic50_A9))      print(p_ic50_A9)
if (!is.null(p_auc_A9))       print(p_auc_A9)

if (!is.null(p_ic50_A7_glob)) print(p_ic50_A7_glob)
if (!is.null(p_auc_A7_glob))  print(p_auc_A7_glob)
if (!is.null(p_ic50_A9_glob)) print(p_ic50_A9_glob)
if (!is.null(p_auc_A9_glob))  print(p_auc_A9_glob)
