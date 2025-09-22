# =====================
# In Silico Validation Plots (read-only; exact same outputs + global Pearson on collapsed drugs)
# =====================

suppressPackageStartupMessages({
  library(readr);    library(dplyr);   library(purrr);  library(tidyr)
  library(ggplot2);  library(ggrepel); library(xml2);   library(jsonlite)
  library(ragg)
  library(patchwork)
})

# ---- Base folders ----
base_dir      <- "~/CESC_Network/6_OCTAD/6_5_Validation"
base_dir_A7   <- file.path(base_dir, "Validation_A7")
base_dir_A9   <- file.path(base_dir, "Validation_A9")
plots_out_dir <- file.path(base_dir, "Plots"); dir.create(plots_out_dir, TRUE, FALSE)

# ---- Constants: exact anchors in the TXT logs ----
ANCHORS <- list(IC50 = "IC50 cortest", AUC = "AUC cortest")

# ---- Small utilities ----
list_cell_lines <- function(base_clade_dir){
  if (!dir.exists(base_clade_dir)) return(character())
  list.dirs(base_clade_dir, full.names = FALSE, recursive = FALSE)
}

# Read OCTAD TSV exactly named "<cellline>_<type>_insilico_data.tsv"
# (points only; we do not recompute anything from panels)
read_insilico_file <- function(cell_line, base_clade_dir, type = c("ic50","auc")){
  type <- match.arg(type)
  f <- file.path(base_clade_dir, cell_line, "CellLineEval",
                 sprintf("%s_%s_insilico_data.tsv", cell_line, type))
  if (!file.exists(f)) return(NULL)
  df <- tryCatch(readr::read_tsv(f, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(df)) return(NULL)
  df$cell_line <- cell_line
  df
}

# ---- Parse r and p from TXT (preferred) ----
# Fixed format:
#   line: "IC50 cortest" or "AUC cortest"
#   nearby: "p-value = ..."
#   then: "sample estimates:" with numeric "cor"
.parse_cor_block <- function(lines, start_idx){
  rng <- seq.int(start_idx, min(length(lines), start_idx + 60))
  t_rel <- grep("p[- ]?value\\s*=", lines[rng], ignore.case = TRUE)
  if (!length(t_rel)) return(NULL)
  tline <- lines[rng[t_rel[1]]]
  p_raw <- sub(".*p[- ]?value\\s*=\\s*([^,\\n]+).*", "\\1", tline, ignore.case = TRUE)
  p_num <- suppressWarnings(as.numeric(gsub("^\\s*[<>]=?\\s*", "", p_raw)))
  est_rel <- grep("sample estimates", lines[rng], ignore.case = TRUE)
  r_val <- NA_real_
  if (length(est_rel)){
    look <- seq.int(rng[est_rel[1]] + 1, min(length(lines), rng[est_rel[1]] + 5))
    for (kk in look){
      got <- regmatches(lines[kk], regexpr("[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?|NaN|NA", lines[kk], perl = TRUE))
      if (length(got) && nzchar(got)){
        r_val <- suppressWarnings(as.numeric(gsub("(NaN|NA)", NA_character_, got)))
        break
      }
    }
  }
  tibble::tibble(r = r_val, p_value_raw = trimws(p_raw), p_value_num = p_num)
}

parse_octad_results_txt <- function(file_path, cell_line){
  if (!file.exists(file_path)) return(NULL)
  ln <- readLines(file_path, warn = FALSE)
  eq_trim <- function(x, val) which(trimws(x) == val)
  i_ic50 <- eq_trim(ln, ANCHORS$IC50)
  i_auc  <- eq_trim(ln, ANCHORS$AUC)
  out <- list()
  if (length(i_ic50)) out <- append(out, list(.parse_cor_block(ln, i_ic50[1]) %>% mutate(cell_line = cell_line, metric = "IC50")))
  if (length(i_auc))  out <- append(out, list(.parse_cor_block(ln, i_auc[1])  %>% mutate(cell_line = cell_line, metric = "AUC")))
  if (!length(out)) return(NULL)
  bind_rows(out)
}

# ---- (Facets) Extract already-drawn smooth line from the HTML (no recompute) ----
extract_smooth_from_html <- function(cell_line, base_clade_dir, metric = c("ic50","auc")){
  metric <- match.arg(metric)
  html_path <- file.path(base_clade_dir, cell_line, "CellLineEval",
                         sprintf("%s_%s_insilico_validation.html", cell_line, metric))
  if (!file.exists(html_path)) return(NULL)
  doc <- tryCatch(xml2::read_html(html_path), error = function(e) NULL); if (is.null(doc)) return(NULL)
  nodes <- xml2::xml_find_all(doc, '//script[@type="application/json"]'); if (!length(nodes)) return(NULL)
  collect_line_traces <- function(obj, acc){
    if (is.list(obj)){
      mode <- if (!is.null(obj$mode)) as.character(obj$mode) else NA_character_
      if (!is.null(obj$x) && !is.null(obj$y)){
        x <- suppressWarnings(as.numeric(unlist(obj$x)))
        y <- suppressWarnings(as.numeric(unlist(obj$y)))
        if (length(x) >= 2 && length(y) >= 2 && isTRUE(grepl("line", mode, TRUE)))
          acc[[length(acc)+1]] <- tibble::tibble(x = x, y = y)
      }
      for (el in obj) acc <- collect_line_traces(el, acc)
    }
    acc
  }
  for (i in seq_along(nodes)){
    js <- tryCatch(jsonlite::fromJSON(xml2::xml_text(nodes[[i]]), simplifyVector = FALSE), error = function(e) NULL)
    if (is.null(js)) next
    acc <- collect_line_traces(js, list())
    if (length(acc)) return(bind_rows(acc) %>% arrange(x, y) %>% mutate(cell_line = cell_line, metric = metric))
  }
  NULL
}

collect_smooth_lines <- function(base_clade_dir, which_metric = c("ic50","auc")){
  which_metric <- match.arg(which_metric)
  cls <- list_cell_lines(base_clade_dir); if (!length(cls)) return(tibble())
  map_dfr(cls, ~ extract_smooth_from_html(.x, base_clade_dir, which_metric))
}

# ---- Load points (TSV) & normalize minimal columns ----
prep_clade_data <- function(base_clade_dir){
  cls <- list_cell_lines(base_clade_dir)
  ic50_tbl <- map_dfr(cls, ~ read_insilico_file(.x, base_clade_dir, "ic50"))
  auc_tbl  <- map_dfr(cls, ~ read_insilico_file(.x, base_clade_dir, "auc"))
  if (nrow(ic50_tbl)){
    if (!"medIC50" %in% names(ic50_tbl)) {
      cand <- intersect(c("medianIC50","MedianIC50","medic50","MEDIC50"), names(ic50_tbl))
      if (length(cand) == 1) ic50_tbl <- ic50_tbl %>% rename(medIC50 = !!cand)
    }
    if (!"StronglyPredicted" %in% names(ic50_tbl)) ic50_tbl <- ic50_tbl %>% mutate(StronglyPredicted = sRGES <= -0.20)
    ic50_tbl <- ic50_tbl %>% filter(is.finite(medIC50), !is.na(sRGES))
  }
  if (nrow(auc_tbl)){
    if (!"medauc" %in% names(auc_tbl)) {
      cand <- intersect(c("medianAUC","MedianAUC","medAUC","MEDAUC","AUC","auc"), names(auc_tbl))
      if (length(cand) == 1) auc_tbl <- auc_tbl %>% rename(medauc = !!cand)
    }
    if (!"StronglyPredicted" %in% names(auc_tbl)) auc_tbl <- auc_tbl %>% mutate(StronglyPredicted = sRGES <= -0.20)
    auc_tbl <- auc_tbl %>% filter(is.finite(medauc), !is.na(sRGES))
  }
  list(ic50 = ic50_tbl, auc = auc_tbl, cells = cls)
}

# ---- Collect r,p labels per cell line (from TXT) ----
collect_octad_metrics <- function(base_clade_dir){
  cls <- list_cell_lines(base_clade_dir); if (!length(cls)) return(tibble())
  out <- list()
  for (cl in cls){
    f_txt <- file.path(base_clade_dir, cl, "CellLineEval", sprintf("%s_drug_sensitivity_insilico_results.txt", cl))
    if (file.exists(f_txt)) {
      got <- parse_octad_results_txt(f_txt, cl)
      if (!is.null(got)) out <- append(out, list(got))
    }
  }
  if (!length(out)) return(tibble())
  bind_rows(out) %>%
    mutate(label = paste0(
      "r = ", if_else(is.na(r), "NA", sprintf("%.2f", r)),
      "; p = ",
      if_else(is.na(p_value_raw),
              if_else(is.na(p_value_num), "NA", format.pval(p_value_num, digits = 2, eps = 1e-16)),
              p_value_raw)))
}

# --- FDA highlight lists (Launched) ---
RGEs_A7 <- readRDS("~/CESC_Network/6_OCTAD/6_3_RES/A7/RES_A7_common3_collapsed_FDA_Launched.rds")
RGEs_A9 <- readRDS("~/CESC_Network/6_OCTAD/6_3_RES/A9/RES_A9_common3_collapsed_FDA_Launched.rds")
fda_set_A7 <- tolower(trimws(RGEs_A7$pert_iname))
fda_set_A9 <- tolower(trimws(RGEs_A9$pert_iname))

add_fda_flag <- function(tbl, fda_set) {
  if (is.null(tbl) || !nrow(tbl)) return(tbl)
  tbl %>% mutate(iname_l = tolower(trimws(pert_iname)),
                 fda_flag = ifelse(iname_l %in% fda_set, "FDA", "Other"))
}

# =====================
# Collapse helpers (1 point per drug) + global Pearson on collapsed data
# =====================
collapse_ic50_by_drug <- function(tbl){
  tbl %>%
    filter(is.finite(sRGES), is.finite(medIC50), !is.na(pert_iname)) %>%
    group_by(pert_iname) %>%
    summarise(
      sRGES  = median(sRGES, na.rm = TRUE),
      y      = median(log10(medIC50), na.rm = TRUE),
      fda_flag = ifelse(any(fda_flag == "FDA", na.rm = TRUE), "FDA", "Other"),
      StronglyPredicted = ifelse(median(sRGES, na.rm = TRUE) < -0.20, "Yes", "No"),
      .groups = "drop"
    )
}

collapse_auc_by_drug <- function(tbl){
  tbl %>%
    filter(is.finite(sRGES), is.finite(medauc), !is.na(pert_iname)) %>%
    group_by(pert_iname) %>%
    summarise(
      sRGES  = median(sRGES, na.rm = TRUE),
      y      = median(medauc,     na.rm = TRUE),
      fda_flag = ifelse(any(fda_flag == "FDA", na.rm = TRUE), "FDA", "Other"),
      StronglyPredicted = ifelse(median(sRGES, na.rm = TRUE) < -0.20, "Yes", "No"),
      .groups = "drop"
    )
}

pearson_label <- function(x, y){
  ct <- suppressWarnings(cor.test(x, y, method = "pearson"))
  paste0("r = ", sprintf("%.2f", as.numeric(ct$estimate)),
         "; p = ", format.pval(ct$p.value, digits = 2, eps = 1e-16))
}

# ---- Load everything ----
data_A7 <- prep_clade_data(base_dir_A7)
data_A9 <- prep_clade_data(base_dir_A9)
metrics_A7 <- collect_octad_metrics(base_dir_A7)
metrics_A9 <- collect_octad_metrics(base_dir_A9)
smooth_A7_ic50 <- collect_smooth_lines(base_dir_A7, "ic50")
smooth_A7_auc  <- collect_smooth_lines(base_dir_A7, "auc")
smooth_A9_ic50 <- collect_smooth_lines(base_dir_A9, "ic50")
smooth_A9_auc  <- collect_smooth_lines(base_dir_A9, "auc")

# Add FDA flag to raw points (used in facets and as input to the collapse)
data_A7$ic50 <- add_fda_flag(data_A7$ic50, fda_set_A7)
data_A7$auc  <- add_fda_flag(data_A7$auc,  fda_set_A7)
data_A9$ic50 <- add_fda_flag(data_A9$ic50, fda_set_A9)
data_A9$auc  <- add_fda_flag(data_A9$auc,  fda_set_A9)

# ---- Plotters (facets) ----
plot_facets_ic50 <- function(tbl, title_prefix, metrics_tbl, smooth_df, fda_set){
  if (!nrow(tbl)) return(NULL)
  ann <- metrics_tbl %>% filter(metric=="IC50") %>% distinct(cell_line, label)
  
  labs_df <- tbl %>%
    mutate(iname_l = tolower(trimws(pert_iname))) %>%
    filter(iname_l %in% fda_set, StronglyPredicted == "Yes",
           is.finite(sRGES), is.finite(medIC50)) %>%
    group_by(cell_line, pert_iname) %>%
    summarise(sRGES = median(sRGES, na.rm = TRUE),
              y     = median(log10(medIC50), na.rm = TRUE),
              .groups = "drop")
  
  # DEFAULT: silence FDA drug names in facets
  labs_df <- labs_df[0, ]
  
  ggplot(tbl, aes(x = sRGES, y = log10(medIC50), color = StronglyPredicted)) +
    geom_point() +
    { if (nrow(labs_df))
      ggrepel::geom_text_repel(
        data = labs_df, inherit.aes = FALSE,
        aes(x = sRGES, y = y, label = pert_iname),
        size = 2.3, segment.color = "grey70", max.overlaps = 5,
        box.padding = 0.5, point.padding = 0.3, direction = "y")
    } +
    { if (nrow(smooth_df))
      geom_path(data = smooth_df, aes(x = x, y = y, group = cell_line),
                inherit.aes = FALSE, color = "black", linewidth = 0.6)
    } +
    facet_wrap(~cell_line, scales = "free") +
    geom_text(data = ann, inherit.aes = FALSE, aes(x = -Inf, y = Inf, label = label),
              hjust = -0.05, vjust = 1.1, size = 3.1) +
    labs(title = sprintf("%s sRGES by cell line", title_prefix),
         y = "log10(medIC50)", x = "sRGES") +
    coord_cartesian(ylim = c(NA, 15)) +
    theme_bw() +
    theme(
      plot.title  = element_text(size = 20, face = "bold"),  # bigger overall title
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16)
    )
}

# ---- Plotters (facets) ----
plot_facets_auc <- function(tbl, title_prefix, metrics_tbl, smooth_df, fda_set){
  if (!nrow(tbl)) return(NULL)
  ann <- metrics_tbl %>% filter(metric=="AUC") %>% distinct(cell_line, label)
  
  labs_df <- tbl %>%
    mutate(iname_l = tolower(trimws(pert_iname))) %>%
    filter(iname_l %in% fda_set, is.finite(sRGES), is.finite(medauc)) %>%
    group_by(cell_line, pert_iname) %>%
    summarise(sRGES = median(sRGES, na.rm = TRUE),
              y     = median(medauc, na.rm = TRUE),
              .groups = "drop")
  
  # DEFAULT: silence FDA drug names in facets
  labs_df <- labs_df[0, ]
  
  ggplot(tbl, aes(x = sRGES, y = medauc, color = StronglyPredicted)) +
    geom_point() +
    { if (nrow(labs_df))
      ggrepel::geom_text_repel(
        data = labs_df, inherit.aes = FALSE,
        aes(x = sRGES, y = y, label = pert_iname),
        size = 2.3, segment.color = "grey70", max.overlaps = 5,
        box.padding = 0.5, point.padding = 0.3, direction = "y")
    } +
    { if (nrow(smooth_df))
      geom_path(data = smooth_df, aes(x = x, y = y, group = cell_line),
                inherit.aes = FALSE, color = "black", linewidth = 0.6)
    } +
    facet_wrap(~cell_line, scales = "free") +
    geom_text(data = ann, inherit.aes = FALSE, aes(x = -Inf, y = Inf, label = label),
              hjust = -0.05, vjust = 1.1, size = 3.1) +
    labs(title = sprintf("%s sRGES by cell line", title_prefix),
         y = "AUC (recomputed)",  # <— updated
         x = "sRGES") +
    theme_bw() +
    theme(
      plot.title  = element_text(size = 20, face = "bold"),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16)
    )
}


# ---- Plotters (globals with collapsed data + Pearson of collapsed points + LM line) ----
plot_global_ic50 <- function(tbl, title_prefix){
  if (is.null(tbl) || !nrow(tbl)) return(NULL)
  
  dd  <- collapse_ic50_by_drug(tbl)          # 1 point per drug
  lab <- pearson_label(dd$sRGES, dd$y)       # Pearson on collapsed points
  
  ggplot(dd, aes(x = sRGES, y = y)) +
    geom_point(aes(color = StronglyPredicted, shape = fda_flag),
               alpha = 0.9, size = 1.9) +
    # Linear fit over collapsed points (same x,y used in Pearson)
    geom_smooth(aes(group = 1), method = "lm", se = TRUE,
                color = "black", fill = "grey80") +
    annotate("text", x = -Inf, y = Inf, label = lab,
             hjust = -0.05, vjust = 1.1, size = 3.1) +
    scale_shape_manual(values = c(FDA = 17, Other = 16),
                       name   = "FDA status",
                       labels = c(FDA = "Launched", Other = "Other")) +
    labs(title = sprintf("%s", title_prefix),
         y = "log10(medIC50)", x = "sRGES") +
    coord_cartesian(ylim = range(dd$y, na.rm = TRUE),
                    xlim = range(dd$sRGES, na.rm = TRUE),
                    clip = "off") +
    theme_minimal(base_size = 16)
}

# ---- Plotters (globals with collapsed data + Pearson of collapsed points + LM line) ----
plot_global_auc <- function(tbl, title_prefix){
  if (is.null(tbl) || !nrow(tbl)) return(NULL)
  
  dd  <- collapse_auc_by_drug(tbl)           # 1 point per drug
  lab <- pearson_label(dd$sRGES, dd$y)       # Pearson on collapsed points
  
  ggplot(dd, aes(x = sRGES, y = y)) +
    geom_point(aes(color = StronglyPredicted, shape = fda_flag),
               alpha = 0.9, size = 1.9) +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE,
                color = "black", fill = "grey80") +
    annotate("text", x = -Inf, y = Inf, label = lab,
             hjust = -0.05, vjust = 1.1, size = 3.1) +
    scale_shape_manual(values = c(FDA = 17, Other = 16),
                       name   = "FDA status",
                       labels = c(FDA = "Launched", Other = "Other")) +
    labs(title = sprintf("%s", title_prefix),
         y = "AUC (recomputed)",  # <— updated
         x = "sRGES") +
    coord_cartesian(ylim = range(dd$y, na.rm = TRUE),
                    xlim = range(dd$sRGES, na.rm = TRUE),
                    clip = "off") +
    theme_minimal(base_size = 16)
}



# ---- Build plots ----
p_ic50_A7 <- plot_facets_ic50(data_A7$ic50, "HPV-A7", metrics_A7, smooth_A7_ic50, fda_set_A7)
p_auc_A7  <- plot_facets_auc (data_A7$auc,  "HPV-A7", metrics_A7, smooth_A7_auc,  fda_set_A7)
p_ic50_A9 <- plot_facets_ic50(data_A9$ic50, "HPV-A9", metrics_A9, smooth_A9_ic50, fda_set_A9)
p_auc_A9  <- plot_facets_auc (data_A9$auc,  "HPV-A9", metrics_A9, smooth_A9_auc,  fda_set_A9)

# Globals (previously suggested titles)
p_ic50_A7_glob  <- plot_global_ic50(data_A7$ic50, "HPV-A7 — Drug response association across evaluated cell lines")
p_auc_A7_glob   <- plot_global_auc (data_A7$auc,  "HPV-A7")
p_ic50_A9_glob  <- plot_global_ic50(data_A9$ic50, "HPV-A9 — Drug response association across evaluated cell lines")
p_auc_A9_glob   <- plot_global_auc (data_A9$auc,  "HPV-A9")

# ---- Save (PDF + 1200 dpi PNG) ----
save_plot <- function(p, name, w = 12, h = 7, dpi = 600, dir = plots_out_dir) {
  if (is.null(p)) return(invisible(NULL))
  ggsave(file.path(dir, paste0(name, ".png")), p,
         device = ragg::agg_png, width = w, height = h, units = "in", dpi = dpi, limitsize = FALSE)
  ggsave(file.path(dir, paste0(name, ".pdf")), p,
         device = "pdf", width = w, height = h, units = "in", useDingbats = FALSE)
}

# Facets at 1200 dpi
save_plot(p_ic50_A7, "A7_facet_ic50", w = 12, h = 7, dpi = 1200)
save_plot(p_auc_A7,  "A7_facet_auc",  w = 12, h = 7, dpi = 1200)
save_plot(p_ic50_A9, "A9_facet_ic50", w = 12, h = 7, dpi = 1200)
save_plot(p_auc_A9,  "A9_facet_auc",  w = 12, h = 7, dpi = 1200)

# Globals at 1200 dpi
save_plot(p_ic50_A7_glob, "A7_global_ic50", w = 9, h = 6, dpi = 1200)
save_plot(p_auc_A7_glob,  "A7_global_auc",  w = 9, h = 6, dpi = 1200)
save_plot(p_ic50_A9_glob, "A9_global_ic50", w = 9, h = 6, dpi = 1200)
save_plot(p_auc_A9_glob,  "A9_global_auc",  w = 9, h = 6, dpi = 1200)

# (Optional interactive display)
for (p in list(p_ic50_A7,p_auc_A7,p_ic50_A9,p_auc_A9,p_ic50_A7_glob,p_auc_A7_glob,p_ic50_A9_glob,p_auc_A9_glob))
  if (!is.null(p)) print(p)

# ==========================
# Grid of the 4 global plots (shared legend; one row title per clade)
# ==========================
p_ic50_A7_glob_n <- p_ic50_A7_glob
p_auc_A7_glob_n  <- p_auc_A7_glob  + theme(plot.title = element_blank())
p_ic50_A9_glob_n <- p_ic50_A9_glob
p_auc_A9_glob_n  <- p_auc_A9_glob  + theme(plot.title = element_blank())

row_A7 <- (p_ic50_A7_glob_n + p_auc_A7_glob_n) + plot_annotation(title = "HPV-A7")
row_A9 <- (p_ic50_A9_glob_n + p_auc_A9_glob_n) + plot_annotation(title = "HPV-A9")

grid_globals <- (row_A7 / row_A9) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
grid_globals

save_plot(grid_globals, "GLOBAL_grid_HPVA7A9", w = 12, h = 10, dpi = 1200)
