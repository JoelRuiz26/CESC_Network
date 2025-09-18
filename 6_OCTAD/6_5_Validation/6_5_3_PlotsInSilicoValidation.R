# =====================
# In Silico Validation Plots (read-only; no recomputation for panels)
# =====================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(purrr); library(stringr)
  library(tidyr); library(ggplot2); library(ggrepel)
  library(jsonlite); library(xml2)
  library(patchwork)   # <-- NUEVO para armar el grid
})


# ---- Base folders ----
base_dir      <- "~/CESC_Network/6_OCTAD/6_5_Validation"
base_dir_A7   <- file.path(base_dir, "Validation_A7")
base_dir_A9   <- file.path(base_dir, "Validation_A9")
plots_out_dir <- file.path(base_dir, "Plots"); dir.create(plots_out_dir, TRUE, FALSE)

# ---- Utilities ----
list_cell_lines <- function(base_clade_dir){
  if (!dir.exists(base_clade_dir)) return(character())
  list.dirs(base_clade_dir, full.names = FALSE, recursive = FALSE)
}

# Read OCTAD-exported TSV (points only; no transforms)
read_insilico_file <- function(cell_line, base_clade_dir, type = c("ic50","auc")){
  type <- match.arg(type)
  cell_dir <- file.path(base_clade_dir, cell_line, "CellLineEval")
  if (!dir.exists(cell_dir)) return(NULL)
  expected  <- file.path(cell_dir, sprintf("%s_%s_insilico_data.tsv", cell_line, type))
  candidate <- if (file.exists(expected)) expected else {
    hits <- list.files(cell_dir, pattern = sprintf("_%s_", type), full.names = TRUE)
    if (length(hits)) hits[[1]] else NULL
  }
  if (is.null(candidate)) return(NULL)
  df <- tryCatch(readr::read_tsv(candidate, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(df)) return(NULL)
  df$cell_line <- cell_line
  df
}

# ---- r,p from OCTAD logs (TXT preferred; HTML fallback) ----
.parse_cor_block <- function(lines, start_idx){
  rng <- seq.int(start_idx, min(length(lines), start_idx + 60))
  t_rel <- grep("t\\s*=.*df\\s*=.*p-value", lines[rng]); if (!length(t_rel)) return(NULL)
  tline <- lines[rng[t_rel[1]]]
  p_raw <- sub(".*p-value\\s*=\\s*([^,\\n]+).*", "\\1", tline)
  p_num <- suppressWarnings(as.numeric(gsub("^\\s*[<>]=?\\s*", "", p_raw)))
  est_rel <- grep("sample estimates", lines[rng]); r_val <- NA_real_
  if (length(est_rel)){
    for (kk in (rng[est_rel[1]]+1):min(length(lines), rng[est_rel[1]]+5)){
      got <- regmatches(lines[kk], regexpr("[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?|NaN|NA", lines[kk], perl=TRUE))
      if (length(got) && nzchar(got)){
        r_val <- suppressWarnings(as.numeric(gsub("(NaN|NA)", NA_character_, got))); break
      }
    }
  }
  tibble::tibble(r = r_val, p_value_raw = trimws(p_raw), p_value_num = p_num)
}

parse_octad_results_txt <- function(file_path, cell_line){
  if (!file.exists(file_path)) return(NULL)
  ln <- readLines(file_path, warn = FALSE)
  ia <- grep("^AUC .*cortest|^AUC cortest|^AUC.*cor\\.test",  ln, ignore.case = TRUE)
  ii <- grep("^IC50 .*cortest|^IC50 cortest|^IC50.*cor\\.test", ln, ignore.case = TRUE)
  out <- list()
  if (length(ii)) out <- append(out, list(.parse_cor_block(ln, ii[1]) %>% mutate(cell_line = cell_line, metric = "IC50")))
  if (length(ia)) out <- append(out, list(.parse_cor_block(ln, ia[1]) %>% mutate(cell_line = cell_line, metric = "AUC")))
  if (!length(out)) return(NULL)
  bind_rows(out)
}

parse_insilico_html_metrics <- function(html_path, cell_line, metric_guess = c("AUC","IC50")){
  metric_guess <- match.arg(metric_guess); if (!file.exists(html_path)) return(NULL)
  ln <- readLines(html_path, warn = FALSE)
  idx <- grep("(Pearson|cor\\.test|correlation|p-value)", ln, ignore.case = TRUE)
  if (!length(idx)) return(NULL)
  seg <- paste(ln[seq.int(idx[1], min(length(ln), idx[1] + 60))], collapse = " ")
  r_val <- {m <- regexpr("r\\s*=\\s*([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)", seg, perl=TRUE, ignore.case=TRUE)
  if (m[1] != -1) as.numeric(sub("r\\s*=\\s*", "", regmatches(seg, m)[1], ignore.case=TRUE)) else NA_real_}
  m_p <- regexpr("p[- ]?value\\s*=\\s*([<>]=?\\s*)?([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)", seg, perl=TRUE, ignore.case=TRUE)
  p_raw <- if (m_p[1] != -1) trimws(sub(".*p[- ]?value\\s*=\\s*", "", regmatches(seg, m_p)[1], ignore.case=TRUE)) else NA_character_
  p_num <- suppressWarnings(as.numeric(gsub("^\\s*[<>]=?\\s*", "", p_raw)))
  if (is.na(r_val) && is.na(p_raw)) return(NULL)
  tibble::tibble(cell_line = cell_line, metric = metric_guess, r = r_val, p_value_raw = p_raw, p_value_num = p_num)
}

collect_octad_metrics <- function(base_clade_dir){
  cls <- list_cell_lines(base_clade_dir); if (!length(cls)) return(tibble())
  out <- list()
  for (cl in cls){
    root <- file.path(base_clade_dir, cl, "CellLineEval"); got <- NULL
    cand_txt <- if (dir.exists(root)) list.files(root, pattern = "_drug_sensitivity_insilico_results\\.txt$", full.names = TRUE, recursive = TRUE) else character()
    if (length(cand_txt)){
      mt <- file.info(cand_txt)$mtime
      got <- parse_octad_results_txt(cand_txt[order(mt, decreasing = TRUE)][1], cl)
    }
    if (is.null(got) || any(!c("IC50","AUC") %in% got$metric)){
      cand_html <- if (dir.exists(root)) list.files(root, pattern = "_(auc|ic50)_insilico_validation\\.html$", full.names = TRUE, ignore.case = TRUE) else character()
      for (hp in cand_html){
        metric_guess <- if (grepl("_auc_", basename(hp), ignore.case = TRUE)) "AUC" else "IC50"
        hdf <- parse_insilico_html_metrics(hp, cl, metric_guess)
        if (!is.null(hdf)) got <- if (is.null(got)) hdf else bind_rows(got, hdf)
      }
    }
    if (!is.null(got)) out <- append(out, list(got))
  }
  if (!length(out)) return(tibble())
  bind_rows(out) %>%
    mutate(label = paste0("r = ", if_else(is.na(r), "NA", sprintf("%.2f", r)),
                          "; p = ", if_else(is.na(p_value_raw),
                                            if_else(is.na(p_value_num), "NA",
                                                    format.pval(p_value_num, digits = 2, eps = 1e-16)),
                                            p_value_raw)))
}

# ---- Extract already-drawn smooth from HTML (facets) ----
extract_smooth_from_html <- function(cell_line, base_clade_dir, metric = c("ic50","auc")){
  metric <- match.arg(metric)
  html_path <- file.path(base_clade_dir, cell_line, "CellLineEval",
                         sprintf("%s_%s_insilico_validation.html", cell_line, metric))
  if (!file.exists(html_path)) return(NULL)
  doc <- tryCatch(xml2::read_html(html_path), error = function(e) NULL); if (is.null(doc)) return(NULL)
  nodes <- xml2::xml_find_all(doc, '//script[@type="application/json"]'); if (!length(nodes)) return(NULL)
  collect_line_traces <- function(obj, acc){
    if (is.list(obj)){
      mode <- tryCatch({ if (!is.null(obj$mode)) as.character(obj$mode) else NA_character_ }, error = function(e) NA_character_)
      if (!is.null(obj$x) && !is.null(obj$y)){
        x <- suppressWarnings(as.numeric(unlist(obj$x))); y <- suppressWarnings(as.numeric(unlist(obj$y)))
        if (length(x) >= 2 && length(y) >= 2 && any(grepl("line", mode %||% "", TRUE)))
          acc[[length(acc)+1]] <- data.frame(x=x, y=y, stringsAsFactors = FALSE)
      }
      for (el in obj) acc <- collect_line_traces(el, acc)
    }
    acc
  }
  for (i in seq_along(nodes)){
    js <- tryCatch(jsonlite::fromJSON(xml2::xml_text(nodes[[i]]), simplifyVector = FALSE), error = function(e) NULL)
    if (is.null(js)) next
    acc <- collect_line_traces(js, list())
    if (length(acc))
      return(bind_rows(acc) %>% arrange(x, y) %>% mutate(cell_line = cell_line, metric = metric))
  }
  NULL
}
collect_smooth_lines <- function(base_clade_dir, which_metric = c("ic50","auc")){
  which_metric <- match.arg(which_metric)
  cls <- list_cell_lines(base_clade_dir); if (!length(cls)) return(tibble())
  map_dfr(cls, ~ extract_smooth_from_html(.x, base_clade_dir, which_metric))
}

# ---- Load points (TSV) & normalize columns ----
prep_clade_data <- function(base_clade_dir){
  cls <- list_cell_lines(base_clade_dir)
  ic50_tbl <- map_dfr(cls, ~ read_insilico_file(.x, base_clade_dir, "ic50"))
  auc_tbl  <- map_dfr(cls, ~ read_insilico_file(.x, base_clade_dir, "auc"))
  if (nrow(ic50_tbl)){
    if (!"medIC50" %in% names(ic50_tbl)){
      cand <- intersect(c("medianIC50","MedianIC50","medic50","MEDIC50"), names(ic50_tbl))
      if (length(cand) == 1) ic50_tbl <- ic50_tbl %>% rename(medIC50 = !!cand)
    }
    if (!"StronglyPredicted" %in% names(ic50_tbl)) ic50_tbl <- ic50_tbl %>% mutate(StronglyPredicted = sRGES <= -0.20)
    ic50_tbl <- ic50_tbl %>% filter(is.finite(medIC50), !is.na(sRGES))
  }
  if (nrow(auc_tbl)){
    if (!"medauc" %in% names(auc_tbl)){
      cand <- intersect(c("medianAUC","MedianAUC","medAUC","MEDAUC","AUC","auc"), names(auc_tbl))
      if (length(cand) == 1) auc_tbl <- auc_tbl %>% rename(medauc = !!cand)
    }
    if (!"StronglyPredicted" %in% names(auc_tbl)) auc_tbl <- auc_tbl %>% mutate(StronglyPredicted = sRGES <= -0.20)
    auc_tbl <- auc_tbl %>% filter(is.finite(medauc), !is.na(sRGES))
  }
  list(ic50 = ic50_tbl, auc = auc_tbl)
}



# --- FDA lists to highlight (only Launched) ---
RGEs_A7 <- readRDS("~/CESC_Network/6_OCTAD/6_3_RES/A7/RES_A7_common3_collapsed_FDA_Launched.rds")
RGEs_A9 <- readRDS("~/CESC_Network/6_OCTAD/6_3_RES/A9/RES_A9_common3_collapsed_FDA_Launched.rds")

fda_set_A7 <- tolower(trimws(RGEs_A7$pert_iname))
fda_set_A9 <- tolower(trimws(RGEs_A9$pert_iname))

add_fda_flag <- function(tbl, fda_set) {
  if (is.null(tbl) || !nrow(tbl)) return(tbl)
  tbl %>%
    mutate(iname_l = tolower(trimws(pert_iname)),
           fda_flag = ifelse(iname_l %in% fda_set, "FDA", "Other"))
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

# ---- Helper: one label per drug (median position) ----
unique_labels_ic50 <- function(tbl){
  tmp <- tbl %>%
    filter(!is.na(pert_iname), is.finite(sRGES), is.finite(medIC50)) %>%
    mutate(y = log10(medIC50))
  if (!nrow(tmp)) return(tmp[0, c("sRGES","y","pert_iname","StronglyPredicted")])
  tmp %>%
    group_by(pert_iname) %>%
    summarise(sRGES = median(sRGES, na.rm = TRUE),
              y     = median(y,     na.rm = TRUE),
              StronglyPredicted = ifelse(mean(StronglyPredicted == "Yes", na.rm = TRUE) >= 0.5, "Yes", "No"),
              .groups = "drop")
}
unique_labels_auc <- function(tbl){
  tmp <- tbl %>%
    filter(!is.na(pert_iname), is.finite(sRGES), is.finite(medauc))
  if (!nrow(tmp)) return(tmp[0, c("sRGES","y","pert_iname","StronglyPredicted")])
  tmp %>%
    group_by(pert_iname) %>%
    summarise(sRGES = median(sRGES, na.rm = TRUE),
              y     = median(medauc, na.rm = TRUE),
              StronglyPredicted = ifelse(mean(StronglyPredicted == "Yes", na.rm = TRUE) >= 0.5, "Yes", "No"),
              .groups = "drop")
}

# ---- Plotters ----
# Facets (sin cambios): points + labels por punto + línea extraída del HTML
plot_facets_ic50 <- function(tbl, title_prefix, metrics_tbl, smooth_df, fda_set){
  if (!nrow(tbl)) return(NULL)
  ann <- metrics_tbl %>% filter(metric=="IC50") %>% dplyr::select(cell_line, label) %>% distinct()
  
  # Etiquetas SOLO FDA y SOLO StronglyPredicted = "Yes" (sRGES <= -0.2)
  labs_df <- tbl %>%
    mutate(iname_l = tolower(trimws(pert_iname))) %>%
    filter(iname_l %in% fda_set,
           StronglyPredicted == "Yes",
           is.finite(sRGES), is.finite(medIC50)) %>%
    group_by(cell_line, pert_iname) %>%
    summarise(sRGES = median(sRGES, na.rm = TRUE),
              y     = median(log10(medIC50), na.rm = TRUE),
              .groups = "drop")
  
  ggplot(tbl, aes(x=sRGES, y=log10(medIC50), color=StronglyPredicted)) +
    geom_point() +
    ggrepel::geom_text_repel(
      data = labs_df,
      inherit.aes = FALSE,
      aes(x = sRGES, y = y, label = pert_iname),
      size=2.3, max.overlaps=5, segment.color='grey70'
    ) +
    { if (nrow(smooth_df)) geom_path(data=smooth_df, aes(x=x,y=y,group=cell_line),
                                     inherit.aes=FALSE, color="black", linewidth=0.6) } +
    facet_wrap(~cell_line, scales="free") +
    geom_text(data=ann, inherit.aes=FALSE, aes(x=-Inf, y=Inf, label=label),
              hjust=-0.05, vjust=1.1, size=3.1) +
    labs(title=sprintf("%s: sRGES vs log10(IC50) by Cell Line", title_prefix),
         y="log10(medIC50)", x="sRGES") +
    coord_cartesian(ylim = c(NA, 15)) +   # <-- recorta arriba en 15
    theme_bw()
}


plot_facets_auc <- function(tbl, title_prefix, metrics_tbl, smooth_df, fda_set){
  if (!nrow(tbl)) return(NULL)
  ann <- metrics_tbl %>% filter(metric=="AUC") %>% dplyr::select(cell_line, label) %>% distinct()
  
  # Etiquetas SOLO FDA
  labs_df <- tbl %>%
    mutate(iname_l = tolower(trimws(pert_iname))) %>%
    filter(iname_l %in% fda_set,
           is.finite(sRGES), is.finite(medauc)) %>%
    group_by(cell_line, pert_iname) %>%
    summarise(sRGES = median(sRGES, na.rm = TRUE),
              y     = median(medauc, na.rm = TRUE),
              .groups = "drop")
  
  ggplot(tbl, aes(x=sRGES, y=medauc, color=StronglyPredicted)) +
    geom_point() +
    # SOLO etiquetar FDA:
    ggrepel::geom_text_repel(
      data = labs_df,
      inherit.aes = FALSE,
      aes(x = sRGES, y = y, label = pert_iname),
      size=2.3, max.overlaps=5, segment.color='grey70'
    ) +
    { if (nrow(smooth_df)) geom_path(data=smooth_df, aes(x=x,y=y,group=cell_line),
                                     inherit.aes=FALSE, color="black", linewidth=0.6) } +
    facet_wrap(~cell_line, scales="free") +
    geom_text(data=ann, inherit.aes=FALSE, aes(x=-Inf, y=Inf, label=label),
              hjust=-0.05, vjust=1.1, size=3.1) +
    labs(title=sprintf("%s: sRGES vs AUC by Cell Line", title_prefix),
         y="AUC", x="sRGES") + theme_bw()
}

# Global (modificado): points + **one label per drug** + single lm line
plot_global_ic50 <- function(tbl, title_prefix){
  if (is.null(tbl) || !nrow(tbl)) return(NULL)
  
  # Etiquetas SOLO FDA y SOLO StronglyPredicted = "Yes"
  labs_df <- tbl %>%
    filter(fda_flag == "FDA",
           StronglyPredicted == "Yes",
           is.finite(sRGES), is.finite(medIC50)) %>%
    group_by(pert_iname) %>%
    summarise(sRGES = median(sRGES, na.rm = TRUE),
              y     = median(log10(medIC50), na.rm = TRUE),
              .groups = "drop")
  
  ggplot(tbl, aes(x = sRGES, y = log10(medIC50))) +
    geom_point(aes(color = StronglyPredicted, shape = fda_flag),
               alpha = ifelse(tbl$fda_flag == "FDA", 0.9, 0.35), size = 1.9) +
    # Recta con TODOS los puntos (no solo FDA)
    geom_smooth(aes(group = 1), method = "lm", se = TRUE,
                color = "black", fill = "grey80") +
    ggrepel::geom_text_repel(
      data = labs_df, inherit.aes = FALSE,
      aes(x = sRGES, y = y, label = pert_iname, color = NULL),
      size = 2.6, segment.color = 'grey80', max.overlaps = Inf
    ) +
    scale_shape_manual(values = c(FDA = 17, Other = 16), name = "FDA status") +
    labs(title = sprintf("%s: Global sRGES vs log10(IC50)", title_prefix),
         y = "log10(medIC50)", x = "sRGES") +
    coord_cartesian(ylim = c(NA, 15)) +   # <-- recorta arriba en 15
    theme_minimal(base_size = 14)
}


plot_global_auc <- function(tbl, title_prefix){
  if (is.null(tbl) || !nrow(tbl)) return(NULL)
  
  labs_df <- tbl %>%
    filter(fda_flag == "FDA", is.finite(sRGES), is.finite(medauc)) %>%
    group_by(pert_iname) %>%
    summarise(sRGES = median(sRGES, na.rm = TRUE),
              y     = median(medauc,     na.rm = TRUE),
              .groups = "drop")
  
  ggplot(tbl, aes(x = sRGES, y = medauc)) +
    geom_point(aes(color = StronglyPredicted, shape = fda_flag),
               alpha = ifelse(tbl$fda_flag == "FDA", 0.9, 0.35), size = 1.9) +
    ## <- CAMBIO: línea con TODOS los puntos
    geom_smooth(aes(group = 1), method = "lm", se = TRUE,
                color = "black", fill = "grey80") +
    ggrepel::geom_text_repel(
      data = labs_df, inherit.aes = FALSE,
      aes(x = sRGES, y = y, label = pert_iname, color = NULL),
      size = 2.6, segment.color = 'grey80', max.overlaps = Inf
    ) +
    scale_shape_manual(values = c(FDA = 17, Other = 16), name = "FDA status") +
    labs(title = sprintf("%s: Global sRGES vs AUC", title_prefix),
         y = "AUC", x = "sRGES") +
    theme_minimal(base_size = 14)
}

# ---- Build plots ----
p_ic50_A7 <- plot_facets_ic50(data_A7$ic50, "A7", metrics_A7, smooth_A7_ic50, fda_set_A7)
p_auc_A7  <- plot_facets_auc (data_A7$auc,  "A7", metrics_A7, smooth_A7_auc,  fda_set_A7)
p_ic50_A9 <- plot_facets_ic50(data_A9$ic50, "A9", metrics_A9, smooth_A9_ic50, fda_set_A9)
p_auc_A9  <- plot_facets_auc (data_A9$auc,  "A9", metrics_A9, smooth_A9_auc,  fda_set_A9)

# Añadir flag a los data.frames usados en los plots globales
data_A7_ic50_fda <- add_fda_flag(data_A7$ic50, fda_set_A7)
data_A7_auc_fda  <- add_fda_flag(data_A7$auc,  fda_set_A7)
data_A9_ic50_fda <- add_fda_flag(data_A9$ic50, fda_set_A9)
data_A9_auc_fda  <- add_fda_flag(data_A9$auc,  fda_set_A9)

# Construir los plots globales con el flag (facet siguen igual que tenías)
p_ic50_A7_glob  <- plot_global_ic50(data_A7_ic50_fda, "A7")
p_auc_A7_glob   <- plot_global_auc (data_A7_auc_fda,  "A7")
p_ic50_A9_glob  <- plot_global_ic50(data_A9_ic50_fda, "A9")
p_auc_A9_glob   <- plot_global_auc (data_A9_auc_fda,  "A9")

# ---- Save (PDF + 1200 dpi PNG) ----
save_plot_all <- function(plot_obj, stem){
  if (is.null(plot_obj)) return(invisible(NULL))
  ggsave(file.path(plots_out_dir, paste0(stem, ".pdf")), plot_obj, device = cairo_pdf)
  ggsave(file.path(plots_out_dir, paste0(stem, ".png")), plot_obj, dpi = 1200)
}

save_plot_all(p_ic50_A7,      "A7_facet_ic50")
save_plot_all(p_auc_A7,       "A7_facet_auc")
save_plot_all(p_ic50_A9,      "A9_facet_ic50")
save_plot_all(p_auc_A9,       "A9_facet_auc")
save_plot_all(p_ic50_A7_glob, "A7_global_ic50")
save_plot_all(p_auc_A7_glob,  "A7_global_auc")
save_plot_all(p_ic50_A9_glob, "A9_global_ic50")
save_plot_all(p_auc_A9_glob,  "A9_global_auc")

# (Optional interactive display)
for (p in list(p_ic50_A7,p_auc_A7,p_ic50_A9,p_auc_A9,p_ic50_A7_glob,p_auc_A7_glob,p_ic50_A9_glob,p_auc_A9_glob))
  if (!is.null(p)) print(p)







length(unique(tolower(data_A7_ic50_fda$pert_iname)))
length(unique(tolower(data_A7_auc_fda$pert_iname)))



