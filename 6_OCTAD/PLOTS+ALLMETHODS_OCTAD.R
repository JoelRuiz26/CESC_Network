# ============================================================
# Drug consensus across methods/clades + averaged sRGES + heatmaps
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(pheatmap)
  # opcional para PNG nítido
  have_ragg <- requireNamespace("ragg", quietly = TRUE)
  # para similitud de líneas (validación visual)
  library(octad)
  library(octad.db)
})

setwd("~/CESC_Network/6_OCTAD/")

# ---------------- Parameters ----------------
methods           <- c("EdgeR","DESeq2","limma")
clades            <- c("A7","A9")
srges_cutoff_tag  <- "0.25"    # coincide con los archivos *_results_0.25.rds
drug_col_candidates <- c("pert_iname","drug","compound","pert","name")
out_dir           <- "6_5_ConsensusHeatmaps"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------- Helpers ----------------
pick_drug_col <- function(df) {
  cand <- intersect(drug_col_candidates, names(df))
  if (length(cand) == 0) stop("No recognizable drug name column in sRGES table.")
  cand[[1]]
}
read_srges_filtered <- function(clade, method, tag = srges_cutoff_tag) {
  f <- sprintf("6_2_1_OCTAD_%s_%s_results_%s.rds", clade, method, tag)
  if (!file.exists(f)) stop("File not found: ", f)
  x <- readRDS(f)
  stopifnot("sRGES" %in% names(x))
  dcol <- pick_drug_col(x)
  x %>%
    dplyr::select(drug = all_of(dcol), sRGES) %>%
    dplyr::mutate(drug = as.character(drug)) %>%
    dplyr::group_by(drug) %>%              # si hay duplicados, promediamos
    dplyr::summarise(sRGES = mean(sRGES, na.rm = TRUE), .groups = "drop")
}

# Construye listas por clado: tablas por método, sets unión e intersección, y media de sRGES
build_clade_lists <- function(clade) {
  tbls <- setNames(lapply(methods, \(m) read_srges_filtered(clade, m)), methods)
  
  # Intersección (comunes dentro del clado)
  common3 <- Reduce(intersect, lapply(tbls, \(d) d$drug))
  # Unión (para definir “únicos del clado” frente al otro)
  union3  <- unique(unlist(lapply(tbls, \(d) d$drug)))
  
  # Matriz con columnas por método y columna 'mean' (promedio)
  mat <- tbls %>%
    purrr::imap(~dplyr::rename(.x, !!paste0(clade,"_",.y) := sRGES)) %>%
    purrr::reduce(full_join, by = "drug")
  
  # mean dentro de clado (sobre columnas de sus 3 métodos)
  clade_cols <- paste0(clade, "_", methods)
  mat[[paste0(clade, "_mean")]] <- mat %>% dplyr::select(all_of(clade_cols)) %>%
    as.matrix() %>% apply(1, function(v) mean(v, na.rm = TRUE))
  
  list(tables = tbls, common = common3, union = union3, mat = mat)
}

# ---------------- Load per-clade data ----------------
clade_data <- setNames(lapply(clades, build_clade_lists), clades)

# Comunes entre clados (intersección de intersecciones)
common_both <- intersect(clade_data$A7$common, clade_data$A9$common)

# Únicos por clado (en unión de clado, excluyendo la unión del otro)
unique_A7 <- setdiff(clade_data$A7$union, clade_data$A9$union)
unique_A9 <- setdiff(clade_data$A9$union, clade_data$A7$union)

# ---------------- Promedios por clado (ya incluidos en mat) ----------------
# Extrae data frames con columnas por método + mean para cada clado
mat_A7 <- clade_data$A7$mat
mat_A9 <- clade_data$A9$mat

# ---------------- Selección para heatmaps ----------------
# Top 20 comunes (por clado), ordenados por sRGES medio más negativo
top20_common_A7 <- mat_A7 %>%
  dplyr::filter(drug %in% clade_data$A7$common) %>%
  dplyr::arrange(!!sym("A7_mean")) %>% dplyr::slice_head(n = 20) %>% dplyr::pull(drug)

top20_common_A9 <- mat_A9 %>%
  dplyr::filter(drug %in% clade_data$A9$common) %>%
  dplyr::arrange(!!sym("A9_mean")) %>% dplyr::slice_head(n = 20) %>% dplyr::pull(drug)

# Top 20 únicos por clado (ordenados por su media de clado, más negativo primero)
top20_unique_A7 <- mat_A7 %>%
  dplyr::filter(drug %in% unique_A7) %>%
  dplyr::arrange(!!sym("A7_mean")) %>% dplyr::slice_head(n = 20) %>% dplyr::pull(drug)

top20_unique_A9 <- mat_A9 %>%
  dplyr::filter(drug %in% unique_A9) %>%
  dplyr::arrange(!!sym("A9_mean")) %>% dplyr::slice_head(n = 20) %>% dplyr::pull(drug)

# Conjunto final de filas a mostrar en el heatmap
rows_heat <- unique(c(top20_common_A7, top20_common_A9, top20_unique_A7, top20_unique_A9))

# ---------------- Construcción de matriz para heatmap ----------------
# Unimos ambas matrices por 'drug' y nos quedamos con columnas de interés
all_mat <- full_join(mat_A7, mat_A9, by = "drug")

# Orden de columnas (métodos + mean por clado)
col_order <- c(paste0("A7_", methods), "A7_mean",
               paste0("A9_", methods), "A9_mean")

# Submatriz y conversión a matrix
hm_df <- all_mat %>%
  dplyr::filter(drug %in% rows_heat) %>%
  dplyr::select(drug, all_of(col_order))

hm_mat <- hm_df %>%
  tibble::column_to_rownames("drug") %>%
  as.matrix()

# Anotación de filas (categoría)
row_cat <- rep(NA_character_, nrow(hm_mat))
names(row_cat) <- rownames(hm_mat)
row_cat[rownames(hm_mat) %in% top20_common_A7] <- "A7_common_top20"
row_cat[rownames(hm_mat) %in% top20_common_A9] <- "A9_common_top20"
row_cat[rownames(hm_mat) %in% top20_unique_A7] <- "A7_unique_top20"
row_cat[rownames(hm_mat) %in% top20_unique_A9] <- "A9_unique_top20"

ann_row <- data.frame(Category = factor(row_cat,
                                        levels = c("A7_common_top20","A9_common_top20",
                                                   "A7_unique_top20","A9_unique_top20")))
rownames(ann_row) <- rownames(hm_mat)

# Paleta divergente centrada en 0 (azul = mejor sRGES negativo)
max_abs <- max(abs(hm_mat), na.rm = TRUE)
bk <- seq(-max_abs, max_abs, length.out = 101)
pal <- colorRampPalette(c("#08306B","#4292C6","#F7FBFF","#FCBBA1","#CB181D"))(100)

# Heatmap
png_fun <- if (have_ragg) ragg::agg_png else png
png_fun(filename = file.path(out_dir, "6_5_1_Heatmap_Common_Unique_A7_A9.png"),
        width = 2200, height = 1600, res = 300)
pheatmap(hm_mat,
         color = pal, breaks = bk, cluster_rows = TRUE, cluster_cols = FALSE,
         annotation_row = ann_row,
         main = "sRGES (negative is better) — A7/A9 by method and mean",
         fontsize = 12, fontsize_row = 8.5, fontsize_col = 11,
         border_color = NA)
dev.off()

# ---------------- Export lists requested ----------------
# 1) Common across the two clade-level intersections (= “en común en las 4 firmas” interpretable como intersección A7∩A9)
#    Si quieres “comunes en las 6 (A7×3 y A9×3)”, usa intersect de las 6 listas en vez de common_both.
writeLines(common_both, file.path(out_dir, "common_drugs_both_clades.txt"))

# 2) Promedios por clado de los fármacos comunes (ya están en hm_mat; aquí exportamos tablas)
common_A7_tbl <- mat_A7 %>%
  dplyr::filter(drug %in% clade_data$A7$common) %>%
  dplyr::select(drug, starts_with("A7_")) %>%
  arrange(A7_mean)
common_A9_tbl <- mat_A9 %>%
  dplyr::filter(drug %in% clade_data$A9$common) %>%
  dplyr::select(drug, starts_with("A9_")) %>%
  arrange(A9_mean)

readr::write_csv(common_A7_tbl, file.path(out_dir, "A7_common_drugs_with_means.csv"))
readr::write_csv(common_A9_tbl, file.path(out_dir, "A9_common_drugs_with_means.csv"))

# ============================================================
# Validation plot with cell lines: barplots of medcor per clade
# (Recomputes similarity just for visualization; comment if not needed)
# ============================================================
phenoDF <- octad.db::get_ExperimentHub_data("EH7274")
sample_clado_A7 <- readRDS("6_0_A7_Samples.rds")
sample_clado_A9 <- readRDS("6_0_A9_Samples.rds")
case_A7 <- intersect(sample_clado_A7, phenoDF$sample.id)
case_A9 <- intersect(sample_clado_A9, phenoDF$sample.id)

get_medcor_tbl <- function(case_ids, thr = 0.31) {
  sim <- computeCellLine(case_id = case_ids, source = "octad.small")
  sim <- sim %>% tibble::rownames_to_column("cell_id")
  sim %>% dplyr::filter(medcor > thr) %>% dplyr::arrange(desc(medcor))
}

med_A7 <- get_medcor_tbl(case_A7, thr = 0.31)
med_A9 <- get_medcor_tbl(case_A9, thr = 0.31)

plot_medcor <- function(df, title, out_png) {
  df$cell_id <- factor(df$cell_id, levels = df$cell_id[order(df$medcor)])
  p <- ggplot(df, aes(x = cell_id, y = medcor)) +
    geom_col(fill = "#3B4CC0", alpha = 0.85) +
    coord_flip() +
    labs(title = title, x = "Cell line", y = "Median correlation (medcor)") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold", size = 16),
          axis.text.y = element_text(size = 10))
  ggsave(out_png, p, width = 7.5, height = 6.5, dpi = 300)
  p
}

plot_medcor(med_A7, "Validation — Similar cell lines used (HPV-A7 clade)",
            file.path(out_dir, "6_5_2_Validation_CellLines_A7.png"))
plot_medcor(med_A9, "Validation — Similar cell lines used (HPV-A9 clade)",
            file.path(out_dir, "6_5_2_Validation_CellLines_A9.png"))

# ---------------- (Optional) print small summaries ----------------
cat("Common within A7 (3 methods):", length(clade_data$A7$common), "\n")
cat("Common within A9 (3 methods):", length(clade_data$A9$common), "\n")
cat("Common in both clades (intersection of intersections):", length(common_both), "\n")
cat("Unique A7 (vs A9):", length(unique_A7), " | Unique A9 (vs A7):", length(unique_A9), "\n")
