############################################################
## 8_7_DGIDB_ClosestProximity_individual_plots.R
############################################################

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

outdir <- "~/CESC_Network/6_DGIDB/"
setwd(outdir)

#-----------------------------------------------------------
# 1) Cargar matrices módulo × fármaco
#-----------------------------------------------------------

mat_A7_df <- read_tsv("8_5_A7_ModuleDrug_dc_matrix.tsv", show_col_types = FALSE)
mat_A9_df <- read_tsv("8_5_A9_ModuleDrug_dc_matrix.tsv", show_col_types = FALSE)

mat_A7 <- mat_A7_df %>% column_to_rownames(var = colnames(mat_A7_df)[1]) %>% as.matrix()
mat_A9 <- mat_A9_df %>% column_to_rownames(var = colnames(mat_A9_df)[1]) %>% as.matrix()

#-----------------------------------------------------------
# 2) Auto–wrap nombres largos (AMBAS redes)
#-----------------------------------------------------------

wrap_names <- function(n, width = 18) {
  sapply(n, function(x) paste(strwrap(x, width = width), collapse = "\n"))
}

colnames(mat_A7) <- wrap_names(colnames(mat_A7))
colnames(mat_A9) <- wrap_names(colnames(mat_A9))

#-----------------------------------------------------------
# 3) Cargar zc full
#-----------------------------------------------------------

results_A7 <- read_tsv("8_5_ClosestProximity_A7_FULL.tsv", show_col_types = FALSE)
results_A9 <- read_tsv("8_5_ClosestProximity_A9_FULL.tsv", show_col_types = FALSE)


#-----------------------------------------------------------
# 4) Función para heatmap individual
#-----------------------------------------------------------

make_single_heatmap <- function(mat, res_df, cohort_label,
                                pdf_file, png_file) {
  
  # Asegurar emparejamiento con nombres originales (pre-wrap)
  original_names <- res_df$drug
  names(res_df$zc) <- wrap_names(original_names)  # match con columnas wrap
  
  common <- intersect(colnames(mat), names(res_df$zc))
  mat <- mat[, common, drop = FALSE]
  zc_vec <- res_df$zc[colnames(mat)]
  
  # ---- Escala zc ----
  z_thresh <- -1.5
  z_min <- min(zc_vec, na.rm = TRUE)
  if (z_min > z_thresh) z_min <- z_thresh
  
  col_zc <- colorRamp2(
    c(z_thresh, z_min),
    c("#e5f5e0", "#00441b")
  )
  
  top_ha <- HeatmapAnnotation(
    zc = zc_vec,
    col = list(zc = col_zc),
    annotation_height = unit(0.8, "cm"),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold")
  )
  
  # ---- Escala dc ----
  dc_min <- min(mat, na.rm = TRUE)
  dc_max <- max(mat, na.rm = TRUE)
  
  col_dc <- colorRamp2(
    c(dc_max, dc_min),
    c("#f7fbff", "#08306B")
  )
  
  ht <- Heatmap(
    mat,
    name = "dc (closest)",
    col = col_dc,
    top_annotation = top_ha,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 11, fontface = "bold"),
    column_names_gp = gpar(fontsize = 11, fontface = "bold"),
    column_title = paste0("HPV ", cohort_label, " | Module-level drug proximity"),
    column_title_gp = gpar(fontsize = 18, fontface = "bold"),
    heatmap_legend_param = list(
      title = "dc (closest)",
      title_gp = gpar(fontsize = 12, fontface = "bold"),
      labels_gp = gpar(fontsize = 10)
    )
  )
  
  # ---- Guardar PDF (más alto) ----
  pdf(pdf_file, width = 16, height = 9)
  draw(ht, merge_legend = TRUE,
       heatmap_legend_side = "right",
       annotation_legend_side = "right")
  dev.off()
  
  # ---- Guardar PNG ----
  png(png_file, width = 3000, height = 1600, res = 180)
  draw(ht, merge_legend = TRUE,
       heatmap_legend_side = "right",
       annotation_legend_side = "right")
  dev.off()
  
  invisible(ht)
}


#-----------------------------------------------------------
# 5) Dibujar A7 y A9
#-----------------------------------------------------------

make_single_heatmap(
  mat_A7, results_A7, "A7",
  file.path(outdir, "8_7_A7_ModuleDrug_Proximity_fixed.pdf"),
  file.path(outdir, "8_7_A7_ModuleDrug_Proximity_fixed.png")
)

make_single_heatmap(
  mat_A9, results_A9, "A9",
  file.path(outdir, "8_7_A9_ModuleDrug_Proximity_fixed.pdf"),
  file.path(outdir, "8_7_A9_ModuleDrug_Proximity_fixed.png")
)

message("Listo — ambos heatmaps individuales generados correctamente.")
