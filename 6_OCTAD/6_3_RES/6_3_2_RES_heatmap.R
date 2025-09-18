

# =========================================
# Heatmap horizontal (A7 vs A9) con leyendas en TOP
# =========================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# ---------------- Config ----------------
base_res   <- "~/CESC_Network/6_OCTAD/6_3_RES"
topN       <- Inf  # usa Inf para TODOS


collapsed_A7 <- readRDS("/home/jjruiz/CESC_Network/6_OCTAD/6_3_RES/A7/RES_A7_common3_collapsed_FDA_Launched.rds")
collapsed_A7 <- collapsed_A7 %>% filter(sRGES <= -0.2)
collapsed_A9 <- readRDS("/home/jjruiz/CESC_Network/6_OCTAD/6_3_RES/A9/RES_A9_common3_collapsed_FDA_Launched.rds")
collapsed_A9 <- collapsed_A9 %>% filter(sRGES <= -0.2)



# --------- Tabla base (A7 vs A9) ----------
wide <- full_join(
  collapsed_A7 %>% dplyr::select(pert_iname, A7 = sRGES),
  collapsed_A9 %>% dplyr::select(pert_iname, A9 = sRGES),
  by = "pert_iname"
) %>%
  mutate(
    membership = case_when(
      !is.na(A7) & !is.na(A9) ~ "Common",
      !is.na(A7) &  is.na(A9) ~ "A7-only",
      is.na(A7) & !is.na(A9) ~ "A9-only",
      TRUE ~ NA_character_
    ),
    min_sRGES = if_else(is.na(A7) & is.na(A9), NA_real_, pmin(A7, A9, na.rm = TRUE))
  )

top_tbl <- wide %>%
  filter(!is.na(min_sRGES)) %>%
  mutate(membership = factor(membership,
                             levels = c("Common","A7-only","A9-only"),
                             labels = c("Common","HPV-A7-only","HPV-A9-only"))) %>%
  arrange(membership, min_sRGES)
if (is.finite(topN)) top_tbl <- slice(top_tbl, 1:min(topN, nrow(top_tbl)))

mat <- top_tbl %>%
  dplyr::select(pert_iname, A7, A9) %>%
  column_to_rownames("pert_iname") %>%
  as.matrix()
colnames(mat) <- c("  HPV-A7","  HPV-A9")

# --------- Rotación a horizontal + “oscuros” a la derecha ----------
# --------- Rotación a horizontal + “oscuros” a la derecha ----------
mat_h <- t(mat)[, rev(rownames(mat)), drop = FALSE]

# --------- Colores (anclado en 0) ----------
vals    <- as.numeric(mat_h)
neg_min <- min(vals, na.rm = TRUE)
pos_max <- max(vals, na.rm = TRUE)
neg_col  <- "#2166AC"; zero_col <- "#FFF3B0"; pos_col <- "#D7301F"
col_fun <- if (pos_max <= 0 || !is.finite(pos_max)) {
  circlize::colorRamp2(c(neg_min, 0),            c(neg_col, zero_col))
} else {
  circlize::colorRamp2(c(neg_min, 0, pos_max),   c(neg_col, zero_col, pos_col))
}

# --------- Anotación arriba (membership) ----------
col_membership_map <- setNames(as.character(top_tbl$membership), rownames(mat))
col_membership     <- col_membership_map[colnames(mat_h)]
membership_cols <- c("Common"="#16A085","HPV-A7-only"="#F39C12","HPV-A9-only"="#8E44AD")
# --------- Anotación arriba (membership) ----------
top_anno <- HeatmapAnnotation(
  membership = col_membership,
  col = list(membership = membership_cols),
  annotation_name_side   = "left",
  annotation_name_gp     = gpar(fontsize = 12, fontface = "bold"),
  annotation_name_offset = unit(2, "mm"),   # ← lo despega de la orilla izq
  height = unit(6, "mm"),                   # ← tira un poquito más de alto a la franja
  show_legend = FALSE
)


# --------- Heatmap ----------
ht <- Heatmap(
  mat_h,
  name = "sRGES",
  col = col_fun,
  na_col = "grey90",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  top_annotation = top_anno,
  show_heatmap_legend = FALSE,
  border = NA,
  
  # filas
  row_names_side = "left",
  row_names_gp   = gpar(fontsize = 14, fontface = "bold"),
  
  # columnas (fármacos)
  column_names_side       = "bottom",
  column_names_rot        = 55,
  column_names_gp         = gpar(fontsize = 10, fontface = "bold"),
  column_names_max_height = unit(26, "mm"),  # ← más espacio para que NO se corten
  
  # cuerpo del heatmap un poquito más alto
  height = unit(2.5, "in"),                  # ← súbelo (2.3–2.7 in) si quieres
  heatmap_legend_param = list(title = "sRGES", direction = "horizontal")
)


# --------- Leyendas arriba, horizontales ----------
# Solo invertimos visualmente la barra (0 a la izq, negativos a la der)
col_fun_leg <- circlize::colorRamp2(c(0, 1), c(zero_col, neg_col))  # amarillo→azul

lgd_sRGES <- Legend(
  title     = "sRGES",
  col_fun   = col_fun_leg,
  direction = "horizontal",
  at        = seq(0, 1, by = 0.25),
  labels    = c("0", "-0.1", "-0.2", "-0.3", "-0.4")
)
lgd_mem <- Legend(
  title = "Membership",
  at = names(membership_cols),
  legend_gp = gpar(fill = membership_cols),
  direction = "horizontal"
)
lgd_packed <- packLegend(lgd_sRGES, lgd_mem, direction = "horizontal", gap = unit(6, "mm"))

draw(
  ht,
  heatmap_legend_list = list(lgd_packed),
  heatmap_legend_side = "top",
  annotation_legend_side = "top",
  merge_legend = FALSE,
  padding = unit(c(2, 2, 2, 2), "mm")
)




# --------------------------
# 6) Save heatmap — mismo ancho, menor alto
# --------------------------
suffix  <- if (is.finite(topN)) paste0("top", topN) else "ALL"
pdf_out <- file.path(base_res, paste0("6_3_3_RES_heatmap_A7_A9_HORIZONTAL_", suffix, ".pdf"))
svg_out <- file.path(base_res, paste0("6_3_3_RES_heatmap_A7_A9_HORIZONTAL_", suffix, ".svg"))
png_out <- file.path(base_res, paste0("6_3_3_RES_heatmap_A7_A9_HORIZONTAL_", suffix, ".png"))

# Reusa el ancho que venías calculando…
n_rows_vert <- nrow(mat)
lab_w_vert_in <- grid::convertWidth(
  max_text_width(rownames(mat), gp = gpar(fontsize = 12, fontface = "bold")) * 1.2,
  "in", valueOnly = TRUE
)
body_w_vert_in   <- 3.8
legend_w_vert_in <- 2.7
pdf_w_vert <- body_w_vert_in + lab_w_vert_in + legend_w_vert_in
pdf_h_vert <- max(5.5, n_rows_vert * 0.165)

# ⬇️ Mantén el ancho rotado, pero fuerza un alto compacto
# --------- Guardado (ajusta la altura total del lienzo) ----------
# Mantén el mismo ancho que ya calculas:
# Mantén el mismo ancho; ajusta un poco la altura total
pdf_w <- pdf_h_vert
FIG_H_IN <- 5.6
pdf_h <- FIG_H_IN

# ---- PDF
pdf(file = pdf_out, width = pdf_w, height = pdf_h, onefile = TRUE)
draw(
  ht,
  heatmap_legend_list   = list(lgd_packed),
  heatmap_legend_side   = "top",
  annotation_legend_side= "top",
  merge_legend          = FALSE,
  padding               = grid::unit(c(4, 6, 26, 16), "mm")  # top, right, bottom, left
)
dev.off()

# ---- SVG
svg(filename = svg_out, width = pdf_w, height = pdf_h)
draw(
  ht,
  heatmap_legend_list   = list(lgd_packed),
  heatmap_legend_side   = "top",
  annotation_legend_side= "top",
  merge_legend          = FALSE,
  padding               = grid::unit(c(4, 6, 26, 16), "mm")
)
dev.off()

# ---- PNG
if (requireNamespace("ragg", quietly = TRUE)) {
  ragg::agg_png(png_out, width = pdf_w, height = pdf_h, units = "in", res = 1200)
} else {
  png(filename = png_out, width = pdf_w, height = pdf_h, units = "in",
      res = 1200, type = "cairo-png")
}
draw(
  ht,
  heatmap_legend_list   = list(lgd_packed),
  heatmap_legend_side   = "top",
  annotation_legend_side= "top",
  merge_legend          = FALSE,
  padding               = grid::unit(c(4, 6, 26, 16), "mm")
)
dev.off()

saveRDS(collapsed_A7, "/home/jjruiz/CESC_Network/6_OCTAD/6_3_RES/A7/RES_A7_common3_collapsed_FDA_Launched_0.20.rds")
saveRDS(collapsed_A9, "/home/jjruiz/CESC_Network/6_OCTAD/6_3_RES/A9/RES_A9_common3_collapsed_FDA_Launched_0.20.rds")



