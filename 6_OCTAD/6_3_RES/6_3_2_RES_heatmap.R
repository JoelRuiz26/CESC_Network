# =========================================
# Simple horizontal heatmap (A7 vs A9)
# Legends merged on top, sRGES 0 -> min, filas un poco más chaparras
# =========================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(ComplexHeatmap); library(circlize); library(grid)
})

# ---------- Inputs ----------
base_res <- "~/CESC_Network/6_OCTAD/6_3_RES"
collapsed_A7 <- readRDS(file.path(base_res, "A7/RES_A7_common3_collapsed_FDA_Launched.rds")) %>% filter(sRGES <= -0.2)
collapsed_A9 <- readRDS(file.path(base_res, "A9/RES_A9_common3_collapsed_FDA_Launched.rds")) %>% filter(sRGES <= -0.2)

# ---------- A7 vs A9 ----------
wide <- full_join(
  collapsed_A7 %>% select(pert_iname, A7 = sRGES),
  collapsed_A9 %>% select(pert_iname, A9 = sRGES),
  by = "pert_iname"
) %>%
  mutate(
    membership = case_when(
      !is.na(A7) & !is.na(A9) ~ "Common",
      !is.na(A7) &  is.na(A9) ~ "HPV-A7-only",
      is.na(A7) & !is.na(A9) ~ "HPV-A9-only",
      TRUE ~ NA_character_
    ),
    min_sRGES = if_else(is.na(A7) & is.na(A9), NA_real_, pmin(A7, A9, na.rm = TRUE))
  ) %>%
  filter(!is.na(min_sRGES)) %>%
  arrange(membership, min_sRGES)

mat <- wide %>% select(pert_iname, A7, A9) %>% column_to_rownames("pert_iname") %>% as.matrix()
colnames(mat) <- c("HPV-A7", "HPV-A9")
mat_h <- t(mat)[, rev(rownames(mat)), drop = FALSE]   # horizontal

# ---------- Colors: 0 -> min (amarillo -> azul) ----------
vals     <- as.numeric(mat_h)
data_min <- min(vals, na.rm = TRUE)
neg_col  <- "#2166AC"
zero_col <- "#FFF3B0"
col_fun  <- colorRamp2(c(-0.1, data_min), c(zero_col, neg_col))

# ---------- Top membership bar ----------
membership_cols <- c("Common"="#16A085","HPV-A7-only"="#F39C12","HPV-A9-only"="#8E44AD")
mem_vec <- setNames(wide$membership, wide$pert_iname)
top_anno <- HeatmapAnnotation(
  membership = mem_vec[colnames(mat_h)],
  col = list(membership = membership_cols),
  annotation_legend_param = list(membership = list(direction = "horizontal", title = "membership")),
  annotation_name_gp = gpar(fontsize = 0)  # ← oculta el texto "membership" del borde derecho
)


# ---------- Heatmap ----------
# --- en el Heatmap() ---
ht <- Heatmap(
  mat_h,
  name = "sRGES",
  col  = col_fun,
  na_col = "grey90",
  cluster_rows = FALSE, cluster_columns = FALSE,
  top_annotation = top_anno,
  
  show_column_names = TRUE,                 # <- asegura que se dibujen
  column_names_side       = "bottom",
  column_names_rot        = 55,
  column_names_gp         = gpar(fontsize = 14, fontface = "bold"),
  column_names_max_height = unit(75, "mm"), # <- más reserva para NO cortar
  
  row_names_side = "left",
  row_names_gp   = gpar(fontsize = 14, fontface = "bold"),
  
  height = unit(3.5, "in"),
  
  heatmap_legend_param = list(direction = "horizontal",
                              at = seq(-0.1, data_min, by = -0.1))
)

# ---------- Save ----------
# -------- AUTOSIZE + SAVE (no corta labels) --------
pdf_out <- file.path(base_res, "heatmap_A7_A9_horizontal.pdf")
svg_out <- file.path(base_res, "heatmap_A7_A9_horizontal.svg")
png_out <- file.path(base_res, "heatmap_A7_A9_horizontal.png")

# 1) Reservas por texto (medidas reales del texto)
lab_gp   <- gpar(fontsize = 14, fontface = "bold")
row_gp   <- gpar(fontsize = 14, fontface = "bold")
angle    <- 55 * pi/180

# Altura necesaria para las etiquetas de columnas rotadas
drug_text_w_in <- grid::convertWidth(
  ComplexHeatmap::max_text_width(colnames(mat_h), gp = lab_gp), "in", valueOnly = TRUE
)
lab_slot_in <- sin(angle) * drug_text_w_in + 0.25   # 0.25in de colchón

# Ancho necesario para las etiquetas de filas (HPV-A7/HPV-A9)
row_text_w_in <- grid::convertWidth(
  ComplexHeatmap::max_text_width(rownames(mat_h), gp = row_gp), "in", valueOnly = TRUE
)

# 2) Cuerpo del heatmap (tú controlas el ancho por columna)
per_col_in <- 0.18
body_w_in  <- ncol(mat_h) * per_col_in
body_h_in  <- 3.2    # filas “ligeramente más chaparras”

# 3) Altura extra para leyendas arriba
#legend_h_in <- 0.7

# 4) Tamaño final del lienzo
fig_w_in <- row_text_w_in + body_w_in + 0.8    # 0.8in margen derecho

#fig_h_in <- body_h_in + lab_slot_in + legend_h_in
legend_h_in  <- 0.9         # más espacio para la franja de leyendas
extra_top_in <- 0.25        # aire extra arriba (en pulgadas)
extra_bot_in <- 0.25        # aire extra abajo (en pulgadas)

fig_h_in <- body_h_in + lab_slot_in + legend_h_in + extra_top_in + extra_bot_in

# márgenes del dispositivo: top, right, bottom, left
pad <- grid::unit(c(10, 6, 10, 6), "mm")   # un poco más arriba y abajo

# 5) Pasa la misma reserva a ComplexHeatmap para que NUNCA recorte
ht <- draw(
  ht,
  heatmap_legend_side    = "top",
  annotation_legend_side = "top",
  merge_legends          = TRUE
)

# Aumenta el “slot” interno para etiquetas (evita cortes dentro del layout)
ht@ht_list[[1]]@column_names_param$max_height <- unit(lab_slot_in, "in")

# 6) Save en tres formatos (márgenes pequeños, ya no hace falta inflarlos)
#pad <- grid::unit(c(4, 6, 6, 6), "mm")  # top, right, bottom, left

# PDF
pdf(pdf_out, width = fig_w_in, height = fig_h_in, onefile = TRUE)
draw(ht, heatmap_legend_side="top", annotation_legend_side="top",
     merge_legends=TRUE, padding=pad)
dev.off()

# SVG
svg(svg_out, width = fig_w_in, height = fig_h_in)
draw(ht, heatmap_legend_side="top", annotation_legend_side="top",
     merge_legends=TRUE, padding=pad)
dev.off()

# PNG
if (requireNamespace("ragg", quietly = TRUE)) {
  ragg::agg_png(png_out, width = fig_w_in, height = fig_h_in, units = "in", res = 1200)
} else {
  png(png_out, width = fig_w_in, height = fig_h_in, units = "in", res = 1200, type = "cairo-png")
}
draw(ht, heatmap_legend_side="top", annotation_legend_side="top",
     merge_legends=TRUE, padding=pad)
dev.off()


saveRDS(collapsed_A7, "/home/jjruiz/CESC_Network/6_OCTAD/6_3_RES/A7/RES_A7_common3_collapsed_FDA_Launched_0.20.rds")
saveRDS(collapsed_A9, "/home/jjruiz/CESC_Network/6_OCTAD/6_3_RES/A9/RES_A9_common3_collapsed_FDA_Launched_0.20.rds")



