# ============================================================
# HEATMAP: Top 100 Genes Contributing to DeltaCon Distance
# A7 vs A9 — Publication Figure (Memory-Efficient Version)
# ============================================================

library(Matrix)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ── 0. Cargar resultados previos ──────────────────────────────
out_dir <- "~/CESC_Network/5_Network_analisis/5_3_DeltaCon/"
# 
# S_A7          <- readRDS(file.path(out_dir, "S_A7_matrix.rds"))
# S_A9          <- readRDS(file.path(out_dir, "S_A9_matrix.rds"))
# delta_con     <- readRDS(file.path(out_dir, "DeltaCon_distance_A7_A9.rds"))
# delta_con_sim <- readRDS(file.path(out_dir, "DeltaCon_similarity_A7_A9.rds"))
# 
# # ── 1. Alinear matrices ───────────────────────────────────────
# genes <- sort(intersect(rownames(S_A7), rownames(S_A9)))
# S1 <- S_A7[genes, genes]
# rm(S_A7); gc()
# 
# S2 <- S_A9[genes, genes]
# rm(S_A9); gc()
# 
# # ── 2. RANKING: contribución cuadrática por gen ───────────────
# sq_diff <- (sqrt(S1) - sqrt(S2))^2
# diag(sq_diff) <- 0
# gene_score <- Matrix::rowSums(sq_diff)
# rm(sq_diff); gc()
# 
# top100 <- names(sort(gene_score, decreasing = TRUE))[1:100]
# 
# # ── 3. HEATMAP: diferencia signada sobre top 100 ─────────────
# delta_mat <- as.matrix(sqrt(S1[top100, top100]) - sqrt(S2[top100, top100]))
# rm(S1, S2); gc()
# 
# diag(delta_mat) <- 0
# 
# # ── 4. Clipping de outliers (percentil 1–99) ──────────────────
# qlim <- quantile(delta_mat[delta_mat != 0], c(0.01, 0.99), na.rm = TRUE)
# delta_clipped <- pmax(pmin(delta_mat, qlim[2]), qlim[1])
# rm(delta_mat, qlim); gc()
# 
# # Escala simétrica alrededor de 0
# max_abs <- max(abs(delta_clipped))
# col_fun <- colorRamp2(
#   c(-max_abs, -max_abs * 0.3, 0, max_abs * 0.3, max_abs),
#   c("#2166AC", "#92C5DE", "white", "#F4A582", "#D6604D")
# )

# rm(max_abs); gc()
# 
# # ── 5. Anotación lateral: magnitud de contribución ───────────
# score_top    <- gene_score[top100]
# score_scaled <- (score_top - min(score_top)) / (max(score_top) - min(score_top))
# rm(gene_score, score_top, top100, genes); gc()
# 
# row_ann <- rowAnnotation(
#   Contribution = anno_barplot(
#     score_scaled,
#     bar_width   = 0.8,
#     gp          = gpar(fill = "#555555", col = NA),
#     axis_param  = list(gp = gpar(fontsize = 6)),
#     width       = unit(1.8, "cm")
#   ),
#   annotation_name_gp  = gpar(fontsize = 7, fontface = "bold"),
#   annotation_name_rot = 90
# )
# rm(score_scaled); gc()
# 
# ── 6. Título y subtítulo ─────────────────────────────────────
main_title <- ""
sub_title  <- sprintf(
  "Distance: %.4f   |   Similarity: %.6f",
  delta_con, delta_con_sim
)
# rm(delta_con, delta_con_sim); gc()
# 
# save.image("/STORAGE/csbig/jruiz/cesc_network/Delta_plot.RData")

load("~/CESC_Network/5_Network_analisis/5_3_DeltaCon/Delta_plot.RData")

# ── 1. Construcción del heatmap ───────────────────────────────
ht <- Heatmap(
  delta_clipped,
  name = "Δ Affinity\n(√S_A7 − √S_A9)",
  col  = col_fun,
  
  # Clustering
  clustering_distance_rows    = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method_rows      = "ward.D2",
  clustering_method_columns   = "ward.D2",
  show_row_dend      = TRUE,
  show_column_dend   = TRUE,
  row_dend_width     = unit(8, "mm"),
  column_dend_height = unit(8, "mm"),
  
  # Etiquetas
  show_row_names    = FALSE,
  show_column_names = FALSE,
  row_names_gp      = gpar(fontsize = 4.5),
  
  # SIN right_annotation
  
  # Título + subtítulo como column_title en dos líneas
  column_title = paste0(
    main_title, "\n",
    sub_title
  ),
  column_title_gp = gpar(
    fontsize  = c(25, 25),        # línea 1 grande, línea 2 igual de visible
    fontface  = c("bold", "bold"),
    col       = c("black", "#333333")
  ),
  
  # Leyenda
  heatmap_legend_param = list(
    title_gp      = gpar(fontsize = 15, fontface = "bold"),
    labels_gp     = gpar(fontsize = 10),
    legend_height = unit(3, "cm")
  ),
  border    = TRUE,
  border_gp = gpar(col = "grey40", lwd = 0.5)
)

# ── 2. Función auxiliar (sin subtítulo al pie ahora) ─────────
draw_heatmap <- function(ht_obj) {
  draw(ht_obj,
       padding                = unit(c(9, 9, 9, 9), "mm"),
       heatmap_legend_side    = "right",
       annotation_legend_side = "right")
}

# ── 3. Guardar figura ─────────────────────────────────────────
out_dir <- "~/CESC_Network/5_Network_analisis/5_3_DeltaCon/"

pdf(file.path(out_dir, "DeltaCon_Top100_Heatmap.pdf"),
    width = 16, height = 10, useDingbats = FALSE)
draw_heatmap(ht)
dev.off()

png(file.path(out_dir, "DeltaCon_Top100_Heatmap.png"),
    width = 5000, height = 4000, res = 600, bg = "white")
draw_heatmap(ht)
dev.off()

cat("✅ Figura guardada.\n")
cat("   → DeltaCon_Top100_Heatmap.pdf\n")
cat("   → DeltaCon_Top100_Heatmap.png\n")
