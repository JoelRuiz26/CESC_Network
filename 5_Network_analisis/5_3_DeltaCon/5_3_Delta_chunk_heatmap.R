# ============================================================
# HEATMAP: Top 100 Genes — DeltaCon A7 vs A9
# Versión óptima: identidad algebraica, sin matriz intermedia n×n
# ============================================================

library(Matrix)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ── 0. Cargar ─────────────────────────────────────────────────
out_dir <- "~/CESC_Network/5_Network_analisis/5_3_DeltaCon/"

delta_con     <- readRDS(file.path(out_dir, "DeltaCon_distance_A7_A9.rds"))
delta_con_sim <- readRDS(file.path(out_dir, "DeltaCon_similarity_A7_A9.rds"))
S1            <- readRDS(file.path(out_dir, "S_A7_matrix.rds"))
S2            <- readRDS(file.path(out_dir, "S_A9_matrix.rds"))

# ── 1. Alinear ────────────────────────────────────────────────
genes <- sort(intersect(rownames(S1), rownames(S2)))
S1 <- S1[genes, genes]
S2 <- S2[genes, genes]

# ── 2. RANKING: identidad algebraica ─────────────────────────
# (√a - √b)² = a + b - 2√(ab)
# → gene_score[i] = Σⱼ S1ᵢⱼ + Σⱼ S2ᵢⱼ - 2·Σⱼ√(S1ᵢⱼ·S2ᵢⱼ)
#
# Ventaja: rowSums(S1) y rowSums(S2) son O(n) sin crear matriz nueva.
# Solo el término cruzado necesita chunks, y solo 1 matriz chunk a la vez.
# Pico de RAM = chunk_size × n (ej: 500×17000×8 bytes ≈ 68 MB)
# vs versión anterior: n×n×8 bytes ≈ 2.3 GB para sq_diff

n  <- length(genes)
rs1 <- Matrix::rowSums(S1)   # O(n²) pero sin copiar — muy rápido
rs2 <- Matrix::rowSums(S2)

# Corrección diagonal (excluir autoafinidad del score)
diag_S1 <- diag(S1)
diag_S2 <- diag(S2)

cat("Calculando término cruzado por chunks...\n")
chunk_size  <- 500            # sube si tienes RAM, baja si no
cross_sums  <- numeric(n)
names(cross_sums) <- genes

for (start in seq(1, n, by = chunk_size)) {
  end <- min(start + chunk_size - 1, n)
  idx <- start:end
  
  # Un solo chunk: √S1[chunk,] × √S2[chunk,]  →  rowSums
  chunk_cross        <- rowSums(sqrt(S1[idx, , drop = FALSE]) *
                                  sqrt(S2[idx, , drop = FALSE]))
  cross_sums[idx]    <- chunk_cross
  rm(chunk_cross)
  # ← sin gc() aquí: R libera automáticamente; gc() explícito es lento
}
gc()   # un solo gc() al terminar el loop

# Ensamble del score con corrección diagonal
# (restamos la contribución i=j que no queremos contar)
diag_correction <- diag_S1 + diag_S2 - 2 * sqrt(diag_S1 * diag_S2)
gene_score      <- (rs1 + rs2 - 2 * cross_sums) - diag_correction
rm(rs1, rs2, cross_sums, diag_S1, diag_S2, diag_correction); gc()

# ── 3. Top 100 ────────────────────────────────────────────────
top100 <- names(sort(gene_score, decreasing = TRUE))[1:100]

# ── 4. Submatriz signada (solo 100×100) ──────────────────────
delta_mat <- as.matrix(
  sqrt(S1[top100, top100]) - sqrt(S2[top100, top100])
)
rm(S1, S2); gc()

diag(delta_mat) <- 0

# ── 5. Clipping ───────────────────────────────────────────────
qlim          <- quantile(delta_mat[delta_mat != 0], c(0.01, 0.99), na.rm = TRUE)
delta_clipped <- pmax(pmin(delta_mat, qlim[2]), qlim[1])
rm(delta_mat, qlim); gc()

max_abs <- max(abs(delta_clipped))
col_fun <- colorRamp2(
  c(-max_abs, -max_abs * 0.3, 0, max_abs * 0.3, max_abs),
  c("#2166AC", "#92C5DE", "white", "#F4A582", "#D6604D")
)
rm(max_abs)

# ── 6. Anotación ──────────────────────────────────────────────
score_top    <- gene_score[top100]
score_scaled <- (score_top - min(score_top)) / (max(score_top) - min(score_top))
rm(gene_score, score_top, top100, genes); gc()

row_ann <- rowAnnotation(
  Contribution = anno_barplot(
    score_scaled,
    bar_width  = 0.8,
    gp         = gpar(fill = "#555555", col = NA),
    axis_param = list(gp = gpar(fontsize = 6)),
    width      = unit(1.8, "cm")
  ),
  annotation_name_gp  = gpar(fontsize = 7, fontface = "bold"),
  annotation_name_rot = 90
)
rm(score_scaled); gc()

# ── 7. Títulos ────────────────────────────────────────────────
main_title <- "DeltaCon Affinity Shift: A7 vs A9 — Top 100 Genes"
sub_title  <- sprintf("Distance: %.4f   |   Similarity: %.6f",
                      delta_con, delta_con_sim)
rm(delta_con, delta_con_sim)

save.image("~/CESC_Network/5_Network_analisis/5_3_DeltaCon/5_3_Image_plots.RData")
#load("~/CESC_Network/5_Network_analisis/5_3_DeltaCon/5_3_Image_plots.RData")

# ── 8. Heatmap ────────────────────────────────────────────────
ht <- Heatmap(
  delta_clipped,
  name = "Δ Affinity\n(√S_A7 − √S_A9)",
  col  = col_fun,
  clustering_distance_rows    = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method_rows      = "ward.D2",
  clustering_method_columns   = "ward.D2",
  show_row_dend      = TRUE,
  show_column_dend   = TRUE,
  row_dend_width     = unit(8, "mm"),
  column_dend_height = unit(8, "mm"),
  show_row_names     = TRUE,
  show_column_names  = FALSE,
  row_names_gp       = gpar(fontsize = 4.5),
  right_annotation   = row_ann,
  column_title       = main_title,
  column_title_gp    = gpar(fontsize = 11, fontface = "bold"),
  heatmap_legend_param = list(
    title_gp      = gpar(fontsize = 8, fontface = "bold"),
    labels_gp     = gpar(fontsize = 7),
    legend_height = unit(3, "cm")
  ),
  border    = TRUE,
  border_gp = gpar(col = "grey40", lwd = 0.5)
)
rm(delta_clipped, col_fun, row_ann, main_title); gc()

# ── 9. Guardar ────────────────────────────────────────────────
draw_with_subtitle <- function(ht_obj, subtitle) {
  draw(ht_obj,
       padding                = unit(c(8, 5, 5, 5), "mm"),
       heatmap_legend_side    = "right",
       annotation_legend_side = "right")
  grid.text(subtitle,
            x = unit(0.5, "npc"), y = unit(3, "mm"),
            just = "bottom",
            gp   = gpar(fontsize = 8, col = "grey30", fontface = "italic"))
}

pdf(file.path(out_dir, "DeltaCon_Top100_Heatmap.pdf"),
    width = 10, height = 11, useDingbats = FALSE)
draw_with_subtitle(ht, sub_title)
dev.off()

png(file.path(out_dir, "DeltaCon_Top100_Heatmap.png"),
    width = 3000, height = 3300, res = 300, bg = "white")
draw_with_subtitle(ht, sub_title)
dev.off()

rm(ht, sub_title, out_dir); gc()
cat("✅ Figura guardada.\n")