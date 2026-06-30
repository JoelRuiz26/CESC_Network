################################################################################
# Louvain Module Comparison: HPV-A7 vs HPV-A9
# Plot B — Jaccard similarity heatmap (todos los módulos, 32 × 41)
################################################################################

suppressPackageStartupMessages({
  library(igraph)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(grid)
  library(scales)
})

DIR <- "~/CESC_Network/7_Modularity/5_1_Modules"

# ── Load graphs ────────────────────────────────────────────────────────────────
g7 <- readRDS(file.path(DIR, "5_1_1_Louvain_A7_graph_max.rds"))
g9 <- readRDS(file.path(DIR, "5_1_1_Louvain_A9_graph_max.rds"))

# ── Gene–module tables ─────────────────────────────────────────────────────────
tbl7 <- tibble(gene  = V(g7)$name,
               mod7  = V(g7)$Community_name,
               size7 = as.integer(V(g7)$Community_size))

tbl9 <- tibble(gene  = V(g9)$name,
               mod9  = V(g9)$Community_name,
               size9 = as.integer(V(g9)$Community_size))

# ── Module sizes (ordenados descendente) ───────────────────────────────────────
sizes7 <- tbl7 %>% distinct(mod7, size7) %>% arrange(desc(size7))
sizes9 <- tbl9 %>% distinct(mod9, size9) %>% arrange(desc(size9))

top7 <- sizes7$mod7   # 32 módulos A7
top9 <- sizes9$mod9   # 41 módulos A9

# ── Jaccard matrix (32 × 41) ──────────────────────────────────────────────────
cat("[1] Calculando matriz Jaccard 32 × 41...\n")

genes7 <- lapply(top7, function(m) tbl7 %>% filter(mod7 == m) %>% pull(gene))
genes9 <- lapply(top9, function(m) tbl9 %>% filter(mod9 == m) %>% pull(gene))
names(genes7) <- top7; names(genes9) <- top9

jaccard_mat <- outer(seq_along(top7), seq_along(top9), Vectorize(function(i, j) {
  a <- genes7[[i]]; b <- genes9[[j]]
  length(intersect(a, b)) / length(union(a, b))
}))
rownames(jaccard_mat) <- top7
colnames(jaccard_mat) <- top9

cat("  Dimensiones:", dim(jaccard_mat), "\n")
cat("  Rango:", round(range(jaccard_mat), 3), "\n")

# ── Annotations ───────────────────────────────────────────────────────────────
ann_row <- data.frame(Size_A7 = sizes7$size7, row.names = sizes7$mod7)
ann_col <- data.frame(Size_A9 = sizes9$size9, row.names = sizes9$mod9)

# ── Helper: recursively make all text grobs bold ───────────────────────────────
make_bold <- function(g) {
  if (inherits(g, "text")) {
    g$gp <- modifyList(if (is.null(g$gp)) gpar() else g$gp,
                       gpar(fontface = "bold"))
  }
  if (!is.null(g$children)) g$children <- lapply(g$children, make_bold)
  if (!is.null(g$grobs))    g$grobs    <- lapply(g$grobs,    make_bold)
  g
}

# ══════════════════════════════════════════════════════════════════════════════
# PLOT B — Jaccard heatmap
# ══════════════════════════════════════════════════════════════════════════════
cat("[2] Generando heatmap...\n")

p <- pheatmap(
  jaccard_mat,
  color           = colorRampPalette(c("white", "#FFF3CD", "#E8602D", "#8B0000"))(100),
  cluster_rows    = TRUE,
  cluster_cols    = TRUE,
  annotation_row  = ann_row,
  annotation_col  = ann_col,
  fontsize        = 22,
  fontsize_row    = 22,
  fontsize_col    = 22,
  cellwidth       = 52,
  cellheight      = 52,
  main            = "Jaccard Similarity between Louvain Modules",
  display_numbers = TRUE,
  number_format   = "%.2f",
  number_color    = "grey20",
  fontsize_number = 16,
  angle_col       = 45,
  margins         = c(14, 14),
  silent          = TRUE
)

# Aplicar negritas a todos los elementos de texto
gt <- make_bold(p$gtable)

png(file.path(DIR, "5_1_3B_jaccard_heatmap.png"),
    width = 11000, height = 8000, res = 300)
grid.draw(gt)
dev.off()

# ── Summary ────────────────────────────────────────────────────────────────────
shared_genes <- intersect(V(g7)$name, V(g9)$name)
jaccard_pair <- function(a, b) length(intersect(a, b)) / length(union(a, b))
best_match <- lapply(top7, function(m) {
  a <- genes7[[m]]
  j <- sapply(genes9, function(b) jaccard_pair(a, b))
  tibble(mod_A7 = m, best_A9 = names(which.max(j)),
         max_jaccard = max(j), size_A7 = length(a))
}) %>% bind_rows() %>% arrange(desc(max_jaccard))

cat("\n══ COMPARATIVE SUMMARY ══════════════════════════════════\n")
cat(sprintf("  A7 modules: %d  |  A9 modules: %d\n", nrow(sizes7), nrow(sizes9)))
cat(sprintf("  Genes shared between networks: %d\n", length(shared_genes)))
cat(sprintf("  Modules with Jaccard > 0.25 (good match): %d / %d\n",
            sum(best_match$max_jaccard > 0.25), nrow(best_match)))
cat(sprintf("  Modules with Jaccard < 0.05 (no clear match): %d / %d\n",
            sum(best_match$max_jaccard < 0.05), nrow(best_match)))
cat("════════════════════════════════════════════════════════\n")
cat("[OK] Saved: 5_1_3B_jaccard_heatmap.png\n")
