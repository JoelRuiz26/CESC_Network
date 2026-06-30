############################################################
## 6_2_DepMap_ECDF.R
##
## Contraste de distribuciones de CRISPR dependency score
## en líneas de cáncer cervicouterino (DepMap):
##
##   Grupo 1 — DEGs significativos (p < 0.01, |LFC| > 1.5)
##             → N_mapped genes tras unir con DepMap
##
##   Grupo 2 — Top connected genes
##             → Se seleccionan los más conectados que
##               TENGAN datos en DepMap, hasta igualar N_mapped
##
## Pruebas estadísticas:
##   - Wilcoxon (Mann-Whitney U)
##   - Kolmogorov-Smirnov
##   - Cliff's delta (effect size)
##
## Output: 6_2_ECDF_combined.pdf / .png
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

out_dir    <- "~/CESC_Network/6_DepMap_cacu"
depmap_dir <- "~/CESC_Network/6_DepMap_cacu"
deg_dir    <- "~/CESC_Network/5_Network_analisis/5_4_DEG_integration"

PAL <- c(
  "Significant DEGs" = "#E63946",
  "Top connected"    = "#457B9D"
)

# ── 1) Cargar tablas ─────────────────────────────────────────
dep_mean <- read_tsv(
  file.path(depmap_dir, "6_1_cervix_dep_mean.tsv"),
  show_col_types = FALSE
)

DE_A7 <- read_tsv(
  file.path(deg_dir, "5_4_1_DEG_A7_Degree.tsv"),
  show_col_types = FALSE
)
DE_A9 <- read_tsv(
  file.path(deg_dir, "5_4_1_DEG_A9_Degree.tsv"),
  show_col_types = FALSE
)

# ── 2) Construir grupos ──────────────────────────────────────
# Top connected: genes ordenados por degree_rank mapeados a DepMap,
# se toman tantos como DEGs mapeados (N = nDEGs).

build_groups <- function(annot_tbl, label) {

  degs <- annot_tbl %>%
    filter(P.Value < 0.01, abs(logFC) > 1.5) %>%
    mutate(gene_u = toupper(Gene))

  g1 <- dep_mean %>%
    filter(gene_u %in% degs$gene_u) %>%
    mutate(grupo = "Significant DEGs")

  N_target <- nrow(g1)

  g2 <- annot_tbl %>%
    arrange(degree_rank) %>%
    mutate(gene_u = toupper(Gene)) %>%
    inner_join(dep_mean, by = "gene_u") %>%
    slice_head(n = N_target) %>%
    select(names(dep_mean)) %>%
    mutate(grupo = "Top connected")

  cat(sprintf("[%s] DEGs: %d totales → %d mapeados\n",
              label, nrow(degs), nrow(g1)))
  cat(sprintf("[%s] Top connected (N = nDEGs, mapeados en DepMap): %d\n",
              label, nrow(g2)))

  list(g1 = g1, g2 = g2,
       n_degs_total = nrow(degs),
       n_mapped     = N_target)
}

sets_A7 <- build_groups(DE_A7, "A7")
sets_A9 <- build_groups(DE_A9, "A9")

# ── 3) Pruebas estadísticas ──────────────────────────────────
sig_stars <- function(p) {
  if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
}

fmt_p <- function(p) {
  if (p < 0.001) "< 0.001" else sprintf("= %.4f", p)
}

cliff_delta <- function(x, y) {
  U  <- as.numeric(wilcox.test(x, y, exact = FALSE)$statistic)
  (2 * U - length(x) * length(y)) / (length(x) * length(y))
}

cliff_magnitude <- function(d) {
  ad <- abs(d)
  if (ad < 0.147) "negligible" else if (ad < 0.330) "small" else if (ad < 0.474) "medium" else "large"
}

run_tests <- function(set, label) {
  x1 <- set$g1$mean_dep
  x2 <- set$g2$mean_dep

  wt <- wilcox.test(x1, x2, exact = FALSE)
  ks <- suppressWarnings(ks.test(x1, x2))
  cd <- cliff_delta(x1, x2)

  cat(sprintf(
    "\n[%s] DEGs n=%d  median=%.4f  |  Top connected n=%d  median=%.4f\n",
    label, length(x1), median(x1), length(x2), median(x2)
  ))
  cat(sprintf("  Wilcoxon W=%.0f  p=%.4g  %s\n",
              wt$statistic, wt$p.value, sig_stars(wt$p.value)))
  cat(sprintf("  KS D=%.4f  p=%.4g\n", ks$statistic, ks$p.value))
  cat(sprintf("  Cliff's delta = %.4f  (%s)\n", cd, cliff_magnitude(cd)))

  list(wilcox = wt, ks = ks, cliff_d = cd)
}

stats_A7 <- run_tests(sets_A7, "A7")
stats_A9 <- run_tests(sets_A9, "A9")

# ── 4) ECDF plots ────────────────────────────────────────────
make_ecdf <- function(set, stats, label) {
  df <- bind_rows(set$g1, set$g2) %>%
    mutate(grupo = factor(grupo, levels = c("Significant DEGs", "Top connected")))

  ann_p <- sprintf("Wilcoxon p %s  %s",
                   fmt_p(stats$wilcox$p.value),
                   sig_stars(stats$wilcox$p.value))

  ggplot(df, aes(x = mean_dep, colour = grupo)) +
    stat_ecdf(linewidth = 1, alpha = 0.9) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    annotate("text",
             x = -Inf, y = 0.12,
             label = ann_p,
             hjust = -0.05, vjust = 0,
             size = 3.5, colour = "grey25") +
    scale_colour_manual(values = PAL) +
    labs(
      title  = sprintf("HPV-%s", label),
      x      = "Mean CRISPR Dependency Score",
      y      = "Cumulative Probability",
      colour = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position  = c(0.28, 0.88),
      legend.text      = element_text(size = 9),
      legend.key.width = unit(1.2, "cm"),
      plot.title       = element_text(face = "bold", size = 11)
    )
}

ecdf_A7 <- make_ecdf(sets_A7, stats_A7, "A7")
ecdf_A9 <- make_ecdf(sets_A9, stats_A9, "A9")

# ── 5) Guardar ───────────────────────────────────────────────
combined <- ecdf_A7 + ecdf_A9 +
  plot_annotation(
    title = "CRISPR Gene Dependency in Cervical Cancer Cell Lines (DepMap)",
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13))
  )

ggsave(file.path(out_dir, "6_2_ECDF_combined.pdf"),
       plot = combined, width = 14, height = 5, device = cairo_pdf)
ggsave(file.path(out_dir, "6_2_ECDF_combined.png"),
       plot = combined, width = 14, height = 5, dpi = 300)

cat("\nPlot guardado: 6_2_ECDF_combined.pdf / .png\n")
cat("✓ Análisis completo.\n")
