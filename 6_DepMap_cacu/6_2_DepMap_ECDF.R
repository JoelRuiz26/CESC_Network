############################################################
## 6_2_DepMap_ECDF.R
##
## Contraste de distribuciones de CRISPR dependency score
## en líneas de cáncer cervicouterino (DepMap):
##
##   Grupo 1 — DEGs significativos (p < 0.01, |LFC| > 1)
##   Grupo 2 — Top 2000 genes por degree (hubs de la red)
##
## Pruebas estadísticas (dos grupos, dos lados):
##   - Wilcoxon (Mann-Whitney U)
##   - Kolmogorov-Smirnov
##
## Outputs (en este mismo directorio):
##   6_2_ECDF_A7.pdf / .png
##   6_2_ECDF_A9.pdf / .png
##   6_2_ECDF_combined.pdf / .png
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

TOP_DEGREE_N <- 2000

out_dir   <- "~/CESC_Network/6_DepMap_cacu"
depmap_dir <- "~/CESC_Network/6_DepMap_cacu"
deg_dir    <- "~/CESC_Network/5_Network_analisis/5_4_DEG_integration"

PAL <- c(
  "DEGs significativos"  = "#E63946",
  "Top 2000 Degree"      = "#457B9D"
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

cat("Genes en dep_mean:", nrow(dep_mean), "\n")
cat("Genes en DE_A7:   ", nrow(DE_A7), "\n")
cat("Genes en DE_A9:   ", nrow(DE_A9), "\n")

# ── 2) Construir grupos y mapear a DepMap ────────────────────
build_groups <- function(annot_tbl, label) {

  # Grupo 1: DEGs significativos
  degs <- annot_tbl %>%
    filter(P.Value < 0.01, abs(logFC) > 1) %>%
    mutate(gene_u = toupper(Gene))

  # Grupo 2: Top 2000 por degree
  top_deg <- annot_tbl %>%
    arrange(degree_rank) %>%
    slice_head(n = TOP_DEGREE_N) %>%
    mutate(gene_u = toupper(Gene))

  cat(sprintf("\n[%s] DEGs (p<0.01, |LFC|>1): %d genes\n", label, nrow(degs)))
  cat(sprintf("[%s] Top %d degree:           %d genes\n", label, TOP_DEGREE_N, nrow(top_deg)))

  # Mapeo a DepMap
  g1 <- dep_mean %>%
    filter(gene_u %in% degs$gene_u) %>%
    mutate(grupo = "DEGs significativos")

  g2 <- dep_mean %>%
    filter(gene_u %in% top_deg$gene_u) %>%
    mutate(grupo = "Top 2000 Degree")

  cat(sprintf("[%s] Mapeados a DepMap — DEGs: %d/%d  |  Degree: %d/%d\n",
              label,
              nrow(g1), nrow(degs),
              nrow(g2), nrow(top_deg)))

  list(g1 = g1, g2 = g2,
       n_degs = nrow(degs), n_deg_mapped = nrow(g1),
       n_top  = nrow(top_deg), n_top_mapped = nrow(g2))
}

sets_A7 <- build_groups(DE_A7, "A7")
sets_A9 <- build_groups(DE_A9, "A9")

# ── 3) Pruebas estadísticas ──────────────────────────────────
run_tests <- function(set, label) {
  x1 <- set$g1$mean_dep   # DEGs
  x2 <- set$g2$mean_dep   # Top degree

  wt <- wilcox.test(x1, x2, exact = FALSE)
  ks <- ks.test(x1, x2)

  cat(sprintf(
    "\n[%s] DEGs  n=%d  median=%.4f\n       Degree n=%d  median=%.4f\n",
    label, length(x1), median(x1), length(x2), median(x2)
  ))
  cat(sprintf("  Wilcoxon W=%.0f  p=%.4g\n", wt$statistic, wt$p.value))
  cat(sprintf("  KS D=%.4f  p=%.4g\n",        ks$statistic, ks$p.value))

  list(wilcox = wt, ks = ks)
}

stats_A7 <- run_tests(sets_A7, "A7")
stats_A9 <- run_tests(sets_A9, "A9")

# ── 4) ECDF plots ────────────────────────────────────────────
fmt_p <- function(p) {
  if (p < 0.001) "< 0.001" else sprintf("= %.4f", p)
}

make_ecdf <- function(set, stats, label) {
  df <- bind_rows(set$g1, set$g2) %>%
    mutate(grupo = factor(grupo, levels = c("DEGs significativos", "Top 2000 Degree")))

  ann <- sprintf(
    "Wilcoxon p %s\nKS D = %.3f, p %s",
    fmt_p(stats$wilcox$p.value),
    stats$ks$statistic,
    fmt_p(stats$ks$p.value)
  )

  sub_txt <- sprintf(
    "DEGs: %d mapeados/%d totales  |  Top 2000 Degree: %d mapeados/%d totales",
    set$n_deg_mapped, set$n_degs,
    set$n_top_mapped, set$n_top
  )

  ggplot(df, aes(x = mean_dep, colour = grupo)) +
    stat_ecdf(linewidth = 1, alpha = 0.9) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    annotate("text",
             x = -2.5, y = 0.90,
             label = ann,
             hjust = 0, size = 3.2,
             colour = "grey25", family = "mono") +
    scale_colour_manual(values = PAL) +
    labs(
      title    = sprintf("ECDF — %s vs SNT", label),
      subtitle = sub_txt,
      x        = "Score de dependencia CRISPR promedio\n(más negativo = mayor esencialidad)",
      y        = "Probabilidad acumulada",
      colour   = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position  = c(0.75, 0.20),
      legend.text      = element_text(size = 9),
      legend.key.width = unit(1.2, "cm"),
      plot.title       = element_text(face = "bold"),
      plot.subtitle    = element_text(colour = "grey40", size = 8)
    )
}

ecdf_A7 <- make_ecdf(sets_A7, stats_A7, "A7")
ecdf_A9 <- make_ecdf(sets_A9, stats_A9, "A9")

# ── 5) Guardar ───────────────────────────────────────────────
save_plot <- function(plot, name, w = 7, h = 5) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")),
         plot = plot, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(out_dir, paste0(name, ".png")),
         plot = plot, width = w, height = h, dpi = 300)
}

save_plot(ecdf_A7, "6_2_ECDF_A7")
save_plot(ecdf_A9, "6_2_ECDF_A9")

combined <- ecdf_A7 + ecdf_A9 +
  plot_annotation(
    title = "ECDF — CRISPR Dependency Score | Líneas Cervicales (DepMap)",
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13))
  )

save_plot(combined, "6_2_ECDF_combined", w = 14, h = 5)

cat("\nPlots guardados en:", out_dir, "\n")
cat("  6_2_ECDF_A7.pdf / .png\n")
cat("  6_2_ECDF_A9.pdf / .png\n")
cat("  6_2_ECDF_combined.pdf / .png\n")
cat("✓ Análisis completo.\n")
