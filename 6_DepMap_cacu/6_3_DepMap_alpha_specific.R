############################################################
## 6_3_DepMap_alpha_specific.R
##
## Alpha-specific and cross-alpha CRISPR dependency analysis
##
## PART 1 — Alpha-specific (matched):
##   HPV-A7 genes  ×  A7 cell lines
##   HPV-A9 genes  ×  A9 cell lines
##
## PART 2 — Cross-alpha:
##   HPV-A7 genes  ×  A9 cell lines
##   HPV-A9 genes  ×  A7 cell lines
##
## Gene sets per comparison:
##   - Significant DEGs (p < 0.01, |LFC| > 1)
##   - Hub genes (top 10% most connected in the network)
##
## Requires outputs from 6_1_DepMap_download.R
##
## Statistical test: Wilcoxon (Mann-Whitney U) + Cliff's delta
##
## Output: 6_3_ECDF_alpha_specific.pdf / .png
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
  "Hub genes"        = "#457B9D"
)

# ── 1) Cargar tablas generadas por 6_1 ───────────────────────
dep_mean_A7 <- read_tsv(
  file.path(depmap_dir, "6_1_cervix_dep_mean_A7.tsv"),
  show_col_types = FALSE
)
dep_mean_A9 <- read_tsv(
  file.path(depmap_dir, "6_1_cervix_dep_mean_A9.tsv"),
  show_col_types = FALSE
)

crispr_lines <- read_tsv(
  file.path(depmap_dir, "6_1_crispr_lines.tsv"),
  show_col_types = FALSE
)

lines_A7 <- crispr_lines %>% filter(hpv_group == "A7") %>% pull(cell_line)
lines_A9 <- crispr_lines %>% filter(hpv_group == "A9") %>% pull(cell_line)

cat("A7 lines con CRISPR (n =", length(lines_A7), "):",
    paste(lines_A7, collapse = ", "), "\n")
cat("A9 lines con CRISPR (n =", length(lines_A9), "):",
    paste(lines_A9, collapse = ", "), "\n")

# ── 2) Load DEG/degree tables ─────────────────────────────────
DE_A7 <- read_tsv(
  file.path(deg_dir, "5_4_1_DEG_A7_Degree.tsv"),
  show_col_types = FALSE
)
DE_A9 <- read_tsv(
  file.path(deg_dir, "5_4_1_DEG_A9_Degree.tsv"),
  show_col_types = FALSE
)

# ── 3) Helpers estadísticos ───────────────────────────────────
sig_stars <- function(p) {
  if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
}

fmt_p <- function(p) {
  if (p < 0.001) "< 0.001" else sprintf("= %.4f", p)
}

cliff_delta <- function(x, y) {
  U <- as.numeric(wilcox.test(x, y, exact = FALSE)$statistic)
  (2 * U - length(x) * length(y)) / (length(x) * length(y))
}

cliff_magnitude <- function(d) {
  ad <- abs(d)
  if (ad < 0.147) "negligible" else if (ad < 0.330) "small" else if (ad < 0.474) "medium" else "large"
}

fmt_lines <- function(cell_lines) {
  paste0(sub("_.*", "", cell_lines), collapse = ", ")
}

# ── 4) Build gene groups ──────────────────────────────────────
# Hub genes = top 10% most connected genes in the network,
# then filtered to those with DepMap CRISPR data.

build_groups <- function(annot_tbl, dep_mean_tbl) {

  degs <- annot_tbl %>%
    filter(P.Value < 0.01, abs(logFC) > 1) %>%
    mutate(gene_u = toupper(Gene))

  g1 <- dep_mean_tbl %>%
    filter(gene_u %in% degs$gene_u) %>%
    mutate(grupo = "Significant DEGs")

  n_hub <- ceiling(0.10 * nrow(annot_tbl))

  g2 <- annot_tbl %>%
    arrange(degree_rank) %>%
    slice_head(n = n_hub) %>%
    mutate(gene_u = toupper(Gene)) %>%
    inner_join(dep_mean_tbl, by = "gene_u") %>%
    select(gene_u, mean_dep) %>%
    mutate(grupo = "Hub genes")

  list(g1 = g1, g2 = g2,
       n_degs = nrow(degs), n_mapped = nrow(g1),
       n_hub = nrow(g2))
}

# ── 5) ECDF plot ──────────────────────────────────────────────
make_ecdf <- function(set, wt, cd, title_txt, lines_txt) {

  ann <- sprintf("Wilcoxon p %s  %s", fmt_p(wt$p.value), sig_stars(wt$p.value))

  df <- bind_rows(set$g1, set$g2) %>%
    mutate(grupo = factor(grupo, levels = c("Significant DEGs", "Hub genes")))

  ggplot(df, aes(x = mean_dep, colour = grupo)) +
    stat_ecdf(linewidth = 1, alpha = 0.9) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    annotate("text",
             x = -Inf, y = 0.12,
             label = ann,
             hjust = -0.05, vjust = 0,
             size = 3.5, colour = "grey25") +
    scale_colour_manual(values = PAL) +
    labs(
      title    = title_txt,
      subtitle = paste0("Cell lines: ", lines_txt),
      x        = "Mean CRISPR Dependency Score",
      y        = "Cumulative Probability",
      colour   = NULL
    ) +
    theme_classic(base_size = 13) +
    theme(
      legend.position  = c(0.30, 0.88),
      legend.text      = element_text(size = 10),
      legend.key.width = unit(1.0, "cm"),
      plot.title       = element_text(face = "bold", size = 13),
      plot.subtitle    = element_text(colour = "grey40", size = 9)
    )
}

# ══════════════════════════════════════════════════════════════
## PART 1 — Alpha-specific (matched)
# ══════════════════════════════════════════════════════════════
cat("\n── PART 1: Alpha-specific ──────────────────────────────\n")

sets_A7_A7 <- build_groups(DE_A7, dep_mean_A7)
wt_A7_A7   <- wilcox.test(sets_A7_A7$g1$mean_dep, sets_A7_A7$g2$mean_dep, exact = FALSE)
cd_A7_A7   <- cliff_delta(sets_A7_A7$g1$mean_dep, sets_A7_A7$g2$mean_dep)
cat(sprintf("[HPV-A7 × A7] DEGs: %d  |  Hub genes: %d\n",
            sets_A7_A7$n_mapped, sets_A7_A7$n_hub))
cat(sprintf("  Wilcoxon W=%.0f  p=%.4g  %s\n",
            wt_A7_A7$statistic, wt_A7_A7$p.value, sig_stars(wt_A7_A7$p.value)))
cat(sprintf("  Cliff's delta = %.4f  (%s)\n", cd_A7_A7, cliff_magnitude(cd_A7_A7)))

sets_A9_A9 <- build_groups(DE_A9, dep_mean_A9)
wt_A9_A9   <- wilcox.test(sets_A9_A9$g1$mean_dep, sets_A9_A9$g2$mean_dep, exact = FALSE)
cd_A9_A9   <- cliff_delta(sets_A9_A9$g1$mean_dep, sets_A9_A9$g2$mean_dep)
cat(sprintf("[HPV-A9 × A9] DEGs: %d  |  Hub genes: %d\n",
            sets_A9_A9$n_mapped, sets_A9_A9$n_hub))
cat(sprintf("  Wilcoxon W=%.0f  p=%.4g  %s\n",
            wt_A9_A9$statistic, wt_A9_A9$p.value, sig_stars(wt_A9_A9$p.value)))
cat(sprintf("  Cliff's delta = %.4f  (%s)\n", cd_A9_A9, cliff_magnitude(cd_A9_A9)))

p1 <- make_ecdf(sets_A7_A7, wt_A7_A7, cd_A7_A7, "HPV-A7", fmt_lines(lines_A7))
p2 <- make_ecdf(sets_A9_A9, wt_A9_A9, cd_A9_A9, "HPV-A9", fmt_lines(lines_A9))

# ══════════════════════════════════════════════════════════════
## PART 2 — Cross-alpha
# ══════════════════════════════════════════════════════════════
cat("\n── PART 2: Cross-alpha ─────────────────────────────────\n")

sets_A7_A9 <- build_groups(DE_A7, dep_mean_A9)
wt_A7_A9   <- wilcox.test(sets_A7_A9$g1$mean_dep, sets_A7_A9$g2$mean_dep, exact = FALSE)
cd_A7_A9   <- cliff_delta(sets_A7_A9$g1$mean_dep, sets_A7_A9$g2$mean_dep)
cat(sprintf("[HPV-A7 genes × A9 lines] DEGs: %d  |  Hub genes: %d\n",
            sets_A7_A9$n_mapped, sets_A7_A9$n_hub))
cat(sprintf("  Wilcoxon W=%.0f  p=%.4g  %s\n",
            wt_A7_A9$statistic, wt_A7_A9$p.value, sig_stars(wt_A7_A9$p.value)))
cat(sprintf("  Cliff's delta = %.4f  (%s)\n", cd_A7_A9, cliff_magnitude(cd_A7_A9)))

sets_A9_A7 <- build_groups(DE_A9, dep_mean_A7)
wt_A9_A7   <- wilcox.test(sets_A9_A7$g1$mean_dep, sets_A9_A7$g2$mean_dep, exact = FALSE)
cd_A9_A7   <- cliff_delta(sets_A9_A7$g1$mean_dep, sets_A9_A7$g2$mean_dep)
cat(sprintf("[HPV-A9 genes × A7 lines] DEGs: %d  |  Hub genes: %d\n",
            sets_A9_A7$n_mapped, sets_A9_A7$n_hub))
cat(sprintf("  Wilcoxon W=%.0f  p=%.4g  %s\n",
            wt_A9_A7$statistic, wt_A9_A7$p.value, sig_stars(wt_A9_A7$p.value)))
cat(sprintf("  Cliff's delta = %.4f  (%s)\n", cd_A9_A7, cliff_magnitude(cd_A9_A7)))

p3 <- make_ecdf(sets_A7_A9, wt_A7_A9, cd_A7_A9,
                "HPV-A7 genes  ×  A9 cell lines", fmt_lines(lines_A9))
p4 <- make_ecdf(sets_A9_A7, wt_A9_A7, cd_A9_A7,
                "HPV-A9 genes  ×  A7 cell lines", fmt_lines(lines_A7))

# ── 6) Output A: alpha-specific only (p1 + p2) ───────────────
specific_only <- p1 + p2 +
  plot_annotation(
    title = "Alpha-Specific CRISPR Dependency | Cervical Cancer Cell Lines (DepMap)",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12)
    )
  )

ggsave(file.path(out_dir, "6_3_ECDF_matched.pdf"),
       plot = specific_only, width = 14, height = 5, device = cairo_pdf)
ggsave(file.path(out_dir, "6_3_ECDF_matched.png"),
       plot = specific_only, width = 14, height = 5, dpi = 300)

# ── 7) Output B: full 2×2 (matched + cross) ──────────────────
combined <- (p1 + p2) / (p3 + p4) +
  plot_annotation(
    title = "Alpha-Specific and Cross-Alpha CRISPR Dependency | Cervical Cancer Cell Lines (DepMap)",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12)
    )
  )

ggsave(file.path(out_dir, "6_3_ECDF_alpha_specific.pdf"),
       plot = combined, width = 14, height = 10, device = cairo_pdf)
ggsave(file.path(out_dir, "6_3_ECDF_alpha_specific.png"),
       plot = combined, width = 14, height = 10, dpi = 300)

cat("\nPlots saved:\n")
cat("  6_3_ECDF_matched.pdf / .png\n")
cat("  6_3_ECDF_alpha_specific.pdf / .png\n")
cat("✓ Analysis complete.\n")

# ══════════════════════════════════════════════════════════════
## PART 3 — Intersection plots (alpha-specific only)
##
## 3 curves per panel:
##   (1) Significant DEGs    — full DEG set (red)
##   (2) Hub genes           — top 10% most connected (blue)
##   (3) DEGs ∩ Hub genes    — DEGs that are also hubs (teal)
##
## Wilcoxon + Cliff's delta: intersection vs full DEGs
# ══════════════════════════════════════════════════════════════
cat("\n── PART 3: Intersection plots ──────────────────────────\n")

PAL3 <- c(
  "Significant DEGs"  = "#E63946",
  "Hub genes"         = "#457B9D",
  "DEGs ∩ Hub genes"  = "#2A9D8F"
)

build_inter_groups <- function(annot_tbl, dep_mean_tbl) {

  degs <- annot_tbl %>%
    filter(P.Value < 0.01, abs(logFC) > 1) %>%
    mutate(gene_u = toupper(Gene))

  g_degs <- dep_mean_tbl %>%
    filter(gene_u %in% degs$gene_u) %>%
    mutate(grupo = "Significant DEGs")

  n_hub <- ceiling(0.10 * nrow(annot_tbl))

  hub_tbl <- annot_tbl %>%
    arrange(degree_rank) %>%
    slice_head(n = n_hub) %>%
    mutate(gene_u = toupper(Gene)) %>%
    inner_join(dep_mean_tbl, by = "gene_u")

  g_hub <- hub_tbl %>%
    select(gene_u, mean_dep) %>%
    mutate(grupo = "Hub genes")

  hub_genes <- hub_tbl %>% pull(gene_u)

  inter <- intersect(degs$gene_u, hub_genes)
  g_inter <- dep_mean_tbl %>%
    filter(gene_u %in% inter) %>%
    mutate(grupo = "DEGs ∩ Hub genes")

  list(g_degs  = g_degs,
       g_hub   = g_hub,
       g_inter = g_inter,
       n_degs  = nrow(degs),
       n_hub   = nrow(g_hub),
       n_inter = nrow(g_inter))
}

make_ecdf_inter <- function(set, title_txt, lines_txt) {

  wt <- wilcox.test(set$g_inter$mean_dep, set$g_degs$mean_dep, exact = FALSE)
  cd <- cliff_delta(set$g_inter$mean_dep, set$g_degs$mean_dep)

  cat(sprintf("  DEGs ∩ Hub genes: n = %d\n", set$n_inter))
  cat(sprintf("  Wilcoxon W=%.0f  p=%.4g  %s\n",
              wt$statistic, wt$p.value, sig_stars(wt$p.value)))
  cat(sprintf("  Cliff's delta = %.4f  (%s)\n", cd, cliff_magnitude(cd)))

  ann <- sprintf("DEGs ∩ Hub genes: Wilcoxon p %s  %s",
                 fmt_p(wt$p.value), sig_stars(wt$p.value))

  df <- bind_rows(set$g_degs, set$g_hub, set$g_inter) %>%
    mutate(grupo = factor(grupo, levels = c(
      "Significant DEGs", "Hub genes", "DEGs ∩ Hub genes"
    )))

  ggplot(df, aes(x = mean_dep, colour = grupo)) +
    stat_ecdf(linewidth = 1, alpha = 0.9) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    annotate("text",
             x = -Inf, y = 0.12,
             label = ann,
             hjust = -0.05, vjust = 0,
             size = 3.3, colour = "grey25") +
    scale_colour_manual(values = PAL3) +
    labs(
      title    = title_txt,
      subtitle = paste0("Cell lines: ", lines_txt),
      x        = "Mean CRISPR Dependency Score",
      y        = "Cumulative Probability",
      colour   = NULL
    ) +
    theme_classic(base_size = 13) +
    theme(
      legend.position  = c(0.32, 0.88),
      legend.text      = element_text(size = 10),
      legend.key.width = unit(1.0, "cm"),
      plot.title       = element_text(face = "bold", size = 13),
      plot.subtitle    = element_text(colour = "grey40", size = 9)
    )
}

cat("\n[HPV-A7 — intersection]\n")
inter_A7 <- build_inter_groups(DE_A7, dep_mean_A7)
pi1 <- make_ecdf_inter(inter_A7, "HPV-A7", fmt_lines(lines_A7))

cat("\n[HPV-A9 — intersection]\n")
inter_A9 <- build_inter_groups(DE_A9, dep_mean_A9)
pi2 <- make_ecdf_inter(inter_A9, "HPV-A9", fmt_lines(lines_A9))

inter_combined <- pi1 + pi2 +
  plot_annotation(
    title = "DEG–Hub Intersection CRISPR Dependency | Cervical Cancer Cell Lines (DepMap)",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12)
    )
  )

ggsave(file.path(out_dir, "6_3_ECDF_intersection.pdf"),
       plot = inter_combined, width = 14, height = 5, device = cairo_pdf)
ggsave(file.path(out_dir, "6_3_ECDF_intersection.png"),
       plot = inter_combined, width = 14, height = 5, dpi = 300)

cat("\n  6_3_ECDF_intersection.pdf / .png saved.\n")
