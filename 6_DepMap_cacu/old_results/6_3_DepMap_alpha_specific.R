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
##   - DE genes (p < 0.01, |LFC| > 1.5)
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
  "DE genes"  = "#E67E22",   # naranja
  "Hub genes" = "#00838F"    # agua marina oscuro
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
    filter(P.Value < 0.01, abs(logFC) > 1.5) %>%
    mutate(gene_u = toupper(Gene))

  g1 <- dep_mean_tbl %>%
    filter(gene_u %in% degs$gene_u) %>%
    mutate(grupo = "DE genes")

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
    mutate(grupo = factor(grupo, levels = c("DE genes", "Hub genes")))

  ggplot(df, aes(x = mean_dep, colour = grupo)) +
    stat_ecdf(linewidth = 1, alpha = 0.9) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    annotate("text",
             x = -Inf, y = 0.22,
             label = ann,
             hjust = -0.05, vjust = 0,
             size = 4.5, colour = "grey25") +
    scale_colour_manual(values = PAL) +
    labs(
      title    = title_txt,
      subtitle = paste0("Cell lines: ", lines_txt),
      x        = "Mean CRISPR Dependency Score",
      y        = "Cumulative Probability",
      colour   = NULL
    ) +
    theme_classic(base_size = 15) +
    theme(
      legend.position  = c(0.30, 0.88),
      legend.text      = element_text(size = 12),
      legend.key.width = unit(1.0, "cm"),
      plot.title       = element_text(face = "bold", size = 15),
      plot.subtitle    = element_text(colour = "grey40", size = 11)
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
      plot.title = element_text(face = "bold", hjust = 0.5, size = 15)
    )
  )

ggsave(file.path(out_dir, "6_3_ECDF_matched.pdf"),
       plot = specific_only, width = 10, height = 5.2, device = cairo_pdf)
ggsave(file.path(out_dir, "6_3_ECDF_matched.png"),
       plot = specific_only, width = 10, height = 5.2, dpi = 300)

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
##   (1) DE genes    — full DEG set (red)
##   (2) Hub genes           — top 10% most connected (blue)
##   (3) DEGs ∩ Hub genes    — DEGs that are also hubs (teal)
##
## Wilcoxon + Cliff's delta: intersection vs full DEGs
# ══════════════════════════════════════════════════════════════
cat("\n── PART 3: Intersection plots ──────────────────────────\n")

PAL3 <- c(
  "DE genes"          = "#E67E22",   # naranja
  "Hub genes"         = "#C2185B",   # fucsia
  "DEGs ∩ Hub genes"  = "#D4AC0D"   # dorado
)

build_inter_groups <- function(annot_tbl, dep_mean_tbl) {

  degs <- annot_tbl %>%
    filter(P.Value < 0.01, abs(logFC) > 1.5) %>%
    mutate(gene_u = toupper(Gene))

  g_degs <- dep_mean_tbl %>%
    filter(gene_u %in% degs$gene_u) %>%
    mutate(grupo = "DE genes")

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
      "DE genes", "Hub genes", "DEGs ∩ Hub genes"
    )))

  ggplot(df, aes(x = mean_dep, colour = grupo)) +
    stat_ecdf(linewidth = 1, alpha = 0.9) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    annotate("text",
             x = -Inf, y = 0.20,
             label = ann,
             hjust = -0.05, vjust = 0,
             size = 4.3, colour = "grey25") +
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

# ══════════════════════════════════════════════════════════════
## PART 4 — Violin plot: DEGs vs Hub genes (alpha-specific matched)
# ══════════════════════════════════════════════════════════════
cat("\n── PART 4: Violin plot (matched) ───────────────────────\n")

suppressPackageStartupMessages(library(ggsignif))

violin_df <- bind_rows(
  bind_rows(sets_A7_A7$g1, sets_A7_A7$g2) %>% mutate(network = "HPV-A7"),
  bind_rows(sets_A9_A9$g1, sets_A9_A9$g2) %>% mutate(network = "HPV-A9")
) %>%
  mutate(
    grupo   = factor(grupo,   levels = c("DE genes", "Hub genes")),
    network = factor(network, levels = c("HPV-A7", "HPV-A9"))
  )

# y position for bracket: just above the 99th percentile of each facet
y_pos <- violin_df %>%
  group_by(network) %>%
  summarise(y = quantile(mean_dep, 0.995) + 0.03, .groups = "drop")

# ggsignif con facetas requiere manual = TRUE y un dataframe por faceta
signif_df <- data.frame(
  network     = factor(c("HPV-A7", "HPV-A9"), levels = c("HPV-A7", "HPV-A9")),
  xmin        = c(1, 1),
  xmax        = c(2, 2),
  annotations = c(sig_stars(wt_A7_A7$p.value), sig_stars(wt_A9_A9$p.value)),
  y_position  = y_pos$y
)

violin_plot <- ggplot(violin_df,
                      aes(x = grupo, y = mean_dep,
                          fill = grupo, colour = grupo)) +
  geom_violin(trim = FALSE, alpha = 0.55, linewidth = 0.4) +
  geom_boxplot(width = 0.10, fill = "white", colour = "grey30",
               outlier.shape = NA, linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  geom_signif(
    data         = signif_df,
    aes(xmin = xmin, xmax = xmax,
        annotations = annotations, y_position = y_position),
    manual       = TRUE,
    inherit.aes  = FALSE,
    tip_length   = 0.01,
    textsize     = 5,
    colour       = "grey20",
    vjust        = 0.4
  ) +
  facet_wrap(~network, scales = "fixed") +
  scale_fill_manual(values   = PAL) +
  scale_colour_manual(values = PAL) +
  labs(
    title = "CRISPR Dependency: DEGs vs Hub genes | Alpha-specific matched",
    x     = NULL,
    y     = "Mean CRISPR Dependency Score",
    fill  = NULL, colour = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position   = "none",
    strip.text        = element_text(face = "bold", size = 12),
    strip.background  = element_blank(),
    axis.text.x       = element_text(size = 11),
    plot.title        = element_text(face = "bold", hjust = 0.5, size = 12)
  )

ggsave(file.path(out_dir, "6_3_violin_matched.pdf"),
       plot = violin_plot, width = 8, height = 6, device = cairo_pdf)
ggsave(file.path(out_dir, "6_3_violin_matched.png"),
       plot = violin_plot, width = 8, height = 6, dpi = 300)

cat("  6_3_violin_matched.pdf / .png saved.\n")

# ══════════════════════════════════════════════════════════════
## PART 5 — Annotated gene tables (ver. 6_3)
##
## Full network gene tables updated with DepMap mean dep score
## and three boolean flags:
##   is_DEG      : P.Value < 0.01 & |logFC| > 1.5
##   is_hub      : top 10% most connected genes (degree_rank)
##   is_interest : is_DEG & is_hub  (DEGs ∩ Hub — intersection set)
##
## Output: 6_3_gene_table_A7.tsv  /  6_3_gene_table_A9.tsv
# ══════════════════════════════════════════════════════════════
cat("\n── PART 5: Annotated gene tables ───────────────────────\n")

build_annot_table <- function(annot_tbl, dep_mean_tbl) {
  n_hub     <- ceiling(0.10 * nrow(annot_tbl))
  hub_genes <- annot_tbl %>%
    arrange(degree_rank) %>%
    slice_head(n = n_hub) %>%
    pull(Gene) %>% toupper()

  deg_genes <- annot_tbl %>%
    filter(P.Value < 0.01, abs(logFC) > 1.5) %>%
    pull(Gene) %>% toupper()

  annot_tbl %>%
    mutate(gene_u = toupper(Gene)) %>%
    left_join(dep_mean_tbl, by = "gene_u") %>%
    mutate(
      is_DEG      = gene_u %in% deg_genes,
      is_hub      = gene_u %in% hub_genes,
      is_interest = is_DEG & is_hub
    ) %>%
    select(Gene, logFC, AveExpr, t, P.Value, adj.P.Val, B,
           degree, pagerank, degree_rank, pagerank_rank,
           mean_dep, is_DEG, is_hub, is_interest) %>%
    arrange(desc(is_interest), desc(is_hub), P.Value)
}

tbl_A7 <- build_annot_table(DE_A7, dep_mean_A7)
tbl_A9 <- build_annot_table(DE_A9, dep_mean_A9)

write_tsv(tbl_A7, file.path(out_dir, "6_3_gene_table_A7.tsv"))
write_tsv(tbl_A9, file.path(out_dir, "6_3_gene_table_A9.tsv"))

cat(sprintf("  A7 — total: %d | DEGs: %d | Hub: %d | Interest: %d\n",
            nrow(tbl_A7),
            sum(tbl_A7$is_DEG),
            sum(tbl_A7$is_hub),
            sum(tbl_A7$is_interest)))
cat(sprintf("  A9 — total: %d | DEGs: %d | Hub: %d | Interest: %d\n",
            nrow(tbl_A9),
            sum(tbl_A9$is_DEG),
            sum(tbl_A9$is_hub),
            sum(tbl_A9$is_interest)))
cat("  Saved: 6_3_gene_table_A7.tsv / 6_3_gene_table_A9.tsv\n")

# ══════════════════════════════════════════════════════════════
## PART 6 — Heatmap: CRISPR dependency × genes of interest
##
## Per-cell-line dep scores for genes flagged is_interest in
## either the A7 or A9 network.
## Rows: genes with dep_score < -0.3 in ≥1 line (essential).
## Cols: cervical cell lines ordered A7 → A9.
##
## Output: 6_3_heatmap_interest.pdf / .png
# ══════════════════════════════════════════════════════════════
cat("\n── PART 6: Heatmap — CRISPR dependency × genes of interest ─\n")

suppressPackageStartupMessages({
  library(depmap)
  library(ExperimentHub)
  library(ComplexHeatmap)
  library(circlize)
  library(tibble)
  library(tidyr)
})

# Union of interest genes from both networks
inter_A7 <- tbl_A7 %>% filter(is_interest) %>% pull(Gene) %>% toupper()
inter_A9 <- tbl_A9 %>% filter(is_interest) %>% pull(Gene) %>% toupper()
all_inter <- union(inter_A7, inter_A9)
cat(sprintf("  Interest genes — A7: %d | A9: %d | union: %d\n",
            length(inter_A7), length(inter_A9), length(all_inter)))

# Per-cell-line CRISPR scores (depmap package; locally cached after first run)
cervix_ids <- crispr_lines$depmap_id

crispr_long <- depmap_crispr() %>%
  filter(depmap_id %in% cervix_ids) %>%
  mutate(gene_u = toupper(gene_name)) %>%
  filter(gene_u %in% all_inter) %>%
  select(gene_u, cell_line, dep_score = dependency)

# gene × cell-line matrix
crispr_mat <- crispr_long %>%
  pivot_wider(names_from = cell_line, values_from = dep_score) %>%
  column_to_rownames("gene_u") %>%
  as.matrix()

# Keep genes with dep_score < -0.3 in at least one line
essential_flag <- apply(crispr_mat, 1, function(x) any(x < -0.3, na.rm = TRUE))
crispr_mat     <- crispr_mat[essential_flag, , drop = FALSE]
cat(sprintf("  Genes after essentiality filter (dep < -0.3): %d × %d cell lines\n",
            nrow(crispr_mat), ncol(crispr_mat)))

if (nrow(crispr_mat) == 0) {
  cat("  No essential interest genes found — skipping heatmap.\n")
} else {

  # Order columns: A7 first, then A9
  col_order  <- c(
    crispr_lines %>% filter(hpv_group == "A7", cell_line %in% colnames(crispr_mat)) %>% pull(cell_line),
    crispr_lines %>% filter(hpv_group == "A9", cell_line %in% colnames(crispr_mat)) %>% pull(cell_line)
  )
  col_order  <- col_order[col_order %in% colnames(crispr_mat)]
  crispr_mat <- crispr_mat[, col_order, drop = FALSE]
  colnames(crispr_mat) <- sub("_CERVIX$", "", colnames(crispr_mat))

  # ── Row clustering: most negative mean dep score at top ───────
  row_hc   <- hclust(dist(crispr_mat), method = "complete")
  row_dend <- reorder(as.dendrogram(row_hc),
                      rowMeans(crispr_mat, na.rm = TRUE),
                      agglo.FUN = mean)
  row_hc_ord <- as.hclust(row_dend)

  # ── Top 20 most essential genes (by mean dep score) ───────────
  top20_essential <- names(sort(rowMeans(crispr_mat, na.rm = TRUE)))[1:20]
  row_index <- which(rownames(crispr_mat) %in% top20_essential)

  # ── Shared HPV group colors ───────────────────────────────────
  hpv_colors <- c("A7 only" = "#27AE60", "A9 only" = "#E74C3C", "A7 & A9" = "#8E44AD")

  # ── Column annotation (one thin bar per cell line) ────────────
  col_hpv <- crispr_lines %>%
    mutate(line_short = sub("_CERVIX$", "", cell_line)) %>%
    filter(line_short %in% colnames(crispr_mat)) %>%
    arrange(match(line_short, colnames(crispr_mat))) %>%
    mutate(`HPV group` = ifelse(hpv_group == "A7", "A7 only", "A9 only")) %>%
    pull(`HPV group`)

  ha_top <- HeatmapAnnotation(
    `HPV group`          = col_hpv,
    col                  = list(`HPV group` = hpv_colors),
    show_legend          = FALSE,
    show_annotation_name = FALSE
  )

  # ── Row annotation (network origin of each gene) ──────────────
  row_hpv <- case_when(
    rownames(crispr_mat) %in% inter_A7 & rownames(crispr_mat) %in% inter_A9 ~ "A7 & A9",
    rownames(crispr_mat) %in% inter_A7 ~ "A7 only",
    rownames(crispr_mat) %in% inter_A9 ~ "A9 only"
  )

  ha_left <- rowAnnotation(
    `HPV group` = row_hpv,
    col         = list(`HPV group` = hpv_colors),
    width       = unit(3, "mm"),
    annotation_name_gp      = gpar(fontsize = 20),
    annotation_legend_param = list(
      title    = "HPV group",
      title_gp = gpar(fontsize = 22, fontface = "bold"),
      labels_gp = gpar(fontsize = 20)
    )
  )

  # ── Color ramp: centered on the actual data median (~-0.4)
  # These hub/DEG genes are genuinely essential: matrix median = -0.4,
  # 75th pct = -0.18, max = +0.19. Anchoring white at the median lets
  # the diverging scale reveal relative essentiality within this set.
  dep_median <- median(crispr_mat, na.rm = TRUE)
  cat(sprintf("  Matrix value range: %.3f … %.3f  median: %.3f\n",
              min(crispr_mat, na.rm=TRUE), max(crispr_mat, na.rm=TRUE), dep_median))
  col_fun <- colorRamp2(
    c(-2.5, -1.5, -0.8, dep_median, 0, 0.2),
    c("#7B241C", "#C0392B", "#E8A99A", "#FFFFFF", "#A9CCE3", "#1A6FB5")
  )

  # ── Heatmap ───────────────────────────────────────────────────
  ht <- Heatmap(
    crispr_mat,
    name              = "Dependency\nScore",
    col               = col_fun,
    cluster_rows      = row_hc_ord,
    cluster_columns   = FALSE,
    show_row_names    = FALSE,
    show_row_dend     = TRUE,
    row_dend_width    = unit(8, "mm"),
    show_column_names = TRUE,
    column_names_gp   = gpar(fontsize = 24, fontface = "bold"),
    column_names_rot  = 45,
    top_annotation    = ha_top,
    left_annotation   = ha_left,
    right_annotation  = rowAnnotation(
      mark = anno_mark(
        at         = row_index,
        labels     = rownames(crispr_mat)[row_index],
        labels_gp  = gpar(fontsize = 20, fontface = "italic"),
        link_width = unit(12, "mm"),
        link_gp    = gpar(lwd = 0.8, col = "grey50")
      )
    ),
    na_col            = "grey90",
    rect_gp           = gpar(col = NA),
    border            = FALSE,
    column_title      = "CRISPR Dependency: Genes of Interest\n× Cervical Cell Lines (DepMap)",
    column_title_gp   = gpar(fontsize = 26, fontface = "bold"),
    heatmap_legend_param = list(
      title_gp      = gpar(fontsize = 22, fontface = "bold"),
      labels_gp     = gpar(fontsize = 20),
      legend_height = unit(5, "cm"),
      direction     = "vertical"
    )
  )

  plot_h <- 14
  plot_w <- max(10, ncol(crispr_mat) * 0.70 + 5)

  cairo_pdf(file.path(out_dir, "6_3_heatmap_interest.pdf"),
            width = plot_w, height = plot_h)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
       merge_legend = TRUE)
  dev.off()

  png(file.path(out_dir, "6_3_heatmap_interest.png"),
      width  = round(plot_w * 150),
      height = round(plot_h * 150),
      res    = 150)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
       merge_legend = TRUE)
  dev.off()

  cat("  Saved: 6_3_heatmap_interest.pdf / .png\n")
}

cat("✓ All parts complete.\n")
