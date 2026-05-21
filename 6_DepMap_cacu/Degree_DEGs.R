############################################################
## DepMap – Cervical Cancer Dependency Score Analysis
##
## Gene sets compared (per list A7 / A9):
##   (A) Significant DEGs (p<0.01, |LFC|>1)   → N_input genes selected
##       mapped to DepMap                       → N_sig_dep with scores
##
##   (B) Top N_input Degree-ranked genes        → N_input genes selected
##       mapped to DepMap                       → N_deg_dep with scores
##       (n may differ from A after mapping –
##        hub genes have higher DepMap coverage)
##
##   (C) Intersection: sig DEG list ∩ top-degree list (pre-mapping)
##       mapped to DepMap                       → N_intersect with scores
##
## ECDF only (no raincloud panels)
##
## Statistical tests (one-sided, H1: group X more negative than DEGs):
##   - Wilcoxon (Mann-Whitney U)
##   - Cliff's delta  (δ < 0 → X more negative)
##   - Kolmogorov-Smirnov
############################################################

suppressPackageStartupMessages({
  library(depmap)
  library(ExperimentHub)
  library(dplyr)
  library(tidyr)
  library(vroom)
  library(ggplot2)
  library(patchwork)
  library(effsize)
})

# ── Colour palette (3 groups) ───────────────────────────────
PAL <- c(
  "Significant DEGs"               = "#E63946",
  "Top Degree-ranked"              = "#457B9D",
  "Intersection (DEG \u2229 Degree)" = "#2A9D8F"
)

############################################################
## 1) Load annotated DEG tables
############################################################
DE_A7_annot <- vroom(
  "~/CESC_Network/5_Network_analisis/5_4_DEG_integration/5_4_1_DEG_A7_Degree.tsv",
  show_col_types = FALSE
)
DE_A9_annot <- vroom(
  "~/CESC_Network/5_Network_analisis/5_4_DEG_integration/5_4_1_DEG_A9_Degree.tsv",
  show_col_types = FALSE
)

############################################################
## 2) DepMap – Cervical Cancer metadata + CRISPR scores
############################################################
meta <- depmap_metadata()

cervix_lines <- meta %>%
  filter(
    sample_collection_site == "cervix",
    grepl("cervical cancer", primary_disease, ignore.case = TRUE)
  ) %>%
  pull(depmap_id)

cat("Cervical cancer cell lines found:", length(cervix_lines), "\n")

crispr_cervix_all <- depmap_crispr() %>%
  filter(depmap_id %in% cervix_lines) %>%
  mutate(gene_u = toupper(gene_name)) %>%
  rename(dep_score = dependency)

cat("CRISPR rows (cervix):", nrow(crispr_cervix_all), "\n")
cat("Cell lines (cervix): ", n_distinct(crispr_cervix_all$cell_line), "\n")

# Average dependency score across cell lines per gene
dep_mean <- crispr_cervix_all %>%
  group_by(gene_u) %>%
  summarise(
    mean_dep   = mean(dep_score, na.rm = TRUE),
    median_dep = median(dep_score, na.rm = TRUE),
    n_lines    = sum(!is.na(dep_score)),
    .groups    = "drop"
  )

############################################################
## 3) Build gene sets
##
##  Selection logic (symmetric before mapping):
##    sig     = top N_input DEGs by significance
##    top_deg = top N_input genes by degree_rank
##  Both lists have the same N_input before DepMap join.
##  After joining dep_mean they may have different n
##  (hub genes are more often in DepMap libraries).
##
##  Intersection = genes in BOTH lists (pre-mapping),
##  then their DepMap scores are retrieved.
############################################################

build_sets <- function(annot_tbl, list_name) {
  
  # ── Significant DEGs (pre-mapping gene list) ───────────
  sig <- annot_tbl %>%
    filter(P.Value < 0.01, abs(logFC) > 1) %>%
    arrange(P.Value)
  
  N_input <- nrow(sig)
  cat(sprintf("[%s] Significant DEGs (p<0.01, |LFC|>1): %d\n", list_name, N_input))
  
  # ── Top N_input Degree-ranked (pre-mapping gene list) ──
  top_deg <- annot_tbl %>%
    arrange(degree_rank) %>%        # rank 1 = highest degree
    slice_head(n = N_input)
  
  cat(sprintf("[%s] Top %d degree-ranked genes selected\n", list_name, N_input))
  
  # ── Intersection of the two gene lists (pre-mapping) ──
  shared_genes <- intersect(toupper(sig$Gene), toupper(top_deg$Gene))
  cat(sprintf("[%s] Gene-list intersection (pre-mapping): %d\n",
              list_name, length(shared_genes)))
  
  # ── Map each list to DepMap independently ──────────────
  dep_sig <- dep_mean %>%
    filter(gene_u %in% toupper(sig$Gene)) %>%
    mutate(group = "Significant DEGs", list = list_name)
  
  dep_degree <- dep_mean %>%
    filter(gene_u %in% toupper(top_deg$Gene)) %>%
    mutate(group = "Top Degree-ranked", list = list_name)
  
  dep_intersect <- dep_mean %>%
    filter(gene_u %in% shared_genes) %>%
    mutate(group = "Intersection (DEG \u2229 Degree)", list = list_name)
  
  cat(sprintf(
    "[%s] DepMap coverage — DEGs: %d/%d (%.1f%%)  |  Degree: %d/%d (%.1f%%)  |  Intersect: %d/%d (%.1f%%)\n",
    list_name,
    nrow(dep_sig),       N_input,             100 * nrow(dep_sig)       / N_input,
    nrow(dep_degree),    N_input,             100 * nrow(dep_degree)    / N_input,
    nrow(dep_intersect), length(shared_genes), 100 * nrow(dep_intersect) / length(shared_genes)
  ))
  
  list(
    dep_sig       = dep_sig,
    dep_degree    = dep_degree,
    dep_intersect = dep_intersect,
    N_input       = N_input,
    N_sig_dep     = nrow(dep_sig),
    N_deg_dep     = nrow(dep_degree),
    N_intersect   = nrow(dep_intersect)
  )
}

sets_A7 <- build_sets(DE_A7_annot, "A7")
sets_A9 <- build_sets(DE_A9_annot, "A9")

############################################################
## 4) Statistical tests
##   Pairwise (one-sided: H1 = left group more negative):
##   (i)   Degree    vs DEGs
##   (ii)  Intersect vs DEGs
##   (iii) Intersect vs Degree
############################################################

run_one_test <- function(x, y, label_x, label_y) {
  wt <- wilcox.test(x, y, alternative = "less", exact = FALSE)
  cd <- cliff.delta(x, y)
  ks <- ks.test(x, y, alternative = "less")
  
  cat(sprintf(
    "\n  -- %s  vs  %s  (H1: %s < %s) --\n",
    label_x, label_y, label_x, label_y
  ))
  cat(sprintf("     n: %d vs %d\n", length(x), length(y)))
  cat(sprintf("     median: %.4f vs %.4f\n", median(x), median(y)))
  cat(sprintf("     Wilcoxon W=%.0f  p=%.4g\n", wt$statistic, wt$p.value))
  cat(sprintf("     Cliff's delta=%.4f [%s]  95%%CI [%.4f, %.4f]\n",
              cd$estimate, cd$magnitude, cd$conf.int[1], cd$conf.int[2]))
  cat(sprintf("     KS D=%.4f  p=%.4g\n", ks$statistic, ks$p.value))
  
  invisible(list(wilcox = wt, cliff = cd, ks = ks))
}

run_stats <- function(set, list_name) {
  x_sig    <- set$dep_sig$mean_dep
  x_degree <- set$dep_degree$mean_dep
  x_inter  <- set$dep_intersect$mean_dep
  
  cat(sprintf(
    "\n══════════════════════════════════════════════════\n %s — Statistical Tests\n══════════════════════════════════════════════════\n",
    list_name
  ))
  cat(sprintf(" Sig DEGs   n=%d  median=%.4f\n", length(x_sig),    median(x_sig)))
  cat(sprintf(" Degree     n=%d  median=%.4f\n", length(x_degree), median(x_degree)))
  cat(sprintf(" Intersect  n=%d  median=%.4f\n", length(x_inter),  median(x_inter)))
  
  r1 <- run_one_test(x_degree, x_sig,    "Degree",    "DEGs")
  r2 <- run_one_test(x_inter,  x_sig,    "Intersect", "DEGs")
  r3 <- run_one_test(x_inter,  x_degree, "Intersect", "Degree")
  
  invisible(list(deg_vs_sig = r1, int_vs_sig = r2, int_vs_deg = r3))
}

stats_A7 <- run_stats(sets_A7, "A7 vs SNT")
stats_A9 <- run_stats(sets_A9, "A9 vs SNT")

############################################################
## 5) ECDF panels  (3 curves per panel)
############################################################

make_ecdf_panel <- function(set, stats_res, list_name) {
  
  df <- bind_rows(set$dep_sig, set$dep_degree, set$dep_intersect) %>%
    mutate(group = factor(group, levels = c(
      "Significant DEGs",
      "Top Degree-ranked",
      "Intersection (DEG \u2229 Degree)"
    )))
  
  # Key stats for in-plot annotation
  delta_deg <- stats_res$deg_vs_sig$cliff$estimate
  p_deg     <- stats_res$deg_vs_sig$wilcox$p.value
  delta_int <- stats_res$int_vs_sig$cliff$estimate
  p_int     <- stats_res$int_vs_sig$wilcox$p.value
  
  fmt_p <- function(p) ifelse(p < 0.001, "< 0.001", sprintf("= %.4f", p))
  
  ann <- sprintf(
    "Degree vs DEGs:      \u03b4 = %.3f [%s], p %s\nIntersect vs DEGs:  \u03b4 = %.3f [%s], p %s",
    delta_deg, stats_res$deg_vs_sig$cliff$magnitude, fmt_p(p_deg),
    delta_int, stats_res$int_vs_sig$cliff$magnitude, fmt_p(p_int)
  )
  
  sub_txt <- sprintf(
    "Input N=%d each  |  DEGs mapped: %d  |  Degree mapped: %d  |  Intersection: %d",
    set$N_input, set$N_sig_dep, set$N_deg_dep, set$N_intersect
  )
  
  ggplot(df, aes(x = mean_dep, colour = group, linewidth = group)) +
    stat_ecdf(alpha = 0.92) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey45", linewidth = 0.4) +
    annotate("text",
             x      = -2.8,
             y      = 0.88,
             label  = ann,
             hjust  = 0,
             size   = 3,
             colour = "grey20",
             family = "mono") +
    scale_colour_manual(values = PAL) +
    scale_linewidth_manual(
      values = c(
        "Significant DEGs"                    = 1.0,
        "Top Degree-ranked"                   = 1.0,
        "Intersection (DEG \u2229 Degree)"    = 1.5
      )
    ) +
    guides(linewidth = "none") +
    labs(
      title    = sprintf("ECDF — %s vs SNT", list_name),
      subtitle = sub_txt,
      x        = "Mean CRISPR Dependency Score\n(more negative = stronger essentiality)",
      y        = "Cumulative probability",
      colour   = NULL,
      caption  = paste(
        "All three gene sets selected with the same N before DepMap mapping.",
        "Intersection = genes in both DEG and top-degree lists.",
        "Left-shifted curve \u2192 more essential. \u03b4 < 0 (Cliff's) \u2192 more negative scores."
      )
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position  = c(0.76, 0.22),
      legend.text      = element_text(size = 9),
      legend.key.width = unit(1.2, "cm"),
      plot.title       = element_text(face = "bold"),
      plot.subtitle    = element_text(colour = "grey40", size = 8.5),
      plot.caption     = element_text(colour = "grey55", size = 7.5)
    )
}

ecdf_A7 <- make_ecdf_panel(sets_A7, stats_A7, "A7")
ecdf_A9 <- make_ecdf_panel(sets_A9, stats_A9, "A9")

############################################################
## 6) Combined ECDF panel – save
############################################################
out_dir <- "~/CESC_Network/5_Network_analisis/5_4_DEG_integration/"

ecdf_combined <- ecdf_A7 + ecdf_A9 +
  plot_annotation(
    title   = "Empirical CDF \u2013 CRISPR Dependency Score | Cervical Cancer Cell Lines",
    caption = paste(
      "DEGs and Degree-ranked sets: top N selected before DepMap mapping (n may differ after join).",
      "Intersection = genes present in both pre-mapping lists.",
      "Left-shifted curve \u2192 stronger CRISPR essentiality."
    ),
    theme = theme(
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.caption = element_text(colour = "grey50", size = 8, hjust = 0.5)
    )
  )

ggsave(
  filename = file.path(out_dir, "5_4_2_DepMap_Dependency_ECDF.pdf"),
  plot     = ecdf_combined,
  width    = 14, height = 6,
  device   = cairo_pdf
)
ggsave(
  filename = file.path(out_dir, "5_4_2_DepMap_Dependency_ECDF.png"),
  plot     = ecdf_combined,
  width    = 14, height = 6,
  dpi      = 300
)

cat("\nECDF plots saved to:", out_dir, "\n")
cat("\u2713 Analysis complete.\n")