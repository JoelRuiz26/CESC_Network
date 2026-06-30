############################################################
## 5_5_Interest_in_modules.R
##
## Distribución de genes de interés (hub ∩ DEG, DepMap)
## dentro de los módulos Louvain de cada red HPV.
##
## Inputs:
##   7_Modularity/5_1_Modules/5_1_1_Louvain_A7_graph_max.rds
##   7_Modularity/5_1_Modules/5_1_1_Louvain_A9_graph_max.rds
##   6_DepMap_cacu/6_3_gene_table_A7.tsv
##   6_DepMap_cacu/6_3_gene_table_A9.tsv
##
## Outputs (en 7_Modularity/5_5_DepMap_interest/):
##   5_5_module_interest_A7.tsv / _A9.tsv   — tabla por módulo
##   5_5_barplot_count.pdf / .png            — nº genes interés × módulo
##   5_5_barplot_proportion.pdf / .png       — % genes interés × módulo
##   5_5_hypergeometric.pdf / .png           — -log10(p.adj) × módulo
############################################################

suppressPackageStartupMessages({
  library(igraph)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(stringr)
})

out_dir  <- "~/CESC_Network/7_Modularity/5_5_DepMap_interest"
mod_dir  <- "~/CESC_Network/7_Modularity/5_1_Modules"
dep_dir  <- "~/CESC_Network/6_DepMap_cacu"

PAL <- c("HPV-A7" = "#1A6FB5", "HPV-A9" = "#C0392B")

# ══════════════════════════════════════════════════════════════
## 1) Cargar datos
# ══════════════════════════════════════════════════════════════
load_module_interest <- function(graph_rds, gene_tsv, label) {
  g   <- readRDS(graph_rds)
  vdf <- igraph::as_data_frame(g, what = "vertices") %>%
    rename(Gene = name) %>%
    mutate(
      Community      = as.integer(Community),
      Community_name = as.character(Community_name),
      Community_size = as.integer(Community_size)
    )

  tbl <- read_tsv(gene_tsv, show_col_types = FALSE)

  dep_genes <- tbl %>%
    filter(is_interest, !is.na(mean_dep), mean_dep <= -0.5) %>%
    pull(Gene)

  vdf %>%
    mutate(
      is_interest = Gene %in% dep_genes,
      group       = label
    )
}

df_A7 <- load_module_interest(
  file.path(mod_dir, "5_1_1_Louvain_A7_graph_max.rds"),
  file.path(dep_dir, "6_3_gene_table_A7.tsv"), "HPV-A7"
)
df_A9 <- load_module_interest(
  file.path(mod_dir, "5_1_1_Louvain_A9_graph_max.rds"),
  file.path(dep_dir, "6_3_gene_table_A9.tsv"), "HPV-A9"
)

cat(sprintf("A7: %d módulos | %d genes de interés\n",
            n_distinct(df_A7$Community), sum(df_A7$is_interest)))
cat(sprintf("A9: %d módulos | %d genes de interés\n",
            n_distinct(df_A9$Community), sum(df_A9$is_interest)))

# ══════════════════════════════════════════════════════════════
## 2) Tabla resumen por módulo + test hipergeométrico
# ══════════════════════════════════════════════════════════════
make_module_stats <- function(df) {
  N <- nrow(df)
  M <- sum(df$is_interest)

  df %>%
    group_by(Community, Community_name, Community_size) %>%
    summarise(
      n_genes    = n(),
      n_interest = sum(is_interest),
      .groups    = "drop"
    ) %>%
    mutate(
      pct_interest = n_interest / n_genes * 100,
      global_pct   = M / N * 100,
      # P(X >= k) hipergeométrica: enriquecimiento
      p_hyper = phyper(n_interest - 1, M, N - M, n_genes, lower.tail = FALSE),
      p_adj   = p.adjust(p_hyper, method = "BH"),
      sig     = case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01  ~ "**",
        p_adj < 0.05  ~ "*",
        TRUE          ~ ""
      ),
      enriched = p_adj < 0.05
    ) %>%
    arrange(desc(n_interest))
}

stats_A7 <- make_module_stats(df_A7)
stats_A9 <- make_module_stats(df_A9)

write_tsv(stats_A7, file.path(out_dir, "5_5_module_interest_A7.tsv"))
write_tsv(stats_A9, file.path(out_dir, "5_5_module_interest_A9.tsv"))

cat("\n── Módulos significativos (BH p < 0.05) ──────────────\n")
cat(sprintf("A7: %d / %d módulos\n",
            sum(stats_A7$enriched), nrow(stats_A7)))
cat(sprintf("A9: %d / %d módulos\n",
            sum(stats_A9$enriched), nrow(stats_A9)))

# ══════════════════════════════════════════════════════════════
## 3) Funciones de plot
# ══════════════════════════════════════════════════════════════

# Etiqueta corta para el módulo (nombre del hub + nº comunidad)
mod_label <- function(name, comm) {
  lbl <- str_trunc(name, 14, ellipsis = "…")
  sprintf("%s\n(M%s)", lbl, comm)
}

# ── 3a) Barplot: nº genes de interés por módulo ──────────────
plot_count <- function(stats, label, accent) {
  df <- stats %>%
    mutate(
      mod_lbl = factor(mod_label(Community_name, Community),
                       levels = mod_label(Community_name[order(desc(n_interest))],
                                         Community[order(desc(n_interest))]))
    )

  ggplot(df, aes(x = mod_lbl, y = n_interest, fill = enriched)) +
    geom_col(width = 0.7, colour = NA) +
    geom_text(aes(label = ifelse(sig != "", sig, "")),
              vjust = -0.3, size = 4.5, colour = "grey20") +
    scale_fill_manual(
      values = c("TRUE" = accent, "FALSE" = "grey75"),
      labels = c("TRUE" = "Enriched (p.adj<0.05)", "FALSE" = "Not significant"),
      name   = NULL
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(
      title = label,
      x     = "Louvain Module",
      y     = "N interest genes"
    ) +
    theme_classic(base_size = 13) +
    theme(
      axis.text.x     = element_text(size = 8,  angle = 45, hjust = 1),
      axis.text.y     = element_text(size = 11, face = "bold"),
      plot.title      = element_text(face = "bold", colour = accent, size = 14),
      legend.position = "top",
      legend.text     = element_text(size = 11),
      panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.3)
    )
}

# ── 3b) Barplot: % interés por módulo (+ línea global) ───────
plot_prop <- function(stats, label, accent) {
  global <- stats$global_pct[1]

  df <- stats %>%
    mutate(
      mod_lbl = factor(mod_label(Community_name, Community),
                       levels = mod_label(Community_name[order(desc(pct_interest))],
                                         Community[order(desc(pct_interest))]))
    )

  ggplot(df, aes(x = mod_lbl, y = pct_interest, fill = enriched)) +
    geom_col(width = 0.7, colour = NA) +
    geom_hline(yintercept = global, linetype = "dashed",
               colour = "grey40", linewidth = 0.6) +
    annotate("text", x = Inf, y = global,
             label = sprintf(" Global %.1f%%", global),
             hjust = 1.05, vjust = -0.5, size = 3.5, colour = "grey40") +
    geom_text(aes(label = ifelse(sig != "", sig, "")),
              vjust = -0.3, size = 4.5, colour = "grey20") +
    scale_fill_manual(
      values = c("TRUE" = accent, "FALSE" = "grey75"),
      labels = c("TRUE" = "Enriched (p.adj<0.05)", "FALSE" = "Not significant"),
      name   = NULL
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)),
                       labels = function(x) paste0(round(x, 1), "%")) +
    labs(
      title = label,
      x     = "Louvain Module",
      y     = "% interest genes in module"
    ) +
    theme_classic(base_size = 13) +
    theme(
      axis.text.x     = element_text(size = 8,  angle = 45, hjust = 1),
      axis.text.y     = element_text(size = 11, face = "bold"),
      plot.title      = element_text(face = "bold", colour = accent, size = 14),
      legend.position = "top",
      legend.text     = element_text(size = 11),
      panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.3)
    )
}

# ── 3c) Lollipop: -log10(p.adj) por módulo ───────────────────
plot_hyper <- function(stats, label, accent) {
  thresh <- -log10(0.05)

  df <- stats %>%
    mutate(
      log_padj = -log10(pmax(p_adj, 1e-10)),
      mod_lbl  = reorder(
        sprintf("%s  (n=%d)", str_trunc(Community_name, 13, ellipsis = "…"), n_genes),
        log_padj
      ),
      n_lab = as.character(n_interest)
    )

  ggplot(df, aes(x = log_padj, y = mod_lbl)) +
    annotate("rect",
             xmin = thresh, xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = accent, alpha = 0.07) +
    geom_vline(xintercept = thresh, linetype = "dashed",
               colour = accent, linewidth = 0.65, alpha = 0.75) +
    annotate("text", x = thresh, y = -Inf,
             label = "  p.adj=0.05", hjust = 0, vjust = 1.6,
             size = 3.0, colour = accent, fontface = "italic") +
    geom_segment(aes(x = 0, xend = log_padj, yend = mod_lbl,
                     colour = log_padj),
                 linewidth = 1.3, lineend = "round") +
    geom_point(aes(colour = log_padj, size = log_padj)) +
    geom_text(aes(label = n_lab),
              hjust = -0.45, size = 3.4, fontface = "bold",
              colour = "grey25") +
    scale_colour_gradient(
      low  = "grey82", high = accent,
      name = expression(-log[10](p[adj]))
    ) +
    scale_size_continuous(range = c(2.5, 7.5), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.22))) +
    labs(
      title    = label,
      subtitle = sprintf("N = %d modules  ·  %d interest genes",
                         nrow(df), sum(stats$n_interest)),
      x        = expression(-log[10](p[adj])~~"(hypergeometric, BH)"),
      y        = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title       = element_text(face = "bold", colour = accent,
                                      size = 15, margin = margin(b = 2)),
      plot.subtitle    = element_text(colour = "grey45", size = 10,
                                      margin = margin(b = 8)),
      axis.text.y      = element_text(size = 10, face = "bold", colour = "grey15"),
      axis.text.x      = element_text(size = 11, face = "bold"),
      axis.title.x     = element_text(size = 11, face = "bold", colour = "grey30"),
      legend.position  = "bottom",
      legend.title     = element_text(size = 10, face = "bold"),
      legend.text      = element_text(size = 9),
      legend.key.width = unit(1.4, "cm"),
      panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.35),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(colour = "grey94", linewidth = 0.25),
      panel.border     = element_rect(colour = "grey70", linewidth = 0.5),
      plot.background  = element_rect(fill = "white", colour = NA)
    )
}

# ══════════════════════════════════════════════════════════════
## 4) Generar y guardar plots
# ══════════════════════════════════════════════════════════════
save_plot <- function(p, fname, w, h) {
  ggsave(file.path(out_dir, paste0(fname, ".pdf")),
         plot = p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(out_dir, paste0(fname, ".png")),
         plot = p, width = w, height = h, dpi = 300)
  cat(sprintf("  Guardado: %s.pdf/png\n", fname))
}

# — Barplot conteo —
p_count <- (plot_count(stats_A7, "HPV-A7", PAL["HPV-A7"]) /
            plot_count(stats_A9, "HPV-A9", PAL["HPV-A9"])) +
  plot_annotation(
    title   = "Interest Genes (Hub ∩ DEG, DepMap dep ≤ -0.5) per Louvain Module",
    caption = "Colored bars = enriched modules (hypergeometric, BH p < 0.05)  ·  * p<0.05  ** p<0.01  *** p<0.001",
    theme   = theme(
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 15),
      plot.caption = element_text(colour = "grey50", size = 10, hjust = 0.5)
    )
  )

n_mod <- max(nrow(stats_A7), nrow(stats_A9))
save_plot(p_count, "5_5_barplot_count",
          w = max(10, n_mod * 0.55 + 2), h = 9)

# — Barplot proporción —
p_prop <- (plot_prop(stats_A7, "HPV-A7", PAL["HPV-A7"]) /
           plot_prop(stats_A9, "HPV-A9", PAL["HPV-A9"])) +
  plot_annotation(
    title   = "Proportion of Interest Genes (dep ≤ -0.5) per Louvain Module",
    caption = "Dashed line = network-wide %  ·  Colored bars = enriched modules (BH p < 0.05)",
    theme   = theme(
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 15),
      plot.caption = element_text(colour = "grey50", size = 10, hjust = 0.5)
    )
  )

save_plot(p_prop, "5_5_barplot_proportion",
          w = max(10, n_mod * 0.55 + 2), h = 9)

# — Lollipop hipergeométrico —
p_hyper <- (plot_hyper(stats_A7, "HPV-A7", PAL["HPV-A7"]) /
            plot_hyper(stats_A9, "HPV-A9", PAL["HPV-A9"])) +
  plot_annotation(
    title   = "Module Enrichment of Interest Genes (dep ≤ -0.5) | Louvain Communities",
    caption = "Hypergeometric test · Benjamini-Hochberg correction · Numbers indicate interest gene count per module",
    theme   = theme(
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 16),
      plot.caption = element_text(colour = "grey50", size = 9.5, hjust = 0.5)
    )
  )

save_plot(p_hyper, "5_5_hypergeometric", w = 6.5, h = 13)

cat("\n✓ 5_5 completo.\n")
