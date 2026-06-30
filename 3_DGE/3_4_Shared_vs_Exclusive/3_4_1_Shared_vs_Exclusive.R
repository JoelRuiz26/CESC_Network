# ============================================================
# 3_4_1_Shared_vs_Exclusive.R
#
# Compara los genes "de interes" (DEG ∩ Hub) segun su categoria:
#   - Exclusivos del clado (presentes solo en A7_int o solo en A9_int)
#   - Compartidos entre ambos clados (A7_int ∩ A9_int, n = 120)
#
# Analisis A: Posicion en la red (degree)
#   → ¿Los genes compartidos son hubs mas conectados en ambas redes?
#
# Analisis B: Magnitud de la expresion diferencial (|logFC|, |t|)
#   → ¿Los genes compartidos estan mas consistentemente desregulados?
#
# Input:
#   - ../../6_DepMap_cacu/6_3_gene_table_A7.tsv
#   - ../../6_DepMap_cacu/6_3_gene_table_A9.tsv
#
# Output:
#   - 3_4_1_A_Connectivity.pdf / .png
#   - 3_4_1_B_DiffExpr.pdf / .png
# ============================================================

setwd("~/CESC_Network/3_DGE/3_4_Shared_vs_Exclusive/")

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

# ── 1) Cargar y clasificar ──────────────────────────────────────

t7 <- read_tsv("../../6_DepMap_cacu/6_3_gene_table_A7.tsv", show_col_types = FALSE)
t9 <- read_tsv("../../6_DepMap_cacu/6_3_gene_table_A9.tsv", show_col_types = FALSE)

a7_int <- t7 %>% filter(is_interest)
a9_int <- t9 %>% filter(is_interest)

shared_genes <- intersect(a7_int$Gene, a9_int$Gene)

a7_int <- a7_int %>%
  mutate(
    abs_logFC = abs(logFC),
    abs_t     = abs(t),
    category  = factor(
      ifelse(Gene %in% shared_genes, "Compartido", "Exclusivo A7"),
      levels = c("Exclusivo A7", "Compartido")
    )
  )

a9_int <- a9_int %>%
  mutate(
    abs_logFC = abs(logFC),
    abs_t     = abs(t),
    category  = factor(
      ifelse(Gene %in% shared_genes, "Compartido", "Exclusivo A9"),
      levels = c("Exclusivo A9", "Compartido")
    )
  )

cat(sprintf(
  "A7_int: %d total  (%d exclusivos, %d compartidos)\n",
  nrow(a7_int),
  sum(a7_int$category == "Exclusivo A7"),
  sum(a7_int$category == "Compartido")
))
cat(sprintf(
  "A9_int: %d total  (%d exclusivos, %d compartidos)\n",
  nrow(a9_int),
  sum(a9_int$category == "Exclusivo A9"),
  sum(a9_int$category == "Compartido")
))

# ── 2) Paleta ───────────────────────────────────────────────────

COL_A7 <- "#2E86AB"
COL_A9 <- "#D35400"
COL_SH <- "#27AE60"

PAL_A7 <- c("Exclusivo A7" = COL_A7, "Compartido" = COL_SH)
PAL_A9 <- c("Exclusivo A9" = COL_A9, "Compartido" = COL_SH)

# ── 3) Helper: boxplot + jitter + Wilcoxon ──────────────────────

make_box <- function(df, yvar, pal, ylab, title) {

  y_vals <- df[[yvar]]
  cats   <- levels(df$category)
  is_c1  <- df$category == cats[1]
  is_c2  <- df$category == cats[2]

  p_val <- wilcox.test(y_vals[is_c1], y_vals[is_c2])$p.value
  p_lbl <- if (p_val < 0.001) "p < 0.001" else sprintf("p = %.3f", p_val)

  y_max  <- max(y_vals, na.rm = TRUE)
  y_seg  <- y_max * 1.10
  y_text <- y_max * 1.18

  # counts per category for x-axis labels
  n_labs <- df %>%
    count(category) %>%
    mutate(label = sprintf("%s\n(n = %d)", category, n))
  lbl_vec <- setNames(n_labs$label, n_labs$category)

  ggplot(df, aes(x = category, y = .data[[yvar]], fill = category)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.65, width = 0.45,
                 linewidth = 0.5) +
    geom_jitter(aes(color = category), width = 0.18, size = 0.9,
                alpha = 0.45) +
    scale_fill_manual(values = pal,  guide = "none") +
    scale_color_manual(values = pal, guide = "none") +
    scale_x_discrete(labels = lbl_vec) +
    coord_cartesian(ylim = c(NA, y_text * 1.10)) +
    annotate("segment",
             x = 1, xend = 2,
             y = y_seg, yend = y_seg,
             color = "grey50", linewidth = 0.45) +
    annotate("segment", x = 1, xend = 1,
             y = y_seg, yend = y_seg * 0.97,
             color = "grey50", linewidth = 0.45) +
    annotate("segment", x = 2, xend = 2,
             y = y_seg, yend = y_seg * 0.97,
             color = "grey50", linewidth = 0.45) +
    annotate("text",
             x = 1.5, y = y_text,
             label = p_lbl, size = 5, color = "grey25") +
    labs(title = title, y = ylab, x = NULL) +
    theme_bw(base_size = 15) +
    theme(
      plot.title   = element_text(face = "bold", size = 15),
      axis.text.x  = element_text(size = 14),
      axis.text.y  = element_text(size = 13),
      axis.title.y = element_text(size = 14),
      panel.grid.major.x = element_blank()
    )
}

# ── 4) Analisis A: Conectividad (degree) ────────────────────────

pA1 <- make_box(
  a7_int, "degree", PAL_A7,
  ylab  = "Degree (red A7)",
  title = "Conectividad en red A7\nExclusivo vs. Compartido"
)

pA2 <- make_box(
  a9_int, "degree", PAL_A9,
  ylab  = "Degree (red A9)",
  title = "Conectividad en red A9\nExclusivo vs. Compartido"
)

fig_A <- pA1 + pA2 +
  plot_annotation(
    title   = "Analisis A: Posicion en la red de co-expresion",
    subtitle = "Genes de interes (DEG x Hub): exclusivos vs. compartidos entre clados",
    theme   = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "grey40"),
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

ggsave("3_4_1_A_Connectivity.pdf", fig_A, width = 9, height = 5)
ggsave("3_4_1_A_Connectivity.png", fig_A, width = 9, height = 5,
       dpi = 300, bg = "white")

# ── 5) Analisis B: Estadistico t ────────────────────────────────

pB1 <- make_box(
  a7_int, "abs_t", PAL_A7,
  ylab  = "|t|",
  title = "Estadistico t A7 (vs. normal)\nExclusivo vs. Compartido"
)

pB2 <- make_box(
  a9_int, "abs_t", PAL_A9,
  ylab  = "|t|",
  title = "Estadistico t A9 (vs. normal)\nExclusivo vs. Compartido"
)

fig_B <- pB1 + pB2 +
  plot_annotation(
    title   = "Analisis B: Consistencia de la expresion diferencial",
    subtitle = "Genes de interes (DEG x Hub): exclusivos vs. compartidos entre clados",
    theme   = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "grey40"),
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

ggsave("3_4_1_B_Tstat.pdf", fig_B, width = 9, height = 5)
ggsave("3_4_1_B_Tstat.png", fig_B, width = 9, height = 5,
       dpi = 300, bg = "white")

cat("Listo: 3_4_1_A_Connectivity y 3_4_1_B_Tstat guardados.\n")
