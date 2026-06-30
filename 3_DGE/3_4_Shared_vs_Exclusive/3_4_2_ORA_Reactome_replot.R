# ============================================================
# 3_4_2_ORA_Reactome_replot.R
#
# Script autónomo: corre el ORA de Reactome y genera el plot.
# No depende de caché. Corre solo este archivo para iterar.
# ============================================================

setwd("~/CESC_Network/3_DGE/3_4_Shared_vs_Exclusive/")

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(msigdbr)
  library(tidyverse)
  library(patchwork)
  library(colorspace)
})

# ── Parámetros (ajustar aquí) ──────────────────────────────────
TOP_N      <- 15
FONT_MULT  <- 6.0
WRAP_WIDTH <- 18
PLOT_W     <- 76
PLOT_H     <- 52
OUT_FILE   <- "3_4_2_ORA_Reactome"

# ── Helpers ────────────────────────────────────────────────────
reorder_within <- function(x, by, within, ...) {
  stats::reorder(paste(x, within, sep = "___"), by, ...)
}
scale_y_reordered <- function(..., sep = "___") {
  scale_y_discrete(labels = function(x) sub("___.*", "", x), ...)
}
clean_description <- function(id_vec) {
  id_vec %>%
    str_replace_all("_", " ") %>%
    str_to_title() %>%
    str_remove_all(regex("^Reactome\\s*", ignore_case = TRUE)) %>%
    str_squish()
}

GROUP_LEVELS <- c("Compartidos (n=120)", "Exclusivos A7 (n=228)", "Exclusivos A9 (n=115)")
GROUP_COL <- c(
  "Compartidos (n=120)"   = "#27AE60",
  "Exclusivos A7 (n=228)" = "#2E86AB",
  "Exclusivos A9 (n=115)" = "#D35400"
)

# ── 1) Universo ────────────────────────────────────────────────
universe_sym <- readRDS(
  "~/CESC_Network/2_Prepro_TCGA_GTEx/2_1_Universe_annotation_table.rds"
) %>% pull(hgnc_symbol) %>% na.omit() %>% unique()
cat(sprintf("Universo: %d genes\n", length(universe_sym)))

# ── 2) Gene sets ───────────────────────────────────────────────
t7 <- read_tsv("../../6_DepMap_cacu/6_3_gene_table_A7.tsv", show_col_types = FALSE)
t9 <- read_tsv("../../6_DepMap_cacu/6_3_gene_table_A9.tsv", show_col_types = FALSE)

gene_sets <- list(
  "Compartidos (n=120)"   = intersect(
    t7 %>% filter(is_interest) %>% pull(Gene),
    t9 %>% filter(is_interest) %>% pull(Gene)
  ),
  "Exclusivos A7 (n=228)" = setdiff(
    t7 %>% filter(is_interest) %>% pull(Gene),
    t9 %>% filter(is_interest) %>% pull(Gene)
  ),
  "Exclusivos A9 (n=115)" = setdiff(
    t9 %>% filter(is_interest) %>% pull(Gene),
    t7 %>% filter(is_interest) %>% pull(Gene)
  )
)
cat(sprintf("Genes — Compartidos: %d | A7: %d | A9: %d\n",
            length(gene_sets[[1]]), length(gene_sets[[2]]), length(gene_sets[[3]])))

# ── 3) Gene sets Reactome ──────────────────────────────────────
t2g <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME") %>%
  select(gs_name, gene_symbol) %>% distinct()
cat(sprintf("Reactome gene sets: %d\n", n_distinct(t2g$gs_name)))

# ── 4) ORA ────────────────────────────────────────────────────
ora_all <- imap_dfr(gene_sets, function(genes, label) {
  cat(sprintf("  [%s] ...", label))
  res <- enricher(
    gene          = genes,
    universe      = universe_sym,
    TERM2GENE     = t2g,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.10,   # BH-adj p < 0.10 para tener suficientes términos
    qvalueCutoff  = 1.00
  )
  if (is.null(res) || nrow(res) == 0) { cat(" sin términos.\n"); return(NULL) }
  cat(sprintf(" %d términos\n", nrow(res)))
  as.data.frame(res) %>%
    filter(p.adjust < 0.10) %>%
    mutate(
      group         = label,
      GeneRatio_num = map_dbl(GeneRatio, ~ {
        p <- as.numeric(strsplit(.x, "/")[[1]]); p[1] / p[2]
      }),
      Description = clean_description(ID)
    )
}) %>%
  mutate(group = factor(group, levels = GROUP_LEVELS))

cat(sprintf("\nTotal términos: %d\n", nrow(ora_all)))
ora_all %>% count(group) %>% print()

# ── 5) Plot ────────────────────────────────────────────────────
make_panel <- function(grp, df_grp) {
  if (nrow(df_grp) == 0) return(NULL)
  accent <- GROUP_COL[grp]

  df_plot <- df_grp %>%
    mutate(
      log_padj  = -log10(p.adjust),
      Desc_wrap = str_wrap(Description, width = WRAP_WIDTH)
    ) %>%
    slice_max(GeneRatio_num, n = TOP_N, with_ties = FALSE) %>%
    mutate(term_ord = reorder_within(Desc_wrap, GeneRatio_num, grp))

  x_max   <- max(df_plot$GeneRatio_num)
  x_breaks <- pretty(c(0, x_max), n = 3)
  x_breaks <- x_breaks[x_breaks <= x_max * 1.01]

  padj_vals <- df_plot$log_padj
  leg_breaks <- round(pretty(padj_vals, n = 3), 2)

  ggplot(df_plot, aes(x = GeneRatio_num, y = term_ord, fill = log_padj)) +
    geom_col(width = 0.72, colour = NA) +
    scale_y_reordered() +
    scale_fill_gradient(
      low    = lighten(accent, 0.55),
      high   = darken(accent, 0.35),
      name   = expression(bold(-log[10](p[adj]))),
      breaks = leg_breaks,
      limits = range(padj_vals),
      guide  = guide_colorbar(
        barheight       = unit(3.5 * FONT_MULT, "cm"),
        barwidth        = unit(0.6  * FONT_MULT, "cm"),
        ticks.linewidth = 1.5,
        label.hjust     = 0.5,
        frame.colour    = "grey40"
      )
    ) +
    scale_x_continuous(
      breaks = x_breaks,
      labels = scales::number_format(accuracy = 0.01),
      expand = expansion(mult = c(0, 0.18))
    ) +
    labs(title = grp, x = "Gene Ratio", y = NULL) +
    theme_classic(base_size = 14 * FONT_MULT) +
    theme(
      plot.title         = element_text(face = "bold", size = 16 * FONT_MULT,
                                        colour = accent, margin = margin(b = 8)),
      axis.text.y        = element_text(size = 11 * FONT_MULT, colour = "grey10", face = "bold"),
      axis.text.x        = element_text(size = 11 * FONT_MULT, face = "bold",
                                        colour = "grey20", angle = 45,
                                        hjust = 1, vjust = 1),
      axis.title.x       = element_text(size = 12 * FONT_MULT, face = "bold",
                                        colour = "grey30", margin = margin(t = 6)),
      axis.line          = element_line(colour = "grey70"),
      legend.title       = element_text(size = 11 * FONT_MULT, face = "bold",
                                        margin = margin(b = 6)),
      legend.text        = element_text(size = 10 * FONT_MULT),
      legend.key.height  = unit(1.2 * FONT_MULT, "cm"),
      panel.grid.major.x = element_line(colour = "grey93", linewidth = 0.4)
    )
}

panels <- imap(
  setNames(GROUP_LEVELS, GROUP_LEVELS),
  function(grp, ...) {
    df_grp <- ora_all %>% filter(as.character(group) == grp)
    make_panel(grp, df_grp)
  }
)
panels <- Filter(Negate(is.null), panels)

combined <- wrap_plots(panels, ncol = length(panels)) +
  plot_annotation(
    title   = "Reactome Pathways — Genes de interés (DEG x Hub) por categoría",
    caption = sprintf(
      "MSigDB · Universo: %d genes · BH-adj p < 0.10 · Top %d términos por Gene Ratio",
      length(universe_sym), TOP_N
    ),
    theme = theme(
      plot.title      = element_text(face = "bold", hjust = 0.5, size = 15 * FONT_MULT),
      plot.caption    = element_text(colour = "grey50", size =  9 * FONT_MULT, hjust = 0.5),
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

ggsave(paste0(OUT_FILE, ".pdf"), combined,
       width = PLOT_W, height = PLOT_H, device = cairo_pdf, limitsize = FALSE)
ggsave(paste0(OUT_FILE, ".png"), combined,
       width = PLOT_W, height = PLOT_H, dpi = 300, bg = "white", limitsize = FALSE)

cat(sprintf("\nGuardado: %s.pdf/png  (%.0f × %.0f in)\n", OUT_FILE, PLOT_W, PLOT_H))
