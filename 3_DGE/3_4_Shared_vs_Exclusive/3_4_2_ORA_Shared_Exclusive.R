# ============================================================
# 3_4_2_ORA_Shared_Exclusive.R
#
# ORA (MSigDB via msigdbr) para los tres subconjuntos de
# genes de interes (DEG x Hub):
#   - Compartidos entre clados  (A7_int ∩ A9_int, n = 120)
#   - Exclusivos A7             (A7_int \ A9_int, n = 228)
#   - Exclusivos A9             (A9_int \ A7_int, n = 115)
#
# Colecciones:
#   - Hallmark Gene Sets      (H)
#   - KEGG Pathways           (C2:CP:KEGG_LEGACY)
#   - Reactome Pathways       (C2:CP:REACTOME)
#   - WikiPathways            (C2:CP:WIKIPATHWAYS)
#   - GO Biological Process   (C5:GO:BP)
#
# Universo: 2_1_Universe_annotation_table.rds
#
# Output (uno por coleccion):
#   3_4_2_ORA_{Hallmark|KEGG|Reactome|WikiPathways|GOBP}.pdf/.png
#   3_4_2_ORA_all_results.tsv
# ============================================================

setwd("~/CESC_Network/3_DGE/3_4_Shared_vs_Exclusive/")

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(msigdbr)
  library(tidyverse)
  library(patchwork)
  library(colorspace)
})

# ── 1) Universo ────────────────────────────────────────────────

universe_sym <- readRDS(
  "~/CESC_Network/2_Prepro_TCGA_GTEx/2_1_Universe_annotation_table.rds"
) %>%
  pull(hgnc_symbol) %>%
  na.omit() %>%
  unique()

cat(sprintf("Universo: %d genes\n", length(universe_sym)))

# ── 2) Gene sets ───────────────────────────────────────────────

t7 <- read_tsv("../../6_DepMap_cacu/6_3_gene_table_A7.tsv", show_col_types = FALSE)
t9 <- read_tsv("../../6_DepMap_cacu/6_3_gene_table_A9.tsv", show_col_types = FALSE)

a7_int <- t7 %>% filter(is_interest) %>% pull(Gene)
a9_int <- t9 %>% filter(is_interest) %>% pull(Gene)

shared_genes <- intersect(a7_int, a9_int)
excl_a7      <- setdiff(a7_int, a9_int)
excl_a9      <- setdiff(a9_int, a7_int)

cat(sprintf(
  "Compartidos: %d | Exclusivos A7: %d | Exclusivos A9: %d\n",
  length(shared_genes), length(excl_a7), length(excl_a9)
))

GROUP_LEVELS <- c("Compartidos (n=120)", "Exclusivos A7 (n=228)", "Exclusivos A9 (n=115)")

gene_sets <- list(
  "Compartidos (n=120)"   = shared_genes,
  "Exclusivos A7 (n=228)" = excl_a7,
  "Exclusivos A9 (n=115)" = excl_a9
)

GROUP_COL <- c(
  "Compartidos (n=120)"   = "#27AE60",
  "Exclusivos A7 (n=228)" = "#2E86AB",
  "Exclusivos A9 (n=115)" = "#D35400"
)

# ── 3) Colecciones a correr ────────────────────────────────────

COLLECTIONS <- list(
  Hallmark = list(
    label      = "Hallmark Gene Sets",
    col        = "H",
    subcol     = NULL,
    top_n      = 20,
    wrap_width = 16,
    font_mult  = 6.0,
    plot_w     = 66,
    plot_h     = 40,
    file       = "3_4_2_ORA_Hallmark"
  ),
  KEGG = list(
    label      = "KEGG Pathways",
    col        = "C2",
    subcol     = "CP:KEGG_LEGACY",
    top_n      = 15,
    wrap_width = 34,
    font_mult  = 1.0,
    file       = "3_4_2_ORA_KEGG"
  ),
  Reactome = list(
    label      = "Reactome Pathways",
    col        = "C2",
    subcol     = "CP:REACTOME",
    top_n      = 15,
    wrap_width = 18,
    font_mult  = 6.0,
    plot_w     = 76,
    plot_h     = 52,
    file       = "3_4_2_ORA_Reactome"
  ),
  WikiPathways = list(
    label      = "WikiPathways",
    col        = "C2",
    subcol     = "CP:WIKIPATHWAYS",
    top_n      = 15,
    wrap_width = 38,
    font_mult  = 1.0,
    file       = "3_4_2_ORA_WikiPathways"
  ),
  GOBP = list(
    label      = "GO Biological Process",
    col        = "C5",
    subcol     = "GO:BP",
    top_n      = 15,
    wrap_width = 36,
    font_mult  = 1.0,
    file       = "3_4_2_ORA_GOBP"
  ),
  GOMF = list(
    label      = "GO Molecular Function",
    col        = "C5",
    subcol     = "GO:MF",
    top_n      = 15,
    wrap_width = 36,
    font_mult  = 1.0,
    file       = "3_4_2_ORA_GOMF"
  )
)

# ── 4) Helpers ─────────────────────────────────────────────────

reorder_within <- function(x, by, within, ...) {
  stats::reorder(paste(x, within, sep = "___"), by, ...)
}
scale_y_reordered <- function(..., sep = "___") {
  scale_y_discrete(labels = function(x) sub(paste0(sep, ".*"), "", x), ...)
}

clean_description <- function(id_vec) {
  id_vec %>%
    str_replace_all("_", " ") %>%
    str_to_title() %>%
    str_remove_all(
      regex("^(Hallmark|Kegg|Reactome|Gobp|Gocc|Gomf|Wp[0-9]+ |Wp )\\s*",
            ignore_case = TRUE)
    ) %>%
    str_squish()
}

run_ora <- function(genes, label, t2g) {
  cat(sprintf("    [%s] ...", label))
  res <- enricher(
    gene          = genes,
    universe      = universe_sym,
    TERM2GENE     = t2g,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 1.00
  )
  if (is.null(res) || nrow(res) == 0) {
    cat(" sin terminos significativos.\n")
    return(data.frame())
  }
  cat(sprintf(" %d terminos.\n", nrow(res)))
  as.data.frame(res) %>%
    mutate(
      group         = label,
      GeneRatio_num = map_dbl(GeneRatio, ~ {
        p <- as.numeric(strsplit(.x, "/")[[1]]); p[1] / p[2]
      }),
      Description = clean_description(ID)
    )
}

make_panel <- function(grp, df_grp, top_n, wrap_width, font_mult = 1) {
  if (is.null(df_grp) || nrow(df_grp) == 0) return(NULL)

  accent <- GROUP_COL[grp]

  df_plot <- df_grp %>%
    mutate(
      log_padj  = -log10(p.adjust),
      Desc_wrap = str_wrap(Description, width = wrap_width)
    ) %>%
    slice_max(GeneRatio_num, n = top_n, with_ties = FALSE) %>%
    mutate(term_ord = reorder_within(Desc_wrap, GeneRatio_num, grp))

  ggplot(df_plot, aes(x = GeneRatio_num, y = term_ord, fill = log_padj)) +
    geom_col(width = 0.7, colour = NA) +
    scale_y_reordered() +
    scale_fill_gradient(
      low    = lighten(accent, 0.55),
      high   = darken(accent, 0.35),
      name   = expression(bold(-log[10](p[adj]))),
      n.breaks = 3,
      guide  = guide_colorbar(
        barheight       = unit(2.2 * font_mult, "cm"),
        barwidth        = unit(0.5  * font_mult, "cm"),
        ticks.linewidth = 1.2,
        label.hjust     = 1
      )
    ) +
    scale_x_continuous(
      expand   = expansion(mult = c(0, 0.14)),
      n.breaks = 3
    ) +
    labs(title = grp, x = "Gene Ratio", y = NULL) +
    theme_classic(base_size = 14 * font_mult) +
    theme(
      plot.title         = element_text(face = "bold", size = 16 * font_mult, colour = accent),
      axis.text.y        = element_text(size = 12 * font_mult, colour = "grey10", face = "bold"),
      axis.text.x        = element_text(size = 12 * font_mult, face = "bold", colour = "grey20",
                                        angle = 45, hjust = 1, vjust = 1),
      axis.title.x       = element_text(size = 13 * font_mult, face = "bold", colour = "grey30"),
      axis.line          = element_line(colour = "grey70"),
      legend.title       = element_text(size = 12 * font_mult, face = "bold",
                                        margin = margin(b = 0.4 * font_mult, unit = "cm")),
      legend.text        = element_text(size = 11 * font_mult),
      legend.spacing.y   = unit(0.3 * font_mult, "cm"),
      panel.grid.major.x = element_line(colour = "grey93", linewidth = 0.35)
    )
}

# ── 5) Loop por coleccion ──────────────────────────────────────

all_results <- list()

for (col_id in names(COLLECTIONS)) {
  meta <- COLLECTIONS[[col_id]]
  cat(sprintf("\n== %s ==\n", meta$label))

  args <- list(species = "Homo sapiens", collection = meta$col)
  if (!is.null(meta$subcol)) args$subcollection <- meta$subcol
  t2g <- do.call(msigdbr, args) %>%
    select(gs_name, gene_symbol) %>%
    distinct()
  cat(sprintf("  Gene sets: %d\n", n_distinct(t2g$gs_name)))

  res_list <- imap(gene_sets, function(genes, label) {
    run_ora(genes, label, t2g)
  })

  ora_all <- bind_rows(res_list) %>%
    filter(!is.na(p.adjust)) %>%
    mutate(
      collection = col_id,
      group      = factor(group, levels = GROUP_LEVELS)
    )

  all_results[[col_id]] <- ora_all
  saveRDS(ora_all, sprintf("3_4_2_%s_oradata.rds", col_id))

  if (nrow(ora_all) == 0) {
    cat("  Sin resultados significativos — omitido.\n"); next
  }

  fm <- if (!is.null(meta$font_mult)) meta$font_mult else 1

  panels <- imap(
    setNames(GROUP_LEVELS, GROUP_LEVELS),
    function(grp, ...) {
      df_grp <- ora_all %>% filter(group == grp)
      make_panel(grp, df_grp, meta$top_n, meta$wrap_width, fm)
    }
  )
  panels <- Filter(Negate(is.null), panels)

  n_rows <- ora_all %>%
    group_by(group) %>%
    slice_max(GeneRatio_num, n = meta$top_n, with_ties = FALSE) %>%
    count() %>% pull(n) %>% max(na.rm = TRUE)

  plot_h <- if (!is.null(meta$plot_h)) meta$plot_h else n_rows * 0.52 * fm + 3.5
  plot_w <- if (!is.null(meta$plot_w)) meta$plot_w else length(panels) * 8 * max(1, fm * 0.75)

  combined <- wrap_plots(panels, ncol = length(panels)) +
    plot_annotation(
      title   = sprintf("%s — Genes de interes (DEG x Hub) por categoria", meta$label),
      caption = sprintf(
        "MSigDB | Universo: %d genes | BH-adj p < 0.05 | Top %d terminos por Gene Ratio",
        length(universe_sym), meta$top_n
      ),
      theme = theme(
        plot.title      = element_text(face = "bold", hjust = 0.5, size = 15 * fm),
        plot.caption    = element_text(colour = "grey50", size = 10 * fm, hjust = 0.5),
        plot.background = element_rect(fill = "white", color = NA)
      )
    )

  ggsave(paste0(meta$file, ".pdf"), combined,
         width = plot_w, height = plot_h, device = cairo_pdf)
  ggsave(paste0(meta$file, ".png"), combined,
         width = plot_w, height = plot_h, dpi = 300, bg = "white")

  cat(sprintf("  Guardado: %s.pdf/png\n", meta$file))
}

# ── 6) TSV combinado ───────────────────────────────────────────

bind_rows(all_results) %>%
  write_tsv("3_4_2_ORA_all_results.tsv")

cat("\nListo: todos los ORA guardados en 3_4_2_ORA_*.pdf/png\n")
