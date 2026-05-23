############################################################
## 6_4_ORA_interest_genes.R
##
## Over-Representation Analysis (ORA) sobre los genes de
## interés (hub ∩ DEG, is_interest = TRUE) de las redes
## HPV-A7 y HPV-A9.
##
## Ontologías : BP, MF, CC
## Universo   : 2_Prepro_TCGA_GTEx/2_1_Universe_annotation_table.rds
## Simplify   : se aplica si > 20 términos significativos
##              por grupo (cutoff semántico = 0.5)
##
## Outputs:
##   6_4_ORA_A7.tsv / 6_4_ORA_A9.tsv   — tablas completas
##   6_4_ORA_barplot.pdf / .png         — gráfico combinado
############################################################

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(forcats)
  library(stringr)
  library(patchwork)
})

# reorder_within / scale_y_reordered (from tidytext, implemented inline)
reorder_within <- function(x, by, within, ...) {
  new_x <- paste(x, within, sep = "___")
  stats::reorder(new_x, by, ...)
}
scale_y_reordered <- function(..., sep = "___") {
  ggplot2::scale_y_discrete(labels = function(x) sub(paste0(sep, ".*"), "", x), ...)
}

out_dir <- "~/CESC_Network/6_DepMap_cacu"

# ══════════════════════════════════════════════════════════════
## 1) Universo
# ══════════════════════════════════════════════════════════════
universe <- readRDS(
  "~/CESC_Network/2_Prepro_TCGA_GTEx/2_1_Universe_annotation_table.rds"
) %>%
  pull(hgnc_symbol) %>%
  na.omit() %>%
  unique()

cat(sprintf("Universo: %d genes\n", length(universe)))

# ══════════════════════════════════════════════════════════════
## 2) Genes de interés (hub ∩ DEG)
# ══════════════════════════════════════════════════════════════
tbl_A7 <- read_tsv(file.path(out_dir, "6_3_gene_table_A7.tsv"), show_col_types = FALSE)
tbl_A9 <- read_tsv(file.path(out_dir, "6_3_gene_table_A9.tsv"), show_col_types = FALSE)

genes_A7 <- tbl_A7 %>% filter(is_interest) %>% pull(Gene)
genes_A9 <- tbl_A9 %>% filter(is_interest) %>% pull(Gene)

cat(sprintf("Genes de interés — A7: %d | A9: %d\n",
            length(genes_A7), length(genes_A9)))

# ══════════════════════════════════════════════════════════════
## 3) Función ORA (acepta cualquier ontología)
# ══════════════════════════════════════════════════════════════
run_ora <- function(genes, universe, label, ont) {
  cat(sprintf("  [%s] %s ... ", label, ont))

  ego <- enrichGO(
    gene          = genes,
    universe      = universe,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.20,
    readable      = TRUE
  )

  if (is.null(ego) || nrow(ego) == 0) {
    cat("sin términos significativos.\n")
    return(data.frame())
  }

  cat(sprintf("%d términos\n", nrow(ego)))

  as.data.frame(ego) %>%
    mutate(
      ontology      = ont,
      group         = label,
      GeneRatio_num = vapply(GeneRatio, function(x) {
        p <- as.numeric(strsplit(x, "/")[[1]])
        p[1] / p[2]
      }, numeric(1))
    )
}

# ══════════════════════════════════════════════════════════════
## 4) Ejecutar ORA para BP, MF, CC
# ══════════════════════════════════════════════════════════════
cat("\n── HPV-A7 ──────────────────────────────────────────────\n")
ora_A7 <- bind_rows(
  run_ora(genes_A7, universe, "HPV-A7", "BP"),
  run_ora(genes_A7, universe, "HPV-A7", "MF"),
  run_ora(genes_A7, universe, "HPV-A7", "CC")
)

cat("\n── HPV-A9 ──────────────────────────────────────────────\n")
ora_A9 <- bind_rows(
  run_ora(genes_A9, universe, "HPV-A9", "BP"),
  run_ora(genes_A9, universe, "HPV-A9", "MF"),
  run_ora(genes_A9, universe, "HPV-A9", "CC")
)

write_tsv(ora_A7, file.path(out_dir, "6_4_ORA_A7.tsv"))
write_tsv(ora_A9, file.path(out_dir, "6_4_ORA_A9.tsv"))
cat(sprintf("\nResultados — A7: %d | A9: %d términos totales\n",
            nrow(ora_A7), nrow(ora_A9)))

# ══════════════════════════════════════════════════════════════
## 5) Función de plot — reutilizable para ambas corridas
# ══════════════════════════════════════════════════════════════
TOP_N     <- 15
GROUP_COL <- c("HPV-A7" = "#1A6FB5", "HPV-A9" = "#C0392B")

ONT_META <- list(
  BP = list(label = "Biological Process",  lo = "#EDE0F5", hi = "#4A235A"),
  MF = list(label = "Molecular Function",  lo = "#D5F5E3", hi = "#1A5C38"),
  CC = list(label = "Cellular Component",  lo = "#FADBD8", hi = "#7B241C")
)

make_panel <- function(grp, accent_col, df_ont, col_lo, col_hi) {
  df <- df_ont %>% filter(group == grp)
  if (nrow(df) == 0) return(NULL)

  ggplot(df, aes(x = GeneRatio_num, y = term_ord, fill = log_padj)) +
    geom_col(width = 0.68, colour = NA) +
    scale_y_reordered() +
    scale_fill_gradient(
      low  = col_lo, high = col_hi,
      name = expression(bold(-log[10](p[adj])))
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.06))) +
    labs(title = grp, x = "Gene Ratio", y = NULL) +
    theme_classic(base_size = 16) +
    theme(
      plot.title         = element_text(face = "bold", size = 19,
                                        colour = accent_col),
      axis.text.y        = element_text(size = 14, colour = "grey10",
                                        face = "bold"),
      axis.text.x        = element_text(size = 13, face = "bold",
                                        colour = "grey20"),
      axis.title.x       = element_text(size = 14, face = "bold",
                                        colour = "grey30"),
      axis.line          = element_line(colour = "grey70"),
      legend.title       = element_text(size = 13, face = "bold"),
      legend.text        = element_text(size = 12, face = "bold"),
      panel.grid.major.x = element_line(colour = "grey93", linewidth = 0.35)
    )
}

ora_all <- bind_rows(ora_A7, ora_A9) %>%
  filter(!is.na(p.adjust)) %>%
  mutate(
    group     = factor(group, levels = c("HPV-A7", "HPV-A9")),
    log_padj  = -log10(p.adjust),
    Desc_wrap = str_wrap(Description, width = 44)
  )

for (ont in c("BP", "MF", "CC")) {
  meta  <- ONT_META[[ont]]
  fname <- sprintf("6_4_ORA_%s_barplot", ont)

  df_ont <- ora_all %>%
    filter(ontology == ont) %>%
    group_by(group) %>%
    slice_max(GeneRatio_num, n = TOP_N, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(term_ord = reorder_within(Desc_wrap, GeneRatio_num, group))

  if (nrow(df_ont) == 0) {
    cat(sprintf("  %s: sin términos, omitido.\n", ont)); next
  }

  p_A7 <- make_panel("HPV-A7", GROUP_COL["HPV-A7"], df_ont, meta$lo, meta$hi)
  p_A9 <- make_panel("HPV-A9", GROUP_COL["HPV-A9"], df_ont, meta$lo, meta$hi)
  panels <- Filter(Negate(is.null), list(p_A7, p_A9))

  combined <- wrap_plots(panels, ncol = 1) +
    plot_annotation(
      title   = sprintf("%s — Hub ∩ DE genes | Cervical Cancer (CESC)", meta$label),
      caption = sprintf(
        "Universe: %d genes · BH-adjusted p < 0.05 · Top %d terms by Gene Ratio",
        length(universe), TOP_N
      ),
      theme = theme(
        plot.title   = element_text(face = "bold", hjust = 0.5, size = 17),
        plot.caption = element_text(colour = "grey50", size = 11, hjust = 0.5)
      )
    )

  n_grp  <- df_ont %>% count(group) %>% pull(n) %>% max(na.rm = TRUE)
  plot_h <- n_grp * 0.33 * 2 + 3
  plot_w <- 9

  ggsave(file.path(out_dir, paste0(fname, ".pdf")),
         plot = combined, width = plot_w, height = plot_h, device = cairo_pdf)
  ggsave(file.path(out_dir, paste0(fname, ".png")),
         plot = combined, width = plot_w, height = plot_h, dpi = 1200)

  cat(sprintf("  Guardado: %s.pdf/png\n", fname))
}

cat("✓ 6_4 completo.\n")
