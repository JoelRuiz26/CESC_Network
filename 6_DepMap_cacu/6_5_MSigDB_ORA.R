############################################################
## 6_5_MSigDB_ORA.R
##
## Over-Representation Analysis con MSigDB (via msigdbr)
## sobre los genes de interés (hub ∩ DEG) de A7 y A9.
##
## Colecciones:
##   H          — Hallmark gene sets (50 firmas de cáncer)
##   C2:KEGG    — KEGG pathways
##   C2:REACTOME— Reactome pathways
##
## Universo: 2_Prepro_TCGA_GTEx/2_1_Universe_annotation_table.rds
##
## Outputs (en 6_DepMap_cacu/):
##   6_5_MSigDB_[H|KEGG|REACTOME]_barplot.pdf / .png
##   6_5_MSigDB_A7.tsv / 6_5_MSigDB_A9.tsv
############################################################

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(msigdbr)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(patchwork)
})

out_dir <- "~/CESC_Network/6_DepMap_cacu"

# Inline reorder_within / scale_y_reordered (sin tidytext)
reorder_within <- function(x, by, within, ...) {
  new_x <- paste(x, within, sep = "___")
  stats::reorder(new_x, by, ...)
}
scale_y_reordered <- function(..., sep = "___") {
  ggplot2::scale_y_discrete(
    labels = function(x) sub(paste0(sep, ".*"), "", x), ...)
}

# ══════════════════════════════════════════════════════════════
## 1) Universo
# ══════════════════════════════════════════════════════════════
universe_sym <- readRDS(
  "~/CESC_Network/2_Prepro_TCGA_GTEx/2_1_Universe_annotation_table.rds"
) %>%
  pull(hgnc_symbol) %>%
  na.omit() %>%
  unique()

cat(sprintf("Universo: %d genes\n", length(universe_sym)))

# ══════════════════════════════════════════════════════════════
## 2) Genes de interés
# ══════════════════════════════════════════════════════════════
tbl_A7 <- read_tsv(file.path(out_dir, "6_3_gene_table_A7.tsv"), show_col_types = FALSE)
tbl_A9 <- read_tsv(file.path(out_dir, "6_3_gene_table_A9.tsv"), show_col_types = FALSE)

genes_A7 <- tbl_A7 %>% filter(is_interest) %>% pull(Gene)
genes_A9 <- tbl_A9 %>% filter(is_interest) %>% pull(Gene)

cat(sprintf("Genes de interés — A7: %d | A9: %d\n",
            length(genes_A7), length(genes_A9)))

# ══════════════════════════════════════════════════════════════
## 3) Colecciones MSigDB
# ══════════════════════════════════════════════════════════════
COLLECTIONS <- list(
  H = list(
    label    = "Hallmark Gene Sets",
    col      = "H",
    subcol   = NULL,
    lo       = "#EDE0F5", hi = "#4A235A",   # morado
    file     = "6_5_MSigDB_H_barplot"
  ),
  KEGG = list(
    label    = "KEGG Pathways",
    col      = "C2",
    subcol   = "CP:KEGG_LEGACY",
    lo       = "#D5F5E3", hi = "#1A5C38",   # verde
    file     = "6_5_MSigDB_KEGG_barplot"
  ),
  REACTOME = list(
    label    = "Reactome Pathways",
    col      = "C2",
    subcol   = "CP:REACTOME",
    lo       = "#EDE0F5", hi = "#4A235A",   # morado
    file     = "6_5_MSigDB_REACTOME_barplot"
  )
)

# Pre-cargar gene sets de MSigDB (especie humana, por símbolo)
get_term2gene <- function(col, subcol = NULL) {
  args <- list(species = "Homo sapiens", collection = col)
  if (!is.null(subcol)) args$subcollection <- subcol
  do.call(msigdbr, args) %>%
    select(gs_name, gene_symbol) %>%
    distinct()
}

# ══════════════════════════════════════════════════════════════
## 4) Función ORA con enricher()
# ══════════════════════════════════════════════════════════════
run_msigdb_ora <- function(genes, universe, label, term2gene) {
  cat(sprintf("  [%s] ... ", label))

  ego <- enricher(
    gene          = genes,
    universe      = universe,
    TERM2GENE     = term2gene,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.20
  )

  if (is.null(ego) || nrow(ego) == 0) {
    cat("sin términos significativos.\n")
    return(data.frame())
  }

  cat(sprintf("%d términos\n", nrow(ego)))

  as.data.frame(ego) %>%
    mutate(
      group         = label,
      GeneRatio_num = vapply(GeneRatio, function(x) {
        p <- as.numeric(strsplit(x, "/")[[1]])
        p[1] / p[2]
      }, numeric(1)),
      Description = str_replace_all(ID, "_", " ") %>%
                    str_to_title() %>%
                    str_remove("^Hallmark |^Kegg |^Reactome ")
    )
}

# ══════════════════════════════════════════════════════════════
## 5) Plot function (igual que 6_4)
# ══════════════════════════════════════════════════════════════
TOP_N     <- 15
GROUP_COL <- c("HPV-A7" = "#1A6FB5", "HPV-A9" = "#C0392B")

make_panel <- function(grp, accent_col, df_grp, col_lo, col_hi) {
  if (nrow(df_grp) == 0) return(NULL)

  ggplot(df_grp, aes(x = GeneRatio_num, y = term_ord, fill = log_padj)) +
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

# ══════════════════════════════════════════════════════════════
## 6) Ejecutar por colección y guardar plots
# ══════════════════════════════════════════════════════════════
all_results <- list()

for (col_id in names(COLLECTIONS)) {
  meta <- COLLECTIONS[[col_id]]
  cat(sprintf("\n══ %s ══════════════════════════════════════════\n", meta$label))

  t2g <- get_term2gene(meta$col, meta$subcol)
  cat(sprintf("  Gene sets en colección: %d\n", n_distinct(t2g$gs_name)))

  ora_A7 <- run_msigdb_ora(genes_A7, universe_sym, "HPV-A7", t2g)
  ora_A9 <- run_msigdb_ora(genes_A9, universe_sym, "HPV-A9", t2g)

  all_results[[col_id]] <- list(A7 = ora_A7, A9 = ora_A9)

  ora_comb <- bind_rows(ora_A7, ora_A9) %>%
    filter(!is.na(p.adjust)) %>%
    mutate(
      group     = factor(group, levels = c("HPV-A7", "HPV-A9")),
      log_padj  = -log10(p.adjust),
      Desc_wrap = str_wrap(Description, width = 44)
    )

  if (nrow(ora_comb) == 0) {
    cat(sprintf("  Sin términos significativos — omitido.\n")); next
  }

  plot_df <- ora_comb %>%
    group_by(group) %>%
    slice_max(GeneRatio_num, n = TOP_N, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(term_ord = reorder_within(Desc_wrap, GeneRatio_num, group))

  p_A7 <- make_panel("HPV-A7", GROUP_COL["HPV-A7"],
                     plot_df %>% filter(group == "HPV-A7"),
                     meta$lo, meta$hi)
  p_A9 <- make_panel("HPV-A9", GROUP_COL["HPV-A9"],
                     plot_df %>% filter(group == "HPV-A9"),
                     meta$lo, meta$hi)
  panels <- Filter(Negate(is.null), list(p_A7, p_A9))

  combined <- wrap_plots(panels, ncol = 1) +
    plot_annotation(
      title   = sprintf("%s — Hub ∩ DE genes | Cervical Cancer (CESC)", meta$label),
      caption = sprintf(
        "MSigDB · Universe: %d genes · BH-adjusted p < 0.05 · Top %d terms by Gene Ratio",
        length(universe_sym), TOP_N
      ),
      theme = theme(
        plot.title   = element_text(face = "bold", hjust = 0.5, size = 17),
        plot.caption = element_text(colour = "grey50", size = 11, hjust = 0.5)
      )
    )

  n_grp  <- plot_df %>% count(group) %>% pull(n) %>% max(na.rm = TRUE)
  plot_h <- n_grp * 0.33 * 2 + 3
  plot_w <- 9

  ggsave(file.path(out_dir, paste0(meta$file, ".pdf")),
         plot = combined, width = plot_w, height = plot_h, device = cairo_pdf)
  ggsave(file.path(out_dir, paste0(meta$file, ".png")),
         plot = combined, width = plot_w, height = plot_h, dpi = 1200)

  cat(sprintf("  Guardado: %s.pdf/png\n", meta$file))
}

# ── Guardar tablas completas ──────────────────────────────────
write_tsv(
  bind_rows(lapply(all_results, `[[`, "A7"), .id = "collection"),
  file.path(out_dir, "6_5_MSigDB_A7.tsv")
)
write_tsv(
  bind_rows(lapply(all_results, `[[`, "A9"), .id = "collection"),
  file.path(out_dir, "6_5_MSigDB_A9.tsv")
)

cat("\n✓ 6_5 completo.\n")
