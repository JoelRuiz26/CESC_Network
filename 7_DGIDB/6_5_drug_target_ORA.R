############################################################
## 6_5_drug_target_ORA.R
##
## Functional enrichment of target genes from proximity-
## significant drugs (z < -1.96) in HPV-A7 and HPV-A9 networks.
##
## Gene groups:
##   Shared     — targets of drugs significant in BOTH cohorts (n=6 drugs)
##   A7 only    — targets of drugs significant only in HPV-A7
##   A9 only    — targets of drugs significant only in HPV-A9
##
## Methods: GO BP  +  KEGG (MSigDB C2:KEGG_LEGACY)
##
## Universe: all 2,800 genes present in the DGIdb bipartite network
##
## Outputs (in 7_DGIDB/):
##   6_5_drugORA_GO_BP.pdf / .png / .tsv
##   6_5_drugORA_KEGG.pdf  / .png / .tsv
############################################################

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(patchwork)
})

setwd("~/CESC_Network/7_DGIDB")
out_dir <- "~/CESC_Network/7_DGIDB"

reorder_within <- function(x, by, within, ...) {
  new_x <- paste(x, within, sep = "___")
  stats::reorder(new_x, by, ...)
}
scale_y_reordered <- function(..., sep = "___") {
  ggplot2::scale_y_discrete(
    labels = function(x) sub(paste0(sep, ".*"), "", x), ...)
}

# ══════════════════════════════════════════════════════════
## 1) Load data & define gene groups
# ══════════════════════════════════════════════════════════
edges   <- read_tsv("6_1_edges_bipartite.tsv", show_col_types = FALSE)
prox_A7 <- read_tsv("6_3_proximity_A7.tsv",    show_col_types = FALSE)
prox_A9 <- read_tsv("6_3_proximity_A9.tsv",    show_col_types = FALSE)

universe_genes <- unique(edges$gene_name)
cat(sprintf("Universe (bipartite network targets): %d genes\n", length(universe_genes)))

sig_A7   <- prox_A7 %>% filter(significant) %>% pull(drug)
sig_A9   <- prox_A9 %>% filter(significant) %>% pull(drug)
both_sig <- intersect(sig_A7, sig_A9)
only_A7  <- setdiff(sig_A7, sig_A9)
only_A9  <- setdiff(sig_A9, sig_A7)

get_targets <- function(drugs) {
  edges %>% filter(drug_name %in% drugs) %>% pull(gene_name) %>% unique()
}

GROUPS <- list(
  "Shared"   = get_targets(both_sig),
  "A7 only"  = get_targets(only_A7),
  "A9 only"  = get_targets(only_A9)
)

cat(sprintf("Drug groups  — Shared: %d drugs (%d targets) | A7 only: %d drugs (%d targets) | A9 only: %d drugs (%d targets)\n",
            length(both_sig), length(GROUPS[["Shared"]]),
            length(only_A7),  length(GROUPS[["A7 only"]]),
            length(only_A9),  length(GROUPS[["A9 only"]])))

# ══════════════════════════════════════════════════════════
## 2) ORA functions
# ══════════════════════════════════════════════════════════
run_go_bp <- function(genes, universe, label) {
  cat(sprintf("  [GO BP] %s ... ", label))
  ego <- enrichGO(
    gene          = genes,
    universe      = universe,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.20,
    readable      = TRUE
  )
  if (is.null(ego) || nrow(ego) == 0) { cat("ns\n"); return(data.frame()) }
  cat(sprintf("%d terms\n", nrow(ego)))
  as.data.frame(ego) %>%
    mutate(
      group         = label,
      GeneRatio_num = vapply(GeneRatio, function(x) {
        p <- as.numeric(strsplit(x, "/")[[1]]); p[1] / p[2]
      }, numeric(1)),
      log_padj      = -log10(p.adjust)
    )
}

run_kegg_msig <- function(genes, universe, label, t2g) {
  cat(sprintf("  [KEGG]  %s ... ", label))
  ego <- enricher(
    gene          = genes,
    universe      = universe,
    TERM2GENE     = t2g,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.20
  )
  if (is.null(ego) || nrow(ego) == 0) { cat("ns\n"); return(data.frame()) }
  cat(sprintf("%d terms\n", nrow(ego)))
  as.data.frame(ego) %>%
    mutate(
      group         = label,
      Description   = str_replace_all(ID, "_", " ") %>%
                      str_to_title() %>%
                      str_remove("^Kegg "),
      GeneRatio_num = vapply(GeneRatio, function(x) {
        p <- as.numeric(strsplit(x, "/")[[1]]); p[1] / p[2]
      }, numeric(1)),
      log_padj      = -log10(p.adjust)
    )
}

# ══════════════════════════════════════════════════════════
## 3) Run enrichments
# ══════════════════════════════════════════════════════════
cat("\n── GO Biological Process ────────────────────────────\n")
go_results <- bind_rows(lapply(names(GROUPS), function(nm)
  run_go_bp(GROUPS[[nm]], universe_genes, nm)
))

cat("\n── KEGG Pathways (MSigDB C2) ────────────────────────\n")
t2g_kegg <- msigdbr(species = "Homo sapiens",
                    collection = "C2", subcollection = "CP:KEGG_LEGACY") %>%
  select(gs_name, gene_symbol) %>% distinct()

kegg_results <- bind_rows(lapply(names(GROUPS), function(nm)
  run_kegg_msig(GROUPS[[nm]], universe_genes, nm, t2g_kegg)
))

write_tsv(go_results,   file.path(out_dir, "6_5_drugORA_GO_BP.tsv"))
write_tsv(kegg_results, file.path(out_dir, "6_5_drugORA_KEGG.tsv"))

# ══════════════════════════════════════════════════════════
## 4) Plot
# ══════════════════════════════════════════════════════════
TOP_N <- 12

GROUP_COLS <- c(
  "Shared"  = "#6C3483",
  "A7 only" = "#1A6FB5",
  "A9 only" = "#C0392B"
)

make_ora_plot <- function(df, main_title, col_lo, col_hi,
                          caption_extra = "") {
  if (nrow(df) == 0) return(NULL)

  df_plot <- df %>%
    filter(!is.na(p.adjust)) %>%
    mutate(group = factor(group, levels = names(GROUP_COLS))) %>%
    group_by(group) %>%
    slice_max(GeneRatio_num, n = TOP_N, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      Desc_wrap = str_wrap(Description, width = 42),
      term_ord  = reorder_within(Desc_wrap, GeneRatio_num, group)
    )

  present_groups <- levels(droplevels(df_plot$group))

  panels <- lapply(present_groups, function(grp) {
    d <- df_plot %>% filter(group == grp)
    if (nrow(d) == 0) return(NULL)
    ggplot(d, aes(x = GeneRatio_num, y = term_ord, fill = log_padj)) +
      geom_col(width = 0.68, colour = NA) +
      scale_y_reordered() +
      scale_fill_gradient(
        low  = col_lo, high = col_hi,
        name = expression(bold(-log[10](p[adj])))
      ) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
      labs(title = grp, x = "Gene Ratio", y = NULL) +
      theme_classic(base_size = 14) +
      theme(
        plot.title         = element_text(face = "bold", size = 15,
                                          colour = GROUP_COLS[grp]),
        axis.text.y        = element_text(size = 12, colour = "grey10", face = "bold"),
        axis.text.x        = element_text(size = 11, face = "bold", colour = "grey20"),
        axis.title.x       = element_text(size = 12, face = "bold", colour = "grey30"),
        axis.line          = element_line(colour = "grey70"),
        legend.title       = element_text(size = 11, face = "bold"),
        legend.text        = element_text(size = 10, face = "bold"),
        panel.grid.major.x = element_line(colour = "grey93", linewidth = 0.35)
      )
  })
  panels <- Filter(Negate(is.null), panels)
  if (length(panels) == 0) return(NULL)

  n_terms <- df_plot %>% count(group) %>% pull(n) %>% max(na.rm = TRUE)
  plot_h  <- n_terms * 0.40 * length(panels) + 3.5

  wrap_plots(panels, ncol = 1) +
    plot_annotation(
      title   = main_title,
      caption = sprintf(
        "Universe: %d drug-target genes (DGIdb bipartite network) · BH-adjusted p < 0.05 · Top %d terms by Gene Ratio%s",
        length(universe_genes), TOP_N, caption_extra
      ),
      theme = theme(
        plot.title   = element_text(face = "bold", hjust = 0.5, size = 16),
        plot.caption = element_text(colour = "grey50", size = 10, hjust = 0.5)
      )
    )
}

save_plot <- function(p, fname, w = 10) {
  if (is.null(p)) { cat(sprintf("  %s: no significant terms\n", fname)); return(invisible()) }
  df_tmp  <- tryCatch(ggplot_build(p)$data[[1]], error = function(e) data.frame())
  plot_h  <- attr(p, "plot_h")
  h <- if (is.null(plot_h)) 14 else plot_h

  ggsave(file.path(out_dir, paste0(fname, ".pdf")),
         plot = p, width = w, height = h, device = cairo_pdf, limitsize = FALSE)
  ggsave(file.path(out_dir, paste0(fname, ".png")),
         plot = p, width = w, height = h, dpi = 300, limitsize = FALSE)
  cat(sprintf("  Saved: %s.pdf/png\n", fname))
}

# ── GO BP ─────────────────────────────────────────────────
cat("\n── Plotting GO BP ───────────────────────────────────\n")
p_go <- make_ora_plot(
  go_results,
  "GO Biological Process — Target Genes of Proximity-Significant Drugs",
  "#D6EAF8", "#1A5276"
)
if (!is.null(p_go)) {
  n_grps <- n_distinct(filter(go_results, !is.na(p.adjust))$group)
  n_t    <- go_results %>% group_by(group) %>% slice_max(GeneRatio_num, n=TOP_N) %>%
            count(group) %>% pull(n) %>% max(na.rm=TRUE)
  h_go   <- n_t * 0.40 * n_grps + 3.5
  ggsave(file.path(out_dir, "6_5_drugORA_GO_BP.pdf"),
         p_go, width = 10, height = h_go, device = cairo_pdf, limitsize = FALSE)
  ggsave(file.path(out_dir, "6_5_drugORA_GO_BP.png"),
         p_go, width = 10, height = h_go, dpi = 300, limitsize = FALSE)
  cat(sprintf("  Saved: 6_5_drugORA_GO_BP.pdf/png\n"))
}

# ── KEGG ──────────────────────────────────────────────────
cat("\n── Plotting KEGG ────────────────────────────────────\n")
p_kegg <- make_ora_plot(
  kegg_results,
  "KEGG Pathways — Target Genes of Proximity-Significant Drugs",
  "#D5F5E3", "#1A5C38"
)
if (!is.null(p_kegg)) {
  n_grps <- n_distinct(filter(kegg_results, !is.na(p.adjust))$group)
  n_t    <- kegg_results %>% group_by(group) %>% slice_max(GeneRatio_num, n=TOP_N) %>%
            count(group) %>% pull(n) %>% max(na.rm=TRUE)
  h_kegg <- n_t * 0.40 * n_grps + 3.5
  ggsave(file.path(out_dir, "6_5_drugORA_KEGG.pdf"),
         p_kegg, width = 10, height = h_kegg, device = cairo_pdf, limitsize = FALSE)
  ggsave(file.path(out_dir, "6_5_drugORA_KEGG.png"),
         p_kegg, width = 10, height = h_kegg, dpi = 300, limitsize = FALSE)
  cat(sprintf("  Saved: 6_5_drugORA_KEGG.pdf/png\n"))
}

cat("\n✓ 6_5 complete.\n")
