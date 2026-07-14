### ============================================================
### 6_3_Hub_essentiality_ORA.R
###
### PRIMARY question: do network hub genes carry stronger DepMap
### dependency than the rest of the network?
###
### Statistical approach (threshold-free on dependency):
###   * Wilcoxon rank-sum (Mann-Whitney U) on the FULL score
###     distribution — no arbitrary essentiality cutoff.
###   * Cliff's delta as the effect size (rank-based, not inflated
###     by the large n; interpretable as P(hub<non) - P(hub>non)).
###
### The ONLY tunable parameter left is the hub definition
### (top 10% by degree). Its robustness is examined separately in
### 6_4_Hub_threshold_sensitivity.R.
###
### A categorical Fisher's exact test at the conventional -0.5
### essentiality threshold is reported ONLY as a confirmatory,
### supplementary anchor (not the headline result).
###
### Hub genes : top 10% by degree k (quantile 0.90)
###
### Outputs:
###   6_3_stats_hub_vs_nonhub.tsv    — Wilcoxon/Cliff (primary) + Fisher (suppl.)
###   6_3_A_dependency_ecdf.png      — ECDF hub vs non-hub (MAIN figure)
###   6_3_ORA_hubs_A7_GO.tsv / A9    — full GO ORA of hub genes
###   6_3_ORA_hubs_A7_GO_simplified.tsv / A9
###   6_3_B_ORA_hub_dotplot.png      — GO ORA dotplot (hub characterization)
### ============================================================
suppressPackageStartupMessages({
  library(igraph)
  library(tidyverse)
  library(vroom)
  library(patchwork)
  library(scales)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(GOSemSim)
  library(enrichplot)
})

dep_dir <- "~/CESC_Network/6_DepMap_cacu"
mod_dir <- "~/CESC_Network/7_Modularity/5_1_Modules"

HUB_QUANTILE     <- 0.90
ESSENTIAL_CUTOFF <- -0.5   # used ONLY for the supplementary Fisher anchor

net_colors <- c("HPV-A7" = "#003366", "HPV-A9" = "#CC2200")

## ── Load data ─────────────────────────────────────────────────────────────────
gene_tbl <- vroom(file.path(dep_dir, "6_2_gene_dep_degree_table.tsv"),
                  show_col_types = FALSE)

g7 <- readRDS(file.path(mod_dir, "5_1_1_Louvain_A7_graph_max.rds"))
g9 <- readRDS(file.path(mod_dir, "5_1_1_Louvain_A9_graph_max.rds"))
universe_A7 <- V(g7)$name
universe_A9 <- V(g9)$name

## ── Classify genes into hub / non-hub ─────────────────────────────────────────
classify <- function(df) {
  thr <- quantile(df$degree, HUB_QUANTILE)
  df %>% mutate(
    gene_class = factor(
      if_else(degree >= thr, "Hub genes", "Non-hub genes"),
      levels = c("Hub genes", "Non-hub genes")
    ),
    hub_thr = thr
  )
}
g7_cl <- classify(filter(gene_tbl, network == "HPV-A7"))
g9_cl <- classify(filter(gene_tbl, network == "HPV-A9"))

## ── Effect size: Cliff's delta via Mann-Whitney U ─────────────────────────────
# delta = 2U/(n1*n2) - 1 ; W from wilcox.test already weights ties by 0.5.
# Two-sided test (conservative / non-directional).
cliff_magnitude <- function(d) case_when(
  abs(d) < 0.147 ~ "negligible",
  abs(d) < 0.330 ~ "small",
  abs(d) < 0.474 ~ "medium",
  TRUE           ~ "large"
)

run_stats <- function(df, label) {
  hub <- df$mean_dep[df$gene_class == "Hub genes"]
  non <- df$mean_dep[df$gene_class == "Non-hub genes"]
  n1 <- length(hub); n2 <- length(non)

  # Primary: Wilcoxon + Cliff's delta (threshold-free)
  w  <- wilcox.test(hub, non, exact = FALSE)
  U  <- as.numeric(w$statistic)
  delta <- (2 * U / (n1 * n2)) - 1

  # Supplementary confirmatory: Fisher at -0.5
  ft <- fisher.test(table(
    Hub       = df$gene_class == "Hub genes",
    Essential = df$mean_dep <= ESSENTIAL_CUTOFF
  ))

  tibble(
    network        = label,
    n_hub          = n1,
    n_non          = n2,
    hub_threshold_k = unique(df$hub_thr),
    median_hub     = median(hub),
    median_non     = median(non),
    # primary
    cliff_delta    = delta,
    cliff_mag      = cliff_magnitude(delta),
    wilcoxon_p     = w$p.value,
    # supplementary (categorical, at -0.5)
    pct_ess_hub    = mean(hub <= ESSENTIAL_CUTOFF) * 100,
    pct_ess_non    = mean(non <= ESSENTIAL_CUTOFF) * 100,
    fisher_OR      = unname(ft$estimate),
    fisher_CI_lo   = ft$conf.int[1],
    fisher_CI_hi   = ft$conf.int[2],
    fisher_p       = ft$p.value
  )
}

stats_res <- bind_rows(run_stats(g7_cl, "HPV-A7"),
                       run_stats(g9_cl, "HPV-A9"))

vroom_write(stats_res, file.path(dep_dir, "6_3_stats_hub_vs_nonhub.tsv"), delim = "\t")
message("[OK] Stats saved.")
print(stats_res %>% dplyr::select(network, n_hub, median_hub, median_non,
                                   cliff_delta, cliff_mag, wilcoxon_p,
                                   fisher_OR, fisher_p))

## ── Theme ─────────────────────────────────────────────────────────────────────
base_theme <- theme_bw(base_family = "Times New Roman") +
  theme(
    axis.text        = element_text(size = 17, face = "bold"),
    axis.title       = element_text(size = 20, face = "bold"),
    plot.title       = element_text(size = 23, face = "bold", hjust = 0.5),
    strip.text       = element_text(size = 21, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", linewidth = 0.8),
    legend.position  = "none"
  )

## ── Plot A (MAIN): ECDF hub vs non-hub ────────────────────────────────────────
# The ECDF is exactly what Wilcoxon / Cliff's delta summarize; the area
# between the two curves is proportional to the effect size.
# Format a p-value for on-plot annotation (Times-friendly, no Cliff's delta).
fmt_p <- function(p) {
  if (p < 2.2e-16) return("p < 2.2e-16")
  if (p < 0.001)   return(sprintf("p = %.1e", p))
  sprintf("p = %.3f", p)
}

make_ecdf_panel <- function(df, net, pval) {
  hub_col <- unname(net_colors[net])
  hub <- df$mean_dep[df$gene_class == "Hub genes"]
  non <- df$mean_dep[df$gene_class == "Non-hub genes"]

  xs <- seq(-3, 0.5, by = 0.01)
  ec <- bind_rows(
    tibble(score = xs, prop = ecdf(hub)(xs) * 100, gene_class = "Hub genes"),
    tibble(score = xs, prop = ecdf(non)(xs) * 100, gene_class = "Non-hub genes")
  )
  ribbon <- tibble(score = xs, hub = ecdf(hub)(xs) * 100, non = ecdf(non)(xs) * 100)

  ggplot() +
    geom_ribbon(data = ribbon, aes(y = score, xmin = non, xmax = hub),
                fill = hub_col, alpha = 0.16, orientation = "y") +
    geom_hline(yintercept = -0.5, linetype = "dashed",
               color = "grey35", linewidth = 0.7) +
    geom_line(data = filter(ec, gene_class == "Non-hub genes"),
              aes(x = prop, y = score), color = "grey55", linewidth = 1.4) +
    geom_line(data = filter(ec, gene_class == "Hub genes"),
              aes(x = prop, y = score), color = hub_col, linewidth = 1.7) +
    annotate("segment", x = 5, xend = 11, y = 0.30, yend = 0.30,
             color = hub_col, linewidth = 1.7) +
    annotate("text", x = 13, y = 0.30, label = "Hub genes", hjust = 0, size = 5.2,
             fontface = "bold", family = "Times New Roman", color = hub_col) +
    annotate("segment", x = 5, xend = 11, y = 0.13, yend = 0.13,
             color = "grey55", linewidth = 1.4) +
    annotate("text", x = 13, y = 0.13, label = "Non-hub genes", hjust = 0, size = 5.2,
             fontface = "bold", family = "Times New Roman", color = "grey45") +
    annotate("text", x = 55, y = -2.7,
             label = paste0("Mann–Whitney U test\n", fmt_p(pval)),
             hjust = 0, vjust = 0, size = 4.8, fontface = "bold",
             family = "Times New Roman", color = "black") +
    scale_y_continuous(limits = c(-3, 0.5), breaks = seq(-3, 0, by = 1)) +
    scale_x_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%"),
                       expand = expansion(mult = c(0.01, 0.03))) +
    labs(title = net, x = "Cumulative % of genes",
         y = "DepMap CRISPR dependency score") +
    base_theme
}

p_a7 <- stats_res$wilcoxon_p[stats_res$network == "HPV-A7"]
p_a9 <- stats_res$wilcoxon_p[stats_res$network == "HPV-A9"]

pA <- make_ecdf_panel(g7_cl, "HPV-A7", p_a7) +
  make_ecdf_panel(g9_cl, "HPV-A9", p_a9) +
  plot_annotation(
    title = "Cumulative distribution of DepMap dependency scores",
    theme = theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5,
                                            family = "Times New Roman",
                                            margin = margin(b = 10)))
  )

ggsave(file.path(dep_dir, "6_3_A_dependency_ecdf.png"),
       pA, width = 16, height = 8, dpi = 300, bg = "white")
message("[OK] Plot A (ECDF) saved.")

## ── GO ORA on hub genes (biological characterization) ─────────────────────────
run_hub_ORA <- function(hub_genes, universe_genes, label) {
  message("GO ORA: ", label, "  (", length(hub_genes), " hub genes)...")
  to_entrez <- function(syms) {
    bitr(syms, fromType = "SYMBOL", toType = "ENTREZID",
         OrgDb = org.Hs.eg.db, drop = TRUE)$ENTREZID
  }
  enrichGO(
    gene          = to_entrez(hub_genes),
    universe      = to_entrez(universe_genes),
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
}

simplify_safe <- function(ora) {
  tryCatch(simplify(ora, cutoff = 0.7, by = "p.adjust", select_fun = min),
           error = function(e) ora)
}

hubs_A7 <- g7_cl$gene[g7_cl$gene_class == "Hub genes"]
hubs_A9 <- g9_cl$gene[g9_cl$gene_class == "Hub genes"]

ora_A7 <- run_hub_ORA(hubs_A7, universe_A7, "HPV-A7")
ora_A9 <- run_hub_ORA(hubs_A9, universe_A9, "HPV-A9")

ora_A7_simp <- simplify_safe(ora_A7)
ora_A9_simp <- simplify_safe(ora_A9)

save_ora <- function(ora, path, label) {
  if (!is.null(ora) && nrow(ora@result) > 0) {
    as.data.frame(ora) %>% vroom_write(path, delim = "\t")
    message(label, ": ", nrow(ora@result), " terms")
    TRUE
  } else { message(label, ": no significant terms"); FALSE }
}

has_A7 <- save_ora(ora_A7,      file.path(dep_dir, "6_3_ORA_hubs_A7_GO.tsv"),            "A7 full")
has_A9 <- save_ora(ora_A9,      file.path(dep_dir, "6_3_ORA_hubs_A9_GO.tsv"),            "A9 full")
          save_ora(ora_A7_simp, file.path(dep_dir, "6_3_ORA_hubs_A7_GO_simplified.tsv"), "A7 simplified")
          save_ora(ora_A9_simp, file.path(dep_dir, "6_3_ORA_hubs_A9_GO_simplified.tsv"), "A9 simplified")

## ── Plot B: ORA dotplot (simplified, top 20) ──────────────────────────────────
make_ora_dotplot <- function(ora_res, color, title, top_n = 20) {
  df <- as.data.frame(ora_res) %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n) %>%
    mutate(
      Description   = str_wrap(Description, width = 48),
      Description   = factor(Description, levels = rev(Description)),
      GeneRatio_num = map_dbl(GeneRatio, ~ {
        p <- strsplit(.x, "/")[[1]]
        as.numeric(p[1]) / as.numeric(p[2])
      })
    )

  ggplot(df, aes(x = GeneRatio_num, y = Description,
                 color = p.adjust, size = Count)) +
    geom_point(alpha = 0.88) +
    scale_color_gradient(low = color, high = "grey78", name = "p.adjust",
                         labels = scientific_format(digits = 2)) +
    scale_size_continuous(name = "Gene count", range = c(4, 12)) +
    scale_x_continuous(labels = number_format(accuracy = 0.01),
                       expand = expansion(mult = c(0.05, 0.15))) +
    labs(title = title, x = "Gene Ratio", y = NULL) +
    base_theme +
    theme(axis.text.y     = element_text(size = 13, face = "bold"),
          legend.position = "right",
          legend.title    = element_text(size = 15, face = "bold"),
          legend.text     = element_text(size = 13, face = "bold"))
}

if (has_A7 || has_A9) {
  plots_b <- list()
  if (has_A7) plots_b[["A7"]] <- make_ora_dotplot(ora_A7_simp, "#003366",
                                                    "HPV-A7  (top 20, simplified)")
  if (has_A9) plots_b[["A9"]] <- make_ora_dotplot(ora_A9_simp, "#CC2200",
                                                    "HPV-A9  (top 20, simplified)")

  pB <- wrap_plots(plots_b, ncol = 2) +
    plot_annotation(
      title = "GO Biological Process enrichment in hub genes",
      theme = theme(plot.title = element_text(size = 30, face = "bold", hjust = 0.5,
                                              family = "Times New Roman"))
    )

  ggsave(file.path(dep_dir, "6_3_B_ORA_hub_dotplot.png"),
         pB, width = 30, height = 16, dpi = 300, bg = "white")
  message("[OK] Plot B (ORA dotplot) saved.")
}

message("✅ All 6_3 outputs saved to: ", dep_dir)
