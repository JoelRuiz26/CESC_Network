### ============================================================
### 6_7_Essential_hub_ORA.R
###
### GO (Biological Process) over-representation analysis of the
### SELECTIVE essential hub genes of each network (pan-essential
### REMOVED), to find clade-relevant, non-housekeeping processes.
###
###   A7 selective essential hubs = hub in A7 network (top 10% degree)
###        AND essential in A7 lines (mean dep <= -0.5) AND NOT pan-essential
###   A9 selective essential hubs = analogous for A9
###
### (The all-essential version was dominated by housekeeping: ~90% of
###  essential hubs are DepMap common-essential. Removing them isolates
###  the selective signal. n is small: ~49 A7 / ~44 A9 — low power.)
###
### Universe = all genes of the respective network.
###
### Outputs:
###   6_7_ORA_selective_hubs_A7_GO.tsv / _simplified.tsv
###   6_7_ORA_selective_hubs_A9_GO.tsv / _simplified.tsv
###   6_7_selective_hub_ORA_dotplot.png   (only if terms are significant)
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
})

dep_dir <- "~/CESC_Network/6_DepMap_cacu"
mod_dir <- "~/CESC_Network/7_Modularity/5_1_Modules"
HUB_QUANTILE     <- 0.90
ESSENTIAL_CUTOFF <- -0.5
col_A7 <- "#003366"; col_A9 <- "#CC2200"

## ── Hub sets + essentiality + pan-essential removal ───────────────────────────
hub_set <- function(g) { d <- igraph::degree(g); names(d)[d >= quantile(d, HUB_QUANTILE)] }
g7 <- readRDS(file.path(mod_dir, "5_1_1_Louvain_A7_graph_max.rds"))
g9 <- readRDS(file.path(mod_dir, "5_1_1_Louvain_A9_graph_max.rds"))
h7 <- hub_set(g7); h9 <- hub_set(g9)
universe_A7 <- V(g7)$name; universe_A9 <- V(g9)$name

d7  <- vroom(file.path(dep_dir, "6_1_cervix_dep_mean_A7.tsv"), show_col_types = FALSE)
d9  <- vroom(file.path(dep_dir, "6_1_cervix_dep_mean_A9.tsv"), show_col_types = FALSE)
pan <- vroom(file.path(dep_dir, "6_5_common_essential_genes.tsv"),
             delim = "\t", show_col_types = FALSE)$gene_u

ess_A7_u <- d7$gene_u[d7$mean_dep <= ESSENTIAL_CUTOFF]
ess_A9_u <- d9$gene_u[d9$mean_dep <= ESSENTIAL_CUTOFF]

# SELECTIVE = essential hub AND NOT pan-essential
sel_hub_A7 <- h7[toupper(h7) %in% ess_A7_u & !(toupper(h7) %in% pan)]
sel_hub_A9 <- h9[toupper(h9) %in% ess_A9_u & !(toupper(h9) %in% pan)]
message(sprintf("Selective essential hubs (pan-essential removed): A7=%d | A9=%d",
                length(sel_hub_A7), length(sel_hub_A9)))

## ── ORA ───────────────────────────────────────────────────────────────────────
to_entrez <- function(syms)
  bitr(syms, "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)$ENTREZID

run_ORA <- function(genes, universe, label) {
  message("GO ORA: ", label, " (", length(genes), " genes)...")
  enrichGO(gene = to_entrez(genes), universe = to_entrez(universe),
           OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
           pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
}
simplify_safe <- function(o) tryCatch(
  simplify(o, cutoff = 0.7, by = "p.adjust", select_fun = min), error = function(e) o)

ora_A7 <- run_ORA(sel_hub_A7, universe_A7, "HPV-A7")
ora_A9 <- run_ORA(sel_hub_A9, universe_A9, "HPV-A9")
ora_A7_s <- simplify_safe(ora_A7)
ora_A9_s <- simplify_safe(ora_A9)

save_ora <- function(o, path, label) {
  df <- if (!is.null(o)) as.data.frame(o) else data.frame()   # significant terms only
  if (nrow(df) > 0) {
    vroom_write(df, path, delim = "\t")
    message(label, ": ", nrow(df), " significant terms"); TRUE
  } else { message(label, ": no significant terms"); FALSE }
}
hA7 <- save_ora(ora_A7,   file.path(dep_dir, "6_7_ORA_selective_hubs_A7_GO.tsv"),            "A7 full")
hA9 <- save_ora(ora_A9,   file.path(dep_dir, "6_7_ORA_selective_hubs_A9_GO.tsv"),            "A9 full")
       save_ora(ora_A7_s, file.path(dep_dir, "6_7_ORA_selective_hubs_A7_GO_simplified.tsv"), "A7 simplified")
       save_ora(ora_A9_s, file.path(dep_dir, "6_7_ORA_selective_hubs_A9_GO_simplified.tsv"), "A9 simplified")

## ── Dotplot ───────────────────────────────────────────────────────────────────
base_theme <- theme_bw(base_family = "Times New Roman") +
  theme(axis.text.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 15, face = "bold"),
        axis.title  = element_text(size = 18, face = "bold"),
        plot.title  = element_text(size = 21, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text  = element_text(size = 12, face = "bold"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 0.8))

make_dotplot <- function(o, color, title, top_n = 20) {
  df <- as.data.frame(o) %>% arrange(p.adjust) %>% slice_head(n = top_n) %>%
    mutate(Description = str_wrap(Description, 46),
           Description = factor(Description, levels = rev(Description)),
           GeneRatio_num = map_dbl(GeneRatio, ~ {
             p <- strsplit(.x, "/")[[1]]; as.numeric(p[1]) / as.numeric(p[2]) }))
  ggplot(df, aes(GeneRatio_num, Description, color = p.adjust, size = Count)) +
    geom_point(alpha = 0.9) +
    scale_color_gradient(low = color, high = "grey80", name = "p.adjust",
                         labels = scientific_format(digits = 2)) +
    scale_size_continuous(name = "Count", range = c(4, 11)) +
    scale_x_continuous(labels = number_format(accuracy = 0.01),
                       expand = expansion(mult = c(0.05, 0.18))) +
    labs(title = title, x = "Gene Ratio", y = NULL) +
    base_theme + theme(legend.position = "right")
}

panels <- list()
if (hA7) panels[["A7"]] <- make_dotplot(ora_A7_s, col_A7,
           sprintf("HPV-A7 selective essential hubs (n=%d)", length(sel_hub_A7)))
if (hA9) panels[["A9"]] <- make_dotplot(ora_A9_s, col_A9,
           sprintf("HPV-A9 selective essential hubs (n=%d)", length(sel_hub_A9)))

if (length(panels) > 0) {
  p <- wrap_plots(panels, ncol = length(panels)) +
    plot_annotation(title = "GO enrichment of SELECTIVE essential hubs (pan-essential removed)",
      theme = theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5,
                    family = "Times New Roman")))
  ggsave(file.path(dep_dir, "6_7_selective_hub_ORA_dotplot.png"), p,
         width = 13 * length(panels), height = 11, dpi = 300, bg = "white")
  message("[OK] Dotplot saved.")
} else {
  message("No significant terms in either network — no dotplot produced.")
}
message("✅ 6_7 outputs saved to: ", dep_dir)
