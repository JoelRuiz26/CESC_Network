### ============================================================
### 6_8_Essential_hub_module_enrichment.R
###
### In which Louvain module(s) are the essential hub genes
### concentrated, per HPV group?
###
### Gene lists (4):
###   A7 complete  = hub in A7 (top 10% degree) & mean dep(A7) <= -0.5
###   A7 selective = A7 complete, pan-essential removed
###   A9 complete  = hub in A9 & mean dep(A9) <= -0.5
###   A9 selective = A9 complete, pan-essential removed
###
### Test: per module, hypergeometric over-representation of the gene
###   list among all network genes (universe = all genes with a
###   community). BH-corrected across modules. This asks whether a
###   module holds MORE essential hubs than expected by chance.
###
### Outputs:
###   6_8_essential_hub_module_enrichment.tsv   — all modules/sets/tests
###   6_8_essential_hub_gene_lists.tsv          — the 4 gene lists w/ module
###   6_8_A_module_enrichment_complete.png      — complete essential hubs
###   6_8_B_module_enrichment_selective.png     — pan-essential removed
### ============================================================
suppressPackageStartupMessages({
  library(igraph)
  library(tidyverse)
  library(vroom)
  library(patchwork)
  library(ggforce)
})

dep_dir <- "~/CESC_Network/6_DepMap_cacu"
mod_dir <- "~/CESC_Network/7_Modularity/5_1_Modules"
HUB_QUANTILE     <- 0.90
ESSENTIAL_CUTOFF <- -0.5
col_A7 <- "#003366"; col_A9 <- "#CC2200"

## ── Data ──────────────────────────────────────────────────────────────────────
hub_set <- function(g) { d <- igraph::degree(g); names(d)[d >= quantile(d, HUB_QUANTILE)] }
mod_tbl <- function(g) tibble(gene = V(g)$name,
                              community_id = V(g)$Community,
                              community_name = V(g)$Community_name)
g7 <- readRDS(file.path(mod_dir, "5_1_1_Louvain_A7_graph_max.rds"))
g9 <- readRDS(file.path(mod_dir, "5_1_1_Louvain_A9_graph_max.rds"))
h7 <- hub_set(g7); h9 <- hub_set(g9)
m7 <- mod_tbl(g7); m9 <- mod_tbl(g9)

d7  <- vroom(file.path(dep_dir, "6_1_cervix_dep_mean_A7.tsv"), show_col_types = FALSE)
d9  <- vroom(file.path(dep_dir, "6_1_cervix_dep_mean_A9.tsv"), show_col_types = FALSE)
pan <- vroom(file.path(dep_dir, "6_5_common_essential_genes.tsv"),
             delim = "\t", show_col_types = FALSE)$gene_u
ess7 <- d7$gene_u[d7$mean_dep <= ESSENTIAL_CUTOFF]
ess9 <- d9$gene_u[d9$mean_dep <= ESSENTIAL_CUTOFF]

## ── The 4 gene lists ──────────────────────────────────────────────────────────
A7_complete  <- h7[toupper(h7) %in% ess7]
A7_selective <- A7_complete[!(toupper(A7_complete) %in% pan)]
A9_complete  <- h9[toupper(h9) %in% ess9]
A9_selective <- A9_complete[!(toupper(A9_complete) %in% pan)]
message(sprintf("A7: complete=%d selective=%d | A9: complete=%d selective=%d",
                length(A7_complete), length(A7_selective),
                length(A9_complete), length(A9_selective)))

# save the 4 gene lists with module annotation
save_lists <- bind_rows(
  m7 %>% filter(gene %in% A7_complete)  %>% mutate(group = "HPV-A7", set = "complete"),
  m7 %>% filter(gene %in% A7_selective) %>% mutate(group = "HPV-A7", set = "selective"),
  m9 %>% filter(gene %in% A9_complete)  %>% mutate(group = "HPV-A9", set = "complete"),
  m9 %>% filter(gene %in% A9_selective) %>% mutate(group = "HPV-A9", set = "selective"))
vroom_write(save_lists, file.path(dep_dir, "6_8_essential_hub_gene_lists.tsv"), delim = "\t")

## ── Hypergeometric per-module enrichment ──────────────────────────────────────
enrich_modules <- function(genes_set, m, group, set) {
  N <- nrow(m)
  K <- sum(m$gene %in% genes_set)
  m %>% group_by(community_id, community_name) %>%
    summarise(n_M = n(), k = sum(gene %in% genes_set), .groups = "drop") %>%
    mutate(expected = n_M * K / N,
           fold     = k / expected,
           p        = phyper(k - 1, K, N - K, n_M, lower.tail = FALSE)) %>%
    group_by(1) %>% mutate(q = p.adjust(p, "BH")) %>% ungroup() %>% select(-`1`) %>%
    mutate(group = group, set = set, K_total = K)
}

enr <- bind_rows(
  enrich_modules(A7_complete,  m7, "HPV-A7", "complete"),
  enrich_modules(A7_selective, m7, "HPV-A7", "selective"),
  enrich_modules(A9_complete,  m9, "HPV-A9", "complete"),
  enrich_modules(A9_selective, m9, "HPV-A9", "selective"))
vroom_write(enr, file.path(dep_dir, "6_8_essential_hub_module_enrichment.tsv"), delim = "\t")

cat("\n=== Significantly enriched modules (BH q<0.05) ===\n")
enr %>% filter(q < 0.05, k > 0) %>%
  arrange(set, group, q) %>%
  select(set, group, module = community_name, n_M, k, expected, fold, q) %>%
  mutate(expected = round(expected, 1), fold = round(fold, 1),
         q = formatC(q, format = "e", digits = 1)) %>%
  as.data.frame() %>% print()

## ── Plot builder ──────────────────────────────────────────────────────────────
base_theme <- theme_bw(base_family = "Times New Roman") +
  theme(axis.text.y = element_text(size = 13, face = "bold"),
        axis.text.x = element_text(size = 15, face = "bold"),
        axis.title  = element_text(size = 19, face = "bold"),
        plot.title  = element_text(size = 22, face = "bold", hjust = 0.5),
        strip.text  = element_text(size = 19, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 15, face = "bold"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 0.8))

make_enrich_plot <- function(set_name, title) {
  df <- enr %>% filter(set == set_name, q < 0.05) %>%   # significant modules only
    mutate(mod   = paste0(community_name, "  (", n_M, ")"),
           group = factor(group, levels = c("HPV-A7", "HPV-A9"))) %>%
    arrange(group, fold) %>%
    mutate(mod = factor(mod, levels = unique(mod)))

  df <- df %>% mutate(neglogq = -log10(q))
  xmax <- max(df$fold) * 1.06
  ggplot(df, aes(x = fold, y = mod)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey55", linewidth = 0.8) +
    geom_segment(aes(x = 1, xend = fold, yend = mod), color = "grey78", linewidth = 1.1) +
    geom_point(aes(size = k, color = neglogq)) +
    facet_col(vars(group), scales = "free_y", space = "free") +
    scale_color_gradient(low = "#F4A582", high = "#67001F",
                         name = expression(bold(-log[10]~italic(q)))) +
    scale_size_continuous(name = "Essential\nhub genes", range = c(5, 13),
                          breaks = c(25, 35, 45)) +
    scale_x_continuous(limits = c(1, xmax), expand = expansion(mult = c(0, 0.02))) +
    labs(title = title,
         subtitle = "Significantly enriched modules (hypergeometric BH q<0.05)",
         x = "Fold enrichment  (observed / expected)", y = NULL) +
    base_theme +
    theme(axis.text.y   = element_text(size = 15, face = "bold"),
          strip.text    = element_text(size = 18, face = "bold"),
          strip.background = element_rect(fill = "grey93", color = "grey60"),
          legend.position = "right",
          plot.subtitle = element_text(size = 13.5, face = "italic", hjust = 0.5, color = "grey30"))
}

pA <- make_enrich_plot("complete",
        "Modules enriched for essential hub genes")
ggsave(file.path(dep_dir, "6_8_A_module_enrichment_complete.png"),
       pA, width = 8.5, height = 12.5, dpi = 300, bg = "white")
message("✅ 6_8 outputs saved to: ", dep_dir)
