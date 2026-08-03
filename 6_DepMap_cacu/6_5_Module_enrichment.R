### ============================================================
### 6_5_Module_enrichment.R
###
### In which Louvain module(s) are the essential hub genes
### concentrated, per HPV group?
###
### Gene lists (4). Thresholds are data-driven, the elbow rule of 6_10 and the
### same cuts as 6_8 / 6_9 (degree elbow for hubs, gene-effect elbow for
### essentiality), applied per network:
###   A7 complete  = hub in A7 (degree >= 35) & gene effect(A7) <= -0.317  -> 395
###   A7 selective = A7 complete, pan-essential removed
###   A9 complete  = hub in A9 (degree >= 42) & gene effect(A9) <= -0.316  -> 436
###   A9 selective = A9 complete, pan-essential removed
###
### Test: per module, hypergeometric over-representation of the gene
###   list among all network genes (universe = all genes with a
###   community). BH-corrected across modules. This asks whether a
###   module holds MORE essential hubs than expected by chance.
###
### Modules tested are RESTRICTED to those with >= 1 significant
### MSigDB Hallmark term in the ORA of 7_Modularity/5_2_ORA
### (5_2_3_Hallmark{7,9}_louv_max.tsv), i.e. functionally annotated
### modules only. Background (universe) stays the full network.
###
### Outputs:
###   6_5_essential_hub_module_enrichment.tsv   — all modules/sets/tests
###   6_5_essential_hub_gene_lists.tsv          — the 4 gene lists w/ module
###   6_5_A_module_enrichment.png               — complete essential hubs
### ============================================================
suppressPackageStartupMessages({
  library(igraph)
  library(tidyverse)
  library(vroom)
  library(patchwork)
  library(ggforce)
})

dep_dir <- "~/CESC_Network/6_DepMap_cacu"
res_dir <- file.path(dep_dir, "results")           # outputs kept apart from scripts
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
mod_dir <- "~/CESC_Network/7_Modularity/5_1_Modules"
col_A7 <- "#003366"; col_A9 <- "#CC2200"

## data-driven thresholds, the elbow rule of 6_10 (identical to 6_8 / 6_9): the
## degree ranking gives the hub cut, the gene-effect ranking gives the essential
## cut, per network. Replaces the former top-10% degree + fixed -0.5 cutoff.
elbow <- function(v) {
  n <- length(v); x <- (seq_len(n) - 1) / (n - 1)
  y <- (v - v[1]) / (v[n] - v[1]); v[which.max(abs(y - x))]
}

## ── Data ──────────────────────────────────────────────────────────────────────
mod_tbl <- function(g) tibble(gene = V(g)$name,
                              community_id = V(g)$Community,
                              community_name = V(g)$Community_name)
g7 <- readRDS(file.path(mod_dir, "5_1_1_Louvain_A7_graph_max.rds"))
g9 <- readRDS(file.path(mod_dir, "5_1_1_Louvain_A9_graph_max.rds"))
m7 <- mod_tbl(g7); m9 <- mod_tbl(g9)

## degree + clade mean gene effect per gene per network (same table and same
## filter as 6_8/6_9/6_10), then the per-network elbow thresholds
gdd <- vroom(file.path(res_dir, "6_2_gene_dep_degree_table.tsv"), show_col_types = FALSE) %>%
  filter(!is.na(degree), !is.na(mean_dep)) %>%
  distinct(network, gene, .keep_all = TRUE)
th <- gdd %>% group_by(network) %>%
  summarise(k_cut = elbow(sort(degree, decreasing = TRUE)),
            d_cut = elbow(sort(mean_dep)), .groups = "drop")
message(sprintf("thresholds  A7: degree>=%g gene effect<=%.3f | A9: degree>=%g gene effect<=%.3f",
                th$k_cut[th$network == "HPV-A7"], th$d_cut[th$network == "HPV-A7"],
                th$k_cut[th$network == "HPV-A9"], th$d_cut[th$network == "HPV-A9"]))

pan <- vroom(file.path(res_dir, "6_1_common_essential_genes.tsv"),
             delim = "\t", show_col_types = FALSE)$gene_u %>% toupper()

## hub = degree >= elbow; essential hub (complete) = hub & gene effect <= elbow;
## selective = essential hub that is NOT pan-essential
flagged <- gdd %>% left_join(th, by = "network") %>%
  mutate(hub           = degree >= k_cut,
         essential     = hub & mean_dep <= d_cut,
         pan_essential = toupper(gene) %in% pan,
         selective     = essential & !pan_essential)

pick <- function(net, col) flagged$gene[flagged$network == net & flagged[[col]]]
h7 <- pick("HPV-A7", "hub");                  h9 <- pick("HPV-A9", "hub")
A7_complete  <- pick("HPV-A7", "essential");  A9_complete  <- pick("HPV-A9", "essential")
A7_selective <- pick("HPV-A7", "selective");  A9_selective <- pick("HPV-A9", "selective")
message(sprintf("A7: complete=%d selective=%d | A9: complete=%d selective=%d",
                length(A7_complete), length(A7_selective),
                length(A9_complete), length(A9_selective)))

## ── Per-hub-gene essentiality annotation table ────────────────────────────────
# One row per (network, hub gene), with the clade mean gene effect and three flags:
#   depmap_dependency : mean CRISPR gene effect in that clade's lines
#   essential         : 1 if gene effect <= dependency elbow (-0.317 / -0.316)
#   pan_essential     : 1 if in DepMap common_essential set
#   important         : 1 if essential & NOT pan-essential (the selective hubs)
hub_table <- flagged %>% filter(hub) %>%
  transmute(network, gene,
            depmap_dependency = mean_dep,
            essential     = as.integer(essential),
            pan_essential = as.integer(pan_essential),
            important     = as.integer(selective))
vroom_write(hub_table, file.path(res_dir, "6_5_hubgene_essentiality_table.tsv"), delim = "\t")

cat("\n=== Hub-gene table: column tallies by network (table()) ===\n")
cat(sprintf("\nTotal hub genes: A7 = %d | A9 = %d\n",
            sum(hub_table$network == "HPV-A7"), sum(hub_table$network == "HPV-A9")))
cat("\n[essential  (gene effect <= dependency elbow: -0.317 A7 / -0.316 A9)]\n")
print(table(network = hub_table$network, essential     = hub_table$essential))
cat("\n[pan_essential]\n")
print(table(network = hub_table$network, pan_essential = hub_table$pan_essential))
cat("\n[important  (essential & not pan-essential)]\n")
print(table(network = hub_table$network, important     = hub_table$important))

# save the 4 gene lists with module annotation
save_lists <- bind_rows(
  m7 %>% filter(gene %in% A7_complete)  %>% mutate(group = "HPV-A7", set = "complete"),
  m7 %>% filter(gene %in% A7_selective) %>% mutate(group = "HPV-A7", set = "selective"),
  m9 %>% filter(gene %in% A9_complete)  %>% mutate(group = "HPV-A9", set = "complete"),
  m9 %>% filter(gene %in% A9_selective) %>% mutate(group = "HPV-A9", set = "selective"))
vroom_write(save_lists, file.path(res_dir, "6_5_essential_hub_gene_lists.tsv"), delim = "\t")

## ── Functional filter: modules with >=1 MSigDB Hallmark term (ORA 5_2_3) ──────
# Keep only Louvain modules that were functionally annotated with at least one
# significant Hallmark term in 7_Modularity/5_2_ORA. Background stays the full
# network; we just restrict WHICH modules are tested + BH-corrected.
ora_dir <- "~/CESC_Network/7_Modularity/5_2_ORA"
hallmark_modules <- function(path) {
  vroom(path, show_col_types = FALSE) %>%
    distinct(Community_Name) %>% pull(Community_Name)
}
hm7 <- hallmark_modules(file.path(ora_dir, "5_2_3_Hallmark7_louv_max.tsv"))
hm9 <- hallmark_modules(file.path(ora_dir, "5_2_3_Hallmark9_louv_max.tsv"))
message(sprintf("Hallmark-annotated modules kept: A7 = %d | A9 = %d", length(hm7), length(hm9)))

## ── Hypergeometric per-module enrichment ──────────────────────────────────────
enrich_modules <- function(genes_set, m, group, set, keep_modules) {
  N <- nrow(m)
  K <- sum(m$gene %in% genes_set)
  m %>% group_by(community_id, community_name) %>%
    summarise(n_M = n(), k = sum(gene %in% genes_set), .groups = "drop") %>%
    filter(community_name %in% keep_modules) %>%    # only Hallmark-annotated modules
    mutate(expected = n_M * K / N,
           fold     = k / expected,
           p        = phyper(k - 1, K, N - K, n_M, lower.tail = FALSE),
           q        = p.adjust(p, "BH")) %>%
    mutate(group = group, set = set, K_total = K)
}

enr <- bind_rows(
  enrich_modules(A7_complete,  m7, "HPV-A7", "complete",  hm7),
  enrich_modules(A7_selective, m7, "HPV-A7", "selective", hm7),
  enrich_modules(A9_complete,  m9, "HPV-A9", "complete",  hm9),
  enrich_modules(A9_selective, m9, "HPV-A9", "selective", hm9))
vroom_write(enr, file.path(res_dir, "6_5_essential_hub_module_enrichment.tsv"), delim = "\t")

cat("\n=== Significantly enriched modules (BH q<0.05) ===\n")
enr %>% filter(q < 0.05, k > 0) %>%
  arrange(set, group, q) %>%
  select(set, group, module = community_name, n_M, k, expected, fold, q) %>%
  mutate(expected = round(expected, 1), fold = round(fold, 1),
         q = formatC(q, format = "e", digits = 1)) %>%
  as.data.frame() %>% print()

## ── Plot builder ──────────────────────────────────────────────────────────────
base_theme <- theme_bw(base_family = "Times New Roman") +
  theme(axis.text.y = element_text(size = 17, face = "bold"),
        axis.text.x = element_text(size = 18, face = "bold"),
        axis.title  = element_text(size = 22, face = "bold"),
        plot.title  = element_text(size = 24, face = "bold", hjust = 0.5),
        strip.text  = element_text(size = 21, face = "bold"),
        legend.text = element_text(size = 17, face = "bold"),
        legend.title = element_text(size = 18, face = "bold"),
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
    scale_size_continuous(name = "Essential\nhub genes", range = c(7, 20),
                          breaks = c(25, 100, 200)) +
    scale_x_continuous(limits = c(1, xmax), expand = expansion(mult = c(0, 0.02))) +
    labs(title = title,
         x = "Fold enrichment  (observed / expected)", y = NULL) +
    base_theme +
    theme(axis.text.y   = element_text(size = 19, face = "bold"),
          strip.text    = element_text(size = 21, face = "bold"),
          strip.background = element_rect(fill = "grey93", color = "grey60"),
          legend.position = "right",
          plot.subtitle = element_text(size = 13.5, face = "italic", hjust = 0.5, color = "grey30"))
}

pA <- make_enrich_plot("complete",
        "Modules enriched for essential hub genes")
ggsave(file.path(res_dir, "6_5_A_module_enrichment.png"),
       pA, width = 10.5, height = 8.2, dpi = 300, bg = "white")
ggsave(file.path(res_dir, "6_5_A_module_enrichment.pdf"),
       pA, width = 10.5, height = 8.2, device = cairo_pdf, bg = "white")
message("✅ 6_4 outputs saved to: ", dep_dir)
