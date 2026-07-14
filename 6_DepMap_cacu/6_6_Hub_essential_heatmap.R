### ============================================================
### 6_6_Hub_essential_heatmap.R
###
### Per-cell-line dependency heatmap of SELECTIVE essential hub genes
### (pan-essential removed), split by hub sharing category.
###
### Shows, in a single figure:
###   * shared hub genes that are also essential (essential across
###     both clade blocks)
###   * clade-specific hub genes (essential mainly in their own clade)
###
### Gene set: hub (top 10% degree in A7 and/or A9), NOT DepMap
###   common-essential, and essential (mean dep <= -0.5) in >=1 clade.
### Rows grouped by: Shared / A7-specific / A9-specific.
### Columns grouped by clade: 6 A7 lines (HPV-18/45) | 4 A9 lines (HPV-16).
###
### Output:
###   6_6_hub_essential_heatmap.png
###   6_6_hub_essential_heatmap.tsv
### ============================================================
suppressPackageStartupMessages({
  library(depmap)
  library(igraph)
  library(tidyverse)
  library(vroom)
  library(scales)
})

dep_dir <- "~/CESC_Network/6_DepMap_cacu"
mod_dir <- "~/CESC_Network/7_Modularity/5_1_Modules"
HUB_QUANTILE     <- 0.90
ESSENTIAL_CUTOFF <- -0.5

## ── 1) Hub classification ─────────────────────────────────────────────────────
hub_set <- function(g) { d <- igraph::degree(g); names(d)[d >= quantile(d, HUB_QUANTILE)] }
g7 <- readRDS(file.path(mod_dir, "5_1_1_Louvain_A7_graph_max.rds"))
g9 <- readRDS(file.path(mod_dir, "5_1_1_Louvain_A9_graph_max.rds"))
h7 <- hub_set(g7); h9 <- hub_set(g9); shared <- intersect(h7, h9)

hub <- tibble(gene = union(h7, h9)) %>%
  mutate(category = factor(case_when(
    gene %in% shared ~ "Shared",
    gene %in% h7     ~ "A7-specific",
    TRUE             ~ "A9-specific"),
    levels = c("Shared", "A7-specific", "A9-specific")),
    gene_u = toupper(gene))

## ── 2) Select selective essential hubs (pan-essential removed) ────────────────
d7  <- vroom(file.path(dep_dir, "6_1_cervix_dep_mean_A7.tsv"), show_col_types = FALSE) %>% rename(dep_A7 = mean_dep)
d9  <- vroom(file.path(dep_dir, "6_1_cervix_dep_mean_A9.tsv"), show_col_types = FALSE) %>% rename(dep_A9 = mean_dep)
pan <- vroom(file.path(dep_dir, "6_5_common_essential_genes.tsv"), delim = "\t", show_col_types = FALSE)$gene_u

sel_ess <- hub %>%
  left_join(d7, by = "gene_u") %>% left_join(d9, by = "gene_u") %>%
  filter(!(gene_u %in% pan), !is.na(dep_A7), !is.na(dep_A9),
         dep_A7 <= ESSENTIAL_CUTOFF | dep_A9 <= ESSENTIAL_CUTOFF) %>%
  mutate(mean_dep = (dep_A7 + dep_A9) / 2)
message("Selective essential hubs: ", nrow(sel_ess))

## ── 3) Per-line dependency for the 10 cervical lines ──────────────────────────
lines <- vroom(file.path(dep_dir, "6_1_crispr_lines.tsv"), show_col_types = FALSE) %>%
  mutate(cell = str_remove(cell_line, "_CERVIX"))
line_order <- lines %>% arrange(hpv_group, cell)

per_line <- depmap_crispr() %>%
  filter(depmap_id %in% lines$depmap_id) %>%
  transmute(gene_u = toupper(gene_name), depmap_id, dep = dependency) %>%
  inner_join(sel_ess %>% select(gene, gene_u, category, mean_dep), by = "gene_u") %>%
  left_join(lines %>% select(depmap_id, cell, hpv_group), by = "depmap_id")

## row order (most essential at top within each category) and column order
per_line <- per_line %>%
  mutate(gene = fct_reorder(gene, -mean_dep),
         cell = factor(cell, levels = line_order$cell),
         clade = factor(paste0(hpv_group, " lines"),
                        levels = c("A7 lines", "A9 lines")))

vroom_write(per_line, file.path(dep_dir, "6_6_hub_essential_heatmap.tsv"), delim = "\t")

## ── 4) Heatmap ────────────────────────────────────────────────────────────────
hm <- ggplot(per_line, aes(x = cell, y = gene, fill = dep)) +
  geom_tile(color = "grey85", linewidth = 0.25) +
  facet_grid(category ~ clade, scales = "free", space = "free", switch = "y") +
  scale_fill_gradient2(
    low = "#08306B", mid = "white", high = "#7F0000", midpoint = 0,
    limits = c(-2, 0.5), oob = squish,
    name = "Dependency\nscore",
    breaks = c(-2, -1, ESSENTIAL_CUTOFF, 0),
    labels = c("-2", "-1", "-0.5", "0")) +
  labs(title = "Selective essential hub genes across cervical cell lines",
       subtitle = "Hub genes (top 10% degree), pan-essential removed, essential in ≥ 1 clade",
       x = NULL, y = NULL) +
  theme_minimal(base_family = "Times New Roman") +
  theme(
    plot.title    = element_text(size = 24, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 15, face = "italic", hjust = 0.5, color = "grey30",
                                 margin = margin(b = 8)),
    axis.text.x   = element_text(size = 13, face = "bold", angle = 45, hjust = 1),
    axis.text.y   = element_text(size = 8,  face = "bold"),
    strip.text.x  = element_text(size = 16, face = "bold"),
    strip.text.y.left = element_text(size = 15, face = "bold", angle = 0),
    strip.background = element_rect(fill = "grey92", color = "grey60"),
    legend.title  = element_text(size = 14, face = "bold"),
    legend.text   = element_text(size = 12, face = "bold"),
    panel.grid    = element_blank(),
    panel.spacing = unit(4, "pt"))

ggsave(file.path(dep_dir, "6_6_hub_essential_heatmap.png"),
       hm, width = 11, height = 20, dpi = 300, bg = "white")
message("[OK] Heatmap saved.")

cat("\n=== genes per category ===\n")
print(as.data.frame(count(sel_ess, category)))
message("✅ 6_6 heatmap saved to: ", dep_dir)
