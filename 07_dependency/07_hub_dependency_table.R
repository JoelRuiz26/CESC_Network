### ============================================================
### 6_7_Hub_dependency_table.R
###
### Builds the per-gene hub-dependency table used by the essential-hub
### figures (6_8 Venn, 6_9 full scatter, 6_10 selective scatter).
### For each hub gene: hub-sharing category (Shared / A7 / A9), mean and
### per-line MEDIAN CRISPR dependency in A7 vs A9 cervical lines, and
### essential-hub class. Pan-essential genes are KEPT here; they are only
### removed downstream (6_10) when the selective view is needed.
###
###   6_7_hub_dependency_medians.tsv — per-gene dependency table (cache)
### ============================================================
suppressPackageStartupMessages({
  library(depmap)
  library(igraph)
  library(tidyverse)
  library(vroom)
})

dep_dir <- "~/CESC_Network/07_dependency"
res_dir <- file.path(dep_dir, "results")           # outputs kept apart from scripts
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
mod_dir <- "~/CESC_Network/06_modules/results"
HUB_QUANTILE     <- 0.90
ESSENTIAL_CUTOFF <- -0.5

col_A7 <- "#003366"; col_A9 <- "#CC2200"; col_shared <- "#6A3D9A"
cat_colors <- c("Shared" = col_shared, "A7-specific" = col_A7, "A9-specific" = col_A9)
cat_levels <- c("Shared", "A7-specific", "A9-specific")

## grey = hub gene that is not an essential hub in either clade (background layer)
ess_levels <- c("Essential hub in both", "Essential hub in A7",
                "Essential hub in A9", "Other network genes")
ess_colors <- c("Essential hub in both" = col_shared,
                "Essential hub in A7"   = col_A7,
                "Essential hub in A9"   = col_A9,
                "Other network genes"   = "grey72")

## ── 1) Hub classification ─────────────────────────────────────────────────────
hub_set <- function(g) { d <- igraph::degree(g); names(d)[d >= quantile(d, HUB_QUANTILE)] }
g7 <- readRDS(file.path(mod_dir, "5_1_1_Louvain_A7_graph_max.rds"))
g9 <- readRDS(file.path(mod_dir, "5_1_1_Louvain_A9_graph_max.rds"))
h7 <- hub_set(g7); h9 <- hub_set(g9); shared <- intersect(h7, h9)
hub <- tibble(gene = union(h7, h9)) %>%
  mutate(category = factor(case_when(
    gene %in% shared ~ "Shared", gene %in% h7 ~ "A7-specific", TRUE ~ "A9-specific"),
    levels = cat_levels), gene_u = toupper(gene))

## ── 2) Per-line MEDIAN dependency (10 cervical lines) ─────────────────────────
lines  <- vroom(file.path(res_dir, "6_1_crispr_lines.tsv"), show_col_types = FALSE)
ids_A7 <- lines$depmap_id[lines$hpv_group == "A7"]
ids_A9 <- lines$depmap_id[lines$hpv_group == "A9"]

per_gene <- depmap_crispr() %>%
  filter(depmap_id %in% lines$depmap_id) %>%
  transmute(gene_u = toupper(gene_name), depmap_id, dep = dependency) %>%
  group_by(gene_u) %>%
  summarise(med_A7  = median(dep[depmap_id %in% ids_A7],  na.rm = TRUE),
            med_A9  = median(dep[depmap_id %in% ids_A9],  na.rm = TRUE),
            med_all = median(dep, na.rm = TRUE), .groups = "drop")

## Essential-hub membership, defined exactly as in 6_5: hub in that clade AND
## clade mean dependency <= -0.5. Pan-essential genes are KEPT.
## The axes use the same clade mean, so the dashed -0.5 lines separate the
## coloured points from the grey ones instead of contradicting them.
d7 <- vroom(file.path(res_dir, "6_1_cervix_dep_mean_A7.tsv"), show_col_types = FALSE)
d9 <- vroom(file.path(res_dir, "6_1_cervix_dep_mean_A9.tsv"), show_col_types = FALSE)
ess7 <- d7$gene_u[d7$mean_dep <= ESSENTIAL_CUTOFF]
ess9 <- d9$gene_u[d9$mean_dep <= ESSENTIAL_CUTOFF]

ess_A7 <- toupper(h7)[toupper(h7) %in% ess7]     # 265
ess_A9 <- toupper(h9)[toupper(h9) %in% ess9]     # 337

## Grey background = EVERY gene of either network, not only the hubs. Restricting
## it to hubs made a grey point in the lower-left quadrant impossible: any hub
## with both dependencies <= -0.5 is essential in the clade where it is a hub,
## so it always got coloured.
all_genes <- tibble(gene = union(V(g7)$name, V(g9)$name)) %>%
  mutate(gene_u = toupper(gene)) %>%
  anti_join(hub, by = "gene_u")

dat <- bind_rows(hub, all_genes) %>%
  left_join(d7 %>% transmute(gene_u, dep_A7 = mean_dep), by = "gene_u") %>%
  left_join(d9 %>% transmute(gene_u, dep_A9 = mean_dep), by = "gene_u") %>%
  left_join(per_gene, by = "gene_u") %>%
  filter(!is.na(dep_A7), !is.na(dep_A9)) %>%
  mutate(ess_class = factor(case_when(
           gene_u %in% ess_A7 & gene_u %in% ess_A9 ~ "Essential hub in both",
           gene_u %in% ess_A7                      ~ "Essential hub in A7",
           gene_u %in% ess_A9                      ~ "Essential hub in A9",
           TRUE                                    ~ "Other network genes"),
           levels = ess_levels),
         clade_bias = dep_A9 - dep_A7)   # <0 more A9-essential, >0 more A7-essential
vroom_write(dat, file.path(res_dir, "6_7_hub_dependency_medians.tsv"), delim = "\t")
message("Hub genes tabulated (pan-essential kept): ", nrow(dat))
print(as.data.frame(count(dat, ess_class)))
message("\u2705 6_7 hub-dependency table saved to: ", res_dir)
