### ============================================================
### 6_2_Network_DepMap_integration.R
### Join Louvain network degree with DepMap dependency scores.
###
### Output:
###   6_2_gene_dep_degree_table.tsv — per-gene joined table (A7 + A9),
###     the input for 6_3_Hub_dependency_CDF.R
### ============================================================
suppressPackageStartupMessages({
  library(igraph)
  library(tidyverse)
  library(vroom)
})

mod_dir <- "~/CESC_Network/7_Modularity/5_1_Modules"
dep_dir <- "~/CESC_Network/6_DepMap_cacu"
res_dir <- file.path(dep_dir, "results")           # outputs kept apart from scripts
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

g7 <- readRDS(file.path(mod_dir, "5_1_1_Louvain_A7_graph_max.rds"))
g9 <- readRDS(file.path(mod_dir, "5_1_1_Louvain_A9_graph_max.rds"))

dep_A7 <- vroom(file.path(res_dir, "6_1_cervix_dep_mean_A7.tsv"), show_col_types = FALSE)
dep_A9 <- vroom(file.path(res_dir, "6_1_cervix_dep_mean_A9.tsv"), show_col_types = FALSE)

## ── Per-gene table ─────────────────────────────────────────────────────────────
# toupper() fixes case mismatches (e.g. "C5orf22" vs "C5ORF22" in DepMap)
extract_gene_tbl <- function(g, dep_tbl, label) {
  deg <- igraph::degree(g)
  igraph::as_data_frame(g, what = "vertices") %>%
    transmute(
      gene           = name,
      gene_join      = toupper(name),
      community_id   = as.integer(Community),
      community_name = Community_name,
      community_size = as.integer(Community_size),
      degree         = as.integer(deg[name])
    ) %>%
    left_join(dep_tbl, by = c("gene_join" = "gene_u")) %>%
    dplyr::select(-gene_join) %>%
    filter(!is.na(mean_dep)) %>%
    mutate(network = label)
}

gene_A7 <- extract_gene_tbl(g7, dep_A7, "HPV-A7")
gene_A9 <- extract_gene_tbl(g9, dep_A9, "HPV-A9")

message("Genes con dato: A7 = ", nrow(gene_A7), " | A9 = ", nrow(gene_A9))

bind_rows(gene_A7, gene_A9) %>%
  vroom_write(file.path(res_dir, "6_2_gene_dep_degree_table.tsv"), delim = "\t")
message("[OK] Table saved.")
