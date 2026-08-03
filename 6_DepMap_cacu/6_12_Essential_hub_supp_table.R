### ============================================================
### 6_12_Essential_hub_supp_table.R
###
### Supplementary table: the complete list of essential hub genes of
### BOTH networks, one row per (network, gene), with the clade degree
### and mean CRISPR gene effect, whether the gene is an essential hub in
### both clades, and whether it is a DepMap pan-essential (common) gene.
###
### Same data-driven thresholds as 6_8 / 6_9 / 6_10 / 6_11 (elbow rule):
###   HPV-A7  degree >= 35   gene effect <= -0.317   -> 395 essential hubs
###   HPV-A9  degree >= 42   gene effect <= -0.316   -> 436 essential hubs
###
###   6_12_essential_hub_supp_table.csv   (Excel-friendly)
### ============================================================
suppressPackageStartupMessages({
  library(tidyverse); library(vroom)
})
res_dir <- "~/CESC_Network/6_DepMap_cacu/results"

elbow <- function(v) {
  n <- length(v); x <- (seq_len(n) - 1) / (n - 1)
  y <- (v - v[1]) / (v[n] - v[1]); v[which.max(abs(y - x))]
}

tab <- vroom(file.path(res_dir, "6_2_gene_dep_degree_table.tsv"), show_col_types = FALSE) %>%
  filter(!is.na(degree), !is.na(mean_dep)) %>%
  distinct(network, gene, .keep_all = TRUE)

th <- tab %>% group_by(network) %>%
  summarise(k_cut = elbow(sort(degree, decreasing = TRUE)),
            d_cut = elbow(sort(mean_dep)), .groups = "drop")

pan <- vroom(file.path(res_dir, "6_1_common_essential_genes.tsv"),
             delim = "\t", show_col_types = FALSE)$gene_u %>% toupper()

ess <- tab %>% left_join(th, by = "network") %>%
  filter(degree >= k_cut, mean_dep <= d_cut) %>%
  mutate(gene = toupper(gene))

## essential in both clades = gene present in the essential set of A7 and A9
both <- intersect(ess$gene[ess$network == "HPV-A7"],
                  ess$gene[ess$network == "HPV-A9"])

supp <- ess %>%
  transmute(network,
            gene,
            degree,
            mean_gene_effect          = round(mean_dep, 4),
            essential_in_both_clades  = if_else(gene %in% both, "yes", "no"),
            pan_essential             = if_else(gene %in% pan,  "yes", "no")) %>%
  arrange(network, mean_gene_effect)                       # strongest dependency first

write.csv(supp, file.path(res_dir, "6_12_essential_hub_supp_table.csv"), row.names = FALSE)

cat(sprintf("essential hubs: A7 = %d | A9 = %d | rows = %d\n",
            sum(supp$network == "HPV-A7"), sum(supp$network == "HPV-A9"), nrow(supp)))
cat("pan-essential among essential hubs:\n")
print(table(network = supp$network, pan_essential = supp$pan_essential))
message("✅ saved 6_12_essential_hub_supp_table.csv to: ", res_dir)
