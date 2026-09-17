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
### ONE ROW PER GENE. A gene that is an essential hub in both clades is a
### single row carrying its degree and gene effect in each network, so the
### class totals are the row counts and nothing is double counted:
###   Both 11 | HPV-A7 only 36 | HPV-A9 only 43 | Common (pan-essential) 573
### `class` is the four-way partition used in the figures; pan-essential takes
### precedence, so the other three classes are the clade-selective dependencies.
###
###   6_12_essential_hub_supp_table.csv
### ============================================================
suppressPackageStartupMessages({
  library(tidyverse); library(vroom)
})
res_dir <- "~/CESC_Network/07_dependency/results"
tab_dir <- "~/CESC_Network/tables"

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

## the four classes used throughout the figures; pan-essential takes precedence,
## so the three remaining classes are exactly the clade-selective dependencies
cls_of <- function(g, net) {
  dplyr::case_when(g %in% pan  ~ "Common (pan-essential)",
                   g %in% both ~ "Both",
                   net == "HPV-A7" ~ "HPV-A7 only",
                   TRUE            ~ "HPV-A9 only")
}

long <- ess %>%
  transmute(network, gene, degree,
            mean_gene_effect = round(mean_dep, 4),
            class            = cls_of(gene, network))

supp <- long %>%
  pivot_wider(names_from = network,
              values_from = c(degree, mean_gene_effect),
              names_glue = "{.value}_{sub('HPV-', '', network)}") %>%
  transmute(gene,
            class,
            essential_hub_A7 = if_else(!is.na(degree_A7), "yes", "no"),
            essential_hub_A9 = if_else(!is.na(degree_A9), "yes", "no"),
            degree_A7, degree_A9,
            mean_gene_effect_A7, mean_gene_effect_A9) %>%
  arrange(factor(class, levels = c("Both", "HPV-A7 only", "HPV-A9 only",
                                   "Common (pan-essential)")), gene)

write.csv(supp, file.path(tab_dir, "Table_S2_essential_hubs.csv"), row.names = FALSE)

cat(sprintf("essential hubs: A7 = %d | A9 = %d | rows (= genes) = %d\n",
            sum(supp$essential_hub_A7 == "yes"),
            sum(supp$essential_hub_A9 == "yes"), nrow(supp)))
cat("genes per class:\n")
print(table(supp$class))
cat(sprintf("clade-selective (not pan-essential): %d genes\n",
            sum(supp$class != "Common (pan-essential)")))
message("saved 6_12_essential_hub_supp_table.csv to: ", res_dir)
