### ============================================================
### 01_consistent_selective_hubs_table.R  (EXPLORATORY / TEMPORARY)
###
### For each cohort (CESC, HTMCP-CC) separately, takes the "selective" hub
### set already used in Figure 4 (hub in HPV-A7 or HPV-A9, essential, with
### pan-essential genes removed) and classifies each gene as A7-selective,
### A9-selective, or essential-in-both. Keeps only genes whose category is
### IDENTICAL in both cohorts ("consistent" across cohorts) and reports them
### as a table grouped by clade (A7-consistent, A9-consistent) and Both.
###
### -> results/01_consistent_selective_hubs.tsv
### ============================================================
suppressPackageStartupMessages({
  library(tidyverse); library(vroom)
})

out_dir <- "~/CESC_Network/08_cross_cohort_validation"
dep_shared_dir <- "~/CESC_Network/07_dependency/results"

pan <- vroom(file.path(dep_shared_dir, "6_1_common_essential_genes.tsv"), delim = "\t", show_col_types = FALSE)$gene_u %>% toupper()

elbow <- function(v) {
  n <- length(v); x <- (seq_len(n) - 1) / (n - 1); y <- (v - v[1]) / (v[n] - v[1])
  list(idx = which.max(abs(y - x)), val = v[which.max(abs(y - x))])
}

## same "selective hub" classification as build_cohort() in figure_4_combined.R
classify_cohort <- function(res_dir) {
  tab <- vroom(file.path(res_dir, "6_2_gene_dep_degree_table.tsv"), show_col_types = FALSE) %>%
    filter(!is.na(degree), !is.na(mean_dep)) %>% distinct(network, gene, .keep_all = TRUE)
  nets <- sort(unique(tab$network))
  info <- map(nets, function(nk) {
    d <- tab %>% filter(network == nk)
    ek <- elbow(sort(d$degree, decreasing = TRUE)); ed <- elbow(sort(d$mean_dep))
    list(net = nk, ek = ek, ed = ed)
  }) %>% set_names(nets)

  flags <- tab %>% left_join(map_dfr(info, ~ tibble(network = .x$net, k_cut = .x$ek$val, d_cut = .x$ed$val)), by = "network") %>%
    mutate(hub = degree >= k_cut, ess = hub & mean_dep <= d_cut, net = str_remove(network, "HPV-")) %>%
    select(gene, net, hub, ess) %>% pivot_wider(names_from = net, values_from = c(hub, ess)) %>%
    mutate(across(starts_with("hub_") | starts_with("ess_"), ~ replace_na(.x, FALSE)))

  dep <- vroom(file.path(res_dir, "6_7_hub_dependency_medians.tsv"), show_col_types = FALSE) %>%
    select(gene, gene_u, dep_A7, dep_A9)

  dat <- dep %>% inner_join(flags, by = "gene") %>%
    filter(hub_A7 | hub_A9, !(toupper(gene_u) %in% pan)) %>%
    mutate(grp = case_when(ess_A7 & ess_A9 ~ "Both", ess_A7 ~ "A7", ess_A9 ~ "A9", TRUE ~ "None"))

  dat
}

cesc  <- classify_cohort("~/CESC_Network/07_dependency/results")
htmcp <- classify_cohort("~/CESC_Network/07_dependency/validation_uganda/results")

merged <- inner_join(
  cesc  %>% select(gene, gene_u, grp_cesc  = grp),
  htmcp %>% select(gene,          grp_htmcp = grp),
  by = "gene"
) %>%
  mutate(
    consistent_group = case_when(
      grp_cesc == "Both" & grp_htmcp == "Both" ~ "Both (A7 & A9)",
      grp_cesc == "A7"   & grp_htmcp == "A7"   ~ "HPV-A7",
      grp_cesc == "A9"   & grp_htmcp == "A9"   ~ "HPV-A9",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(consistent_group)) %>%
  select(gene, gene_u, consistent_group) %>%
  arrange(factor(consistent_group, levels = c("HPV-A7", "HPV-A9", "Both (A7 & A9)")), gene)

cat("Consistent selective hubs (same category in CESC AND HTMCP-CC):\n")
print(table(merged$consistent_group))
cat("\n")
merged %>% split(.$consistent_group) %>% iwalk(~ {
  cat(sprintf("== %s (n=%d) ==\n", .y, nrow(.x)))
  cat(paste(.x$gene_u, collapse = ", "), "\n\n")
})

write_tsv(merged, file.path(out_dir, "results/01_consistent_selective_hubs.tsv"))
message("Saved: results/01_consistent_selective_hubs.tsv")
