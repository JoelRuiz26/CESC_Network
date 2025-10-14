library(vroom)
library(igraph)
library(tibble)
library(dplyr)
library(fgsea)

# ----------------------------
# 1. Load networks and communities
# ----------------------------
filtered_graph_A7 <- readRDS("~/CESC_Network/7_Modularity/7_1_Louvain/7_1_1_Louvain_A7_graph.rds")
filtered_graph_A9 <- readRDS("~/CESC_Network/7_Modularity/7_1_Louvain/7_1_1_Louvain_A9_graph.rds")

create_Infomap_df <- function(graph) {
  tibble(
    Community = as.integer(V(graph)$Community),
    Community_name = as.character(V(graph)$Community_name),
    size = as.integer(V(graph)$Community_size),
    gene = as.character(V(graph)$name)
  )
}

df_A7 <- create_Infomap_df(filtered_graph_A7)
df_A9 <- create_Infomap_df(filtered_graph_A9)

# ----------------------------
# 2. Load full logFC lists
# ----------------------------
dge_A7 <- vroom("~/CESC_Network/3_DGE/3_1_DGE_A7_vs_SNT_AllResults.tsv")
dge_A7 <- dge_A7 %>% filter(Gene %in% df_A7$gene)

dge_A9 <- vroom("~/CESC_Network/3_DGE/3_3_DGE_A9_vs_SNT_AllResults.tsv")
dge_A9 <- dge_A9 %>% filter(Gene %in% df_A9$gene)

# ----------------------------
# 3. Create rankings for GSEA
# ----------------------------
ranking_A7 <- dge_A7 %>%
  mutate(score = logFC * -log10(adj.P.Val)) %>%
  select(Gene, score) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  arrange(desc(score)) %>%
  deframe()

ranking_A9 <- dge_A9 %>%
  mutate(score = logFC * -log10(adj.P.Val)) %>%
  select(Gene, score) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  arrange(desc(score)) %>%
  deframe()

# ----------------------------
# 4. Create gene sets by community
# ----------------------------
pathways_A7 <- df_A7 %>%
  group_by(Community_name) %>%
  summarise(genes = list(unique(gene))) %>%
  deframe()

pathways_A9 <- df_A9 %>%
  group_by(Community_name) %>%
  summarise(genes = list(unique(gene))) %>%
  deframe()

# ----------------------------
# 5. Run GSEA
# ----------------------------
set.seed(123)
fgsea_res_A7 <- fgsea(
  pathways = pathways_A7,
  stats = ranking_A7)

fgsea_res_A9 <- fgsea(
  pathways = pathways_A9,
  stats = ranking_A9)

# ----------------------------
# 6. Save results
# ----------------------------
vroom_write(fgsea_res_A7, file = "~/CESC_Network/7_Modularity/7_4_GSEA_module/7_4_1_GSEA_Communities_A7_full.tsv")
vroom_write(fgsea_res_A9, file = "~/CESC_Network/7_Modularity/7_4_GSEA_module/7_4_1_GSEA_Communities_A9_full.tsv")

# ----------------------------
# 7. Quick visualization
# ----------------------------
library(ggplot2)
# ========================================
# 8. Combined Plot (A7 + A9 in one figure)
# ========================================
# Combine significant results from both fgsea results
combined_top <- bind_rows(
  fgsea_res_A7 %>% filter(padj < 0.05) %>% mutate(Condition = "A7"),
  fgsea_res_A9 %>% filter(padj < 0.05) %>% mutate(Condition = "A9")
)

# Check if data is available
if (nrow(combined_top) == 0) {
  message("⚠️ No significant communities (padj < 0.005) found in either A7 or A9!")
} else {
  ggplot(combined_top, aes(x = reorder(pathway, NES), y = NES, fill = Condition)) +
    geom_col(position = position_dodge()) +
    coord_flip() +
    labs(
      title = "Top Enriched Communities (A7 and A9 Combined)",
      x = "Community",
      y = "Normalized Enrichment Score (NES)"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.text.y = element_text(size = 8)
    )
}







# ========================================
# 9. Fisher's Exact Test (Over-representation Analysis)
# ========================================
library(purrr)

# Flag DE genes (adj.P.Val < 0.05)
dge_A7 <- dge_A7 %>%
  mutate(isDE = adj.P.Val < 0.05)
dge_A9 <- dge_A9 %>%
  mutate(isDE = adj.P.Val < 0.05)

# Join DE flag into community dataframes
df_A7 <- df_A7 %>%
  left_join(dge_A7 %>% select(Gene, isDE), by = c("gene" = "Gene"))

df_A9 <- df_A9 %>%
  left_join(dge_A9 %>% select(Gene, isDE), by = c("gene" = "Gene"))

# Define universe counts
all_genes_A7 <- unique(df_A7$gene)
all_DE_A7 <- df_A7 %>% filter(isDE) %>% pull(gene)

all_genes_A9 <- unique(df_A9$gene)
all_DE_A9 <- df_A9 %>% filter(isDE) %>% pull(gene)

# ----------------------------------------
# A7: Fisher Test with Monte Carlo simulation
# ----------------------------------------
fisher_tables_A7 <- df_A7 %>%
  group_by(Community_name) %>%
  summarise(
    DE_in_module = sum(isDE, na.rm = TRUE),
    NonDE_in_module = sum(!isDE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    DE_outside = length(all_DE_A7) - DE_in_module,
    NonDE_outside = length(all_genes_A7) - length(all_DE_A7) - NonDE_in_module
  )

set.seed(123)
fisher_res_A7 <- fisher_tables_A7 %>%
  mutate(
    pval = pmap_dbl(
      list(DE_in_module, NonDE_in_module, DE_outside, NonDE_outside),
      ~ {
        mat <- matrix(c(..1, ..2, ..3, ..4), nrow = 2)
        fisher.test(mat, simulate.p.value = TRUE, B = 1e5)$p.value
      }
    )
  ) %>%
  mutate(padj = p.adjust(pval, method = "BH"))

sig_A7 <- fisher_res_A7 %>% filter(padj < 0.05)
cat("✅ Number of significant modules in A7 (FDR < 0.05):", nrow(sig_A7), "\n")

# ----------------------------------------
# A9: Fisher Test with Monte Carlo simulation
# ----------------------------------------
fisher_tables_A9 <- df_A9 %>%
  group_by(Community_name) %>%
  summarise(
    DE_in_module = sum(isDE, na.rm = TRUE),
    NonDE_in_module = sum(!isDE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    DE_outside = length(all_DE_A9) - DE_in_module,
    NonDE_outside = length(all_genes_A9) - length(all_DE_A9) - NonDE_in_module
  )

set.seed(123)
fisher_res_A9 <- fisher_tables_A9 %>%
  mutate(
    pval = pmap_dbl(
      list(DE_in_module, NonDE_in_module, DE_outside, NonDE_outside),
      ~ {
        mat <- matrix(c(..1, ..2, ..3, ..4), nrow = 2)
        fisher.test(mat, simulate.p.value = TRUE, B = 1e5)$p.value
      }
    )
  ) %>%
  mutate(padj = p.adjust(pval, method = "BH"))

sig_A9 <- fisher_res_A9 %>% filter(padj < 0.05)
cat("✅ Number of significant modules in A9 (FDR < 0.05):", nrow(sig_A9), "\n")

# ----------------------------------------
# Save Fisher results
# ----------------------------------------
vroom_write(fisher_res_A7, "~/CESC_Network/7_Modularity/7_4_GSEA_module/7_4_2_Fisher_A7.tsv")
vroom_write(fisher_res_A9, "~/CESC_Network/7_Modularity/7_4_GSEA_module/7_4_2_Fisher_A9.tsv")

