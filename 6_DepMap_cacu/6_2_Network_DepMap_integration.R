### ============================================================
### 6_2_Network_DepMap_integration.R
### Join Louvain network degree with DepMap dependency scores.
### Shows weak continuous correlation (motivates categorical approach).
###
### Outputs:
###   6_2_A_gene_degree_vs_dep.png  — scatter: degree vs dep (gene level)
###   6_2_gene_dep_degree_table.tsv — per-gene joined table (A7 + A9)
### ============================================================
suppressPackageStartupMessages({
  library(igraph)
  library(tidyverse)
  library(vroom)
  library(patchwork)
  library(scales)
})

mod_dir <- "~/CESC_Network/7_Modularity/5_1_Modules"
dep_dir <- "~/CESC_Network/6_DepMap_cacu"

g7 <- readRDS(file.path(mod_dir, "5_1_1_Louvain_A7_graph_max.rds"))
g9 <- readRDS(file.path(mod_dir, "5_1_1_Louvain_A9_graph_max.rds"))

dep_A7 <- vroom(file.path(dep_dir, "6_1_cervix_dep_mean_A7.tsv"), show_col_types = FALSE)
dep_A9 <- vroom(file.path(dep_dir, "6_1_cervix_dep_mean_A9.tsv"), show_col_types = FALSE)

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
  vroom_write(file.path(dep_dir, "6_2_gene_dep_degree_table.tsv"), delim = "\t")
message("[OK] Table saved.")

## ── Theme ─────────────────────────────────────────────────────────────────────
base_theme <- theme_bw(base_family = "Times New Roman") +
  theme(
    axis.text        = element_text(size = 18, face = "bold"),
    axis.title       = element_text(size = 22, face = "bold"),
    plot.title       = element_text(size = 26, face = "bold", hjust = 0.5,
                                    margin = margin(b = 8)),
    legend.title     = element_text(size = 17, face = "bold"),
    legend.text      = element_text(size = 15, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", linewidth = 0.8)
  )

## ── Plot A: Gene-level scatter (degree vs dep) ────────────────────────────────
make_gene_scatter <- function(gene_tbl, color, title) {
  r  <- cor(gene_tbl$degree, gene_tbl$mean_dep, method = "spearman")
  pv <- cor.test(gene_tbl$degree, gene_tbl$mean_dep, method = "spearman")$p.value
  ann <- sprintf("ρ = %.3f\np < 0.001\nn = %s", r, format(nrow(gene_tbl), big.mark = ","))

  ggplot(gene_tbl, aes(x = degree, y = mean_dep)) +
    geom_hex(bins = 60, aes(fill = after_stat(log10(count)))) +
    scale_fill_gradient(low = "grey90", high = color,
                        name = "log₁₀(n genes)",
                        labels = function(x) round(10^x)) +
    geom_smooth(method = "lm", color = "black", se = TRUE,
                linewidth = 0.9, alpha = 0.2) +
    annotate("text", x = Inf, y = Inf, label = ann,
             hjust = 1.05, vjust = 1.2, size = 5.5, fontface = "bold",
             family = "Times New Roman") +
    scale_x_continuous(breaks = pretty_breaks(5)) +
    labs(title = title, x = "Degree (k)", y = "Mean dependency score") +
    base_theme +
    theme(legend.position = "right")
}

pA_A7 <- make_gene_scatter(gene_A7, "#003366", "HPV-A7")
pA_A9 <- make_gene_scatter(gene_A9, "#CC2200", "HPV-A9")

pA <- pA_A7 + pA_A9 +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "right")

pA <- pA + plot_annotation(
  title = "Gene degree vs. DepMap dependency score",
  theme = theme(plot.title = element_text(size = 30, face = "bold", hjust = 0.5,
                                          family = "Times New Roman",
                                          margin = margin(b = 10)))
)

ggsave(file.path(dep_dir, "6_2_A_gene_degree_vs_dep.png"),
       pA, width = 22, height = 9, dpi = 300, bg = "white")
message("✅ Plot A saved.")
