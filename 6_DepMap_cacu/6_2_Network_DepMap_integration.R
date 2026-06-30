### ============================================================
### 6_2_Network_DepMap_integration.R
### Integrate Louvain network degree with DepMap dependency scores
### per gene and per module, for HPV-A7 and HPV-A9 networks.
###
### Outputs:
###   6_2_A_gene_degree_vs_dep.png      — scatter: degree vs dep (gene level)
###   6_2_B_module_degree_vs_dep.png    — scatter: mean degree vs mean dep (module)
###   6_2_C_module_dep_boxplot.png      — boxplot: dep distribution per module
###   6_2_gene_dep_degree_table.tsv     — per-gene joined table (A7 + A9)
###   6_2_module_summary_table.tsv      — per-module summary (A7 + A9)
### ============================================================
suppressPackageStartupMessages({
  library(igraph)
  library(tidyverse)
  library(vroom)
  library(patchwork)
  library(scales)
  library(ggrepel)
})

## ── Paths ─────────────────────────────────────────────────────────────────────
mod_dir <- "~/CESC_Network/7_Modularity/5_1_Modules"
dep_dir <- "~/CESC_Network/6_DepMap_cacu"

## ── Load networks and DepMap scores ───────────────────────────────────────────
g7 <- readRDS(file.path(mod_dir, "5_1_1_Louvain_A7_graph_max.rds"))
g9 <- readRDS(file.path(mod_dir, "5_1_1_Louvain_A9_graph_max.rds"))

dep_A7 <- vroom(file.path(dep_dir, "6_1_cervix_dep_mean_A7.tsv"), show_col_types = FALSE)
dep_A9 <- vroom(file.path(dep_dir, "6_1_cervix_dep_mean_A9.tsv"), show_col_types = FALSE)

## ── Per-gene table ─────────────────────────────────────────────────────────────
extract_gene_tbl <- function(g, dep_tbl, label) {
  deg <- igraph::degree(g)
  igraph::as_data_frame(g, what = "vertices") %>%
    transmute(
      gene           = name,
      community_id   = as.integer(Community),
      community_name = Community_name,
      community_size = as.integer(Community_size),
      degree         = as.integer(deg[name])
    ) %>%
    left_join(dep_tbl, by = c("gene" = "gene_u")) %>%
    filter(!is.na(mean_dep)) %>%
    mutate(network = label)
}

gene_A7 <- extract_gene_tbl(g7, dep_A7, "HPV-A7")
gene_A9 <- extract_gene_tbl(g9, dep_A9, "HPV-A9")

message("Genes con dato: A7 = ", nrow(gene_A7), " | A9 = ", nrow(gene_A9))

## ── Per-module summary ─────────────────────────────────────────────────────────
module_summary <- function(gene_tbl) {
  gene_tbl %>%
    group_by(network, community_id, community_name, community_size) %>%
    summarise(
      n_dep        = n(),
      mean_degree  = mean(degree),
      median_degree = median(degree),
      mean_dep     = mean(mean_dep),
      median_dep   = median(mean_dep),
      sd_dep       = sd(mean_dep),
      .groups      = "drop"
    )
}

mod_A7 <- module_summary(gene_A7)
mod_A9 <- module_summary(gene_A9)

## ── Save tables ───────────────────────────────────────────────────────────────
bind_rows(gene_A7, gene_A9) %>%
  vroom_write(file.path(dep_dir, "6_2_gene_dep_degree_table.tsv"), delim = "\t")

bind_rows(mod_A7, mod_A9) %>%
  arrange(network, median_dep) %>%
  vroom_write(file.path(dep_dir, "6_2_module_summary_table.tsv"), delim = "\t")

message("[OK] Tables saved.")

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
    panel.border     = element_rect(color = "black", linewidth = 0.8),
    strip.text       = element_text(size = 18, face = "bold")
  )

spearman_label <- function(x, y) {
  r   <- cor(x, y, method = "spearman")
  pv  <- cor.test(x, y, method = "spearman")$p.value
  pvl <- if (pv < 0.001) "p < 0.001" else sprintf("p = %.3f", pv)
  sprintf("ρ = %.3f\n%s\nn = %s", r, pvl, format(length(x), big.mark = ","))
}

## ── Plot A: Gene-level scatter (degree vs dep) ────────────────────────────────
make_gene_scatter <- function(gene_tbl, color, title) {
  ann <- spearman_label(gene_tbl$degree, gene_tbl$mean_dep)

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
  theme = theme(
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5,
                              family = "Times New Roman", margin = margin(b = 10))
  )
)

ggsave(file.path(dep_dir, "6_2_A_gene_degree_vs_dep.png"),
       pA, width = 22, height = 9, dpi = 300, bg = "white")
message("[OK] Plot A saved.")

## ── Plot B: Module-level scatter (mean degree vs mean dep) ────────────────────
make_module_scatter <- function(mod_tbl, color, title) {
  ann <- spearman_label(mod_tbl$mean_degree, mod_tbl$mean_dep)

  ggplot(mod_tbl, aes(x = mean_degree, y = mean_dep)) +
    geom_smooth(method = "lm", color = "grey40", se = TRUE,
                linewidth = 0.8, alpha = 0.15) +
    geom_point(aes(size = community_size), color = color, alpha = 0.75) +
    geom_text_repel(
      aes(label = community_name),
      size = 4.8, fontface = "bold", family = "Times New Roman",
      max.overlaps = 40, box.padding = 0.45, min.segment.length = 0.2,
      segment.color = "grey50", segment.size = 0.4
    ) +
    annotate("text", x = Inf, y = Inf, label = ann,
             hjust = 1.05, vjust = 1.2, size = 5.5, fontface = "bold",
             family = "Times New Roman") +
    scale_size_continuous(name = "Module size", range = c(3, 13),
                          labels = comma) +
    labs(title = title, x = "Mean degree (k)", y = "Mean dependency score") +
    base_theme +
    theme(legend.position = "right")
}

pB_A7 <- make_module_scatter(mod_A7, "#003366", "HPV-A7")
pB_A9 <- make_module_scatter(mod_A9, "#CC2200", "HPV-A9")

pB <- pB_A7 + pB_A9 +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "Module mean degree vs. mean DepMap dependency score",
    theme = theme(
      plot.title = element_text(size = 30, face = "bold", hjust = 0.5,
                                family = "Times New Roman", margin = margin(b = 10))
    )
  )

ggsave(file.path(dep_dir, "6_2_B_module_degree_vs_dep.png"),
       pB, width = 24, height = 10, dpi = 300, bg = "white")
message("[OK] Plot B saved.")

## ── Plot C: Boxplot per module (dep distribution, sorted by median) ───────────
make_module_boxplot <- function(gene_tbl, mod_tbl, color, title) {
  # Sort modules by median DepMap (most essential first = most negative)
  lvls <- mod_tbl %>% arrange(desc(median_dep)) %>% pull(community_name)

  # Join module mean_degree for fill color
  gene_plot <- gene_tbl %>%
    left_join(mod_tbl %>% select(community_name, mean_degree, median_dep),
              by = "community_name") %>%
    mutate(community_name = factor(community_name, levels = lvls))

  # Annotate n genes per module (right side)
  n_ann <- mod_tbl %>%
    mutate(community_name = factor(community_name, levels = lvls),
           label = paste0("n=", n_dep))

  ggplot(gene_plot, aes(x = community_name, y = mean_dep, fill = mean_degree)) +
    geom_hline(yintercept = -0.5, linetype = "dashed",
               color = "firebrick3", linewidth = 0.7, alpha = 0.8) +
    geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.25, linewidth = 0.45) +
    geom_text(data = n_ann, aes(x = community_name, y = Inf, label = label),
              inherit.aes = FALSE, hjust = -0.05, vjust = 0.5,
              size = 3.8, fontface = "bold", family = "Times New Roman") +
    scale_fill_gradient(low = "#D6EAF8", high = color,
                        name = "Mean\ndegree (k)") +
    scale_y_continuous(labels = number_format(accuracy = 0.1)) +
    coord_flip(clip = "off") +
    labs(title = title, x = NULL, y = "Mean dependency score") +
    base_theme +
    theme(
      axis.text.y      = element_text(size = 14, face = "bold"),
      axis.text.x      = element_text(size = 16, face = "bold"),
      legend.position  = "right",
      plot.margin      = margin(5, 55, 5, 5)
    )
}

pC_A7 <- make_module_boxplot(gene_A7, mod_A7, "#003366", "HPV-A7  (32 modules)")
pC_A9 <- make_module_boxplot(gene_A9, mod_A9, "#CC2200", "HPV-A9  (41 modules)")

pC <- pC_A7 + pC_A9 +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "right")

pC <- pC + plot_annotation(
  title = "DepMap dependency score distribution per Louvain module",
  theme = theme(
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5,
                              family = "Times New Roman", margin = margin(b = 10))
  )
)

ggsave(file.path(dep_dir, "6_2_C_module_dep_boxplot.png"),
       pC, width = 30, height = 16, dpi = 300, bg = "white")
message("[OK] Plot C saved.")

message("✅ All outputs saved to: ", dep_dir)
