### ============================================================
### 6_6_Hub_dependency_plots.R
###
### Four alternative views of hub-gene dependency (pan-essential
### removed), by hub-sharing category (Shared / A7-specific /
### A9-specific). Uses per-line MEDIAN dependency.
###
###   6_6_1_scatter_A7_vs_A9.png   — median dep A7 vs A9, by category
###   6_6_2_MA_plot.png            — overall dependency vs clade bias
###   6_6_3_waterfall.png          — genes ranked by median dependency
###   6_6_4_scatter_marginal.png   — scatter + marginal densities
### ============================================================
suppressPackageStartupMessages({
  library(depmap)
  library(igraph)
  library(tidyverse)
  library(vroom)
  library(patchwork)
  library(ggrepel)
})

dep_dir <- "~/CESC_Network/6_DepMap_cacu"
mod_dir <- "~/CESC_Network/7_Modularity/5_1_Modules"
HUB_QUANTILE     <- 0.90
ESSENTIAL_CUTOFF <- -0.5

col_A7 <- "#003366"; col_A9 <- "#CC2200"; col_shared <- "#6A3D9A"
cat_colors <- c("Shared" = col_shared, "A7-specific" = col_A7, "A9-specific" = col_A9)
cat_levels <- c("Shared", "A7-specific", "A9-specific")

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
lines  <- vroom(file.path(dep_dir, "6_1_crispr_lines.tsv"), show_col_types = FALSE)
ids_A7 <- lines$depmap_id[lines$hpv_group == "A7"]
ids_A9 <- lines$depmap_id[lines$hpv_group == "A9"]

per_gene <- depmap_crispr() %>%
  filter(depmap_id %in% lines$depmap_id) %>%
  transmute(gene_u = toupper(gene_name), depmap_id, dep = dependency) %>%
  group_by(gene_u) %>%
  summarise(med_A7  = median(dep[depmap_id %in% ids_A7],  na.rm = TRUE),
            med_A9  = median(dep[depmap_id %in% ids_A9],  na.rm = TRUE),
            med_all = median(dep, na.rm = TRUE), .groups = "drop")

pan <- vroom(file.path(dep_dir, "6_5_common_essential_genes.tsv"),
             delim = "\t", show_col_types = FALSE)$gene_u

dat <- hub %>%
  left_join(per_gene, by = "gene_u") %>%
  filter(!(gene_u %in% pan), !is.na(med_A7), !is.na(med_A9)) %>%
  mutate(clade_bias = med_A9 - med_A7)   # <0 more A9-essential, >0 more A7-essential
vroom_write(dat, file.path(dep_dir, "6_6_hub_dependency_medians.tsv"), delim = "\t")
message("Hub genes (pan-essential removed): ", nrow(dat))
print(as.data.frame(count(dat, category)))

## ── Theme ─────────────────────────────────────────────────────────────────────
base_theme <- theme_bw(base_family = "Times New Roman") +
  theme(axis.text   = element_text(size = 17, face = "bold"),
        axis.title  = element_text(size = 21, face = "bold"),
        plot.title  = element_text(size = 24, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 15, face = "italic", hjust = 0.5, color = "grey30"),
        legend.text = element_text(size = 15, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 0.8))

## ── Plot 1: scatter med A7 vs A9 ──────────────────────────────────────────────
# Label logic (by quadrant):
#   both essential (med_A7<=-0.5 & med_A9<=-0.5): only top 10 by combined dep
#   A7-only quadrant (med_A7<=-0.5 & med_A9>-0.5): ALL A7-specific genes
#   A9-only quadrant (med_A9<=-0.5 & med_A7>-0.5): ALL A9-specific genes
lab_both <- dat %>% filter(med_A7 <= ESSENTIAL_CUTOFF, med_A9 <= ESSENTIAL_CUTOFF) %>%
  mutate(s = med_A7 + med_A9) %>% slice_min(s, n = 10)
lab_A7 <- dat %>% filter(category == "A7-specific",
                         med_A7 <= ESSENTIAL_CUTOFF, med_A9 > ESSENTIAL_CUTOFF)
lab_A9 <- dat %>% filter(category == "A9-specific",
                         med_A9 <= ESSENTIAL_CUTOFF, med_A7 > ESSENTIAL_CUTOFF)
lab1 <- bind_rows(lab_both, lab_A7, lab_A9) %>% distinct(gene, .keep_all = TRUE)

# axis limits: only the actual data range (square, small margin)
rng <- range(c(dat$med_A7, dat$med_A9), na.rm = TRUE)
lim <- c(rng[1] - 0.08, rng[2] + 0.08)

p1 <- ggplot(dat, aes(med_A7, med_A9, color = category)) +
  annotate("rect", xmin = -Inf, xmax = ESSENTIAL_CUTOFF, ymin = -Inf, ymax = ESSENTIAL_CUTOFF,
           fill = "#2E8B57", alpha = 0.06) +
  geom_hline(yintercept = ESSENTIAL_CUTOFF, linetype = "dashed", linewidth = 0.7) +
  geom_vline(xintercept = ESSENTIAL_CUTOFF, linetype = "dashed", linewidth = 0.7) +
  geom_abline(slope = 1, intercept = 0, color = "grey65", linewidth = 0.5) +
  geom_point(size = 2.2, alpha = 0.35) +
  geom_text_repel(data = lab1, aes(label = gene), size = 4, fontface = "bold",
                  family = "Times New Roman", max.overlaps = 40, min.segment.length = 0,
                  box.padding = 0.45, seed = 1, show.legend = FALSE,
                  bg.color = "white", bg.r = 0.2) +
  scale_color_manual(values = cat_colors, name = "Hub category") +
  scale_x_continuous(limits = lim) + scale_y_continuous(limits = lim) +
  labs(title = "Hub gene dependency: A7 vs A9 cervical lines",
       subtitle = "Median CRISPR dependency, pan-essential removed",
       x = "Median dependency — A7 lines", y = "Median dependency — A9 lines") +
  base_theme + theme(legend.position = c(0.86, 0.16),
                     legend.background = element_rect(fill = alpha("white", 0.85), color = "grey70"))
ggsave(file.path(dep_dir, "6_6_1_scatter_A7_vs_A9.png"), p1,
       width = 10, height = 11, dpi = 300, bg = "white")
message("[OK] Plot 1 saved.")

## ── Plot 2: MA-plot (overall dependency vs clade bias) ────────────────────────
lab2 <- bind_rows(slice_max(dat, clade_bias, n = 8), slice_min(dat, clade_bias, n = 8))
p2 <- ggplot(dat, aes(med_all, clade_bias, color = category)) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.7) +
  geom_vline(xintercept = ESSENTIAL_CUTOFF, linetype = "dashed", linewidth = 0.7) +
  geom_point(size = 2.4, alpha = 0.7) +
  geom_text_repel(data = lab2, aes(label = gene), size = 3.8, fontface = "bold",
                  family = "Times New Roman", max.overlaps = 30, min.segment.length = 0,
                  show.legend = FALSE, seed = 1) +
  annotate("text", x = 0.35, y = max(dat$clade_bias) * 0.9, label = "more A7-essential",
           hjust = 1, size = 4.4, fontface = "bold.italic", family = "Times New Roman", color = col_A7) +
  annotate("text", x = 0.35, y = min(dat$clade_bias) * 0.9, label = "more A9-essential",
           hjust = 1, size = 4.4, fontface = "bold.italic", family = "Times New Roman", color = col_A9) +
  scale_color_manual(values = cat_colors, name = "Hub category") +
  labs(title = "Hub gene essentiality vs clade bias",
       subtitle = "x: overall dependency  |  y: A9 − A7 difference",
       x = "Median dependency (all cervical lines)",
       y = "Clade bias  (A9 − A7)") +
  base_theme + theme(legend.position = c(0.14, 0.83),
                     legend.background = element_rect(fill = alpha("white", 0.85), color = "grey70"))
ggsave(file.path(dep_dir, "6_6_2_MA_plot.png"), p2,
       width = 11, height = 9, dpi = 300, bg = "white")
message("[OK] Plot 2 saved.")

## ── Plot 3: waterfall (ranked by median dependency) ───────────────────────────
wf <- dat %>% arrange(med_all) %>% mutate(rank = row_number())
p3 <- ggplot(wf, aes(rank, med_all, color = category)) +
  geom_hline(yintercept = ESSENTIAL_CUTOFF, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_point(size = 1.6, alpha = 0.8) +
  annotate("text", x = nrow(wf) * 0.5, y = ESSENTIAL_CUTOFF + 0.08,
           label = "essential threshold (−0.5)", size = 4.4, fontface = "bold.italic",
           family = "Times New Roman", color = "grey30") +
  scale_color_manual(values = cat_colors, name = "Hub category") +
  labs(title = "Hub genes ranked by dependency",
       subtitle = "Each gene ordered by its median dependency (pan-essential removed)",
       x = "Gene rank (most essential → least)", y = "Median dependency") +
  base_theme + theme(legend.position = c(0.14, 0.83),
                     legend.background = element_rect(fill = alpha("white", 0.85), color = "grey70")) +
  guides(color = guide_legend(override.aes = list(size = 4)))
ggsave(file.path(dep_dir, "6_6_3_waterfall.png"), p3,
       width = 12, height = 8, dpi = 300, bg = "white")
message("[OK] Plot 3 saved.")

## ── Plot 4: scatter + marginal densities ──────────────────────────────────────
main4 <- ggplot(dat, aes(med_A7, med_A9, color = category)) +
  geom_hline(yintercept = ESSENTIAL_CUTOFF, linetype = "dashed", linewidth = 0.6) +
  geom_vline(xintercept = ESSENTIAL_CUTOFF, linetype = "dashed", linewidth = 0.6) +
  geom_abline(slope = 1, intercept = 0, color = "grey65", linewidth = 0.5) +
  geom_point(size = 2.1, alpha = 0.65) +
  scale_color_manual(values = cat_colors, name = "Hub category") +
  scale_x_continuous(limits = c(-2.8, 0.5)) + scale_y_continuous(limits = c(-2.8, 0.5)) +
  labs(x = "Median dependency — A7 lines", y = "Median dependency — A9 lines") +
  base_theme + theme(legend.position = c(0.85, 0.16),
                     legend.background = element_rect(fill = alpha("white", 0.85), color = "grey70"))

top_d <- ggplot(dat, aes(med_A7, color = category, fill = category)) +
  geom_density(alpha = 0.25, linewidth = 1) +
  scale_color_manual(values = cat_colors, guide = "none") +
  scale_fill_manual(values = cat_colors, guide = "none") +
  scale_x_continuous(limits = c(-2.8, 0.5)) +
  theme_void(base_family = "Times New Roman")

right_d <- ggplot(dat, aes(med_A9, color = category, fill = category)) +
  geom_density(alpha = 0.25, linewidth = 1) +
  scale_color_manual(values = cat_colors, guide = "none") +
  scale_fill_manual(values = cat_colors, guide = "none") +
  scale_x_continuous(limits = c(-2.8, 0.5)) +
  coord_flip() + theme_void(base_family = "Times New Roman")

p4 <- (top_d + plot_spacer() + main4 + right_d +
       plot_layout(ncol = 2, widths = c(4, 1), heights = c(1, 4))) +
  plot_annotation(
    title = "Hub gene dependency with marginal distributions",
    subtitle = "Median CRISPR dependency by hub category, pan-essential removed",
    theme = theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5,
                    family = "Times New Roman"),
                  plot.subtitle = element_text(size = 15, face = "italic", hjust = 0.5,
                    color = "grey30", family = "Times New Roman")))
ggsave(file.path(dep_dir, "6_6_4_scatter_marginal.png"), p4,
       width = 11, height = 10, dpi = 300, bg = "white")
message("✅ All 4 plots saved to: ", dep_dir)
