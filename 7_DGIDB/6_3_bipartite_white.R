# ====================================================
# 6_3_bipartite_white.R
# Bipartite Drug–Gene Network  ·  Dark Arc Layout  ·  White background
# Input : 6_1_edges_bipartite.tsv  (built by 6_1_DGIBD_metrics.R)
# Output: 6_3_bipartite_white.{pdf,png}
# ====================================================

setwd("~/CESC_Network/7_DGIDB/")

library(tidyverse)
library(ggplot2)
library(ggrepel)
library(scales)

set.seed(42)

# ── 1. Load pre-computed edge table ───────────────────────────────────────────
edges_all <- read_tsv("6_1_edges_bipartite.tsv", show_col_types = FALSE) %>%
  mutate(interaction_group = factor(interaction_group,
                                    levels = c("inhibits", "activates", "regulates")))

drug_stats <- edges_all %>%
  group_by(drug_name) %>%
  summarise(n_genes = n_distinct(gene_name), .groups = "drop") %>%
  arrange(desc(n_genes))

gene_stats <- edges_all %>%
  group_by(gene_name) %>%
  summarise(n_drugs = n_distinct(drug_name), .groups = "drop") %>%
  arrange(desc(n_drugs))

n_drugs    <- nrow(drug_stats)
n_genes    <- nrow(gene_stats)
total_ints <- nrow(edges_all)

cat(sprintf("Drugs: %d  |  Genes (A7+A9 network with drug links): %d  |  Interactions: %s\n",
            n_drugs, n_genes, format(total_ints, big.mark = ",")))

# ── 2. Layout parameters ──────────────────────────────────────────────────────
R_ARC <- 6.5
X_GAP <- 1.6

# ── 3. Node positions on facing semicircles ───────────────────────────────────
drug_nodes <- drug_stats %>%
  mutate(
    rank    = row_number(),
    angle   = seq(-pi / 2, pi / 2, length.out = n_drugs),
    x       = -(X_GAP + R_ARC * cos(angle)),
    y       = R_ARC * sin(angle),
    sz_core = rescale(log1p(n_genes), to = c(0.3, 3.5)),
    sz_glow = sz_core * 5
  )

gene_nodes <- gene_stats %>%
  mutate(
    rank    = row_number(),
    angle   = seq(pi / 2, -pi / 2, length.out = n_genes),
    x       = X_GAP + R_ARC * cos(angle),
    y       = R_ARC * sin(angle),
    sz_core = rescale(log1p(n_drugs), to = c(0.3, 3.5)),
    sz_glow = sz_core * 5
  )

# ── 4. Top nodes for label columns ────────────────────────────────────────────
LABEL_X_DRUG <- -(X_GAP + R_ARC + 3.2)
LABEL_X_GENE <-  (X_GAP + R_ARC + 3.2)

top_drug_nodes <- drug_nodes %>%
  filter(rank <= 10) %>%
  mutate(
    label   = paste0(drug_name, "  (", n_genes, ")"),
    nudge_x = LABEL_X_DRUG - x
  )

top_gene_nodes <- gene_nodes %>%
  filter(rank <= 10) %>%
  mutate(
    label   = paste0(gene_name, "  (", n_drugs, ")"),
    nudge_x = LABEL_X_GENE - x
  )

# ── 5. Edge coordinates ───────────────────────────────────────────────────────
edges_coord <- edges_all %>%
  left_join(drug_nodes %>% select(drug_name, x_drug = x, y_drug = y),
            by = "drug_name") %>%
  left_join(gene_nodes %>% select(gene_name, x_gene = x, y_gene = y),
            by = "gene_name")

# ── 6. Arc guide paths ────────────────────────────────────────────────────────
arc_angles    <- seq(-pi / 2, pi / 2, length.out = 400)
arc_drug_path <- data.frame(
  x = -(X_GAP + R_ARC * cos(arc_angles)),
  y = R_ARC * sin(arc_angles)
)
arc_gene_path <- data.frame(
  x = X_GAP + R_ARC * cos(arc_angles),
  y = R_ARC * sin(arc_angles)
)

# ── 7. Colors (adapted for white background) ──────────────────────────────────
BG_COL   <- "#FFFFFF"
DRUG_COL <- "#D4770A"   # dark amber-orange  (drugs)
GENE_COL <- "#0B7A8F"   # deep teal          (genes)
TEXT_COL <- "#1A1A1A"   # near-black for labels and title

group_colors <- c(
  inhibits  = "#E8735A",   # soft coral-red
  activates = "#52C78A",   # soft mint-green
  regulates = "#5BAED6"    # soft sky-blue
)

# ── 8. Build plot ──────────────────────────────────────────────────────────────
p <- ggplot() +

  # Subtle arc guide lines
  geom_path(data = arc_drug_path, aes(x = x, y = y),
            colour = DRUG_COL, alpha = 0.30, linewidth = 0.55) +
  geom_path(data = arc_gene_path, aes(x = x, y = y),
            colour = GENE_COL, alpha = 0.30, linewidth = 0.55) +

  # All edges — overlapping transparency creates density
  geom_segment(data = edges_coord,
               aes(x = x_drug, y = y_drug,
                   xend = x_gene, yend = y_gene,
                   colour = interaction_group),
               alpha = 0.025, linewidth = 0.16) +

  # Drug nodes — outer glow
  geom_point(data = drug_nodes, aes(x = x, y = y, size = sz_glow),
             colour = DRUG_COL, alpha = 0.06, shape = 16) +
  # Drug nodes — mid halo
  geom_point(data = drug_nodes, aes(x = x, y = y, size = sz_core * 2.8),
             colour = DRUG_COL, alpha = 0.14, shape = 16) +
  # Drug nodes — core
  geom_point(data = drug_nodes, aes(x = x, y = y, size = sz_core),
             colour = DRUG_COL, alpha = 0.90, shape = 16) +

  # Gene nodes — outer glow
  geom_point(data = gene_nodes, aes(x = x, y = y, size = sz_glow),
             colour = GENE_COL, alpha = 0.06, shape = 16) +
  # Gene nodes — mid halo
  geom_point(data = gene_nodes, aes(x = x, y = y, size = sz_core * 2.8),
             colour = GENE_COL, alpha = 0.14, shape = 16) +
  # Gene nodes — core
  geom_point(data = gene_nodes, aes(x = x, y = y, size = sz_core),
             colour = GENE_COL, alpha = 0.90, shape = 16) +

  # Drug labels
  geom_text_repel(
    data      = top_drug_nodes,
    aes(x = x, y = y, label = label),
    nudge_x   = top_drug_nodes$nudge_x,
    colour    = TEXT_COL,
    size      = 4.2, fontface = "bold",
    direction = "y", hjust = 1,
    segment.colour = DRUG_COL, segment.alpha = 0.60, segment.size = 0.30,
    box.padding = 0.12, point.padding = 0.08,
    max.overlaps = 40, seed = 42, min.segment.length = 0.05
  ) +

  # Gene labels
  geom_text_repel(
    data      = top_gene_nodes,
    aes(x = x, y = y, label = label),
    nudge_x   = top_gene_nodes$nudge_x,
    colour    = TEXT_COL,
    size      = 4.2, fontface = "bold",
    direction = "y", hjust = 0,
    segment.colour = GENE_COL, segment.alpha = 0.60, segment.size = 0.30,
    box.padding = 0.12, point.padding = 0.08,
    max.overlaps = 35, seed = 123, min.segment.length = 0.05
  ) +

  # Rotated side annotations
  annotate("text", x = LABEL_X_DRUG - 0.7, y = 0,
           label = "FDA DRUGS", colour = DRUG_COL,
           size = 8.5, fontface = "bold", angle = 90, hjust = 0.5) +
  annotate("text", x = LABEL_X_GENE + 0.7, y = 0,
           label = "NETWORK GENES", colour = GENE_COL,
           size = 8.5, fontface = "bold", angle = -90, hjust = 0.5) +

  # Title block
  annotate("text", x = 0, y = R_ARC + 2.2,
           label = "FDA-Approved Drug · Gene Interaction Network",
           colour = TEXT_COL, size = 9.0, fontface = "bold", hjust = 0.5) +
  annotate("text", x = 0, y = R_ARC + 1.5,
           label = sprintf("%s drugs  ·  %s genes  ·  %s interactions  (DGIdb)",
                           format(n_drugs,    big.mark = ","),
                           format(n_genes,    big.mark = ","),
                           format(total_ints, big.mark = ",")),
           colour = "#555555", size = 7.5, fontface = "bold", hjust = 0.5) +

  # Scales
  scale_colour_manual(
    name   = "Interaction type",
    values = group_colors,
    guide  = guide_legend(
      direction      = "horizontal",
      title.position = "top",
      override.aes   = list(linewidth = 2.5, alpha = 0.85)
    )
  ) +
  scale_size_identity(guide  = "none") +
  scale_alpha_identity(guide = "none") +

  coord_fixed(
    xlim = c(LABEL_X_DRUG - 2.0, LABEL_X_GENE + 2.0),
    ylim = c(-R_ARC - 2.5, R_ARC + 3.2)
  ) +

  theme_void(base_size = 14) +
  theme(
    plot.background   = element_rect(fill = BG_COL, colour = NA),
    panel.background  = element_rect(fill = BG_COL, colour = NA),
    legend.position   = c(0.5, 0.03),
    legend.direction  = "horizontal",
    legend.background = element_rect(fill     = "#F5F5F5",
                                     colour   = "#CCCCCC",
                                     linewidth = 0.5),
    legend.margin     = margin(6, 16, 6, 16),
    legend.title      = element_text(colour = TEXT_COL, face = "bold",
                                     size = 17, hjust = 0.5,
                                     margin = margin(b = 5)),
    legend.text       = element_text(colour = "#333333", size = 16),
    legend.key.width  = unit(2.0, "cm"),
    legend.key.height = unit(0.3, "cm"),
    plot.margin       = margin(10, 15, 40, 15)
  )

# ── 9. Save ───────────────────────────────────────────────────────────────────
ggsave("6_3_bipartite_white.pdf", p,
       width = 20, height = 13, bg = BG_COL)
ggsave("6_3_bipartite_white.png", p,
       width = 20, height = 13, dpi = 300, bg = BG_COL)

cat(sprintf(
  "\n✓ Saved: 6_3_bipartite_white.pdf / .png  (%d drugs · %d genes · %s interactions)\n",
  n_drugs, n_genes, format(total_ints, big.mark = ",")
))
