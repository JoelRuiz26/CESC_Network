### ============================================================
### 05_figure_panel.R
###
### Final combined figure, CESC A7 vs A9 latent-space comparison.
###   A) Scree plot and the embedding dimension (profile-likelihood elbow)
###   B) Aligned latent space, A7 vs A9 overall (density contours)
###   C) All significant genes (FDR<0.01) as a displacement flow field,
###      top 10 labelled
###   D) Degree vs. latent distance, coloured by significance -- shows why
###      the test had to control for degree (hubs displace more by
###      construction) instead of repeating C's gene labels
###
### -> figures/Figure_latent_embedding.png / .pdf
### ============================================================
suppressPackageStartupMessages({ library(igraph); library(dplyr); library(tidyr); library(readr); library(ggplot2); library(ggrepel); library(patchwork) })

res_dir  <- "~/CESC_Network/11_latent_embedding/results"
topo_dir <- "~/CESC_Network/05_topology/results"
fig_dir  <- "~/CESC_Network/11_latent_embedding/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
FF <- "Times New Roman"
COL_A7 <- "#003366"; COL_A9 <- "#CC0000"; COL_SIG <- "#CC0000"; COL_NS <- "grey75"
P_CUT <- 0.01

X7         <- readRDS(file.path(res_dir, "01_ase_X7.rds"))
X9_aligned <- readRDS(file.path(res_dir, "01_ase_X9_aligned.rds"))
genes      <- read_tsv(file.path(res_dir, "01_shared_gene_order.tsv"), show_col_types = FALSE)$gene
scree      <- read_tsv(file.path(res_dir, "01_scree_singular_values.tsv"), show_col_types = FALSE)
sig        <- read_tsv(file.path(res_dir, "04_gene_significance.tsv"), show_col_types = FALSE)

base <- theme_bw(base_family = FF) +
  theme(axis.text = element_text(size = 10, colour = "grey15"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 13, face = "bold"),
        legend.title = element_blank(), legend.text = element_text(size = 9),
        legend.key.size = unit(0.35, "cm"),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2),
        panel.border = element_rect(colour = "grey30", linewidth = 0.6))

## -- A) Scree / dimension selection --------------------------------------
d <- scree$elbow_d[1]
pA <- ggplot(scree, aes(k, singular_value)) +
  geom_line(colour = "grey40", linewidth = 0.4) +
  geom_point(size = 1.4, colour = "grey30") +
  geom_point(data = scree |> filter(k == d), colour = COL_A9, size = 2.6) +
  geom_vline(xintercept = d, linetype = "dashed", colour = COL_A9, linewidth = 0.4) +
  annotate("text", x = d + 1, y = max(scree$singular_value), label = sprintf("d = %d", d),
           family = FF, fontface = "bold", colour = COL_A9, size = 3.6, hjust = 0) +
  labs(title = "Dimension selection", x = "k", y = "Singular value") + base

## -- B) Overall latent space, A7 vs A9 (density contours) ----------------
pca <- prcomp(rbind(X7, X9_aligned), center = TRUE, scale. = FALSE)
pts <- tibble(gene = rep(genes, 2), network = rep(c("A7", "A9"), each = length(genes)),
              PC1 = pca$x[, 1], PC2 = pca$x[, 2])
pB <- ggplot(pts, aes(PC1, PC2, colour = network)) +
  geom_density_2d(linewidth = 0.55, bins = 7) +
  scale_colour_manual(values = c(A7 = COL_A7, A9 = COL_A9)) +
  labs(title = "Latent space (A7 vs A9)", x = "PC1", y = "PC2") + base +
  theme(legend.position = c(0.88, 0.88))

## -- C) ALL significant genes as a displacement flow field ----------------
sig_genes <- sig |> filter(p_adj < P_CUT) |> pull(gene)
top10 <- sig |> filter(p_adj < P_CUT) |> slice_min(p_adj, n = 10)

seg_all <- pts |> filter(gene %in% sig_genes) |>
  pivot_wider(id_cols = gene, names_from = network, values_from = c(PC1, PC2))
seg_top <- seg_all |> filter(gene %in% top10$gene)

pC <- ggplot() +
  geom_point(data = pts, aes(PC1, PC2), colour = "grey85", size = 0.25, alpha = 0.35) +
  geom_segment(data = seg_all, aes(x = PC1_A7, y = PC2_A7, xend = PC1_A9, yend = PC2_A9),
               colour = COL_SIG, linewidth = 0.15, alpha = 0.1) +
  geom_segment(data = seg_top, aes(x = PC1_A7, y = PC2_A7, xend = PC1_A9, yend = PC2_A9),
               colour = "grey20", linewidth = 0.5, arrow = arrow(length = unit(0.07, "inches"))) +
  geom_point(data = pts |> filter(gene %in% top10$gene), aes(PC1, PC2, colour = network), size = 2.4) +
  geom_text_repel(data = pts |> filter(gene %in% top10$gene, network == "A9"), aes(PC1, PC2, label = gene),
                   family = FF, size = 2.9, fontface = "bold.italic", max.overlaps = Inf,
                   box.padding = 0.8, point.padding = 0.2, force = 8, force_pull = 0.3, max.iter = 20000,
                   min.segment.length = 0, segment.size = 0.25, segment.alpha = 0.6, seed = 1) +
  scale_colour_manual(values = c(A7 = COL_A7, A9 = COL_A9), guide = "none") +
  labs(title = sprintf("Significant genes (FDR<%.2f, n=%d)", P_CUT, length(sig_genes)), x = "PC1", y = "PC2") + base

## -- D) Degree vs. latent distance, coloured by significance ---------------
pD <- ggplot(sig |> arrange(p_adj < P_CUT), aes(min_degree, distance, colour = p_adj < P_CUT)) +
  geom_point(size = 0.7, alpha = 0.55) +
  scale_colour_manual(values = c(`TRUE` = COL_SIG, `FALSE` = COL_NS), labels = c(`TRUE` = "Significant", `FALSE` = "Not significant")) +
  scale_x_log10() +
  labs(title = "Latent distance vs. degree", x = "min(degree A7, A9), log scale", y = "Latent distance") + base +
  theme(legend.position = c(0.78, 0.18))

## -- Assemble --------------------------------------------------------------
fig <- (pA | pB) / (pC | pD) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = FF, size = 14, face = "bold"))

ggsave(file.path(fig_dir, "Figure_latent_embedding.png"), fig, width = 11, height = 10, dpi = 400, bg = "white")
ggsave(file.path(fig_dir, "Figure_latent_embedding.pdf"), fig, width = 11, height = 10, device = cairo_pdf, bg = "white")
message("[OK] Saved: figures/Figure_latent_embedding.png / .pdf")
