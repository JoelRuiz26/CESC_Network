################################################################################
# 5_2_6_Topology_panel_patchwork.R
#
# Assembles the three network-topology panels into a single horizontal,
# paper-sized figure:
#
#   A) Complementary Cumulative Degree Distribution (CCDF)   [left, full height]
#   B) Gene degree correlation                               [top right]
#   C) Edge overlap (Jaccard index)                          [bottom right]
#
# Panels are rebuilt here from the same graphs used by 5_2_2 / 5_2_3_scatter /
# 5_2_Network_topology, with font sizes retuned for a combined figure. Only the
# three panel titles are shown (no technical subtitles).
################################################################################

suppressPackageStartupMessages({
  library(igraph)
  library(ggplot2)
  library(ggrepel)
  library(ggvenn)
  library(scales)
  library(dplyr)
  library(tibble)
  library(grid)
  library(patchwork)
})

sub_dir <- "~/CESC_Network/5_Network_analisis/5_1_Subnetwork"
out_dir <- "~/CESC_Network/5_Network_analisis/5_2_Network_topology"

A7_HPV <- readRDS(file.path(sub_dir, "5_2_4_Subnetwork_A7_max.rds"))
A9_HPV <- readRDS(file.path(sub_dir, "5_2_4_Subnetwork_A9_max.rds"))

col_A7 <- "#003366"; col_A9 <- "#CC0000"

## font sizes for the COMBINED figure (smaller than the standalone exports)
ttl  <- 24   # panel title
axt  <- 18   # axis title
axx  <- 18   # axis text
lgt  <- 18   # legend text
lgtt <- 20   # legend title

base_thm <- theme_bw(base_family = "Times New Roman") +
  theme(plot.title       = element_text(size = ttl, face = "bold", hjust = 0.5),
        axis.title       = element_text(size = axt, face = "bold"),
        axis.text        = element_text(size = axx, color = "black"),
        panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.9),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        legend.title     = element_text(size = lgtt, face = "bold"),
        legend.text      = element_text(size = lgt))

################################################################################
# A) CCDF
################################################################################
deg_A7 <- as.integer(igraph::degree(A7_HPV, mode = "all"))
deg_A9 <- as.integer(igraph::degree(A9_HPV, mode = "all"))
ks_result <- suppressWarnings(ks.test(deg_A7, deg_A9))

build_ccdf_df <- function(graph, label) {
  deg <- as.integer(igraph::degree(graph, mode = "all")); deg <- deg[deg > 0]
  tab <- table(deg); k <- as.integer(names(tab)); ord <- order(k)
  k <- k[ord]; counts <- as.numeric(tab[ord])
  ccdf <- rev(cumsum(rev(counts))) / sum(counts)
  tibble(Degree = k, CCDF = ccdf,
         Network = factor(label, levels = c("HPV_A7", "HPV_A9")))
}
ccdf_df <- bind_rows(build_ccdf_df(A7_HPV, "HPV_A7"),
                     build_ccdf_df(A9_HPV, "HPV_A9"))
ann_box <- sprintf("KS test: D = %.3f, p < 2.2×10⁻¹⁶",
                   unname(ks_result$statistic))

pA <- ggplot(ccdf_df, aes(Degree, CCDF, color = Network, linetype = Network)) +
  geom_step(linewidth = 0.9, alpha = 0.95, direction = "vh") +
  geom_point(size = 1.5, alpha = 0.6) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                expand = expansion(mult = c(0.02, 0.05))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                expand = expansion(mult = c(0.05, 0.02))) +
  annotation_logticks(sides = "bl") +
  scale_color_manual(values = c(HPV_A7 = col_A7, HPV_A9 = col_A9),
                     labels = c("HPV-A7", "HPV-A9")) +
  scale_linetype_manual(values = c(HPV_A7 = "solid", HPV_A9 = "dashed"),
                        labels = c("HPV-A7", "HPV-A9")) +
  annotate("label", x = min(ccdf_df$Degree) * 1.25, y = min(ccdf_df$CCDF) * 1.8,
           label = ann_box, hjust = 0, vjust = 0, family = "Times New Roman",
           fontface = "bold", size = 6.0, label.size = 0.3,
           label.padding = unit(0.4, "lines"), fill = "white", color = "grey10") +
  labs(title = "Complementary Cumulative Degree\nDistribution (CCDF)",
       x = "Degree (k)", y = expression(P(K >= k)),
       color = "Network", linetype = "Network") +
  base_thm +
  theme(legend.position      = c(0.96, 0.96),
        legend.justification = c("right", "top"),
        legend.background    = element_rect(fill = alpha("white", 0.9), color = "grey30"),
        legend.key           = element_blank(),
        legend.key.width     = unit(1.1, "cm"))

################################################################################
# B) Gene degree correlation
################################################################################
common_nodes <- intersect(V(A7_HPV)$name, V(A9_HPV)$name)
degree_df <- data.frame(
  Gene   = common_nodes,
  deg_A7 = as.numeric(igraph::degree(A7_HPV, v = common_nodes)),
  deg_A9 = as.numeric(igraph::degree(A9_HPV, v = common_nodes))) |>
  filter(deg_A7 > 0, deg_A9 > 0) |>
  mutate(color_group = factor(ifelse(deg_A7 - deg_A9 >= 0, "HPV-A7", "HPV-A9"),
                              levels = c("HPV-A7", "HPV-A9")))
cor_test <- suppressWarnings(cor.test(degree_df$deg_A7, degree_df$deg_A9, method = "spearman"))
rho_val  <- round(unname(cor_test$estimate), 3)
ann_label <- sprintf("Spearman ρ = %.3f\n(moderate)\np < 2.2 × 10⁻¹⁶", rho_val)

## both axes are the same quantity (degree), so they share one common range and
## an equal (square) aspect: the identity line is then a true 45 degrees
deg_lim <- range(c(degree_df$deg_A7, degree_df$deg_A9))

pB <- ggplot(degree_df, aes(deg_A7, deg_A9, color = color_group)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.5, color = "grey50") +
  geom_point(alpha = 0.28, size = 1.5, key_glyph = draw_key_point) +
  annotate("label", x = deg_lim[1], y = deg_lim[2],
           hjust = 0, vjust = 1, label = ann_label, fontface = "bold", size = 6.0,
           fill = alpha("white", 0.85), color = "grey15", label.size = 0.35,
           label.padding = unit(0.35, "lines")) +
  scale_x_log10(limits = deg_lim,
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                expand = expansion(mult = c(0.05, 0.08))) +
  scale_y_log10(limits = deg_lim,
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                expand = expansion(mult = c(0.05, 0.08))) +
  coord_fixed(ratio = 1) +
  annotation_logticks(sides = "bl", color = "grey50", linewidth = 0.35) +
  scale_color_manual(values = c("HPV-A7" = col_A7, "HPV-A9" = "#CC2200"),
                     guide = guide_legend(override.aes = list(alpha = 0.85, size = 5))) +
  labs(title = "Gene degree correlation",
       x = "Degree in HPV-A7  (log scale)",
       y = "Degree in HPV-A9  (log scale)", color = NULL) +
  base_thm +
  theme(legend.position      = c(0.97, 0.05),
        legend.justification = c("right", "bottom"),
        legend.background    = element_rect(fill = alpha("white", 0.88),
                                            color = "grey60", linewidth = 0.4),
        legend.key           = element_blank())

################################################################################
# C) Edge overlap (Jaccard index)
################################################################################
edges_A7 <- apply(as_edgelist(A7_HPV), 1, function(e) paste(sort(e), collapse = "-"))
edges_A9 <- apply(as_edgelist(A9_HPV), 1, function(e) paste(sort(e), collapse = "-"))
list_edges <- list("HPV-A7" = edges_A7, "HPV-A9" = edges_A9)

pC <- ggvenn(list_edges,
             fill_color = c(col_A7, "#FF0000"), fill_alpha = 0.5,
             stroke_size = 0.7, set_name_size = 6.5, text_size = 5.2,
             show_percentage = TRUE, label_sep = "\n", digits = 1) +
  labs(title = "Edge overlap (Jaccard index)") +
  coord_fixed(clip = "off") +
  theme(text       = element_text(family = "Times New Roman"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = ttl,
                                  family = "Times New Roman"))

################################################################################
# Assemble: left column A (top) + C Venn (bottom); B scatter on the right
################################################################################
## A (CCDF) is the main panel and takes the wider left column; the Venn (C) sits
## directly beneath it. B is a square degree-vs-degree scatter given the whole
## right column at full height, but the narrower right column keeps it medium
## sized rather than dominating the figure.
## A four-row grid gives finer vertical control. A takes the taller top two rows
## and C the shorter bottom two on the left; B spans rows 1-3 on the right, so it
## sits high (near its title) with an empty fourth row beneath it. C keeps its own
## cell, so nudging B up no longer shrinks the Venn.
design <- c(area(1, 1, 2, 1),   # A: top left  (rows 1-2, taller)
            area(1, 2, 3, 2),   # B: right, rows 1-3 -> sits high near its title
            area(3, 1, 4, 1))   # C: bottom left (Venn, rows 3-4, shorter)

p_all <- pA + pB + pC +
  plot_layout(design = design, widths = c(1.2, 1.12),
              heights = c(1.35, 1.35, 0.85, 0.85)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 28, face = "bold", family = "Times New Roman"))

ggsave(file.path(out_dir, "5_2_6_Topology_panel.png"),
       p_all, width = 18, height = 11, dpi = 400, bg = "white")
ggsave(file.path(out_dir, "5_2_6_Topology_panel.pdf"),
       p_all, width = 18, height = 11, device = cairo_pdf, bg = "white")
message("Saved: 5_2_6_Topology_panel.png / .pdf")
