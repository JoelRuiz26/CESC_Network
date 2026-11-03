################################################################################
# figure_2_combined.R
#
# CESC/TCGA + Uganda/HTMCP-CC (independent validation cohort) combined view of
# Figure 2. Does NOT touch the original CESC-only Figure 2
# (05_topology/04_figure_2_topology_panel.R -> figures/Figure_2_topology.png,
# left untouched). This script only ADDS a combined version, grouped by
# ANALYSIS TYPE rather than by cohort: each row is one analysis (CCDF, degree
# correlation, edge overlap), with CESC and Uganda side by side under the
# same tag letter.
#
#   A) CCDF                    [CESC | Uganda]
#   B) Gene degree correlation [CESC | Uganda]
#   C) Edge overlap (Jaccard)  [CESC | Uganda]
#
# -> figures/validation_uganda/Figure_2_topology_uganda.png / .pdf
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

fig_dir <- "~/CESC_Network/figures/validation_uganda"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

ttl  <- 25; axt  <- 19; axx  <- 18; lgt  <- 17; lgtt <- 18

base_thm <- theme_bw(base_family = "Times New Roman") +
  theme(plot.title       = element_text(size = ttl, face = "bold", hjust = 0.5),
        axis.title       = element_text(size = axt, face = "bold"),
        axis.text        = element_text(size = axx, color = "black"),
        panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.9),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        legend.title     = element_text(size = lgtt, face = "bold"),
        legend.text      = element_text(size = lgt))

row_title <- function(txt) ggplot() + theme_void() +
  labs(title = txt) +
  theme(plot.title = element_text(family = "Times New Roman", size = ttl + 3, face = "bold",
                                   hjust = 0.5, margin = margin(b = 4)))

build_panels <- function(sub_dir, cohort_label, col_A7, col_A9) {
  A7_HPV <- readRDS(file.path(sub_dir, "5_2_4_Subnetwork_A7_max.rds"))
  A9_HPV <- readRDS(file.path(sub_dir, "5_2_4_Subnetwork_A9_max.rds"))

  deg_A7 <- as.integer(igraph::degree(A7_HPV, mode = "all"))
  deg_A9 <- as.integer(igraph::degree(A9_HPV, mode = "all"))
  ks_result <- suppressWarnings(ks.test(deg_A7, deg_A9))

  build_ccdf_df <- function(graph, label) {
    deg <- as.integer(igraph::degree(graph, mode = "all")); deg <- deg[deg > 0]
    tab <- table(deg); k <- as.integer(names(tab)); ord <- order(k)
    k <- k[ord]; counts <- as.numeric(tab[ord])
    ccdf <- rev(cumsum(rev(counts))) / sum(counts)
    tibble(Degree = k, CCDF = ccdf, Network = factor(label, levels = c("HPV_A7", "HPV_A9")))
  }
  ccdf_df <- bind_rows(build_ccdf_df(A7_HPV, "HPV_A7"), build_ccdf_df(A9_HPV, "HPV_A9"))
  ann_box <- sprintf("KS test: D = %.3f, p < 2.2×10⁻¹⁶", unname(ks_result$statistic))

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
    scale_color_manual(values = c(HPV_A7 = col_A7, HPV_A9 = col_A9), labels = c("HPV-A7", "HPV-A9")) +
    scale_linetype_manual(values = c(HPV_A7 = "solid", HPV_A9 = "dashed"), labels = c("HPV-A7", "HPV-A9")) +
    annotate("label", x = min(ccdf_df$Degree) * 1.25, y = min(ccdf_df$CCDF) * 1.8,
             label = ann_box, hjust = 0, vjust = 0, family = "Times New Roman",
             fontface = "bold", size = 5.6, label.padding = unit(0.35, "lines"),
             fill = "white", color = "grey10") +
    labs(title = cohort_label, x = "Degree (k)", y = expression(P(K >= k)),
         color = "Network", linetype = "Network") +
    base_thm +
    theme(legend.position = c(0.04, 0.96), legend.justification = c("left", "top"),
          legend.background = element_rect(fill = alpha("white", 0.9), color = "grey30"),
          legend.key = element_blank(), legend.key.width = unit(1.0, "cm"))

  common_nodes <- intersect(V(A7_HPV)$name, V(A9_HPV)$name)
  degree_df <- data.frame(
    Gene = common_nodes,
    deg_A7 = as.numeric(igraph::degree(A7_HPV, v = common_nodes)),
    deg_A9 = as.numeric(igraph::degree(A9_HPV, v = common_nodes))) |>
    filter(deg_A7 > 0, deg_A9 > 0) |>
    mutate(color_group = factor(ifelse(deg_A7 - deg_A9 >= 0, "HPV-A7", "HPV-A9"), levels = c("HPV-A7", "HPV-A9")))
  cor_test <- suppressWarnings(cor.test(degree_df$deg_A7, degree_df$deg_A9, method = "spearman"))
  rho_val  <- round(unname(cor_test$estimate), 3)
  ann_label <- sprintf("Spearman ρ = %.3f\np < 2.2 × 10⁻¹⁶", rho_val)
  deg_lim <- range(c(degree_df$deg_A7, degree_df$deg_A9))

  pB <- ggplot(degree_df, aes(deg_A7, deg_A9, color = color_group)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey50") +
    geom_point(alpha = 0.30, size = 1.9, key_glyph = draw_key_point) +
    annotate("label", x = deg_lim[1], y = deg_lim[2], hjust = 0, vjust = 1, label = ann_label,
             fontface = "bold", size = 5.6, fill = alpha("white", 0.85), color = "grey15",
             label.padding = unit(0.3, "lines")) +
    scale_x_log10(limits = deg_lim, breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)), expand = expansion(mult = c(0.05, 0.08))) +
    scale_y_log10(limits = deg_lim, breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)), expand = expansion(mult = c(0.05, 0.08))) +
    coord_fixed(ratio = 1) +
    annotation_logticks(sides = "bl", color = "grey50", linewidth = 0.35) +
    scale_color_manual(values = c("HPV-A7" = col_A7, "HPV-A9" = col_A9),
                       guide = guide_legend(override.aes = list(alpha = 0.85, size = 4))) +
    labs(title = cohort_label,
         x = "Degree in HPV-A7  (log scale)", y = "Degree in HPV-A9  (log scale)", color = NULL) +
    base_thm +
    theme(legend.position = c(0.97, 0.05), legend.justification = c("right", "bottom"),
          legend.background = element_rect(fill = alpha("white", 0.88), color = "grey60", linewidth = 0.4),
          legend.key = element_blank())

  edges_A7 <- apply(as_edgelist(A7_HPV), 1, function(e) paste(sort(e), collapse = "-"))
  edges_A9 <- apply(as_edgelist(A9_HPV), 1, function(e) paste(sort(e), collapse = "-"))
  list_edges <- list("HPV-A7" = edges_A7, "HPV-A9" = edges_A9)

  pC <- ggvenn(list_edges, fill_color = c(col_A7, col_A9), fill_alpha = 0.5,
               stroke_size = 0.7, set_name_size = 6.5, text_size = 5.8,
               show_percentage = TRUE, label_sep = "\n", digits = 1) +
    labs(title = cohort_label) +
    coord_fixed(clip = "off") +
    theme(text = element_text(family = "Times New Roman"),
          plot.title = element_text(hjust = 0.5, face = "bold", size = ttl, family = "Times New Roman"))

  list(pA = pA, pB = pB, pC = pC)
}

cesc   <- build_panels("~/CESC_Network/05_topology/results", "CESC", "#003366", "#CC0000")
uganda <- build_panels("~/CESC_Network/05_topology/validation_uganda/results", "HTMCP-CC", "#D55E00", "#009E73")

## ── grouped by analysis type: each row is CESC | Uganda, one shared tag letter.
## the analysis-type heading (CCDF / Gene degree correlation / Edge overlap) is
## written ONCE per row, above both panels; each panel's own title is just its
## cohort name (CESC / HTMCP-CC — kept homogeneous, neither gets a "TCGA-"
## prefix since HTMCP-CC is not a TCGA project).
## row B (scatter) given more height than A; row C (Venn) given less, so the
## Venn circles render smaller while their internal text (set above) reads bigger
rowA <- wrap_elements(full = (row_title("CCDF") / (cesc$pA | uganda$pA) + plot_layout(heights = c(0.09, 1))))
rowB <- wrap_elements(full = (row_title("Gene degree correlation") / (cesc$pB | uganda$pB) + plot_layout(heights = c(0.09, 1))))
rowC <- wrap_elements(full = (row_title("Edge overlap (Jaccard index)") / (cesc$pC | uganda$pC) + plot_layout(heights = c(0.09, 1))))

p_all <- rowA / rowB / rowC +
  plot_layout(heights = c(1, 1.3, 0.9)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 28, face = "bold", family = "Times New Roman"))

ggsave(file.path(fig_dir, "Figure_2_topology_uganda.png"), p_all, width = 14, height = 16, dpi = 400, bg = "white")
ggsave(file.path(fig_dir, "Figure_2_topology_uganda.pdf"), p_all, width = 14, height = 16, device = cairo_pdf, bg = "white")
message("Saved: Figure_2_topology_uganda.png / .pdf  (A=CCDF, B=degree corr, C=edge overlap; each CESC | Uganda)")
