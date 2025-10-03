# ============================================================
# Bubble plots: MeSH & ChEMBL + 1×2 grid with equal panel width
# - Keep inner grid lines, remove panel frames
# - Bigger "Signature" legend keys and extra spacing from title
# - Same X range and same export size for both panels
# - One centered caption for both plots, one unified legend
# - High-res export (PDF + PNG @ 1200 dpi)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(viridis)   # color palette
  library(cowplot)   # equal panel width + composing caption/legend
})

# -------------------- Paths & data --------------------
out_base <- "~/CESC_Network/6_OCTAD/6_4_EnrichedDrugs"
rds_in   <- file.path(out_base, "6_4_1_enrichment_sig_with_n_present.rds")
enrich   <- readRDS(rds_in)

# Consistent labels
enrich <- enrich %>%
  mutate(signature = dplyr::recode(signature, "A7" = "HPV-A7", "A9" = "HPV-A9"))

# Global filters (adjust if needed)
alpha_cut   <- 0.01
min_present <- 3
score_min   <- 0.33

# ---------- Common X range for both panels ----------
xmax_global <- enrich %>%
  filter(target_type %in% c("mesh", "chembl_targets"),
         signature %in% c("HPV-A7","HPV-A9"),
         padj <= alpha_cut, n_present >= min_present, score >= score_min) %>%
  summarise(xmax = ceiling(max(score, na.rm = TRUE) * 20) / 20) %>% pull(xmax)
x_limits <- c(0.30, xmax_global)
x_breaks <- seq(x_limits[1], x_limits[2], by = 0.05)

# ---------- Common theme ----------
theme_bubbles <- theme_minimal(base_size = 16) +
  theme(
    # keep inner grid (slightly thin)
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.major.y = element_line(color = "grey88", linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    
    # remove facet strips and any panel frame/box
    strip.text.y       = element_blank(),
    strip.background   = element_blank(),
    panel.border       = element_blank(),
    
    # white backgrounds (avoid any dark backdrop)
    panel.background   = element_rect(fill = "white", color = NA),
    plot.background    = element_rect(fill = "white", color = NA),
    
    # legend tweaks: bigger keys and extra spacing after title
    legend.position    = "right",
    legend.key.size    = grid::unit(7, "mm"),
    legend.text        = element_text(size = 13),
    legend.title       = element_text(size = 13, face = "bold", margin = margin(b = 6)),
    legend.spacing.x   = grid::unit(8, "pt"),
    
    plot.margin        = margin(12, 20, 22, 28)
  )

# ---------- Build MeSH plot (no caption here; we'll add one globally) ----------
mesh_df <- enrich %>%
  filter(target_type == "mesh",
         signature %in% c("HPV-A7","HPV-A9"),
         padj <= alpha_cut, n_present >= min_present, score >= score_min) %>%
  arrange(signature, padj, dplyr::desc(score)) %>%
  group_by(signature) %>%
  mutate(target = factor(target, levels = rev(unique(target[order(-score)])))) %>%
  ungroup()

p_mesh <- ggplot(mesh_df, aes(score, target)) +
  geom_point(aes(size = n_present, fill = signature),
             shape = 21, color = "grey20", stroke = 0.25, alpha = 0.95) +
  scale_size_continuous(name = "n   ", range = c(3.2, 10), breaks = c(4, 7, 10),    
                        labels = c("4","7","10")) +
  scale_fill_viridis_d(name = "Signature    ", option = "D", end = 0.85) +
  scale_x_continuous(limits = x_limits, breaks = x_breaks,
                     labels = scales::label_number(accuracy = 0.01),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(expand = expansion(add = 0.8)) +
  labs(title = "MeSH", x = "Enrichment score", y = NULL, caption = NULL) +
  facet_grid(signature ~ ., scales = "free_y", space = "free_y") +
  theme_bubbles +
  guides(fill = guide_legend(override.aes = list(size = 6, alpha = 1, shape = 21, color = "grey20")),
         size = guide_legend(order = 2))

# ---------- Build ChEMBL plot (no caption here) ----------
chembl_df <- enrich %>%
  filter(target_type == "chembl_targets",
         signature %in% c("HPV-A7","HPV-A9"),
         padj <= alpha_cut, n_present >= min_present, score >= score_min) %>%
  arrange(signature, padj, dplyr::desc(score)) %>%
  group_by(signature) %>%
  mutate(target = factor(target, levels = rev(unique(target[order(-score)])))) %>%
  ungroup()

p_chembl <- ggplot(chembl_df, aes(score, target)) +
  geom_point(aes(size = n_present, fill = signature),
             shape = 21, color = "grey20", stroke = 0.25, alpha = 0.95) +
  scale_size_continuous(name = "n", range = c(3.2, 10)) +
  scale_fill_viridis_d(name = "Signature    ", option = "D", end = 0.85) +
  scale_x_continuous(limits = x_limits, breaks = x_breaks,
                     labels = scales::label_number(accuracy = 0.01),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(expand = expansion(add = 0.8)) +
  labs(title = "ChEMBL targets", x = "Enrichment score", y = NULL, caption = NULL) +
  facet_grid(signature ~ ., scales = "free_y", space = "free_y") +
  theme_bubbles +
  guides(
    fill = guide_legend(override.aes = list(size = 6, alpha = 1, shape = 21, color = "grey20")),
    size = guide_legend(order = 2)
  )

# ---------- Save individual figures (same size for both) ----------
w_in <- 15; h_in <- 9
ggsave(file.path(out_base, "6_4_plot_bubble_mesh_fixed.pdf"),
       p_mesh + labs(caption = sprintf("All points are significant (FDR \u2264 %.2f)", alpha_cut)),
       width = w_in, height = h_in, units = "in", dpi = 1200, useDingbats = FALSE)
ggsave(file.path(out_base, "6_4_plot_bubble_mesh_fixed.png"),
       p_mesh + labs(caption = sprintf("All points are significant (FDR \u2264 %.2f)", alpha_cut)),
       width = w_in, height = h_in, units = "in", dpi = 1200)

ggsave(file.path(out_base, "6_4_plot_bubble_chembl.pdf"),
       p_chembl + labs(caption = sprintf("All points are significant (FDR \u2264 %.2f)", alpha_cut)),
       width = w_in, height = h_in, units = "in", dpi = 1200, useDingbats = FALSE)
ggsave(file.path(out_base, "6_4_plot_bubble_chembl.png"),
       p_chembl + labs(caption = sprintf("All points are significant (FDR \u2264 %.2f)", alpha_cut)),
       width = w_in, height = h_in, units = "in", dpi = 1200)

# ---------- Grid with equal panel width + one centered caption ----------
# 1) Build a single legend
legend_combined <- cowplot::get_legend(p_mesh + theme(legend.position = "bottom"))

# 2) Remove legends from each panel
p_mesh_nl   <- p_mesh   + theme(legend.position = "none")
p_chembl_nl <- p_chembl + theme(legend.position = "none")

# 3) Align panels so their *inner panel area* has the same width
aligned <- cowplot::align_plots(p_mesh_nl, p_chembl_nl, align = "h", axis = "tb")
grid_panels <- cowplot::plot_grid(aligned[[1]], aligned[[2]],
                                  ncol = 2, rel_widths = c(1, 1))

# 4) Create a centered caption as its own grob
cap_text <- sprintf("All points are significant (FDR \u2264 %.2f)", alpha_cut)
cap_grob <- cowplot::ggdraw() + cowplot::draw_label(cap_text, x = 0.5, y = 0.5,
                                                    hjust = 0.5, vjust = 0.5, size = 12)

# 5) Stack: panels (top) + caption (middle) + legend (bottom)
p_grid <- cowplot::plot_grid(grid_panels, cap_grob, legend_combined,
                             ncol = 1, rel_heights = c(1, 0.07, 0.12))
plot(p_grid)

# 6) Save grid
ggsave(file.path(out_base, "6_4_plot_bubble_grid_mesh_chembl_equalwidth.pdf"),
       p_grid, width = 18, height = 9.5, units = "in", dpi = 1200, useDingbats = FALSE)
ggsave(file.path(out_base, "6_4_plot_bubble_grid_mesh_chembl_equalwidth.png"),
       p_grid, width = 18, height = 9.5, units = "in", dpi = 1200)
