# ============================================================
# HUBnDEG_plots.R
# Publication-ready figures for drug proximity to HUBnDEG genes
# (intersect of DEGs with top-N most connected nodes, N = |DEGs|)
# ============================================================
# Loads: proximity_HUBnDEG_A7.tsv, proximity_HUBnDEG_A9.tsv
# Output: HUBnDEG_P*.pdf / HUBnDEG_P*.png
# ============================================================

setwd("~/CESC_Network/7_DGIDB/HUBnDEG/")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(patchwork)
  library(scales)
  library(ggridges)
  library(RColorBrewer)
})

# ── Palette & theme ───────────────────────────────────────────

PAL <- list(
  A7        = "#1B6CA8",
  A9        = "#C0392B",
  sig       = "#E74C3C",
  nonsig    = "#BDC3C7",
  gradient  = c("#08306B", "#2171B5", "#6BAED6", "#C6DBEF", "#FFF5EB",
                "#FDAE6B", "#E6550D", "#A63603"),
  highlight = "#F39C12"
)

theme_pub <- function(base = 12) {
  theme_minimal(base_size = base) +
    theme(
      plot.title        = element_text(face = "bold", size = base + 2, hjust = 0),
      plot.subtitle     = element_text(size = base - 1, color = "#555555", hjust = 0,
                                       margin = margin(b = 8)),
      axis.title        = element_text(face = "bold", size = base),
      axis.text         = element_text(color = "#333333"),
      panel.grid.major  = element_line(color = "#EBEBEB"),
      panel.grid.minor  = element_blank(),
      legend.title      = element_text(face = "bold", size = base - 1),
      legend.text       = element_text(size = base - 2),
      strip.text        = element_text(face = "bold", size = base),
      plot.background   = element_rect(fill = "white", color = NA),
      panel.background  = element_rect(fill = "white", color = NA),
      plot.margin       = margin(14, 14, 14, 14)
    )
}

save_fig <- function(p, name, w = 10, h = 7, dpi = 300) {
  ggsave(paste0(name, ".pdf"), p, width = w, height = h)
  ggsave(paste0(name, ".png"), p, width = w, height = h, dpi = dpi, bg = "white")
  cat("  Saved:", name, "\n")
}

# ── Load data ─────────────────────────────────────────────────

cat("[1] Loading HUBnDEG proximity results...\n")
prox_A7   <- read_tsv("proximity_HUBnDEG_A7.tsv", show_col_types = FALSE)
prox_A9   <- read_tsv("proximity_HUBnDEG_A9.tsv", show_col_types = FALSE)
prox_both <- bind_rows(prox_A7, prox_A9) %>%
  mutate(
    cohort    = factor(cohort, levels = c("A7", "A9")),
    sig_label = case_when(
      significant & z_score < -3  ~ "Significant (z < -3)",
      significant & z_score < -2  ~ "Significant (z < -2)",
      significant                  ~ "Significant (z < -1.96)",
      TRUE                         ~ "Not significant"
    ),
    sig_label = factor(sig_label,
                       levels = c("Significant (z < -3)", "Significant (z < -2)",
                                  "Significant (z < -1.96)", "Not significant"))
  )

# Gene-set sizes from data (inferred from proximity objects summary)
n_A7 <- prox_A7 %>% filter(!is.na(d_c_HUBnDEG)) %>% nrow()  # proxy count
n_A9 <- prox_A9 %>% filter(!is.na(d_c_HUBnDEG)) %>% nrow()
cat(sprintf("  A7: %d drugs | A9: %d drugs\n", nrow(prox_A7), nrow(prox_A9)))

top_label <- function(df, n = 12) {
  df %>%
    group_by(cohort) %>%
    slice_min(z_score, n = n, with_ties = FALSE) %>%
    ungroup()
}

# ══════════════════════════════════════════════════════════════
# PLOT 1 · Volcano: z-score vs -log10(p_value)
# ══════════════════════════════════════════════════════════════

cat("[P1] Volcano plot...\n")

p1_data <- prox_both %>%
  filter(!is.na(z_score), !is.na(p_value)) %>%
  mutate(log10p = -log10(p_value))

label_df <- top_label(p1_data %>% filter(significant), n = 10)

P1_XLIM <- c(-4.8, 4.2)
clip_annot <- p1_data %>%
  group_by(cohort) %>%
  summarise(n_clip = sum(z_score > 4.2, na.rm = TRUE), .groups = "drop") %>%
  mutate(z_score = 4.15, log10p = Inf,
         label = paste0(n_clip, " drugs\nz > 4.2\nnot shown\n(n ≤3 targets)"))

thr_196 <- -log10(pnorm(-1.96))
thr_2   <- -log10(pnorm(-2))
thr_3   <- -log10(pnorm(-3))

p1 <- ggplot(p1_data, aes(x = z_score, y = log10p, color = sig_label, size = n_targets)) +
  annotate("rect",
           xmin = P1_XLIM[1], xmax = -1.96, ymin = thr_196, ymax = Inf,
           fill = "#EAF4FB", alpha = 0.55) +
  geom_hline(yintercept = thr_196, linetype = "dashed", color = "#2980B9", linewidth = 0.7) +
  geom_hline(yintercept = thr_2,   linetype = "dotted", color = "#E74C3C", linewidth = 0.5) +
  geom_hline(yintercept = thr_3,   linetype = "dotted", color = "#C0392B", linewidth = 0.5) +
  geom_vline(xintercept = 0,     linetype = "solid",  color = "#BBBBBB", linewidth = 0.4) +
  geom_vline(xintercept = -1.96, linetype = "dashed", color = "#2980B9", linewidth = 0.7) +
  geom_point(alpha = 0.65) +
  geom_text_repel(
    data     = label_df,
    aes(label = drug),
    size     = 2.9, fontface = "bold",
    box.padding   = 0.45, point.padding = 0.3,
    max.overlaps  = 25, min.segment.length = 0.2,
    segment.size  = 0.3, segment.color = "#777777",
    xlim     = c(P1_XLIM[1], -0.5),
    show.legend   = FALSE
  ) +
  annotate("text", x = (P1_XLIM[1] + -1.96) / 2, y = Inf,
           label = "← significant candidates",
           hjust = 0.5, vjust = 1.6, size = 3.0, color = "#2471A3", fontface = "italic") +
  annotate("text", x = P1_XLIM[2], y = thr_196 + 0.05,
           label = "z = -1.96", hjust = 1, vjust = 0, size = 2.5, color = "#2980B9") +
  annotate("text", x = P1_XLIM[2], y = thr_3 + 0.05,
           label = "z = -3",    hjust = 1, vjust = 0, size = 2.5, color = "#C0392B") +
  geom_text(data = clip_annot,
            aes(x = z_score, y = log10p, label = label),
            hjust = 1, vjust = 1.4, size = 2.5, color = "#999999",
            lineheight = 1.1, inherit.aes = FALSE) +
  scale_color_manual(
    name   = "Significance (Guney z-threshold)",
    values = c("Significant (z < -3)"    = "#C0392B",
               "Significant (z < -2)"    = PAL$sig,
               "Significant (z < -1.96)" = PAL$highlight,
               "Not significant"         = PAL$nonsig)
  ) +
  scale_size_continuous(name = "# Targets in network", range = c(1.5, 6),
                        breaks = c(1, 5, 20, 50)) +
  scale_x_continuous(breaks = c(-4, -3, -2, -1.96, 0, 2, 4),
                     labels = c("-4", "-3", "-2", "-1.96", "0", "2", "4")) +
  scale_y_continuous(name = expression(-log[10](p[permutation]))) +
  coord_cartesian(xlim = P1_XLIM) +
  facet_wrap(~cohort, labeller = labeller(cohort = c(A7 = "HPV-A7", A9 = "HPV-A9"))) +
  labs(
    title    = "Drug proximity to HUBnDEG genes (DEG ∩ top-N connected)",
    subtitle = paste0(
      "x: permutation z-score (Guney et al. 2016); gene set: DEGs intersected with top-N most connected nodes\n",
      "y: -log10(one-sided permutation p-value); shaded quadrant: z < -1.96 and p < 0.025"
    ),
    x = "Permutation z-score"
  ) +
  theme_pub() +
  theme(panel.grid.major.x = element_line(color = "#F0F0F0"))

save_fig(p1, "HUBnDEG_P1_volcano", w = 14, h = 7)

# ══════════════════════════════════════════════════════════════
# PLOT 2 · Lollipop: top 20 drugs per cohort by z-score
# ══════════════════════════════════════════════════════════════

cat("[P2] Lollipop top drugs...\n")

top20 <- prox_both %>%
  filter(!is.na(z_score)) %>%
  group_by(cohort) %>%
  slice_min(z_score, n = 20, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    drug = str_trunc(drug, 30),
    drug = fct_reorder2(drug, cohort, z_score, .fun = function(c, z) -mean(z))
  )

p2 <- ggplot(top20, aes(x = z_score, y = drug, color = cohort)) +
  geom_vline(xintercept = c(-1.96, 0), linetype = "dotted",
             color = "#CCCCCC", linewidth = 0.5) +
  geom_segment(aes(xend = 0, yend = drug), linewidth = 0.8, alpha = 0.7) +
  geom_point(aes(size = n_targets), alpha = 0.9) +
  geom_text(aes(label = ifelse(p_adj < 0.05, "*", "")),
            hjust = 0, nudge_x = 0.05, size = 4, color = "#E74C3C", show.legend = FALSE) +
  scale_color_manual(name   = "HPV clade",
                     values = c(A7 = PAL$A7, A9 = PAL$A9)) +
  scale_size_continuous(name = "# Targets", range = c(2, 8),
                        breaks = c(1, 5, 20, 50)) +
  facet_wrap(~cohort, scales = "free_y",
             labeller = labeller(cohort = c(A7 = "HPV-A7", A9 = "HPV-A9"))) +
  labs(
    title    = "Top 20 candidate drugs by proximity z-score (HUBnDEG gene set)",
    subtitle = "Ranked by permutation z-score (more negative = closer to HUBnDEG genes than chance)\n* = BH-adjusted p < 0.05",
    x        = "Permutation z-score",
    y        = NULL
  ) +
  theme_pub() +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 9))

save_fig(p2, "HUBnDEG_P2_lollipop_top20", w = 14, h = 10)

# ══════════════════════════════════════════════════════════════
# PLOT 3 · Scatter: d_c_HUBnDEG vs d_c_nonHUBnDEG
# ══════════════════════════════════════════════════════════════

cat("[P3] Distance scatter...\n")

scatter_df <- prox_both %>%
  filter(!is.na(d_c_HUBnDEG), !is.na(d_c_nonHUBnDEG), !is.na(z_score))

P3_XLIM <- c(-0.02, quantile(scatter_df$d_c_HUBnDEG,    0.98, na.rm = TRUE) + 0.05)
P3_YLIM <- c(-0.02, quantile(scatter_df$d_c_nonHUBnDEG, 0.98, na.rm = TRUE) + 0.05)
shared_max <- max(P3_XLIM[2], P3_YLIM[2])
P3_XLIM[2] <- shared_max
P3_YLIM[2] <- shared_max

lab_scatter <- scatter_df %>%
  filter(significant) %>%
  group_by(cohort) %>%
  slice_min(z_score, n = 8) %>%
  ungroup()

set.seed(42)

p3 <- ggplot(scatter_df, aes(x = d_c_HUBnDEG, y = d_c_nonHUBnDEG)) +
  annotate("polygon",
           x = c(0, P3_XLIM[2], P3_XLIM[2], 0),
           y = c(0, P3_XLIM[2], P3_YLIM[2], P3_YLIM[2]),
           fill = "#EBF5FB", alpha = 0.6) +
  annotate("text",
           x = P3_XLIM[2] - 0.02, y = P3_YLIM[2] - 0.02,
           label = "Absolute preference\nfor HUBnDEG genes",
           hjust = 1, vjust = 1, size = 2.8, color = "#2471A3",
           fontface = "italic", lineheight = 1.1) +
  geom_abline(slope = 1, intercept = 0,
              color = "#1A252F", linewidth = 1.1, linetype = "solid") +
  geom_jitter(aes(color = z_score, size = n_targets),
              width = 0.008, height = 0.005,
              alpha = 0.55, shape = 16) +
  geom_text_repel(
    data        = lab_scatter %>% filter(d_c_HUBnDEG    <= P3_XLIM[2],
                                         d_c_nonHUBnDEG <= P3_YLIM[2]),
    aes(label   = drug),
    size        = 2.8, fontface = "bold",
    box.padding = 0.45, point.padding = 0.25,
    max.overlaps = 18,
    segment.size = 0.3, segment.color = "#555555",
    show.legend = FALSE
  ) +
  scale_color_gradientn(
    name   = "Permutation\nz-score",
    colors = rev(brewer.pal(9, "RdYlBu")),
    limits = c(-4, 4), oob = squish,
    breaks = c(-4, -2, 0, 2, 4),
    guide  = guide_colorbar(barwidth = 0.8, barheight = 5,
                            ticks.colour = "grey40", frame.colour = "grey40")
  ) +
  scale_size_continuous(name = "# Targets", range = c(0.8, 5.5),
                        breaks = c(1, 5, 20, 50)) +
  coord_fixed(ratio = 1, xlim = P3_XLIM, ylim = P3_YLIM, expand = FALSE) +
  facet_wrap(~cohort, labeller = labeller(cohort = c(A7 = "HPV-A7", A9 = "HPV-A9"))) +
  labs(
    title    = "Drug proximity: HUBnDEG genes vs non-HUBnDEG genes",
    subtitle = paste0(
      "Each point is a drug. Solid diagonal (y = x): equal distance to both gene sets.\n",
      "Shaded region (above diagonal): drug targets are absolutely closer to HUBnDEG genes; ",
      "color encodes significance (blue = stronger candidate)."
    ),
    x = expression(italic(d)[c] * " to HUBnDEG genes  (closer →  left)"),
    y = expression(italic(d)[c] * " to non-HUBnDEG genes  (closer →  bottom)")
  ) +
  theme_pub() +
  theme(panel.border = element_rect(color = "#CCCCCC", fill = NA, linewidth = 0.5))

save_fig(p3, "HUBnDEG_P3_scatter_dist", w = 14, h = 7)

# ══════════════════════════════════════════════════════════════
# PLOT 4 · Z-score distribution across all drugs (density + rug)
# ══════════════════════════════════════════════════════════════

cat("[P4] Z-score distribution...\n")

dist_df <- prox_both %>% filter(!is.na(z_score))

sig_counts <- dist_df %>%
  group_by(cohort) %>%
  summarise(
    n_total = n(),
    n_sig   = sum(significant,        na.rm = TRUE),
    n_sig2  = sum(sig_z2,             na.rm = TRUE),
    n_sig3  = sum(sig_tier == "z<-3", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    lbl = paste0(
      "z < -1.96: ", n_sig,  " drugs\n",
      "z < -2:     ", n_sig2, " drugs\n",
      "z < -3:     ", n_sig3, " drugs"
    )
  )

p4 <- ggplot(dist_df, aes(x = z_score, fill = cohort, color = cohort)) +
  annotate("rect", xmin = -Inf, xmax = -1.96, ymin = -Inf, ymax = Inf,
           fill = "#EAF4FB", alpha = 0.50) +
  geom_density(alpha = 0.30, linewidth = 0.9, trim = FALSE) +
  geom_rug(
    data  = dist_df %>% filter(significant),
    aes(x = z_score, color = cohort),
    sides = "b", length = unit(0.04, "npc"), linewidth = 0.4, alpha = 0.7
  ) +
  geom_vline(xintercept = 0,     linetype = "solid",  color = "#BBBBBB", linewidth = 0.5) +
  geom_vline(xintercept = -1.96, linetype = "dashed", color = "#2980B9", linewidth = 0.8) +
  geom_vline(xintercept = -2,    linetype = "dotted", color = "#E74C3C", linewidth = 0.6) +
  geom_vline(xintercept = -3,    linetype = "dotted", color = "#C0392B", linewidth = 0.6) +
  geom_text(
    data = sig_counts %>% mutate(xpos = ifelse(cohort == "A7", -4.6, -2.2)),
    aes(x = xpos, y = Inf, label = lbl, color = cohort),
    hjust = 0, vjust = 1.5, size = 2.8, lineheight = 1.3,
    fontface = "plain", inherit.aes = FALSE
  ) +
  scale_fill_manual(values  = c(A7 = PAL$A7, A9 = PAL$A9),
                    name    = "HPV clade",
                    labels  = c(A7 = "HPV-A7", A9 = "HPV-A9")) +
  scale_color_manual(values = c(A7 = PAL$A7, A9 = PAL$A9),
                     name   = "HPV clade",
                     labels = c(A7 = "HPV-A7", A9 = "HPV-A9")) +
  coord_cartesian(xlim = c(-4.8, 4.2)) +
  scale_x_continuous(breaks = c(-4, -3, -2, -1.96, 0, 2, 4),
                     labels = c("-4", "-3", "-2", "-1.96", "0", "2", "4")) +
  labs(
    title    = "Genome-wide drug proximity z-score distribution (HUBnDEG gene set)",
    subtitle = paste0(
      "Density of permutation z-scores across all tested drugs; shaded region: z < -1.96\n",
      "Rug marks: significant drugs (z < -1.96) per cohort"
    ),
    x = "Permutation z-score",
    y = "Density"
  ) +
  theme_pub() +
  theme(panel.grid.major.x = element_line(color = "#F0F0F0"))

save_fig(p4, "HUBnDEG_P4_zscore_distribution", w = 10, h = 6)

# ══════════════════════════════════════════════════════════════
# PLOT 5 · Cross-cohort comparison (A7 vs A9 common drugs)
# ══════════════════════════════════════════════════════════════

cat("[P5] A7 vs A9 cross-cohort scatter...\n")

common_drugs <- intersect(prox_A7$drug, prox_A9$drug)

cross <- prox_A7 %>%
  filter(drug %in% common_drugs) %>%
  select(drug, n_targets_A7 = n_targets,
         z_A7 = z_score, sig_A7 = significant) %>%
  inner_join(
    prox_A9 %>%
      filter(drug %in% common_drugs) %>%
      select(drug, n_targets_A9 = n_targets,
             z_A9 = z_score, sig_A9 = significant),
    by = "drug"
  ) %>%
  mutate(
    sig_both = sig_A7 & sig_A9,
    sig_any  = sig_A7 | sig_A9,
    grp = case_when(
      sig_A7 & sig_A9    ~ "Both",
      sig_A7 & !sig_A9   ~ "A7 only",
      !sig_A7 & sig_A9   ~ "A9 only",
      TRUE               ~ "Neither"
    ),
    grp = factor(grp, levels = c("Both", "A7 only", "A9 only", "Neither"))
  )

cross_label <- cross %>%
  filter(sig_any) %>%
  slice_min(z_A7 + z_A9, n = 10, with_ties = FALSE)

p5 <- ggplot(cross, aes(x = z_A7, y = z_A9, color = grp)) +
  geom_hline(yintercept = c(-1.96, 0), linetype = "dotted", color = "#CCCCCC") +
  geom_vline(xintercept = c(-1.96, 0), linetype = "dotted", color = "#CCCCCC") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "#AAAAAA", linewidth = 0.5) +
  geom_point(aes(size = (n_targets_A7 + n_targets_A9) / 2), alpha = 0.7) +
  geom_text_repel(
    data        = cross_label,
    aes(label   = drug),
    size        = 3, fontface = "bold",
    box.padding = 0.5, max.overlaps = 15,
    segment.size = 0.25, segment.color = "#666666",
    show.legend = FALSE
  ) +
  scale_color_manual(
    name   = "Significant in",
    values = c("Both"    = "#8E44AD",
               "A7 only" = PAL$A7,
               "A9 only" = PAL$A9,
               "Neither" = PAL$nonsig)
  ) +
  scale_size_continuous(name = "Avg # targets", range = c(1.5, 7),
                        breaks = c(1, 5, 20, 50), guide = "none") +
  annotate("rect", xmin = -Inf, xmax = -1.96, ymin = -Inf, ymax = -1.96,
           fill = "#8E44AD", alpha = 0.04) +
  labs(
    title    = "Cross-cohort drug proximity comparison (HUBnDEG gene set)",
    subtitle = paste0(length(common_drugs),
                      " drugs present in both A7 and A9 networks\n",
                      "Shaded quadrant: significant in both cohorts"),
    x = "Permutation z-score (HPV-A7)",
    y = "Permutation z-score (HPV-A9)"
  ) +
  theme_pub()

save_fig(p5, "HUBnDEG_P5_cross_cohort", w = 9, h = 8)

# ══════════════════════════════════════════════════════════════
# PLOT 6 · Ranked bubble plot: top 30 per cohort
# ══════════════════════════════════════════════════════════════

cat("[P6] Ranked bubble plot...\n")

top30 <- prox_both %>%
  filter(!is.na(z_score)) %>%
  group_by(cohort) %>%
  slice_min(z_score, n = 30, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    drug  = str_trunc(drug, 28),
    drug  = fct_reorder(drug, z_score, .desc = FALSE),
    stars = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ ""
    )
  )

p6 <- ggplot(top30, aes(x = z_score, y = drug)) +
  geom_vline(xintercept = c(-1.96, 0), linetype = "dotted",
             color = "#CCCCCC", linewidth = 0.5) +
  geom_segment(aes(xend = 0, yend = drug, color = -log10(p_adj + 1e-10)),
               linewidth = 1.2, alpha = 0.8) +
  geom_point(aes(size = n_targets, fill = -log10(p_adj + 1e-10)),
             shape = 21, color = "white", stroke = 0.4, alpha = 0.95) +
  geom_text(aes(label = stars, x = z_score - 0.1),
            hjust = 1, vjust = 0.5, size = 3.5, color = "#333333") +
  scale_fill_gradientn(
    name   = expression(-log[10](p[adj])),
    colors = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
    limits = c(0, 3), oob = squish
  ) +
  scale_color_gradientn(
    colors = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
    limits = c(0, 3), oob = squish,
    guide  = "none"
  ) +
  scale_size_continuous(name = "# Targets", range = c(2, 8),
                        breaks = c(1, 5, 20, 50)) +
  facet_wrap(~cohort, scales = "free_y",
             labeller = labeller(cohort = c(A7 = "HPV-A7", A9 = "HPV-A9"))) +
  labs(
    title    = "Top 30 candidate drugs by network proximity (HUBnDEG gene set)",
    subtitle = expression(paste(
      "Color/segment: BH-adjusted significance (-log"[10], "(p"[adj], "))  |  Size: # targets  |  *, **, *** = p"[adj], " < 0.05/0.01/0.001"
    )),
    x = "Permutation z-score",
    y = NULL
  ) +
  theme_pub(base = 11) +
  theme(axis.text.y = element_text(size = 8.5))

save_fig(p6, "HUBnDEG_P6_bubble_top30", w = 16, h = 12)

# ══════════════════════════════════════════════════════════════
# PLOT 7 · Null distribution for top 5 drugs per cohort
# ══════════════════════════════════════════════════════════════

cat("[P7] Null distribution ridge plot...\n")

top5_A7 <- prox_A7 %>%
  filter(significant) %>%
  slice_min(z_score, n = 5, with_ties = FALSE)

top5_A9 <- prox_A9 %>%
  filter(significant) %>%
  slice_min(z_score, n = 5, with_ties = FALSE)

make_ridge_data <- function(df, cohort_lbl) {
  purrr::map_dfr(seq_len(nrow(df)), function(i) {
    mu  <- df$mean_null[i]
    sig <- df$sd_null[i]
    obs <- df$d_c_HUBnDEG[i]
    tibble(
      drug   = df$drug[i],
      x      = seq(mu - 4 * sig, mu + 4 * sig, length.out = 300),
      y      = dnorm(seq(mu - 4 * sig, mu + 4 * sig, length.out = 300), mu, sig),
      obs    = obs,
      z      = df$z_score[i],
      cohort = cohort_lbl
    )
  })
}

ridge_df <- bind_rows(
  make_ridge_data(top5_A7, "A7"),
  make_ridge_data(top5_A9, "A9")
) %>%
  mutate(drug = str_trunc(drug, 28),
         drug = factor(drug))

obs_pts <- ridge_df %>%
  distinct(drug, obs, z, cohort) %>%
  left_join(ridge_df %>%
              group_by(drug) %>%
              filter(x == x[which.min(abs(x - obs))]) %>%
              summarise(y_obs = max(y), .groups = "drop"),
            by = "drug")

p7 <- ggplot(ridge_df, aes(x = x, y = drug, height = y, group = drug)) +
  geom_density_ridges(
    stat = "identity", scale = 1.4,
    fill = "#AED6F1", color = "#2980B9", alpha = 0.7, linewidth = 0.6
  ) +
  geom_segment(
    data = obs_pts,
    aes(x = obs, xend = obs,
        y = as.numeric(drug), yend = as.numeric(drug) + 0.9),
    color = PAL$sig, linewidth = 1.2, inherit.aes = FALSE
  ) +
  geom_text(
    data = obs_pts,
    aes(x = obs, y = as.numeric(drug) + 0.95,
        label = sprintf("z=%.2f", z)),
    color = PAL$sig, size = 3, hjust = 0.5, fontface = "bold", inherit.aes = FALSE
  ) +
  facet_wrap(~cohort, scales = "free",
             labeller = labeller(cohort = c(A7 = "HPV-A7", A9 = "HPV-A9"))) +
  labs(
    title    = "Null distributions for top significant drugs (HUBnDEG gene set)",
    subtitle = "Blue curve: permutation null (Gaussian approximation from 500 permutations)\nRed line: observed d_c(drug targets -> HUBnDEG genes)",
    x        = expression(italic(d)[c] * " (closest proximity to HUBnDEG genes)"),
    y        = NULL
  ) +
  theme_pub()

save_fig(p7, "HUBnDEG_P7_null_ridges", w = 14, h = 7)

# ══════════════════════════════════════════════════════════════
# PLOT 8 · Summary table: significant drugs both cohorts
# ══════════════════════════════════════════════════════════════

cat("[P8] Summary stats table (significant drugs)...\n")

sig_summary <- prox_both %>%
  filter(significant) %>%
  select(cohort, drug, n_targets, d_c_HUBnDEG, d_c_nonHUBnDEG, z_score, p_adj) %>%
  arrange(cohort, z_score)

write_tsv(sig_summary, "HUBnDEG_significant_drugs_summary.tsv")
cat("  Saved: HUBnDEG_significant_drugs_summary.tsv\n")

cat("\nSignificant drugs (z < -1.96):\n")
print(sig_summary, n = 40)

cat("\n[Done] All HUBnDEG plots saved.\n")
cat("Files:\n")
for (f in list.files(".", pattern = "^HUBnDEG_P[0-9]_.*\\.(png|pdf)$")) cat("  ", f, "\n")
