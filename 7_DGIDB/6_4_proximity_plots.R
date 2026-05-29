# ============================================================
# 6_4_proximity_plots.R
# Publication-ready figures — drug proximity to interest genes
# ============================================================
# Method: Closest proximity d_c (Guney et al. 2016)
#
# Two independent dimensions per drug:
#   d_c_interest  → actual network distance (hops); lower = physically closer
#   z_score       → how much closer than expected by chance (permutation)
#                   these are NOT redundant: a drug can be close but not
#                   statistically significant (many targets), or significant
#                   but at moderate distance (few very well-placed targets)
#   p_value = pnorm(z_score) → same info as z, NOT plotted separately
#
# Input:  6_3_proximity_A7.tsv, 6_3_proximity_A9.tsv
# Output: 6_4_P*.pdf / 6_4_P*.png
#         6_4_significant_drugs_summary.tsv
#
# Plots
#   P1 · Z-score distribution  — landscape of all tested drugs
#   P2 · d_c vs z-score        — raw distance × statistical significance
#   P3 · Ranked by z-score     — top 20 candidates per cohort
#   P4 · Ranked by d_c         — top 20 closest (by actual distance) per cohort
#   P5 · Cross-cohort          — HPV-A7 vs HPV-A9 z-scores
# ============================================================

setwd("~/CESC_Network/7_DGIDB/")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(patchwork)
  library(scales)
})

# ── Palette & helpers ──────────────────────────────────────────

PAL_A7 <- "#1B6CA8"
PAL_A9 <- "#A93226"
GRAY   <- "#B2BABB"

TIER_COLORS <- c(
  "z < -3"          = "#6C3483",
  "z < -2"          = "#A93226",
  "z < -1.96"       = "#CA6F1E",
  "Not significant"  = "#B2BABB"
)

COHORT_LABELS <- c(A7 = "HPV-A7", A9 = "HPV-A9")

tier_label <- function(z) {
  case_when(
    z < -3    ~ "z < -3",
    z < -2    ~ "z < -2",
    z < -1.96 ~ "z < -1.96",
    TRUE      ~ "Not significant"
  )
}

# Asymmetric y-scale: identity below thresh, compressed above
make_zsqueeze_trans <- function(thresh = -1.96, k = 0.22) {
  scales::trans_new(
    "z_squeeze",
    transform = function(x) ifelse(x <= thresh, x, thresh + (x - thresh) * k),
    inverse   = function(x) ifelse(x <= thresh, x, thresh + (x - thresh) / k)
  )
}

# Asymmetric x-scale: identity below thresh (interesting region), compressed above
make_xsqueeze_trans <- function(thresh = 1.5, k = 0.28) {
  scales::trans_new(
    "x_squeeze",
    transform = function(x) ifelse(x <= thresh, x, thresh + (x - thresh) * k),
    inverse   = function(x) ifelse(x <= thresh, x, thresh + (x - thresh) / k)
  )
}

theme_pub <- function(base = 12) {
  theme_minimal(base_size = base) +
    theme(
      plot.title        = element_text(face = "bold", size = base + 3, hjust = 0,
                                       margin = margin(b = 5)),
      plot.subtitle     = element_text(size = base - 1.5, color = "#555555", hjust = 0,
                                       margin = margin(b = 10), lineheight = 1.35),
      axis.title        = element_text(face = "bold", size = base),
      axis.text         = element_text(color = "#333333", size = base - 1),
      panel.grid.major  = element_line(color = "#EFEFEF", linewidth = 0.4),
      panel.grid.minor  = element_blank(),
      legend.title      = element_text(face = "bold", size = base - 1),
      legend.text       = element_text(size = base - 2),
      strip.text        = element_text(face = "bold", size = base + 1, color = "#222222"),
      strip.background  = element_rect(fill = "#F4F4F4", color = NA),
      panel.border      = element_rect(color = "#DDDDDD", fill = NA, linewidth = 0.5),
      plot.background   = element_rect(fill = "white", color = NA),
      panel.background  = element_rect(fill = "white", color = NA),
      plot.margin       = margin(16, 18, 12, 16)
    )
}

save_fig <- function(p, name, w = 10, h = 7, dpi = 300) {
  ggsave(paste0(name, ".pdf"), p, width = w, height = h)
  ggsave(paste0(name, ".png"), p, width = w, height = h, dpi = dpi, bg = "white")
  message("  Saved: ", name)
}

# ── Load data ─────────────────────────────────────────────────

message("[1] Loading data...")
prox_A7 <- read_tsv("6_3_proximity_A7.tsv", show_col_types = FALSE)
prox_A9 <- read_tsv("6_3_proximity_A9.tsv", show_col_types = FALSE)

prox_both <- bind_rows(prox_A7, prox_A9) %>%
  filter(!is.na(z_score), !is.na(d_c_interest)) %>%
  mutate(
    cohort = factor(cohort, levels = c("A7", "A9")),
    tier   = factor(tier_label(z_score), levels = names(TIER_COLORS))
  )

Z_LIM <- c(-5.5, 5.0)

# ══════════════════════════════════════════════════════════════
# PLOT 1 · Z-score distribution — landscape of all drugs
# ══════════════════════════════════════════════════════════════

message("[P1] Z-score distribution...")

n_sig_tbl <- prox_both %>%
  group_by(cohort) %>%
  summarise(n = n(), n_sig = sum(significant), .groups = "drop") %>%
  mutate(pct = round(100 * n_sig / n, 1))

sub_p1 <- paste(
  sprintf("%s: %d / %d drugs significant at z < -1.96 (%.1f%%)",
          COHORT_LABELS[as.character(n_sig_tbl$cohort)],
          n_sig_tbl$n_sig, n_sig_tbl$n, n_sig_tbl$pct),
  collapse = "   |   "
)

p1 <- ggplot(prox_both, aes(x = z_score, color = cohort, fill = cohort)) +
  annotate("rect",
           xmin = -Inf, xmax = -1.96, ymin = 0, ymax = Inf,
           fill = "#D6EAF8", alpha = 0.40) +
  geom_density(alpha = 0.18, linewidth = 0.95) +
  geom_rug(
    data = filter(prox_both, significant),
    sides = "b", length = unit(0.04, "npc"),
    linewidth = 0.35, alpha = 0.70
  ) +
  geom_vline(xintercept = 0,     color = "#CCCCCC", linewidth = 0.45) +
  geom_vline(xintercept = -1.96, color = "#2471A3", linewidth = 0.85, linetype = "dashed") +
  geom_vline(xintercept = -3,    color = "#6C3483", linewidth = 0.65, linetype = "dotted") +
  annotate("text", x = -2.10, y = Inf, label = "z = -1.96",
           angle = 90, hjust = 1.15, vjust = -0.3, size = 2.8, color = "#2471A3") +
  annotate("text", x = -3.12, y = Inf, label = "z = -3",
           angle = 90, hjust = 1.15, vjust = -0.3, size = 2.8, color = "#6C3483") +
  scale_fill_manual(values = c(A7 = PAL_A7, A9 = PAL_A9),
                    labels = COHORT_LABELS, name = "HPV clade") +
  scale_color_manual(values = c(A7 = PAL_A7, A9 = PAL_A9),
                     labels = COHORT_LABELS, name = "HPV clade") +
  coord_cartesian(xlim = Z_LIM) +
  scale_x_continuous(breaks = c(-4, -3, -2, -1.96, 0, 2, 4),
                     labels = c("-4", "-3", "-2", "-1.96", "0", "2", "4")) +
  labs(
    title    = "Distribution of permutation z-scores across all tested drugs",
    subtitle = sub_p1,
    x        = "Permutation z-score",
    y        = "Density"
  ) +
  theme_pub() +
  theme(legend.position = c(0.88, 0.78))

save_fig(p1, "6_4_P1_distribution", w = 10, h = 6)

# ══════════════════════════════════════════════════════════════
# PLOT 2 · d_c_interest vs z_score — two independent axes
# ══════════════════════════════════════════════════════════════
# d_c: actual network distance (lower = physically closer)
# z  : normalized by permutation null (more negative = rarer)
# These are NOT the same: a drug can be close (low d_c) but
# not significant (z > -1.96) if it has many network targets
# and the null distribution shifts accordingly.
# Best candidates: lower-left quadrant (close AND significant)
# ══════════════════════════════════════════════════════════════

message("[P2] d_c vs z-score scatter...")

# Reference: median d_c of significant drugs per cohort
dc_med <- prox_both %>%
  filter(significant) %>%
  group_by(cohort) %>%
  summarise(med_dc = median(d_c_interest), .groups = "drop")

labels_p2 <- prox_both %>%
  filter(significant) %>%
  group_by(cohort) %>%
  slice_min(d_c_interest + z_score, n = 10, with_ties = FALSE) %>%  # closest AND most significant
  ungroup()

p2 <- ggplot(prox_both, aes(x = d_c_interest, y = z_score)) +
  # Shading: candidate quadrant (significant z AND below median d_c of significant drugs)
  geom_rect(
    data = dc_med,
    aes(xmin = -Inf, xmax = med_dc, ymin = -Inf, ymax = -1.96),
    fill = "#D5F5E3", alpha = 0.45, inherit.aes = FALSE
  ) +
  annotate("text", x = -Inf, y = -Inf,
           label = "Closest + significant\n(best candidates)",
           hjust = -0.05, vjust = -0.2, size = 2.7,
           color = "#1E8449", fontface = "italic") +
  # Significance threshold
  geom_hline(yintercept = -1.96, color = "#2471A3", linewidth = 0.80, linetype = "dashed") +
  geom_hline(yintercept = -3,    color = "#6C3483", linewidth = 0.60, linetype = "dotted") +
  geom_hline(yintercept = 0,     color = "#CCCCCC", linewidth = 0.45) +
  # Median d_c of significant drugs (per-facet reference)
  geom_vline(
    data = dc_med,
    aes(xintercept = med_dc),
    color = "#1E8449", linewidth = 0.65, linetype = "dashed"
  ) +
  # Not-significant background
  geom_point(
    data  = function(d) filter(d, tier == "Not significant"),
    color = GRAY, size = 1.1, alpha = 0.30, shape = 16
  ) +
  # Significant points on top
  geom_point(
    data = function(d) filter(d, tier != "Not significant"),
    aes(color = tier, size = n_targets), alpha = 0.85, shape = 16
  ) +
  geom_text_repel(
    data          = labels_p2,
    aes(label = drug, color = tier),
    size          = 2.7, fontface = "bold",
    box.padding   = 0.40, point.padding = 0.25,
    max.overlaps  = 20, min.segment.length = 0.15,
    segment.size  = 0.28, segment.color = "#BBBBBB",
    show.legend   = FALSE
  ) +
  annotate("text", x = Inf, y = -2.05,
           label = "z = -1.96", hjust = 1.05, size = 2.5, color = "#2471A3") +
  scale_color_manual(name = "Significance", values = TIER_COLORS, drop = FALSE) +
  scale_size_continuous(name = "# Network\ntargets",
                        range = c(1.2, 6), breaks = c(1, 5, 20, 50)) +
  scale_y_continuous(
    trans  = make_zsqueeze_trans(),
    limits = c(-6.5, 5),
    breaks = c(-5, -4, -3, -2, -1.96, 0, 2, 4),
    labels = c("-5", "-4", "-3", "-2", "-1.96", "0", "2", "4")
  ) +
  scale_x_continuous(
    trans  = make_xsqueeze_trans(),
    breaks = c(0.5, 1.0, 1.2, 1.4, 2.0, 3.0),
    labels = c("0.5", "1.0", "1.2", "1.4", "2.0", "3.0"),
    expand = expansion(mult = c(0.03, 0.04))
  ) +
  facet_wrap(~cohort, labeller = labeller(cohort = COHORT_LABELS)) +
  labs(
    title    = "Network distance vs. statistical significance: all drugs",
    subtitle = paste0(
      "x-axis: d_c(interest genes) -- actual mean-minimum network distance (hops); lower = physically closer.\n",
      "y-axis: permutation z-score -- how much closer than random expectation; more negative = rarer.\n",
      "Green zone: d_c below median of significant drugs AND z < -1.96. ",
      "Dashed green line: median d_c of significant drugs per cohort."
    ),
    x = expression(italic(d)[c] * " to interest genes  (network hops)"),
    y = "Permutation z-score"
  ) +
  theme_pub() +
  guides(color = guide_legend(override.aes = list(size = 3.2)))

save_fig(p2, "6_4_P2_dc_vs_zscore", w = 14, h = 7)

# ══════════════════════════════════════════════════════════════
# PLOT 3 · Top 20 by z-score — statistical significance ranking
# ══════════════════════════════════════════════════════════════

message("[P3] Top candidates by z-score...")

TOP_N <- 20

make_lollipop_z <- function(df, panel_title) {
  dat <- df %>%
    filter(!is.na(z_score)) %>%
    slice_min(z_score, n = TOP_N, with_ties = FALSE) %>%
    mutate(
      drug = str_trunc(drug, 32),
      drug = fct_reorder(drug, z_score, .desc = TRUE),
      tier = factor(tier_label(z_score), levels = names(TIER_COLORS)[1:3])
    )

  ggplot(dat, aes(x = z_score, y = drug)) +
    geom_vline(xintercept = -1.96, linetype = "dashed",
               color = "#2471A3", linewidth = 0.65) +
    geom_vline(xintercept = 0, color = "#DDDDDD", linewidth = 0.40) +
    geom_segment(aes(xend = 0, yend = drug, color = tier),
                 linewidth = 1.1, alpha = 0.80) +
    geom_point(aes(color = tier, size = n_targets), alpha = 0.95) +
    scale_color_manual(name = "z-score tier",
                       values = TIER_COLORS[1:3], drop = FALSE) +
    scale_size_continuous(name = "# Network\ntargets",
                          range = c(2.5, 8), breaks = c(1, 5, 20, 50)) +
    scale_x_continuous(breaks = c(-4, -3, -2, -1.96, 0),
                       labels = c("-4", "-3", "-2", "-1.96", "0")) +
    labs(title = panel_title, x = "Permutation z-score", y = NULL) +
    theme_pub() +
    theme(axis.text.y = element_text(size = 9))
}

p3 <- (make_lollipop_z(prox_A7, "HPV-A7") |
       make_lollipop_z(prox_A9, "HPV-A9")) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = paste0("Top ", TOP_N, " drugs by permutation z-score (statistical proximity)"),
    subtitle = paste0(
      "Ranking by how much closer drug targets are to interest genes than random expectation.\n",
      "Point size = number of drug targets in the network. All drugs shown have z < -1.96."
    ),
    theme = theme_pub() + theme(
      plot.title    = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10.5, color = "#555555", lineheight = 1.35)
    )
  )

save_fig(p3, "6_4_P3_top_zscore", w = 14, h = 10)

# ══════════════════════════════════════════════════════════════
# PLOT 4 · Top 20 by d_c — physical distance ranking
# ══════════════════════════════════════════════════════════════
# Only among significant drugs (z < -1.96).
# This answers: "which candidates are literally closest to
# the interest genes in the network?"
# ══════════════════════════════════════════════════════════════

message("[P4] Top candidates by d_c (physical distance)...")

make_lollipop_dc <- function(df, panel_title) {
  dat <- df %>%
    filter(significant, !is.na(d_c_interest)) %>%
    slice_min(d_c_interest, n = TOP_N, with_ties = FALSE) %>%
    mutate(
      drug = str_trunc(drug, 32),
      drug = fct_reorder(drug, d_c_interest),   # smallest d_c at top
      tier = factor(tier_label(z_score), levels = names(TIER_COLORS)[1:3])
    )

  ggplot(dat, aes(x = d_c_interest, y = drug)) +
    geom_vline(xintercept = 0, color = "#DDDDDD", linewidth = 0.40) +
    geom_segment(aes(xend = 0, yend = drug, color = tier),
                 linewidth = 1.1, alpha = 0.80) +
    geom_point(aes(color = tier, size = n_targets), alpha = 0.95) +
    geom_text(aes(label = sprintf("z = %.2f", z_score)),
              hjust = -0.15, size = 2.5, color = "#555555") +
    scale_color_manual(name = "z-score tier",
                       values = TIER_COLORS[1:3], drop = FALSE) +
    scale_size_continuous(name = "# Network\ntargets",
                          range = c(2.5, 8), breaks = c(1, 5, 20, 50)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(
      title = panel_title,
      x     = expression(italic(d)[c] * " to interest genes  (network hops)"),
      y     = NULL
    ) +
    theme_pub() +
    theme(axis.text.y = element_text(size = 9))
}

p4 <- (make_lollipop_dc(prox_A7, "HPV-A7") |
       make_lollipop_dc(prox_A9, "HPV-A9")) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = paste0("Top ", TOP_N, " significant drugs by physical network distance"),
    subtitle = paste0(
      "Ranked by d_c(interest genes): the mean-minimum distance from drug targets to interest genes.\n",
      "Only drugs significant at z < -1.96 are included. ",
      "Labels show the permutation z-score for each drug."
    ),
    theme = theme_pub() + theme(
      plot.title    = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10.5, color = "#555555", lineheight = 1.35)
    )
  )

save_fig(p4, "6_4_P4_top_dc", w = 14, h = 10)

# ══════════════════════════════════════════════════════════════
# PLOT 5 · Cross-cohort — HPV-A7 vs HPV-A9 z-scores
# ══════════════════════════════════════════════════════════════

message("[P5] Cross-cohort...")

cross <- inner_join(
  prox_A7 %>% select(drug, z_A7 = z_score, n_A7 = n_targets, sig_A7 = significant),
  prox_A9 %>% select(drug, z_A9 = z_score, n_A9 = n_targets, sig_A9 = significant),
  by = "drug"
) %>%
  filter(!is.na(z_A7), !is.na(z_A9)) %>%
  mutate(
    grp = case_when(
      sig_A7 & sig_A9 ~ "Both cohorts",
      sig_A7          ~ "HPV-A7 only",
      sig_A9          ~ "HPV-A9 only",
      TRUE            ~ "Neither"
    ),
    grp = factor(grp, levels = c("Both cohorts", "HPV-A7 only", "HPV-A9 only", "Neither"))
  )

# Balanced labeling: top N per group, ranked by the axis that is significant
cross_lbl <- bind_rows(
  cross %>% filter(grp == "Both cohorts") %>%
    slice_min(z_A7 + z_A9, n = 6,  with_ties = FALSE),
  cross %>% filter(grp == "HPV-A7 only") %>%
    slice_min(z_A7,         n = 9,  with_ties = FALSE),
  cross %>% filter(grp == "HPV-A9 only") %>%
    slice_min(z_A9,         n = 9,  with_ties = FALSE)
)

n_both   <- sum(cross$grp == "Both cohorts")
n_common <- nrow(cross)

p5 <- ggplot(cross, aes(x = z_A7, y = z_A9)) +
  annotate("rect",
           xmin = -Inf, xmax = -1.96, ymin = -Inf, ymax = -1.96,
           fill = "#E8DAEF", alpha = 0.55) +
  annotate("label", x = -5.3, y = -2.1,
           label = paste0("Both cohorts\nn = ", n_both),
           hjust = 0, vjust = 1, size = 5.0,
           color = "#6C3483", fontface = "italic",
           fill = "#F5EEF8", label.size = 0.3, label.padding = unit(0.3, "lines")) +
  geom_vline(xintercept = 0,     color = "#DDDDDD", linewidth = 0.40) +
  geom_vline(xintercept = -1.96, color = "#2471A3", linewidth = 0.78, linetype = "dashed") +
  geom_hline(yintercept = 0,     color = "#DDDDDD", linewidth = 0.40) +
  geom_hline(yintercept = -1.96, color = "#2471A3", linewidth = 0.78, linetype = "dashed") +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "#CCCCCC", linewidth = 0.50) +
  geom_point(
    data  = function(d) filter(d, grp == "Neither"),
    color = GRAY, size = 1.2, alpha = 0.28, shape = 16
  ) +
  geom_point(
    data = function(d) filter(d, grp != "Neither"),
    aes(color = grp, size = (n_A7 + n_A9) / 2),
    alpha = 0.85, shape = 16
  ) +
  geom_text_repel(
    data              = cross_lbl,
    aes(label = drug, color = grp),
    size              = 4.5, fontface = "bold",
    box.padding       = 0.70, point.padding = 0.40,
    force             = 7,   force_pull = 0.3,
    max.overlaps      = Inf, min.segment.length = 0.08,
    segment.size      = 0.35, segment.color = "#AAAAAA",
    show.legend       = FALSE
  ) +
  scale_color_manual(
    name   = "Significant in",
    values = c("Both cohorts" = "#6C3483",
               "HPV-A7 only" = PAL_A7,
               "HPV-A9 only" = PAL_A9,
               "Neither"     = GRAY)
  ) +
  scale_size_continuous(range = c(1.5, 6.5), guide = "none") +
  coord_cartesian(xlim = c(-5.5, 1.5), ylim = c(-5.5, 1.5)) +
  labs(
    title    = "Cross-cohort drug proximity: HPV-A7 vs HPV-A9",
    subtitle = paste0(
      n_common, " drugs tested in both networks | ",
      "Shaded quadrant: significant (z < -1.96) in both HPV clades (n = ", n_both, ") | ",
      "Dashed diagonal: equal z-score in both cohorts"
    ),
    x = "Permutation z-score  (HPV-A7)",
    y = "Permutation z-score  (HPV-A9)"
  ) +
  theme_pub(base = 15) +
  theme(
    axis.text       = element_text(size = 14),
    axis.title      = element_text(size = 15, face = "bold"),
    legend.text     = element_text(size = 13),
    legend.title    = element_text(size = 14, face = "bold"),
    plot.title      = element_text(size = 18, face = "bold"),
    plot.subtitle   = element_text(size = 12.5, color = "#555555", lineheight = 1.35)
  ) +
  guides(color = guide_legend(override.aes = list(size = 5.0)))

save_fig(p5, "6_4_P5_cross_cohort", w = 14, h = 13)

# ── Summary table ─────────────────────────────────────────────

message("[Table] Saving summary...")

sig_summary <- prox_both %>%
  filter(significant) %>%
  select(cohort, drug, n_targets, d_c_interest, z_score, p_adj, sig_tier) %>%
  arrange(cohort, z_score)

write_tsv(sig_summary, "6_4_significant_drugs_summary.tsv")
message("  Saved: 6_4_significant_drugs_summary.tsv")

message("\n[Done] All plots saved.")
for (f in sort(list.files(".", pattern = "^6_4_P[0-9]_.*\\.(png|pdf)$")))
  message("  ", f)
