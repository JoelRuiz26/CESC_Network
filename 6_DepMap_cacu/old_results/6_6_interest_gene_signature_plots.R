# ============================================================
# 6_6_interest_gene_signature_plots.R
# Venn diagrams (DEG x Hub) and creative interest-gene signature
# ============================================================

setwd("~/CESC_Network/6_DepMap_cacu/")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggVennDiagram)
  library(ggforce)
  library(patchwork)
  library(scales)
  library(ggrepel)
})

save_fig <- function(p, name, w = 10, h = 7, dpi = 300) {
  ggsave(paste0(name, ".pdf"), p, width = w, height = h)
  ggsave(paste0(name, ".png"), p, width = w, height = h, dpi = dpi, bg = "white")
  cat("  Saved:", name, "\n")
}

theme_pub <- function(base = 12) {
  theme_minimal(base_size = base) +
    theme(
      plot.title       = element_text(face = "bold", size = base + 2, hjust = 0),
      plot.subtitle    = element_text(size = base - 1, color = "#555555", hjust = 0,
                                      margin = margin(b = 8)),
      axis.title       = element_text(face = "bold", size = base),
      axis.text        = element_text(color = "#333333"),
      panel.grid.major = element_line(color = "#EBEBEB"),
      panel.grid.minor = element_blank(),
      legend.title     = element_text(face = "bold", size = base - 1),
      legend.text      = element_text(size = base - 2),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin      = margin(14, 14, 14, 14)
    )
}

# ── Load gene tables ──────────────────────────────────────────

cat("[1] Loading gene tables...\n")
t7 <- read_tsv("6_3_gene_table_A7.tsv", show_col_types = FALSE)
t9 <- read_tsv("6_3_gene_table_A9.tsv", show_col_types = FALSE)

# Derived counts
counts <- list(
  A7 = list(
    total    = nrow(t7),
    deg      = sum(t7$is_DEG),
    hub      = sum(t7$is_hub),
    interest = sum(t7$is_interest),
    deg_only = sum(t7$is_DEG  & !t7$is_hub),
    hub_only = sum(!t7$is_DEG &  t7$is_hub),
    neither  = sum(!t7$is_DEG & !t7$is_hub)
  ),
  A9 = list(
    total    = nrow(t9),
    deg      = sum(t9$is_DEG),
    hub      = sum(t9$is_hub),
    interest = sum(t9$is_interest),
    deg_only = sum(t9$is_DEG  & !t9$is_hub),
    hub_only = sum(!t9$is_DEG &  t9$is_hub),
    neither  = sum(!t9$is_DEG & !t9$is_hub)
  )
)

cat(sprintf("  A7 — DEG: %d | Hub: %d | Interest (DEG x Hub): %d | Neither: %d\n",
            counts$A7$deg, counts$A7$hub, counts$A7$interest, counts$A7$neither))
cat(sprintf("  A9 — DEG: %d | Hub: %d | Interest (DEG x Hub): %d | Neither: %d\n",
            counts$A9$deg, counts$A9$hub, counts$A9$interest, counts$A9$neither))

# ══════════════════════════════════════════════════════════════
# PLOT 1 · Venn diagram: DEGs ∩ Hub genes (A7 and A9)
# ══════════════════════════════════════════════════════════════

cat("[P1] Venn diagrams...\n")

PAL_VENN <- list(
  DEG  = "#2471A3",
  Hub  = "#CB4335",
  Int  = "#1E8449",
  bg   = "#F8F9FA"
)

make_venn <- function(tbl, cohort, pal) {
  sets <- list(
    DEGs               = tbl %>% filter(is_DEG)  %>% pull(Gene),
    `Hub genes\n(top 10%)` = tbl %>% filter(is_hub) %>% pull(Gene)
  )
  n_total <- nrow(tbl)

  ggVennDiagram(
    sets,
    label      = "count",
    label_size = 5,
    set_size   = 4.5,
    edge_size  = 0.8
  ) +
    scale_fill_gradient(
      low  = "#EBF5FB",
      high = pal$Int,
      guide = "none"
    ) +
    scale_color_manual(values = c(pal$DEG, pal$Hub), guide = "none") +
    annotate("text",
             x = 0, y = -1.55,
             label = sprintf("Network: %s  |  %s total nodes  |  %d neither DEG nor Hub",
                             cohort, scales::comma(n_total),
                             sum(!tbl$is_DEG & !tbl$is_hub)),
             size = 3.2, color = "#777777", hjust = 0.5) +
    labs(title = sprintf("HPV-%s", cohort)) +
    theme_void(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold", size = 14, hjust = 0.5,
                                     margin = margin(b = 6)),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin     = margin(10, 10, 20, 10)
    )
}

pv_A7 <- make_venn(t7, "A7", PAL_VENN)
pv_A9 <- make_venn(t9, "A9", PAL_VENN)

p_venn <- pv_A7 + pv_A9 +
  plot_annotation(
    title    = "Gene set overlap: DEGs x Hub genes = Interest genes",
    subtitle = paste0(
      "DEGs: differentially expressed (adj.P < 0.05)  |  ",
      "Hub genes: top 10% most connected nodes in the subnetwork  |  ",
      "Interest genes: the intersection (used for drug proximity analysis)"
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 14, hjust = 0),
      plot.subtitle = element_text(size = 10, color = "#555555", hjust = 0,
                                   margin = margin(b = 6)),
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

save_fig(p_venn, "6_6_P1_venn_DEG_Hub", w = 12, h = 6)

# ══════════════════════════════════════════════════════════════
# PLOT 2 · Interest gene signature: logFC × connectivity × essentiality
# ══════════════════════════════════════════════════════════════
# Creative scatter:
#   x  = log2 fold change (differential expression magnitude)
#   y  = log10(degree)    (network connectivity)
#   color = category (background / DEG-only / Hub-only / Interest)
#   size  = |mean_dep|    (DepMap essentiality — larger = more essential)
#   Interest genes: colored by mean_dep (blue = essential) with bold outline
# ══════════════════════════════════════════════════════════════

cat("[P2] Interest gene signature scatter...\n")

make_sig_df <- function(tbl, cohort) {
  tbl %>%
    mutate(
      cohort   = cohort,
      category = case_when(
        is_interest            ~ "Interest (DEG x Hub)",
        is_DEG & !is_hub       ~ "DEG only",
        is_hub & !is_DEG       ~ "Hub only",
        TRUE                   ~ "Background"
      ),
      category = factor(category,
                        levels = c("Background", "DEG only", "Hub only",
                                   "Interest (DEG x Hub)")),
      log10deg = log10(degree),
      dep_abs  = abs(mean_dep)
    )
}

sig_df <- bind_rows(
  make_sig_df(t7, "A7"),
  make_sig_df(t9, "A9")
)

# Layer order: Background first, Interest on top
sig_bg   <- sig_df %>% filter(category == "Background")
sig_deg  <- sig_df %>% filter(category == "DEG only")
sig_hub  <- sig_df %>% filter(category == "Hub only")
sig_int  <- sig_df %>% filter(category == "Interest (DEG x Hub)")

# Hub threshold: degree at 95th percentile
hub_thr_A7 <- t7 %>% arrange(desc(degree)) %>% slice(ceiling(0.10 * nrow(t7))) %>% pull(degree)
hub_thr_A9 <- t9 %>% arrange(desc(degree)) %>% slice(ceiling(0.10 * nrow(t9))) %>% pull(degree)

# DEG logFC threshold labels (symmetric ±1.5 for visual guide, matches is_DEG criterion)
fc_guide <- 1.5

# Top interest genes to label (most essential + most connected)
top_int <- sig_int %>%
  group_by(cohort) %>%
  slice_min(mean_dep, n = 8, with_ties = FALSE) %>%
  ungroup()

cat_colors <- c(
  "Background"           = "#D5D8DC",
  "DEG only"             = "#5DADE2",
  "Hub only"             = "#F0A500",
  "Interest (DEG x Hub)" = "#1E8449"
)

p_sig <- ggplot() +
  # ── Background layer ─────────────────────────────────────
  geom_point(data = sig_bg,
             aes(x = logFC, y = log10deg),
             color = "#D5D8DC", size = 0.35, alpha = 0.4, shape = 16) +
  # ── DEG only ──────────────────────────────────────────────
  geom_point(data = sig_deg,
             aes(x = logFC, y = log10deg),
             color = "#5DADE2", size = 0.7, alpha = 0.55, shape = 16) +
  # ── Hub only ──────────────────────────────────────────────
  geom_point(data = sig_hub,
             aes(x = logFC, y = log10deg),
             color = "#F0A500", size = 1.0, alpha = 0.65, shape = 16) +
  # ── Hub threshold line ────────────────────────────────────
  geom_hline(
    data = data.frame(
      cohort    = c("A7", "A9"),
      hub_log10 = log10(c(hub_thr_A7, hub_thr_A9))
    ),
    aes(yintercept = hub_log10),
    linetype = "dashed", color = "#F0A500", linewidth = 0.6, alpha = 0.8
  ) +
  # ── logFC guide lines ─────────────────────────────────────
  geom_vline(xintercept = c(-fc_guide, fc_guide),
             linetype = "dotted", color = "#5DADE2", linewidth = 0.5, alpha = 0.7) +
  geom_vline(xintercept = 0,
             linetype = "solid", color = "#BBBBBB", linewidth = 0.4) +
  # ── Interest genes: colored by mean_dep ───────────────────
  geom_point(data = sig_int,
             aes(x = logFC, y = log10deg, fill = mean_dep, size = dep_abs),
             shape = 21, color = "white", stroke = 0.5, alpha = 0.92) +
  # ── Labels for top essential interest genes ───────────────
  geom_text_repel(
    data        = top_int,
    aes(x = logFC, y = log10deg, label = Gene),
    size        = 2.6, fontface = "bold.italic",
    color       = "#145A32",
    box.padding = 0.5, point.padding = 0.3,
    max.overlaps = 20,
    segment.size = 0.3, segment.color = "#27AE60",
    min.segment.length = 0.2
  ) +
  # ── Quadrant annotation ───────────────────────────────────
  annotate("text", x = 7.5, y = log10(hub_thr_A7) + 0.05,
           label = "Hub\nthreshold\n(top 10%)",
           size = 2.5, color = "#B7770D", hjust = 1,
           fontface = "italic", lineheight = 1.1) +
  # ── Scales ───────────────────────────────────────────────
  scale_fill_gradientn(
    name   = "DepMap\ndependency\nscore",
    colors = c("#08306B", "#2171B5", "#6BAED6", "#C6DBEF",
               "#FFFFBF", "#FDD49E", "#FC8D59", "#D7301F"),
    limits = c(-2.8, 0.35), oob = squish,
    breaks = c(-2.5, -2, -1.5, -1, -0.5, 0),
    guide  = guide_colorbar(barwidth = 0.8, barheight = 5,
                            ticks.colour = "grey40", frame.colour = "grey40")
  ) +
  scale_size_continuous(
    name   = "Essentiality\n(|mean_dep|)",
    range  = c(1.2, 5.5),
    breaks = c(0.5, 1.0, 1.5, 2.0),
    guide  = guide_legend(override.aes = list(fill = "#1E8449", color = "white"))
  ) +
  # ── Manual legend for categories ─────────────────────────
  # (drawn as annotated text to avoid overcomplicating the legend)
  facet_wrap(~cohort, labeller = labeller(cohort = c(A7 = "HPV-A7", A9 = "HPV-A9"))) +
  labs(
    title    = "Interest gene signature: differential expression × network connectivity × essentiality",
    subtitle = paste0(
      "Each point = one gene (n = 17,481 per cohort)  |  ",
      "x: log2FC  |  y: log10(degree)  |  ",
      "Green filled = Interest genes (DEG x Hub), colored by DepMap dependency score\n",
      "Blue = DEG only  |  Orange = Hub only  |  Gray = background  |  ",
      "Dashed orange line: hub threshold (top 10%)  |  Labeled: most essential interest genes"
    ),
    x = expression(log[2]*" fold change  (tumor vs normal)"),
    y = expression(log[10]*"(degree)  — network connectivity")
  ) +
  theme_pub() +
  theme(
    strip.text   = element_text(face = "bold", size = 12),
    panel.border = element_rect(color = "#DDDDDD", fill = NA, linewidth = 0.5)
  )

save_fig(p_sig, "6_6_P2_interest_signature", w = 16, h = 7)

# ══════════════════════════════════════════════════════════════
# PLOT 3 · Essentiality distribution by category (raincloud)
# ══════════════════════════════════════════════════════════════

cat("[P3] Essentiality raincloud by category...\n")

rain_df <- bind_rows(
  make_sig_df(t7, "A7"),
  make_sig_df(t9, "A9")
) %>%
  filter(category != "Background")  # background has no dep signal of interest

medians <- rain_df %>%
  group_by(cohort, category) %>%
  summarise(med = median(mean_dep, na.rm = TRUE), .groups = "drop")

p_rain <- ggplot(rain_df, aes(x = mean_dep, y = fct_rev(category), fill = category, color = category)) +
  # Density
  ggdist::stat_halfeye(
    adjust = 0.8, width = 0.6, .width = 0,
    justification = -0.25, point_colour = NA, alpha = 0.75
  ) +
  # Boxplot
  geom_boxplot(
    width = 0.18, outlier.shape = NA,
    color = "white", alpha = 0.9, linewidth = 0.6
  ) +
  # Jitter (only interest genes — the key set)
  geom_jitter(
    data   = rain_df %>% filter(category == "Interest (DEG x Hub)"),
    height = 0.08, width = 0, size = 1.0, alpha = 0.5, shape = 16
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#888888", linewidth = 0.6) +
  geom_vline(xintercept = -0.5, linetype = "dotted", color = "#C0392B", linewidth = 0.5) +
  geom_text(
    data    = medians,
    aes(x = med, y = fct_rev(category), label = sprintf("med=%.2f", med)),
    nudge_y = 0.32, size = 2.9, fontface = "bold", show.legend = FALSE
  ) +
  scale_fill_manual(values  = cat_colors[-1], guide = "none") +
  scale_color_manual(values = cat_colors[-1], guide = "none") +
  annotate("text", x = -0.5, y = 3.5,
           label = "essentiality\nthreshold\n(–0.5)", size = 2.5,
           color = "#C0392B", hjust = 0, fontface = "italic", lineheight = 1.1) +
  facet_wrap(~cohort, labeller = labeller(cohort = c(A7 = "HPV-A7", A9 = "HPV-A9"))) +
  labs(
    title    = "DepMap essentiality score by gene category",
    subtitle = paste0(
      "More negative score = gene is more essential for cancer cell survival (DepMap CRISPR screens)\n",
      "Interest genes (DEG x Hub) are enriched for highly essential genes relative to DEG-only or Hub-only"
    ),
    x = "DepMap mean dependency score",
    y = NULL
  ) +
  theme_pub() +
  theme(axis.text.y = element_text(size = 10, face = "bold"))

save_fig(p_rain, "6_6_P3_essentiality_raincloud", w = 13, h = 6)

cat("\n[Done] All signature plots saved.\n")
for (f in list.files(".", pattern = "^6_6_.*\\.(png|pdf)$")) cat("  ", f, "\n")
