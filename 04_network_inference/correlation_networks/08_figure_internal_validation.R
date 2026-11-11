#!/usr/bin/env Rscript
############################################################
# 08_figure_internal_validation.R
#
# Publication figure: internal validation of the 200 Spearman leave-x-out
# networks. Four metrics, three pair groups, violin + box, with the two
# contrasts the design actually tests drawn as brackets.
#
# Why two brackets and not one statistic. The earlier version reported a
# single QAP value per metric, inter minus the pooled intra, which invites the
# reading that A7:A7 is being compared with A9:A9. That comparison is not
# admissible: two A7 subsamples share 55-58 of 60 patients (A7 has only 65),
# two A9 subsamples share 8-29. The admissible contrasts are each intra group
# against inter, and those are what the brackets show.
#
# The dagger on the A7 bracket is not decoration. Inter-clade pairs share zero
# patients by construction, so every intra-vs-inter contrast starts with an
# overlap advantage. Within A9 the shared count varies enough (8-29) to
# regress it out: extrapolated to zero shared patients intra-A9 still beats
# inter on all four metrics (45_overlap_control.R). Within A7 it does not vary
# enough to extrapolate, so that contrast cannot be cleared the same way and
# the figure says so.
#
# Statistics: 44_qap_contrasts.R. Permutation of the clade label of each of
# the 200 networks, B = 5000, 100/100 split preserved. The permutation unit is
# the network, not the pair.
#
# Design: palette #2E6DA4 / #CC0000 / #8C8C8C -- validated for the lightness
# band, CVD separation (worst adjacent deuteranopia dE 14.5, tritan 24.4),
# normal-vision dE 24.4 and contrast. The third slot is a deliberate neutral:
# A7:A9 is the reference group, and identity is carried by the x axis, never
# by colour alone. No legend, because the groups ARE the x categories. The
# direction of each metric rides in the panel title as an arrow instead of a
# subtitle line repeating "higher = more similar".

options(bitmapType = "cairo")
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

REPO <- path.expand("~/CESC_Network")
BASE <- file.path(REPO, "04_network_inference", "correlation_networks")
OUT  <- file.path(BASE, "figure"); dir.create(OUT, showWarnings = FALSE)
METHOD <- "spearman"

PD <- fread(file.path(BASE, "pairwise_metrics", "pairwise_metrics_long.csv.gz"))[method_key == METHOD]
CT <- fread(file.path(OUT, "qap_contrasts.csv"))[method == METHOD]

lev <- c("intra-A7", "intra-A9", "inter A7-A9")
PD[, g2 := factor(group, levels = lev, labels = c("A7 : A7", "A9 : A9", "A7 : A9"))]
col <- c("A7 : A7" = "#2E6DA4", "A9 : A9" = "#CC0000", "A7 : A9" = "#8C8C8C")

metrics <- data.table(
  mkey = c("Jaccard", "Spearman", "ARI", "KS"),
  lab  = c("Edge Jaccard index  ↑",
           "Degree-sequence correlation ρ  ↑",
           "Module agreement, ARI  ↑",
           "Degree-distribution distance, KS D  ↓"),
  tag  = LETTERS[1:4])

fmt <- function(p, d, dag = FALSE)
  sprintf("paste(italic(P) %s, '  ', italic(d) == '%s'%s)",
          if (p <= 2.1e-4) "< 10^-3" else sprintf("== %.3g", p),
          sub("-", "−", sprintf("%.2f", d)),
          if (dag) ", '  †'" else "")

thm <- theme_classic(base_family = "Times New Roman", base_size = 9) +
  theme(
    plot.title   = element_text(size = 9.5, face = "bold", hjust = 0, margin = margin(b = 4)),
    plot.tag     = element_text(size = 11, face = "bold", family = "Times New Roman"),
    plot.tag.position = c(0.005, 0.985),
    axis.title   = element_blank(),
    axis.text    = element_text(size = 8.5, colour = "black"),
    axis.line    = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks   = element_line(linewidth = 0.3, colour = "black"),
    panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey90"),
    legend.position = "none",
    plot.margin  = margin(6, 7, 2, 5))

panels <- lapply(seq_len(nrow(metrics)), function(i) {
  m <- metrics[i]
  d <- PD[, .(g2, val = get(m$mkey))][!is.na(val)]
  c7 <- CT[metric == m$mkey & grepl("^A7", contrast)]
  c9 <- CT[metric == m$mkey & grepl("^A9", contrast)]
  lo <- min(d$val); hi <- max(d$val); R <- hi - lo
  y9 <- hi + 0.07 * R; y7 <- hi + 0.20 * R; tick <- 0.025 * R
  ## ticks stay inside the data range: the brackets must not stretch the axis
  ## into territory the metric never reaches (Jaccard above 1, say)
  brk <- pretty(c(lo, hi), n = 4); brk <- brk[brk >= lo - 0.02 * R & brk <= hi + 0.02 * R]

  ggplot(d, aes(g2, val, fill = g2)) +
    geom_violin(trim = TRUE, alpha = 0.5, linewidth = 0.25, colour = "grey25", width = 0.85) +
    geom_boxplot(width = 0.13, outlier.shape = NA, linewidth = 0.3,
                 fill = "white", colour = "grey15") +
    ## bracket A9:A9 vs A7:A9
    annotate("segment", x = 2, xend = 3, y = y9, yend = y9, linewidth = 0.3) +
    annotate("segment", x = 2, xend = 2, y = y9, yend = y9 - tick, linewidth = 0.3) +
    annotate("segment", x = 3, xend = 3, y = y9, yend = y9 - tick, linewidth = 0.3) +
    annotate("text", x = 2.5, y = y9 + 0.022 * R, size = 2.65, parse = TRUE,
             family = "Times New Roman", colour = "grey15", vjust = 0,
             label = fmt(c9$QAP_p, c9$cohen_d)) +
    ## bracket A7:A7 vs A7:A9
    annotate("segment", x = 1, xend = 3, y = y7, yend = y7, linewidth = 0.3) +
    annotate("segment", x = 1, xend = 1, y = y7, yend = y7 - tick, linewidth = 0.3) +
    annotate("segment", x = 3, xend = 3, y = y7, yend = y7 - tick, linewidth = 0.3) +
    annotate("text", x = 2, y = y7 + 0.022 * R, size = 2.65, parse = TRUE,
             family = "Times New Roman", colour = "grey15", vjust = 0,
             label = fmt(c7$QAP_p, c7$cohen_d, dag = TRUE)) +
    scale_fill_manual(values = col) +
    scale_y_continuous(breaks = brk, expand = expansion(mult = c(0.05, 0.30))) +
    labs(title = m$lab, tag = m$tag) + thm
})

fig <- wrap_plots(panels, nrow = 2, ncol = 2)
ggsave(file.path(OUT, "Figure_Spearman_internal_validation.pdf"), fig,
       width = 7.2, height = 6.0, device = cairo_pdf)
ggsave(file.path(OUT, "Figure_Spearman_internal_validation.png"), fig,
       width = 7.2, height = 6.0, dpi = 600, bg = "white")
cat("listo ->", OUT, "\n")
