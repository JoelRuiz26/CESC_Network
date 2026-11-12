#!/usr/bin/env Rscript
############################################################
# 43_figure_spearman_paper.R
#
# Publication figure: internal validation of the Spearman leave-x-out
# networks. Four metrics, three pair groups, violin + box, with the two
# contrasts the design actually tests drawn as brackets.
#
# THE STATISTIC. Not the QAP of 44_qap_contrasts.R. That test permuted the
# clade label of the 200 finished networks, which assumes the networks are
# exchangeable under H0; they are not, because A7 has only 65 patients and two
# A7 subsamples of 60 therefore share 55-58 of them, while cross-clade pairs
# share none. Permuting network labels breaks that overlap structure, so its
# null was far too flat and every P value hit the 1/(B+1) floor.
#
# The test behind these brackets permutes the clade label of the PATIENTS
# (46_patient_label_null.R): the 267 tumours are split at random into a
# pseudo-A7 of 65 and a pseudo-A9 of 202, the leave-x-out subsampling and the
# networks are rebuilt, and the four metrics are recomputed. B = 100, so the
# floor is P = 1/101 = 0.0099. Overlap structure is preserved exactly; only
# clade identity is destroyed.
#
# ONE FIGURE, ONE QUESTION. This figure answers the descriptive one: how alike
# are networks within a clade and between clades. The inferential one -- is
# that gap bigger than chance -- has its own figure, 48_figure_null_distribution.R,
# which draws the permutation null with the observed statistic on top. An
# earlier version drew the null as a band inside these violins; two questions
# in one panel, and it read as a second boxplot.
#
# Design: palette #2E6DA4 / #CC0000 / #8C8C8C -- validated for the lightness
# band, CVD separation (worst adjacent deuteranopia dE 14.5, tritan 24.4),
# normal-vision dE 24.4 and contrast. The third slot is a deliberate neutral:
# A7:A9 is the reference group, and identity is carried by the x axis, never
# by colour alone. No legend, because the groups ARE the x categories. The
# box takes the hue of its violin with a white median bar.
#
# Note on n. The violins show every pair among the 100 + 100 networks. The
# test is run on a 20 + 20 design, because each permutation has to rebuild its
# networks from scratch and 100 permutations at 100 + 100 would be weeks. The
# observed statistic is computed on the same 20 + 20 so that it is compared
# with a null of matching variance; the group means barely move (Jaccard
# intra-A7 0.6953 at 20 networks, 0.6987 at 100).

options(bitmapType = "cairo")
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

BASE <- "/STORAGE/csbig/jruiz/cesc_network/ARACNE_AP/Leave_x_out"
OUT  <- file.path(BASE, "43_figure_paper"); dir.create(OUT, showWarnings = FALSE)
METHOD <- "spearman"

PD <- fread(file.path(BASE, "29_full_panel", "pairwise_long_allmethods.csv"))[method_key == METHOD]
PV <- fread(file.path(BASE, "46_patient_null", "patient_null_pvalues.csv"))

lev <- c("intra-A7", "intra-A9", "inter A7-A9")
PD[, g2 := factor(group, levels = lev, labels = c("A7 : A7", "A9 : A9", "A7 : A9"))]
col <- c("A7 : A7" = "#2E6DA4", "A9 : A9" = "#CC0000", "A7 : A9" = "#8C8C8C")

metrics <- data.table(
  mkey = c("Jaccard", "Spearman", "ARI", "KS"),
  lab  = c("Edge Jaccard index",
           "Degree-sequence correlation ρ",
           "Module agreement, ARI",
           "Degree-distribution distance, KS D"),
  tag  = LETTERS[1:4])

## P as the test returns it. At the 1/(B+1) floor the exact value is printed
## and the caption says no permutation reached the observation; above it, two
## significant digits. No "< 10^-3": that hides the resolution of the test.
fmt <- function(p) sprintf("paste(italic(P) == %s)", formatC(p, format = "fg", digits = 2))

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
  pv <- PV[metric == m$mkey]
  lo <- min(d$val); hi <- max(d$val); R <- hi - lo
  y9 <- hi + 0.05 * R; y7 <- hi + 0.16 * R; tick <- 0.022 * R
  brk <- pretty(c(lo, hi), n = 4); brk <- brk[brk >= lo - 0.02 * R & brk <= hi + 0.02 * R]

  ggplot(d, aes(g2, val, fill = g2)) +
    geom_violin(trim = TRUE, alpha = 0.5, linewidth = 0.25, colour = "grey25", width = 0.85) +
    geom_boxplot(width = 0.14, outlier.shape = NA, linewidth = 0.3, colour = "grey15") +
    stat_summary(fun = median, fun.min = median, fun.max = median,
                 geom = "crossbar", width = 0.14, linewidth = 0.3,
                 colour = "white", fatten = 0) +
    ## bracket A9:A9 vs A7:A9
    annotate("segment", x = 2, xend = 3, y = y9, yend = y9, linewidth = 0.3) +
    annotate("segment", x = 2, xend = 2, y = y9, yend = y9 - tick, linewidth = 0.3) +
    annotate("segment", x = 3, xend = 3, y = y9, yend = y9 - tick, linewidth = 0.3) +
    annotate("text", x = 2.5, y = y9 + 0.022 * R, size = 2.65, parse = TRUE,
             family = "Times New Roman", colour = "grey15", vjust = 0,
             label = fmt(pv$p_g9)) +
    ## bracket A7:A7 vs A7:A9
    annotate("segment", x = 1, xend = 3, y = y7, yend = y7, linewidth = 0.3) +
    annotate("segment", x = 1, xend = 1, y = y7, yend = y7 - tick, linewidth = 0.3) +
    annotate("segment", x = 3, xend = 3, y = y7, yend = y7 - tick, linewidth = 0.3) +
    annotate("text", x = 2, y = y7 + 0.022 * R, size = 2.65, parse = TRUE,
             family = "Times New Roman", colour = "grey15", vjust = 0,
             label = fmt(pv$p_g7)) +
    scale_fill_manual(values = col) +
    scale_y_continuous(breaks = brk, expand = expansion(mult = c(0.05, 0.24))) +
    labs(title = m$lab, tag = m$tag) + thm
})

## The band needs naming inside the figure: an unexplained graphical element
## is worse than none. One line under the grid, not a legend box.
fig <- wrap_plots(panels, nrow = 2, ncol = 2) +
  plot_annotation(
    caption = paste0(
      "Each point is a pair of networks: all pairs among the 100 + 100 ",
      "leave-x-out networks.  Brackets give the permutation P for each intra ",
      "group against inter,\nfrom the patient-label null of Figure S; ",
      "0.0099 is its floor, 0 of 100 permutations."),
    theme = theme(plot.caption = element_text(size = 7, hjust = 0, colour = "grey25",
                                              family = "Times New Roman",
                                              margin = margin(t = 6))))
ggsave(file.path(OUT, "Figure_Spearman_internal_validation.pdf"), fig,
       width = 7.2, height = 6.25, device = cairo_pdf)
ggsave(file.path(OUT, "Figure_Spearman_internal_validation.png"), fig,
       width = 7.2, height = 6.25, dpi = 600, bg = "white")
cat("listo ->", OUT, "\n")
