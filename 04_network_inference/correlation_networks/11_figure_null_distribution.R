#!/usr/bin/env Rscript
############################################################
# 48_figure_null_distribution.R
#
# The inferential half of the validation, as its own figure.
#
# 43_figure_spearman_paper.R answers the descriptive question -- how alike are
# networks within a clade and between clades. This one answers the other one:
# is that gap bigger than the gap random clades produce. They were in one
# panel and that was a mistake; a null drawn inside the violins reads as a
# second boxplot and makes a reader hunt for which comparison is which.
#
# Statistic, the same one the brackets of the other figure report:
#     T = mean(inter pairs) - mean(intra pairs)
# so a clade effect pushes T away from zero -- negative for the three
# similarities, positive for KS D, which is a distance.
#
# The null is 46_patient_label_null.R: the 267 tumours are split at random
# into a pseudo-A7 of 65 and a pseudo-A9 of 202, the leave-x-out subsampling
# and the networks are rebuilt, and T is recomputed. 100 permutations. The
# histogram IS that null; the vertical line is the observed T.
#
# What the figure shows at a glance: for the Jaccard the observed line sits
# inside its own null, for the degree metrics it sits clear of it.

options(bitmapType = "cairo")
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

BASE <- "/STORAGE/csbig/jruiz/cesc_network/ARACNE_AP/Leave_x_out"
OUT  <- file.path(BASE, "43_figure_paper"); dir.create(OUT, showWarnings = FALSE)

NS <- fread(file.path(BASE, "46_patient_null", "null_stats.csv"))
OB <- fread(file.path(BASE, "46_patient_null", "observed.csv"))
PV <- fread(file.path(BASE, "46_patient_null", "patient_null_pvalues.csv"))

CON <- c("A7 : A7  vs  A7 : A9", "A9 : A9  vs  A7 : A9")
col <- setNames(c("#2E6DA4", "#CC0000"), CON)

metrics <- data.table(
  mkey = c("Jaccard", "Spearman", "ARI", "KS"),
  lab  = c("Edge Jaccard index",
           "Degree-sequence correlation ρ",
           "Module agreement, ARI",
           "Degree-distribution distance, KS D"),
  tag  = LETTERS[1:4])

thm <- theme_classic(base_family = "Times New Roman", base_size = 9) +
  theme(
    plot.title   = element_text(size = 9.5, face = "bold", hjust = 0, margin = margin(b = 4)),
    plot.tag     = element_text(size = 11, face = "bold", family = "Times New Roman"),
    plot.tag.position = c(0.005, 0.985),
    axis.title.x = element_text(size = 8.5, margin = margin(t = 3)),
    axis.title.y = element_text(size = 8.5),
    axis.text    = element_text(size = 8, colour = "black"),
    axis.line    = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks   = element_line(linewidth = 0.3, colour = "black"),
    panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey90"),
    legend.position = "none",
    plot.margin  = margin(6, 8, 2, 5))

panels <- lapply(seq_len(nrow(metrics)), function(i) {
  m  <- metrics[i]
  ns <- NS[metric == m$mkey]
  ob <- OB[metric == m$mkey]
  pv <- PV[metric == m$mkey]
  d  <- rbind(data.table(con = CON[1], T = ns$T_g7),
              data.table(con = CON[2], T = ns$T_g9))
  d[, con := factor(con, levels = CON)]
  ov <- data.table(con = factor(CON, levels = CON),
                   T = c(ob$T_g7, ob$T_g9),
                   p = c(pv$p_g7, pv$p_g9))
  rng <- range(c(d$T, ov$T)); pad <- 0.12 * diff(rng)
  bw  <- diff(rng) / 26
  ## a label on a line near the right edge has to hang to the left of it, or
  ## it runs off the panel and collides with the other contrast's label
  ov[, hj := ifelse(T > mean(rng), 1.08, -0.08)]
  ## and stagger them vertically: with both lines near the middle the two
  ## labels meet head-on whatever the hjust
  ov[, vj := c(1.4, 3.1)]

  ggplot(d, aes(T, fill = con)) +
    geom_histogram(binwidth = bw, alpha = 0.55, colour = "grey30",
                   linewidth = 0.15, position = "identity") +
    geom_vline(data = ov, aes(xintercept = T, colour = con),
               linewidth = 0.55) +
    geom_text(data = ov, aes(x = T, colour = con, hjust = hj, vjust = vj,
                             label = sprintf("italic(P) == %s",
                                             formatC(p, format = "fg", digits = 2))),
              y = Inf, size = 2.5, parse = TRUE,
              family = "Times New Roman", show.legend = FALSE) +
    geom_blank() +
    scale_fill_manual(values = col) + scale_colour_manual(values = col) +
    scale_x_continuous(limits = c(rng[1] - pad, rng[2] + pad)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.34))) +
    labs(title = m$lab, tag = m$tag,
         x = expression(italic(T) == mean(inter) - mean(intra)),
         y = "permutations") + thm
})

## one legend strip, built by hand so the panels keep none
leg <- ggplot(data.table(x = 1:2, con = factor(CON, levels = CON)),
              aes(x, 0, fill = con)) +
  geom_point(shape = 22, size = 2.6, colour = "grey30", stroke = 0.3) +
  geom_text(aes(label = con), hjust = -0.18, size = 2.7,
            family = "Times New Roman", colour = "grey15") +
  scale_fill_manual(values = col) +
  scale_x_continuous(limits = c(0.75, 3.4)) +
  theme_void() + theme(legend.position = "none")

## caption as its own row, not plot_annotation: patchwork here refuses the
## theme object and the text comes out unstyled and clipped.
cap <- paste0(
  "Histograms: the null. 100 random splits of the 267 patients into a pseudo-A7 of 65\n",
  "and a pseudo-A9 of 202, rebuilding the subsamples, the networks and the metrics\n",
  "each time. Vertical lines: the observed statistic. The split keeps the patient-overlap\n",
  "structure of the real design and destroys only clade identity, so a line inside its own\n",
  "histogram means the intra\u2013inter gap is what the subsampling produces by itself.\n",
  "P = 0.0099 is the floor: 0 of 100 permutations.")
capg <- ggplot() + annotate("text", x = 0, y = 0, label = cap, hjust = 0, vjust = 1,
                            size = 2.45, family = "Times New Roman", colour = "grey25",
                            lineheight = 1.25) +
  scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(-1, 0)) +
  theme_void()

## the 2x2 grid counts as ONE row of the outer layout, so three heights
fig <- wrap_plots(panels, nrow = 2, ncol = 2) / leg / capg +
  plot_layout(heights = c(1, 0.055, 0.22))

ggsave(file.path(OUT, "Figure_Spearman_permutation_null.pdf"), fig,
       width = 7.2, height = 6.5, device = cairo_pdf)
ggsave(file.path(OUT, "Figure_Spearman_permutation_null.png"), fig,
       width = 7.2, height = 6.5, dpi = 600, bg = "white")
cat("listo ->", OUT, "\n")
