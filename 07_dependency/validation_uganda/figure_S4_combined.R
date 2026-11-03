### ============================================================
### figure_S4_combined.R
###
### CESC + Uganda combined view of Figure S4. Does NOT touch the original
### CESC-only figure. Grouped by ANALYSIS TYPE: each row is one analysis
### (hub-vs-nonhub CDF, permutation null of the mean), CESC and Uganda side
### by side under the same tag letter. Nothing is recomputed; the panels
### are exactly the ones 03_hub_dependency_tests.R already produced per
### cohort.
###
###   A) Hub vs non-hub dependency, empirical CDF   [CESC | Uganda]
###   B) Permutation null of the mean dependency    [CESC | Uganda]
###
### -> figures/validation_uganda/Figure_S4_hub_essentiality_uganda.png/.pdf
### ============================================================
suppressPackageStartupMessages({
  library(png); library(grid); library(patchwork); library(ggplot2)
})

fig_dir <- "~/CESC_Network/figures/validation_uganda"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
FF <- "Times New Roman"

cesc_dir   <- "~/CESC_Network/07_dependency/results"
uganda_dir <- "~/CESC_Network/07_dependency/validation_uganda/results"

read_grob <- function(path) {
  img <- readPNG(path)
  list(grob = wrap_elements(full = rasterGrob(img, interpolate = TRUE)),
       aspect = dim(img)[2] / dim(img)[1])
}

cA <- read_grob(file.path(cesc_dir,   "6_3_A_hub_vs_nonhub_cdf.png"))
uA <- read_grob(file.path(uganda_dir, "6_3_A_hub_vs_nonhub_cdf.png"))
cB <- read_grob(file.path(cesc_dir,   "6_3_C_mean_dependency_permutation_null.png"))
uB <- read_grob(file.path(uganda_dir, "6_3_C_mean_dependency_permutation_null.png"))

rowA <- wrap_elements(full = (cA$grob | uA$grob))
rowB <- wrap_elements(full = (cB$grob | uB$grob))

fig <- rowA / rowB +
  plot_layout(heights = c(mean(1 / c(cA$aspect, uA$aspect)), mean(1 / c(cB$aspect, uB$aspect)))) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = FF, size = 26, face = "bold"))

W <- 20
H <- W / 2 * (mean(1 / c(cA$aspect, uA$aspect)) + mean(1 / c(cB$aspect, uB$aspect)))

ggsave(file.path(fig_dir, "Figure_S4_hub_essentiality_uganda.png"), fig, width = W, height = H, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(file.path(fig_dir, "Figure_S4_hub_essentiality_uganda.pdf"), fig, width = W, height = H, device = cairo_pdf, bg = "white", limitsize = FALSE)
message("Saved: Figure_S4_hub_essentiality_uganda.png / .pdf  (A=CDF, B=permutation null; each CESC | Uganda)")
