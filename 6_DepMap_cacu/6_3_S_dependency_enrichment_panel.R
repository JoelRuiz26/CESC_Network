### ============================================================
### 6_3_S_dependency_enrichment_panel.R
###
### Assembles a single supplementary figure from two existing 6_3 plots,
### stacked as two rows and tagged (A) and (B):
###   A  hub vs non-hub dependency, empirical CDF   (6_3_A_hub_vs_nonhub_cdf)
###   B  permutation null of the mean dependency    (6_3_C_mean_dependency_null)
### Each source plot already carries its own HPV-A7 / HPV-A9 sub-titles, so only
### the outer A / B tags are added here. Row heights are set from the real
### aspect ratio of each image so the rasterised plots are not distorted.
###
###   6_3_S_dependency_enrichment_panel.png/.pdf
### Nothing is recomputed; the panels are exactly the ones 6_3 produced.
### ============================================================
suppressPackageStartupMessages({
  library(png)
  library(grid)
  library(patchwork)
  library(ggplot2)
})

res_dir <- "~/CESC_Network/6_DepMap_cacu/results"
FF <- "Times New Roman"

files <- c(A = "6_3_A_hub_vs_nonhub_cdf.png",
           B = "6_3_C_mean_dependency_permutation_null.png")

## read each PNG, keep its native aspect (width / height)
imgs   <- lapply(file.path(res_dir, files), readPNG)
aspect <- vapply(imgs, function(r) dim(r)[2] / dim(r)[1], numeric(1))

grobs <- lapply(imgs, function(r) wrap_elements(full = rasterGrob(r, interpolate = TRUE)))

fig <- grobs[[1]] / grobs[[2]] +
  ## row height proportional to 1/aspect keeps both images undistorted at equal width
  plot_layout(heights = 1 / aspect) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = FF, size = 30, face = "bold"))

W <- 14                       # canvas width (in)
H <- sum(W / aspect)          # total height that renders both rows at native aspect

ggsave(file.path(res_dir, "6_3_S_dependency_enrichment_panel.png"), fig,
       width = W, height = H, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(file.path(res_dir, "6_3_S_dependency_enrichment_panel.pdf"), fig,
       width = W, height = H, device = cairo_pdf, bg = "white", limitsize = FALSE)
message("✅ saved 6_3_S_dependency_enrichment_panel.png / .pdf")
