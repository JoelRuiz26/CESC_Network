### Joel Ruiz
### Assembles the Uganda-cohort counterpart of Supplementary Figure S1 from
### the individual QC plots 01_normalise_rnaseq.R already wrote to
### plots_prepro_data/ -- same layout, same method as
### 03_preprocessing/02_figure_S1_qc_panel.R for CESC (rasterise + trim the
### PDFs, nothing recomputed). See that script's header for the exact
### pdftoppm/convert commands used to produce S1_panels/trim_*.png.
suppressPackageStartupMessages({
  library(png); library(grid); library(patchwork); library(ggplot2)
})

base_dir <- "~/CESC_Network/03_preprocessing/validation_uganda"
tmp_dir  <- file.path(base_dir, "S1_panels")
fig_dir  <- "~/CESC_Network/figures/validation_uganda"   # shared figures/ folder, alongside CESC's own
FF <- "Times New Roman"

img <- function(name) {
  r <- readPNG(file.path(tmp_dir, paste0("trim_", name, ".png")))
  wrap_elements(full = rasterGrob(r, interpolate = TRUE))
}

panel_A <- (img("QC_boxplot_counts_before") | img("QC_boxplot_counts_after"))
panel_B <- ((img("QC_length_bias_A7_before") | img("QC_length_bias_A7_after")) /
            (img("QC_length_bias_A9_before") | img("QC_length_bias_A9_after")))
panel_C <- ((img("QC_GC_bias_A7_before") | img("QC_GC_bias_A7_after")) /
            (img("QC_GC_bias_A9_before") | img("QC_GC_bias_A9_after")))
panel_D <- (img("QC_Mvalues_mycd_before") | img("QC_Mvalues_mycd_after"))

design <- c(area(1,  1, 1,  50), area(1, 51, 1, 100),
            area(2,  1, 2,  50), area(2, 51, 2, 100))

fig <- wrap_elements(full = panel_A) + wrap_elements(full = panel_B) +
       wrap_elements(full = panel_C) + wrap_elements(full = panel_D) +
  plot_layout(design = design, heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = FF, size = 30, face = "bold"))

ggsave(file.path(fig_dir, "Figure_S1_qc_normalisation_uganda.png"), fig,
       width = 22, height = 14, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(file.path(fig_dir, "Figure_S1_qc_normalisation_uganda.pdf"), fig,
       width = 22, height = 14, device = cairo_pdf, bg = "white", limitsize = FALSE)
message("saved Figure_S1_qc_normalisation_uganda.png / .pdf")
