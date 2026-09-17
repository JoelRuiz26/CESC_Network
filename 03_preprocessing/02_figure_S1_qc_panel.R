### ============================================================
### 2_9_QC_supplementary_panel.R
###
### Assembles Supplementary Figure S1 from the individual QC plots that
### 2_Prepro_TCGA_GTEX_RNA.R already wrote to plots_prepro_data/. It
### rasterises those PDFs and lays them out, so the panels are exactly
### the ones produced by the pipeline; nothing is recomputed and no
### pipeline output is touched.
###
### Layout: A and B on the top row, C and D on the bottom row. Inside every
### panel the LEFT column is before normalisation and the RIGHT column is
### after; in B and C the top row is clade A7 and the bottom row is A9.
### Before and after were stacked vertically in an earlier version, which made
### the whole panel nearly square; side by side halves the height of A and D
### and gives the horizontal format.
###   A  global expression distribution per clade      (1 unit)
###   B  transcript-length bias, A7 and A9             (2 units)
###   C  GC-content bias, A7 and A9                    (2 units)
###   D  M-value distributions                         (1 unit)
###
###   2_9_QC_supplementary_panel.png/.pdf
###
### The PNGs are produced beforehand with
###   pdftoppm -r 300 -png -singlefile <plot>.pdf S1_panels/<plot>
###   convert <plot>.png -trim +repage -bordercolor white -border 12 trim_<plot>.png
### The trim step removes the wide empty margins that the base-R devices leave
### around each plot, which were the main source of white space in the panel.
### The two M-value plots of panel D additionally get their title chopped off:
###   convert <plot>.png -trim +repage -gravity North -chop 0x170 ...
### that device was written at a different size, so its title rendered far
### larger than the titles of the other panels. The reference sample it named
### (TCGA-4J-AA1J-01) belongs in the figure caption instead.
### ============================================================
suppressPackageStartupMessages({
  library(png)
  library(grid)
  library(patchwork)
  library(ggplot2)
})

base_dir <- "~/CESC_Network/03_preprocessing"
fig_dir  <- "~/CESC_Network/figures"
tmp_dir  <- file.path(base_dir, "S1_panels")
FF <- "Times New Roman"

## a rasterised plot, keeping its own aspect ratio
img <- function(name) {
  r <- readPNG(file.path(tmp_dir, paste0("trim_", name, ".png")))
  wrap_elements(full = rasterGrob(r, interpolate = TRUE))
}

## a small caption strip above a row of plots
strip <- function(txt)
  wrap_elements(full = textGrob(txt, gp = gpar(fontfamily = FF, fontsize = 15,
                                               fontface = "italic", col = "grey25")))

## ── Panels ────────────────────────────────────────────────────────────────────
## left = before normalisation, right = after
panel_A <- (img("QC_boxplot_counts_before") | img("QC_boxplot_counts_after"))

panel_B <- ((img("QC_length_bias_A7_before") | img("QC_length_bias_A7_after")) /
            (img("QC_length_bias_A9_before") | img("QC_length_bias_A9_after")))

panel_C <- ((img("QC_GC_bias_A7_before") | img("QC_GC_bias_A7_after")) /
            (img("QC_GC_bias_A9_before") | img("QC_GC_bias_A9_after")))

panel_D <- (img("QC_Mvalues_mycd_before") | img("QC_Mvalues_mycd_after"))

## ── Composite ────────────────────────────────────────────────────────────────
## Column widths are a compromise. Solving for equal row heights gave A 60% and
## D 63% of their rows, which dwarfed B and C. Both are capped at half their row
## instead, so the four panels are comparable in size; A and D end up shorter
## than their slot and carry some vertical air, which is the price of it.
design <- c(area(1,  1, 1,  50),   # A
            area(1, 51, 1, 100),   # B
            area(2,  1, 2,  50),   # C
            area(2, 51, 2, 100))   # D

fig <- wrap_elements(full = panel_A) + wrap_elements(full = panel_B) +
       wrap_elements(full = panel_C) + wrap_elements(full = panel_D) +
  plot_layout(design = design, heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = FF, size = 30, face = "bold"))

ggsave(file.path(fig_dir, "Figure_S1_qc_normalisation.png"), fig,
       width = 22, height = 14, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(file.path(fig_dir, "Figure_S1_qc_normalisation.pdf"), fig,
       width = 22, height = 14, device = cairo_pdf, bg = "white", limitsize = FALSE)
message("saved 2_9_QC_supplementary_panel.png / .pdf")
