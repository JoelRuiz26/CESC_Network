### ============================================================
### 6_16_Essential_hub_module_circos.R
###
### Grouped chord (circos) schematic of WHERE the non pan-essential
### essential hubs sit (which network module) and WHICH cancer hallmark
### they respond to (leading-edge membership), one panel per clade.
###
###   * outer labels  = essential-hub gene names + hallmark process names
###   * coloured band = the module each gene belongs to
###   * chords        = gene -> hallmark (leading-edge), coloured by hallmark
###   * gene node col = purple if essential in BOTH clades, else clade colour
###
### Genes/mapping come from 6_9_B_nonpan_essential_hub_locations.tsv crossed
### with the module hallmark ORA leading edge (5_2_3_Hallmark{7,9}_louv_max).
###
###   -> results/6_16_Essential_hub_module_circos.png / .pdf
### ============================================================
suppressPackageStartupMessages({
  library(tidyverse); library(vroom); library(circlize)
})

res  <- "~/CESC_Network/6_DepMap_cacu/results"
ora  <- "~/CESC_Network/7_Modularity/5_2_ORA"

col_A7 <- "#4E79A7"; col_A9 <- "#CB6A5C"; col_both <- "#5E3C99"

## refined, publication-grade palette (muted / harmonious, still separable)
hm_col <- c(
  OXIDATIVE_PHOSPHORYLATION = "#2E7D6F", MYC_TARGETS_V1 = "#C46A34",
  DNA_REPAIR = "#5B6BA8", GLYCOLYSIS = "#C24E78", HYPOXIA = "#4C9CC0",
  INFLAMMATORY_RESPONSE = "#D6A339", MYOGENESIS = "#8A6D3B",
  TNFA_SIGNALING_VIA_NFKB = "#7C7F86", ANDROGEN_RESPONSE = "#9E4C6E",
  APICAL_JUNCTION = "#7BA05B", E2F_TARGETS = "#DE8079",
  ESTROGEN_RESPONSE_EARLY = "#B4A0CE", P53_PATHWAY = "#E0B36A")

## gently desaturate the process palette toward its own luminance -> a calmer,
## print-friendly "Nature" feel while keeping hues separable
desat <- function(cols, keep = 0.78) {
  m <- col2rgb(cols) / 255
  g <- 0.30 * m[1, ] + 0.59 * m[2, ] + 0.11 * m[3, ]
  rgb(keep * m[1, ] + (1 - keep) * g,
      keep * m[2, ] + (1 - keep) * g,
      keep * m[3, ] + (1 - keep) * g)
}
hm_col <- setNames(desat(hm_col), names(hm_col))

FAM <- "Nimbus Sans"   # Helvetica-equivalent (URW), the journal-standard sans
hm_pretty <- c(
  OXIDATIVE_PHOSPHORYLATION = "Oxidative phosphorylation", MYC_TARGETS_V1 = "MYC targets V1",
  DNA_REPAIR = "DNA repair", GLYCOLYSIS = "Glycolysis", HYPOXIA = "Hypoxia",
  INFLAMMATORY_RESPONSE = "Inflammatory response", MYOGENESIS = "Myogenesis",
  TNFA_SIGNALING_VIA_NFKB = "TNFα signaling (NF-κB)", ANDROGEN_RESPONSE = "Androgen response",
  APICAL_JUNCTION = "Apical junction", E2F_TARGETS = "E2F targets",
  ESTROGEN_RESPONSE_EARLY = "Estrogen response (early)", P53_PATHWAY = "p53 pathway")

map <- vroom(file.path(ora, "5_2_ess_hub_module_hallmark_map.tsv"), show_col_types = FALSE)

## ---- one circos panel -------------------------------------------------------
panel <- function(cl, clade_col, tag) {
  d <- map |> filter(clade == cl)
  genes <- d |> distinct(gene, module, grp) |> arrange(module, gene)
  hmord <- d |> count(HM) |> arrange(desc(n), HM) |> pull(HM)
  gene_order <- genes$gene
  sectors <- c(gene_order, hmord)

  grid.col <- c(setNames(if_else(genes$grp == "Both", col_both, clade_col), genes$gene),
                setNames(hm_col[hmord], hmord))
  chord_col <- setNames(adjustcolor(hm_col[d$HM], 0.55), NULL)

  ## gaps: small within a module, big at module borders and between blocks
  gaps <- setNames(rep(2.6, length(sectors)), sectors)
  last_of_mod <- genes |> group_by(module) |> slice_tail(n = 1) |> pull(gene)
  gaps[last_of_mod]        <- 6.5
  gaps[tail(gene_order, 1)] <- 20      # gene block -> hallmark block
  gaps[tail(hmord, 1)]      <- 20      # hallmark block -> gene block

  circos.clear()
  circos.par(start.degree = 90, gap.after = unname(gaps[sectors]),
             canvas.xlim = c(-1.26, 1.26), canvas.ylim = c(-1.26, 1.26),
             points.overflow.warning = FALSE, cell.padding = c(0, 0, 0, 0))

  chordDiagram(
    d[, c("gene", "HM")], order = sectors,
    grid.col = grid.col, col = chord_col, grid.border = "white",
    transparency = 0.45, link.lwd = 0.3, link.border = adjustcolor("white", 0.4),
    link.sort = TRUE, link.decreasing = TRUE,
    annotationTrack = "grid", annotationTrackHeight = mm_h(3.0),
    preAllocateTracks = list(list(track.height = mm_h(22)),   # 1: outer labels
                             list(track.height = mm_h(3.8)))) # 2: module band

  ## outer labels: gene symbols in italic (journal convention), process names
  ## in regular weight; single dark-grey ink for a calm, typographic hierarchy
  circos.track(track.index = 1, bg.border = NA, panel.fun = function(x, y) {
    si  <- CELL_META$sector.index
    isH <- si %in% hmord
    circos.text(CELL_META$xcenter, CELL_META$ylim[1] + mm_y(0.7),
                if (isH) hm_pretty[si] else si,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),
                cex = if (isH) 1.14 else 1.0,
                font = if (isH) 1 else 3,
                family = FAM, col = "grey15")
  })
  ## module bands: soft clade tint, module label inked in the clade colour
  ## (no saturated chip) — reads as a quiet grouping annotation
  band_fill <- adjustcolor(clade_col, 0.20)
  for (m in unique(genes$module)) {
    gs <- genes$gene[genes$module == m]
    highlight.sector(gs, track.index = 2, col = band_fill, border = "white",
                     text = m, facing = "bending.inside", niceFacing = TRUE,
                     cex = 0.92, font = 2, family = FAM,
                     text.col = adjustcolor(clade_col, 0.95))
  }
  highlight.sector(hmord, track.index = 2, col = "#EFEFF1", border = "white",
                   text = "Hallmark processes", facing = "bending.inside",
                   niceFacing = TRUE, cex = 0.9, font = 1, family = FAM,
                   text.col = "grey45")

  title(sprintf("HPV-%s essential hubs", cl), family = FAM,
        cex.main = 2.0, font.main = 2, col.main = clade_col, line = -1.0)
  ## main-figure panel tag (a / b), top-left
  mtext(tag, side = 3, line = -0.4, adj = 0.02, cex = 2.2, font = 2,
        col = "grey15", family = FAM)
  circos.clear()
}

## ---- combined figure: A7 | A9  +  node legend -------------------------------
draw_all <- function() {
  par(family = FAM)
  layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(6.7, 0.62))
  par(mar = c(0.2, 0.4, 1.6, 0.4)); panel("A7", col_A7, "a")
  par(mar = c(0.2, 0.4, 1.6, 0.4)); panel("A9", col_A9, "b")

  par(mar = c(0, 1, 0.2, 1)); plot.new()
  legend("top", horiz = TRUE,
         legend = c("essential in both clades", "HPV-A7 specific", "HPV-A9 specific"),
         pt.bg = c(col_both, col_A7, col_A9), pch = 21, pt.cex = 2.3, col = "white",
         bty = "n", cex = 1.15, text.col = "grey25", text.width = 0.135)
}

out_png <- file.path(res, "6_16_Essential_hub_module_circos.png")
out_pdf <- file.path(res, "6_16_Essential_hub_module_circos.pdf")
png(out_png, width = 16, height = 8.8, units = "in", res = 340, bg = "white", type = "cairo"); draw_all(); dev.off()
cairo_pdf(out_pdf, width = 16, height = 8.8, bg = "white"); draw_all(); dev.off()
message("[OK] saved ", out_png)
