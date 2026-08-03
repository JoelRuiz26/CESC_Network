### ============================================================
### 6_19_All_essential_hub_leadingedge_circos.R
###
### Same design as 6_16, but WITHOUT the pan-essential filter: every
### essential hub (the "complete" set in 6_5, 395 A7 / 436 A9) is shown,
### not only the clade-selective / non pan-essential ones. Each gene is
### linked to a cancer hallmark by LEADING-EDGE membership (it sits in the
### leading edge of a hallmark that its network module enriches), one panel
### per clade.
###
###   * outer labels  = essential-hub gene names + hallmark process names
###   * coloured band = the module each gene belongs to (clade-tinted)
###   * chords        = gene -> hallmark (leading-edge), coloured by hallmark
###   * gene node col = purple if essential hub in BOTH clades, else clade col
###
### Inputs:
###   6_DepMap_cacu/results/6_5_essential_hub_gene_lists.tsv  (set == complete)
###   7_Modularity/5_2_ORA/5_2_3_Hallmark{7,9}_louv_max.tsv   (leading edge = geneID)
###   -> results/6_19_All_essential_hub_leadingedge_circos.png / .pdf
### ============================================================
suppressPackageStartupMessages({
  library(tidyverse); library(vroom); library(circlize)
})

res  <- "~/CESC_Network/6_DepMap_cacu/results"
ora  <- "~/CESC_Network/7_Modularity/5_2_ORA"

col_A7 <- "#4E79A7"; col_A9 <- "#CB6A5C"; col_both <- "#5E3C99"
FAM <- "Nimbus Sans"

## ---- stable colour + pretty label per hallmark (shared across panels) -------
hm_col <- c(
  OXIDATIVE_PHOSPHORYLATION="#2E7D6F", MYC_TARGETS_V1="#C46A34", MYC_TARGETS_V2="#E08A4A",
  DNA_REPAIR="#5B6BA8", GLYCOLYSIS="#C24E78", HYPOXIA="#4C9CC0",
  INFLAMMATORY_RESPONSE="#D6A339", MYOGENESIS="#8A6D3B", TNFA_SIGNALING_VIA_NFKB="#7C7F86",
  ANDROGEN_RESPONSE="#9E4C6E", APICAL_JUNCTION="#7BA05B", E2F_TARGETS="#DE8079",
  ESTROGEN_RESPONSE_EARLY="#B4A0CE", ESTROGEN_RESPONSE_LATE="#8C6BB1", P53_PATHWAY="#E0B36A",
  G2M_CHECKPOINT="#3E7C4F", MITOTIC_SPINDLE="#5AB4A6", MTORC1_SIGNALING="#C9556E",
  PI3K_AKT_MTOR_SIGNALING="#B23A2E", PROTEIN_SECRETION="#6BAED6",
  UNFOLDED_PROTEIN_RESPONSE="#546E7A", UV_RESPONSE_DN="#A9886B", ADIPOGENESIS="#C6A15B",
  IL2_STAT5_SIGNALING="#9C7CB8", INTERFERON_ALPHA_RESPONSE="#7FB0C9",
  INTERFERON_GAMMA_RESPONSE="#4F8FA6")

hm_pretty <- c(
  OXIDATIVE_PHOSPHORYLATION="Oxidative phosphorylation", MYC_TARGETS_V1="MYC targets V1",
  MYC_TARGETS_V2="MYC targets V2", DNA_REPAIR="DNA repair", GLYCOLYSIS="Glycolysis",
  HYPOXIA="Hypoxia", INFLAMMATORY_RESPONSE="Inflammatory response", MYOGENESIS="Myogenesis",
  TNFA_SIGNALING_VIA_NFKB="TNFα signaling (NF-κB)", ANDROGEN_RESPONSE="Androgen response",
  APICAL_JUNCTION="Apical junction", E2F_TARGETS="E2F targets",
  ESTROGEN_RESPONSE_EARLY="Estrogen resp. (early)", ESTROGEN_RESPONSE_LATE="Estrogen resp. (late)",
  P53_PATHWAY="p53 pathway", G2M_CHECKPOINT="G2M checkpoint", MITOTIC_SPINDLE="Mitotic spindle",
  MTORC1_SIGNALING="mTORC1 signaling", PI3K_AKT_MTOR_SIGNALING="PI3K–AKT–mTOR",
  PROTEIN_SECRETION="Protein secretion", UNFOLDED_PROTEIN_RESPONSE="Unfolded protein response",
  UV_RESPONSE_DN="UV response DN", ADIPOGENESIS="Adipogenesis", IL2_STAT5_SIGNALING="IL2–STAT5",
  INTERFERON_ALPHA_RESPONSE="IFN-α response", INTERFERON_GAMMA_RESPONSE="IFN-γ response")
pretty <- function(x) ifelse(x %in% names(hm_pretty), hm_pretty[x],
                             str_to_sentence(gsub("_"," ",x)))

## ---- build the leading-edge link table (ALL essential hubs) ------------------
hubs <- vroom(file.path(res, "6_5_essential_hub_gene_lists.tsv"), show_col_types = FALSE) |>
  filter(set == "complete") |>
  transmute(clade = sub("HPV-", "", group), gene, module = community_name)

## grp: Both = essential hub in the complete set of BOTH clades
grp <- hubs |> distinct(clade, gene) |> count(gene) |>
  transmute(gene, grp = if_else(n == 2, "Both", "one"))

le <- bind_rows(
  vroom(file.path(ora,"5_2_3_Hallmark7_louv_max.tsv"), show_col_types = FALSE) |>
    transmute(clade="A7", module=Community_Name, HM=sub("^HALLMARK_","",ID), geneID),
  vroom(file.path(ora,"5_2_3_Hallmark9_louv_max.tsv"), show_col_types = FALSE) |>
    transmute(clade="A9", module=Community_Name, HM=sub("^HALLMARK_","",ID), geneID)) |>
  separate_rows(geneID, sep="/") |> rename(gene = geneID)

map <- hubs |> inner_join(le, by = c("clade","module","gene")) |>
  left_join(grp, by = "gene") |> distinct()

message("[info] links A7/A9: ",
        paste(map |> count(clade) |> mutate(s=paste0(clade,"=",n)) |> pull(s), collapse="  "))

## ---- one circos panel -------------------------------------------------------
panel <- function(cl, clade_col) {
  d <- map |> filter(clade == cl)
  genes <- d |> distinct(gene, module, grp) |> arrange(module, gene)
  hmord <- d |> count(HM) |> arrange(desc(n), HM) |> pull(HM)
  gene_order <- genes$gene
  sectors <- c(gene_order, hmord)

  grid.col <- c(setNames(if_else(genes$grp == "Both", col_both, clade_col), genes$gene),
                setNames(hm_col[hmord], hmord))
  chord_col <- setNames(adjustcolor(hm_col[d$HM], 0.55), NULL)

  gaps <- setNames(rep(1.4, length(sectors)), sectors)
  last_of_mod <- genes |> group_by(module) |> slice_tail(n = 1) |> pull(gene)
  gaps[last_of_mod]         <- 3.5
  gaps[tail(gene_order, 1)] <- 12      # gene block -> hallmark block
  gaps[tail(hmord, 1)]      <- 12      # hallmark block -> gene block

  circos.clear()
  circos.par(start.degree = 90, gap.after = unname(gaps[sectors]),
             canvas.xlim = c(-1.40, 1.40), canvas.ylim = c(-1.40, 1.40),
             points.overflow.warning = FALSE, cell.padding = c(0, 0, 0, 0))

  chordDiagram(
    d[, c("gene", "HM")], order = sectors,
    grid.col = grid.col, col = chord_col, grid.border = "white",
    transparency = 0.4, link.lwd = 0.3, link.border = adjustcolor("white", 0.4),
    link.sort = TRUE, link.decreasing = TRUE,
    annotationTrack = "grid", annotationTrackHeight = mm_h(2.6),
    preAllocateTracks = list(list(track.height = mm_h(19)),   # 1: outer labels
                             list(track.height = mm_h(3.2)))) # 2: module band

  circos.track(track.index = 1, bg.border = NA, panel.fun = function(x, y) {
    si  <- CELL_META$sector.index
    isH <- si %in% hmord
    circos.text(CELL_META$xcenter, CELL_META$ylim[1] + mm_y(0.5),
                if (isH) pretty(si) else si,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),
                cex = if (isH) 1.02 else 0.62, font = 2,
                family = FAM, col = "grey12")
  })
  band <- adjustcolor(clade_col, 0.9)
  for (m in unique(genes$module)) {
    gs <- genes$gene[genes$module == m]
    highlight.sector(gs, track.index = 2, col = band, border = "white",
                     text = m, facing = "bending.inside", niceFacing = TRUE,
                     cex = 0.9, font = 2, family = FAM, text.col = "white")
  }
  highlight.sector(hmord, track.index = 2, col = "#EAEAEE", border = "white",
                   text = "Hallmark processes", facing = "bending.inside",
                   niceFacing = TRUE, cex = 0.98, font = 2, family = FAM,
                   text.col = "grey35")

  title(sprintf("HPV-%s essential hubs (all)", cl), family = FAM,
        cex.main = 2.4, font.main = 2, col.main = clade_col, line = -1.2)
  circos.clear()
}

## ---- combined figure: A7 | A9  +  node legend -------------------------------
draw_all <- function() {
  par(family = FAM)
  layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(7.2, 0.7))
  par(mar = c(0.2, 0.2, 1.8, 0.2)); panel("A7", col_A7)
  par(mar = c(0.2, 0.2, 1.8, 0.2)); panel("A9", col_A9)

  par(mar = c(0, 1, 0.2, 1)); plot.new()
  legend("top", title = "Essential-hub node", horiz = TRUE,
         legend = c("Both clades", "HPV-A7 specific", "HPV-A9 specific"),
         pt.bg = c(col_both, col_A7, col_A9), pch = 21, pt.cex = 2.6, col = "white",
         bty = "n", cex = 1.35, title.font = 2, text.width = 0.14)
}

out_png <- file.path(res, "6_19_All_essential_hub_leadingedge_circos.png")
out_pdf <- file.path(res, "6_19_All_essential_hub_leadingedge_circos.pdf")
png(out_png, width = 20, height = 10.6, units = "in", res = 320, bg = "white", type = "cairo"); draw_all(); dev.off()
cairo_pdf(out_pdf, width = 20, height = 10.6, bg = "white"); draw_all(); dev.off()
message("[OK] saved ", out_png)
