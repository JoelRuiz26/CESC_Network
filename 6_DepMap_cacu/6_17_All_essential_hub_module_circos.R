### ============================================================
### 6_17_All_essential_hub_module_circos.R
###
### Comprehensive circos: EVERY non pan-essential essential hub (47 A7,
### 54 A9) linked to the cancer hallmarks enriched in ITS network module
### -- a gene is linked to a process purely by sitting in a module that
### enriches it (leading-edge membership NOT required).
###
### Encoding (module-driven, clade-tinted):
###   * colour = network MODULE (cool palette for A7, warm for A9)
###   * chords = gene -> hallmark (every hallmark its module enriches)
###   * hallmark side = neutral grey blocks, black bold names
###
### Genes whose module enriches NO hallmark are omitted (nothing to link);
### they are printed to the console.
###
### Inputs:
###   6_DepMap_cacu/results/6_9_B_nonpan_essential_hub_locations.tsv
###   7_Modularity/5_2_ORA/5_2_3_Hallmark{7,9}_louv_max.tsv
###   -> results/6_17_All_essential_hub_module_circos.png / .pdf
### ============================================================
suppressPackageStartupMessages({
  library(tidyverse); library(vroom); library(circlize)
})

res <- "~/CESC_Network/6_DepMap_cacu/results"
ora <- "~/CESC_Network/7_Modularity/5_2_ORA"
col_A7 <- "#4E79A7"; col_A9 <- "#CB6A5C"; col_both <- "#5E3C99"
FAM <- "Nimbus Sans"

## stable colour per hallmark process (same colour in both panels).
## Chords are coloured by PROCESS; module bands are clade-tinted (see panel()).
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
  INTERFERON_GAMMA_RESPONSE="#4F8FA6", ALLOGRAFT_REJECTION="#B76E79",
  ANGIOGENESIS="#6FA8A0", APICAL_SURFACE="#9DB86F", CHOLESTEROL_HOMEOSTASIS="#D0946B",
  COAGULATION="#7E9CB8", COMPLEMENT="#8FA35B", EPITHELIAL_MESENCHYMAL_TRANSITION="#A05C7B",
  IL6_JAK_STAT3_SIGNALING="#B58AC0", KRAS_SIGNALING_UP="#C77E5A", NOTCH_SIGNALING="#6B8E6B",
  REACTIVE_OXYGEN_SPECIES_PATHWAY="#5FA0B0", XENOBIOTIC_METABOLISM="#9A8C5C")

hm_pretty <- c(
  ADIPOGENESIS="Adipogenesis", ALLOGRAFT_REJECTION="Allograft rejection",
  ANDROGEN_RESPONSE="Androgen response", ANGIOGENESIS="Angiogenesis",
  APICAL_JUNCTION="Apical junction", APICAL_SURFACE="Apical surface",
  CHOLESTEROL_HOMEOSTASIS="Cholesterol homeostasis", COAGULATION="Coagulation",
  COMPLEMENT="Complement", DNA_REPAIR="DNA repair", E2F_TARGETS="E2F targets",
  EPITHELIAL_MESENCHYMAL_TRANSITION="EMT",
  ESTROGEN_RESPONSE_EARLY="Estrogen resp. (early)",
  ESTROGEN_RESPONSE_LATE="Estrogen resp. (late)", G2M_CHECKPOINT="G2M checkpoint",
  GLYCOLYSIS="Glycolysis", HYPOXIA="Hypoxia", IL2_STAT5_SIGNALING="IL2–STAT5",
  IL6_JAK_STAT3_SIGNALING="IL6–JAK–STAT3", INFLAMMATORY_RESPONSE="Inflammatory response",
  INTERFERON_ALPHA_RESPONSE="IFN-α response", INTERFERON_GAMMA_RESPONSE="IFN-γ response",
  KRAS_SIGNALING_UP="KRAS signaling up", MITOTIC_SPINDLE="Mitotic spindle",
  MTORC1_SIGNALING="mTORC1 signaling", MYC_TARGETS_V1="MYC targets V1",
  MYC_TARGETS_V2="MYC targets V2", MYOGENESIS="Myogenesis",
  NOTCH_SIGNALING="Notch signaling", OXIDATIVE_PHOSPHORYLATION="Oxidative phosphorylation",
  P53_PATHWAY="p53 pathway", PI3K_AKT_MTOR_SIGNALING="PI3K–AKT–mTOR",
  PROTEIN_SECRETION="Protein secretion", REACTIVE_OXYGEN_SPECIES_PATHWAY="ROS pathway",
  TNFA_SIGNALING_VIA_NFKB="TNFα signaling (NF-κB)",
  UNFOLDED_PROTEIN_RESPONSE="Unfolded protein response", UV_RESPONSE_DN="UV response DN",
  XENOBIOTIC_METABOLISM="Xenobiotic metabolism")
pretty <- function(x) ifelse(x %in% names(hm_pretty), hm_pretty[x],
                             str_to_sentence(gsub("_"," ",x)))

## ---- data -------------------------------------------------------------------
loc <- vroom(file.path(res, "6_9_B_nonpan_essential_hub_locations.tsv"), show_col_types = FALSE)
mh <- bind_rows(
  vroom(file.path(ora,"5_2_3_Hallmark7_louv_max.tsv"), show_col_types = FALSE) |>
    transmute(clade="A7", module=Community_Name, HM=sub("^HALLMARK_","",ID)),
  vroom(file.path(ora,"5_2_3_Hallmark9_louv_max.tsv"), show_col_types = FALSE) |>
    transmute(clade="A9", module=Community_Name, HM=sub("^HALLMARK_","",ID))) |> distinct()

genes_all <- bind_rows(
  loc |> filter(grp %in% c("A7 only","Both")) |> transmute(clade="A7", gene, grp, module=module_A7),
  loc |> filter(grp %in% c("A9 only","Both")) |> transmute(clade="A9", gene, grp, module=module_A9))

## which genes drop out (module enriches nothing)
dropped <- genes_all |> anti_join(mh, by=c("clade","module"))
message("[info] omitted (module has no enriched hallmark): ",
        dropped |> count(clade) |> mutate(s=paste0(clade,"=",n)) |> pull(s) |> paste(collapse="  "))
message("       ", paste(sort(dropped$gene), collapse=", "))

## final link table: gene -> every hallmark enriched in its module
LINKS <- genes_all |> inner_join(mh, by=c("clade","module"))

## ---- one circos panel -------------------------------------------------------
panel <- function(cl, clade_col) {
  d <- LINKS |> filter(clade == cl)
  ## modules ordered by gene count desc -> palette assignment
  gcount <- d |> distinct(gene, module) |> count(module, name="ng") |> arrange(desc(ng), module)
  modord <- gcount$module

  genes <- d |> distinct(gene, module, grp) |>
    mutate(module = factor(module, levels = modord)) |> arrange(module, gene)
  gene_order <- genes$gene
  hmord <- d |> count(HM, name="n") |> arrange(desc(n), HM) |> pull(HM)
  sectors <- c(gene_order, hmord)

  ## node colour = biology status (purple = essential in BOTH clades, else clade)
  grid.col <- c(setNames(if_else(genes$grp == "Both", col_both, clade_col), gene_order),
                setNames(hm_col[hmord], hmord))
  ## chord colour = biological PROCESS (hallmark)
  chord_col <- setNames(adjustcolor(hm_col[d$HM], 0.55), NULL)

  gaps <- setNames(rep(1.6, length(sectors)), sectors)
  last_of_mod <- genes |> group_by(module) |> slice_tail(n=1) |> pull(gene)
  gaps[last_of_mod]         <- 4
  gaps[tail(gene_order,1)]  <- 11
  gaps[tail(hmord,1)]       <- 11

  circos.clear()
  circos.par(start.degree = 90, gap.after = unname(gaps[sectors]),
             canvas.xlim = c(-1.38,1.38), canvas.ylim = c(-1.38,1.38),
             points.overflow.warning = FALSE, cell.padding = c(0,0,0,0))

  chordDiagram(
    d[, c("gene","HM")], order = sectors,
    grid.col = grid.col, col = chord_col, grid.border = "white",
    transparency = 0.4, link.lwd = 0.3, link.border = adjustcolor("white",0.4),
    link.sort = TRUE, link.decreasing = TRUE,
    annotationTrack = "grid", annotationTrackHeight = mm_h(2.6),
    preAllocateTracks = list(list(track.height = mm_h(20)),
                             list(track.height = mm_h(3.6))))

  circos.track(track.index = 1, bg.border = NA, panel.fun = function(x,y){
    si <- CELL_META$sector.index; isH <- si %in% hmord
    circos.text(CELL_META$xcenter, CELL_META$ylim[1] + mm_y(0.5),
                if (isH) pretty(si) else si, facing="clockwise", niceFacing=TRUE,
                adj=c(0,0.5), cex=if(isH) 1.02 else 0.98, font=2,
                family=FAM, col="grey12")
  })
  ## module bands: clade-tinted, module name in white
  band <- adjustcolor(clade_col, 0.9)
  for (m in modord) {
    gs <- genes$gene[genes$module == m]
    highlight.sector(gs, track.index=2, col=band, border="white",
                     text=m, facing="bending.inside", niceFacing=TRUE,
                     cex=1.0, font=2, family=FAM, text.col="white")
  }
  highlight.sector(hmord, track.index=2, col="#EAEAEE", border="white",
                   text="Hallmark processes", facing="bending.inside", niceFacing=TRUE,
                   cex=0.98, font=2, family=FAM, text.col="grey35")

  title(sprintf("HPV-%s essential hubs", cl), family=FAM,
        cex.main=2.4, font.main=2, col.main=clade_col, line=-1.0)
  circos.clear()
}

## ---- combined figure --------------------------------------------------------
draw_all <- function() {
  par(family = FAM)
  layout(matrix(c(1,2,3,3), nrow=2, byrow=TRUE), heights=c(7.2, 0.6))
  par(mar=c(0.2,0.2,1.6,0.2)); panel("A7", col_A7)
  par(mar=c(0.2,0.2,1.6,0.2)); panel("A9", col_A9)
  par(mar=c(0,1,0,1)); plot.new()
  legend("top", title = "Essential-hub node", horiz=TRUE,
         legend = c("Both clades", "HPV-A7 specific", "HPV-A9 specific"),
         pt.bg = c(col_both, col_A7, col_A9), pch = 21, pt.cex = 2.6, col = "white",
         bty="n", cex=1.3, title.font=2, text.width=0.14)
  mtext("Each gene links to the hallmarks enriched in its module  ·  colour = biological process (hallmark)",
        side=1, line=-1.2, cex=0.95, col="grey35", font=2)
}

out_png <- file.path(res, "6_17_All_essential_hub_module_circos.png")
out_pdf <- file.path(res, "6_17_All_essential_hub_module_circos.pdf")
png(out_png, width=19, height=10, units="in", res=320, bg="white", type="cairo"); draw_all(); dev.off()
cairo_pdf(out_pdf, width=19, height=10, bg="white"); draw_all(); dev.off()
message("[OK] saved ", out_png)
