### ============================================================
### 6_18_Both_essential_hub_module_circos.R
###
### Circos restricted to the 11 essential hubs that are dependencies in
### BOTH clades, each linked to the cancer hallmarks enriched in ITS
### network module (leading-edge membership NOT required). Because a
### "both" gene lives in a different module per clade, the two panels
### show how the SAME shared dependency wires into different processes.
###
### Encoding = module-driven, clade-tinted (cool = A7, warm = A9);
### hallmark side neutral grey with black bold names.
### Both-genes whose module enriches no hallmark are omitted (printed).
###
### Inputs identical to 6_17.
###   -> results/6_18_Both_essential_hub_module_circos.png / .pdf
### ============================================================
suppressPackageStartupMessages({
  library(tidyverse); library(vroom); library(circlize)
})

res <- "~/CESC_Network/6_DepMap_cacu/results"
ora <- "~/CESC_Network/7_Modularity/5_2_ORA"
col_A7 <- "#4E79A7"; col_A9 <- "#CB6A5C"
FAM <- "Nimbus Sans"

pal_cool <- c("#2C5F8A","#4E9CC0","#5AB4A6","#3E7C4F","#6BAED6",
              "#8C6BB1","#5C8AA6","#7FB77E","#A9C4DE","#546E7A")
pal_warm <- c("#B23A2E","#E1743C","#E1A33C","#C9556E","#CF7FA6",
              "#9E5B33","#D98F6B","#E0B24A","#B5654A","#8A5A5A")

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

loc <- vroom(file.path(res, "6_9_B_nonpan_essential_hub_locations.tsv"), show_col_types = FALSE)
mh <- bind_rows(
  vroom(file.path(ora,"5_2_3_Hallmark7_louv_max.tsv"), show_col_types = FALSE) |>
    transmute(clade="A7", module=Community_Name, HM=sub("^HALLMARK_","",ID)),
  vroom(file.path(ora,"5_2_3_Hallmark9_louv_max.tsv"), show_col_types = FALSE) |>
    transmute(clade="A9", module=Community_Name, HM=sub("^HALLMARK_","",ID))) |> distinct()

## only the 11 shared (Both) dependencies
both <- bind_rows(
  loc |> filter(grp=="Both") |> transmute(clade="A7", gene, module=module_A7),
  loc |> filter(grp=="Both") |> transmute(clade="A9", gene, module=module_A9))

dropped <- both |> anti_join(mh, by=c("clade","module"))
message("[info] both-genes omitted (module has no enriched hallmark):")
message("       ", paste(sprintf("%s:%s", dropped$clade, dropped$gene), collapse=", "))

LINKS <- both |> inner_join(mh, by=c("clade","module"))

panel <- function(cl, clade_col, pal) {
  d <- LINKS |> filter(clade == cl)
  gcount <- d |> distinct(gene, module) |> count(module, name="ng") |> arrange(desc(ng), module)
  modord <- gcount$module
  mod_col <- setNames(pal[seq_along(modord)], modord)

  genes <- d |> distinct(gene, module) |>
    mutate(module = factor(module, levels = modord)) |> arrange(module, gene)
  gene_order <- genes$gene
  hmord <- d |> count(HM, name="n") |> arrange(desc(n), HM) |> pull(HM)
  sectors <- c(gene_order, hmord)

  g2mod <- setNames(as.character(genes$module), genes$gene)
  grid.col <- c(setNames(mod_col[g2mod[gene_order]], gene_order),
                setNames(rep("#C7C7CD", length(hmord)), hmord))
  chord_col <- adjustcolor(unname(mod_col[g2mod[d$gene]]), 0.55)

  gaps <- setNames(rep(2.4, length(sectors)), sectors)
  last_of_mod <- genes |> group_by(module) |> slice_tail(n=1) |> pull(gene)
  gaps[last_of_mod]        <- 6
  gaps[tail(gene_order,1)] <- 16
  gaps[tail(hmord,1)]      <- 16

  circos.clear()
  circos.par(start.degree = 90, gap.after = unname(gaps[sectors]),
             canvas.xlim = c(-1.34,1.34), canvas.ylim = c(-1.34,1.34),
             points.overflow.warning = FALSE, cell.padding = c(0,0,0,0))

  chordDiagram(
    d[, c("gene","HM")], order = sectors,
    grid.col = grid.col, col = chord_col, grid.border = "white",
    transparency = 0.4, link.lwd = 0.4, link.border = adjustcolor("white",0.5),
    link.sort = TRUE, link.decreasing = TRUE,
    annotationTrack = "grid", annotationTrackHeight = mm_h(3.0),
    preAllocateTracks = list(list(track.height = mm_h(22)),
                             list(track.height = mm_h(3.8))))

  circos.track(track.index = 1, bg.border = NA, panel.fun = function(x,y){
    si <- CELL_META$sector.index; isH <- si %in% hmord
    circos.text(CELL_META$xcenter, CELL_META$ylim[1] + mm_y(0.6),
                if (isH) pretty(si) else si, facing="clockwise", niceFacing=TRUE,
                adj=c(0,0.5), cex=if(isH) 1.02 else 0.98, font=2,
                family=FAM, col="grey12")
  })
  for (m in modord) {
    gs <- genes$gene[genes$module == m]
    highlight.sector(gs, track.index=2, col=mod_col[m], border="white",
                     text=m, facing="bending.inside", niceFacing=TRUE,
                     cex=0.98, font=2, family=FAM, text.col="white")
  }
  highlight.sector(hmord, track.index=2, col="#EAEAEE", border="white",
                   text="Hallmark processes", facing="bending.inside", niceFacing=TRUE,
                   cex=0.96, font=2, family=FAM, text.col="grey35")

  title(sprintf("HPV-%s  ·  shared essential hubs", cl), family=FAM,
        cex.main=2.3, font.main=2, col.main=clade_col, line=-1.0)
  circos.clear()
}

draw_all <- function() {
  par(family = FAM)
  layout(matrix(c(1,2,3,3), nrow=2, byrow=TRUE), heights=c(7.2, 0.5))
  par(mar=c(0.2,0.2,1.8,0.2)); panel("A7", col_A7, pal_cool)
  par(mar=c(0.2,0.2,1.8,0.2)); panel("A9", col_A9, pal_warm)
  par(mar=c(0,1,0,1)); plot.new()
  legend("center", horiz=TRUE,
         legend="11 hubs essential in BOTH clades  ·  linked to the hallmarks enriched in each clade's module  ·  colour = module",
         text.col="grey30", bty="n", cex=1.15, text.font=2)
}

out_png <- file.path(res, "6_18_Both_essential_hub_module_circos.png")
out_pdf <- file.path(res, "6_18_Both_essential_hub_module_circos.pdf")
png(out_png, width=18, height=9.4, units="in", res=320, bg="white", type="cairo"); draw_all(); dev.off()
cairo_pdf(out_pdf, width=18, height=9.4, bg="white"); draw_all(); dev.off()
message("[OK] saved ", out_png)
