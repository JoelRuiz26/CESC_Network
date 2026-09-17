### ============================================================
### 6_19_All_essential_hub_leadingedge_circos.R
###
### Same design as 6_16, but WITHOUT the pan-essential filter: every
### essential hub (the "complete" set in 6_5, 395 A7 / 436 A9) is eligible,
### not only the clade-selective / non pan-essential ones. Each gene is
### linked to a cancer hallmark by LEADING-EDGE membership (it sits in the
### leading edge of a hallmark that its network module enriches), one panel
### per clade.
###
### Every essential hub with at least MIN_HM leading-edge links is drawn, so
### the default MIN_HM = 1 keeps all of them.
###
###   * outer labels  = essential-hub gene names + hallmark process names
###   * coloured band = the module each gene belongs to (clade-tinted)
###   * chords        = gene -> hallmark (leading-edge), coloured by hallmark
###   * gene node col = purple if essential hub in BOTH clades, else clade col
###
### Inputs:
###   07_dependency/results/6_5_essential_hub_gene_lists.tsv  (set == complete)
###   06_modules/5_2_3_Hallmark{7,9}_louv_max.tsv   (leading edge = geneID)
###   -> results/6_19_All_essential_hub_leadingedge_circos.png / .pdf
### ============================================================
suppressPackageStartupMessages({
  library(tidyverse); library(vroom); library(circlize)
})

res  <- "~/CESC_Network/07_dependency/results"
fig_dir <- "~/CESC_Network/figures"
ora  <- "~/CESC_Network/06_modules/results"

col_A7 <- "#4E79A7"; col_A9 <- "#CB6A5C"; col_both <- "#5E3C99"
FAM  <- "Nimbus Sans"
GAP  <- 1.7     # degrees between two neighbouring process sectors
GAP_G <- 2.3    # ... and between two gene sectors, which buys the gene names
                # the arc they need to be set larger than the process names
LEAD <- 1.09    # minimum white space between two neighbouring labels
## gene names are italic, as gene symbols are set in print; process names stay
## upright bold. The mark on a gene name says one thing only: the gene is NOT
## pan-essential (it is outside DepMap's common-essential set). Which clade a
## gene is essential in is already carried by the colour of its node, so the
## mark does not repeat it.
FNT_GENE <- 3; FNT_HM <- 2
MARK <- "\u25C6"    # the selectivity mark, set on the name's own baseline and
                    # always on the end of the symbol that faces the circle, so
                    # the marks read as one ring and never as part of the name
CX_MOD <- 3.75   # largest type the module names are allowed to take
MIN_HM <- 1     # a gene is drawn only if it is leading edge of >= MIN_HM hallmarks

darker <- function(x, k = 0.7) rgb(t(col2rgb(x) / 255 * k))

## pull the process palette toward its own luminance: print-friendly, still
## separable (same treatment as 6_16)
desat <- function(cols, keep = 0.78) {
  m <- col2rgb(cols) / 255
  g <- 0.30 * m[1, ] + 0.59 * m[2, ] + 0.11 * m[3, ]
  rgb(keep * m[1, ] + (1 - keep) * g,
      keep * m[2, ] + (1 - keep) * g,
      keep * m[3, ] + (1 - keep) * g)
}

## Sector labels are the MSigDB gene-set names exactly as they come out of the
## ORA output: the HALLMARK_ prefix is dropped (on read) and underscores become
## spaces. One line, nothing renamed, abbreviated, truncated or wrapped.
pretty <- function(x) gsub("_", " ", x)

## ---- build the leading-edge link table (ALL essential hubs) ------------------
hub_tbl <- vroom(file.path(res, "6_5_essential_hub_gene_lists.tsv"), show_col_types = FALSE) |>
  mutate(clade = sub("HPV-", "", group))
hubs <- hub_tbl |> filter(set == "complete") |>
  transmute(clade, gene, module = community_name) |>
  ## clade-selective = essential here and not pan-essential; it is the "selective"
  ## set of 6_5, so the face of a gene name is read off the data, never typed in
  left_join(hub_tbl |> filter(set == "selective") |>
              transmute(clade, gene, selective = TRUE), by = c("clade", "gene")) |>
  mutate(selective = replace_na(selective, FALSE))

## grp: Both = essential hub in the complete set of BOTH clades
grp <- hubs |> distinct(clade, gene) |> count(gene) |>
  transmute(gene, grp = if_else(n == 2, "Both", "one"))

le <- bind_rows(
  vroom(file.path(ora,"5_2_3_Hallmark7_louv_max.tsv"), show_col_types = FALSE) |>
    transmute(clade="A7", module=Community_Name, HM=sub("^HALLMARK_","",ID), geneID),
  vroom(file.path(ora,"5_2_3_Hallmark9_louv_max.tsv"), show_col_types = FALSE) |>
    transmute(clade="A9", module=Community_Name, HM=sub("^HALLMARK_","",ID), geneID)) |>
  separate_rows(geneID, sep="/") |> rename(gene = geneID)

map_all <- hubs |> inner_join(le, by = c("clade","module","gene")) |>
  left_join(grp, by = "gene") |> distinct()

## keep the multi-process hubs (see header): degree counted within the clade
deg <- map_all |> distinct(clade, gene, HM) |> count(clade, gene, name = "k")
map <- map_all |> inner_join(deg |> filter(k >= MIN_HM), by = c("clade","gene"))

message("[info] genes kept (k >= ", MIN_HM, ") / dropped: ",
        paste(deg |> group_by(clade) |>
                summarise(s = paste0(cur_group()$clade, "=", sum(k >= MIN_HM),
                                     "/", sum(k < MIN_HM)), .groups = "drop") |>
                pull(s), collapse = "  "))
message("[info] links A7/A9: ",
        paste(map |> count(clade) |> mutate(s=paste0(clade,"=",n)) |> pull(s), collapse="  "))

## ---- colour per hallmark, derived from the gene sets actually present -------
## Keyed on the sorted set of hallmarks in the data, so a hallmark can never be
## missing from the palette and both panels always agree on the colour.
hm_all <- sort(unique(map$HM))
hm_col <- setNames(desat(hcl.colors(length(hm_all), palette = "Dark 3")), hm_all)

## verification: every drawn label must be its gene-set name, spaces for "_"
stopifnot(identical(pretty(hm_all), gsub("_", " ", hm_all)), !any(grepl("\n", pretty(hm_all))))
## a marked gene label is its symbol, verbatim, with the mark and one space on
## one end of it -- whichever end the label's own rotation puts next to the circle
stopifnot(identical(gsub(paste0("^", MARK, " | ", MARK, "$"), "",
                         c(paste(MARK, "TP53"), paste("TP53", MARK))),
                    c("TP53", "TP53")))
message("[check] ", length(hm_all), " hallmark labels verbatim, single line")

## ---- the sectors of one panel, before anything is drawn ---------------------
RLAB <- 0.97   # radius the outer labels start at, in circle radii
LBL  <- 6      # thickness of the track the outer labels sit in, mm
EXP  <- 1.08   # the 4 % axis expansion R adds at both ends of a plot range

geom <- function(cl, clade_col) {
  d <- map |> filter(clade == cl)
  genes <- d |> distinct(gene, module, grp, selective) |> arrange(module, gene)
  hmord <- d |> count(HM) |> arrange(desc(n), HM) |> pull(HM)
  gene_order <- genes$gene
  sectors <- c(gene_order, hmord)

  gaps <- setNames(c(rep(GAP_G, length(gene_order)), rep(GAP, length(hmord))), sectors)
  last_of_mod <- genes |> group_by(module) |> slice_tail(n = 1) |> pull(gene)
  gaps[last_of_mod]         <- GAP_G * 2.5   # module boundary
  gaps[tail(gene_order, 1)] <- 8             # gene block -> hallmark block
  gaps[tail(hmord, 1)]      <- 8             # hallmark block -> gene block

  span <- c(table(d$gene)[gene_order], table(d$HM)[hmord]) *
            (360 - sum(gaps)) / (2 * nrow(d))               # sector width, deg
  ## sector edges and mid-angles: circlize starts at 90 deg and runs clockwise
  sdeg <- 90 - (cumsum(span + gaps[sectors]) - gaps[sectors])   # leading edge
  edeg <- sdeg - span                                           # trailing edge
  mid  <- (sdeg - span / 2) * pi / 180

  ## a radial label whose sector sits on the left half of the circle is turned
  ## round, which moves the first character of the string from the inner end of
  ## the label to the outer end. The mark is therefore put on whichever end of
  ## the symbol comes to lie next to the circle, and every mark sits on the
  ## same ring, at the same distance from the band, all round both panels.
  degm  <- (sdeg - span / 2) %% 360
  flip  <- degm >= 90 & degm < 270
  fg    <- flip[seq_len(nrow(genes))]
  gdisp <- if_else(genes$selective,
                   if_else(fg, paste(gene_order, MARK), paste(MARK, gene_order)),
                   gene_order)

  list(d = d, genes = genes, hmord = hmord, sectors = sectors, gaps = gaps,
       labs = c(gdisp, pretty(hmord)),
       ## face of each sector's own label, and which block it belongs to
       fnt  = c(rep(FNT_GENE, nrow(genes)), rep(FNT_HM, length(hmord))),
       disp = setNames(c(gdisp, pretty(hmord)), sectors),
       sdeg = setNames(sdeg, sectors), edeg = setNames(edeg, sectors),
       is_g = c(rep(TRUE, nrow(genes)), rep(FALSE, length(hmord))),
       sel  = setNames(genes$selective, gene_order),
       arc  = (span + gaps[sectors]) * pi / 180,            # label pitch, rad
       mid = mid, ux = cos(mid), uy = sin(mid),
       grid.col = c(setNames(if_else(genes$grp == "Both", col_both, clade_col), genes$gene),
                    setNames(hm_col[hmord], hmord)),
       chord_col = setNames(adjustcolor(hm_col[d$HM], 0.55), NULL))
}

## ---- solve type size and circle together, on the region already open --------
## Type size is set by the arc each sector owns and the circle by the labels, so
## the two are solved by iteration instead of being typed in. Each sector
## reserves only what its own label needs in the direction that label runs, so
## the short gene names on one side of the circle stop paying for the long
## process names on the other and the circle grows into the space they leave.
## Nothing here depends on which gene sets happen to be present.
fit <- function(g, cex_fix = NA) {
  pw <- par("pin")[1] * 25.4                   # the box to fill, mm
  ph <- par("pin")[2] * 25.4
  ch <- strheight("M", "inches", cex = 1, font = 2) * 25.4
  ## one type size per block -- gene names, process names -- each as large as
  ## the arcs of its own block allow
  cx <- if (length(cex_fix) == 2) cex_fix else c(1, 1)
  R  <- min(pw, ph) / 2
  for (i in 1:40) {                            # damped, converges in a few steps
    cl <- ifelse(g$is_g, cx[1], cx[2])
    w  <- mapply(function(l, m, ft) strwidth(l, "inches", cex = m, font = ft),
                 g$labs, cl, g$fnt) * 25.4 + ch * cl
    e  <- RLAB + w / R                         # how far each label reaches, radii
    xr <- max(c(e * g$ux,  1.02)); xl <- max(c(-e * g$ux, 1.02))
    yt <- max(c(e * g$uy,  1.02)); yb <- max(c(-e * g$uy, 1.02))
    R  <- min(pw / (xr + xl), ph / (yt + yb)) / EXP
    if (length(cex_fix) != 2) {
      K  <- R * RLAB / (ch * LEAD)
      cx <- (cx + c(min(g$arc[g$is_g]), min(g$arc[!g$is_g])) * K) / 2
    }
  }
  list(cex = cx, R = R, ch = ch, xl = xl, xr = xr, yb = yb, yt = yt,
       binds = if (pw / (xr + xl) < ph / (yt + yb)) "x" else "y")
}

## ---- one circos panel -------------------------------------------------------
panel <- function(cl, clade_col, cex_lab = NA, R_fix = NA, fresh = TRUE) {
  g <- geom(cl, clade_col)
  if (fresh) plot.new()                        # establish this panel's region
  f <- fit(g, cex_lab)
  cex_lab <- f$cex                             # c(genes, processes)
  d <- g$d; genes <- g$genes; hmord <- g$hmord; sectors <- g$sectors
  ## both circles are drawn at the same radius: the panel that needs the least
  ## room simply keeps the margin it does not use, so the two clades are read
  ## at one scale
  k <- if (is.na(R_fix)) 1 else f$R / R_fix

  par(new = TRUE)                              # keep circlize on that region
  circos.clear()
  circos.par(start.degree = 90, gap.after = unname(g$gaps[sectors]),
             canvas.xlim = c(-f$xl, f$xr) * k, canvas.ylim = c(-f$yb, f$yt) * k,
             points.overflow.warning = FALSE, cell.padding = c(0, 0, 0, 0))

  chordDiagram(
    d[, c("gene", "HM")], order = sectors,
    grid.col = g$grid.col, col = g$chord_col, grid.border = "white",
    transparency = 0.4, link.lwd = 0.3, link.border = adjustcolor("white", 0.4),
    link.sort = TRUE, link.decreasing = TRUE,
    annotationTrack = "grid", annotationTrackHeight = mm_h(4),
    preAllocateTracks = list(list(track.height = mm_h(LBL)),  # 1: outer labels
                             list(track.height = mm_h(9))))   # 2: module band

  ## the names, each preceded by the selectivity mark when it carries one
  circos.track(track.index = 1, bg.border = NA, panel.fun = function(x, y) {
    si <- CELL_META$sector.index
    hm <- si %in% hmord
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], g$disp[si],
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),
                cex  = if (hm) cex_lab[2] else cex_lab[1],
                font = if (hm) FNT_HM else FNT_GENE,
                family = FAM, col = "grey12")
  })

  ## module bands: every module is labelled, whatever its width. The name is
  ## sized on the band PLUS the empty gap that flanks it on either side, which
  ## belongs to no other label, so a module holding a single gene still gets a
  ## readable name instead of a sliver of type. If running along the arc still
  ## does not pay, the name is turned radially, where the limit is the band
  ## thickness and not the arc.
  r_top <- get.cell.meta.data("cell.top.radius", sector.index = sectors[1],
                              track.index = 2)
  r_bot <- get.cell.meta.data("cell.bottom.radius", sector.index = sectors[1],
                              track.index = 2)
  r_mid  <- (r_top + r_bot) / 2
  mm_deg <- r_mid / mm_h(1) * pi / 180          # mm of arc per degree, at the band
  BAND <- 9                                     # module band thickness, mm
  gp   <- g$gaps[sectors]                       # gap that follows each sector
  for (m in unique(genes$module)) {
    gs   <- genes$gene[genes$module == m]
    i1   <- match(gs[1], sectors); i2 <- match(tail(gs, 1), sectors)
    bs   <- g$sdeg[i1]; be <- g$edeg[i2]                             # the band
    free <- (bs - be) + 0.85 * ((if (i1 > 1) gp[i1 - 1] else gp[length(gp)]) +
                                gp[i2])                              # + its gaps
    w_mm <- strwidth(m, "inches", cex = 1, font = 2) * 25.4
    cx   <- min(min(free * mm_deg / w_mm, BAND / (f$ch * 1.4)) * 0.9, CX_MOD)
    fits <- (bs - be) * mm_deg >= w_mm * cx * 1.08
    highlight.sector(gs, track.index = 2, col = adjustcolor(clade_col, 0.92),
                     border = "white", text = if (fits) m else "",
                     facing = "bending.inside", niceFacing = TRUE, cex = cx,
                     font = 2, family = FAM, text.col = "white")
    if (!fits) {
      ## a band too narrow to hold its own name is widened locally into the
      ## empty gaps beside it, which belong to no other label, and the name is
      ## set flat on that tab -- white on the clade's colour, the same ink and
      ## the same ground as every name that does fit inside its band
      half <- min(w_mm * cx * 1.18 / 2 / mm_deg, free / 2)
      th   <- (bs + be) / 2
      clip(par("usr")[1], par("usr")[2], par("usr")[3], par("usr")[4])
      draw.sector(th + half, th - half, rou1 = r_top, rou2 = r_bot,
                  col = adjustcolor(clade_col, 0.92), border = "white", lwd = 1.2)
      srt <- (th + 90) %% 360
      if (srt > 90 && srt < 270) srt <- srt - 180
      text(r_mid * cos(th * pi / 180), r_mid * sin(th * pi / 180), m,
           srt = srt, adj = c(0.5, 0.5), cex = cx, font = 2, family = FAM,
           col = "white")
    }
  }
  highlight.sector(hmord, track.index = 2, col = "#EAEAEE", border = "white",
                   text = "Hallmark processes", facing = "bending.inside",
                   niceFacing = TRUE, cex = 2.9, font = 2, family = FAM,
                   text.col = "grey35")
  at18 <- 12 * 180 / (25.4 * par("din")[1])
  message(sprintf(
    "[fit] %s: %d sectors | genes cex %.2f (%.1f pt at 18 cm) | processes cex %.2f (%.1f pt) | radius %.0f mm | box %.2f x %.2f, binds %s",
    cl, length(sectors), cex_lab[1], cex_lab[1] * at18, cex_lab[2], cex_lab[2] * at18,
    f$R / k, f$xl + f$xr, f$yb + f$yt, f$binds))

  title(sprintf("HPV-%s essential hubs", cl), family = FAM,
        cex.main = 7.0, font.main = 2, col.main = clade_col, line = -2.4)
  circos.clear()
}

## ---- combined figure: A7 | A9  +  node legend -------------------------------
draw_all <- function() {
  par(family = FAM)
  ## The two panels sit side by side and the frame is kept close to a landscape
  ## 16:9, because that is what makes the type as large as it can be for a
  ## reader who has the whole figure in front of them: a label is only as big as
  ## the smaller of (frame width / figure width) and (frame height / figure
  ## height) makes it, and a tall figure throws the first of those away. At this
  ## aspect the width of a panel and the height of the row run out at the same
  ## moment, so no radius is left unused in either direction.
  layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(9.1, 0.9))
  ## both panels are labelled at the same size: the region is opened once, each
  ## panel is solved on it, and the smaller of the two type sizes is the one
  ## drawn, so nothing in the figure is set in two different sizes
  par(mar = c(0.2, 0.2, 4.0, 0.2)); plot.new()
  g7 <- geom("A7", col_A7); g9 <- geom("A9", col_A9)
  cex_c <- pmin(fit(g7)$cex, fit(g9)$cex)     # c(genes, processes), both panels
  R_c   <- min(fit(g7, cex_c)$R, fit(g9, cex_c)$R)
  panel("A7", col_A7, cex_c, R_c, fresh = FALSE)
  par(mar = c(0.2, 0.2, 4.0, 0.2)); panel("A9", col_A9, cex_c, R_c)

  par(mar = c(0, 1, 0, 1)); plot.new()
  legend(x = 0.02, y = 0.98, title = "Essential-hub node", horiz = TRUE, xpd = NA,
         legend = c("Both clades", "HPV-A7 only", "HPV-A9 only"),
         pt.bg = c(col_both, col_A7, col_A9), pch = 21, pt.cex = 9.2, col = "white",
         bty = "n", cex = 5.0, title.font = 2, text.width = 0.158)
  legend(x = 0.56, y = 0.98, title = "Gene mark", xpd = NA,
         legend = paste(MARK, " Not pan-essential in DepMap"),
         bty = "n", cex = 5.0, title.font = 2, text.col = "grey12")
}

out_png <- file.path(fig_dir, "Figure_5_leading_edge_circos.png")
out_pdf <- file.path(fig_dir, "Figure_5_leading_edge_circos.pdf")
W <- 40; H <- 25.5                              # the frame, inches
png(out_png, width = W, height = H, units = "in", res = 400, bg = "white", type = "cairo"); draw_all(); dev.off()
cairo_pdf(out_pdf, width = W, height = H, bg = "white"); draw_all(); dev.off()
message("[OK] saved ", out_png)
