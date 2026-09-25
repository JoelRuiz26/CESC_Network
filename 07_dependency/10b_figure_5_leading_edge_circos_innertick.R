### ============================================================
### 6_19b_All_essential_hub_leadingedge_circos_innertick.R
###
### Variant of 10_figure_5_leading_edge_circos.R: the same essential-hub ->
### hallmark chord diagram, but the two flags (shared between clades / not
### pan-essential) are signalled INSIDE the circle instead of with a diamond
### after the gene name outside it. chordDiagram already draws one small
### rectangle per sector where its chords meet it (the "grid" annotation
### track, right between the chords and the module band); for a gene sector
### this rectangle is repurposed from a constant clade colour into the flag.
### Every gene gets exactly ONE rectangle, one solid colour -- never split or
### doubled -- one of four:
###   - purple = essential hub in BOTH clades
###   - amber  = NOT a DepMap pan-essential gene
###   - teal   = BOTH of the above at once
###   - grey   = neither -- the plain, most common case, still drawn so no
###              gene's sector is left blank
### An earlier version split a dual-flag gene's rectangle into two stacked
### halves; at a glance that read as two separate marks rather than one, so
### the four states are now four solid colours instead. The amber and teal
### were chosen by a search, not by eye (see INNER_SEL / INNER_COMBO, below):
### this same inner rectangle carries a Hallmark's own colour on the OTHER
### side of the ring, so the two flag colours also had to clear all 26 of
### those. A first pass found two hues that cleared the Hallmark palette but
### sat too close to EACH OTHER to tell apart; this pair keeps one warm hue
### and one cool hue so the two flags are as distinct from one another as
### they are from everything else. Gene names go back to plain, unmarked
### italic text; only this inner rectangle carries the flags now. Everything
### else -- modules, chords, hallmark palette, layout, sizing -- is unchanged
### from 10_figure_5_leading_edge_circos.R.
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
### Inputs:
###   07_dependency/results/6_5_essential_hub_gene_lists.tsv  (set == complete)
###   06_modules/5_2_3_Hallmark{7,9}_louv_max.tsv   (leading edge = geneID)
###   -> figures/Figure_5_leading_edge_circos_innertick.png / .pdf
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
LEAD <- 1.04    # minimum white space between two neighbouring labels
## gene names are italic, as gene symbols are set in print, and plain: no
## colour, no mark. The flags live in the inner grid rectangle instead (see
## INNER_BOTH / INNER_SEL / INNER_COMBO / INNER_PLAIN, drawn per gene in
## panel(), below) -- one gene, one rectangle, one solid colour, always.
FNT_GENE <- 3; FNT_HM <- 2
INNER_BOTH  <- col_both     # purple -- essential hub in both clades
INNER_SEL   <- "#A97200"    # amber -- not pan-essential in DepMap
INNER_COMBO <- "#007F7F"    # teal  -- essential hub in both clades AND not
                            # pan-essential, one gene at once
## INNER_SEL and INNER_COMBO were picked by search, not by eye, then checked
## against each other: each of the desaturated Hallmark colours (hm_col,
## below) and the two clade colours were treated as "forbidden", and
## candidates were scored by their minimum RGB distance to every one of them,
## restricted to a hue band that stays clear of the red/blue/purple already
## spoken for. A first pass picked two hues from the same, cooler half of
## that band (teal and green); both cleared the Hallmark/clade colours, but
## sat close enough to one another to be mistaken for the same colour. This
## pair is one warm hue and one cool hue instead -- as far apart from each
## other as they are from everything else -- so the two flags read apart
## from one another as clearly as from the rest of the figure. An even
## warmer first choice for INNER_SEL sat only 20 RGB units from the Hallmark
## APOPTOSIS colour -- about as close as two colours can be before they read
## as the same one -- which is what sent this search looking in the first
## place; amber keeps some of that warmth while clearing every Hallmark
## colour by at least as much as the paper's own blue and red already do.
INNER_PLAIN <- "grey80"     # neither -- still drawn, a plain grey tick, so no
                            # gene's sector is left blank
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
message("[check] ", length(hm_all), " hallmark labels verbatim, single line")

## ---- the sectors of one panel, before anything is drawn ---------------------
## the whole figure is drawn ~8% larger than the original frame (in, mm and pt
## all scale together) purely for legibility; SCALE is the single dial for it
SCALE <- 1.08
RLAB <- 0.97   # radius the outer labels start at, in circle radii
LBL  <- 6 * SCALE     # thickness of the track the outer labels sit in, mm
EXP  <- 1.02   # the 4 % axis expansion R adds at both ends of a plot range

geom <- function(cl, clade_col) {
  d <- map |> filter(clade == cl)
  ## within a module, genes are grouped by which flag they carry -- both
  ## flags, shared-only, selective-only, then plain -- and alphabetical only
  ## within that group, so every flag colour forms one run around the ring
  ## instead of being scattered through the module by name
  genes <- d |> distinct(gene, module, grp, selective) |>
    mutate(flag_cat = dplyr::case_when(
      grp == "Both" & selective  ~ 1L,
      grp == "Both" & !selective ~ 2L,
      grp != "Both" & selective  ~ 3L,
      TRUE                       ~ 4L)) |>
    arrange(module, flag_cat, gene) |>
    dplyr::select(-flag_cat)
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

  ## gene names carry no mark, so fit() only ever measures the bare symbol
  gdisp <- gene_order

  ## the flag colour, one solid value per gene, set here and only here: this
  ## IS the grid.col chordDiagram itself draws, not a layer painted on top of
  ## it afterwards. An earlier version drew a grey base rectangle here and a
  ## second, coloured one on top in panel() -- the two never lined up exactly
  ## at their shared edge, leaving a sliver of the grey base showing beside
  ## the real colour, i.e. every gene appeared to carry two bars instead of
  ## one. Setting the true colour directly, in the one call chordDiagram uses
  ## to draw the rectangle, removes the second layer entirely.
  inner_col <- dplyr::case_when(
    genes$grp == "Both" & genes$selective ~ INNER_COMBO,
    genes$grp == "Both"                   ~ INNER_BOTH,
    genes$selective                       ~ INNER_SEL,
    TRUE                                   ~ INNER_PLAIN)

  list(d = d, genes = genes, hmord = hmord, sectors = sectors, gaps = gaps,
       labs = c(gdisp, pretty(hmord)),
       ## one face for every gene label; hallmark names keep their own
       fnt  = c(rep(FNT_GENE, nrow(genes)), rep(FNT_HM, length(hmord))),
       hm_disp = setNames(pretty(hmord), hmord),
       sdeg = setNames(sdeg, sectors), edeg = setNames(edeg, sectors),
       is_g = c(rep(TRUE, nrow(genes)), rep(FALSE, length(hmord))),
       arc  = (span + gaps[sectors]) * pi / 180,            # label pitch, rad
       mid = mid, ux = cos(mid), uy = sin(mid),
       grid.col = c(setNames(inner_col, genes$gene),
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
    ## this track is where the flag actually lives now (see the loop below
    ## panel.fun), so it is given real thickness -- three times what a plain
    ## chord-anchor strip needs -- rather than the thin 4 mm sliver a chord
    ## diagram uses by default, so the mark reads clearly at this figure size
    annotationTrack = "grid", annotationTrackHeight = mm_h(12 * SCALE),
    preAllocateTracks = list(list(track.height = mm_h(LBL)),  # 1: outer labels
                             list(track.height = mm_h(9 * SCALE))))   # 2: module band

  ## every label is now plain: the flags are carried by the inner grid
  ## rectangle (grid.col, set in geom()), not by anything drawn here
  circos.track(track.index = 1, bg.border = NA, panel.fun = function(x, y) {
    si <- CELL_META$sector.index
    hm <- si %in% hmord
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], if (hm) g$hm_disp[si] else si,
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
  BAND <- 9 * SCALE                             # module band thickness, mm
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
        cex.main = 11.5, font.main = 2, col.main = clade_col, line = -3.0)
  circos.clear()
}

## ---- combined figure: A7 | A9  +  node legend -------------------------------
draw_all <- function() {
  par(family = FAM)
  ## The two panels are stacked over a thin band that carries the key. The solved
  ## drawing is about 2.65 radii wide and 3.09 high, so a stacked pair reaches
  ## its largest circle when the frame is 2 x 3.09 / 2.65 = 2.33 times as tall as
  ## it is wide. Beyond that ratio the extra height buys nothing and only
  ## stretches the figure, so the frame is set at it and the band is added on top.
  layout(matrix(c(1, 2, 3), nrow = 3), heights = c(4.85, 4.85, 0.62))
  par(mar = c(0.2, 0.2, 4.0, 0.2)); plot.new()
  g7 <- geom("A7", col_A7); g9 <- geom("A9", col_A9)
  cex_c <- pmin(fit(g7)$cex, fit(g9)$cex)
  R_c   <- min(fit(g7, cex_c)$R, fit(g9, cex_c)$R)
  message(sprintf("[check] type solved: genes %.2f | processes %.2f | radius %.0f mm | panel %.1f x %.1f in",
                  cex_c[1], cex_c[2], R_c, par("pin")[1], par("pin")[2]))
  panel("A7", col_A7, cex_c, R_c, fresh = FALSE)
  par(mar = c(0.2, 0.2, 4.0, 0.2)); panel("A9", col_A9, cex_c, R_c)

  ## the key runs along the foot of the figure, clear of both circles: each
  ## entry is a small filled square in the colour the inner rectangle takes,
  ## and the sentence it means. The scale is solved from the actual string
  ## widths (same principle fit() uses for the circles) rather than picked
  ## by hand, so it can never overflow the row or clip a caption.
  par(mar = c(0, 1, 0, 1)); plot.new()
  key_y      <- 0.52
  SWATCH     <- "■"  # a small filled square -- the same shape the inner
                          # grid rectangle draws, just as a key glyph here
  MARK_RATIO <- 1.05     # the swatch, about the size of its caption
  TARGET_W   <- 0.94     # the key fills almost the whole row width, centred
  GAP        <- 0.016    # small breathing room between pieces
  entries <- list(
    list(col = INNER_BOTH,  caption = "essential hub in both clades"),
    list(col = INNER_SEL,   caption = "not pan-essential in DepMap"),
    list(col = INNER_COMBO, caption = "essential hub in both clades and not pan-essential"))

  ## every entry as an ordered (swatch, caption) pair, one swatch each -- the
  ## same one rectangle, one colour rule the figure itself follows
  pieces <- function(e, lead) list(
    list(txt = paste0(lead, SWATCH),    col = e$col,    font = 1, rel = MARK_RATIO),
    list(txt = paste0("  ", e$caption), col = "grey20", font = 1, rel = 1))
  all_pieces <- c(pieces(entries[[1]], ""), pieces(entries[[2]], "        "),
                  pieces(entries[[3]], "        "))

  ## the scale that makes the whole key exactly TARGET_W wide, solved from the
  ## actual string widths (same principle fit() uses for the circles) rather
  ## than picked by hand, then centred in the row. GAP itself is already a
  ## plot-fraction constant, not something that scales with the solved cex, so
  ## it is subtracted from TARGET_W up front rather than folded into w_at_1 --
  ## otherwise, close to TARGET_W = 1, the unscaled gaps push the last caption
  ## past the edge of the row.
  w_at_1  <- sum(vapply(all_pieces, function(p)
    strwidth(p$txt, cex = p$rel, font = p$font), numeric(1)))
  key_cex <- (TARGET_W - GAP * (length(all_pieces) - 1)) / w_at_1
  message(sprintf("[check] legend solved: caption cex %.2f, swatch cex %.2f",
                  key_cex, key_cex * MARK_RATIO))

  x0 <- (1 - TARGET_W) / 2
  for (p in all_pieces) {
    text(x0, key_y, p$txt, font = p$font, col = p$col, family = FAM,
         cex = key_cex * p$rel, adj = c(0, 0.5), xpd = NA)
    x0 <- x0 + strwidth(p$txt, cex = key_cex * p$rel, font = p$font) + GAP
  }
}

out_png <- file.path(fig_dir, "Figure_5_leading_edge_circos_innertick.png")
out_pdf <- file.path(fig_dir, "Figure_5_leading_edge_circos_innertick.pdf")
W <- 25 * SCALE; H <- 60 * SCALE                # the frame, inches (see SCALE above)
png(out_png, width = W, height = H, units = "in", res = 400, bg = "white", type = "cairo"); draw_all(); dev.off()
cairo_pdf(out_pdf, width = W, height = H, bg = "white"); draw_all(); dev.off()
message("[OK] saved ", out_png)
