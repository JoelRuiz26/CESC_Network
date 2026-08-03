### ============================================================
### 6_14_Selective_hub_targets_circos.R
###
### Per-network circos of the SELECTIVE essential hubs (essential hubs that are
### NOT DepMap pan-essential) and their co-expression connectivity.
###   HPV-A7: 47 selective hubs  (11 shared + 36 A7-unique)
###   HPV-A9: 54 selective hubs  (11 shared + 43 A9-unique)
###
### Layout (same for both panels):
###   * top block  = all selective hubs, labelled. 11 shared-in-both genes in
###                  purple; network-unique hubs in the clade colour. An outer
###                  bar gives each hub's TOTAL number of co-expression targets.
###   * the 11 shared (purple) genes fan out to ALL of their individual targets
###     (bottom block), emphasised: this is the point, that the SAME genes reach
###     largely DIFFERENT partners in each clade.
###   * chords across the top block = co-expression edges among the selective
###     hubs themselves (each hub's link to other essential hubs).
###
###   6_14_selective_hub_target_circos.png / .pdf   two circos (A7 | A9)
### ============================================================
suppressPackageStartupMessages({
  library(tidyverse); library(vroom); library(igraph); library(circlize)
})

sub_dir <- "~/CESC_Network/5_Network_analisis/5_1_Subnetwork"
res_dir <- "~/CESC_Network/6_DepMap_cacu/results"

col_common <- "#8C79B5"   # shared-in-both selective hubs (purple)
col_A7     <- "#4E79A7"   # HPV-A7-unique selective hubs
col_A9     <- "#CB6A5C"   # HPV-A9-unique selective hubs
col_link   <- "#9AA0A6"   # hub-hub (inter-essential-hub) chord
col_share  <- "#C9BAE6"   # light purple = target of a purple gene shared by BOTH clades
                          # divergent targets are drawn in each network's own colour (col_A7 / col_A9)

A7 <- readRDS(file.path(sub_dir, "5_2_4_Subnetwork_A7_max.rds"))
A9 <- readRDS(file.path(sub_dir, "5_2_4_Subnetwork_A9_max.rds"))
V(A7)$name <- toupper(V(A7)$name); V(A9)$name <- toupper(V(A9)$name)

elbow <- function(v) { n <- length(v); x <- (seq_len(n) - 1) / (n - 1)
  y <- (v - v[1]) / (v[n] - v[1]); v[which.max(abs(y - x))] }

tab <- vroom(file.path(res_dir, "6_2_gene_dep_degree_table.tsv"), show_col_types = FALSE) %>%
  filter(!is.na(degree), !is.na(mean_dep)) %>% distinct(network, gene, .keep_all = TRUE)
th <- tab %>% group_by(network) %>%
  summarise(k = elbow(sort(degree, decreasing = TRUE)),
            d = elbow(sort(mean_dep)), .groups = "drop")
ess <- tab %>% left_join(th, by = "network") %>%
  filter(degree >= k, mean_dep <= d) %>% mutate(gene = toupper(gene))
pan <- toupper(vroom(file.path(res_dir, "6_1_common_essential_genes.tsv"),
                     delim = "\t", show_col_types = FALSE)$gene_u)
sel <- ess %>% filter(!gene %in% pan)
sel7 <- sel$gene[sel$network == "HPV-A7"]; sel9 <- sel$gene[sel$network == "HPV-A9"]
common <- intersect(sel7, sel9)

common_order <- c("NDUFA13", "NDUFA2", "NDUFA9", "NDUFS3",
                  "MRPS15", "MRPS2", "SSBP1", "ZNHIT1", "ZWINT", "EGFR", "FOSL1")
common_order <- common_order[common_order %in% common]

## for each shared gene, the targets it connects to in BOTH clades (shared);
## any other target is specific to one clade -> used to colour the purple fans
nb_in <- function(G, g) if (g %in% V(G)$name) neighbors(G, g)$name else character(0)
share_by_gene <- setNames(
  lapply(common, function(g) intersect(nb_in(A7, g), nb_in(A9, g))), common)

## global max degree among selective hubs (comparable bar heights) ---------
gmax <- max(igraph::degree(A7, v = intersect(sel7, V(A7)$name)),
            igraph::degree(A9, v = intersect(sel9, V(A9)$name)))

WID_HUB <- 2.2; WID_TGT <- 0.5

draw_panel <- function(G, netlab, sel_set, col_unique) {
  present <- intersect(sel_set, V(G)$name)
  purple  <- common_order[common_order %in% present]
  deg     <- setNames(as.numeric(igraph::degree(G, v = present)), present)
  uniq    <- setdiff(present, purple)
  uniq    <- uniq[order(deg[uniq], decreasing = TRUE)]
  hubs    <- c(purple, uniq)
  hub_class <- setNames(ifelse(hubs %in% purple, "common", "unique"), hubs)
  hub_col   <- ifelse(hub_class[hubs] == "common", col_common, col_unique)
  names(hub_col) <- hubs

  nb <- function(g) neighbors(G, g)$name
  purple_nb <- setNames(lapply(purple, nb), purple)

  ## bottom block = individual targets of the purple genes that are NOT hubs
  targetsB <- setdiff(unique(unlist(purple_nb)), hubs)
  prim <- sapply(targetsB, function(t)
    which(vapply(purple, function(g) t %in% purple_nb[[g]], logical(1)))[1])
  targetsB <- targetsB[order(prim, targetsB)]
  prim <- prim[order(prim, names(prim))]

  ## inter-hub edges among the selective hubs (each hub's link to other hubs)
  eH <- as_edgelist(induced_subgraph(G, hubs))

  sectors <- c(hubs, targetsB)
  xlim <- cbind(0, c(rep(WID_HUB, length(hubs)), rep(WID_TGT, length(targetsB))))
  rownames(xlim) <- sectors

  gap <- c(rep(2.2, length(hubs) - 1), 7,
           rep(0.1, length(targetsB) - 1), 7)

  circos.clear()
  circos.par(start.degree = 80, gap.after = gap, cell.padding = c(0, 0, 0, 0),
             points.overflow.warning = FALSE,
             canvas.xlim = c(-1.4, 1.4), canvas.ylim = c(-1.4, 1.4))
  circos.initialize(factors = factor(sectors, levels = sectors), xlim = xlim)

  ## track 1: degree bar (total targets, hubs only) + hub labels outside ----
  circos.track(ylim = c(0, 1), track.height = 0.14, bg.border = NA,
    panel.fun = function(x, y) {
      s <- CELL_META$sector.index; w <- CELL_META$xrange
      if (s %in% hubs) {
        circos.rect(0.15 * w, 0, 0.85 * w, deg[s] / gmax,
                    col = hub_col[s], border = NA)
        cl <- hub_col[s]; fw <- if (hub_class[s] == "common") 2 else 1
        circos.text(w / 2, 1.35, s, facing = "clockwise", niceFacing = TRUE,
                    adj = c(0, 0.5), cex = 0.5, font = fw, col = cl)
      }
    })
  ## track 2: coloured identity ring ---------------------------------------
  circos.track(ylim = c(0, 1), track.height = 0.03, bg.border = NA,
    panel.fun = function(x, y) {
      s <- CELL_META$sector.index; w <- CELL_META$xrange
      if (s %in% hubs) circos.rect(0, 0, w, 1, col = hub_col[s], border = NA)
      else             circos.rect(0, 0.2, w, 0.8, col = "#DCDCDC", border = NA)
    })

  ## inter-hub chords (grey, behind) ---------------------------------------
  if (nrow(eH) > 0) for (i in seq_len(nrow(eH))) {
    circos.link(eH[i, 1], WID_HUB / 2, eH[i, 2], WID_HUB / 2,
                col = adjustcolor(col_link, 0.5), lwd = 0.8, h.ratio = 0.5)
  }
  ## purple fans: each shared gene -> all its individual targets, coloured by
  ## whether the gene reaches that target in BOTH clades (shared, light) or only
  ## here (clade-specific, dark). Draw shared first so the dark divergent links
  ## sit on top and the pattern is legible.
  fan <- list()
  for (g in purple) {
    ts <- intersect(purple_nb[[g]], targetsB)
    if (!length(ts)) next
    fan[[g]] <- list(ts = ts,
                     xs = seq(0.18 * WID_HUB, 0.82 * WID_HUB, length.out = length(ts)),
                     shared = ts %in% share_by_gene[[g]])
  }
  for (g in names(fan)) { f <- fan[[g]]; for (j in which(f$shared))
    circos.link(g, f$xs[j], f$ts[j], WID_TGT / 2,
                col = adjustcolor(col_share, 0.85), lwd = 0.45, h.ratio = 0.62) }
  for (g in names(fan)) { f <- fan[[g]]; for (j in which(!f$shared))
    circos.link(g, f$xs[j], f$ts[j], WID_TGT / 2,
                col = adjustcolor(col_unique, 0.6), lwd = 0.55, h.ratio = 0.62) }
  text(-0.92, 1.28, netlab, cex = 1.4, font = 2, adj = 0)
  circos.clear()
}

make_fig <- function(open_dev) {
  open_dev()
  layout(matrix(1:2, 1, 2)); par(mar = c(2.2, 1, 3, 1), family = "sans")
  draw_panel(A7, "HPV-A7", sel7, col_A7)
  draw_panel(A9, "HPV-A9", sel9, col_A9)
  par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0)); plot.new()
  text(0.5, 0.975, "Bar height = total number of co-expression targets of each hub.  Ribbons trace every target of the 11 shared genes:",
       cex = 0.82, col = "grey30")
  text(0.5, 0.945, "light purple = target reached in both clades;  clade colour (blue in HPV-A7, red in HPV-A9) = target reached only in this clade.",
       cex = 0.82, col = "grey30")
  legend("bottom", horiz = TRUE, bty = "n", cex = 0.92, inset = 0.004,
         seg.len = 1.3, lwd = c(NA, NA, NA, 3), pch = c(15, 15, 15, NA),
         pt.cex = 1.7,
         col = c(col_common, col_A7, col_A9, col_link),
         legend = c("shared selective hub (both clades)", "HPV-A7-unique selective hub",
                    "HPV-A9-unique selective hub", "co-expression edge between hubs"))
  dev.off()
}

make_fig(function() png(file.path(res_dir, "6_14_selective_hub_target_circos.png"),
                        width = 17, height = 9.6, units = "in", res = 400))
make_fig(function() cairo_pdf(file.path(res_dir, "6_14_selective_hub_target_circos.pdf"),
                              width = 17, height = 9.6))

cat(sprintf("selective hubs  A7=%d (11 shared + %d unique)  A9=%d (11 shared + %d unique)\n",
            length(sel7), length(sel7) - 11, length(sel9), length(sel9) - 11))
message("Saved 6_14_selective_hub_target_circos.png / .pdf to ", res_dir)
