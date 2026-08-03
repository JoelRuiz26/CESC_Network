### ============================================================
### 6_13_Essential_hub_circos.R
###
### Circos / chord diagrams built on the 168 genes that are ESSENTIAL HUBS
### in BOTH the HPV-A7 and HPV-A9 co-expression networks (same data-driven
### elbow thresholds as 6_8 / 6_9 / 6_10 / 6_11 / 6_12).
###
### The point: although these genes are shared essential hubs, their
### connectivity (their co-expression partners) is largely different between
### clades. Median neighbour-set Jaccard ~ 0.10; only 127 of the ~700 hub-hub
### edges are shared (edge Jaccard 0.17).
###
###   6_13_A_hub_wiring_circos.png/.pdf   two circos (A7 | A9), same 168
###                                       sectors, hub-hub edges per network
###   6_13_B_jaccard_circos.png/.pdf      one circos, edges coloured
###                                       A7-only / shared / A9-only, with an
###                                       outer neighbourhood-Jaccard track
### ============================================================
suppressPackageStartupMessages({
  library(tidyverse); library(vroom); library(igraph); library(circlize)
})

sub_dir <- "~/CESC_Network/5_Network_analisis/5_1_Subnetwork"
res_dir <- "~/CESC_Network/6_DepMap_cacu/results"

## ---- palette (paper convention) ----------------------------------------
col_A7   <- "#4E79A7"   # A7-unique edges
col_A9   <- "#CB6A5C"   # A9-unique edges
col_both <- "#8C79B5"   # edges present in both networks
grid_grey <- "#B7B7B7"

## ---- networks & essential-hub definition -------------------------------
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
both <- intersect(ess$gene[ess$network == "HPV-A7"], ess$gene[ess$network == "HPV-A9"])

## ---- sector order: by A7 module, then A7 degree (descending) ------------
a7info <- tab %>% filter(network == "HPV-A7") %>% mutate(gene = toupper(gene)) %>%
  filter(gene %in% both) %>%
  transmute(gene, a7_mod = community_name, a7_deg = degree)
a9info <- tab %>% filter(network == "HPV-A9") %>% mutate(gene = toupper(gene)) %>%
  filter(gene %in% both) %>%
  transmute(gene, a9_mod = community_name, a9_deg = degree)

modnum <- function(x) as.integer(sub("^M", "", x))
gi <- a7info %>% left_join(a9info, by = "gene") %>%
  arrange(modnum(a7_mod), desc(a7_deg))
gene_order <- gi$gene
stopifnot(length(gene_order) == length(both))

## module colours (order by frequency; M2 is the dominant A7 module) -------
mod_levels <- gi %>% count(a7_mod) %>% arrange(desc(n)) %>% pull(a7_mod)
## module palette deliberately AVOIDS the edge colours (A7 blue / A9 red /
## shared purple) so the identity ring is never confused with the chords
mod_pal <- setNames(
  colorRampPalette(c("#59A14F", "#F28E2B", "#EDC948", "#76B7B2", "#9C755F",
                     "#2E6E4E", "#FF9DA7", "#B8860B", "#86BCB6", "#BAB0AC",
                     "#7F9C3A", "#D4A56A"))(length(mod_levels)),
  mod_levels)
gi$mod_col <- mod_pal[gi$a7_mod]
sector_col <- setNames(gi$mod_col, gi$gene)

## ---- hub-hub edge sets among the 168 shared essential hubs --------------
edge_key <- function(g) {
  s <- induced_subgraph(g, both)
  el <- as_edgelist(s)
  if (nrow(el) == 0) return(character(0))
  apply(el, 1, function(x) paste(sort(x), collapse = "\r"))
}
e7 <- edge_key(A7); e9 <- edge_key(A9)
shared_e <- intersect(e7, e9)
only7    <- setdiff(e7, e9)
only9    <- setdiff(e9, e7)
split_edge <- function(k) strsplit(k, "\r", fixed = TRUE)[[1]]

## ---- per-hub full-neighbourhood Jaccard (A7 vs A9 partners) -------------
jac <- sapply(gene_order, function(g) {
  n7 <- neighbors(A7, g)$name; n9 <- neighbors(A9, g)$name
  length(intersect(n7, n9)) / length(union(n7, n9))
})

## gaps: small within a module, larger between modules -------------------
mod_of  <- setNames(gi$a7_mod, gi$gene)
gaps <- sapply(seq_along(gene_order), function(i) {
  j <- if (i == length(gene_order)) 1 else i + 1
  if (mod_of[gene_order[i]] != mod_of[gene_order[j]]) 3 else 0.6
})

## ------------------------------------------------------------------------
## helper: draw one circos of hub-hub edges
##   edge_lists : named list edgeset -> colour
## ------------------------------------------------------------------------
draw_circos <- function(edge_lists, title, label_genes = NULL,
                        outer_track = NULL, outer_col = grid_grey,
                        link_alpha = 0.5, label_cex = 0.30) {
  circos.clear()
  circos.par(gap.after = gaps, cell.padding = c(0, 0, 0, 0),
             start.degree = 90, track.margin = c(0.004, 0.004),
             points.overflow.warning = FALSE)
  circos.initialize(factors = factor(gene_order, levels = gene_order),
                    xlim = c(0, 1))

  ## optional outer numeric track (e.g. neighbourhood Jaccard)
  if (!is.null(outer_track)) {
    ymax <- max(outer_track)
    circos.track(ylim = c(0, ymax), track.height = 0.10, bg.border = NA,
      panel.fun = function(x, y) {
        g <- CELL_META$sector.index
        circos.rect(0, 0, 1, outer_track[[g]], col = outer_col, border = NA)
      })
    circos.yaxis(side = "left", at = pretty(c(0, ymax), 2),
                 labels.cex = 0.4, sector.index = gene_order[1], track.index = 1)
  }

  ## module ring + optional gene labels
  circos.track(ylim = c(0, 1), track.height = 0.055, bg.border = NA,
    panel.fun = function(x, y) {
      g <- CELL_META$sector.index
      circos.rect(0, 0, 1, 1, col = sector_col[[g]], border = NA)
      if (!is.null(label_genes) && g %in% label_genes)
        circos.text(0.5, 1.4, g, facing = "clockwise", niceFacing = TRUE,
                    adj = c(0, 0.5), cex = label_cex, font = 2)
    })

  ## links (draw shared last so the contrast is visible; unique on top)
  draw_set <- function(keys, col) {
    if (!length(keys)) return(invisible())
    for (k in keys) { p <- split_edge(k)
      circos.link(p[1], 0.5, p[2], 0.5,
                  col = adjustcolor(col, link_alpha), lwd = 0.6, h.ratio = 0.6) }
  }
  for (nm in names(edge_lists)) draw_set(edge_lists[[nm]]$keys, edge_lists[[nm]]$col)

  text(0, 1.04, title, cex = 1.15, font = 2)
  circos.clear()
}

## label the most rewired shared hubs (lowest neighbourhood Jaccard) + anchors
anchor <- c("MYC", "EGFR")
lab <- unique(c(names(sort(jac))[1:12], intersect(anchor, gene_order)))

## ========================================================================
## Panel A: two circos, hub-hub wiring per network (same 168 sectors)
## ========================================================================
png(file.path(res_dir, "6_13_A_hub_wiring_circos.png"),
    width = 15, height = 8.4, units = "in", res = 400)
layout(matrix(1:2, 1, 2)); par(mar = c(1, 1, 2.5, 1), family = "sans")

draw_circos(list(only  = list(keys = only7, col = col_A7),
                 share = list(keys = shared_e, col = col_both)),
            title = "HPV-A7  hub–hub wiring", label_genes = lab)
draw_circos(list(only  = list(keys = only9, col = col_A9),
                 share = list(keys = shared_e, col = col_both)),
            title = "HPV-A9  hub–hub wiring", label_genes = lab)
## shared legend
par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0))
plot.new()
legend("bottom", horiz = TRUE, bty = "n", inset = 0.005, cex = 0.95,
       legend = c("HPV-A7-unique edge", "shared edge", "HPV-A9-unique edge"),
       col = c(col_A7, col_both, col_A9), lwd = 3, seg.len = 1.4)
dev.off()

## ========================================================================
## Panel B: one "Jaccard" circos, edges coloured A7-only/shared/A9-only,
##          outer track = per-hub neighbourhood Jaccard
## ========================================================================
png(file.path(res_dir, "6_13_B_jaccard_circos.png"),
    width = 9.5, height = 9.5, units = "in", res = 400)
par(mar = c(1, 1, 2, 1), family = "sans")

draw_circos(list(o7    = list(keys = only7,    col = col_A7),
                 o9    = list(keys = only9,    col = col_A9),
                 share = list(keys = shared_e, col = col_both)),
            title = "Divergent connectivity of the shared essential hubs",
            label_genes = lab, outer_track = jac, outer_col = "#6E6E6E",
            link_alpha = 0.55, label_cex = 0.32)
legend("bottomright", bty = "n", cex = 0.85, inset = 0.02,
       title = "co-expression edge",
       legend = c("HPV-A7 only", "shared (both)", "HPV-A9 only"),
       col = c(col_A7, col_both, col_A9), lwd = 3, seg.len = 1.4)
legend("bottomleft", bty = "n", cex = 0.8, inset = 0.02,
       legend = "outer bars: neighbourhood\nJaccard (A7 vs A9 partners)")
dev.off()

## also PDFs (vector) ------------------------------------------------------
for (dev_spec in list(
  list(f = "6_13_A_hub_wiring_circos.pdf", w = 15, h = 8.4, two = TRUE),
  list(f = "6_13_B_jaccard_circos.pdf",   w = 9.5, h = 9.5, two = FALSE))) {
  cairo_pdf(file.path(res_dir, dev_spec$f), width = dev_spec$w, height = dev_spec$h)
  if (dev_spec$two) {
    layout(matrix(1:2, 1, 2)); par(mar = c(1, 1, 2.5, 1), family = "sans")
    draw_circos(list(only = list(keys = only7, col = col_A7),
                     share = list(keys = shared_e, col = col_both)),
                "HPV-A7  hub–hub wiring", label_genes = lab)
    draw_circos(list(only = list(keys = only9, col = col_A9),
                     share = list(keys = shared_e, col = col_both)),
                "HPV-A9  hub–hub wiring", label_genes = lab)
  } else {
    par(mar = c(1, 1, 2, 1), family = "sans")
    draw_circos(list(o7 = list(keys = only7, col = col_A7),
                     o9 = list(keys = only9, col = col_A9),
                     share = list(keys = shared_e, col = col_both)),
                "Divergent connectivity of the shared essential hubs",
                label_genes = lab, outer_track = jac, outer_col = "#6E6E6E")
  }
  dev.off()
}

cat(sprintf("shared essential hubs: %d\n", length(both)))
cat(sprintf("hub-hub edges  A7=%d  A9=%d  shared=%d  (edge Jaccard %.3f)\n",
            length(e7), length(e9), length(shared_e),
            length(shared_e) / length(union(e7, e9))))
cat(sprintf("neighbourhood Jaccard  median=%.3f  mean=%.3f\n",
            median(jac), mean(jac)))
message("Saved 6_13_A / 6_13_B circos plots to ", res_dir)
