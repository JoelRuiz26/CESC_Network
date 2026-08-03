### ============================================================
### 6_15_Shared_hub_fusion_circos.R
###
### FUSED DIVERGENT connectivity of the 11 genes that are selective essential
### hubs in BOTH clades (the "purple" genes). Both networks are merged into one
### diagram, but ONLY the clade-divergent edges are drawn:
###   * blue  = target reached by the gene only in HPV-A7
###   * red   = target reached by the gene only in HPV-A9
### Targets reached in both clades (convergent) are removed, as are targets that
### are reached through HPV-A7 for some hubs and HPV-A9 for others (mixed), so
### every drawn target is unambiguously HPV-A7-divergent or HPV-A9-divergent.
### Each target is one sector (merged), sized by the number of shared hubs that
### reach it, so genes co-targeted by several hubs stand out.
###
###   6_15_shared_hub_fusion_circos.png / .pdf
### ============================================================
suppressPackageStartupMessages({
  library(tidyverse); library(vroom); library(igraph); library(circlize)
})

sub_dir <- "~/CESC_Network/5_Network_analisis/5_1_Subnetwork"
res_dir <- "~/CESC_Network/6_DepMap_cacu/results"

col_hub  <- "#8C79B5"                 # shared (purple) hubs
col_a7   <- "#4E79A7"; col_a9 <- "#CB6A5C"; col_both <- "#7E5FB0"
tcol     <- c(A7 = "#9DB6D3", both = "#B7A6D9", A9 = "#E0AEA6")

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
common <- intersect(sel$gene[sel$network == "HPV-A7"], sel$gene[sel$network == "HPV-A9"])
co <- c("NDUFA13","NDUFA2","NDUFA9","NDUFS3","MRPS15","MRPS2","SSBP1","ZNHIT1","ZWINT","EGFR","FOSL1")
hubs <- c(intersect(co, common), setdiff(common, co))

nb <- function(G, g) if (g %in% V(G)$name) neighbors(G, g)$name else character(0)

## edges of each shared hub, tagged: reached only in A7, only in A9, or in both
E <- map_dfr(hubs, function(h) {
  n7 <- nb(A7, h); n9 <- nb(A9, h)
  bind_rows(tibble(hub = h, target = setdiff(n7, n9), tag = "A7"),
            tibble(hub = h, target = setdiff(n9, n7), tag = "A9"),
            tibble(hub = h, target = intersect(n7, n9), tag = "both"))
}) %>% filter(!target %in% hubs)

## classify targets: only in HPV-A7 (blue), only in HPV-A9 (red), or in both
## clades (purple: a shared hub reaches the same target in both networks, or
## different shared hubs reach it in each clade)
tc <- E %>% group_by(target) %>%
  summarise(a7 = any(tag %in% c("A7", "both")), a9 = any(tag %in% c("A9", "both")),
            nhub = n_distinct(hub), .groups = "drop") %>%
  mutate(cls = if_else(a7 & a9, "both", if_else(a7, "A7", "A9"))) %>%
  arrange(factor(cls, levels = c("A7", "both", "A9")), desc(nhub), target)
targets <- tc$target
tcls  <- setNames(tc$cls, targets); tnhub <- setNames(tc$nhub, targets)

## geometry -----------------------------------------------------------------
WID_HUB <- 5
tw <- 0.3 + 0.5 * (tnhub - 1)                         # bigger target = more hubs share it
sectors <- c(hubs, targets)
xlim <- cbind(0, c(rep(WID_HUB, length(hubs)), tw)); rownames(xlim) <- sectors

gap <- c(rep(1.4, length(hubs) - 1), 7,
         rep(0.12, length(targets) - 1), 7)

draw <- function() {
  circos.clear()
  circos.par(start.degree = 90, gap.after = gap, cell.padding = c(0, 0, 0, 0),
             points.overflow.warning = FALSE,
             canvas.xlim = c(-1.4, 1.4), canvas.ylim = c(-1.4, 1.4))
  circos.initialize(factors = factor(sectors, levels = sectors), xlim = xlim)

  circos.track(ylim = c(0, 1), track.height = 0.05, bg.border = NA,
    panel.fun = function(x, y) {
      s <- CELL_META$sector.index; w <- CELL_META$xrange
      if (s %in% hubs) {
        circos.rect(0, 0, w, 1, col = col_hub, border = NA)
        circos.text(w / 2, 1.45, s, facing = "clockwise", niceFacing = TRUE,
                    adj = c(0, 0.5), cex = 0.72, font = 2, col = "#54338C")
      } else circos.rect(0, 0, w, 1, col = tcol[tcls[s]], border = NA)
    })

  ## draw divergent (A7/A9) links first, then the shared (both) links on top
  El <- E %>% arrange(match(tag, c("A7", "A9", "both")))
  lcol <- c(A7 = col_a7, A9 = col_a9, both = col_both)
  lalp <- c(A7 = 0.4, A9 = 0.4, both = 0.7); llwd <- c(A7 = 0.5, A9 = 0.5, both = 0.8)
  hx <- setNames(rep(WID_HUB / 2, length(hubs)), hubs)
  for (i in seq_len(nrow(El))) {
    tg <- El$tag[i]
    circos.link(El$hub[i], hx[[El$hub[i]]], El$target[i], tw[[El$target[i]]] / 2,
                col = adjustcolor(lcol[tg], lalp[tg]), lwd = llwd[tg], h.ratio = 0.6)
  }
  circos.clear()
}

make_fig <- function(open_dev) {
  open_dev(); par(mar = c(1, 1, 2.5, 1), family = "sans"); draw()
  text(0, 1.34, "Divergent connectivity of the 11 shared selective essential hubs",
       cex = 1.15, font = 2)
  legend("bottomright", bty = "n", cex = 0.9, inset = 0.02, lwd = 3, seg.len = 1.3,
         col = c(col_a7, col_a9, col_both),
         legend = c("edge only in HPV-A7", "edge only in HPV-A9",
                    "edge in both clades (same target)"))
  legend("bottomleft", bty = "n", cex = 0.9, inset = 0.02, pch = 15, pt.cex = 1.7,
         col = col_hub, legend = "shared selective hub (both clades)")
  text(0, -1.36, "Blue / red = clade-divergent edges; purple = the same target reached in both clades.  Target sector size grows with the number of shared hubs reaching it.",
       cex = 0.66, col = "grey35")
  dev.off()
}

make_fig(function() png(file.path(res_dir, "6_15_shared_hub_fusion_circos.png"),
                        width = 12, height = 12.5, units = "in", res = 400))
make_fig(function() cairo_pdf(file.path(res_dir, "6_15_shared_hub_fusion_circos.pdf"),
                              width = 12, height = 12.5))

cat(sprintf("11 shared hubs | targets A7=%d both=%d A9=%d | edges: %s\n",
            sum(tcls == "A7"), sum(tcls == "both"), sum(tcls == "A9"),
            paste(names(table(E$tag)), table(E$tag), collapse = " ")))
message("Saved 6_15_shared_hub_fusion_circos.png / .pdf to ", res_dir)
