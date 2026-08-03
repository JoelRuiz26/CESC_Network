### ============================================================
### 5_2_9_Module_pair_cellcycle_organic.R
###
### Organic-network variant (companion to arc diagram 5_2_8) for the cell-cycle
### example: HPV-A7 module M8 vs HPV-A9 module M13, both enriched for the same
### hallmark set (G2M checkpoint, E2F targets, mitotic spindle).
###
### Clean version: only the hallmark-DRIVER genes are drawn (drivers + the
### drivers shared by both modules, in purple) — no grey context, no non-driver
### shared genes. A force-directed cloud per side, links to the central hubs.
###
###   -> 5_2_9_Module_pair_cellcycle_organic.png / .pdf
### ============================================================
suppressPackageStartupMessages({
  library(igraph); library(tidyverse); library(vroom); library(ggrepel)
})

ora    <- "~/CESC_Network/7_Modularity/5_2_ORA"
moddir <- "~/CESC_Network/7_Modularity/5_1_Modules"

col_A7 <- "#4E79A7"; col_A9 <- "#CB6A5C"; col_shared_drv <- "#5E3C99"

SHARED_HM <- c("G2M_CHECKPOINT", "E2F_TARGETS", "MITOTIC_SPINDLE")
hm_col <- c(G2M_CHECKPOINT = "#E69F00", E2F_TARGETS = "#009E73", MITOTIC_SPINDLE = "#D55E00")
hm_lab <- c(G2M_CHECKPOINT = "G2M checkpoint", E2F_TARGETS = "E2F targets",
            MITOTIC_SPINDLE = "Mitotic spindle")
hm_y   <- c(G2M_CHECKPOINT = 4.4, E2F_TARGETS = 0, MITOTIC_SPINDLE = -4.4)

## enforce a homogeneous minimum spacing between nodes (de-overlap)
declump <- function(L, min_d, iters = 400, step = 0.65) {
  for (i in seq_len(iters)) {
    D <- as.matrix(dist(L)); diag(D) <- Inf
    if (min(D) >= min_d) break
    n <- nrow(L); disp <- matrix(0, n, 2)
    idx <- which(D < min_d & upper.tri(D), arr.ind = TRUE)
    for (k in seq_len(nrow(idx))) {
      a <- idx[k, 1]; b <- idx[k, 2]
      d <- L[a, ] - L[b, ]; dd <- sqrt(sum(d^2)) + 1e-9
      push <- (min_d - dd) / 2 * step * d / dd
      disp[a, ] <- disp[a, ] + push; disp[b, ] <- disp[b, ] - push
    }
    L <- L + disp
  }
  L
}

g7 <- readRDS(file.path(moddir, "5_1_1_Louvain_A7_graph_max.rds"))
g9 <- readRDS(file.path(moddir, "5_1_1_Louvain_A9_graph_max.rds"))

hall_genes <- function(tsv, mod) {
  vroom(tsv, show_col_types = FALSE) |>
    mutate(HM = sub("^HALLMARK_", "", Description)) |>
    filter(Community_Name == mod, HM %in% SHARED_HM) |>
    transmute(HM, gene = strsplit(geneID, "/")) |> unnest(gene)
}
hg7 <- hall_genes(file.path(ora, "5_2_3_Hallmark7_louv_max.tsv"), "M8")
hg9 <- hall_genes(file.path(ora, "5_2_3_Hallmark9_louv_max.tsv"), "M13")
drv7 <- unique(hg7$gene); drv9 <- unique(hg9$gene)

## ---- one side: hallmark-driver sub-network only -----------------------------
make_side <- function(g_full, mod, hg, net, other_drv, xc, seed = 11,
                      Sx = 3.25, Sy = 4.75) {
  drv <- unique(hg$gene)
  sub <- induced_subgraph(g_full, V(g_full)[V(g_full)$Community_name == mod &
                                            V(g_full)$name %in% drv])
  set.seed(seed)
  L <- layout_with_fr(sub, niter = 3000)
  L <- sweep(L, 2, colMeans(L)); rad <- sqrt(rowSums(L^2)); L <- L / quantile(rad, 0.95)
  rad <- sqrt(rowSums(L^2)); L <- L * pmin(1, 1.0 / rad)
  L[, 1] <- L[, 1] * Sx; L[, 2] <- L[, 2] * Sy
  nm <- V(sub)$name
  ty <- tapply(hm_y[hg$HM], hg$gene, mean); yi <- match(nm, names(ty)); has <- !is.na(yi)
  L[has, 2] <- 0.74 * L[has, 2] + 0.26 * ty[yi[has]]
  L <- declump(L, min_d = 0.82)                      # homogeneous spacing, no overlaps
  shared <- intersect(nm, other_drv)
  ## degree within the FULL module (all genes, before dropping non-drivers)
  mdeg <- degree(induced_subgraph(g_full, V(g_full)[V(g_full)$Community_name == mod]))
  nodes <- tibble(gene = nm, x = L[, 1] + xc, y = L[, 2],
                  deg = degree(sub), mdeg = mdeg[nm], net = net) |>
    mutate(role = if_else(gene %in% shared, "Shared driver", paste0(net, " driver")))
  P <- setNames(nodes$x, nodes$gene); Q <- setNames(nodes$y, nodes$gene)
  el <- as_edgelist(sub) |> as_tibble(.name_repair = ~c("from", "to")) |>
    mutate(x = P[from], y = Q[from], xend = P[to], yend = Q[to])
  links <- hg |> filter(gene %in% nm) |>
    mutate(x = P[gene], y = Q[gene], xend = 0, yend = hm_y[HM])
  list(nodes = nodes, el = el, links = links, n_drv = length(drv), n_shared = length(shared),
       rx = max(abs(L[, 1])), ry = max(abs(L[, 2])))
}

XC <- 4.7; SX <- 3.25; SY <- 4.75
A <- make_side(g7, "M8",  hg7, "HPV-A7", drv9, xc = -XC, Sx = SX, Sy = SY)
B <- make_side(g9, "M13", hg9, "HPV-A9", drv7, xc =  XC, Sx = SX, Sy = SY)

nodes <- bind_rows(A$nodes, B$nodes)
el    <- bind_rows(A$el, B$el)
links <- bind_rows(A$links, B$links)
hubs  <- tibble(HM = SHARED_HM, x = 0, y = hm_y[SHARED_HM])

ell <- function(cx, rx, ry) { t <- seq(0, 2 * pi, length.out = 160); tibble(x = cx + rx * cos(t), y = ry * sin(t)) }
bgA <- ell(-XC, A$rx * 1.14, A$ry * 1.08); bgB <- ell(XC, B$rx * 1.14, B$ry * 1.08)
XLIM <- max(A$rx, B$rx) + XC + 0.8; YLIM <- max(A$ry, B$ry) + 0.6

## labels: TOP 10 most connected genes per module (by full-module degree)
topA <- A$nodes |> slice_max(mdeg, n = 10, with_ties = FALSE)
topB <- B$nodes |> slice_max(mdeg, n = 10, with_ties = FALSE)
lab_drvA  <- topA |> filter(role == "HPV-A7 driver")
lab_drvB  <- topB |> filter(role == "HPV-A9 driver")
lab_shared <- bind_rows(topA, topB) |> filter(role == "Shared driver")
## single label set (one repel layer -> ggrepel separates ALL of them, no cross-layer overlaps)
labs <- bind_rows(
  transmute(lab_drvA,  gene, x, y, col = col_A7,  ff = "bold"),
  transmute(lab_drvB,  gene, x, y, col = col_A9,  ff = "bold"),
  transmute(lab_shared, gene, x, y, col = "grey20", ff = "italic"))

p <- ggplot() +
  geom_polygon(data = bgA, aes(x, y), fill = col_A7, alpha = 0.05, colour = NA) +
  geom_polygon(data = bgB, aes(x, y), fill = col_A9, alpha = 0.05, colour = NA) +
  geom_segment(data = links, aes(x, y, xend = xend, yend = yend, colour = HM),
               linewidth = 0.35, alpha = 0.20) +
  geom_curve(data = el, aes(x, y, xend = xend, yend = yend),
             curvature = 0.15, colour = "grey60", linewidth = 0.45, alpha = 0.5,
             lineend = "round") +
  geom_point(data = nodes, aes(x, y, fill = role), shape = 21, size = 7.4,
             stroke = 1.1, colour = "white", alpha = 0.88) +
  geom_point(data = hubs, aes(x, y, colour = HM), shape = 18, size = 20) +
  geom_label(data = hubs, aes(x, y, label = hm_lab[HM]), colour = "grey15",
             size = 4.6, fontface = "bold", label.size = 0, lineheight = 0.9,
             fill = alpha("white", 0.9), label.padding = unit(0.16, "lines")) +
  geom_text_repel(data = labs, aes(x, y, label = gene),
                  colour = labs$col, fontface = labs$ff, size = 4.9,
                  max.overlaps = Inf, min.segment.length = 0, segment.size = 0.3,
                  segment.colour = "grey55", box.padding = 1.3, point.padding = 0.7,
                  force = 3, seed = 3, bg.color = "white", bg.r = 0.16) +
  annotate("text", x = -XC, y = 5.6, label = "HPV-A7  ·  module M8",
           fontface = "bold", size = 7.6, colour = col_A7) +
  annotate("text", x = -XC, y = 5.22, label = sprintf("%d hallmark drivers · %d shared", A$n_drv, A$n_shared),
           size = 4.6, colour = "grey40") +
  annotate("text", x =  XC, y = 5.6, label = "HPV-A9  ·  module M13",
           fontface = "bold", size = 7.6, colour = col_A9) +
  annotate("text", x =  XC, y = 5.22, label = sprintf("%d hallmark drivers · %d shared", B$n_drv, B$n_shared),
           size = 4.6, colour = "grey40") +
  annotate("text", x = 0, y = 5.45, label = "Shared cancer\nhallmarks",
           fontface = "bold", size = 5.6, colour = "grey20", lineheight = 0.9) +
  scale_fill_manual(values = c("HPV-A7 driver" = col_A7, "HPV-A9 driver" = col_A9,
                               "Shared driver" = col_shared_drv),
                    breaks = c("HPV-A7 driver", "HPV-A9 driver", "Shared driver"), name = NULL,
                    guide = guide_legend(override.aes = list(size = 6.5, stroke = 1.1, colour = "white"))) +
  scale_colour_manual(values = hm_col,
    labels = c(G2M_CHECKPOINT = "G2M checkpoint", E2F_TARGETS = "E2F targets",
               MITOTIC_SPINDLE = "Mitotic spindle"), name = "Hallmark process",
    guide = guide_legend(override.aes = list(shape = 18, size = 6, linewidth = 0))) +
  coord_cartesian(xlim = c(-XLIM, XLIM), ylim = c(-YLIM, max(YLIM + 0.3, 6.0)), clip = "off") +
  labs(title = "The same cell-cycle programme assembled from largely different genes",
       subtitle = "HPV-A7 M8 and HPV-A9 M13 both enrich G2M checkpoint, E2F targets and mitotic spindle (gene Jaccard = 0.199; 28 shared drivers, in purple)") +
  theme_void(base_size = 13) +
  theme(legend.position = "bottom", legend.box = "horizontal",
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = margin(b = 4)),
        plot.subtitle = element_text(size = 12, colour = "grey35", hjust = 0.5, margin = margin(b = 4)),
        legend.title = element_text(face = "bold", size = 11),
        legend.text = element_text(size = 10), plot.margin = margin(12, 14, 8, 14))

out_png <- file.path(ora, "5_2_9_Module_pair_cellcycle_organic.png")
out_pdf <- file.path(ora, "5_2_9_Module_pair_cellcycle_organic.pdf")
ggsave(out_png, p, width = 17, height = 13.5, dpi = 350, bg = "white", limitsize = FALSE)
ggsave(out_pdf, p, width = 17, height = 13.5, device = cairo_pdf, bg = "white", limitsize = FALSE)
message(sprintf("[OK] saved %s  (A7 %d nodes/%d edges; A9 %d nodes/%d edges; %d shared drivers)",
                out_png, nrow(A$nodes), nrow(A$el), nrow(B$nodes), nrow(B$el), A$n_shared))
