### ============================================================
### 5_2_7_Module_pair_organic_network.R
###
### Organic-network variant of the same-programme/different-genes example
### (companion to the tidy arc diagram 5_2_6). Each network is a force-directed
### cloud on its side (A7 left / A9 right); the 3 shared hallmarks sit in the
### central column and the driver genes link inward. Gentle grouping pulls each
### driver toward its hallmark band. Only a FEW representative genes are labelled
### (the shared drivers + the best-connected drivers per side).
###
###   -> 5_2_7_Module_pair_organic_network.png / .pdf
### ============================================================
suppressPackageStartupMessages({
  library(igraph); library(tidyverse); library(vroom); library(ggrepel)
})

ora    <- "~/CESC_Network/7_Modularity/5_2_ORA"
moddir <- "~/CESC_Network/7_Modularity/5_1_Modules"

col_A7 <- "#4E79A7"; col_A9 <- "#CB6A5C"
col_shared <- "#8C79B5"; col_shared_drv <- "#5E3C99"; grey_node <- "grey78"

SHARED_HM <- c("MYC_TARGETS_V1", "MYC_TARGETS_V2", "OXIDATIVE_PHOSPHORYLATION")
hm_col <- c(MYC_TARGETS_V1 = "#E69F00", MYC_TARGETS_V2 = "#D55E00",
           OXIDATIVE_PHOSPHORYLATION = "#009E73")
hm_lab <- c(MYC_TARGETS_V1 = "MYC targets V1", MYC_TARGETS_V2 = "MYC targets V2",
            OXIDATIVE_PHOSPHORYLATION = "Oxidative\nphosphorylation")
hm_y   <- c(MYC_TARGETS_V1 = 4.4, OXIDATIVE_PHOSPHORYLATION = 0, MYC_TARGETS_V2 = -4.4)

g7 <- readRDS(file.path(moddir, "5_1_1_Louvain_A7_graph_max.rds"))
g9 <- readRDS(file.path(moddir, "5_1_1_Louvain_A9_graph_max.rds"))

hall_genes <- function(tsv, mod) {
  vroom(tsv, show_col_types = FALSE) |>
    mutate(HM = sub("^HALLMARK_", "", Description)) |>
    filter(Community_Name == mod, HM %in% SHARED_HM) |>
    transmute(HM, gene = strsplit(geneID, "/")) |> unnest(gene)
}
hg7 <- hall_genes(file.path(ora, "5_2_3_Hallmark7_louv_max.tsv"), "M7")
hg9 <- hall_genes(file.path(ora, "5_2_3_Hallmark9_louv_max.tsv"), "M8")

## ---- build one side: functional sub-network (drivers + shared) + grey context
make_side <- function(g_full, mod, hg, net, other_genes, xc, seed = 11,
                      Sx = 3.1, Sy = 4.55) {
  drv    <- unique(hg$gene)
  own    <- V(g_full)$name[V(g_full)$Community_name == mod]
  shared <- intersect(own, other_genes)
  keep   <- union(drv, shared)
  sub <- induced_subgraph(g_full, V(g_full)[V(g_full)$Community_name == mod &
                                            V(g_full)$name %in% keep])
  set.seed(seed)
  L <- layout_with_fr(sub, niter = 3000)
  L <- sweep(L, 2, colMeans(L)); rad <- sqrt(rowSums(L^2)); L <- L / quantile(rad, 0.95)
  rad <- sqrt(rowSums(L^2)); L <- L * pmin(1, 1.0 / rad)
  L[, 1] <- L[, 1] * Sx; L[, 2] <- L[, 2] * Sy
  nm <- V(sub)$name
  ## gentle grouping toward each driver's hallmark band
  ty <- tapply(hm_y[hg$HM], hg$gene, mean); yi <- match(nm, names(ty)); has <- !is.na(yi)
  L[has, 2] <- 0.68 * L[has, 2] + 0.32 * ty[yi[has]]
  nrr <- sqrt((L[, 1] / Sx)^2 + (L[, 2] / Sy)^2); s <- pmin(1, 1 / nrr)
  L[, 1] <- L[, 1] * s; L[, 2] <- L[, 2] * s
  nodes <- tibble(gene = nm, x = L[, 1] + xc, y = L[, 2], deg = degree(sub), net = net) |>
    mutate(role = case_when(
      gene %in% shared & gene %in% drv ~ "Shared driver",
      gene %in% shared                 ~ "Shared gene",
      TRUE                             ~ paste0(net, " driver")))
  P <- setNames(nodes$x, nodes$gene); Q <- setNames(nodes$y, nodes$gene)
  el <- as_edgelist(sub) |> as_tibble(.name_repair = ~c("from", "to")) |>
    mutate(x = P[from], y = Q[from], xend = P[to], yend = Q[to])
  links <- hg |> filter(gene %in% nm) |>
    mutate(x = P[gene], y = Q[gene], xend = 0, yend = hm_y[HM])
  ## grey context: all OTHER module genes near their subnet neighbours, bounded by ellipse
  subfull <- induced_subgraph(g_full, V(g_full)[V(g_full)$Community_name == mod])
  others  <- setdiff(own, keep); set.seed(seed + 1)
  opos <- t(vapply(others, function(gn) {
    fnb <- intersect(names(neighbors(subfull, gn)), keep)
    if (length(fnb) > 0) c(mean(P[fnb]) + rnorm(1, 0, 0.40), mean(Q[fnb]) + rnorm(1, 0, 0.60))
    else { ang <- runif(1, 0, 2 * pi); rr <- sqrt(runif(1, 0.35, 1)); c(xc + cos(ang) * Sx * rr, sin(ang) * Sy * rr) }
  }, numeric(2)))
  enx <- Sx * 1.05; eny <- Sy * 1.05
  nr <- sqrt(((opos[, 1] - xc) / enx)^2 + (opos[, 2] / eny)^2); sc <- pmin(1, 1 / nr)
  opos[, 1] <- xc + (opos[, 1] - xc) * sc; opos[, 2] <- opos[, 2] * sc
  others_nodes <- tibble(gene = others, x = opos[, 1], y = opos[, 2],
                         deg = NA_integer_, net = net, role = "Other module gene")
  ce <- lapply(others, function(gn) {
    fnb <- intersect(names(neighbors(subfull, gn)), keep); if (length(fnb)) tibble(gene = gn, to = fnb) else NULL
  }) |> bind_rows()
  Ox <- setNames(others_nodes$x, others_nodes$gene); Oy <- setNames(others_nodes$y, others_nodes$gene)
  ctx_edges <- if (nrow(ce)) ce |> mutate(x = Ox[gene], y = Oy[gene], xend = P[to], yend = Q[to]) else ce
  list(nodes = nodes, el = el, links = links, others = others_nodes, ctx_edges = ctx_edges,
       n_drv = length(drv), n_shared = length(shared))
}

XC <- 4.1; SX <- 3.1; SY <- 4.55
A <- make_side(g7, "M7", hg7, "HPV-A7", V(g9)$name[V(g9)$Community_name == "M8"], xc = -XC, Sx = SX, Sy = SY)
B <- make_side(g9, "M8", hg9, "HPV-A9", V(g7)$name[V(g7)$Community_name == "M7"], xc =  XC, Sx = SX, Sy = SY)

nodes     <- bind_rows(A$nodes, B$nodes)
el        <- bind_rows(A$el, B$el)
links     <- bind_rows(A$links, B$links)
others    <- bind_rows(A$others, B$others)
ctx_edges <- bind_rows(A$ctx_edges, B$ctx_edges)
hubs      <- tibble(HM = SHARED_HM, x = 0, y = hm_y[SHARED_HM])

## soft clade-tinted background ellipse per side (delimits each network)
ell <- function(cx, rx, ry) { t <- seq(0, 2 * pi, length.out = 160); tibble(x = cx + rx * cos(t), y = ry * sin(t)) }
bgA <- ell(-XC, SX * 1.14, SY * 1.07); bgB <- ell(XC, SX * 1.14, SY * 1.07)

big_roles  <- c("HPV-A7 driver", "HPV-A9 driver", "Shared driver")
## label only a FEW representative genes: shared drivers + best-connected drivers per side
lab_shared <- nodes |> filter(role == "Shared driver")
lab_drvA <- A$nodes |> filter(role == "HPV-A7 driver") |> slice_max(deg, n = 5)
lab_drvB <- B$nodes |> filter(role == "HPV-A9 driver") |> slice_max(deg, n = 5)

## ---- plot ------------------------------------------------------------------
p <- ggplot() +
  ## soft background per side
  geom_polygon(data = bgA, aes(x, y), fill = col_A7, alpha = 0.05, colour = NA) +
  geom_polygon(data = bgB, aes(x, y), fill = col_A9, alpha = 0.05, colour = NA) +
  ## faint module context (dots only, no edges -> much cleaner)
  geom_point(data = others, aes(x, y, fill = role),
             shape = 21, size = 1.8, stroke = 0, colour = "white", alpha = 0.55) +
  ## gene -> hallmark links (soft, underneath)
  geom_segment(data = links, aes(x, y, xend = xend, yend = yend, colour = HM),
               linewidth = 0.4, alpha = 0.25) +
  ## co-expression edges, curved for a cleaner network look
  geom_curve(data = el, aes(x, y, xend = xend, yend = yend),
             curvature = 0.16, colour = "grey48", linewidth = 0.5, alpha = 0.55,
             lineend = "round") +
  geom_point(data = filter(nodes, role == "Shared gene"), aes(x, y, fill = role),
             shape = 21, size = 5.6, stroke = 0.6, colour = "white") +
  geom_point(data = filter(nodes, role %in% big_roles), aes(x, y, fill = role),
             shape = 21, size = 7.2, stroke = 1.1, colour = "white") +
  geom_point(data = hubs, aes(x, y, colour = HM), shape = 18, size = 19) +
  geom_label(data = hubs, aes(x, y, label = hm_lab[HM]), colour = "grey15",
             size = 4.4, fontface = "bold", label.size = 0, lineheight = 0.9,
             fill = alpha("white", 0.9), label.padding = unit(0.16, "lines")) +
  geom_text_repel(data = lab_drvA, aes(x, y, label = gene), size = 5.0, fontface = "bold",
                  colour = col_A7, max.overlaps = 50, min.segment.length = 0,
                  segment.size = 0.3, segment.colour = "grey55", box.padding = 0.9,
                  point.padding = 0.5, seed = 3) +
  geom_text_repel(data = lab_drvB, aes(x, y, label = gene), size = 5.0, fontface = "bold",
                  colour = col_A9, max.overlaps = 50, min.segment.length = 0,
                  segment.size = 0.3, segment.colour = "grey55", box.padding = 0.9,
                  point.padding = 0.5, seed = 3) +
  geom_text_repel(data = lab_shared, aes(x, y, label = gene), size = 4.7, fontface = "italic",
                  colour = "grey20", max.overlaps = 50, min.segment.length = 0,
                  segment.size = 0.3, segment.colour = "grey60", box.padding = 0.8,
                  point.padding = 0.5, seed = 3) +
  annotate("text", x = -XC, y = 5.5, label = "HPV-A7  ·  module M7",
           fontface = "bold", size = 7.6, colour = col_A7) +
  annotate("text", x = -XC, y = 5.13, label = sprintf("%d hallmark drivers · %d shared genes", A$n_drv, A$n_shared),
           size = 4.6, colour = "grey40") +
  annotate("text", x =  XC, y = 5.5, label = "HPV-A9  ·  module M8",
           fontface = "bold", size = 7.6, colour = col_A9) +
  annotate("text", x =  XC, y = 5.13, label = sprintf("%d hallmark drivers · %d shared genes", B$n_drv, B$n_shared),
           size = 4.6, colour = "grey40") +
  annotate("text", x = 0, y = 5.35, label = "Shared cancer\nhallmarks",
           fontface = "bold", size = 5.6, colour = "grey20", lineheight = 0.9) +
  scale_fill_manual(
    values = c("HPV-A7 driver" = col_A7, "HPV-A9 driver" = col_A9,
               "Shared gene" = col_shared, "Shared driver" = col_shared_drv,
               "Other module gene" = grey_node),
    breaks = c("HPV-A7 driver", "HPV-A9 driver", "Shared gene", "Shared driver", "Other module gene"),
    name = NULL,
    guide = guide_legend(nrow = 1, override.aes = list(
      size = c(6, 6, 6, 6, 2.4), stroke = c(1.1, 1.1, 0.6, 1.1, 0), colour = "white"))) +
  scale_colour_manual(values = hm_col, labels = c("MYC targets V1", "MYC targets V2",
                      "Oxidative phosphorylation"), name = "Hallmark process",
                      guide = guide_legend(override.aes = list(shape = 18, size = 6, linewidth = 0))) +
  coord_cartesian(xlim = c(-(XC + 3.1), XC + 3.1), ylim = c(-5.1, 5.85), clip = "off") +
  labs(title = "A shared oncogenic programme assembled from largely different genes",
       subtitle = "HPV-A7 M7 and HPV-A9 M8 both feed MYC targets and oxidative phosphorylation, yet share only 30 genes (gene Jaccard = 0.024; 6 shared drivers)") +
  theme_void(base_size = 13) +
  theme(legend.position = "bottom", legend.box = "vertical",
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = margin(b = 4)),
        plot.subtitle = element_text(size = 12, colour = "grey35", hjust = 0.5, margin = margin(b = 4)),
        legend.title = element_text(face = "bold", size = 11),
        legend.text = element_text(size = 10),
        plot.margin = margin(12, 14, 8, 14))

out_png <- file.path(ora, "5_2_7_Module_pair_organic_network.png")
out_pdf <- file.path(ora, "5_2_7_Module_pair_organic_network.pdf")
ggsave(out_png, p, width = 17, height = 13.5, dpi = 350, bg = "white", limitsize = FALSE)
ggsave(out_pdf, p, width = 17, height = 13.5, device = cairo_pdf, bg = "white", limitsize = FALSE)
message(sprintf("[OK] saved %s  (A7 %d nodes/%d edges; A9 %d nodes/%d edges)",
                out_png, nrow(A$nodes), nrow(A$el), nrow(B$nodes), nrow(B$el)))
