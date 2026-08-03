### ============================================================
### 5_2_5_Module_pair_MYC_OXPHOS_network.R
###
### Two-panel network figure illustrating the "same programme, different genes"
### example: HPV-A7 module M7 vs HPV-A9 module M8. Both modules are strongly &
### almost identically enriched for MYC TARGETS V1/V2 and OXIDATIVE
### PHOSPHORYLATION, yet their gene composition is nearly disjoint
### (gene Jaccard = 0.024; 30 shared genes, only 6 of them programme drivers).
###
### Per panel:
###   * every module gene = tiny light-grey node, co-expression edges faint grey
###   * hallmark-driver genes (leading-edge of the 3 shared hallmarks) = large,
###     coloured by clade; their mutual co-expression edges drawn darker
###   * genes SHARED between the two modules = purple (ring = shared AND driver)
###   * 3 hallmark "process" nodes on the periphery, linked to the genes that
###     enriched them (link colour = hallmark)
###
###   -> 5_2_5_Module_pair_MYC_OXPHOS_network.png / .pdf
### ============================================================
suppressPackageStartupMessages({
  library(igraph); library(tidyverse); library(vroom)
  library(ggrepel); library(patchwork)
})

ora   <- "~/CESC_Network/7_Modularity/5_2_ORA"
moddir<- "~/CESC_Network/7_Modularity/5_1_Modules"

col_A7 <- "#4E79A7"; col_A9 <- "#CB6A5C"
col_shared <- "#8C79B5"; col_shared_drv <- "#5E3C99"
grey_node  <- "grey80"

## shared hallmarks to illustrate + their colours (Okabe-Ito, clade/purple-safe)
SHARED_HM <- c("MYC_TARGETS_V1", "MYC_TARGETS_V2", "OXIDATIVE_PHOSPHORYLATION")
hm_col <- c(MYC_TARGETS_V1 = "#E69F00", MYC_TARGETS_V2 = "#D55E00",
           OXIDATIVE_PHOSPHORYLATION = "#009E73")
hm_lab <- c(MYC_TARGETS_V1 = "MYC targets V1", MYC_TARGETS_V2 = "MYC targets V2",
            OXIDATIVE_PHOSPHORYLATION = "Oxidative phosphorylation")
hm_ang <- c(MYC_TARGETS_V1 = 90, MYC_TARGETS_V2 = 210, OXIDATIVE_PHOSPHORYLATION = 330)

g7 <- readRDS(file.path(moddir, "5_1_1_Louvain_A7_graph_max.rds"))
g9 <- readRDS(file.path(moddir, "5_1_1_Louvain_A9_graph_max.rds"))

hall_genes <- function(tsv, mod) {
  vroom(tsv, show_col_types = FALSE) |>
    mutate(HM = sub("^HALLMARK_", "", Description)) |>
    filter(Community_Name == mod, HM %in% SHARED_HM) |>
    transmute(HM, gene = strsplit(geneID, "/")) |>
    unnest(gene)
}
hg7 <- hall_genes(file.path(ora, "5_2_3_Hallmark7_louv_max.tsv"), "M7")
hg9 <- hall_genes(file.path(ora, "5_2_3_Hallmark9_louv_max.tsv"), "M8")

## ---- panel builder --------------------------------------------------------
build_panel <- function(g_full, mod, hg, accent, other_genes, seed = 7) {
  sub <- induced_subgraph(g_full, V(g_full)[V(g_full)$Community_name == mod])
  set.seed(seed)
  L  <- layout_with_fr(sub)
  ct <- colMeans(L); L <- sweep(L, 2, ct)                    # centre
  rad <- sqrt(rowSums(L^2)); L <- L / quantile(rad, 0.95)    # scale
  nm <- V(sub)$name
  nodes <- tibble(gene = nm, x = L[, 1], y = L[, 2], deg = degree(sub))

  drv    <- unique(hg$gene)
  shared <- intersect(nm, other_genes)
  nodes <- nodes |>
    mutate(is_drv = gene %in% drv, is_shared = gene %in% shared,
           role = case_when(
             is_shared &  is_drv ~ "Shared driver",
             is_shared & !is_drv ~ "Shared gene",
             !is_shared & is_drv ~ "Hallmark driver",
             TRUE                ~ "Other gene"))

  ## co-expression edges
  el <- as_edgelist(sub) |> as_tibble(.name_repair = ~c("from", "to"))
  P  <- setNames(nodes$x, nodes$gene); Q <- setNames(nodes$y, nodes$gene)
  hot <- union(drv, shared)
  el <- el |> mutate(x = P[from], y = Q[from], xend = P[to], yend = Q[to],
                     core = from %in% hot & to %in% hot)

  ## hallmark hub positions + membership links
  hubs <- tibble(HM = SHARED_HM, ang = hm_ang[SHARED_HM] * pi / 180) |>
    mutate(x = 1.55 * cos(ang), y = 1.55 * sin(ang))
  Hx <- setNames(hubs$x, hubs$HM); Hy <- setNames(hubs$y, hubs$HM)
  links <- hg |> filter(gene %in% nodes$gene) |>
    mutate(x = P[gene], y = Q[gene], xend = Hx[HM], yend = Hy[HM])

  lab_df  <- nodes |> filter(role == "Shared driver")
  lab_drv <- nodes |> filter(role == "Hallmark driver") |> slice_max(deg, n = 6)

  ggplot() +
    ## edges
    geom_segment(data = el, aes(x, y, xend = xend, yend = yend),
                 colour = "grey78", linewidth = 0.10, alpha = 0.35) +
    geom_segment(data = filter(el, core), aes(x, y, xend = xend, yend = yend),
                 colour = "grey45", linewidth = 0.28, alpha = 0.55) +
    ## gene -> hallmark membership links
    geom_segment(data = links, aes(x, y, xend = xend, yend = yend, colour = HM),
                 linewidth = 0.30, alpha = 0.30) +
    ## nodes, drawn small -> large
    geom_point(data = filter(nodes, role == "Other gene"), aes(x, y, fill = role),
               shape = 21, size = 0.6, stroke = 0, colour = "white") +
    geom_point(data = filter(nodes, role == "Shared gene"), aes(x, y, fill = role),
               shape = 21, size = 2.3, stroke = 0.3, colour = "white") +
    geom_point(data = filter(nodes, role == "Hallmark driver"), aes(x, y, fill = role),
               shape = 21, size = 3.1, stroke = 0.7, colour = "white") +
    geom_point(data = filter(nodes, role == "Shared driver"), aes(x, y, fill = role),
               shape = 21, size = 3.1, stroke = 0.7, colour = "white") +
    ## hallmark hubs + labels
    geom_point(data = hubs, aes(x, y, colour = HM), shape = 18, size = 10) +
    geom_label(data = hubs, aes(x, y, label = hm_lab[HM]), colour = "grey15",
               size = 3.2, fontface = "bold", label.size = 0,
               fill = alpha("white", 0.85), label.padding = unit(0.12, "lines"),
               vjust = ifelse(hubs$y >= 0, -0.9, 1.9)) +
    ## label a few clade-specific (non-shared) hallmark drivers
    geom_text_repel(data = lab_drv, aes(x, y, label = gene), size = 2.7,
                    fontface = "bold", colour = accent, max.overlaps = 30,
                    min.segment.length = 0, segment.size = 0.2,
                    segment.colour = "grey55", box.padding = 0.5, seed = 3) +
    ## label the shared driver genes
    geom_text_repel(data = lab_df, aes(x, y, label = gene), size = 2.7,
                    fontface = "italic", colour = "grey20", max.overlaps = 20,
                    min.segment.length = 0, segment.size = 0.2,
                    segment.colour = "grey60", box.padding = 0.4, seed = 3) +
    scale_fill_manual(
      values = c("Other gene" = grey_node, "Hallmark driver" = accent,
                 "Shared gene" = col_shared, "Shared driver" = col_shared_drv),
      breaks = c("Hallmark driver", "Shared gene", "Shared driver", "Other gene"),
      name = NULL,
      guide = guide_legend(override.aes = list(
        size = c(3.4, 2.6, 3.4, 1.4), stroke = c(0.7, 0.3, 0.7, 0),
        colour = c("white", "white", "white", "white")))) +
    scale_colour_manual(values = hm_col, labels = hm_lab, name = "Hallmark process",
                        guide = guide_legend(override.aes = list(
                          shape = 18, size = 5, linewidth = 0))) +
    coord_equal(clip = "off") +
    theme_void(base_size = 13)
}

genes9_all <- V(g9)$name[V(g9)$Community_name == "M8"]
genes7_all <- V(g7)$name[V(g7)$Community_name == "M7"]

pA <- build_panel(g7, "M7", hg7, col_A7, genes9_all) +
  ggtitle("HPV-A7  ·  module M7", "867 genes  ·  70 hallmark drivers")
pB <- build_panel(g9, "M8", hg9, col_A9, genes7_all) +
  ggtitle("HPV-A9  ·  module M8", "395 genes  ·  44 hallmark drivers")

fig <- (pA | pB) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "A shared oncogenic programme assembled from largely different genes",
    subtitle = expression("Modules HPV-A7 M7 and HPV-A9 M8 are both enriched for MYC targets and oxidative phosphorylation, yet share only 30 genes (gene Jaccard = 0.024; 6 shared drivers)"),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 18, hjust = 0),
      plot.subtitle = element_text(size = 12, colour = "grey30", hjust = 0),
      plot.margin   = margin(12, 12, 12, 12))) &
  theme(legend.position = "bottom", legend.box = "vertical",
        plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
        plot.subtitle = element_text(size = 11, colour = "grey35", hjust = 0.5),
        legend.title = element_text(face = "bold", size = 11),
        legend.text  = element_text(size = 10))

out_png <- file.path(ora, "5_2_5_Module_pair_MYC_OXPHOS_network.png")
out_pdf <- file.path(ora, "5_2_5_Module_pair_MYC_OXPHOS_network.pdf")
ggsave(out_png, fig, width = 15, height = 9, dpi = 350, bg = "white", limitsize = FALSE)
ggsave(out_pdf, fig, width = 15, height = 9, device = cairo_pdf, bg = "white", limitsize = FALSE)
message("[OK] saved ", out_png)
