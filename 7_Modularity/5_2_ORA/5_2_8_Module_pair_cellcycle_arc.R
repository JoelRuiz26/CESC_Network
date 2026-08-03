### ============================================================
### 5_2_8_Module_pair_cellcycle_arc.R
###
### Same structured "arc" schematic as 5_2_6, for the cell-cycle example:
### HPV-A7 module M8 vs HPV-A9 module M13. Both enrich the SAME set of hallmarks
### (G2M checkpoint, E2F targets, mitotic spindle) yet are built from largely
### different genes (gene Jaccard = 0.199; 28 shared drivers, in purple).
###
###   -> 5_2_8_Module_pair_cellcycle_arc.png / .pdf
### ============================================================
suppressPackageStartupMessages({
  library(tidyverse); library(vroom); library(ggforce)
})

ora <- "~/CESC_Network/7_Modularity/5_2_ORA"
col_A7 <- "#4E79A7"; col_A9 <- "#CB6A5C"; col_shared_drv <- "#5E3C99"

SHARED_HM <- c("G2M_CHECKPOINT", "E2F_TARGETS", "MITOTIC_SPINDLE")   # top->bottom
hm_col <- c(G2M_CHECKPOINT = "#E69F00", E2F_TARGETS = "#009E73", MITOTIC_SPINDLE = "#D55E00")
hm_lab <- c(G2M_CHECKPOINT = "G2M checkpoint", E2F_TARGETS = "E2F targets",
            MITOTIC_SPINDLE = "Mitotic spindle")
hm_y   <- c(G2M_CHECKPOINT = 3.7, E2F_TARGETS = 0, MITOTIC_SPINDLE = -3.7)

hall_long <- function(tsv, mod) {
  vroom(tsv, show_col_types = FALSE) |>
    mutate(HM = sub("^HALLMARK_", "", Description)) |>
    filter(Community_Name == mod, HM %in% SHARED_HM) |>
    transmute(HM, padj = p.adjust, gene = strsplit(geneID, "/")) |> unnest(gene)
}
L7 <- hall_long(file.path(ora, "5_2_3_Hallmark7_louv_max.tsv"), "M8")
L9 <- hall_long(file.path(ora, "5_2_3_Hallmark9_louv_max.tsv"), "M13")
shared_drv <- intersect(unique(L7$gene), unique(L9$gene))

YTOP <- 5.4; YBOT <- -5.4; GAP <- 2
make_col <- function(Lg, xc) {
  prim <- Lg |> group_by(gene) |> slice_min(padj, n = 1, with_ties = FALSE) |> ungroup() |>
    mutate(band = factor(HM, levels = SHARED_HM)) |> arrange(band, padj)
  ord <- prim$gene
  jump <- c(0, cumsum(as.integer(prim$band)[-1] != as.integer(prim$band)[-length(prim$band)]))
  slot <- seq_along(ord) + jump * GAP
  y <- YTOP - (slot - 1) / (max(slot) - 1) * (YTOP - YBOT)
  nodes <- tibble(gene = ord, x = xc, y = y,
                  role = if_else(gene %in% shared_drv, "Shared driver",
                                 if_else(xc < 0, "HPV-A7 driver", "HPV-A9 driver")))
  P <- setNames(nodes$x, nodes$gene); Q <- setNames(nodes$y, nodes$gene)
  links <- Lg |> mutate(x = P[gene], y = Q[gene], xend = 0, yend = hm_y[HM],
                        id = paste(gene, HM, xc))
  list(nodes = nodes, links = links)
}
XC <- 4.6
A <- make_col(L7, -XC); B <- make_col(L9, XC)
nodes <- bind_rows(A$nodes, B$nodes); links <- bind_rows(A$links, B$links)
hubs  <- tibble(HM = SHARED_HM, x = 0, y = hm_y[SHARED_HM])

p <- ggplot() +
  ggforce::geom_diagonal(data = links,
    aes(x = x, y = y, xend = xend, yend = yend, colour = HM, group = id),
    linewidth = 0.5, alpha = 0.5, strength = 0.55, n = 80) +
  geom_point(data = nodes, aes(x, y, fill = role), shape = 21, size = 5.0,
             stroke = 0.9, colour = "white") +
  geom_text(data = filter(nodes, x < 0), aes(x, y, label = gene),
            hjust = 1, nudge_x = -0.22, size = 2.6, colour = "grey15") +
  geom_text(data = filter(nodes, x > 0), aes(x, y, label = gene),
            hjust = 0, nudge_x = 0.22, size = 2.6, colour = "grey15") +
  geom_point(data = hubs, aes(x, y, colour = HM), shape = 18, size = 20) +
  geom_label(data = hubs, aes(x, y, label = hm_lab[HM]), colour = "grey15",
             size = 4.4, fontface = "bold", label.size = 0, lineheight = 0.9,
             fill = alpha("white", 0.92), label.padding = unit(0.18, "lines")) +
  annotate("text", x = -XC, y = 6.2, label = "HPV-A7 · M8", fontface = "bold",
           size = 7, colour = col_A7) +
  annotate("text", x = -XC, y = 5.85, label = "61 hallmark drivers", size = 4.2, colour = "grey40") +
  annotate("text", x =  XC, y = 6.2, label = "HPV-A9 · M13", fontface = "bold",
           size = 7, colour = col_A9) +
  annotate("text", x =  XC, y = 5.85, label = "41 hallmark drivers", size = 4.2, colour = "grey40") +
  scale_colour_manual(values = hm_col,
    labels = c(G2M_CHECKPOINT = "G2M checkpoint", E2F_TARGETS = "E2F targets",
               MITOTIC_SPINDLE = "Mitotic spindle"), name = "Hallmark process",
    guide = guide_legend(order = 2, override.aes = list(shape = 18, size = 6, linewidth = 0))) +
  scale_fill_manual(values = c("HPV-A7 driver" = col_A7, "HPV-A9 driver" = col_A9,
                               "Shared driver" = col_shared_drv),
                    breaks = c("HPV-A7 driver", "HPV-A9 driver", "Shared driver"), name = NULL,
                    guide = guide_legend(order = 1, override.aes = list(size = 6, stroke = 0.9, colour = "white"))) +
  coord_cartesian(xlim = c(-(XC + 2.6), XC + 2.6), ylim = c(-5.9, 6.5), clip = "off") +
  labs(title = "The same cell-cycle programme assembled from largely different genes",
       subtitle = "HPV-A7 M8 and HPV-A9 M13 enrich the same hallmark set (G2M checkpoint, E2F targets, mitotic spindle) through largely different genes (gene Jaccard = 0.199; 28 shared drivers)") +
  theme_void(base_size = 13) +
  theme(legend.position = "bottom", legend.box = "horizontal",
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = margin(b = 4)),
        plot.subtitle = element_text(size = 12, colour = "grey35", hjust = 0.5, margin = margin(b = 6)),
        legend.title = element_text(face = "bold", size = 11),
        legend.text = element_text(size = 10), plot.margin = margin(12, 16, 8, 16))

out_png <- file.path(ora, "5_2_8_Module_pair_cellcycle_arc.png")
out_pdf <- file.path(ora, "5_2_8_Module_pair_cellcycle_arc.pdf")
ggsave(out_png, p, width = 16, height = 13, dpi = 350, bg = "white", limitsize = FALSE)
ggsave(out_pdf, p, width = 16, height = 13, device = cairo_pdf, bg = "white", limitsize = FALSE)
message(sprintf("[OK] saved %s  (A7 %d genes, A9 %d genes, %d shared drivers)",
                out_png, nrow(A$nodes), nrow(B$nodes), length(shared_drv)))
