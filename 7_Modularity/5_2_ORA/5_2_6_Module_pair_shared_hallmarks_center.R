### ============================================================
### 5_2_6_Module_pair_shared_hallmarks_center.R
###
### Structured "arc" (bipartite) schematic of the same-programme/different-genes
### example (HPV-A7 M7 vs HPV-A9 M8). The 3 shared hallmarks sit in a central
### column; each network's hallmark-driver genes form a tidy ordered column on
### its side (grouped top->bottom: MYC targets V1 / oxidative phosphorylation /
### MYC targets V2), linked to their hallmark by curved ribbons. Genes shared as
### drivers by BOTH modules are drawn in purple and joined across the middle.
###
###   -> 5_2_6_Module_pair_shared_hallmarks_center.png / .pdf
### ============================================================
suppressPackageStartupMessages({
  library(tidyverse); library(vroom); library(ggforce); library(ggrepel)
})

ora <- "~/CESC_Network/7_Modularity/5_2_ORA"
col_A7 <- "#4E79A7"; col_A9 <- "#CB6A5C"; col_shared_drv <- "#5E3C99"

SHARED_HM <- c("MYC_TARGETS_V1", "OXIDATIVE_PHOSPHORYLATION", "MYC_TARGETS_V2")  # top->bottom
hm_col <- c(MYC_TARGETS_V1 = "#E69F00", MYC_TARGETS_V2 = "#D55E00",
           OXIDATIVE_PHOSPHORYLATION = "#009E73")
hm_lab <- c(MYC_TARGETS_V1 = "MYC targets V1", MYC_TARGETS_V2 = "MYC targets V2",
            OXIDATIVE_PHOSPHORYLATION = "Oxidative\nphosphorylation")
hm_y   <- c(MYC_TARGETS_V1 = 3.7, OXIDATIVE_PHOSPHORYLATION = 0, MYC_TARGETS_V2 = -3.7)

hall_long <- function(tsv, mod) {
  vroom(tsv, show_col_types = FALSE) |>
    mutate(HM = sub("^HALLMARK_", "", Description)) |>
    filter(Community_Name == mod, HM %in% SHARED_HM) |>
    transmute(HM, padj = p.adjust, gene = strsplit(geneID, "/")) |> unnest(gene)
}
L7 <- hall_long(file.path(ora, "5_2_3_Hallmark7_louv_max.tsv"), "M7")
L9 <- hall_long(file.path(ora, "5_2_3_Hallmark9_louv_max.tsv"), "M8")
shared_drv <- intersect(unique(L7$gene), unique(L9$gene))          # driver in BOTH modules

## ---- one ordered gene column ------------------------------------------------
YTOP <- 5.4; YBOT <- -5.4; GAP <- 2
make_col <- function(Lg, xc) {
  prim <- Lg |> group_by(gene) |> slice_min(padj, n = 1, with_ties = FALSE) |> ungroup() |>
    mutate(band = factor(HM, levels = SHARED_HM)) |>
    arrange(band, padj)
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
nodes <- bind_rows(A$nodes, B$nodes)
links <- bind_rows(A$links, B$links)
hubs  <- tibble(HM = SHARED_HM, x = 0, y = hm_y[SHARED_HM])

## connectors joining the SAME shared-driver gene across the two columns
sd <- nodes |> filter(role == "Shared driver")
sd_link <- sd |> filter(x < 0) |> select(gene, x1 = x, y1 = y) |>
  inner_join(sd |> filter(x > 0) |> select(gene, x2 = x, y2 = y), by = "gene")

## ---- plot ------------------------------------------------------------------
p <- ggplot() +
  ## gene -> hallmark curved ribbons
  ggforce::geom_diagonal(data = links,
    aes(x = x, y = y, xend = xend, yend = yend, colour = HM, group = id),
    linewidth = 0.5, alpha = 0.5, strength = 0.55, n = 80) +
  ## gene nodes
  geom_point(data = nodes, aes(x, y, fill = role), shape = 21, size = 5.4,
             stroke = 0.9, colour = "white") +
  ## gene labels (outer side of each column)
  geom_text(data = filter(nodes, x < 0), aes(x, y, label = gene),
            hjust = 1, nudge_x = -0.22, size = 2.75, colour = "grey15") +
  geom_text(data = filter(nodes, x > 0), aes(x, y, label = gene),
            hjust = 0, nudge_x = 0.22, size = 2.75, colour = "grey15") +
  ## hallmark hubs + labels
  geom_point(data = hubs, aes(x, y, colour = HM), shape = 18, size = 20) +
  geom_label(data = hubs, aes(x, y, label = hm_lab[HM]), colour = "grey15",
             size = 4.4, fontface = "bold", label.size = 0, lineheight = 0.9,
             fill = alpha("white", 0.92), label.padding = unit(0.18, "lines")) +
  ## headers
  annotate("text", x = -XC, y = 6.2, label = "HPV-A7 · M7", fontface = "bold",
           size = 7, colour = col_A7) +
  annotate("text", x = -XC, y = 5.85, label = "59 hallmark drivers", size = 4.2, colour = "grey40") +
  annotate("text", x =  XC, y = 6.2, label = "HPV-A9 · M8", fontface = "bold",
           size = 7, colour = col_A9) +
  annotate("text", x =  XC, y = 5.85, label = "41 hallmark drivers", size = 4.2, colour = "grey40") +
  scale_colour_manual(values = hm_col, labels = c("MYC targets V1", "MYC targets V2",
                      "Oxidative phosphorylation"), name = "Hallmark process",
                      guide = guide_legend(order = 2, override.aes = list(shape = 18, size = 6, linewidth = 0))) +
  scale_fill_manual(values = c("HPV-A7 driver" = col_A7, "HPV-A9 driver" = col_A9,
                               "Shared driver" = col_shared_drv),
                    breaks = c("HPV-A7 driver", "HPV-A9 driver", "Shared driver"), name = NULL,
                    guide = guide_legend(order = 1, override.aes = list(size = 6, stroke = 0.9, colour = "white"))) +
  coord_cartesian(xlim = c(-(XC + 2.6), XC + 2.6), ylim = c(-5.9, 6.5), clip = "off") +
  labs(title = "A shared oncogenic programme assembled from largely different genes",
       subtitle = "HPV-A7 M7 and HPV-A9 M8 drive the same hallmarks (MYC targets, oxidative phosphorylation) through almost disjoint gene sets (gene Jaccard = 0.024; 6 shared drivers, in purple)") +
  theme_void(base_size = 13) +
  theme(legend.position = "bottom", legend.box = "horizontal",
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = margin(b = 4)),
        plot.subtitle = element_text(size = 12, colour = "grey35", hjust = 0.5, margin = margin(b = 6)),
        legend.title = element_text(face = "bold", size = 11),
        legend.text = element_text(size = 10),
        plot.margin = margin(12, 16, 8, 16))

out_png <- file.path(ora, "5_2_6_Module_pair_shared_hallmarks_center.png")
out_pdf <- file.path(ora, "5_2_6_Module_pair_shared_hallmarks_center.pdf")
ggsave(out_png, p, width = 15, height = 13, dpi = 350, bg = "white", limitsize = FALSE)
ggsave(out_pdf, p, width = 15, height = 13, device = cairo_pdf, bg = "white", limitsize = FALSE)
message(sprintf("[OK] saved %s  (A7 %d genes, A9 %d genes, %d shared drivers)",
                out_png, nrow(A$nodes), nrow(B$nodes), length(shared_drv)))
