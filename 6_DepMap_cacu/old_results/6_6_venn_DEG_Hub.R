# ============================================================
# 6_6_venn_DEG_Hub.R
# Venn diagrams: DEGs x Hub genes (intra-grupo) +
#                genes de interes A7 vs A9 (inter-grupo)
# ============================================================

setwd("~/CESC_Network/6_DepMap_cacu/")

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
  library(scales)
})

# ── Counts ────────────────────────────────────────────────────

t7 <- read_tsv("6_3_gene_table_A7.tsv", show_col_types = FALSE)
t9 <- read_tsv("6_3_gene_table_A9.tsv", show_col_types = FALSE)

a7_int <- t7$Gene[t7$is_interest]
a9_int <- t9$Gene[t9$is_interest]

n <- list(
  A7 = list(
    deg_only = sum( t7$is_DEG & !t7$is_hub),
    hub_only = sum(!t7$is_DEG &  t7$is_hub),
    int      = sum( t7$is_interest)
  ),
  A9 = list(
    deg_only = sum( t9$is_DEG & !t9$is_hub),
    hub_only = sum(!t9$is_DEG &  t9$is_hub),
    int      = sum( t9$is_interest)
  ),
  inter = list(
    shared  = length(intersect(a7_int, a9_int)),
    excl_a7 = length(setdiff(a7_int, a9_int)),
    excl_a9 = length(setdiff(a9_int, a7_int))
  )
)

# ── Colors ────────────────────────────────────────────────────

COL_DEG <- "#6C3483"
COL_HUB <- "#D35400"
COL_INT <- "#1C2833"
COL_A7  <- "#2E86AB"
COL_A9  <- "#D35400"
COL_SH  <- "#27AE60"

# ── Helper: circle polygon ────────────────────────────────────

circle_poly <- function(cx, cy, r, n = 300) {
  theta <- seq(0, 2 * pi, length.out = n + 1)
  data.frame(x = cx + r * cos(theta), y = cy + r * sin(theta))
}

# ── Intra-grupo: DEG x Hub por clado ─────────────────────────

make_venn <- function(deg_only, hub_only, int_n, cohort) {

  r    <- 1.50
  cx_l <- -0.82
  cx_r <-  0.82

  cl <- circle_poly(cx_l, 0, r)
  cr <- circle_poly(cx_r, 0, r)

  ggplot() +
    geom_polygon(data = cl, aes(x, y),
                 fill = COL_DEG, color = NA, alpha = 0.18) +
    geom_polygon(data = cr, aes(x, y),
                 fill = COL_HUB, color = NA, alpha = 0.18) +
    geom_path(data = cl, aes(x, y),
              color = COL_DEG, linewidth = 1.6) +
    geom_path(data = cr, aes(x, y),
              color = COL_HUB, linewidth = 1.6) +
    annotate("text",
             x = cx_l - 0.55, y = 1.72,
             label = "DEGs", hjust = 0.5,
             size = 6.5, fontface = "bold", color = COL_DEG) +
    annotate("text",
             x = cx_r + 0.55, y = 1.72,
             label = "Hub genes", hjust = 0.5,
             size = 6.5, fontface = "bold", color = COL_HUB) +
    annotate("text",
             x = cx_l - 0.62, y = 0,
             label = format(deg_only, big.mark = ","),
             size = 10, fontface = "bold", color = COL_DEG, hjust = 0.5) +
    annotate("text",
             x = 0, y = 0,
             label = format(int_n, big.mark = ","),
             size = 9, fontface = "bold", color = COL_INT, hjust = 0.5) +
    annotate("text",
             x = cx_r + 0.62, y = 0,
             label = format(hub_only, big.mark = ","),
             size = 10, fontface = "bold", color = COL_HUB, hjust = 0.5) +
    annotate("text",
             x = 0, y = 2.25,
             label = paste0("HPV-", cohort),
             size = 8, fontface = "bold", color = "#1C2833", hjust = 0.5) +
    coord_fixed(xlim = c(-3.0, 3.0), ylim = c(-2.1, 2.55),
                expand = FALSE, clip = "off") +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA),
          plot.margin     = margin(8, 8, 8, 8))
}

# ── Inter-grupo: genes de interes A7 vs A9 ───────────────────

make_venn_inter <- function(excl_a7, excl_a9, shared_n) {

  r    <- 1.50
  cx_l <- -0.82
  cx_r <-  0.82

  cl <- circle_poly(cx_l, 0, r)
  cr <- circle_poly(cx_r, 0, r)

  ggplot() +
    geom_polygon(data = cl, aes(x, y),
                 fill = COL_A7, color = NA, alpha = 0.18) +
    geom_polygon(data = cr, aes(x, y),
                 fill = COL_A9, color = NA, alpha = 0.18) +
    geom_path(data = cl, aes(x, y),
              color = COL_A7, linewidth = 1.6) +
    geom_path(data = cr, aes(x, y),
              color = COL_A9, linewidth = 1.6) +
    annotate("text",
             x = cx_l - 0.62, y = 0,
             label = format(excl_a7, big.mark = ","),
             size = 10, fontface = "bold", color = COL_A7, hjust = 0.5) +
    annotate("text",
             x = 0, y = 0,
             label = format(shared_n, big.mark = ","),
             size = 9, fontface = "bold", color = COL_SH, hjust = 0.5) +
    annotate("text",
             x = cx_r + 0.62, y = 0,
             label = format(excl_a9, big.mark = ","),
             size = 10, fontface = "bold", color = COL_A9, hjust = 0.5) +
    coord_fixed(xlim = c(-3.0, 3.0), ylim = c(-2.1, 2.55),
                expand = FALSE, clip = "off") +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA),
          plot.margin     = margin(8, 8, 8, 8))
}

# ── Construir paneles ─────────────────────────────────────────

pA7    <- make_venn(n$A7$deg_only, n$A7$hub_only, n$A7$int, "A7")
pA9    <- make_venn(n$A9$deg_only, n$A9$hub_only, n$A9$int, "A9")
pInter <- make_venn_inter(n$inter$excl_a7, n$inter$excl_a9, n$inter$shared)

# Layout: los dos intra arriba, el inter centrado abajo
top_row <- plot_grid(pA7, pA9, ncol = 2)
bot_row <- plot_grid(NULL, pInter, NULL, ncol = 3, rel_widths = c(0.5, 1, 0.5))

caption_txt <- paste0(
  "Top row — Purple = DEGs only  |  Center = DEG x Hub (interest genes)  |  Orange = Hub genes only\n",
  "Bottom row — Blue = HPV-A7 exclusive  |  Green = shared between clades  |  Orange = HPV-A9 exclusive"
)
caption_grob <- ggdraw() +
  draw_label(caption_txt, x = 0.5, y = 0.5,
             hjust = 0.5, vjust = 0.5,
             size = 9, colour = "#777777")

p_final <- plot_grid(
  top_row, bot_row, caption_grob,
  ncol = 1,
  rel_heights = c(1, 1, 0.12)
)

ggsave("6_6_P1_venn_DEG_Hub.pdf", p_final, width = 12, height = 12,
       bg = "white")
ggsave("6_6_P1_venn_DEG_Hub.png", p_final, width = 12, height = 12,
       dpi = 300, bg = "white")
cat("Saved: 6_6_P1_venn_DEG_Hub\n")
