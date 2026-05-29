# ============================================================
# 6_6_venn_DEG_Hub.R
# Clean Venn diagram: DEGs x Hub genes — A7 and A9
# ============================================================

setwd("~/CESC_Network/6_DepMap_cacu/")

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
})

# ── Counts ────────────────────────────────────────────────────

t7 <- read_tsv("6_3_gene_table_A7.tsv", show_col_types = FALSE)
t9 <- read_tsv("6_3_gene_table_A9.tsv", show_col_types = FALSE)

n <- list(
  A7 = list(
    total    = nrow(t7),
    deg_only = sum( t7$is_DEG & !t7$is_hub),
    hub_only = sum(!t7$is_DEG &  t7$is_hub),
    int      = sum( t7$is_interest)
  ),
  A9 = list(
    total    = nrow(t9),
    deg_only = sum( t9$is_DEG & !t9$is_hub),
    hub_only = sum(!t9$is_DEG &  t9$is_hub),
    int      = sum( t9$is_interest)
  )
)

# ── Venn builder ─────────────────────────────────────────────

COL_DEG <- "#D35400"   # amber-burnt orange
COL_HUB <- "#6C3483"   # deep purple
COL_INT <- "#1C2833"   # near-black for intersection number

# Build a circle as a closed polygon in data coordinates
circle_poly <- function(cx, cy, r, n = 300) {
  theta <- seq(0, 2 * pi, length.out = n + 1)
  data.frame(x = cx + r * cos(theta), y = cy + r * sin(theta))
}

make_venn <- function(deg_only, hub_only, int_n, cohort) {

  r    <- 1.50
  cx_l <- -0.82
  cx_r <-  0.82

  cl <- circle_poly(cx_l, 0, r)
  cr <- circle_poly(cx_r, 0, r)

  ggplot() +

    # ── Filled circles ────────────────────────────────────────
    geom_polygon(data = cl, aes(x, y),
                 fill = COL_DEG, color = NA, alpha = 0.18) +
    geom_polygon(data = cr, aes(x, y),
                 fill = COL_HUB, color = NA, alpha = 0.18) +

    # ── Outlines ──────────────────────────────────────────────
    geom_path(data = cl, aes(x, y),
              color = COL_DEG, linewidth = 1.6) +
    geom_path(data = cr, aes(x, y),
              color = COL_HUB, linewidth = 1.6) +

    # ── Set name labels ───────────────────────────────────────
    annotate("text",
             x = cx_l - 0.55, y = 1.72,
             label = "DEGs", hjust = 0.5,
             size = 6.5, fontface = "bold", color = COL_DEG) +
    annotate("text",
             x = cx_r + 0.55, y = 1.72,
             label = "Hub genes", hjust = 0.5,
             size = 6.5, fontface = "bold", color = COL_HUB) +

    # ── Counts ────────────────────────────────────────────────
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

    # ── Cohort title ─────────────────────────────────────────
    annotate("text",
             x = 0, y = 2.25,
             label = paste0("HPV-", cohort),
             size = 8, fontface = "bold", color = "#1C2833", hjust = 0.5) +

    coord_fixed(xlim = c(-3.0, 3.0), ylim = c(-2.1, 2.55),
                expand = FALSE, clip = "off") +
    theme_void() +
    theme(plot.background  = element_rect(fill = "white", color = NA),
          plot.margin      = margin(8, 8, 8, 8))
}

pA7 <- make_venn(n$A7$deg_only, n$A7$hub_only, n$A7$int, "A7")
pA9 <- make_venn(n$A9$deg_only, n$A9$hub_only, n$A9$int, "A9")

# ── Combine ───────────────────────────────────────────────────

p_final <- pA7 + pA9 +
  plot_annotation(
    caption = "Orange = DEGs only   |   Center = intersection (interest genes)   |   Purple = Hub genes only",
    theme   = theme(
      plot.caption    = element_text(size = 9.5, color = "#777777", hjust = 0.5,
                                     margin = margin(t = 2, b = 8)),
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

ggsave("6_6_P1_venn_DEG_Hub.pdf", p_final, width = 12, height = 6.5)
ggsave("6_6_P1_venn_DEG_Hub.png", p_final, width = 12, height = 6.5, dpi = 300, bg = "white")
cat("Saved: 6_6_P1_venn_DEG_Hub\n")
