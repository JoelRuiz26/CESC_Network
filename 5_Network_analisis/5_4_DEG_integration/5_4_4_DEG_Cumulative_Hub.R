# ============================================================
# 5_4_4_DEG_Cumulative_Hub.R
# Cumulative DEG count across the top 1%-20% most-connected genes
#
# For each network, as the degree-percentile threshold p grows from
# 1% to 20% (i.e., considering an increasingly larger set of the
# most-connected genes), this plot tracks how many DEGs have been
# accumulated by that point — per network (HPV-A7, HPV-A9) and the
# DEGs shared between both networks within that hub region.
#
# DEGs: adj.P.Val <= 0.01, |logFC| >= 1.5
# ============================================================

setwd("~/CESC_Network/5_Network_analisis/5_4_DEG_integration/")

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

PADJ_THR <- 0.01
LFC_THR  <- 1.5
ZOOM_MAX <- 20

# ── Load DEG + degree tables ──────────────────────────────────
deg_A7 <- read_tsv("5_4_1_DEG_A7_Degree.tsv", show_col_types = FALSE)
deg_A9 <- read_tsv("5_4_1_DEG_A9_Degree.tsv", show_col_types = FALSE)

prep <- function(tbl, network_name) {
  N <- nrow(tbl)
  tbl |>
    mutate(
      is_DEG   = adj.P.Val <= PADJ_THR & abs(logFC) >= LFC_THR,
      pct_rank = 100 * degree_rank / N,
      Network  = network_name
    )
}

a7 <- prep(deg_A7, "HPV-A7")
a9 <- prep(deg_A9, "HPV-A9")

df <- bind_rows(a7, a9) |>
  mutate(Network = factor(Network, levels = c("HPV-A7", "HPV-A9")))

# ============================================================
# Cumulative DEG counts across the 1%-20% hub region
# ============================================================
PCT_SEQ <- seq(1, ZOOM_MAX, by = 1)

cum_df <- map_dfr(PCT_SEQ, function(p) {
  df |>
    filter(pct_rank <= p) |>
    group_by(Network) |>
    summarise(n_deg = sum(is_DEG), .groups = "drop") |>
    mutate(percentile = p)
}) |>
  rename(Series = Network)

shared_df <- map_dfr(PCT_SEQ, function(p) {
  genes_A7 <- a7 |> filter(pct_rank <= p, is_DEG) |> pull(Gene)
  genes_A9 <- a9 |> filter(pct_rank <= p, is_DEG) |> pull(Gene)
  tibble(percentile = p, n_deg = length(intersect(genes_A7, genes_A9)))
}) |>
  mutate(Series = "Shared (A7 ∩ A9)")

plot_df <- bind_rows(cum_df, shared_df) |>
  mutate(Series = factor(Series, levels = c("HPV-A7", "HPV-A9", "Shared (A7 ∩ A9)")))

cat("\nCumulative DEGs across percentile 1-20% (every 5%):\n")
print(plot_df |> filter(percentile %% 5 == 0) |> arrange(percentile, Series))

pal <- c("HPV-A7" = "#003366", "HPV-A9" = "#CC0000", "Shared (A7 ∩ A9)" = "#6C3483")

# ============================================================
# Plot
# ============================================================
p_cum <- ggplot(plot_df, aes(x = percentile, y = n_deg, color = Series)) +

  geom_line(linewidth = 1.2) +
  geom_point(size = 2.2) +

  scale_color_manual(values = pal, name = NULL) +

  scale_x_continuous(breaks = seq(1, ZOOM_MAX, by = 1),
                      minor_breaks = NULL,
                      expand = expansion(mult = c(0.015, 0.05))) +
  scale_y_continuous(breaks = pretty_breaks(8),
                      expand = expansion(mult = c(0.02, 0.12))) +

  labs(
    title = "Differentially Expressed Genes Among the Most Connected Genes",
    x     = "Most connected genes (%)",
    y     = "Number of differentially expressed genes"
  ) +

  theme_bw(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 14,
                                     margin = margin(b = 10)),
    axis.title       = element_text(face = "bold", size = 11),
    axis.text        = element_text(size = 9.5, color = "black"),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    legend.position  = "top",
    legend.text      = element_text(size = 10),
    plot.margin      = margin(12, 16, 10, 12)
  )

ggsave("5_4_4_DEG_Cumulative_Hub.pdf", p_cum,
       width = 10, height = 7, device = cairo_pdf)
ggsave("5_4_4_DEG_Cumulative_Hub.png", p_cum,
       width = 10, height = 7, dpi = 300, bg = "white")

cat("\nSaved: 5_4_4_DEG_Cumulative_Hub.pdf / .png\n")
