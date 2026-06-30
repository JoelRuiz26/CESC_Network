# ============================================================
# 5_4_2_Scatter_Degree_SharedDEGs.R
# Degree scatter (A7 vs A9) — restricted to shared DEGs only
# Spearman correlation on log-log scale
# ============================================================

setwd("~/CESC_Network/5_Network_analisis/5_4_DEG_integration/")

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

PADJ_THR <- 0.01
LFC_THR  <- 1.5

# ── Load DEG + degree tables ──────────────────────────────────
deg_A7 <- read_tsv("5_4_1_DEG_A7_Degree.tsv", show_col_types = FALSE)
deg_A9 <- read_tsv("5_4_1_DEG_A9_Degree.tsv", show_col_types = FALSE)

# ── Filter to significant DEGs in each cohort ─────────────────
sig_A7 <- deg_A7 |> filter(adj.P.Val <= PADJ_THR, abs(logFC) >= LFC_THR)
sig_A9 <- deg_A9 |> filter(adj.P.Val <= PADJ_THR, abs(logFC) >= LFC_THR)

# ── Keep only shared DEGs ────────────────────────────────────
shared_genes <- intersect(sig_A7$Gene, sig_A9$Gene)

df <- sig_A7 |>
  filter(Gene %in% shared_genes) |>
  select(Gene, logFC_A7 = logFC, degree_A7 = degree) |>
  inner_join(
    sig_A9 |>
      filter(Gene %in% shared_genes) |>
      select(Gene, logFC_A9 = logFC, degree_A9 = degree),
    by = "Gene"
  ) |>
  filter(degree_A7 > 0, degree_A9 > 0) |>
  mutate(
    direction = case_when(
      logFC_A7 > 0 & logFC_A9 > 0 ~ "Both Up",
      logFC_A7 < 0 & logFC_A9 < 0 ~ "Both Down",
      TRUE                          ~ "Discordant"
    ),
    direction = factor(direction, levels = c("Both Up", "Both Down", "Discordant"))
  )

cat(sprintf("Shared DEGs in scatter: %d\n", nrow(df)))
cat(sprintf("  Both Up: %d  |  Both Down: %d  |  Discordant: %d\n",
            sum(df$direction == "Both Up"),
            sum(df$direction == "Both Down"),
            sum(df$direction == "Discordant")))

# ── Spearman correlation ──────────────────────────────────────
sp       <- cor.test(df$degree_A7, df$degree_A9, method = "spearman")
rho_val  <- round(sp$estimate, 3)
rho_str  <- case_when(
  abs(rho_val) >= 0.70 ~ "strong",
  abs(rho_val) >= 0.50 ~ "moderate",
  abs(rho_val) >= 0.30 ~ "low-moderate",
  TRUE                  ~ "weak"
)

ann_label <- sprintf("Spearman rho = %.3f  (%s)\np < 2.2e-16\nn = %d shared DEGs",
                     rho_val, rho_str, nrow(df))

# ── Palette ───────────────────────────────────────────────────
pal <- c(
  "Both Up"    = "#C0392B",
  "Both Down"  = "#1A5276",
  "Discordant" = "#888888"
)

# ── Plot ──────────────────────────────────────────────────────
p <- ggplot(df, aes(x = degree_A7, y = degree_A9, color = direction)) +

  # diagonal: perfect concordance
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", linewidth = 0.55, color = "grey50") +

  # points
  geom_point(alpha = 0.35, size = 1.2, shape = 16) +

  # stats box
  annotate("label",
           x = min(df$degree_A7)^1.05,
           y = max(df$degree_A9) * 0.88,
           hjust = 0, vjust = 1,
           label = ann_label,
           size  = 4.2, fontface = "bold",
           fill  = alpha("white", 0.90), color = "grey15",
           label.padding = unit(0.4, "lines")) +

  # log-log scales
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x)),
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x)),
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  annotation_logticks(sides = "bl", color = "grey50", linewidth = 0.3) +

  scale_color_manual(
    values = pal,
    name   = NULL,
    guide  = guide_legend(override.aes = list(alpha = 0.9, size = 3.5))
  ) +

  labs(
    title    = "Network Degree — Shared DEGs (HPV-A7 vs HPV-A9)",
    subtitle = sprintf(
      "Log-log scale  |  Diagonal = perfect concordance  |  adj.P.Val <= %.2f, |logFC| >= %g",
      PADJ_THR, LFC_THR
    ),
    x     = "Degree in HPV-A7 network  (log scale)",
    y     = "Degree in HPV-A9 network  (log scale)",
    color = NULL
  ) +

  theme_bw(base_size = 13) +
  theme(
    plot.title        = element_text(face = "bold", hjust = 0.5, size = 14,
                                     margin = margin(b = 4)),
    plot.subtitle     = element_text(hjust = 0.5, size = 9.5, color = "grey35",
                                     margin = margin(b = 8)),
    axis.title        = element_text(face = "bold", size = 12),
    axis.text         = element_text(size = 10, color = "black"),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor  = element_blank(),
    legend.position   = c(0.88, 0.12),
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = alpha("white", 0.88),
                                     color = "grey60", linewidth = 0.4),
    legend.text       = element_text(size = 11),
    legend.key        = element_blank(),
    plot.background   = element_rect(fill = "white", color = NA),
    plot.margin       = margin(12, 15, 10, 12)
  )

ggsave("5_4_2_Scatter_Degree_SharedDEGs.pdf", p,
       width = 8, height = 7, device = cairo_pdf)
ggsave("5_4_2_Scatter_Degree_SharedDEGs.png", p,
       width = 8, height = 7, dpi = 300, bg = "white")

cat("Saved: 5_4_2_Scatter_Degree_SharedDEGs.pdf / .png\n")
