# ============================================================
# 3_2_volcano_A7_A9_vs_CTL.R
# Volcano plots: HPV-A7 vs Normal and HPV-A9 vs Normal
# ============================================================

setwd("~/CESC_Network/3_DGE/")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(patchwork)
})

PADJ_THR <- 0.01
LFC_THR  <- 1
Y_CAP    <- 50      # genes with -log10(adj.P) > 50 drawn as triangles at cap line

COL_UP <- "#C0392B"   # red  — up-regulated in tumor
COL_DN <- "#1A5276"   # blue — down-regulated in tumor
COL_NS <- "#BBBBBB"   # grey — not significant

# ── Load data ─────────────────────────────────────────────────
de_A7 <- read_tsv("3_1_limma_A7_vs_SNT_allGenes.tsv", show_col_types = FALSE)
de_A9 <- read_tsv("3_1_limma_A9_vs_SNT_allGenes.tsv", show_col_types = FALSE)

# ── Volcano builder ───────────────────────────────────────────
make_volcano <- function(de, cohort, col_up = COL_UP, col_dn = COL_DN) {

  df <- de %>%
    mutate(
      status = case_when(
        adj.P.Val <= PADJ_THR & logFC >=  LFC_THR ~ "Up",
        adj.P.Val <= PADJ_THR & logFC <= -LFC_THR ~ "Down",
        TRUE ~ "NS"
      ),
      status      = factor(status, levels = c("Up", "Down", "NS")),
      neg_log10_p = -log10(adj.P.Val),
      y_plot      = pmin(neg_log10_p, Y_CAP),
      capped      = neg_log10_p > Y_CAP
    )

  n_up <- sum(df$status == "Up")
  n_dn <- sum(df$status == "Down")

  pal <- c("Up" = col_up, "Down" = col_dn, "NS" = COL_NS)

  ggplot(df, aes(x = logFC, y = y_plot, color = status)) +

    # ── NS layer (bottom) ────────────────────────────────────
    geom_point(data = filter(df, status == "NS"),
               size = 0.65, alpha = 0.25, shape = 16) +

    # ── Significant — regular points ─────────────────────────
    geom_point(data = filter(df, status != "NS", !capped),
               size = 1.2, alpha = 0.70, shape = 16) +

    # ── Significant — capped (triangles at Y_CAP) ────────────
    geom_point(data = filter(df, status != "NS", capped),
               size = 2.0, alpha = 0.90, shape = 17) +

    # ── Reference lines ──────────────────────────────────────
    geom_hline(yintercept = -log10(PADJ_THR),
               linetype = "dashed", color = "grey40", linewidth = 0.45) +
    geom_vline(xintercept = c(-LFC_THR, LFC_THR),
               linetype = "dashed", color = "grey40", linewidth = 0.45) +
    geom_hline(yintercept = Y_CAP,
               linetype = "dotted", color = "grey60", linewidth = 0.35) +

    # ── Scales ────────────────────────────────────────────────
    scale_color_manual(
      values = pal,
      labels = c(
        "Up"   = paste0("Up in tumor (n = ", n_up, ")"),
        "Down" = paste0("Down in tumor (n = ", n_dn, ")"),
        "NS"   = "Not significant"
      ),
      name = NULL
    ) +
    scale_y_continuous(
      breaks = c(0, 10, 20, 30, 40, 50),
      labels = c("0", "10", "20", "30", "40", paste0(">", Y_CAP)),
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    scale_x_continuous(expand = expansion(mult = 0.04)) +

    labs(
      title    = paste0("HPV-", cohort, " vs Normal"),
      subtitle = paste0("adj.P.Val <= ", PADJ_THR,
                        "  |  |logFC| >= ", LFC_THR,
                        "  |  triangles: -log10(adj.p) > ", Y_CAP),
      x = expression(log[2]~Fold~Change~~(Tumor~vs~Normal)),
      y = expression(-log[10]~(adj.~p~value))
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title        = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle     = element_text(hjust = 0.5, size = 9, color = "grey35"),
      legend.position   = c(0.82, 0.12),
      legend.background = element_rect(fill = alpha("white", 0.80), color = NA),
      legend.text       = element_text(size = 9),
      panel.border      = element_rect(color = "#CCCCCC", fill = NA, linewidth = 0.5),
      plot.background   = element_rect(fill = "white", color = NA)
    )
}

# ── Build plots ───────────────────────────────────────────────
pA7 <- make_volcano(de_A7, "A7")
pA9 <- make_volcano(de_A9, "A9")

# ── Combined figure ───────────────────────────────────────────
p_combined <- pA7 + pA9 +
  plot_annotation(
    title = "Differential Gene Expression — HPV A7 and A9 vs Normal Tissue",
    theme = theme(
      plot.title      = element_text(face = "bold", hjust = 0.5, size = 15),
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

ggsave("3_2_volcano_A7_A9_vs_CTL.pdf", p_combined,
       width = 16, height = 7, device = cairo_pdf)
ggsave("3_2_volcano_A7_A9_vs_CTL.png", p_combined,
       width = 16, height = 7, dpi = 300, bg = "white")

# ── Individual saves ──────────────────────────────────────────
ggsave("3_2_volcano_A7_vs_CTL.pdf", pA7,
       width = 8, height = 7, device = cairo_pdf)
ggsave("3_2_volcano_A7_vs_CTL.png", pA7,
       width = 8, height = 7, dpi = 300, bg = "white")

ggsave("3_2_volcano_A9_vs_CTL.pdf", pA9,
       width = 8, height = 7, device = cairo_pdf)
ggsave("3_2_volcano_A9_vs_CTL.png", pA9,
       width = 8, height = 7, dpi = 300, bg = "white")

cat("Saved:\n")
cat("  3_2_volcano_A7_A9_vs_CTL.pdf / .png  (combined)\n")
cat("  3_2_volcano_A7_vs_CTL.pdf / .png\n")
cat("  3_2_volcano_A9_vs_CTL.pdf / .png\n")

# ── Console summary ───────────────────────────────────────────
for (nm in c("A7", "A9")) {
  de  <- if (nm == "A7") de_A7 else de_A9
  up  <- sum(de$adj.P.Val <= PADJ_THR & de$logFC >=  LFC_THR, na.rm = TRUE)
  dn  <- sum(de$adj.P.Val <= PADJ_THR & de$logFC <= -LFC_THR, na.rm = TRUE)
  cat(sprintf("  HPV-%s: %d up, %d down  (total DEG = %d)\n", nm, up, dn, up + dn))
}
