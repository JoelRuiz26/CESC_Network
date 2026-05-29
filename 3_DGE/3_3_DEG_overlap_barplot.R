# ============================================================
# 3_3_DEG_overlap_barplot.R
# DEG set overlap — horizontal diverging bar chart
# Replaces Venn: shows set membership + regulation direction
# ============================================================

setwd("~/CESC_Network/3_DGE/")

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

PADJ_THR <- 0.01
LFC_THR  <- 1

# ── Load & classify ───────────────────────────────────────────
de_A7 <- read_tsv("3_1_limma_A7_vs_SNT_allGenes.tsv", show_col_types = FALSE)
de_A9 <- read_tsv("3_1_limma_A9_vs_SNT_allGenes.tsv", show_col_types = FALSE)

sig_A7 <- de_A7 |> filter(adj.P.Val <= PADJ_THR, abs(logFC) >= LFC_THR)
sig_A9 <- de_A9 |> filter(adj.P.Val <= PADJ_THR, abs(logFC) >= LFC_THR)

shared  <- intersect(sig_A7$Gene, sig_A9$Gene)
only_A7 <- setdiff(sig_A7$Gene, shared)
only_A9 <- setdiff(sig_A9$Gene, shared)

sh_lfc <- sig_A7 |>
  filter(Gene %in% shared) |>
  select(Gene, lfc_A7 = logFC) |>
  left_join(sig_A9 |> select(Gene, lfc_A9 = logFC), by = "Gene")

# ── Count table ───────────────────────────────────────────────
counts <- tibble(
  set = factor(
    rep(c("HPV-A7 only", "Shared  (A7 ∩ A9)", "HPV-A9 only"), each = 2),
    levels = c("HPV-A9 only", "Shared  (A7 ∩ A9)", "HPV-A7 only")
  ),
  direction = rep(c("Up", "Down"), 3),
  n = c(
    sum(filter(sig_A7, Gene %in% only_A7)$logFC > 0),   # A7-only Up
    sum(filter(sig_A7, Gene %in% only_A7)$logFC < 0),   # A7-only Down
    sum(sh_lfc$lfc_A7 > 0),                              # Shared Up
    sum(sh_lfc$lfc_A7 < 0),                              # Shared Down
    sum(filter(sig_A9, Gene %in% only_A9)$logFC > 0),   # A9-only Up
    sum(filter(sig_A9, Gene %in% only_A9)$logFC < 0)    # A9-only Down
  )
) |>
  mutate(
    n_plot    = ifelse(direction == "Down", -n, n),
    fill_id   = paste(set, direction),
    label_x   = ifelse(direction == "Up", n + 60, -n - 60),
    hjust_lab = ifelse(direction == "Up", 0, 1)
  )

# ── Colors ────────────────────────────────────────────────────
pal <- c(
  "HPV-A7 only Up"                  = "#003366",
  "HPV-A7 only Down"                = "#5B9BD5",
  "Shared  (A7 ∩ A9) Up"      = "#5B2C6F",
  "Shared  (A7 ∩ A9) Down"    = "#A569BD",
  "HPV-A9 only Up"                  = "#7B241C",
  "HPV-A9 only Down"                = "#E57373"
)

# ── Shared % annotation ───────────────────────────────────────
pct_A7 <- round(100 * length(shared) / nrow(sig_A7))
pct_A9 <- round(100 * length(shared) / nrow(sig_A9))

# ── Plot ──────────────────────────────────────────────────────
p <- ggplot(counts, aes(y = set, x = n_plot, fill = fill_id)) +
  geom_col(width = 0.52, color = "white", linewidth = 0.5) +

  # center line
  geom_vline(xintercept = 0, color = "#333333", linewidth = 0.6) +

  # count labels outside bars
  geom_text(
    aes(x = label_x, label = comma(n), hjust = hjust_lab),
    size = 5.5, fontface = "bold", color = "grey20"
  ) +

  # axis direction labels
  annotate("text", x = -3600, y = 3.55,
           label = "<- Down-regulated", hjust = 0,
           size = 5, color = "grey50", fontface = "italic") +
  annotate("text", x =  3600, y = 3.55,
           label = "Up-regulated ->", hjust = 1,
           size = 5, color = "grey50", fontface = "italic") +

  scale_fill_manual(
    values = pal,
    breaks = c("HPV-A7 only Up", "HPV-A7 only Down",
               "Shared  (A7 ∩ A9) Up", "Shared  (A7 ∩ A9) Down",
               "HPV-A9 only Up", "HPV-A9 only Down"),
    labels = c("HPV-A7 only  Up", "HPV-A7 only  Down",
               "Shared  Up", "Shared  Down",
               "HPV-A9 only  Up", "HPV-A9 only  Down"),
    name = NULL
  ) +
  scale_x_continuous(
    labels = function(x) comma(abs(x)),
    limits = c(-4200, 4200),
    expand = c(0, 0)
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = FALSE)) +

  labs(
    title    = "Differentially Expressed Genes — HPV-A7 and HPV-A9 vs Normal",
    subtitle = sprintf(
      "adj.P.Val <= %.2f  |  |logFC| >= %g  |  n(A7 DEGs) = %s  |  n(A9 DEGs) = %s  |  n(shared) = %s",
      PADJ_THR, LFC_THR,
      comma(nrow(sig_A7)), comma(nrow(sig_A9)), comma(length(shared))
    ),
    x       = "Number of DEGs",
    y       = NULL,
    caption = sprintf(
      "%d%% of HPV-A7 DEGs and %d%% of HPV-A9 DEGs are shared between both groups",
      pct_A7, pct_A9
    )
  ) +

  theme_classic(base_size = 17) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 19),
    plot.subtitle   = element_text(hjust = 0.5, size = 12, color = "grey35",
                                   margin = margin(b = 10)),
    axis.text.y     = element_text(size = 16, face = "bold", color = "#222222"),
    axis.text.x     = element_text(size = 14),
    axis.title.x    = element_text(size = 15, face = "bold"),
    axis.line.y     = element_blank(),
    axis.ticks.y    = element_blank(),
    legend.position = "top",
    legend.text     = element_text(size = 13),
    legend.key.size = unit(0.55, "cm"),
    legend.spacing.x = unit(0.4, "cm"),
    panel.grid.major.x = element_line(color = "#EBEBEB", linewidth = 0.35),
    plot.background = element_rect(fill = "white", color = NA),
    plot.caption    = element_text(size = 13, color = "grey35", hjust = 0.5,
                                   face = "italic", margin = margin(t = 6)),
    plot.margin     = margin(12, 20, 12, 12)
  )

ggsave("3_3_DEG_overlap_barplot.pdf", p,
       width = 9, height = 6.5, device = cairo_pdf)
ggsave("3_3_DEG_overlap_barplot.png", p,
       width = 9, height = 6.5, dpi = 300, bg = "white")

cat("Saved: 3_3_DEG_overlap_barplot.pdf / .png\n")
cat(sprintf("  HPV-A7 only : %d up  /  %d down\n",
            counts$n[counts$set=="HPV-A7 only" & counts$direction=="Up"],
            counts$n[counts$set=="HPV-A7 only" & counts$direction=="Down"]))
cat(sprintf("  Shared      : %d up  /  %d down\n",
            counts$n[grepl("Shared", counts$set) & counts$direction=="Up"],
            counts$n[grepl("Shared", counts$set) & counts$direction=="Down"]))
cat(sprintf("  HPV-A9 only : %d up  /  %d down\n",
            counts$n[counts$set=="HPV-A9 only" & counts$direction=="Up"],
            counts$n[counts$set=="HPV-A9 only" & counts$direction=="Down"]))
