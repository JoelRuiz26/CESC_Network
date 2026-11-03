### ============================================================
### 04_summary_panel.R
###
### Publication-style multi-panel summary of the A7-vs-A9 DEG analysis:
###   A) Volcano plot -- CESC            (n up / n down annotated)
###   B) Volcano plot -- HTMCP-CC         (n up / n down annotated)
###   C) logFC correlation scatter -- shared DEGs, raw values, with a real
###      AXIS BREAK (ggbreak) over the empty |logFC|<1 band
###   D) GSEA (Hallmark) barplot -- union of terms significant in either
###      cohort; bars for a term NOT significant in that specific cohort
###      (shown only for comparison) are rendered at reduced opacity
###
### Colour convention: Up (A7) = red, Down / Up-A9 = blue (inverted from
### the earlier draft, per explicit request).
###
### Reuses the tables already written by 01/03/05, no recomputation here.
###
### -> figures/04_summary_panel.png / .pdf
### ============================================================
suppressPackageStartupMessages({
  library(tidyverse); library(vroom); library(ggrepel); library(patchwork); library(ggbreak); library(ggplotify)
})

out_dir <- "~/CESC_Network/09_functional_similarity_degs"
FF <- "Times New Roman"
LFC_CUT <- 1; P_CUT <- 0.05

## ── palette: each cohort's volcano uses its own project-wide identity
## colours (CESC = navy/red, HTMCP-CC = orange/green); the scatter (C) and
## GSEA (D) panels share one additional purple/teal pair, since they're not
## tied to a single cohort
CESC_UP   <- "#CC0000"; CESC_DN   <- "#003366"
HTMCP_UP  <- "#D55E00"; HTMCP_DN  <- "#009E73"
COL_NS <- "#C9C9C9"
ALT_A7 <- "#7A4EAB"; ALT_A9 <- "#1D9A8C"
COL_BORDER <- "#2B2B2B"; COL_GRID <- "#ECECEC"; COL_STRIP <- "#262626"

full_cesc  <- vroom(file.path(out_dir, "results/01_deg_full_cesc.tsv"),  show_col_types = FALSE)
full_htmcp <- vroom(file.path(out_dir, "results/01_deg_full_htmcp.tsv"), show_col_types = FALSE)
shared     <- vroom(file.path(out_dir, "results/03_shared_degs.tsv"), show_col_types = FALSE)

base_thm <- theme_bw(base_family = FF) +
  theme(axis.title = element_text(size = 18, face = "bold"), axis.text = element_text(size = 16),
        plot.title = element_text(size = 21, face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = COL_GRID, linewidth = 0.35),
        panel.border = element_rect(colour = COL_BORDER, linewidth = 0.8),
        legend.position = "none")

## ── A / B: volcano plots, top 5 up / top 5 down labelled (by adj.P.Val),
## with the n of DEGs on each side annotated directly on the plot
make_volcano <- function(tt, title, col_up, col_dn) {
  d <- tt %>% mutate(grp = case_when(adj.P.Val < P_CUT & logFC > LFC_CUT ~ "Up",
                                     adj.P.Val < P_CUT & logFC < -LFC_CUT ~ "Down",
                                     TRUE ~ "NS"),
                     neglogp = -log10(pmax(adj.P.Val, 1e-300)))
  lab_up <- d %>% filter(grp == "Up")   %>% slice_min(adj.P.Val, n = 5)
  lab_dn <- d %>% filter(grp == "Down") %>% slice_min(adj.P.Val, n = 5)
  n_up <- sum(d$grp == "Up"); n_dn <- sum(d$grp == "Down")
  xrng <- range(d$logFC, na.rm = TRUE)
  ggplot(d, aes(logFC, neglogp)) +
    geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", colour = "grey55", linewidth = 0.4) +
    geom_hline(yintercept = -log10(P_CUT), linetype = "dashed", colour = "grey55", linewidth = 0.4) +
    geom_point(aes(colour = grp), size = 1.8, alpha = 0.5) +
    geom_text_repel(data = lab_dn, aes(label = gene), family = FF, fontface = "bold", size = 5.6,
                    max.overlaps = Inf, seed = 1, box.padding = 0.4, force = 2, force_pull = 0.6,
                    segment.size = 0.25, segment.colour = "grey50", min.segment.length = 0.1,
                    xlim = c(xrng[1], -LFC_CUT * 0.7)) +
    geom_text_repel(data = lab_up, aes(label = gene), family = FF, fontface = "bold", size = 5.6,
                    max.overlaps = Inf, seed = 1, box.padding = 0.4, force = 2, force_pull = 0.6,
                    segment.size = 0.25, segment.colour = "grey50", min.segment.length = 0.1,
                    xlim = c(LFC_CUT * 0.7, xrng[2])) +
    annotate("label", x = -Inf, y = -Inf, hjust = -0.08, vjust = -0.5, family = FF, fontface = "bold", size = 5.6,
             colour = col_dn, fill = "white", label.r = unit(0.2, "lines"),
             label.padding = unit(0.35, "lines"), label = sprintf("Down (A9): n = %d", n_dn)) +
    annotate("label", x = Inf, y = -Inf, hjust = 1.08, vjust = -0.5, family = FF, fontface = "bold", size = 5.6,
             colour = col_up, fill = "white", label.r = unit(0.2, "lines"),
             label.padding = unit(0.35, "lines"), label = sprintf("Up (A7): n = %d", n_up)) +
    scale_colour_manual(values = c("Up" = col_up, "Down" = col_dn, "NS" = COL_NS)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    labs(title = title, x = "LFC", y = expression(-log[10]~italic(p)[adj])) +
    base_thm
}

pA <- make_volcano(full_cesc,  "CESC — HPV-A7 vs HPV-A9",     CESC_UP,  CESC_DN)
pB <- make_volcano(full_htmcp, "HTMCP-CC — HPV-A7 vs HPV-A9", HTMCP_UP, HTMCP_DN)

## ── C: logFC correlation scatter, shared DEGs, raw values -- with a real
## axis break (ggbreak) collapsing the empty |logFC|<1 band on BOTH axes
ct <- cor.test(shared$logFC_cesc, shared$logFC_htmcp, method = "spearman")
shared2 <- shared %>% mutate(direction = if_else(logFC_cesc > 0, "Up (A7 vs A9)", "Down (A7 vs A9)"),
                             mean_lfc = (logFC_cesc + logFC_htmcp) / 2)
n_up <- sum(shared2$direction == "Up (A7 vs A9)"); n_dn <- sum(shared2$direction == "Down (A7 vs A9)")
lab_c <- bind_rows(shared2 %>% slice_max(mean_lfc, n = 8), shared2 %>% slice_min(mean_lfc, n = 8))
rng0 <- range(c(shared2$logFC_cesc, shared2$logFC_htmcp))
pad <- diff(rng0) * 0.09
rng <- c(rng0[1] - pad, rng0[2] + pad)
rho_lab <- sprintf("Spearman ρ = %.3f\np < 2.2×10⁻¹⁶\nn = %d shared genes", unname(ct$estimate), nrow(shared))

pC_base <- ggplot(shared2, aes(logFC_cesc, logFC_htmcp)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.55) +
  geom_point(aes(colour = direction), size = 3.0, alpha = 0.55) +
  geom_text_repel(data = lab_c, aes(label = gene), family = FF, fontface = "bold", size = 5.6,
                  max.overlaps = Inf, seed = 1, box.padding = 0.35, force = 1.5, force_pull = 0.8,
                  max.iter = 20000, max.time = 2, min.segment.length = 0.1,
                  segment.size = 0.3, segment.colour = "grey45") +
  annotate("label", x = rng[1], y = rng[2], hjust = 0, vjust = 1, family = FF, fontface = "bold",
           size = 6.4, colour = "grey10", fill = "#F7F7F5", label.r = unit(0.3, "lines"),
           label.padding = unit(0.6, "lines"), label = rho_lab) +
  ## NOTE: finite coordinates on purpose -- with BOTH axes broken (ggbreak),
  ## any label placed at -Inf/Inf gets duplicated across every sub-panel;
  ## a finite point that falls inside only one quadrant's own sub-range
  ## renders exactly once. Count boxes pinned to the padded plot limits
  ## (`rng`, just outside the data's own extent `rng0`) at the bottom-left
  ## / top-right corners, clearing the most extreme gene points/labels.
  annotate("label", x = rng[1], y = rng[1], hjust = 0, vjust = 0, family = FF, fontface = "bold",
           size = 5.6, colour = ALT_A9, fill = "white", label.r = unit(0.2, "lines"), label.padding = unit(0.35, "lines"),
           label = sprintf("Down (A9): n = %d", n_dn)) +
  annotate("label", x = rng[2], y = rng[2], hjust = 1, vjust = 1, family = FF, fontface = "bold",
           size = 5.6, colour = ALT_A7, fill = "white", label.r = unit(0.2, "lines"), label.padding = unit(0.35, "lines"),
           label = sprintf("Up (A7): n = %d", n_up)) +
  scale_colour_manual(values = c("Up (A7 vs A9)" = ALT_A7, "Down (A7 vs A9)" = ALT_A9), name = NULL) +
  scale_x_continuous(limits = rng, breaks = scales::breaks_width(1)) +
  scale_y_continuous(limits = rng, breaks = scales::breaks_width(1)) +
  labs(title = "LFC Correlation — CESC vs HTMCP-CC",
       x = "LFC (A7 vs A9) — CESC", y = "LFC (A7 vs A9) — HTMCP-CC") +
  base_thm + theme(legend.position = "none", plot.margin = margin(16, 16, 10, 10))

pC_break <- pC_base +
  scale_x_break(c(-LFC_CUT, LFC_CUT), scales = "fixed", space = 0.04, ticklabels = NULL) +
  scale_y_break(c(-LFC_CUT, LFC_CUT), scales = "fixed", space = 0.04, ticklabels = NULL)
pC <- wrap_elements(full = as.ggplot(pC_break))

## ── D: GSEA (Hallmark) barplot, reusing results already written by
## 05_gsea_hallmark.R. Union of terms significant (adj.P<0.05) in EITHER
## cohort. Direction (purple = A7, teal = A9) is self-evident from which
## side of zero the bar extends and from the panel title, so it gets NO
## legend at all. The only thing the legend explains is significance:
## solid colour = significant, pale tint of the SAME hue = not significant
## (shown only for comparison) -- one simple 2-entry legend, generic grey
## swatches standing in for "solid vs pale" regardless of hue.
tint <- function(hex, amt = 0.6) {
  rgb_mix <- (col2rgb(hex) * (1 - amt) + 255 * amt) / 255
  rgb(rgb_mix[1], rgb_mix[2], rgb_mix[3])
}
gsea_cesc  <- vroom(file.path(out_dir, "results/05_gsea_hallmark_cesc.tsv"),  show_col_types = FALSE)
gsea_htmcp <- vroom(file.path(out_dir, "results/05_gsea_hallmark_htmcp.tsv"), show_col_types = FALSE)
sig_terms <- union(gsea_cesc$HM[gsea_cesc$p.adjust < 0.05], gsea_htmcp$HM[gsea_htmcp$p.adjust < 0.05])
gsea_both <- bind_rows(gsea_cesc, gsea_htmcp) %>% filter(HM %in% sig_terms) %>%
  mutate(cohort = factor(cohort, c("CESC", "HTMCP-CC")),
         sig_here = p.adjust < 0.05,
         fill_cat = factor(
           case_when(NES > 0 & sig_here  ~ "Enriched in A7 — significant",
                    NES > 0 & !sig_here ~ "Enriched in A7 — not significant",
                    NES < 0 & sig_here  ~ "Enriched in A9 — significant",
                    TRUE                ~ "Enriched in A9 — not significant"),
           levels = c("Enriched in A7 — significant", "Enriched in A7 — not significant",
                     "Enriched in A9 — significant", "Enriched in A9 — not significant")),
         significance = factor(if_else(sig_here, "Significant", "Not significant"),
                               levels = c("Significant", "Not significant")))
term_order <- gsea_cesc %>% filter(HM %in% sig_terms) %>% arrange(NES) %>% pull(HM)
term_order <- c(setdiff(sig_terms, term_order), term_order)
gsea_both <- gsea_both %>% mutate(HM = factor(HM, term_order))

fill_cats <- c("Enriched in A7 — significant" = ALT_A7, "Enriched in A7 — not significant" = tint(ALT_A7),
              "Enriched in A9 — significant" = ALT_A9, "Enriched in A9 — not significant" = tint(ALT_A9))
pD <- ggplot(gsea_both, aes(NES, HM, fill = fill_cat)) +
  geom_col(width = 0.72, colour = "grey25", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.4) +
  ## invisible layer purely to drive a clean, generic 2-entry legend
  geom_col(aes(alpha = significance), fill = NA, colour = NA, width = 0.72, show.legend = TRUE) +
  facet_wrap(~cohort, nrow = 1) +
  scale_fill_manual(values = fill_cats, guide = "none") +
  scale_alpha_manual(values = c("Significant" = 1, "Not significant" = 1), name = NULL,
                     guide = guide_legend(override.aes = list(fill = c("grey25", "grey85"), colour = "grey40"))) +
  labs(title = "GSEA (Hallmark) — HPV-A7 vs HPV-A9",
       x = "Normalized Enrichment Score (NES)", y = NULL) +
  theme_bw(base_family = FF) +
  theme(strip.background = element_rect(fill = COL_STRIP, colour = NA),
        strip.text = element_text(colour = "white", size = 19, face = "bold"),
        axis.text.y = element_text(size = 17, face = "bold", colour = "grey10"),
        axis.text.x = element_text(size = 16), axis.title.x = element_text(size = 18, face = "bold"),
        plot.title = element_text(size = 21, face = "bold", hjust = 0.5),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(colour = COL_GRID, linewidth = 0.35),
        panel.border = element_rect(colour = COL_BORDER, linewidth = 0.7),
        panel.spacing = unit(1.2, "lines"),
        legend.position = "bottom", legend.text = element_text(size = 14))

## ── assemble: row 1 = two volcanoes; row 2 = scatter (narrower) beside the
## GSEA barplot (wider) -- a more compact, balanced layout than stacking
## everything in one column. Tags are applied INSIDE each row (A/B, then
## C/D) before wrap_elements() makes that row opaque, so all four panels
## get their own letter instead of collapsing to one tag per row.
rowAB <- wrap_elements(full = (
  (pA | pB) + plot_annotation(tag_levels = list(c("A", "B"))) &
    theme(plot.tag = element_text(family = FF, size = 24, face = "bold"))
))
rowCD <- wrap_elements(full = (
  (pC + plot_layout(widths = 1) | pD + plot_layout(widths = 1.7)) +
    plot_annotation(tag_levels = list(c("C", "D"))) &
    theme(plot.tag = element_text(family = FF, size = 24, face = "bold"))
))

fig <- (rowAB / rowCD) +
  plot_layout(heights = c(1, 1.7)) &
  theme(plot.tag = element_text(family = FF, size = 26, face = "bold"))

ggsave(file.path(out_dir, "figures/04_summary_panel.png"), fig, width = 24, height = 18, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(file.path(out_dir, "figures/04_summary_panel.pdf"), fig, width = 24, height = 18, device = cairo_pdf, bg = "white", limitsize = FALSE)
message("Saved: figures/04_summary_panel.png / .pdf")
