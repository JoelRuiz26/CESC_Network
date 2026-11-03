### ============================================================
### figure_S3_combined.R
###
### CESC + Uganda combined view of Figure S3 (elbow thresholds). Does NOT
### touch the original CESC-only figure. Grouped by ANALYSIS TYPE: each row
### is one measure x clade (A7 connectivity, A9 connectivity, A7 dependency,
### A9 dependency), CESC and Uganda side by side under the same tag letter.
###
###   A) HPV-A7, connectivity  [CESC | Uganda]
###   B) HPV-A9, connectivity  [CESC | Uganda]
###   C) HPV-A7, dependency    [CESC | Uganda]
###   D) HPV-A9, dependency    [CESC | Uganda]
###
### -> figures/validation_uganda/Figure_S3_elbow_thresholds_uganda.png/.pdf
### ============================================================
suppressPackageStartupMessages({
  library(tidyverse); library(vroom); library(patchwork); library(ggrepel)
})

fig_dir <- "~/CESC_Network/figures/validation_uganda"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
FF <- "Times New Roman"

elbow <- function(v) {
  n <- length(v); x <- (seq_len(n) - 1) / (n - 1); y <- (v - v[1]) / (v[n] - v[1])
  i <- which.max(abs(y - x)); list(idx = i, val = v[i])
}
base <- theme_bw(base_family = FF) +
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"), legend.position = "none",
        panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey95", linewidth = .3),
        panel.border = element_rect(colour = "grey35", linewidth = .7), plot.margin = margin(8, 10, 8, 8))

rank_panel <- function(v, e, col, ttl, lab, corner) {
  n <- length(v); chord_y <- v[1] + (v[n] - v[1]) * (e$idx - 1) / (n - 1)
  ggplot(tibble(rank = seq_len(n), y = v), aes(rank, y)) +
    geom_segment(x = 1, y = v[1], xend = n, yend = v[n], colour = "grey62", linetype = "22", linewidth = .5) +
    geom_line(colour = col, linewidth = 1.0) +
    geom_segment(x = e$idx, y = v[e$idx], xend = e$idx, yend = chord_y, colour = "grey30", linewidth = .55) +
    geom_point(x = e$idx, y = e$val, colour = "grey10", size = 3) +
    annotate("text", x = Inf, y = if (corner == "top") Inf else -Inf, family = FF, fontface = "bold",
             size = 5.2, label = lab, hjust = 1.06, vjust = if (corner == "top") 1.6 else -1.0) +
    labs(title = ttl, x = "Gene rank", y = if (corner == "top") "Degree" else "Mean CRISPR gene effect") + base
}

analyse_cohort <- function(res_dir) {
  tab <- vroom(file.path(res_dir, "6_2_gene_dep_degree_table.tsv"), show_col_types = FALSE) %>%
    filter(!is.na(degree), !is.na(mean_dep)) %>% distinct(network, gene, .keep_all = TRUE)
  nets <- sort(unique(tab$network))
  map(nets, function(nk) {
    d <- tab %>% filter(network == nk)
    kv <- sort(d$degree, decreasing = TRUE); dv <- sort(d$mean_dep)
    list(net = nk, kv = kv, dv = dv, ek = elbow(kv), ed = elbow(dv))
  }) %>% set_names(nets)
}

cesc   <- analyse_cohort("~/CESC_Network/07_dependency/results")
uganda <- analyse_cohort("~/CESC_Network/07_dependency/validation_uganda/results")

mk <- function(i, col, cohort_label, measure) {
  if (measure == "k") rank_panel(i$kv, i$ek, col, paste0(i$net, " connectivity — ", cohort_label),
                                paste0("degree ≥ ", i$ek$val), "top")
  else rank_panel(i$dv, i$ed, col, paste0(i$net, " dependency — ", cohort_label),
                  paste0("gene effect ≤ ", sprintf("%.3f", i$ed$val)), "bottom")
}

rowA <- wrap_elements(full = (mk(cesc[["HPV-A7"]], "#4E79A7", "CESC", "k") | mk(uganda[["HPV-A7"]], "#D55E00", "Uganda", "k")))
rowB <- wrap_elements(full = (mk(cesc[["HPV-A9"]], "#CB6A5C", "CESC", "k") | mk(uganda[["HPV-A9"]], "#009E73", "Uganda", "k")))
rowC <- wrap_elements(full = (mk(cesc[["HPV-A7"]], "#4E79A7", "CESC", "d") | mk(uganda[["HPV-A7"]], "#D55E00", "Uganda", "d")))
rowD <- wrap_elements(full = (mk(cesc[["HPV-A9"]], "#CB6A5C", "CESC", "d") | mk(uganda[["HPV-A9"]], "#009E73", "Uganda", "d")))

fig_supp <- (rowA / rowB / rowC / rowD) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = FF, size = 22, face = "bold"))

ggsave(file.path(fig_dir, "Figure_S3_elbow_thresholds_uganda.png"), fig_supp, width = 14, height = 18, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(file.path(fig_dir, "Figure_S3_elbow_thresholds_uganda.pdf"), fig_supp, width = 14, height = 18, device = cairo_pdf, bg = "white", limitsize = FALSE)
message("Saved: Figure_S3_elbow_thresholds_uganda.png / .pdf  (A=A7 connectivity, B=A9 connectivity, C=A7 dependency, D=A9 dependency; each CESC | Uganda)")
