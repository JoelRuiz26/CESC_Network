### ============================================================
### 6_4_Hub_threshold_sensitivity.R
###
### Robustness / stability analysis for the ONLY tunable parameter
### in the hub-dependency result: the hub definition (top X% by
### degree). The dependency side is threshold-free (Wilcoxon +
### Cliff's delta on the full score distribution), so it needs no
### sensitivity sweep.
###
### For each hub quantile in [0.70, 0.95] we recompute:
###   * Cliff's delta (effect size, hub vs non-hub)
###   * Wilcoxon rank-sum p-value
### and show the effect size is stable/monotone and always highly
### significant — i.e. the conclusion does not hinge on the 10% cut.
###
### Outputs:
###   6_4_hub_threshold_sensitivity.tsv
###   6_4_A_hub_threshold_sensitivity.png
### ============================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(vroom)
  library(patchwork)
})

dep_dir <- "~/CESC_Network/6_DepMap_cacu"
net_colors <- c("HPV-A7" = "#003366", "HPV-A9" = "#CC2200")

gene_tbl <- vroom(file.path(dep_dir, "6_2_gene_dep_degree_table.tsv"),
                  show_col_types = FALSE)

## ── Cliff's delta via Mann-Whitney U (two-sided) ──────────────────────────────
cliff_wilcox <- function(hub, non) {
  n1 <- length(hub); n2 <- length(non)
  w  <- wilcox.test(hub, non, exact = FALSE)
  U  <- as.numeric(w$statistic)
  list(delta = (2 * U / (n1 * n2)) - 1, p = w$p.value, n1 = n1)
}

QUANTILES <- seq(0.10, 0.95, by = 0.05)

sens <- map_dfr(names(net_colors), function(net) {
  d <- filter(gene_tbl, network == net)
  map_dfr(QUANTILES, function(q) {
    thr <- quantile(d$degree, q)
    r <- cliff_wilcox(d$mean_dep[d$degree >= thr], d$mean_dep[d$degree < thr])
    tibble(network = net, quantile = q, top_pct = (1 - q) * 100,
           k_thr = as.numeric(thr), n_hub = r$n1,
           cliff_delta = r$delta, wilcoxon_p = r$p)
  })
}) %>%
  mutate(
    cliff_mag = case_when(
      abs(cliff_delta) < 0.147 ~ "negligible",
      abs(cliff_delta) < 0.330 ~ "small",
      abs(cliff_delta) < 0.474 ~ "medium",
      TRUE                     ~ "large"
    ),
    network = factor(network, levels = c("HPV-A7", "HPV-A9"))
  )

vroom_write(sens, file.path(dep_dir, "6_4_hub_threshold_sensitivity.tsv"), delim = "\t")
message("[OK] Sensitivity table saved.")

max_p   <- max(sens$wilcoxon_p, na.rm = TRUE)
p_label <- if (max_p < 0.001) "all comparisons  p < 0.001" else
  sprintf("max p = %.3g", max_p)
print(sens %>% dplyr::select(network, top_pct, k_thr, n_hub, cliff_delta, cliff_mag, wilcoxon_p),
      n = Inf)

## ── Plot ──────────────────────────────────────────────────────────────────────
base_theme <- theme_bw(base_family = "Times New Roman") +
  theme(
    axis.text        = element_text(size = 18, face = "bold"),
    axis.title       = element_text(size = 21, face = "bold"),
    plot.title       = element_text(size = 24, face = "bold", hjust = 0.5),
    legend.text      = element_text(size = 17, face = "bold"),
    legend.title     = element_blank(),
    legend.position  = c(0.83, 0.88),
    legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", linewidth = 0.8)
  )

mag_lines <- tibble(y = c(-0.147, -0.330, -0.474),
                    lab = c("small", "medium", "large"))

p_sens <- ggplot(sens, aes(x = top_pct, y = cliff_delta, color = network)) +
  geom_hline(data = mag_lines, aes(yintercept = y),
             linetype = "dotted", color = "grey60", linewidth = 0.7) +
  geom_text(data = mag_lines, aes(x = 90, y = y, label = lab), inherit.aes = FALSE,
            hjust = 0, vjust = -0.4, size = 4.6, fontface = "italic",
            family = "Times New Roman", color = "grey45") +
  geom_vline(xintercept = 10, linetype = "dashed", color = "black", linewidth = 0.8) +
  annotate("text", x = 10, y = -0.02, label = "chosen\n(top 10%)", size = 4.6,
           fontface = "bold.italic", family = "Times New Roman",
           color = "grey25", lineheight = 0.85, hjust = -0.08) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 4) +
  scale_color_manual(values = net_colors) +
  scale_x_reverse(breaks = seq(90, 10, by = -10), labels = function(x) paste0(x, "%")) +
  scale_y_continuous(limits = c(-0.5, 0), breaks = seq(-0.5, 0, by = 0.1)) +
  labs(title = "Effect size is robust across hub definitions",
       x = "Hub threshold (top % by degree)  →  more stringent",
       y = "Cliff's δ  (hub vs non-hub)") +
  annotate("text", x = 50, y = -0.49, label = p_label,
           size = 5, fontface = "bold.italic", family = "Times New Roman") +
  base_theme

ggsave(file.path(dep_dir, "6_4_A_hub_threshold_sensitivity.png"),
       p_sens, width = 12, height = 8, dpi = 300, bg = "white")
message("✅ 6_4 outputs saved to: ", dep_dir)
