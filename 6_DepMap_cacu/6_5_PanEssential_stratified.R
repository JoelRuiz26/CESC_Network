### ============================================================
### 6_5_PanEssential_stratified.R
###
### Robustness check for the hub -> dependency result (6_3 / 6_4):
### is the signal driven by pan-essential (housekeeping) genes, or
### does it persist among context-selective dependencies?
###
### Pan-essential = DepMap official "common_essential" flag from
### depmap_gene_summary(), dataset "Chronos_Achilles" (same Achilles
### CRISPR screen used in 6_1). Release 22Q2. NOT re-derived by us.
###
### The hub definition is IDENTICAL to 6_3 (top 10% by degree on the
### full dependency-data pool). We then recompute Wilcoxon + Cliff's
### delta after REMOVING pan-essential genes from BOTH groups.
###
### Two outcomes, both valid:
###   * |delta| persists  -> hubs capture SELECTIVE dependency (strong)
###   * |delta| collapses -> signal was housekeeping (weaker story)
###
### Does NOT modify 6_3 / 6_4.
###
### Outputs:
###   6_5_common_essential_genes.tsv     — pan-essential set used (provenance)
###   6_5_hub_panessential_composition.tsv — are hubs enriched for pan-ess?
###   6_5_stratified_stats.tsv           — delta: all vs excl pan-essential
###   6_5_sensitivity_stratified.tsv     — hub-threshold sweep, both strata
###   6_5_A_delta_comparison.png         — dumbbell all vs excl pan-essential
###   6_5_B_sensitivity_stratified.png   — sensitivity curves, both strata
### ============================================================
suppressPackageStartupMessages({
  library(depmap)
  library(tidyverse)
  library(vroom)
  library(patchwork)
})

dep_dir <- "~/CESC_Network/6_DepMap_cacu"
HUB_QUANTILE <- 0.90
net_colors <- c("HPV-A7" = "#003366", "HPV-A9" = "#CC2200")

## ── 1) Pan-essential set (DepMap official flag) ───────────────────────────────
gs <- depmap_gene_summary()
common_ess <- gs %>%
  filter(dataset == "Chronos_Achilles", common_essential) %>%
  transmute(gene_u = toupper(gene_name)) %>%
  distinct()

vroom_write(common_ess, file.path(dep_dir, "6_5_common_essential_genes.tsv"), delim = "\t")
message("Pan-essential genes (Chronos_Achilles): ", nrow(common_ess))

pan_set <- common_ess$gene_u

## ── 2) Network gene table + pan-essential label ───────────────────────────────
gene_tbl <- vroom(file.path(dep_dir, "6_2_gene_dep_degree_table.tsv"),
                  show_col_types = FALSE) %>%
  mutate(pan_essential = toupper(gene) %in% pan_set)

## ── Helper: Cliff's delta via Mann-Whitney U (two-sided) ──────────────────────
cliff_wilcox <- function(hub, non) {
  n1 <- length(hub); n2 <- length(non)
  if (n1 < 3 || n2 < 3) return(list(delta = NA_real_, p = NA_real_, n1 = n1, n2 = n2))
  w <- wilcox.test(hub, non, exact = FALSE)
  U <- as.numeric(w$statistic)
  list(delta = (2 * U / (n1 * n2)) - 1, p = w$p.value, n1 = n1, n2 = n2)
}
cliff_mag <- function(d) case_when(
  is.na(d)       ~ NA_character_,
  abs(d) < 0.147 ~ "negligible",
  abs(d) < 0.330 ~ "small",
  abs(d) < 0.474 ~ "medium",
  TRUE           ~ "large"
)

## ── 3) Classify hubs (same as 6_3), then stratify ─────────────────────────────
classify <- function(df) {
  thr <- quantile(df$degree, HUB_QUANTILE)
  df %>% mutate(gene_class = if_else(degree >= thr, "Hub", "Non-hub"), hub_thr = thr)
}

## Composition: are hubs enriched for pan-essential genes? (Fisher)
composition <- map_dfr(names(net_colors), function(net) {
  d  <- classify(filter(gene_tbl, network == net))
  ft <- fisher.test(table(Hub = d$gene_class == "Hub", Pan = d$pan_essential))
  tibble(
    network         = net,
    pct_pan_hub     = mean(d$pan_essential[d$gene_class == "Hub"]) * 100,
    pct_pan_nonhub  = mean(d$pan_essential[d$gene_class == "Non-hub"]) * 100,
    OR_hub_pan      = unname(ft$estimate),
    p               = ft$p.value
  )
})
vroom_write(composition, file.path(dep_dir, "6_5_hub_panessential_composition.tsv"), delim = "\t")
message("\n[Composition] fraction of hubs that are pan-essential:")
print(composition)

## Stratified delta: all genes vs excluding pan-essential
stratified <- map_dfr(names(net_colors), function(net) {
  d <- classify(filter(gene_tbl, network == net))
  map_dfr(c("All genes", "Excl. pan-essential"), function(stratum) {
    dd <- if (stratum == "All genes") d else filter(d, !pan_essential)
    r  <- cliff_wilcox(dd$mean_dep[dd$gene_class == "Hub"],
                       dd$mean_dep[dd$gene_class == "Non-hub"])
    tibble(network = net, stratum = stratum,
           n_hub = r$n1, n_non = r$n2,
           cliff_delta = r$delta, cliff_mag = cliff_mag(r$delta), wilcoxon_p = r$p)
  })
}) %>% mutate(network = factor(network, levels = names(net_colors)),
              stratum = factor(stratum, levels = c("All genes", "Excl. pan-essential")))

vroom_write(stratified, file.path(dep_dir, "6_5_stratified_stats.tsv"), delim = "\t")
message("\n[Stratified] Cliff's delta — all vs excl pan-essential:")
print(stratified)

## ── 4) Sensitivity sweep, both strata ─────────────────────────────────────────
QUANTILES <- seq(0.10, 0.95, by = 0.05)
sens <- map_dfr(names(net_colors), function(net) {
  d0 <- filter(gene_tbl, network == net)
  map_dfr(QUANTILES, function(q) {
    thr <- quantile(d0$degree, q)
    d   <- d0 %>% mutate(gene_class = if_else(degree >= thr, "Hub", "Non-hub"))
    map_dfr(c("All genes", "Excl. pan-essential"), function(stratum) {
      dd <- if (stratum == "All genes") d else filter(d, !pan_essential)
      r  <- cliff_wilcox(dd$mean_dep[dd$gene_class == "Hub"],
                         dd$mean_dep[dd$gene_class == "Non-hub"])
      tibble(network = net, stratum = stratum, top_pct = (1 - q) * 100,
             cliff_delta = r$delta, wilcoxon_p = r$p)
    })
  })
}) %>% mutate(network = factor(network, levels = names(net_colors)),
              stratum = factor(stratum, levels = c("All genes", "Excl. pan-essential")))

vroom_write(sens, file.path(dep_dir, "6_5_sensitivity_stratified.tsv"), delim = "\t")

## ── Theme ─────────────────────────────────────────────────────────────────────
base_theme <- theme_bw(base_family = "Times New Roman") +
  theme(
    axis.text        = element_text(size = 17, face = "bold"),
    axis.title       = element_text(size = 20, face = "bold"),
    plot.title       = element_text(size = 23, face = "bold", hjust = 0.5),
    legend.text      = element_text(size = 16, face = "bold"),
    legend.title     = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", linewidth = 0.8)
  )
mag_lines <- tibble(y = c(-0.147, -0.330, -0.474), lab = c("small", "medium", "large"))

## ── Plot A: dumbbell all vs excl pan-essential ────────────────────────────────
wide <- stratified %>%
  select(network, stratum, cliff_delta) %>%
  pivot_wider(names_from = stratum, values_from = cliff_delta) %>%
  rename(all = "All genes", excl = "Excl. pan-essential") %>%
  mutate(network = factor(network, levels = rev(names(net_colors))))

pA <- ggplot(wide) +
  geom_hline(data = mag_lines, aes(yintercept = y), linetype = "dotted",
             color = "grey65", linewidth = 0.7) +
  geom_segment(aes(x = network, xend = network, y = all, yend = excl),
               color = "grey70", linewidth = 2.5) +
  geom_point(aes(x = network, y = all, color = network), size = 9) +
  geom_point(aes(x = network, y = excl, color = network), size = 9,
             shape = 21, fill = "white", stroke = 2.2) +
  geom_text(aes(x = network, y = all,  label = sprintf("%.2f", all)),
            vjust = -1.6, size = 5, fontface = "bold", family = "Times New Roman") +
  geom_text(aes(x = network, y = excl, label = sprintf("%.2f", excl)),
            vjust = -1.6, size = 5, fontface = "bold", family = "Times New Roman") +
  scale_color_manual(values = net_colors, guide = "none") +
  scale_y_continuous(limits = c(-0.5, 0), breaks = seq(-0.5, 0, 0.1)) +
  coord_flip() +
  labs(title = "Effect persists after removing pan-essential genes",
       subtitle = "Filled = all genes   |   Hollow = excl. pan-essential",
       x = NULL, y = "Cliff's δ  (hub vs non-hub)") +
  base_theme +
  theme(plot.subtitle = element_text(size = 15, face = "bold.italic",
                                      hjust = 0.5, color = "grey30"))

ggsave(file.path(dep_dir, "6_5_A_delta_comparison.png"),
       pA, width = 12, height = 6, dpi = 300, bg = "white")
message("[OK] Plot A saved.")

## ── Plot B: sensitivity curves, both strata ───────────────────────────────────
pB <- ggplot(sens, aes(x = top_pct, y = cliff_delta,
                       color = network, linetype = stratum)) +
  geom_hline(data = mag_lines, aes(yintercept = y), linetype = "dotted",
             color = "grey65", linewidth = 0.6) +
  geom_text(data = mag_lines, aes(x = 90, y = y, label = lab), inherit.aes = FALSE,
            hjust = 0, vjust = -0.4, size = 4.2, fontface = "italic",
            family = "Times New Roman", color = "grey45") +
  geom_vline(xintercept = 10, linetype = "dashed", color = "black", linewidth = 0.7) +
  geom_line(linewidth = 1.4) +
  geom_point(size = 3.5) +
  scale_color_manual(values = net_colors) +
  scale_linetype_manual(values = c("All genes" = "solid", "Excl. pan-essential" = "22")) +
  scale_x_reverse(breaks = seq(90, 10, by = -10), labels = function(x) paste0(x, "%")) +
  labs(title = "Robustness across hub thresholds and strata",
       x = "Hub threshold (top % by degree)  →  more stringent",
       y = "Cliff's δ  (hub vs non-hub)") +
  base_theme +
  theme(legend.position = c(0.80, 0.20),
        legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
        legend.spacing.y = unit(2, "pt"))

ggsave(file.path(dep_dir, "6_5_B_sensitivity_stratified.png"),
       pB, width = 12, height = 8, dpi = 300, bg = "white")
message("✅ 6_5 outputs saved to: ", dep_dir)
