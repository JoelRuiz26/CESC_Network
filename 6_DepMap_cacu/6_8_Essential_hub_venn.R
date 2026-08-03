### ============================================================
### 6_8_Essential_hub_venn.R
###
### Two-set Venn of the ESSENTIAL HUB genes of each clade.
###
### essential hub = hub in that clade's network AND strongly required in
### that clade's cell lines. Both thresholds are taken from the data with
### the rule of 6_10: genes are ranked by degree and by mean gene effect
### within each network and the cut is the elbow of the ranked curve,
### the point of maximum distance to the chord joining its two ends.
###   HPV-A7  degree >= 35   gene effect <= -0.317   ->  395 genes
###   HPV-A9  degree >= 42   gene effect <= -0.316   ->  436 genes
###
### The set totals are printed in the set names because a Venn draws the
### region counts, and 227 + 168 = 395 is not obvious from the regions.
###
### Same palette as 6_9, 6_10 and 6_11.
###   6_8_Essential_hub_venn.png/.pdf
###   6_8_essential_hub_venn_lists.tsv
### ============================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(vroom)
  library(ggvenn)
})

res_dir <- "~/CESC_Network/6_DepMap_cacu/results"
FF     <- "Times New Roman"
col_A7 <- "#4E79A7"
col_A9 <- "#CB6A5C"

## ── Data-driven thresholds: elbow of each ranked curve ───────────────────────
elbow <- function(v) {
  n <- length(v)
  x <- (seq_len(n) - 1) / (n - 1)
  y <- (v - v[1]) / (v[n] - v[1])
  v[which.max(abs(y - x))]
}

tab <- vroom(file.path(res_dir, "6_2_gene_dep_degree_table.tsv"), show_col_types = FALSE) %>%
  filter(!is.na(degree), !is.na(mean_dep)) %>%
  distinct(network, gene, .keep_all = TRUE)

th <- tab %>%
  group_by(network) %>%
  summarise(k_cut = elbow(sort(degree, decreasing = TRUE)),
            d_cut = elbow(sort(mean_dep)), .groups = "drop")
print(as.data.frame(th))

ess <- tab %>%
  left_join(th, by = "network") %>%
  filter(degree >= k_cut, mean_dep <= d_cut) %>%
  mutate(gene = toupper(gene))

ess_hub_A7 <- sort(unique(ess$gene[ess$network == "HPV-A7"]))
ess_hub_A9 <- sort(unique(ess$gene[ess$network == "HPV-A9"]))

shared <- intersect(ess_hub_A7, ess_hub_A9)
onlyA7 <- setdiff(ess_hub_A7, ess_hub_A9)
onlyA9 <- setdiff(ess_hub_A9, ess_hub_A7)

cat(sprintf("\nEssential hub genes (data-driven thresholds)\n  A7 total : %d\n  A9 total : %d\n",
            length(ess_hub_A7), length(ess_hub_A9)))
cat(sprintf("  Shared (intersection) : %d\n  A7-only               : %d\n  A9-only               : %d\n",
            length(shared), length(onlyA7), length(onlyA9)))

## the three gene lists, padded to equal length for a tidy TSV
n   <- max(length(shared), length(onlyA7), length(onlyA9))
pad <- function(x) c(x, rep(NA_character_, n - length(x)))
vroom_write(tibble(shared = pad(shared), A7_only = pad(onlyA7), A9_only = pad(onlyA9)),
            file.path(res_dir, "6_8_essential_hub_venn_lists.tsv"), delim = "\t")

## ── Venn ──────────────────────────────────────────────────────────────────────
venn_list <- list(ess_hub_A7, ess_hub_A9)
names(venn_list) <- c(sprintf("HPV-A7\n(n = %d)", length(ess_hub_A7)),
                      sprintf("HPV-A9\n(n = %d)", length(ess_hub_A9)))

p <- ggvenn(venn_list,
            fill_color   = c(col_A7, col_A9), fill_alpha = 0.5,
            stroke_size  = 0.7, set_name_size = 11, text_size = 9,
            show_percentage = FALSE, digits = 0) +
  coord_fixed(clip = "off") +
  theme(text = element_text(family = FF),
        plot.margin = margin(20, 20, 20, 20))

ggsave(file.path(res_dir, "6_8_Essential_hub_venn.png"), p,
       width = 9, height = 8, dpi = 400, bg = "white")
ggsave(file.path(res_dir, "6_8_Essential_hub_venn.pdf"), p,
       width = 9, height = 8, device = cairo_pdf, bg = "white")
message("✅ saved 6_8 Venn + gene lists to: ", res_dir)
