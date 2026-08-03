### ============================================================
### 6_6_Module_enrichment_stability.R
###
### Sensitivity of the 6_5 module enrichment to its two arbitrary
### parameters, swept jointly:
###   hub definition   top 20 / 15 / 10 / 5 % of degree
###   essentiality     dependency <= -0.3 ... -1.3
### 44 parameter combinations per network.
###
### Everything else is held exactly as in 6_5: pan-essential genes
### kept, only Hallmark-annotated modules tested, hypergeometric
### against the full network, BH within each (network, hub, cutoff).
###
### Outputs:
###   6_6_module_enrichment_stability.tsv — every module x cell test
###   6_6_A_module_stability.png/.pdf     — evidence across the grid
### ============================================================
suppressPackageStartupMessages({
  library(igraph)
  library(tidyverse)
  library(vroom)
})

dep_dir <- "~/CESC_Network/6_DepMap_cacu"
res_dir <- file.path(dep_dir, "results")
mod_dir <- "~/CESC_Network/7_Modularity/5_1_Modules"
ora_dir <- "~/CESC_Network/7_Modularity/5_2_ORA"

HUB_QUANTILES <- c(0.80, 0.85, 0.90, 0.95)   # top 20 / 15 / 10 / 5 %
CUTOFFS       <- seq(-0.3, -1.3, by = -0.1)
REF_HUB       <- 0.90                        # the values used in 6_5
REF_CUTOFF    <- -0.5
ALPHA         <- 0.05

## ── Data ──────────────────────────────────────────────────────────────────────
mod_tbl <- function(g) tibble(gene = V(g)$name, community_name = V(g)$Community_name)

g7 <- readRDS(file.path(mod_dir, "5_1_1_Louvain_A7_graph_max.rds"))
g9 <- readRDS(file.path(mod_dir, "5_1_1_Louvain_A9_graph_max.rds"))
deg7 <- igraph::degree(g7); deg9 <- igraph::degree(g9)
m7 <- mod_tbl(g7); m9 <- mod_tbl(g9)

d7 <- vroom(file.path(res_dir, "6_1_cervix_dep_mean_A7.tsv"), show_col_types = FALSE)
d9 <- vroom(file.path(res_dir, "6_1_cervix_dep_mean_A9.tsv"), show_col_types = FALSE)

hallmark_modules <- function(path)
  vroom(path, show_col_types = FALSE) %>% distinct(Community_Name) %>% pull(Community_Name)
hm7 <- hallmark_modules(file.path(ora_dir, "5_2_3_Hallmark7_louv_max.tsv"))
hm9 <- hallmark_modules(file.path(ora_dir, "5_2_3_Hallmark9_louv_max.tsv"))

## ── One cell of the grid ──────────────────────────────────────────────────────
# identical to enrich_modules() in 6_5; hub quantile and cutoff are the free args
enrich_cell <- function(deg, m, dep, keep_modules, group, hq, cut) {
  hubs  <- names(deg)[deg >= quantile(deg, hq)]
  ess   <- dep$gene_u[dep$mean_dep <= cut]
  genes <- hubs[toupper(hubs) %in% ess]
  N <- nrow(m); K <- sum(m$gene %in% genes)
  m %>% group_by(community_name) %>%
    summarise(n_M = n(), k = sum(gene %in% genes), .groups = "drop") %>%
    filter(community_name %in% keep_modules) %>%
    mutate(expected = n_M * K / N,
           fold     = k / expected,
           p        = phyper(k - 1, K, N - K, n_M, lower.tail = FALSE),
           q        = p.adjust(p, "BH"),
           network  = group, hub_q = hq, cutoff = cut,
           n_hub = length(hubs), K_total = K)
}

grid <- expand_grid(hq = HUB_QUANTILES, cut = CUTOFFS)
sweep <- bind_rows(lapply(seq_len(nrow(grid)), function(i) bind_rows(
  enrich_cell(deg7, m7, d7, hm7, "HPV-A7", grid$hq[i], grid$cut[i]),
  enrich_cell(deg9, m9, d9, hm9, "HPV-A9", grid$hq[i], grid$cut[i]))))
vroom_write(sweep, file.path(res_dir, "6_6_module_enrichment_stability.tsv"), delim = "\t")

## ── Stability summary ─────────────────────────────────────────────────────────
stability <- sweep %>% group_by(network, community_name) %>%
  summarise(n_sig = sum(q < ALPHA), n_cells = n(), frac_sig = n_sig / n_cells,
            sig_at_ref = any(q < ALPHA & hub_q == REF_HUB & cutoff == REF_CUTOFF),
            fold_min = min(fold), fold_max = max(fold), .groups = "drop") %>%
  arrange(network, desc(frac_sig))
cat("\n=== Modules significant in at least one of the 44 cells ===\n")
stability %>% filter(n_sig > 0) %>%
  mutate(across(c(frac_sig, fold_min, fold_max), ~round(.x, 2))) %>%
  as.data.frame() %>% print()

## ── Figure ────────────────────────────────────────────────────────────────────
## Lines, not tiles: the distance to the q = 0.05 line carries how much evidence
## a module has, which a significant/not-significant tile throws away.
base_theme <- theme_bw(base_family = "Times New Roman") +
  theme(axis.text    = element_text(size = 14, face = "bold"),
        axis.title   = element_text(size = 19, face = "bold"),
        plot.title   = element_text(size = 23, face = "bold", hjust = 0.5),
        plot.subtitle= element_text(size = 14, face = "italic", hjust = 0.5, color = "grey30"),
        strip.text   = element_text(size = 16, face = "bold"),
        strip.background = element_rect(fill = "grey93", color = "grey60"),
        legend.text  = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 15, face = "bold"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 0.8))

mod_pal <- c("M2" = "#67001F", "M6" = "#D6604D", "M8" = "#F4A582",
             "M7" = "#4393C3", "M12" = "#2166AC", "M13" = "#7F7F7F",
             "M9" = "#B8B8B8", "M18" = "#9C9C9C", "M22" = "#5F5F5F")

# only modules that reach significance somewhere in the grid
keep <- stability %>% filter(n_sig > 0) %>% distinct(network, community_name)
plt <- sweep %>% inner_join(keep, by = c("network", "community_name")) %>%
  mutate(hub_lab = factor(sprintf("Hub = top %g%%", round((1 - hub_q) * 100)),
                          levels = sprintf("Hub = top %g%%",
                                           round((1 - sort(HUB_QUANTILES, TRUE)) * 100))),
         module = factor(community_name, levels = names(mod_pal)),
         neglogq = -log10(q))

pA <- ggplot(plt, aes(abs(cutoff), neglogq, color = module)) +
  geom_hline(yintercept = -log10(ALPHA), linetype = "dashed",
             color = "grey30", linewidth = 0.8) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.0) +
  facet_grid(network ~ hub_lab, scales = "free_y") +
  scale_color_manual(values = mod_pal, name = "Module", drop = TRUE) +
  scale_x_continuous(breaks = seq(0.3, 1.3, by = 0.2)) +
  # sqrt on an already-log axis: without it a module at -log10 q = 60 squashes
  # the 0-5 band, which is exactly where modules cross q = 0.05
  scale_y_sqrt(breaks = c(1.3, 5, 10, 20, 40, 60),
               labels = c("1.3\n(q=0.05)", "5", "10", "20", "40", "60")) +
  labs(title = "Sensitivity of the module enrichment to both thresholds",
       subtitle = "Dashed line: q = 0.05. Points above it are significant.",
       x = "Essentiality cutoff  (|dependency|, stricter to the right)",
       y = expression(bold(-log[10]~italic(q)))) +
  base_theme
ggsave(file.path(res_dir, "6_6_A_module_stability.png"), pA,
       width = 16, height = 10, dpi = 300, bg = "white")
ggsave(file.path(res_dir, "6_6_A_module_stability.pdf"), pA,
       width = 16, height = 10, device = cairo_pdf, bg = "white")

## ── Calibrating the stability score ───────────────────────────────────────────
## "Significant in >= 80% of cells" would just move the arbitrariness elsewhere.
## Instead the module labels are permuted (module SIZES and every dependency
## value are left untouched, only the gene -> module assignment is broken) and
## the whole 44-cell grid is recomputed. That gives, for a module of each size,
## the retention frequency the grid produces by chance. The cut comes from the
## null, not from us.
B <- 2000
set.seed(1)

null_retention <- function(deg, m, dep, keep_modules) {
  N <- nrow(m); mods <- keep_modules; nm <- length(mods)
  n_M <- as.integer(table(factor(m$community_name, levels = mods)))
  cells <- lapply(seq_len(nrow(grid)), function(i) {
    hubs  <- names(deg)[deg >= quantile(deg, grid$hq[i])]
    ess   <- dep$gene_u[dep$mean_dep <= grid$cut[i]]
    which(m$gene %in% hubs[toupper(hubs) %in% ess])
  })
  Ks <- lengths(cells)
  out <- matrix(0, B, nm, dimnames = list(NULL, mods))
  for (b in seq_len(B)) {
    labi <- match(sample(m$community_name), mods)   # sizes preserved
    sig  <- matrix(FALSE, length(cells), nm)
    for (i in seq_along(cells)) {
      k <- tabulate(labi[cells[[i]]], nbins = nm)
      p <- phyper(k - 1L, Ks[i], N - Ks[i], n_M, lower.tail = FALSE)
      sig[i, ] <- p.adjust(p, "BH") < ALPHA
    }
    out[b, ] <- colMeans(sig)
  }
  out
}

null7 <- null_retention(deg7, m7, d7, hm7)
null9 <- null_retention(deg9, m9, d9, hm9)

# pivot each network BEFORE binding: the two networks have different module
# names, so binding the wide matrices would pad the gaps with NA
long_null <- function(mat, net) as_tibble(mat) %>%
  pivot_longer(everything(), names_to = "community_name", values_to = "frac_null") %>%
  mutate(network = net)
null_tbl <- bind_rows(long_null(null7, "HPV-A7"), long_null(null9, "HPV-A9"))

stab <- stability %>%
  left_join(null_tbl %>% group_by(network, community_name) %>%
              summarise(null_max = max(frac_null),
                        null_q95 = quantile(frac_null, 0.95),
                        p_perm   = NA_real_, .groups = "drop"),
            by = c("network", "community_name")) %>%
  left_join(null_tbl %>% inner_join(stability %>% select(network, community_name, frac_sig),
                                    by = c("network", "community_name")) %>%
              group_by(network, community_name) %>%
              summarise(p_perm = (1 + sum(frac_null >= first(frac_sig))) / (B + 1),
                        .groups = "drop"),
            by = c("network", "community_name"), suffix = c("_drop", "")) %>%
  select(-p_perm_drop) %>%
  # q95, not the max: the max of B draws is a single extreme realisation and it
  # reordered modules that the retention frequency itself ranks the other way
  mutate(retained = frac_sig > null_q95)
vroom_write(stab, file.path(res_dir, "6_6_module_stability_summary.tsv"), delim = "\t")

cat("\n=== Retention frequency vs permuted null (", B, " label shuffles) ===\n")
stab %>% filter(n_sig > 0) %>%
  transmute(network, module = community_name, frac_sig = round(frac_sig, 2),
            null_max = round(null_max, 2), p_perm, retained) %>%
  as.data.frame() %>% print()

## ── Figure B: retention frequency ─────────────────────────────────────────────
## The permutation null above is reported in the TSV but is NOT used to draw a
## line here: its 95th percentile is 0 for every module (a randomly labelled
## module is essentially never significant in any cell), so it would retain a
## module that survives 2% of the grid. Retention across analyst choices is a
## descriptive question, not an inferential one. The two criteria used are:
##   - the largest gap in the observed retention values (data-driven, not chosen)
##   - significance at the pre-specified reference point (-0.5, top 10%)
gap_line <- stab %>% filter(n_sig > 0) %>% group_by(network) %>%
  arrange(desc(frac_sig), .by_group = TRUE) %>%
  mutate(gap = frac_sig - lead(frac_sig)) %>%
  summarise(cut_at = { i <- which.max(gap); mean(c(frac_sig[i], frac_sig[i + 1])) },
            .groups = "drop")

plt_b <- stab %>% filter(n_sig > 0) %>%
  mutate(module = factor(community_name, levels = names(mod_pal)))

pB <- ggplot(plt_b, aes(frac_sig, fct_reorder(community_name, frac_sig))) +
  geom_vline(data = gap_line, aes(xintercept = cut_at), linetype = "dashed",
             color = "grey35", linewidth = 0.9) +
  geom_segment(aes(x = 0, xend = frac_sig, yend = fct_reorder(community_name, frac_sig)),
               color = "grey75", linewidth = 1.2) +
  geom_point(aes(color = module, shape = sig_at_ref), size = 6.5, stroke = 1.5) +
  facet_wrap(~network, nrow = 1, scales = "free_y") +
  scale_color_manual(values = mod_pal, guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 19, `FALSE` = 1),
                     labels = c(`TRUE` = "Significant at the reference parameters",
                                `FALSE` = "Not significant at the reference parameters"),
                     name = NULL) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1.04),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Which modules survive the parameter grid",
       subtitle = "Dashed line: largest gap in the observed retention values",
       x = "Cells of the 44-point grid where the module stays significant",
       y = NULL) +
  base_theme + theme(legend.position = "bottom")
ggsave(file.path(res_dir, "6_6_B_stability_vs_null.png"), pB,
       width = 14, height = 7, dpi = 300, bg = "white")
ggsave(file.path(res_dir, "6_6_B_stability_vs_null.pdf"), pB,
       width = 14, height = 7, device = cairo_pdf, bg = "white")

message("✅ 6_6 outputs saved to: ", res_dir)
