### ============================================================
### 6_9_Essential_hub_scatter.R
###
### Gene effect scatter of ALL hub genes (union of the A7 and A9 hub
### sets, n = 2534) shown as a grey cloud, with the ESSENTIAL hubs
### highlighted in colour:
###   grey   1871   hub, not an essential hub in either clade
###   blue    227   essential hub in A7 only
###   red     268   essential hub in A9 only
###   purple  168   essential hub in BOTH clades
###
### Thresholds are the data-driven ones of 6_10: genes are ranked by
### degree and by mean gene effect within each network and the cut is
### the elbow of each ranked curve.
###   A7  hub: degree >= 35   essential: gene effect <= -0.317
###   A9  hub: degree >= 42   essential: gene effect <= -0.316
###
### axes = mean CRISPR gene effect in the A7 (x) and A9 (y) lines.
### The grey cloud is the full hub universe; the colour marks which
### clade(s) the gene is an essential hub in. Clade-unique hubs (blue /
### red) are NOT in the essential-hub set of the other network.
###
### Gene NAMES are shown only for the non pan-essential essential hubs;
### pan-essential housekeeping genes stay coloured but unlabelled.
### Pan-essential set = 6_1_common_essential_genes.tsv (Chronos_Combined).
###
###   6_9_A_essential_hub_scatter.png/.pdf
### ============================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(vroom)
  library(ggrepel)
})

res_dir <- "~/CESC_Network/6_DepMap_cacu/results"

lvl  <- c("Hub, not essential", "Essential hub in A7 only",
          "Essential hub in A9 only", "Essential hub in both")
cols <- c("Hub, not essential"       = "grey87",
          "Essential hub in A7 only" = "#4E79A7",   # muted blue       (A7)
          "Essential hub in A9 only" = "#CB6A5C",   # muted terracotta (A9)
          "Essential hub in both"    = "#8C79B5")   # muted purple     (both)

## ── Data-driven thresholds: elbow of each ranked curve (same rule as 6_10) ─────
elbow <- function(v) {
  n <- length(v)
  x <- (seq_len(n) - 1) / (n - 1)
  y <- (v - v[1]) / (v[n] - v[1])
  v[which.max(abs(y - x))]
}

net <- vroom(file.path(res_dir, "6_2_gene_dep_degree_table.tsv"), show_col_types = FALSE) %>%
  filter(!is.na(degree), !is.na(mean_dep)) %>%
  distinct(network, gene, .keep_all = TRUE)

th <- net %>%
  group_by(network) %>%
  summarise(k_cut = elbow(sort(degree, decreasing = TRUE)),
            d_cut = elbow(sort(mean_dep)), .groups = "drop")
print(as.data.frame(th))

flags <- net %>%
  left_join(th, by = "network") %>%
  mutate(hub = degree >= k_cut, ess = hub & mean_dep <= d_cut,
         net = str_remove(network, "HPV-")) %>%
  select(gene, net, hub, ess) %>%
  pivot_wider(names_from = net, values_from = c(hub, ess)) %>%
  mutate(across(starts_with("hub_") | starts_with("ess_"), ~ replace_na(.x, FALSE)))

## ── Clade mean gene effect per gene, for the two axes ─────────────────────────
dat <- vroom(file.path(res_dir, "6_7_hub_dependency_medians.tsv"), show_col_types = FALSE) %>%
  select(gene, gene_u, dep_A7, dep_A9) %>%
  inner_join(flags, by = "gene") %>%
  filter(hub_A7 | hub_A9) %>%                      # the hub universe
  mutate(grp = factor(case_when(
    ess_A7 & ess_A9  ~ "Essential hub in both",
    ess_A7           ~ "Essential hub in A7 only",
    ess_A9           ~ "Essential hub in A9 only",
    TRUE             ~ "Hub, not essential"), levels = lvl))

x_cut <- th$d_cut[th$network == "HPV-A7"]
y_cut <- th$d_cut[th$network == "HPV-A9"]

## the legend title already says "Essential hub", so the keys only carry the clade
short <- c("Essential hub in both"    = "Both",
           "Essential hub in A7 only" = "HPV-A7 only",
           "Essential hub in A9 only" = "HPV-A9 only")

grp_n     <- dat %>% count(grp)
grp_labs  <- setNames(sprintf("%s  (n = %d)", short[as.character(grp_n$grp)], grp_n$n), grp_n$grp)
leg_order <- c("Essential hub in both", "Essential hub in A7 only", "Essential hub in A9 only")

## Names are shown ONLY for non pan-essential essential hubs, so the labels
## flag the biologically selective dependencies. Colouring is not affected.
pan <- vroom(file.path(res_dir, "6_1_common_essential_genes.tsv"),
             delim = "\t", show_col_types = FALSE)$gene_u %>% toupper()
sel <- dat %>% filter(grp != "Hub, not essential", !(toupper(gene_u) %in% pan))
TOP_N <- 5
hi <- bind_rows(
  sel %>% filter(grp == "Essential hub in both")    %>% slice_min(dep_A7 + dep_A9, n = TOP_N),
  sel %>% filter(grp == "Essential hub in A7 only") %>% slice_min(dep_A7, n = TOP_N),
  sel %>% filter(grp == "Essential hub in A9 only") %>% slice_min(dep_A9, n = TOP_N)
) %>% distinct(gene, .keep_all = TRUE)

rng <- range(c(dat$dep_A7, dat$dep_A9), na.rm = TRUE)
lim <- c(rng[1] - 0.08, rng[2] + 0.1)

## ── Theme ─────────────────────────────────────────────────────────────────────
base_theme <- theme_bw(base_family = "Times New Roman") +
  theme(axis.text     = element_text(size = 17, face = "bold"),
        axis.title    = element_text(size = 21, face = "bold"),
        plot.title    = element_text(size = 24, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, face = "italic", hjust = 0.5, color = "grey30"),
        legend.text   = element_text(size = 14, face = "bold"),
        legend.title  = element_text(size = 16, face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey93", linewidth = 0.3),
        panel.border  = element_rect(color = "grey35", linewidth = 0.6))

## ── Scatter: grey cloud first, coloured essential hubs on top ─────────────────
grey_lyr <- dat %>% filter(grp == "Hub, not essential")
col_lyr  <- dat %>% filter(grp != "Hub, not essential") %>% droplevels()

p <- ggplot(mapping = aes(dep_A7, dep_A9)) +
  geom_abline(slope = 1, intercept = 0, color = "grey87", linewidth = 0.4) +
  geom_hline(yintercept = y_cut, linetype = "dashed", linewidth = 0.4, color = "grey62") +
  geom_vline(xintercept = x_cut, linetype = "dashed", linewidth = 0.4, color = "grey62") +
  geom_point(data = grey_lyr, color = "grey78", size = 2.0, alpha = 0.55, stroke = 0) +
  geom_point(data = col_lyr, aes(fill = grp), shape = 21, color = "white",
             size = 2.0, stroke = 0.20, alpha = 0.82) +
  geom_text_repel(data = hi, aes(label = gene), color = "grey20",
                  size = 3.2, fontface = "bold", family = "Times New Roman",
                  max.overlaps = Inf, min.segment.length = 0, box.padding = 0.75,
                  point.padding = 0.35, seed = 1, force = 6, force_pull = 0.4,
                  segment.color = "grey60", segment.size = 0.22,
                  bg.color = "white", bg.r = 0.12) +
  scale_fill_manual(values = cols, name = "Essential hub",
                    breaks = leg_order,
                    labels = grp_labs[leg_order]) +
  guides(fill = guide_legend(override.aes = list(size = 4.6, alpha = 1))) +
  scale_x_continuous(limits = lim) + scale_y_continuous(limits = lim) +
  coord_equal() +
  labs(title = "Hub gene dependencies across clades",
       x = "Mean CRISPR gene effect — HPV-A7 lines",
       y = "Mean CRISPR gene effect — HPV-A9 lines") +
  base_theme +
  theme(legend.position = c(0.985, 0.02),
        legend.justification = c("right", "bottom"),
        legend.background = element_rect(fill = alpha("white", 0.92), color = "grey75"),
        legend.key = element_blank())

ggsave(file.path(res_dir, "6_9_A_essential_hub_scatter.png"), p,
       width = 10, height = 10, dpi = 300, bg = "white")
ggsave(file.path(res_dir, "6_9_A_essential_hub_scatter.pdf"), p,
       width = 10, height = 10, device = cairo_pdf, bg = "white")
message("✅ hub-universe scatter with essential hubs highlighted saved to: ", res_dir)
print(as.data.frame(grp_n))

### ============================================================
### 6_9_B: the same scatter with every pan-essential gene removed,
### from the grey cloud as well as from the coloured points, so that
### what remains are only the clade-selective dependencies. Axes are
### fitted to the remaining data rather than to the full range.
###   6_9_B_essential_hub_scatter_nonpan.png/.pdf
### ============================================================
dat_np <- dat %>% filter(!(toupper(gene_u) %in% pan))

grp_n_np   <- dat_np %>% count(grp)
grp_labs_np <- setNames(sprintf("%s  (n = %d)", short[as.character(grp_n_np$grp)], grp_n_np$n),
                        grp_n_np$grp)

## labels: within each quadrant, the genes furthest from the intersection of the
## two thresholds, by euclidean distance. Both axes are gene effect, so the
## distance is in the units of the plot and matches what the eye reads. A summed
## percentile rank was tried first and rejected: it favours genes that are mildly
## extreme on both axes over genes that are strongly extreme on one, so it left
## CCNC and NHLRC2 unlabelled while labelling points closer to the crossing.
## The upper-right quadrant holds 1829 grey genes and no highlighted one, so it
## contributes no labels.
TOP_N <- 5
hi_np <- dat_np %>%
  filter(grp != "Hub, not essential") %>%
  mutate(quad = paste0(if_else(dep_A7 <= x_cut, "L", "R"),
                       if_else(dep_A9 <= y_cut, "B", "T")),
         dist = sqrt((dep_A7 - x_cut)^2 + (dep_A9 - y_cut)^2)) %>%
  group_by(quad) %>%
  slice_max(dist, n = TOP_N, with_ties = FALSE) %>%
  ungroup()

rng_np <- range(c(dat_np$dep_A7, dat_np$dep_A9), na.rm = TRUE)
lim_np <- c(rng_np[1] - 0.05, rng_np[2] + 0.05)

grey_np <- dat_np %>% filter(grp == "Hub, not essential")
col_np  <- dat_np %>% filter(grp != "Hub, not essential") %>% droplevels()

pB <- ggplot(mapping = aes(dep_A7, dep_A9)) +
  geom_abline(slope = 1, intercept = 0, color = "grey87", linewidth = 0.4) +
  geom_hline(yintercept = y_cut, linetype = "dashed", linewidth = 0.4, color = "grey62") +
  geom_vline(xintercept = x_cut, linetype = "dashed", linewidth = 0.4, color = "grey62") +
  geom_point(data = grey_np, color = "grey78", size = 2.4, alpha = 0.55, stroke = 0) +
  geom_point(data = col_np, aes(fill = grp), shape = 21, color = "white",
             size = 2.4, stroke = 0.22, alpha = 0.88) +
  geom_text_repel(data = hi_np, aes(label = gene), color = "grey20",
                  size = 3.4, fontface = "bold", family = "Times New Roman",
                  max.overlaps = Inf, min.segment.length = 0, box.padding = 0.75,
                  point.padding = 0.35, seed = 1, force = 6, force_pull = 0.4,
                  segment.color = "grey60", segment.size = 0.22,
                  bg.color = "white", bg.r = 0.12) +
  scale_fill_manual(values = cols, name = "Essential hub",
                    breaks = leg_order,
                    labels = grp_labs_np[leg_order]) +
  guides(fill = guide_legend(override.aes = list(size = 4.6, alpha = 1))) +
  scale_x_continuous(limits = lim_np) + scale_y_continuous(limits = lim_np) +
  coord_equal() +
  labs(title = "Non pan-essential hub dependencies",
       x = "Mean CRISPR gene effect — HPV-A7 lines",
       y = "Mean CRISPR gene effect — HPV-A9 lines") +
  base_theme +
  theme(legend.position = c(0.985, 0.02),
        legend.justification = c("right", "bottom"),
        legend.background = element_rect(fill = alpha("white", 0.92), color = "grey75"),
        legend.key = element_blank())

ggsave(file.path(res_dir, "6_9_B_essential_hub_scatter_nonpan.png"), pB,
       width = 10, height = 10, dpi = 300, bg = "white")
ggsave(file.path(res_dir, "6_9_B_essential_hub_scatter_nonpan.pdf"), pB,
       width = 10, height = 10, device = cairo_pdf, bg = "white")
message("✅ pan-essential-free scatter saved to: ", res_dir)
print(as.data.frame(grp_n_np))
