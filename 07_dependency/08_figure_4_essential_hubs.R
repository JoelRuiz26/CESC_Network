### ============================================================
### 6_11_Main_figure.R
###
### Composite main figure assembled with patchwork.
###   A  degree vs gene effect, one panel per network, split into four
###      quadrants by the data-driven thresholds (as in 6_10_A)
###   B  Venn of the essential hub sets of the two clades
###   C  gene effect of the essential hubs with pan-essentials removed
###      (as in 6_9_B)
###
### Thresholds are the elbow of each ranked curve, computed within each
### network (see 6_10 and its supplementary figure):
###   HPV-A7  degree >= 35   gene effect <= -0.317
###   HPV-A9  degree >= 42   gene effect <= -0.316
###
### The three panels share one palette: #4E79A7 for HPV-A7, #CB6A5C for
### HPV-A9 and #8C79B5 for genes essential in both.
###
###   6_11_Main_figure.png/.pdf
### ============================================================
suppressPackageStartupMessages({
  library(tidyverse); library(vroom); library(patchwork)
  library(ggrepel);   library(ggvenn)
})

res_dir <- "~/CESC_Network/07_dependency/results"
fig_dir <- "~/CESC_Network/figures"
FF   <- "Times New Roman"
A7c  <- "#4E79A7"
A9c  <- "#CB6A5C"
BOTH <- "#8C79B5"
cols_net <- c("HPV-A7" = A7c, "HPV-A9" = A9c)

## ── Thresholds: elbow of each ranked curve ───────────────────────────────────
elbow <- function(v) {
  n <- length(v)
  x <- (seq_len(n) - 1) / (n - 1)
  y <- (v - v[1]) / (v[n] - v[1])
  list(idx = which.max(abs(y - x)), val = v[which.max(abs(y - x))])
}

tab <- vroom(file.path(res_dir, "6_2_gene_dep_degree_table.tsv"), show_col_types = FALSE) %>%
  filter(!is.na(degree), !is.na(mean_dep)) %>%
  distinct(network, gene, .keep_all = TRUE)

nets <- sort(unique(tab$network))
info <- map(nets, function(nk) {
  d  <- tab %>% filter(network == nk)
  ek <- elbow(sort(d$degree, decreasing = TRUE))
  ed <- elbow(sort(d$mean_dep))
  N  <- nrow(d)
  q  <- d %>% mutate(hi_k = degree >= ek$val, hi_d = mean_dep <= ed$val) %>%
    count(hi_k, hi_d) %>% mutate(pct = 100 * n / N)
  getq <- function(k, dd) {
    r <- q %>% filter(hi_k == k, hi_d == dd)
    if (!nrow(r)) list(n = 0L, pct = 0) else list(n = r$n, pct = r$pct)
  }
  list(net = nk, d = d, N = N, ek = ek, ed = ed,
       qTL = getq(FALSE, FALSE), qTR = getq(TRUE, FALSE),
       qBL = getq(FALSE, TRUE),  qBR = getq(TRUE, TRUE),
       sel = d %>% filter(degree >= ek$val, mean_dep <= ed$val))
}) %>% set_names(nets)

ess_A7 <- toupper(info[["HPV-A7"]]$sel$gene)
ess_A9 <- toupper(info[["HPV-A9"]]$sel$gene)

cat("essential hubs  A7:", length(ess_A7), " A9:", length(ess_A9),
    " shared:", length(intersect(ess_A7, ess_A9)), "\n")

## ── Theme shared by the composite panels ─────────────────────────────────────
base <- theme_bw(base_family = FF) +
  theme(axis.text        = element_text(size = 13),
        axis.title       = element_text(size = 15, face = "bold"),
        plot.title       = element_text(size = 21, face = "bold"),
        legend.position  = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey95", linewidth = .3),
        panel.border     = element_rect(colour = "grey35", linewidth = .7),
        plot.margin      = margin(8, 12, 8, 8))

## ══════════════════════════════════════════════════════════════════════════════
## PANEL A: the four quadrants, one plot per network
## ══════════════════════════════════════════════════════════════════════════════
ylim_all <- range(tab$mean_dep, na.rm = TRUE) + c(-0.05, 0.05)

## gene labels, quadrant boxes and legend keys all share one size across panels
ANNO <- 5.0

qlab <- function(x, y, hj, vj, txt)
  annotate("label", x = x, y = y, label = txt, family = FF, fontface = "bold",
           size = ANNO, hjust = hj, vjust = vj, colour = "grey15",
           fill = alpha("white", .82), label.padding = unit(.24, "lines"))

pQuad <- map(info, function(i) {
  lab <- i$sel %>% slice_max(degree, n = 8)
  ggplot(i$d, aes(degree, mean_dep)) +
    annotate("rect", xmin = i$ek$val, xmax = Inf, ymin = -Inf, ymax = i$ed$val,
             fill = cols_net[i$net], alpha = .09) +
    geom_point(colour = "grey70", size = .85, alpha = .9, stroke = 0) +
    geom_point(data = i$sel, colour = cols_net[i$net], size = 1.7, alpha = .85, stroke = 0) +
    geom_vline(xintercept = i$ek$val, colour = "grey20", linetype = "22", linewidth = .55) +
    geom_hline(yintercept = i$ed$val, colour = "grey20", linetype = "22", linewidth = .55) +
    qlab(-Inf,  Inf, -0.09,  1.35, sprintf("%.1f%%", i$qTL$pct)) +
    qlab( Inf,  Inf,  1.09,  1.35, sprintf("%.1f%%", i$qTR$pct)) +
    qlab(-Inf, -Inf, -0.09, -0.45, sprintf("%.1f%%", i$qBL$pct)) +
    qlab( Inf, -Inf,  1.09, -0.45, sprintf("%.1f%%  (n = %d)", i$qBR$pct, i$qBR$n)) +
    geom_text_repel(data = lab, aes(label = gene), family = FF, fontface = "bold",
                    size = ANNO, colour = "grey12", max.overlaps = Inf, seed = 1,
                    box.padding = .5, force = 6, min.segment.length = 0,
                    segment.colour = "grey55", segment.size = .25) +
    scale_y_continuous(limits = ylim_all) +
    labs(title = i$net, x = "Degree", y = "Mean CRISPR gene effect") + base
})

panel_A <- wrap_elements(full = (
  (pQuad[[1]] | pQuad[[2]]) +
    plot_annotation(title = "Connectivity and dependency across networks",
                    theme = theme(plot.title = element_text(family = FF, size = 21,
                                                            face = "bold", hjust = 0,
                                                            margin = margin(l = 48, b = 4))))
))

## ══════════════════════════════════════════════════════════════════════════════
## PANEL B: Venn of the two essential hub sets
## ══════════════════════════════════════════════════════════════════════════════
## ggvenn only draws the two sets side by side, which wastes the height of this
## slot, so the Venn is drawn directly: two circles stacked vertically, with the
## region counts inside each. Set totals are not repeated here; they are already
## given in panel A as n = 395 and n = 436 inside each shaded quadrant.
circle_xy <- function(cx, cy, r, n = 300) {
  t <- seq(0, 2 * pi, length.out = n)
  tibble(x = cx + r * cos(t), y = cy + r * sin(t))
}

n_A7 <- length(setdiff(ess_A7, ess_A9))
n_sh <- length(intersect(ess_A7, ess_A9))
n_A9 <- length(setdiff(ess_A9, ess_A7))

panel_B <- ggplot() +
  geom_polygon(data = circle_xy(0,  0.60, 1), aes(x, y),
               fill = A7c, alpha = .5, colour = "grey20", linewidth = .55) +
  geom_polygon(data = circle_xy(0, -0.60, 1), aes(x, y),
               fill = A9c, alpha = .5, colour = "grey20", linewidth = .55) +
  annotate("text", x = 0, y =  1.05, label = n_A7, family = FF, size = 6.5) +
  annotate("text", x = 0, y =  0.00, label = n_sh, family = FF, size = 6.5) +
  annotate("text", x = 0, y = -1.05, label = n_A9, family = FF, size = 6.5) +
  annotate("text", x = 0, y =  1.95, label = "HPV-A7", family = FF,
           fontface = "bold", size = 5.8) +
  annotate("text", x = 0, y = -1.95, label = "HPV-A9", family = FF,
           fontface = "bold", size = 5.8) +
  coord_fixed(xlim = c(-1.45, 1.45), ylim = c(-2.75, 2.75), expand = FALSE) +
  labs(title = "Essential hubs per clade") +
  theme_void(base_family = FF) +
  theme(plot.title  = element_text(hjust = 0.5, face = "bold", size = 21, family = FF),
        plot.margin = margin(8, 8, 8, 8))

## save the three lists alongside the figure
shared <- sort(intersect(ess_A7, ess_A9))
onlyA7 <- sort(setdiff(ess_A7, ess_A9))
onlyA9 <- sort(setdiff(ess_A9, ess_A7))
nmax <- max(length(shared), length(onlyA7), length(onlyA9))
pad  <- function(x) c(x, rep(NA_character_, nmax - length(x)))
vroom_write(tibble(shared = pad(shared), A7_only = pad(onlyA7), A9_only = pad(onlyA9)),
            file.path(res_dir, "6_11_essential_hub_venn_lists.tsv"), delim = "\t")

## ══════════════════════════════════════════════════════════════════════════════
## PANEL C: gene effect scatter with pan-essentials removed
## ══════════════════════════════════════════════════════════════════════════════
lvl  <- c("Hub, not essential", "Essential hub in A7 only",
          "Essential hub in A9 only", "Essential hub in both")
cols <- c("Hub, not essential"       = "grey87",
          "Essential hub in A7 only" = A7c,
          "Essential hub in A9 only" = A9c,
          "Essential hub in both"    = BOTH)
short <- c("Essential hub in both"    = "Both",
           "Essential hub in A7 only" = "HPV-A7 only",
           "Essential hub in A9 only" = "HPV-A9 only")

flags <- tab %>%
  left_join(map_dfr(info, ~ tibble(network = .x$net, k_cut = .x$ek$val, d_cut = .x$ed$val)),
            by = "network") %>%
  mutate(hub = degree >= k_cut, ess = hub & mean_dep <= d_cut,
         net = str_remove(network, "HPV-")) %>%
  select(gene, net, hub, ess) %>%
  pivot_wider(names_from = net, values_from = c(hub, ess)) %>%
  mutate(across(starts_with("hub_") | starts_with("ess_"), ~ replace_na(.x, FALSE)))

pan <- vroom(file.path(res_dir, "6_1_common_essential_genes.tsv"),
             delim = "\t", show_col_types = FALSE)$gene_u %>% toupper()

dat_np <- vroom(file.path(res_dir, "6_7_hub_dependency_medians.tsv"), show_col_types = FALSE) %>%
  select(gene, gene_u, dep_A7, dep_A9) %>%
  inner_join(flags, by = "gene") %>%
  filter(hub_A7 | hub_A9, !(toupper(gene_u) %in% pan)) %>%
  mutate(grp = factor(case_when(
    ess_A7 & ess_A9 ~ "Essential hub in both",
    ess_A7          ~ "Essential hub in A7 only",
    ess_A9          ~ "Essential hub in A9 only",
    TRUE            ~ "Hub, not essential"), levels = lvl))

x_cut <- info[["HPV-A7"]]$ed$val
y_cut <- info[["HPV-A9"]]$ed$val

gn        <- dat_np %>% count(grp)
grp_labs  <- setNames(sprintf("%s  (n = %d)", short[as.character(gn$grp)], gn$n), gn$grp)
leg_order <- c("Essential hub in both", "Essential hub in A7 only", "Essential hub in A9 only")

## labels: within each quadrant, the genes furthest from the intersection of the
## two thresholds, by euclidean distance. Both axes are gene effect, so the
## distance is in the units of the plot and matches what the eye reads. A summed
## percentile rank was tried first and rejected: it favours genes that are mildly
## extreme on both axes over genes that are strongly extreme on one, so it left
## CCNC and NHLRC2 unlabelled while labelling points closer to the crossing.
## The upper-right quadrant holds 1829 grey genes and no highlighted one, so it
## contributes no labels.
TOP_N <- 5
hi <- dat_np %>%
  filter(grp != "Hub, not essential") %>%
  mutate(quad = paste0(if_else(dep_A7 <= x_cut, "L", "R"),
                       if_else(dep_A9 <= y_cut, "B", "T")),
         dist = sqrt((dep_A7 - x_cut)^2 + (dep_A9 - y_cut)^2)) %>%
  group_by(quad) %>%
  slice_max(dist, n = TOP_N, with_ties = FALSE) %>%
  ungroup()

## the axes are cropped to the range of the highlighted genes, not of the whole
## grey cloud: the cloud runs up to +0.2 in both axes and pushed the coloured
## points into a corner. 200 of the 1861 grey genes (10.7%) fall outside the
## crop; all 90 coloured ones are inside it.
hl  <- dat_np %>% filter(grp != "Hub, not essential")
rng <- range(c(hl$dep_A7, hl$dep_A9), na.rm = TRUE)
lim <- c(rng[1] - 0.05, rng[2] + 0.05)

grey_np <- dat_np %>% filter(grp == "Hub, not essential")
col_np  <- dat_np %>% filter(grp != "Hub, not essential") %>% droplevels()

panel_C <- ggplot(mapping = aes(dep_A7, dep_A9)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey87", linewidth = .4) +
  geom_hline(yintercept = y_cut, linetype = "dashed", linewidth = .4, colour = "grey62") +
  geom_vline(xintercept = x_cut, linetype = "dashed", linewidth = .4, colour = "grey62") +
  geom_point(data = grey_np, colour = "grey75", size = 2.2, alpha = .6, stroke = 0) +
  geom_point(data = col_np, aes(fill = grp), shape = 21, colour = "white",
             size = 3.4, stroke = .3, alpha = .92) +
  geom_text_repel(data = hi, aes(label = gene), colour = "grey15",
                  size = ANNO, fontface = "bold", family = FF,
                  max.overlaps = Inf, min.segment.length = 0, box.padding = .6,
                  point.padding = .3, seed = 1, force = 6, force_pull = .4,
                  segment.colour = "grey60", segment.size = .22,
                  bg.color = "white", bg.r = .12) +
  scale_fill_manual(values = cols, name = "Essential hub",
                    breaks = leg_order, labels = grp_labs[leg_order]) +
  guides(fill = guide_legend(override.aes = list(size = 5.2, alpha = 1))) +
  coord_fixed(xlim = lim, ylim = lim, expand = FALSE) +
  labs(title = "Selective hub dependencies in both networks",
       x = "Mean CRISPR gene effect — HPV-A7 lines",
       y = "Mean CRISPR gene effect — HPV-A9 lines") +
  base +
  theme(legend.position = c(0.98, 0.02),
        legend.justification = c("right", "bottom"),
        legend.text  = element_text(size = 15, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.background = element_rect(fill = alpha("white", .92), colour = "grey75"),
        legend.key = element_blank())

## ══════════════════════════════════════════════════════════════════════════════
## Composite
## ══════════════════════════════════════════════════════════════════════════════
fig <- panel_A / (panel_B + panel_C + plot_layout(widths = c(0.48, 1.52))) +
  plot_layout(heights = c(1, 1.95)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = FF, size = 24, face = "bold"))

ggsave(file.path(fig_dir, "Figure_4_essential_hubs.png"), fig,
       width = 16, height = 16, dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "Figure_4_essential_hubs.pdf"), fig,
       width = 16, height = 16, device = cairo_pdf, bg = "white")

message("main figure saved to: ", res_dir)
print(as.data.frame(gn))
