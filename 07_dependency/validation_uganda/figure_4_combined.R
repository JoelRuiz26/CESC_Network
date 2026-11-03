### ============================================================
### figure_4_combined.R
###
### CESC + HTMCP-CC combined view of Figure 4. Does NOT touch the original
### CESC-only figure. TWO separate output figures:
###
### Figure 4 (essential hubs): grouped by COHORT -- each cohort's two
### quadrant scatters (A7, A9) are stacked vertically, with the essential-hub
### Venn they share to the right; the two cohort blocks sit side by side.
###   [ (HPV-A7 quad / HPV-A9 quad) | shared Venn ]   CESC | HTMCP-CC
### -> figures/validation_uganda/Figure_4_essential_hubs_uganda.png/.pdf
###
### Figure 4b (selective hubs): the selective-hub dependency scatter on its
### own, CESC | HTMCP-CC side by side.
### -> figures/validation_uganda/Figure_4_selective_hub_dependency_uganda.png/.pdf
### ============================================================
suppressPackageStartupMessages({
  library(tidyverse); library(vroom); library(patchwork)
  library(ggrepel);   library(ggvenn)
})

fig_dir <- "~/CESC_Network/figures/validation_uganda"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dep_shared_dir <- "~/CESC_Network/07_dependency/results"
FF   <- "Times New Roman"

elbow <- function(v) {
  n <- length(v); x <- (seq_len(n) - 1) / (n - 1); y <- (v - v[1]) / (v[n] - v[1])
  list(idx = which.max(abs(y - x)), val = v[which.max(abs(y - x))])
}
base <- theme_bw(base_family = FF) +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 19, face = "bold"), legend.position = "none",
        panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey95", linewidth = .3),
        panel.border = element_rect(colour = "grey35", linewidth = .7), plot.margin = margin(8, 12, 8, 8))
ANNO <- 5.2

build_cohort <- function(res_dir, cohort_label, A7c, A9c, BOTH) {
  cols_net <- c("HPV-A7" = A7c, "HPV-A9" = A9c)
  tab <- vroom(file.path(res_dir, "6_2_gene_dep_degree_table.tsv"), show_col_types = FALSE) %>%
    filter(!is.na(degree), !is.na(mean_dep)) %>% distinct(network, gene, .keep_all = TRUE)
  nets <- sort(unique(tab$network))
  info <- map(nets, function(nk) {
    d <- tab %>% filter(network == nk)
    ek <- elbow(sort(d$degree, decreasing = TRUE)); ed <- elbow(sort(d$mean_dep)); N <- nrow(d)
    q <- d %>% mutate(hi_k = degree >= ek$val, hi_d = mean_dep <= ed$val) %>% count(hi_k, hi_d) %>% mutate(pct = 100 * n / N)
    getq <- function(k, dd) { r <- q %>% filter(hi_k == k, hi_d == dd); if (!nrow(r)) list(n = 0L, pct = 0) else list(n = r$n, pct = r$pct) }
    list(net = nk, d = d, N = N, ek = ek, ed = ed, qTL = getq(FALSE, FALSE), qTR = getq(TRUE, FALSE),
         qBL = getq(FALSE, TRUE), qBR = getq(TRUE, TRUE), sel = d %>% filter(degree >= ek$val, mean_dep <= ed$val))
  }) %>% set_names(nets)

  ess_A7 <- toupper(info[["HPV-A7"]]$sel$gene); ess_A9 <- toupper(info[["HPV-A9"]]$sel$gene)
  cat(sprintf("[%s] essential hubs A7:%d A9:%d shared:%d\n", cohort_label, length(ess_A7), length(ess_A9), length(intersect(ess_A7, ess_A9))))

  ylim_all <- range(tab$mean_dep, na.rm = TRUE) + c(-0.05, 0.05)
  qlab <- function(x, y, hj, vj, txt) annotate("label", x = x, y = y, label = txt, family = FF, fontface = "bold",
    size = ANNO, hjust = hj, vjust = vj, colour = "grey15", fill = alpha("white", .82), label.padding = unit(.22, "lines"))
  pQuad <- map(info, function(i) {
    lab <- i$sel %>% slice_max(degree, n = 8)
    ggplot(i$d, aes(degree, mean_dep)) +
      annotate("rect", xmin = i$ek$val, xmax = Inf, ymin = -Inf, ymax = i$ed$val, fill = cols_net[i$net], alpha = .09) +
      geom_point(colour = "grey70", size = .8, alpha = .9, stroke = 0) +
      geom_point(data = i$sel, colour = cols_net[i$net], size = 1.5, alpha = .85, stroke = 0) +
      geom_vline(xintercept = i$ek$val, colour = "grey20", linetype = "22", linewidth = .5) +
      geom_hline(yintercept = i$ed$val, colour = "grey20", linetype = "22", linewidth = .5) +
      qlab(-Inf, Inf, -0.09, 1.35, sprintf("%.1f%%", i$qTL$pct)) + qlab(Inf, Inf, 1.09, 1.35, sprintf("%.1f%%", i$qTR$pct)) +
      qlab(-Inf, -Inf, -0.09, -0.45, sprintf("%.1f%%", i$qBL$pct)) + qlab(Inf, -Inf, 1.09, -0.45, sprintf("%.1f%%  (n=%d)", i$qBR$pct, i$qBR$n)) +
      geom_text_repel(data = lab, aes(label = gene), family = FF, fontface = "bold", size = ANNO, colour = "grey12",
                      max.overlaps = Inf, seed = 1, box.padding = .5, force = 6, min.segment.length = 0, segment.colour = "grey55", segment.size = .22) +
      scale_y_continuous(limits = ylim_all) +
      labs(title = paste0(i$net, " — ", cohort_label), x = "Degree", y = "Mean CRISPR gene effect") + base
  })

  circle_xy <- function(cx, cy, r, n = 300) { t <- seq(0, 2 * pi, length.out = n); tibble(x = cx + r * cos(t), y = cy + r * sin(t)) }
  n_A7 <- length(setdiff(ess_A7, ess_A9)); n_sh <- length(intersect(ess_A7, ess_A9)); n_A9 <- length(setdiff(ess_A9, ess_A7))
  n_tot <- n_A7 + n_sh + n_A9
  lab_A7 <- sprintf("%d\n(%.1f%%)", n_A7, 100 * n_A7 / n_tot)
  lab_sh <- sprintf("%d\n(%.1f%%)", n_sh, 100 * n_sh / n_tot)
  lab_A9 <- sprintf("%d\n(%.1f%%)", n_A9, 100 * n_A9 / n_tot)
  panel_B <- ggplot() +
    geom_polygon(data = circle_xy(0, 0.60, 1), aes(x, y), fill = A7c, alpha = .5, colour = "grey20", linewidth = .55) +
    geom_polygon(data = circle_xy(0, -0.60, 1), aes(x, y), fill = A9c, alpha = .5, colour = "grey20", linewidth = .55) +
    annotate("text", x = 0, y = 1.05, label = lab_A7, family = FF, size = 5.6, lineheight = 0.9) +
    annotate("text", x = 0, y = 0.00, label = lab_sh, family = FF, size = 5.6, lineheight = 0.9) +
    annotate("text", x = 0, y = -1.05, label = lab_A9, family = FF, size = 5.6, lineheight = 0.9) +
    annotate("text", x = 0, y = 1.95, label = "HPV-A7", family = FF, fontface = "bold", size = 5.8) +
    annotate("text", x = 0, y = -1.95, label = "HPV-A9", family = FF, fontface = "bold", size = 5.8) +
    coord_fixed(xlim = c(-1.45, 1.45), ylim = c(-2.75, 2.75), expand = FALSE) +
    labs(title = cohort_label) + theme_void(base_family = FF) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 17, family = FF), plot.margin = margin(8, 8, 8, 8))

  lvl <- c("Hub, not essential", "Essential hub in A7 only", "Essential hub in A9 only", "Essential hub in both")
  cols <- c("Hub, not essential" = "grey87", "Essential hub in A7 only" = A7c, "Essential hub in A9 only" = A9c, "Essential hub in both" = BOTH)
  short <- c("Essential hub in both" = "Both", "Essential hub in A7 only" = "HPV-A7 only", "Essential hub in A9 only" = "HPV-A9 only")
  flags <- tab %>% left_join(map_dfr(info, ~ tibble(network = .x$net, k_cut = .x$ek$val, d_cut = .x$ed$val)), by = "network") %>%
    mutate(hub = degree >= k_cut, ess = hub & mean_dep <= d_cut, net = str_remove(network, "HPV-")) %>%
    select(gene, net, hub, ess) %>% pivot_wider(names_from = net, values_from = c(hub, ess)) %>%
    mutate(across(starts_with("hub_") | starts_with("ess_"), ~ replace_na(.x, FALSE)))
  pan <- vroom(file.path(dep_shared_dir, "6_1_common_essential_genes.tsv"), delim = "\t", show_col_types = FALSE)$gene_u %>% toupper()
  dat_np <- vroom(file.path(res_dir, "6_7_hub_dependency_medians.tsv"), show_col_types = FALSE) %>%
    select(gene, gene_u, dep_A7, dep_A9) %>% inner_join(flags, by = "gene") %>%
    filter(hub_A7 | hub_A9, !(toupper(gene_u) %in% pan)) %>%
    mutate(grp = factor(case_when(ess_A7 & ess_A9 ~ "Essential hub in both", ess_A7 ~ "Essential hub in A7 only",
                                  ess_A9 ~ "Essential hub in A9 only", TRUE ~ "Hub, not essential"), levels = lvl))
  x_cut <- info[["HPV-A7"]]$ed$val; y_cut <- info[["HPV-A9"]]$ed$val
  gn <- dat_np %>% count(grp); grp_labs <- setNames(sprintf("%s (n=%d)", short[as.character(gn$grp)], gn$n), gn$grp)
  leg_order <- c("Essential hub in both", "Essential hub in A7 only", "Essential hub in A9 only")
  hi <- dat_np %>% filter(grp != "Hub, not essential") %>%
    mutate(quad = paste0(if_else(dep_A7 <= x_cut, "L", "R"), if_else(dep_A9 <= y_cut, "B", "T")),
           dist = sqrt((dep_A7 - x_cut)^2 + (dep_A9 - y_cut)^2)) %>%
    group_by(quad) %>% slice_max(dist, n = 5, with_ties = FALSE) %>% ungroup()
  hl <- dat_np %>% filter(grp != "Hub, not essential"); rng <- range(c(hl$dep_A7, hl$dep_A9), na.rm = TRUE); lim <- c(rng[1] - 0.05, rng[2] + 0.05)
  grey_np <- dat_np %>% filter(grp == "Hub, not essential"); col_np <- dat_np %>% filter(grp != "Hub, not essential") %>% droplevels()

  panel_C <- ggplot(mapping = aes(dep_A7, dep_A9)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey87", linewidth = .4) +
    geom_hline(yintercept = y_cut, linetype = "dashed", linewidth = .4, colour = "grey62") +
    geom_vline(xintercept = x_cut, linetype = "dashed", linewidth = .4, colour = "grey62") +
    geom_point(data = grey_np, colour = "grey75", size = 1.8, alpha = .6, stroke = 0) +
    geom_point(data = col_np, aes(fill = grp), shape = 21, colour = "white", size = 2.8, stroke = .3, alpha = .92) +
    geom_text_repel(data = hi, aes(label = gene), colour = "grey15", size = ANNO, fontface = "bold", family = FF,
                    max.overlaps = Inf, min.segment.length = 0, box.padding = .6, point.padding = .3, seed = 1, force = 6,
                    force_pull = .4, segment.colour = "grey60", segment.size = .2, bg.color = "white", bg.r = .12) +
    scale_fill_manual(values = cols, name = "Essential hub", breaks = leg_order, labels = grp_labs[leg_order]) +
    guides(fill = guide_legend(override.aes = list(size = 4.2, alpha = 1))) +
    coord_fixed(xlim = lim, ylim = lim, expand = FALSE) +
    labs(title = cohort_label, x = "Mean CRISPR gene effect — HPV-A7 lines", y = "Mean CRISPR gene effect — HPV-A9 lines") +
    base + theme(legend.position = c(0.98, 0.02), legend.justification = c("right", "bottom"),
                legend.text = element_text(size = 13, face = "bold"), legend.title = element_text(size = 14, face = "bold"),
                legend.background = element_rect(fill = alpha("white", .92), colour = "grey75"), legend.key = element_blank())

  print(as.data.frame(gn))
  list(pA7 = pQuad[["HPV-A7"]], pA9 = pQuad[["HPV-A9"]], panel_B = panel_B, panel_C = panel_C)
}

cesc   <- build_cohort("~/CESC_Network/07_dependency/results", "CESC", "#4E79A7", "#CB6A5C", "#8C79B5")
uganda <- build_cohort("~/CESC_Network/07_dependency/validation_uganda/results", "HTMCP-CC", "#D55E00", "#009E73", "#6A7E39")

## ── Figure 4: per-cohort ROW -- its two quadrant scatters side by side
## (shorter/less-squished than stacking them), its shared essential-hub Venn
## at the end; the two cohort rows stacked vertically (CESC on top, HTMCP-CC
## below). Tags are applied INSIDE each row before wrap_elements(full=...)
## makes it opaque, with each row given its OWN three letters (A-C, then
## D-F) so the full six-panel figure reads as a single unique A-F sequence,
## as a published figure requires, instead of repeating A-C twice.
cohort_row <- function(pA7, pA9, venn, tags) {
  row <- (pA7 | pA9 | venn) +
    plot_layout(widths = c(1, 1, 0.55)) +
    plot_annotation(tag_levels = list(tags)) &
    theme(plot.tag = element_text(family = FF, size = 22, face = "bold"))
  wrap_elements(full = row)
}
rowCESC  <- cohort_row(cesc$pA7, cesc$pA9, cesc$panel_B, c("A", "B", "C"))
rowHTMCP <- cohort_row(uganda$pA7, uganda$pA9, uganda$panel_B, c("D", "E", "F"))

fig4 <- wrap_elements(full = (rowCESC / rowHTMCP))

ggsave("~/CESC_Network/figures/validation_uganda/Figure_4_essential_hubs_uganda.png", fig4, width = 20, height = 11, dpi = 300, bg = "white", limitsize = FALSE)
ggsave("~/CESC_Network/figures/validation_uganda/Figure_4_essential_hubs_uganda.pdf", fig4, width = 20, height = 11, device = cairo_pdf, bg = "white", limitsize = FALSE)
message("Saved: Figure_4_essential_hubs_uganda.png / .pdf  ([A7/A9 quad | shared Venn] per cohort; CESC | HTMCP-CC)")

## ── separate figure: selective-hub dependency scatter (the "bottom" scatters)
fig4b <- wrap_elements(full = (cesc$panel_C | uganda$panel_C)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = FF, size = 24, face = "bold"))

ggsave("~/CESC_Network/figures/validation_uganda/Figure_4_selective_hub_dependency_uganda.png", fig4b, width = 18, height = 9, dpi = 300, bg = "white", limitsize = FALSE)
ggsave("~/CESC_Network/figures/validation_uganda/Figure_4_selective_hub_dependency_uganda.pdf", fig4b, width = 18, height = 9, device = cairo_pdf, bg = "white", limitsize = FALSE)
message("Saved: Figure_4_selective_hub_dependency_uganda.png / .pdf  (selective hub dependency scatter; CESC | HTMCP-CC)")
