### ============================================================
### 6_10_Elbow_quadrants.R
###
### Thresholds for degree and dependency are taken from the data: the
### genes of each network are ranked by each measure and the threshold
### is the elbow of the ranked curve, defined as the point of maximum
### distance to the chord joining the two ends of the curve, with both
### axes normalised to [0, 1]. The same rule is applied to both measures
### and to both networks.
###
### Outputs
###   6_10_A_Degree_dependency_quadrants.png/.pdf   main figure
###       the two thresholds split each network into four quadrants,
###       with the percentage of genes falling in each
###   6_10_S_Elbow_threshold_curves.png/.pdf        supplementary
###       the four ranked curves and the elbow on each, i.e. where the
###       thresholds of the main figure come from
###
### Reads 6_2_gene_dep_degree_table.tsv (degree and clade mean
### dependency per gene per network). No DepMap download.
### ============================================================
suppressPackageStartupMessages({
  library(tidyverse); library(vroom); library(patchwork); library(ggrepel)
})
res_dir <- "~/CESC_Network/6_DepMap_cacu/results"
FF   <- "Times New Roman"
cols <- c("HPV-A7" = "#4E79A7", "HPV-A9" = "#CB6A5C")
labs_net <- c("HPV-A7" = "A7 network", "HPV-A9" = "A9 network")

tab <- vroom(file.path(res_dir, "6_2_gene_dep_degree_table.tsv"), show_col_types = FALSE) %>%
  filter(!is.na(degree), !is.na(mean_dep)) %>%
  distinct(network, gene, .keep_all = TRUE)

## ── elbow of a ranked curve: maximum distance to the chord ────────────────────
elbow <- function(v) {
  n <- length(v)
  x <- (seq_len(n) - 1) / (n - 1)
  y <- (v - v[1]) / (v[n] - v[1])
  i <- which.max(abs(y - x))
  list(idx = i, val = v[i])
}

nets <- sort(unique(tab$network))
info <- map(nets, function(nk) {
  d  <- tab %>% filter(network == nk)
  kv <- sort(d$degree, decreasing = TRUE)
  dv <- sort(d$mean_dep)
  ek <- elbow(kv); ed <- elbow(dv)
  N  <- nrow(d)
  q  <- d %>%
    mutate(hi_k = degree   >= ek$val,
           hi_d = mean_dep <= ed$val) %>%
    count(hi_k, hi_d) %>%
    mutate(pct = 100 * n / N)
  getq <- function(k, dd) {
    r <- q %>% filter(hi_k == k, hi_d == dd)
    if (!nrow(r)) list(n = 0L, pct = 0) else list(n = r$n, pct = r$pct)
  }
  list(net = nk, d = d, N = N, kv = kv, dv = dv, ek = ek, ed = ed,
       qTL = getq(FALSE, FALSE), qTR = getq(TRUE, FALSE),
       qBL = getq(FALSE, TRUE),  qBR = getq(TRUE, TRUE),
       sel = d %>% filter(degree >= ek$val, mean_dep <= ed$val))
}) %>% set_names(nets)

cat("\n================ ELBOW THRESHOLDS ================\n")
walk(info, function(i) {
  cat("\n", i$net, "  (N = ", i$N, ")\n", sep = "")
  cat("  degree     : >= ", i$ek$val, "   rank ", i$ek$idx,
      " (top ", round(100 * i$ek$idx / i$N, 1), "%)\n", sep = "")
  cat("  dependency : <= ", round(i$ed$val, 3), "  rank ", i$ed$idx,
      " (", round(100 * i$ed$idx / i$N, 1), "%)\n", sep = "")
  cat("  quadrants  : high degree + strong dependency ", i$qBR$n,
      " (", round(i$qBR$pct, 1), "%)\n", sep = "")
  cat("               low degree  + strong dependency ", i$qBL$n,
      " (", round(i$qBL$pct, 1), "%)\n", sep = "")
  cat("               high degree + weak dependency   ", i$qTR$n,
      " (", round(i$qTR$pct, 1), "%)\n", sep = "")
  cat("               low degree  + weak dependency   ", i$qTL$n,
      " (", round(i$qTL$pct, 1), "%)\n", sep = "")
})

## ── Shared theme, sized for a full-page figure ────────────────────────────────
base <- theme_bw(base_family = FF) +
  theme(axis.text        = element_text(size = 17),
        axis.title       = element_text(size = 20, face = "bold"),
        plot.title       = element_text(size = 22, face = "bold"),
        plot.subtitle    = element_text(size = 15, face = "italic", colour = "grey35"),
        legend.position  = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey95", linewidth = .3),
        panel.border     = element_rect(colour = "grey35", linewidth = .7),
        plot.margin      = margin(10, 14, 10, 10))

## ══════════════════════════════════════════════════════════════════════════════
## MAIN FIGURE: the four quadrants
## ══════════════════════════════════════════════════════════════════════════════
## the dependency scale is the same measurement in both networks, so the y axis
## is shared; degree is not, because it scales with each network's density
ylim_all <- range(tab$mean_dep, na.rm = TRUE) + c(-0.05, 0.05)

qlab <- function(x, y, hj, vj, txt)
  annotate("label", x = x, y = y, label = txt, family = FF, fontface = "bold",
           size = 6.2, hjust = hj, vjust = vj, colour = "grey15",
           fill = alpha("white", .82), label.size = 0, label.padding = unit(.28, "lines"))

pMain <- map(info, function(i) {
  lab <- i$sel %>% slice_max(degree, n = 10)
  ggplot(i$d, aes(degree, mean_dep)) +
    annotate("rect", xmin = i$ek$val, xmax = Inf, ymin = -Inf, ymax = i$ed$val,
             fill = cols[i$net], alpha = .09) +
    geom_point(colour = "grey78", size = .9, alpha = .72, stroke = 0) +
    geom_point(data = i$sel, colour = cols[i$net], size = 2, alpha = .85, stroke = 0) +
    geom_vline(xintercept = i$ek$val, colour = "grey20", linetype = "22", linewidth = .6) +
    geom_hline(yintercept = i$ed$val, colour = "grey20", linetype = "22", linewidth = .6) +
    qlab(-Inf,  Inf, -0.08,  1.3, sprintf("%.1f%%", i$qTL$pct)) +
    qlab( Inf,  Inf,  1.08,  1.3, sprintf("%.1f%%", i$qTR$pct)) +
    qlab(-Inf, -Inf, -0.08, -0.4, sprintf("%.1f%%", i$qBL$pct)) +
    qlab( Inf, -Inf,  1.08, -0.4, sprintf("%.1f%%  (n = %d)", i$qBR$pct, i$qBR$n)) +
    geom_text_repel(data = lab, aes(label = gene), family = FF, fontface = "bold",
                    size = 4.6, colour = "grey12", max.overlaps = Inf, seed = 1,
                    box.padding = .6, force = 6, min.segment.length = 0,
                    segment.colour = "grey55", segment.size = .3) +
    scale_y_continuous(limits = ylim_all) +
    labs(title = i$net,
         x = "Degree",
         y = "Mean CRISPR gene effect") + base
})

fig_main <- (pMain[[1]] | pMain[[2]]) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = FF, size = 26, face = "bold"))

ggsave(file.path(res_dir, "6_10_A_Degree_dependency_quadrants.png"), fig_main,
       width = 17, height = 9, dpi = 300, bg = "white")
ggsave(file.path(res_dir, "6_10_A_Degree_dependency_quadrants.pdf"), fig_main,
       width = 17, height = 9, device = cairo_pdf, bg = "white")

## ══════════════════════════════════════════════════════════════════════════════
## SUPPLEMENTARY FIGURE: where each threshold comes from
## ══════════════════════════════════════════════════════════════════════════════
## the label used to sit on the elbow point, where it crossed the curve; it now
## goes to the corner of the panel that the curve leaves empty: above the curve
## for the decaying degree ranking, below it for the rising dependency ranking
rank_panel <- function(i, v, e, ttl, sub, ylab, lab, corner = "top") {
  n <- length(v)
  chord_y <- v[1] + (v[n] - v[1]) * (e$idx - 1) / (n - 1)
  ggplot(tibble(rank = seq_len(n), y = v), aes(rank, y)) +
    geom_segment(x = 1, y = v[1], xend = n, yend = v[n],
                 colour = "grey62", linetype = "22", linewidth = .55) +
    geom_line(colour = cols[i$net], linewidth = 1.2) +
    geom_segment(x = e$idx, y = v[e$idx], xend = e$idx, yend = chord_y,
                 colour = "grey30", linewidth = .6) +
    geom_point(x = e$idx, y = e$val, colour = "grey10", size = 4) +
    annotate("text", x = Inf, y = if (corner == "top") Inf else -Inf,
             family = FF, fontface = "bold", size = 7.5, label = lab,
             hjust = 1.08, vjust = if (corner == "top") 1.7 else -1.1) +
    labs(title = ttl, subtitle = sub, x = "Gene rank", y = ylab) + base +
    ## larger fonts for the standalone supplementary figure only
    theme(plot.title    = element_text(size = 28, face = "bold"),
          plot.subtitle = element_text(size = 20, face = "italic", colour = "grey35"),
          axis.title    = element_text(size = 26, face = "bold"),
          axis.text     = element_text(size = 22))
}

## Hub size is the tie-safe count of genes AT OR ABOVE the elbow degree, not the
## elbow rank index: many genes tie at the threshold degree (171 at degree 35 in
## A7, 120 at degree 42 in A9), so the index undercounts the set. This value cut
## is the hub set used by the quadrant figure and by 6_3.
nhub <- function(i) sum(i$d$degree >= i$ek$val)

pS <- list(
  rank_panel(info[[1]], info[[1]]$kv, info[[1]]$ek,
             paste0(labs_net[nets[1]], " — ranked by connectivity"),
             paste0("threshold at degree ≥ ", info[[1]]$ek$val, ", the top ",
                    round(100 * nhub(info[[1]]) / info[[1]]$N, 1), "% (",
                    nhub(info[[1]]), " of ", info[[1]]$N, " genes)"),
             "Degree", paste0("degree ≥ ", info[[1]]$ek$val)),
  rank_panel(info[[1]], info[[1]]$dv, info[[1]]$ed,
             paste0(labs_net[nets[1]], " — ranked by dependency"),
             paste0("threshold at rank ", info[[1]]$ed$idx, " of ", info[[1]]$N,
                    ", the strongest ", round(100 * info[[1]]$ed$idx / info[[1]]$N, 1), "%"),
             "Mean CRISPR gene effect", paste0("gene effect ≤ ", sprintf("%.3f", info[[1]]$ed$val)),
             corner = "bottom"),
  rank_panel(info[[2]], info[[2]]$kv, info[[2]]$ek,
             paste0(labs_net[nets[2]], " — ranked by connectivity"),
             paste0("threshold at degree ≥ ", info[[2]]$ek$val, ", the top ",
                    round(100 * nhub(info[[2]]) / info[[2]]$N, 1), "% (",
                    nhub(info[[2]]), " of ", info[[2]]$N, " genes)"),
             "Degree", paste0("degree ≥ ", info[[2]]$ek$val)),
  rank_panel(info[[2]], info[[2]]$dv, info[[2]]$ed,
             paste0(labs_net[nets[2]], " — ranked by dependency"),
             paste0("threshold at rank ", info[[2]]$ed$idx, " of ", info[[2]]$N,
                    ", the strongest ", round(100 * info[[2]]$ed$idx / info[[2]]$N, 1), "%"),
             "Mean CRISPR gene effect", paste0("gene effect ≤ ", sprintf("%.3f", info[[2]]$ed$val)),
             corner = "bottom")
)

fig_supp <- (pS[[1]] | pS[[2]]) / (pS[[3]] | pS[[4]])

ggsave(file.path(res_dir, "6_10_S_Elbow_threshold_curves.png"), fig_supp,
       width = 17, height = 13, dpi = 300, bg = "white")
ggsave(file.path(res_dir, "6_10_S_Elbow_threshold_curves.pdf"), fig_supp,
       width = 17, height = 13, device = cairo_pdf, bg = "white")

walk(info, function(i)
  vroom_write(i$sel %>% select(gene, degree, mean_dep) %>% arrange(mean_dep),
              file.path(res_dir, paste0("6_10_selected_genes_",
                                        str_remove(i$net, "HPV-"), ".tsv"))))
message("✅ 6_10 main and supplementary figures saved to: ", res_dir)
