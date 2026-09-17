### ============================================================
### 6_3_Hub_dependency_CDF.R
###
### QUESTION: are hub genes shifted toward CRISPR essentiality, and
### does ANY random gene set of the SAME size do the same?
###
### Three figures, three jobs:
###   6_3_A_hub_vs_nonhub_cdf   hub vs non-hub CDF of dependency, with the
###                             two-sided KS (D, asymptotic p). Shows the
###                             effect in interpretable dependency units.
###   6_3_B_KS_permutation_null permutation null of D: does any random
###                             same-size set separate from its rest as much
###                             as the hub set? Shows the separation is not
###                             chance. D is non-directional, hence C.
###   6_3_C_mean_dependency_permutation_null  permutation null of the MEAN
###                             dependency: is the hub set shifted toward
###                             essentiality more than any random same-size
###                             set? The x axis is in dependency units, so
###                             the direction is read off the figure while
###                             the test stays two-sided.
###
### A and B run the SAME test: the TWO-SIDED two-sample KS,
###   D = max |F_set - F_non-set|.
### A one-sided D+ was considered and rejected. One-sided testing requires both
### that the direction be pre-specified AND that the opposite result be of no
### interest (Bland & Altman, BMJ 1994;309:248). The first holds here
### (centrality-lethality rule, Jeong et al. 2001 Nature); the second does NOT:
### hubs turning out LESS essential would refute that rule and would certainly
### be reported. One-sided would only halve an already overwhelming p
### (7e-10 -> 4e-10; 3e-42 -> 1e-42) at the cost of the standard objection.
### Two-sided is also the conservative choice for figure B: its null is wider
### (mean D 0.027 vs 0.012 one-sided) yet the verdict is unchanged (0/10,000).
###
### DIRECTION is therefore not carried by the test but MEASURED and reported
### separately, as KS_Dplus / KS_Dminus in the stats file:
###   D+ = max(F_hub - F_non-hub)      lean toward essential
###   D- = max(F_non-hub - F_hub)      lean toward non-essential
### D- is 0.001 in A7 and exactly 0 in A9 (complete stochastic dominance): the
### hub CDF never runs below the non-hub CDF. Dependency is negative-is-
### essential, so F_hub above F_non-hub means hubs pile up on the essential side.
### Figure C carries the direction as a FIGURE rather than a footnote, by
### permuting a signed statistic (the set's mean dependency) on a dependency axis.
###
### A gives the asymptotic p; B and C give permutation p against a size-matched
### null: B random sets of the hub's size drawn from the NON-HUB genes, each
### compared with its own complement.
###
### No z-score / sigma is reported anywhere. The D null is right-skewed by
### construction (a maximum of gaps, bounded at 0), so a sigma would invite an
### invalid sigma -> p reading; sigma is dropped from C too, for consistency and
### because the shift in dependency units is the interpretable effect size.
### The permutation P is the honest statement; it cannot resolve below 1/(B+1).
###
### D is evaluated at the observed values (exact; matches ks.test to ~1e-17),
### never on an interpolation grid, which would underestimate a maximum.
### ============================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(vroom)
  library(patchwork)
})

dep_dir <- "~/CESC_Network/07_dependency"
res_dir <- file.path(dep_dir, "results")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

B            <- 10000
set.seed(123)
net_colors <- c("HPV-A7" = "#003366", "HPV-A9" = "#CC2200")

gene_tbl <- vroom(file.path(res_dir, "6_2_gene_dep_degree_table.tsv"),
                  show_col_types = FALSE) %>%
  filter(!is.na(degree), !is.na(mean_dep)) %>%
  distinct(network, gene, .keep_all = TRUE)

## ── Data-driven hub threshold ─────────────────────────────────────────────────
## Hubs are no longer the top decile of degree. The cut-off is the elbow of the
## degree ranking, the SAME rule and the same code used to define the hubs in
## 6_10: rank the genes by decreasing degree, normalise both axes to [0, 1] and
## take the point of maximum distance to the chord joining the two ends of the
## curve. Computed per network, so each clade keeps its own connectivity cut-off.
elbow_val <- function(v) {                 # v sorted so the ranked curve bends once
  n <- length(v)
  x <- (seq_len(n) - 1) / (n - 1)
  y <- (v - v[1]) / (v[n] - v[1])
  v[which.max(abs(y - x))]
}
hub_thr <- vapply(names(net_colors), function(net) {
  d <- filter(gene_tbl, network == net)
  elbow_val(sort(d$degree, decreasing = TRUE))
}, numeric(1))
cat("\n=== Data-driven hub thresholds (degree elbow) ===\n")
print(hub_thr)

## ── Permutation ───────────────────────────────────────────────────────────────
dist_perm <- function(net) {
  x   <- filter(gene_tbl, network == net)
  dep <- x$mean_dep; N <- length(dep)
  thr <- hub_thr[[net]]; hub <- x$degree >= thr
  hub_idx <- which(hub); nonhub_idx <- which(!hub); nh <- length(hub_idx)

  ## In sorted order, set membership -> cumsum gives F_set at every observed
  ## value in O(N); the complement follows from
  ##   F_rest = (N*F_all - n*F_set)/(N - n).
  ## Under ties the ECDF is only defined at the last index of each tied run.
  ord <- order(dep); dep_s <- dep[ord]
  pos <- integer(N); pos[ord] <- seq_len(N)     # gene -> its position in dep_s
  Fall_s  <- seq_len(N) / N
  eval_at <- !duplicated(dep_s, fromLast = TRUE)
  D_exact <- function(idx) {
    n <- length(idx)
    m <- integer(N); m[pos[idx]] <- 1L
    Fset  <- cumsum(m) / n
    Frest <- (N * Fall_s - n * Fset) / (N - n)
    g <- (Fset - Frest)[eval_at]
    # two = the tested statistic; plus/minus only document the direction
    c(two = max(abs(g)), plus = max(g), minus = max(-g))
  }
  obs <- D_exact(hub_idx)
  Dnull <- numeric(B)
  ## Mean dependency of the set: the SIGNED companion to D. D says the hub set
  ## separates from its rest; it cannot say in which direction (it is a maximum
  ## of ABSOLUTE gaps). The mean is signed and lives in dependency units, so its
  ## permutation null answers "is the hub set shifted toward essentiality more
  ## than any random same-size set" while staying two-sided: we ask whether the
  ## hub is extreme in EITHER tail of the null, and report which one it lands in.
  ## Mean rather than median: the effect is tail-driven (an excess of strongly
  ## essential genes), and the median is blind to the tail (shift -0.026 vs
  ## -0.121 in A7). The median is still reported in the stats file as a
  ## robustness check, since a mean can be driven by outliers.
  mean_hub <- mean(dep[hub_idx]); med_hub <- median(dep[hub_idx])
  mean_null <- numeric(B); med_null <- numeric(B)
  for (b in seq_len(B)) {
    s <- sample(nonhub_idx, nh)                            # random set from NON-hubs
    Dnull[b]     <- D_exact(s)[["two"]]
    mean_null[b] <- mean(dep[s])
    med_null[b]  <- median(dep[s])
  }
  p_two <- function(obs, null) {                           # two-sided permutation p
    mu <- mean(null); (1 + sum(abs(null - mu) >= abs(obs - mu))) / (B + 1)
  }

  ## ECDF curves on a grid: for DRAWING figure A only, never for computing D.
  grid    <- seq(min(dep), max(dep), length.out = 400)
  Fall    <- findInterval(grid, dep_s) / N
  Fhub    <- findInterval(grid, sort(dep[hub_idx])) / nh
  Fnonhub <- (N * Fall - nh * Fhub) / (N - nh)

  list(net = net, nh = nh, grid = grid, Fhub = Fhub, Fnonhub = Fnonhub,
       D_hub = obs[["two"]], Dplus_hub = obs[["plus"]], Dminus_hub = obs[["minus"]],
       Dnull = Dnull,
       p_perm = (1 + sum(Dnull >= obs[["two"]])) / (B + 1),
       n_ge   = sum(Dnull >= obs[["two"]]),
       mean_hub = mean_hub, mean_null = mean_null,
       p_mean   = p_two(mean_hub, mean_null),
       n_ge_mean = sum(abs(mean_null - mean(mean_null)) >= abs(mean_hub - mean(mean_null))),
       med_hub = med_hub, p_med = p_two(med_hub, med_null))
}

## Cached: the permutations are the slow part, so editing only the plots
## re-runs in seconds. Set TRUE (or delete the .rds) to recompute.
FORCE_RECOMPUTE <- FALSE
cache_file <- file.path(res_dir, "6_3_perm_cache.rds")
if (!FORCE_RECOMPUTE && file.exists(cache_file)) {
  message("Loading cached permutations from: ", cache_file)
  DR <- readRDS(cache_file)
} else {
  message("Computing permutations (slow part)...")
  DR <- setNames(lapply(names(net_colors), dist_perm), names(net_colors))
  saveRDS(DR, cache_file)
}

stats <- bind_rows(lapply(DR, function(d) tibble(
  network = d$net, n_hub = d$nh, B = B,
  ## separation (figures A / B)
  KS_D = d$D_hub,                                   # the tested statistic
  KS_Dplus = d$Dplus_hub, KS_Dminus = d$Dminus_hub, # direction, documented
  null_D_mean = mean(d$Dnull), null_D_max = max(d$Dnull),
  n_null_reaching_hub = d$n_ge, p_perm_D = d$p_perm,
  ## essentiality shift (figure C)
  mean_dep_hub = d$mean_hub, mean_dep_null = mean(d$mean_null),
  mean_dep_shift = d$mean_hub - mean(d$mean_null),
  null_mean_min = min(d$mean_null), null_mean_max = max(d$mean_null),
  n_null_beyond_hub = d$n_ge_mean, p_perm_mean = d$p_mean,
  median_dep_hub = d$med_hub, p_perm_median = d$p_med)))
vroom_write(stats, file.path(res_dir, "6_3_hub_KS_permutation_stats.tsv"), delim = "\t")
cat("\n=== Hub set vs", B, "random same-size sets ===\n")
print(as.data.frame(stats), row.names = FALSE)

## ── Theme ─────────────────────────────────────────────────────────────────────
base_theme <- theme_bw(base_family = "Times New Roman") +
  theme(axis.text        = element_text(size = 23, color = "black"),
        axis.title       = element_text(size = 27),
        axis.title.x     = element_text(margin = margin(t = 10)),
        axis.title.y     = element_text(margin = margin(r = 10)),
        axis.ticks       = element_line(linewidth = 0.8, color = "black"),
        plot.title       = element_text(size = 30, hjust = 0.5, margin = margin(b = 12)),
        panel.border     = element_rect(color = "black", linewidth = 1.0, fill = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
        plot.margin      = margin(16, 20, 12, 16))

## ── Figure A: hub vs non-hub CDF (classic axes) ───────────────────────────────
## TWO-SIDED KS, the same test figure B permutes (see header for why not D+).
##
## ks.test()'s C routine underflows to exactly 0 once D is large (A9), so fall
## back to the leading term of the asymptotic Kolmogorov series,
##   Q(t) = 2*sum_k (-1)^(k-1) exp(-2 k^2 t^2),  t = sqrt(ne)*D, ne = n1*n2/(n1+n2).
## Higher terms are O(exp(-8t^2)), negligible here. Verified: on A7 this
## reproduces ks.test's own two-sided p exactly (7.03e-10).
ks_p_twosided <- function(a, b) {
  p <- suppressWarnings(ks.test(a, b)$p.value)
  if (p > 0) return(p)
  n1 <- length(a); n2 <- length(b); ne <- n1 * n2 / (n1 + n2)
  D  <- suppressWarnings(unname(ks.test(a, b)$statistic))
  exp(log(2) - 2 * ne * D^2)
}

make_cdf_panel <- function(d) {
  col <- unname(net_colors[d$net])
  xx  <- filter(gene_tbl, network == d$net)
  hh  <- xx$degree >= hub_thr[[d$net]]
  p_ks  <- ks_p_twosided(xx$mean_dep[hh], xx$mean_dep[!hh])
  e     <- floor(log10(p_ks))
  D_lab <- sprintf("KS~italic(D) == %.3f", d$D_hub)
  p_lab <- sprintf("italic(p) == %.0f %%*%% 10^%d", p_ks / 10^e, e)
  df <- tibble(score = d$grid, hub = d$Fhub, rest = d$Fnonhub)
  x0 <- min(d$grid); w <- max(d$grid) - x0
  ggplot() +
    geom_line(data = df, aes(score, rest), color = "grey55", linewidth = 1.8) +
    geom_line(data = df, aes(score, hub),  color = col, linewidth = 2.2) +
    annotate("rect", xmin = x0 + 0.03 * w, xmax = x0 + 0.50 * w, ymin = 0.60, ymax = 0.985,
             fill = "white", color = "grey75", linewidth = 0.4, alpha = 0.92) +
    annotate("segment", x = x0 + 0.06 * w, xend = x0 + 0.14 * w, y = 0.925, yend = 0.925,
             color = col, linewidth = 2.2) +
    annotate("text", x = x0 + 0.16 * w, y = 0.925, label = "Hub set", hjust = 0,
             size = 7.6, fontface = "bold", family = "Times New Roman", color = col) +
    annotate("segment", x = x0 + 0.06 * w, xend = x0 + 0.14 * w, y = 0.840, yend = 0.840,
             color = "grey55", linewidth = 1.8) +
    annotate("text", x = x0 + 0.16 * w, y = 0.840, label = "non-hub set", hjust = 0,
             size = 7.6, family = "Times New Roman", color = "grey30") +
    annotate("text", x = x0 + 0.06 * w, y = 0.745, label = D_lab, parse = TRUE, hjust = 0,
             size = 7, family = "Times New Roman", color = "black") +
    annotate("text", x = x0 + 0.06 * w, y = 0.665, label = p_lab, parse = TRUE, hjust = 0,
             size = 7, family = "Times New Roman", color = "black") +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                       expand = expansion(mult = c(0.01, 0.02))) +
    coord_cartesian(clip = "off") +
    labs(title = d$net, x = "CRISPR gene effect",
         y = "Cumulative proportion of genes") +
    base_theme +
    theme(axis.text  = element_text(size = 27, color = "black"),
          axis.title = element_text(size = 31),
          plot.title = element_text(size = 34, hjust = 0.5, margin = margin(b = 12)))
}
pA <- make_cdf_panel(DR[["HPV-A7"]]) + make_cdf_panel(DR[["HPV-A9"]])
ggsave(file.path(res_dir, "6_3_A_hub_vs_nonhub_cdf.png"),
       pA, width = 16, height = 7.6, dpi = 300, bg = "white")
ggsave(file.path(res_dir, "6_3_A_hub_vs_nonhub_cdf.pdf"),
       pA, width = 16, height = 7.6, device = cairo_pdf, bg = "white")

## ── Figure B: permutation null of the KS D ────────────────────────────────────
make_Dnull_panel <- function(d) {
  col  <- unname(net_colors[d$net])
  dens <- density(d$Dnull); dd <- tibble(x = dens$x, y = dens$y); ymax <- max(dd$y)
  hub_lab <- sprintf("atop(bold('Hub set'), italic(D) == %.3f)", d$D_hub)
  p_lab   <- if (d$n_ge == 0) "italic(P) < 10^-4"
             else sprintf("italic(P) == %.3g", d$p_perm)
  xr <- range(c(d$Dnull, d$D_hub)); pad <- diff(xr) * 0.08
  xlim <- c(xr[1] - pad, xr[2] + pad)
  y_top <- ymax * 0.9
  ggplot() +
    geom_area(data = dd, aes(x, y), fill = "grey86", color = NA) +
    geom_line(data = dd, aes(x, y), color = "grey45", linewidth = 0.9) +
    annotate("text", x = mean(d$Dnull), y = ymax * 0.45,
             label = sprintf("%s random sets", format(B, big.mark = ",")),
             hjust = -0.05, size = 5.6, fontface = "italic",
             family = "Times New Roman", color = "grey40") +
    annotate("text", x = xlim[1], y = ymax * 1.16, label = p_lab, parse = TRUE,
             hjust = 0, vjust = 1, size = 6, family = "Times New Roman", color = "black") +
    annotate("segment", x = d$D_hub, xend = d$D_hub, y = 0, yend = y_top,
             color = col, linewidth = 1.0) +
    annotate("point", x = d$D_hub, y = y_top, shape = 18, size = 6, color = col) +
    annotate("text", x = d$D_hub, y = y_top * 0.70, label = hub_lab, parse = TRUE,
             hjust = 1.14, size = 6.2, family = "Times New Roman", color = col) +
    scale_x_continuous(limits = xlim) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
    labs(title = d$net, x = "Kolmogorov–Smirnov D", y = "Density") +
    base_theme
}
pB <- make_Dnull_panel(DR[["HPV-A7"]]) + make_Dnull_panel(DR[["HPV-A9"]])
ggsave(file.path(res_dir, "6_3_B_KS_permutation_null.png"),
       pB, width = 15, height = 6.8, dpi = 300, bg = "white")
ggsave(file.path(res_dir, "6_3_B_KS_permutation_null.pdf"),
       pB, width = 15, height = 6.8, device = cairo_pdf, bg = "white")

## ── Figure C: permutation null of the MEAN dependency ─────────────────────────
## What B cannot say: B shows the hub set separates more than any random set,
## but D is non-directional. Here the x axis IS the dependency score, so the
## direction is read off the figure: the hub set sits to the LEFT of the null
## (more negative = more essential). The test stays two-sided.
make_mean_panel <- function(d) {
  col  <- unname(net_colors[d$net])
  dens <- density(d$mean_null); dd <- tibble(x = dens$x, y = dens$y); ymax <- max(dd$y)
  hub_lab <- sprintf("atop(bold('Hub set'), bar(italic(x)) == %.3f)", d$mean_hub)
  p_lab   <- if (d$n_ge_mean == 0) "italic(P) < 10^-4"
             else sprintf("italic(P) == %.3g", d$p_mean)
  xr <- range(c(d$mean_null, d$mean_hub)); pad <- diff(xr) * 0.10
  xlim <- c(xr[1] - pad, xr[2] + pad)
  y_top <- ymax * 0.9
  ggplot() +
    geom_area(data = dd, aes(x, y), fill = "grey86", color = NA) +
    geom_line(data = dd, aes(x, y), color = "grey45", linewidth = 0.9) +
    # null label sits in the empty gap left of the cloud (the hub is further left,
    # at a different height), so it never runs off the panel edge
    annotate("text", x = min(d$mean_null), y = ymax * 0.30,
             label = sprintf("%s random sets", format(B, big.mark = ",")),
             hjust = 1.06, size = 5.6, fontface = "italic",
             family = "Times New Roman", color = "grey40") +
    annotate("text", x = xlim[2], y = ymax * 1.16, label = p_lab, parse = TRUE,
             hjust = 1, vjust = 1, size = 6, family = "Times New Roman", color = "black") +
    annotate("segment", x = d$mean_hub, xend = d$mean_hub, y = 0, yend = y_top,
             color = col, linewidth = 1.0) +
    annotate("point", x = d$mean_hub, y = y_top, shape = 18, size = 6, color = col) +
    annotate("text", x = d$mean_hub, y = y_top * 0.70, label = hub_lab, parse = TRUE,
             hjust = -0.14, size = 6.2, family = "Times New Roman", color = col) +
    scale_x_continuous(limits = xlim) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
    labs(title = d$net, x = "Mean CRISPR gene effect per gene set", y = "Density") +
    base_theme
}
pC <- make_mean_panel(DR[["HPV-A7"]]) + make_mean_panel(DR[["HPV-A9"]])
ggsave(file.path(res_dir, "6_3_C_mean_dependency_permutation_null.png"),
       pC, width = 15, height = 6.8, dpi = 300, bg = "white")
ggsave(file.path(res_dir, "6_3_C_mean_dependency_permutation_null.pdf"),
       pC, width = 15, height = 6.8, device = cairo_pdf, bg = "white")
## ── Figure D: B and C combined, one supplementary panel ───────────────────────
## Same panels as B and C, stacked. Row 1 asks how far the hub set separates
## from random sets; row 2 asks how essential it is compared with them.
row_title <- function(txt)
  ggplot() + theme_void() +
    annotate("text", x = 0, y = 0, label = txt, hjust = 0.5,
             size = 8, fontface = "bold", family = "Times New Roman")

# network name kept on the top row only; the columns are aligned, so repeating
# it on the bottom row would just be four titles competing with the two real ones
strip_title <- function(p) p + labs(title = NULL)

pD <- row_title("Separation from the rest of the network") /
      (make_Dnull_panel(DR[["HPV-A7"]]) + make_Dnull_panel(DR[["HPV-A9"]])) /
      row_title("Essentiality compared with random gene sets") /
      (strip_title(make_mean_panel(DR[["HPV-A7"]])) +
       strip_title(make_mean_panel(DR[["HPV-A9"]]))) +
      plot_layout(heights = c(0.10, 1, 0.10, 1))

ggsave(file.path(res_dir, "6_3_D_permutation_nulls_combined.png"),
       pD, width = 15, height = 10.5, dpi = 300, bg = "white")
ggsave(file.path(res_dir, "6_3_D_permutation_nulls_combined.pdf"),
       pD, width = 15, height = 10.5, device = cairo_pdf, bg = "white")

message("6_5 outputs saved to: ", res_dir)
