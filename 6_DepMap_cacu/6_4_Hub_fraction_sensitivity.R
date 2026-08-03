### ============================================================
### 6_4_Hub_fraction_sensitivity.R
###
### QUESTION: 6_3 fixed the hub set at the top 10% by degree and showed
### (A/B) that this set SEPARATES from the rest of the network (KS D) and
### (C) that it is SHIFTED toward CRISPR essentiality more than any random
### same-size set. Both effects were overwhelming at 10%.
###
### Here we ask the sensitivity question: as we LOOSEN the hub definition
### (label a larger and larger fraction of genes as "hubs", from 2% up to
### 95% by degree), at what fraction do these significant effects vanish?
###
### Two effects tracked at every fraction f (hub = top f by degree):
###
###   1) SEPARATION   two-sided two-sample KS between the hub set and its
###      complement (the same test as 6_3 A/B). p is the asymptotic
###      two-sided KS p, matching the value printed on figure 6_3_A.
###
###   2) ESSENTIALITY SHIFT   mean CRISPR dependency of the hub set vs a
###      size-matched random gene set (the same idea as 6_3 C). Two-sided
###      permutation p, B resamples.
###
### NULL UNIVERSE. 6_3 C draws the random sets from the NON-hub genes. That
### null degenerates once f approaches 50% (the non-hub pool becomes as small
### as the hub set, so there is essentially one possible random set and the
### null collapses). To keep a valid null across the WHOLE sweep, the random
### sets here are drawn from ALL network genes. At f = 10% this is
### indistinguishable from 6_3 (hubs are only 10% of the pool), and it stays
### valid as f grows. This is stated so the 10% row reconciles with 6_3.
###
### Outputs:
###   6_4_hub_fraction_sensitivity.tsv      per-fraction D, shift, p-values
###   6_4_hub_fraction_crossings.tsv        the fraction where each effect
###                                         first loses significance (p>=0.05)
###   6_4_A_hub_fraction_significance.png    -log10(p) vs hub fraction
###   6_4_B_hub_fraction_effectsize.png      KS D and mean shift vs fraction
### ============================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(vroom)
  library(parallel)
})

dep_dir <- "~/CESC_Network/6_DepMap_cacu"
res_dir <- file.path(dep_dir, "results")

REF_FRAC <- 0.10          # the hub definition used everywhere else (top 10%)
ALPHA    <- 0.05
B        <- 2000          # permutation resamples for the mean-dependency null
set.seed(123)
net_colors <- c("HPV-A7" = "#003366", "HPV-A9" = "#CC2200")
FRACS <- seq(0.02, 0.95, by = 0.01)

gene_tbl <- vroom(file.path(res_dir, "6_2_gene_dep_degree_table.tsv"),
                  show_col_types = FALSE)

## Asymptotic two-sided KS p, with the same underflow fallback as 6_3_A:
## ks.test()'s C routine returns exactly 0 once D is large, so drop to the
## leading term of the Kolmogorov series. Verified to reproduce ks.test on A7.
ks_p_twosided <- function(a, b) {
  kt <- suppressWarnings(ks.test(a, b))
  p  <- kt$p.value
  D  <- unname(kt$statistic)
  if (p > 0) return(list(D = D, p = p))
  n1 <- length(a); n2 <- length(b); ne <- n1 * n2 / (n1 + n2)
  list(D = D, p = exp(log(2) - 2 * ne * D^2))
}

## ── Sweep: one job per (network, fraction), run across cores ──────────────────
## Each of the 94 fractions x 2 networks is an independent permutation job, so
## the sweep is embarrassingly parallel. mclapply spreads the jobs over the
## available cores (single-core it would take several minutes; parallel ~1 min).
netdata <- lapply(names(net_colors), function(net) {
  x <- filter(gene_tbl, network == net)
  list(net = net, dep = x$mean_dep, deg = x$degree)
})
names(netdata) <- names(net_colors)
jobs <- expand.grid(net = names(net_colors), f = FRACS, stringsAsFactors = FALSE)

run_job <- function(i) {
  net <- jobs$net[i]; f <- jobs$f[i]
  nd  <- netdata[[net]]; dep <- nd$dep; deg <- nd$deg
  N   <- length(dep); grand <- mean(dep)
  thr <- quantile(deg, 1 - f, type = 7)
  hub <- deg >= thr; nh <- sum(hub)
  if (min(nh, N - nh) < 30) return(NULL)            # need both sides populated

  ## 1) separation: hub vs complement (two-sided KS, asymptotic p)
  ks <- ks_p_twosided(dep[hub], dep[!hub])

  ## 2) essentiality shift: hub mean vs random same-size sets from ALL genes
  mean_hub <- mean(dep[hub])
  null <- replicate(B, mean(dep[sample.int(N, nh)]))
  mu   <- mean(null)
  p_ess <- (1 + sum(abs(null - mu) >= abs(mean_hub - mu))) / (B + 1)

  tibble(network = net, frac = f, n_hub = nh, frac_real = nh / N,
         KS_D = ks$D, p_sep = ks$p,
         mean_dep_hub = mean_hub, grand_mean = grand,
         mean_shift = mean_hub - grand, p_ess = p_ess)
}

## Cached: the permutation sweep is the slow part. Set FORCE_RECOMPUTE <- TRUE
## (or delete the tsv) to recompute; editing only the plots re-runs in seconds.
FORCE_RECOMPUTE <- FALSE
cache_tsv <- file.path(res_dir, "6_4_hub_fraction_sensitivity.tsv")
if (!FORCE_RECOMPUTE && file.exists(cache_tsv)) {
  message("Loading cached sweep from: ", cache_tsv)
  sweep <- vroom(cache_tsv, show_col_types = FALSE)
} else {
  n_cores <- max(1, detectCores() - 1)
  message(sprintf("Sweeping %d jobs over %d cores...", nrow(jobs), n_cores))
  ## L'Ecuyer-CMRG gives reproducible independent RNG streams across workers
  RNGkind("L'Ecuyer-CMRG"); set.seed(123)
  sweep <- bind_rows(mclapply(seq_len(nrow(jobs)), run_job, mc.cores = n_cores))
  vroom_write(sweep, cache_tsv, delim = "\t")
}

## ── Crossing points: first fraction (> reference) where p >= ALPHA ────────────
first_cross <- function(d, pcol) {
  d <- d %>% arrange(frac_real) %>% filter(frac_real >= REF_FRAC)
  hit <- which(d[[pcol]] >= ALPHA)
  if (length(hit) == 0) return(NA_real_)
  d$frac_real[hit[1]]
}
crossings <- sweep %>% group_by(network) %>%
  group_modify(~ tibble(sep_lost_at = first_cross(.x, "p_sep"),
                        ess_lost_at = first_cross(.x, "p_ess"))) %>%
  ungroup()
vroom_write(crossings, file.path(res_dir, "6_4_hub_fraction_crossings.tsv"), delim = "\t")

cat("\n=== Fraction at which each effect first loses significance (p >=", ALPHA, ") ===\n")
crossings %>%
  mutate(across(c(sep_lost_at, ess_lost_at),
                ~ ifelse(is.na(.), "never in range", sprintf("%.0f%%", . * 100)))) %>%
  as.data.frame() %>% print(row.names = FALSE)

## ── Long table for plotting significance, capped for readability ──────────────
CAP <- 8   # display cap on -log10(p); the significant plateau runs far above it
sig_long <- sweep %>%
  transmute(network, frac_real,
            `Separation (KS, hub vs rest)`                  = p_sep,
            `Essentiality shift (hub vs random gene sets)`  = p_ess) %>%
  pivot_longer(-c(network, frac_real), names_to = "effect", values_to = "p") %>%
  mutate(neglogp = pmin(-log10(p), CAP),
         effect  = factor(effect, levels = c("Separation (KS, hub vs rest)",
                                             "Essentiality shift (hub vs random gene sets)")))

cross_long <- crossings %>%
  pivot_longer(-network, names_to = "which", values_to = "frac_real") %>%
  mutate(effect = ifelse(which == "sep_lost_at",
                         "Separation (KS, hub vs rest)",
                         "Essentiality shift (hub vs random gene sets)"),
         effect = factor(effect, levels = levels(sig_long$effect))) %>%
  filter(!is.na(frac_real))

base_theme <- theme_bw(base_family = "Times New Roman") +
  theme(axis.text        = element_text(size = 16, color = "black"),
        axis.title       = element_text(size = 19, face = "bold"),
        plot.title       = element_text(size = 22, face = "bold", hjust = 0.5),
        plot.subtitle    = element_text(size = 14, hjust = 0.5, color = "grey30"),
        strip.text       = element_text(size = 15, face = "bold"),
        strip.background  = element_rect(fill = "grey93", color = "grey60"),
        legend.position  = "none",
        panel.grid.minor = element_blank(),
        panel.border     = element_rect(color = "black", linewidth = 0.9, fill = NA))

## ── Figure A: significance vs hub fraction ────────────────────────────────────
pA <- ggplot(sig_long, aes(frac_real * 100, neglogp, color = network)) +
  geom_hline(yintercept = -log10(ALPHA), linetype = "dashed",
             color = "grey35", linewidth = 0.8) +
  geom_vline(xintercept = REF_FRAC * 100, linetype = "dotted",
             color = "grey55", linewidth = 0.7) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 1.6, alpha = 0.8) +
  geom_vline(data = cross_long, aes(xintercept = frac_real * 100, color = network),
             linetype = "solid", linewidth = 0.8, alpha = 0.5) +
  geom_text(data = cross_long,
            aes(x = frac_real * 100, y = CAP * 0.9,
                label = sprintf("lost at %.0f%%", frac_real * 100)),
            inherit.aes = FALSE, family = "Times New Roman", fontface = "bold",
            size = 4.6, hjust = -0.06, color = "grey15") +
  annotate("text", x = REF_FRAC * 100, y = -log10(ALPHA) + 0.35,
           label = "p = 0.05", family = "Times New Roman", size = 4.2,
           color = "grey35", hjust = -0.05) +
  facet_grid(effect ~ network) +
  scale_color_manual(values = net_colors) +
  scale_x_continuous(breaks = seq(0, 90, 20)) +
  labs(title = "At what hub fraction do the essentiality effects vanish?",
       subtitle = sprintf("Dotted line = the top-10%% hub definition used elsewhere; -log10(p) capped at %d for display", CAP),
       x = "Genes labelled as hubs (top % by degree)",
       y = expression(-log[10](p))) +
  base_theme
ggsave(file.path(res_dir, "6_4_A_hub_fraction_significance.png"),
       pA, width = 14, height = 9, dpi = 300, bg = "white")
ggsave(file.path(res_dir, "6_4_A_hub_fraction_significance.pdf"),
       pA, width = 14, height = 9, device = cairo_pdf, bg = "white")

## ── Figure B: effect size vs hub fraction ─────────────────────────────────────
eff_long <- sweep %>%
  transmute(network, frac_real,
            `Separation  (KS D)`                       = KS_D,
            `Essentiality shift  (mean dependency)`    = mean_shift) %>%
  pivot_longer(-c(network, frac_real), names_to = "effect", values_to = "value") %>%
  mutate(effect = factor(effect, levels = c("Separation  (KS D)",
                                           "Essentiality shift  (mean dependency)")))

pB <- ggplot(eff_long, aes(frac_real * 100, value, color = network)) +
  geom_hline(data = tibble(effect = factor("Essentiality shift  (mean dependency)",
                                           levels = levels(eff_long$effect)), y = 0),
             aes(yintercept = y), linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = REF_FRAC * 100, linetype = "dotted",
             color = "grey55", linewidth = 0.7) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 1.6, alpha = 0.8) +
  facet_grid(effect ~ network, scales = "free_y") +
  scale_color_manual(values = net_colors) +
  scale_x_continuous(breaks = seq(0, 90, 20)) +
  labs(title = "Effect size decays as the hub definition is loosened",
       subtitle = "Dotted line = the top-10% hub definition used elsewhere",
       x = "Genes labelled as hubs (top % by degree)", y = NULL) +
  base_theme
ggsave(file.path(res_dir, "6_4_B_hub_fraction_effectsize.png"),
       pB, width = 14, height = 9, dpi = 300, bg = "white")
ggsave(file.path(res_dir, "6_4_B_hub_fraction_effectsize.pdf"),
       pB, width = 14, height = 9, device = cairo_pdf, bg = "white")

## ── Figure C: mean-dependency section only (A7 | A9) ──────────────────────────
pC <- ggplot(sweep, aes(frac_real * 100, mean_shift, color = network)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = REF_FRAC * 100, linetype = "dotted",
             color = "grey55", linewidth = 0.7) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 1.8, alpha = 0.85) +
  facet_wrap(~ network) +
  scale_color_manual(values = net_colors) +
  scale_x_continuous(breaks = seq(0, 90, 20)) +
  labs(title = "Mean CRISPR dependency of the hub set across hub definitions",
       subtitle = "Dotted line = the top-10% hub definition used elsewhere",
       x = "Genes labelled as hubs (top % by degree)",
       y = "Mean dependency shift vs the whole network") +
  base_theme
ggsave(file.path(res_dir, "6_4_C_mean_dependency_vs_fraction.png"),
       pC, width = 14, height = 6.2, dpi = 300, bg = "white")
ggsave(file.path(res_dir, "6_4_C_mean_dependency_vs_fraction.pdf"),
       pC, width = 14, height = 6.2, device = cairo_pdf, bg = "white")

message("✅ 6_7 outputs saved to: ", res_dir)
