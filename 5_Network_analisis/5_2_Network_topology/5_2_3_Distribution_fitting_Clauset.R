################################################################################
# 5_2_3_Distribution_fitting_Clauset.R
#
# GOAL (relative only): identify the BEST RELATIVE DESCRIPTOR of each network's
# degree distribution among three candidate laws. We do NOT test whether any law
# is the *true* generating distribution (no absolute goodness-of-fit test) — only
# which one fits comparatively better.
#
# Method (tools from Clauset, Shalizi & Newman 2009; reference implementation
# poweRlaw, Gillespie 2015, J Stat Softw 64(2)):
#   For each model, applied even-handedly (not power-law-centric):
#     1. Maximum-likelihood fit; the tail start x_min is chosen by minimizing the
#        Kolmogorov-Smirnov distance (estimate_xmin).
#     2. KS distance (data vs fit) = how close each model is to the data.
#     3. Likelihood-ratio test (Vuong) between models on a common support = which
#        model wins head-to-head, and whether the win is significant.
#   No bootstrap GOF (that answers "is it the true law?", which is out of scope).
#
# Candidate set = the three laws implemented exactly by poweRlaw:
#   Power-law, Log-normal, Exponential.
#
# Output: 5_2_3_Distribution_benchmark.png / .pdf   (CCDF + fitted models, per network)
#         5_2_3_fit_table.tsv                (x_min, params, KS per model)
#         5_2_3_LR_table.tsv                 (pairwise likelihood-ratio tests)
################################################################################

suppressPackageStartupMessages({
  library(igraph)
  library(poweRlaw)
  library(ggplot2)
  library(scales)
  library(dplyr)
  library(tibble)
  library(patchwork)
})

set.seed(123)
base_dir <- "~/CESC_Network/5_Network_analisis/5_2_Network_topology"

################################################################################
# 1) Load networks and degree sequences
################################################################################

A7_HPV <- readRDS("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_2_4_Subnetwork_A7_max.rds")
A9_HPV <- readRDS("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_2_4_Subnetwork_A9_max.rds")

deg <- list(
  `HPV-A7` = as.integer(igraph::degree(A7_HPV, mode = "all")),
  `HPV-A9` = as.integer(igraph::degree(A9_HPV, mode = "all"))
)

models  <- list(`Power-law` = displ, `Log-normal` = dislnorm, `Exponential` = disexp)
mcolors <- c(`Power-law` = "#7B1FA2", `Log-normal` = "#003366", `Exponential` = "#009900")
mlines  <- c(`Power-law` = "longdash", `Log-normal` = "solid", `Exponential` = "dotted")

################################################################################
# 2) Fit each model (own KS-optimal x_min) + pairwise likelihood-ratio
################################################################################

fit_one <- function(d, ctor) {
  d <- as.integer(d[d > 0])
  m <- ctor$new(d); est <- estimate_xmin(m); m$setXmin(est)
  list(m = m, xmin = m$getXmin(), pars = m$getPars(), D = est$gof,
       ntail = sum(d >= m$getXmin()), n = length(d))
}

analyse <- function(net) {
  d <- as.integer(deg[[net]][deg[[net]] > 0])
  fits <- lapply(models, function(ctor) fit_one(d, ctor))
  nms  <- names(fits)

  fit_tbl <- bind_rows(lapply(nms, function(nm) {
    f <- fits[[nm]]
    tibble(Network = net, Model = nm, xmin = f$xmin,
           Params = paste(sprintf("%.4g", f$pars), collapse = ", "),
           KS_D = round(f$D, 4), N_tail = f$ntail,
           Pct_tail = round(100 * f$ntail / f$n, 1))
  })) |> arrange(KS_D)

  ## Pairwise likelihood-ratio (Vuong) at a common x_min = max of the pair
  lr_tbl <- bind_rows(lapply(seq_len(nrow(t(combn(nms, 2)))), function(k) {
    pair <- combn(nms, 2)[, k]
    xm <- max(fits[[pair[1]]]$xmin, fits[[pair[2]]]$xmin)
    a <- models[[pair[1]]]$new(d); a$setXmin(xm); a$setPars(estimate_pars(a)$pars)
    b <- models[[pair[2]]]$new(d); b$setXmin(xm); b$setPars(estimate_pars(b)$pars)
    cmp <- compare_distributions(a, b)
    tibble(Network = net, Comparison = paste(pair[1], "vs", pair[2]),
           common_xmin = xm, Vuong_stat = round(cmp$test_statistic, 2),
           p_value = signif(cmp$p_two_sided, 3),
           Favored = ifelse(abs(cmp$test_statistic) < 1.96 |
                              cmp$p_two_sided > 0.05, "tie (n.s.)",
                            ifelse(cmp$test_statistic > 0, pair[1], pair[2])))
  }))

  ## Best relative descriptor = smallest KS, not significantly beaten in any LR
  best <- fit_tbl$Model[1]
  beaten <- lr_tbl |> filter(grepl(best, Comparison), Favored != "tie (n.s.)", Favored != best)
  best_note <- if (nrow(beaten) == 0)
    sprintf("%s (smallest KS; wins or ties every likelihood-ratio test)", best)
  else sprintf("%s smallest KS, but likelihood-ratio favors %s",
               best, paste(unique(beaten$Favored), collapse = "/"))

  list(net = net, d = d, fits = fits, fit_tbl = fit_tbl, lr_tbl = lr_tbl,
       best = best, best_note = best_note)
}

res <- lapply(names(deg), analyse); names(res) <- names(deg)

for (net in names(res)) {
  cat("\n=====================", net, "=====================\n")
  print(as.data.frame(res[[net]]$fit_tbl), row.names = FALSE)
  cat("--- likelihood-ratio (Vuong) ---\n")
  print(as.data.frame(res[[net]]$lr_tbl), row.names = FALSE)
  cat("Best relative descriptor:", res[[net]]$best_note, "\n")
}

write.table(bind_rows(lapply(res, `[[`, "fit_tbl")),
            file.path(base_dir, "5_2_3_fit_table.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(bind_rows(lapply(res, `[[`, "lr_tbl")),
            file.path(base_dir, "5_2_3_LR_table.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

################################################################################
# 3) Benchmark plot — empirical CCDF + each fitted model (drawn from its x_min)
#    Model CCDF is scaled by the empirical fraction >= x_min so it overlays the
#    data at x_min (standard Clauset survival-plot convention).
################################################################################

emp_ccdf <- function(d) {
  tab <- table(d); k <- as.integer(names(tab)); ord <- order(k)
  k <- k[ord]; cnt <- as.numeric(tab[ord])
  tibble(k = k, CCDF = rev(cumsum(rev(cnt))) / sum(cnt))
}

make_panel <- function(r) {
  d <- r$d
  emp <- emp_ccdf(d)

  curves <- bind_rows(lapply(names(models), function(nm) {
    f <- r$fits[[nm]]; ks <- seq(f$xmin, max(d))
    scale <- mean(d >= f$xmin)                       # empirical P(K >= xmin)
    Smod  <- dist_cdf(f$m, ks, lower_tail = FALSE)   # model survival | X >= xmin
    tibble(k = ks, CCDF = scale * Smod, Model = nm)
  })) |> mutate(Model = factor(Model, levels = names(models)))

  fb  <- r$fits[[r$best]]
  pfmt <- switch(r$best,
    "Log-normal"  = sprintf("μ = %.2f, σ = %.2f", fb$pars[1], fb$pars[2]),
    "Power-law"   = sprintf("α = %.2f", fb$pars[1]),
    "Exponential" = sprintf("λ = %.3f", fb$pars[1]))
  ann <- sprintf("Best-supported model: %s\n%s\nx_min = 10^%.2f (k = %d)",
                 r$best, pfmt, log10(fb$xmin), fb$xmin)

  ggplot() +
    geom_point(data = emp, aes(k, CCDF), color = "grey55", alpha = 0.55, size = 2.4) +
    geom_line(data = curves, aes(k, CCDF, color = Model, linetype = Model), linewidth = 1.6) +
    annotate("label", x = min(emp$k), y = min(emp$CCDF), label = ann,
             hjust = 0, vjust = 0, family = "Times New Roman", size = 8,
             label.size = 0.4, label.padding = unit(0.35, "lines"),
             fill = "white", color = "grey15", lineheight = 1.15) +
    scale_color_manual(values = mcolors) +
    scale_linetype_manual(values = mlines) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides = "bl") +
    labs(title = paste("Degree distribution —", r$net),
         x = "Degree (k)", y = expression(P(K >= k)),
         color = "Model", linetype = "Model") +
    theme_bw(base_family = "Times New Roman") +
    theme(
      plot.title       = element_text(size = 33, face = "bold", hjust = 0.5),
      plot.subtitle    = element_text(size = 22, color = "grey20", hjust = 0.5,
                                      lineheight = 1.3, margin = margin(t = 6, b = 8)),
      axis.title       = element_text(size = 34, face = "bold"),
      axis.text        = element_text(size = 30),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_line(color = "grey95"),
      legend.position  = "top",
      legend.title     = element_text(size = 36, face = "bold"),
      legend.text      = element_text(size = 34),
      legend.key.width = unit(2.6, "cm"),
      legend.key.size  = unit(1.4, "cm")
    )
}

p_bench <- make_panel(res[["HPV-A7"]]) + make_panel(res[["HPV-A9"]]) +
  plot_layout(ncol = 2, guides = "collect") & theme(legend.position = "bottom")

ggsave(file.path(base_dir, "5_2_3_Distribution_benchmark.png"),
       p_bench, width = 17, height = 8.8, dpi = 300, bg = "white")
ggsave(file.path(base_dir, "5_2_3_Distribution_benchmark.pdf"),
       p_bench, width = 17, height = 8.8, device = cairo_pdf, bg = "white")
message("Saved: 5_2_3_Distribution_benchmark.png / .pdf")
