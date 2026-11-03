################################################################################
# figure_S2_combined.R
#
# CESC + Uganda combined view of Figure S2 (degree-distribution model fit).
# Does NOT touch the original CESC-only figure. Grouped by ANALYSIS TYPE:
# each row is one clade (A7, A9), CESC and Uganda side by side under the
# same tag letter.
#
#   A) HPV-A7  [CESC | Uganda]
#   B) HPV-A9  [CESC | Uganda]
#
# -> figures/validation_uganda/Figure_S2_degree_distribution_fit_uganda.png/.pdf
################################################################################

suppressPackageStartupMessages({
  library(igraph); library(poweRlaw); library(ggplot2); library(scales)
  library(dplyr); library(tibble); library(patchwork)
})

set.seed(123)
fig_dir <- "~/CESC_Network/figures/validation_uganda"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

models  <- list(`Power-law` = displ, `Log-normal` = dislnorm, `Exponential` = disexp)
mcolors <- c(`Power-law` = "#7B1FA2", `Log-normal` = "#003366", `Exponential` = "#009900")
mlines  <- c(`Power-law` = "longdash", `Log-normal` = "solid", `Exponential` = "dotted")

fit_one <- function(d, ctor) {
  d <- as.integer(d[d > 0]); m <- ctor$new(d); est <- estimate_xmin(m); m$setXmin(est)
  list(m = m, xmin = m$getXmin(), pars = m$getPars(), D = est$gof, ntail = sum(d >= m$getXmin()), n = length(d))
}
analyse <- function(d, net) {
  d <- as.integer(d[d > 0]); fits <- lapply(models, function(ctor) fit_one(d, ctor)); nms <- names(fits)
  fit_tbl <- bind_rows(lapply(nms, function(nm) {
    f <- fits[[nm]]
    tibble(Network = net, Model = nm, xmin = f$xmin, KS_D = round(f$D, 4))
  })) |> arrange(KS_D)
  best <- fit_tbl$Model[1]
  list(net = net, d = d, fits = fits, best = best)
}
emp_ccdf <- function(d) {
  tab <- table(d); k <- as.integer(names(tab)); ord <- order(k); k <- k[ord]; cnt <- as.numeric(tab[ord])
  tibble(k = k, CCDF = rev(cumsum(rev(cnt))) / sum(cnt))
}
make_panel <- function(r, title_tag) {
  d <- r$d; emp <- emp_ccdf(d)
  curves <- bind_rows(lapply(names(models), function(nm) {
    f <- r$fits[[nm]]; ks <- seq(f$xmin, max(d)); scale <- mean(d >= f$xmin)
    Smod <- dist_cdf(f$m, ks, lower_tail = FALSE)
    tibble(k = ks, CCDF = scale * Smod, Model = nm)
  })) |> mutate(Model = factor(Model, levels = names(models)))
  fb <- r$fits[[r$best]]
  pfmt <- switch(r$best, "Log-normal" = sprintf("μ = %.2f, σ = %.2f", fb$pars[1], fb$pars[2]),
                "Power-law" = sprintf("α = %.2f", fb$pars[1]), "Exponential" = sprintf("λ = %.3f", fb$pars[1]))
  ann <- sprintf("Best: %s\n%s\nx_min = %d", r$best, pfmt, fb$xmin)
  ggplot() +
    geom_point(data = emp, aes(k, CCDF), color = "grey55", alpha = 0.55, size = 1.6) +
    geom_line(data = curves, aes(k, CCDF, color = Model, linetype = Model), linewidth = 1.1) +
    annotate("label", x = min(emp$k), y = min(emp$CCDF), label = ann, hjust = 0, vjust = 0,
             family = "Times New Roman", size = 4.6, label.padding = unit(0.28, "lines"),
             fill = "white", color = "grey15", lineheight = 1.05) +
    scale_color_manual(values = mcolors) + scale_linetype_manual(values = mlines) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides = "bl") +
    labs(title = title_tag, x = "Degree (k)", y = expression(P(K >= k)), color = "Model", linetype = "Model") +
    theme_bw(base_family = "Times New Roman") +
    theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 16, face = "bold"), axis.text = element_text(size = 13),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          legend.position = "top", legend.title = element_text(size = 15, face = "bold"),
          legend.text = element_text(size = 13), legend.key.width = unit(1.4, "cm"))
}

run_cohort <- function(base_dir) {
  A7_HPV <- readRDS(file.path(base_dir, "5_2_4_Subnetwork_A7_max.rds"))
  A9_HPV <- readRDS(file.path(base_dir, "5_2_4_Subnetwork_A9_max.rds"))
  deg <- list(`HPV-A7` = as.integer(igraph::degree(A7_HPV, mode = "all")),
             `HPV-A9` = as.integer(igraph::degree(A9_HPV, mode = "all")))
  list(A7 = analyse(deg[["HPV-A7"]], "HPV-A7"), A9 = analyse(deg[["HPV-A9"]], "HPV-A9"))
}

cesc   <- run_cohort("~/CESC_Network/05_topology/results")
uganda <- run_cohort("~/CESC_Network/05_topology/validation_uganda/results")

rowA <- wrap_elements(full = (make_panel(cesc$A7, "HPV-A7 — CESC") | make_panel(uganda$A7, "HPV-A7 — Uganda")))
rowB <- wrap_elements(full = (make_panel(cesc$A9, "HPV-A9 — CESC") | make_panel(uganda$A9, "HPV-A9 — Uganda")))

p_bench <- rowA / rowB +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 22, face = "bold", family = "Times New Roman"))

ggsave(file.path(fig_dir, "Figure_S2_degree_distribution_fit_uganda.png"), p_bench, width = 15, height = 13, dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "Figure_S2_degree_distribution_fit_uganda.pdf"), p_bench, width = 15, height = 13, device = cairo_pdf, bg = "white")
message("Saved: Figure_S2_degree_distribution_fit_uganda.png / .pdf  (A=HPV-A7, B=HPV-A9; each CESC | Uganda)")
