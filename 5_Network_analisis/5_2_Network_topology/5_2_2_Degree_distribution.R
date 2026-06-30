
################################################################################
# Degree Distribution Analysis for Cervical Cancer Co-expression Networks
# Author: Joel Ruiz Hernandez
#
# Model menu (discrete MLE, xmin = 1, AIC/BIC comparable across all models):
#   1) Power-law                 (Barabasi & Albert 1999; Clauset et al. 2009)
#   2) Log-normal                (Limpert et al. 2001; Khanin & Wit 2006;
#                                  Broido & Clauset 2019)
#   3) Discrete Weibull          (Amaral et al. 2000; Clauset et al. 2009)
#   4) Exponential               (Amaral et al. 2000; Clauset et al. 2009)
#
# Power-law, log-normal and exponential are fit with poweRlaw (Gillespie 2015,
# J Stat Softw 64(2)) at xmin = 1. Discrete Weibull is fit by direct
# discrete-likelihood maximization (no closed form in poweRlaw).
# A Vuong likelihood-ratio test (poweRlaw::compare_distributions) contrasts
# Power-law vs Log-normal directly.
################################################################################

suppressPackageStartupMessages({
  library(igraph)
  library(ggplot2)
  library(scales)
  library(dplyr)
  library(tibble)
  library(poweRlaw)
  library(patchwork)
})

################################################################################
# 1) Load graphs
################################################################################

A7_HPV <- readRDS("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_2_4_Subnetwork_A7_max.rds")

# A9 recortada al mismo número de edges que A7 max (consenso = ecount A7)
CONSENSUS_N <- igraph::ecount(A7_HPV)   # 117 184
el_A9 <- readRDS("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_1_2_A9_edgelist_individual_allnodes.rds") |>
  dplyr::arrange(dplyr::desc(MI)) |> dplyr::slice(1:CONSENSUS_N)
A9_HPV <- igraph::graph_from_data_frame(el_A9[, c("GenA", "GenB")], directed = FALSE) |> igraph::simplify()

################################################################################
# 2) Extract degree sequences
################################################################################

deg_A7 <- as.integer(igraph::degree(A7_HPV, mode = "all"))
deg_A9 <- as.integer(igraph::degree(A9_HPV, mode = "all"))

################################################################################
# 3) Statistical comparison: KS + Wilcoxon
################################################################################

ks_result <- ks.test(deg_A7, deg_A9)
wx_result  <- wilcox.test(deg_A7, deg_A9, alternative = "two.sided", conf.int = TRUE)

print(ks_result)
print(wx_result)

tibble(
  Network  = c("HPV_A7", "HPV_A9"),
  N_nodes  = c(length(deg_A7), length(deg_A9)),
  Min      = c(min(deg_A7),    min(deg_A9)),
  Median   = c(median(deg_A7), median(deg_A9)),
  Mean     = round(c(mean(deg_A7), mean(deg_A9)), 2),
  Max      = c(max(deg_A7),    max(deg_A9)),
  SD       = round(c(sd(deg_A7), sd(deg_A9)), 2)
) |> print()

################################################################################
# 4) Distribution fitting — discrete MLE + AIC, 4-model menu
################################################################################

# ---- Discrete Stretched Exponential (Clauset et al. 2009, Section 4.3) ------
# PMF: P(k) = exp(-lambda * k^beta) / Z,   lambda > 0, beta > 0
# Z = sum_{j=1}^{inf} exp(-lambda * j^beta)  (computed by truncated sum)
# Reduces to discrete exponential when beta = 1

log_Z_strexp <- function(lambda, beta, k_max) {
  lv <- -lambda * seq_len(k_max)^beta
  m  <- max(lv)
  m + log(sum(exp(lv - m)))
}

negll_strexp <- function(par, x) {
  lambda <- exp(par[1])
  beta   <- exp(par[2])
  if (!is.finite(lambda) || !is.finite(beta)) return(1e10)
  k_max <- max(max(x) * 5L, 5000L)
  lZ <- log_Z_strexp(lambda, beta, k_max)
  ll <- -lambda * sum(x^beta) - length(x) * lZ
  -ll
}

fit_strexp <- function(x) {
  best <- NULL
  for (init_l in c(0.005, 0.05, 0.5)) {
    for (init_b in c(0.5, 1.0, 1.5)) {
      opt <- try(optim(c(log(init_l), log(init_b)), negll_strexp, x = x,
                       method = "Nelder-Mead", control = list(maxit = 5000)),
                 silent = TRUE)
      if (!inherits(opt, "try-error") && (is.null(best) || opt$value < best$value))
        best <- opt
    }
  }
  lambda_hat <- exp(best$par[1])
  beta_hat   <- exp(best$par[2])
  ll <- -best$value
  list(lambda = lambda_hat, beta = beta_hat, loglik = ll,
       aic = 2 * 2 - 2 * ll, bic = 2 * log(length(x)) - 2 * ll)
}

pmf_strexp <- function(k, lambda, beta, k_max = NULL) {
  if (is.null(k_max)) k_max <- max(k) * 5L
  lZ <- log_Z_strexp(lambda, beta, k_max)
  exp(-lambda * k^beta - lZ)
}

# ---- Main fitting routine ---------------------------------------------------
fit_network <- function(deg_vec, network_name) {
  message("\n  Fitting: ", network_name)
  d <- as.integer(deg_vec[deg_vec > 0])
  n <- length(d)

  # 1) Power-law (discrete, xmin = 1)
  m_pl <- displ$new(d); m_pl$setXmin(1); m_pl$setPars(estimate_pars(m_pl)$pars)
  ll_pl <- dist_ll(m_pl)

  # 2) Log-normal (discrete, xmin = 1)
  m_ln <- dislnorm$new(d); m_ln$setXmin(1); m_ln$setPars(estimate_pars(m_ln)$pars)
  ll_ln <- dist_ll(m_ln)

  # 3) Exponential (discrete, xmin = 1)
  m_exp <- disexp$new(d); m_exp$setXmin(1); m_exp$setPars(estimate_pars(m_exp)$pars)
  ll_exp <- dist_ll(m_exp)

  # 4) Stretched Exponential (Clauset et al. 2009)
  strexp <- fit_strexp(d)

  aic_tbl <- tibble(
    Distribution = c("Power-law", "Log-normal", "Stretched Exponential", "Exponential"),
    LogLik = c(ll_pl, ll_ln, strexp$loglik, ll_exp),
    AIC    = c(2 * 1 - 2 * ll_pl, 2 * 2 - 2 * ll_ln,
               strexp$aic, 2 * 1 - 2 * ll_exp),
    BIC    = c(1 * log(n) - 2 * ll_pl, 2 * log(n) - 2 * ll_ln,
               strexp$bic, 1 * log(n) - 2 * ll_exp)
  ) |>
    arrange(AIC) |>
    mutate(DeltaAIC = round(AIC - min(AIC), 1))

  print(aic_tbl)

  # Vuong likelihood-ratio test: Power-law vs Log-normal
  vuong <- compare_distributions(m_pl, m_ln)
  message(sprintf(
    "  Vuong test (Power-law vs Log-normal): statistic = %.3f, p (two-sided) = %.3g",
    vuong$test_statistic, vuong$p_two_sided
  ))

  invisible(list(
    m_pl = m_pl, m_ln = m_ln, m_exp = m_exp, strexp = strexp,
    aic = aic_tbl, vuong = vuong,
    deg_int = d, n = n,
    lnorm_pars = m_ln$pars,
    mu    = round(m_ln$pars[1], 3),
    sigma = round(m_ln$pars[2], 3)
  ))
}

fits_A7 <- fit_network(deg_A7, "HPV_A7")
fits_A9 <- fit_network(deg_A9, "HPV_A9")

# Save AIC/BIC comparison table for both networks
aic_out <- bind_rows(
  fits_A7$aic |> mutate(Network = "HPV-A7", .before = 1),
  fits_A9$aic |> mutate(Network = "HPV-A9", .before = 1)
)
write.table(
  aic_out,
  "~/CESC_Network/5_Network_analisis/5_2_Network_topology/5_2_2_Distribution_AIC_comparison.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE
)

################################################################################
# 5) Plot 1 — Distribution benchmark
#    Empirical PMF + 5 fitted discrete PMFs per network, log-log scale
################################################################################

make_benchmark_panel <- function(fits, deg_vec, network_name) {

  d  <- as.integer(deg_vec[deg_vec > 0])
  ks <- sort(unique(d))

  emp_df <- as.data.frame(table(d))
  names(emp_df) <- c("k", "count")
  emp_df <- emp_df |>
    mutate(k = as.integer(as.character(k)),
           Probability = count / sum(count))

  curves_df <- bind_rows(
    tibble(k = ks, Probability = dist_pdf(fits$m_pl, ks),  Model = "Power-law"),
    tibble(k = ks, Probability = dist_pdf(fits$m_ln, ks),  Model = "Log-normal"),
    tibble(k = ks, Probability = pmf_strexp(ks, fits$strexp$lambda, fits$strexp$beta,
                                             k_max = max(ks) * 5L),
           Model = "Stretched Exponential"),
    tibble(k = ks, Probability = dist_pdf(fits$m_exp, ks), Model = "Exponential")
  ) |>
    mutate(Model = factor(Model, levels = c(
      "Log-normal", "Stretched Exponential", "Power-law", "Exponential"
    )))

  colors <- c(
    "Log-normal"           = "#003366",
    "Stretched Exponential"= "#CC0000",
    "Power-law"            = "#7B1FA2",
    "Exponential"          = "#009900"
  )
  linetypes <- c(
    "Log-normal"           = "solid",
    "Stretched Exponential"= "dashed",
    "Power-law"            = "longdash",
    "Exponential"          = "dotted"
  )

  # Wrap Delta-AIC list into 2 entries per line for readability at large font size
  aic_entries <- sprintf("%s: ΔAIC = %s",
            fits$aic$Distribution,
            formatC(fits$aic$DeltaAIC, format = "f", digits = 1, big.mark = ","))
  aic_groups  <- split(aic_entries, ceiling(seq_along(aic_entries) / 2))
  aic_str <- paste(vapply(aic_groups, paste, character(1), collapse = "   |   "),
                    collapse = "\n")

  ggplot() +
    geom_point(data = emp_df, aes(x = k, y = Probability),
               color = "grey50", alpha = 0.5, size = 1.8) +
    geom_line(data = curves_df,
              aes(x = k, y = Probability, color = Model, linetype = Model),
              linewidth = 1.3) +
    scale_color_manual(values = colors) +
    scale_linetype_manual(values = linetypes) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                   labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                   labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides = "bl") +
    labs(
      title    = paste("Fitted distributions —", network_name),
      subtitle = aic_str,
      x = "Degree (k)", y = "P(K = k)",
      color = "Model", linetype = "Model"
    ) +
    theme_bw(base_family = "Times New Roman") +
    theme(
      plot.title       = element_text(size = 26, face = "bold", hjust = 0.5),
      plot.subtitle    = element_text(size = 20, color = "grey20", hjust = 0.5,
                                       lineheight = 1.4,
                                       margin = margin(t = 6, b = 8)),
      axis.title       = element_text(size = 22, face = "bold"),
      axis.text        = element_text(size = 19),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_line(color = "grey95"),
      legend.position  = "top",
      legend.title     = element_text(size = 22, face = "bold"),
      legend.text      = element_text(size = 20),
      legend.key.width = unit(2.2, "cm"),
      legend.key.size  = unit(1.3, "cm")
    )
}

p_bench <- make_benchmark_panel(fits_A7, deg_A7, "HPV-A7") +
  make_benchmark_panel(fits_A9, deg_A9, "HPV-A9") +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(
    title = "Distribution Benchmark — CC Co-expression Networks",
    theme = theme(
      plot.title = element_text(size = 30, face = "bold", hjust = 0.5,
                                 family = "Times New Roman",
                                 margin = margin(b = 10))
    )
  ) &
  theme(legend.position = "bottom")

ggsave(
  "~/CESC_Network/5_Network_analisis/5_2_Network_topology/5_2_3_Distribution_benchmark.png",
  p_bench, width = 22, height = 12, dpi = 300, bg = "white"
)
message("Saved: Distribution benchmark")

################################################################################
# 6) Plot 2 — Degree distribution PDF, log-log
################################################################################

get_pdf_df <- function(graph, label) {
  deg  <- as.numeric(igraph::degree(graph, mode = "all"))
  dist <- table(deg) / length(deg)
  data.frame(Degree = as.numeric(names(dist)),
             Probability = as.numeric(dist),
             Network = label)
}

df_pdf <- rbind(get_pdf_df(A7_HPV, "HPV_A7"), get_pdf_df(A9_HPV, "HPV_A9"))

p_pdf <- ggplot(df_pdf, aes(Degree, Probability, color = Network, shape = Network)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c(HPV_A7 = "#003366", HPV_A9 = "#CC0000")) +
  scale_shape_manual(values = c(HPV_A7 = 16, HPV_A9 = 17)) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "bl") +
  labs(title = "Degree Distribution — CC Co-expression Networks",
       x = "Degree (k)", y = "Probability p(k)",
       color = "Network", shape = "Network") +
  theme_bw(base_family = "Times New Roman") +
  theme(
    plot.title           = element_text(size = 24, face = "bold", hjust = 0.5),
    axis.title           = element_text(size = 20, face = "bold"),
    axis.text            = element_text(size = 17),
    panel.border         = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major     = element_line(color = "grey90"),
    panel.grid.minor     = element_line(color = "grey95"),
    legend.position      = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background    = element_rect(fill = alpha("white", 0.8), color = "black"),
    legend.key           = element_blank(),
    legend.title         = element_text(size = 19, face = "bold"),
    legend.text          = element_text(size = 17)
  )

ggsave(
  "~/CESC_Network/5_Network_analisis/5_2_Network_topology/5_2_2_Degree_distribution_max.png",
  p_pdf, width = 12, height = 9, dpi = 600, bg = "white"
)
message("Saved: PDF degree distribution")

################################################################################
# 7) Plot 3 — CCDF of degree distribution (publication-ready)
#    Empirical survival function: P(K >= k)
#    Log-log axes, discrete step representation, clean annotation box
################################################################################

library(grid)

# ------------------------------------------------------------------------------
# Helper: build empirical CCDF from positive degrees only
# ------------------------------------------------------------------------------
build_ccdf_df <- function(graph, label) {
  deg <- as.integer(igraph::degree(graph, mode = "all"))
  deg <- deg[deg > 0]  # log-scale and log-normal fit require positive values

  if (length(deg) == 0) {
    stop("No positive degrees found for: ", label)
  }

  tab <- table(deg)
  k   <- as.integer(names(tab))
  ord <- order(k)

  k      <- k[ord]
  counts <- as.numeric(tab[ord])
  ccdf   <- rev(cumsum(rev(counts))) / sum(counts)

  tibble(
    Degree = k,
    CCDF   = ccdf,
    Network = factor(label, levels = c("HPV_A7", "HPV_A9"))
  )
}

# ------------------------------------------------------------------------------
# Data frame
# ------------------------------------------------------------------------------
ccdf_df <- bind_rows(
  build_ccdf_df(A7_HPV, "HPV_A7"),
  build_ccdf_df(A9_HPV, "HPV_A9")
)

# ------------------------------------------------------------------------------
# Dynamic model labels from AIC tables
# ------------------------------------------------------------------------------
best_A7 <- fits_A7$aic$Distribution[1]
best_A9 <- fits_A9$aic$Distribution[1]

ann_box <- sprintf(
  paste0(
    "KS test: D = %.3f, p < 2.2×10^-16\n",
    "Best fit: %s (A7), %s (A9)\n",
    "Log-normal fit\n",
    "HPV-A7: μ = %.3f, σ = %.3f\n",
    "HPV-A9: μ = %.3f, σ = %.3f"
  ),
  unname(ks_result$statistic),
  best_A7, best_A9,
  fits_A7$mu, fits_A7$sigma,
  fits_A9$mu, fits_A9$sigma
)
# ------------------------------------------------------------------------------
# Annotation placement
# ------------------------------------------------------------------------------
ccdf_label_x <- min(ccdf_df$Degree) * 1.25
ccdf_label_y <- min(ccdf_df$CCDF) * 1.8

# ------------------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------------------
p_ccdf <- ggplot(ccdf_df, aes(x = Degree, y = CCDF,
                              color = Network, linetype = Network)) +
  geom_step(linewidth = 1.05, alpha = 0.95, direction = "vh") +
  geom_point(size = 1.4, alpha = 0.60) +

  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x)),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x)),
    expand = expansion(mult = c(0.05, 0.02))
  ) +
  annotation_logticks(sides = "bl") +
  scale_color_manual(values = c(HPV_A7 = "#003366", HPV_A9 = "#CC0000")) +
  scale_linetype_manual(values = c(HPV_A7 = "solid", HPV_A9 = "dashed")) +

  annotate(
    "label",
    x = ccdf_label_x, y = ccdf_label_y,
    label = ann_box,
    hjust = 0, vjust = 0,
    family = "Times New Roman",
    fontface = "bold",
    size = 6,
    label.size = 0.35,
    label.padding = grid::unit(0.45, "lines"),
    fill = "white",
    color = "grey10"
  ) +

  labs(
    title   = "CCDF of Degree Distribution",
    x       = "Degree (k)",
    y       = expression(P(K >= k)),
    color   = "Network",
    linetype = "Network"
  ) +
  theme_bw(base_family = "Times New Roman") +
  theme(
    plot.title       = element_text(size = 24, face = "bold", hjust = 0.5),
    axis.title       = element_text(size = 20, face = "bold"),
    axis.text        = element_text(size = 17, color = "black"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_line(color = "grey95", linewidth = 0.2),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1),

    legend.position      = c(0.96, 0.96),
    legend.justification = c("right", "top"),
    legend.background    = element_rect(fill = alpha("white", 0.9),
                                        color = "grey30"),
    legend.key           = element_blank(),
    legend.key.width     = unit(1.5, "cm"),
    legend.key.height    = unit(0.8, "cm"),
    legend.title         = element_text(face = "bold", size = 19),
    legend.text          = element_text(size = 17),
    legend.spacing.y     = unit(0.3, "cm")
  )

ggsave(
  "~/CESC_Network/5_Network_analisis/5_2_Network_topology/5_1_2_Degree_CCDF_distribution_max.png",
  p_ccdf,
  width = 12, height = 8, dpi = 600, bg = "white"
)

################################################################################
# 8) Plot 4 — Log-normal fit validation
#    QQ plot against theoretical log-normal quantiles
################################################################################

build_lognorm_qq <- function(deg_vec, label, lnorm_pars) {

  d <- sort(as.numeric(deg_vec[deg_vec > 0]))

  n <- length(d)
  p <- ppoints(n)

  theoretical <- qlnorm(
    p,
    meanlog = lnorm_pars[1],
    sdlog   = lnorm_pars[2]
  )

  tibble(
    Theoretical = theoretical,
    Observed    = d,
    Network     = label
  )
}

qq_df <- bind_rows(
  build_lognorm_qq(deg_A7, "HPV_A7", fits_A7$lnorm_pars),
  build_lognorm_qq(deg_A9, "HPV_A9", fits_A9$lnorm_pars)
)

p_lognorm <- ggplot(
  qq_df,
  aes(Theoretical, Observed, color = Network)
) +

  geom_point(size = 1.7, alpha = 0.65) +

  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.9,
    color = "black"
  ) +

  scale_color_manual(
    values = c(
      HPV_A7 = "#003366",
      HPV_A9 = "#CC0000"
    )
  ) +

  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +

  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +

  annotation_logticks(sides = "bl") +

  labs(
    title = "Log-normal Fit Validation",
    subtitle = "Observed vs theoretical quantiles",
    x = "Theoretical quantiles",
    y = "Observed degree quantiles",
    color = "Network"
  ) +

  theme_bw(base_family = "Times New Roman") +
  theme(
    plot.title    = element_text(size = 24, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 18, hjust = 0.5, color = "grey30"),
    axis.title    = element_text(size = 20, face = "bold"),
    axis.text     = element_text(size = 17),
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey95"),
    legend.position      = c(0.97, 0.03),
    legend.justification = c("right", "bottom"),
    legend.background    = element_rect(fill = alpha("white", 0.9), color = "grey40"),
    legend.key           = element_blank(),
    legend.title         = element_text(size = 19, face = "bold"),
    legend.text          = element_text(size = 17)
  )

ggsave(
  "~/CESC_Network/5_Network_analisis/5_2_Network_topology/5_2_4_Lognormal_QQ_fit.png",
  p_lognorm,
  width = 12,
  height = 9,
  dpi = 600,
  bg = "white"
)

message("Saved: Log-normal QQ fit")
