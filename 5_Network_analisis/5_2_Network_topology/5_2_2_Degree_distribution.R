
################################################################################
# Degree Distribution Analysis for Cervical Cancer Co-expression Networks
# Author: Joel Ruiz Hernandez
################################################################################

suppressPackageStartupMessages({
  library(igraph)
  library(ggplot2)
  library(scales)
  library(dplyr)
  library(tibble)
  library(fitdistrplus)
  library(patchwork)
})

################################################################################
# 1) Load graphs
################################################################################

A7_HPV <- readRDS("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_2_4_Subnetwork_A7_max.rds")
A9_HPV <- readRDS("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_2_4_Subnetwork_A9_max.rds")

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
# 4) Distribution fitting — MLE + AIC
################################################################################

fit_network <- function(deg_vec, network_name) {
  message("\n  Fitting: ", network_name)
  d <- as.numeric(deg_vec[deg_vec > 0])
  
  fit_ln <- fitdist(d, "lnorm",   method = "mle")
  fit_wb <- fitdist(d, "weibull", method = "mle")
  fit_gm <- fitdist(d, "gamma",   method = "mle")
  fit_ex <- fitdist(d, "exp",     method = "mle")
  
  aic_tbl <- tibble(
    Distribution = c("Log-normal", "Weibull", "Gamma", "Exponential"),
    AIC          = c(fit_ln$aic, fit_wb$aic, fit_gm$aic, fit_ex$aic),
    BIC          = c(fit_ln$bic, fit_wb$bic, fit_gm$bic, fit_ex$bic)
  ) |>
    arrange(AIC) |>
    mutate(DeltaAIC = round(AIC - min(AIC), 1))
  
  print(aic_tbl)
  
  invisible(list(
    lnorm   = fit_ln, weibull = fit_wb,
    gamma   = fit_gm, exp     = fit_ex,
    aic     = aic_tbl,
    deg_num = d,
    mu      = round(fit_ln$estimate["meanlog"], 3),
    sigma   = round(fit_ln$estimate["sdlog"],   3)
  ))
}

fits_A7 <- fit_network(deg_A7, "HPV_A7")
fits_A9 <- fit_network(deg_A9, "HPV_A9")

################################################################################
# 5) Plot 1 — Distribution benchmark
#    Empirical histogram + four theoretical PDFs per network, side by side.
################################################################################

make_benchmark_panel <- function(fits, network_name) {
  
  d     <- fits$deg_num
  x_max <- quantile(d, 0.99)
  xs    <- seq(1, x_max, length.out = 600)
  
  lines_df <- tibble(
    x = rep(xs, 4),
    Density = c(
      dlnorm(xs,   fits$lnorm$estimate["meanlog"],  fits$lnorm$estimate["sdlog"]),
      dweibull(xs, fits$weibull$estimate["shape"],  fits$weibull$estimate["scale"]),
      dgamma(xs,   fits$gamma$estimate["shape"],    fits$gamma$estimate["rate"]),
      dexp(xs,     fits$exp$estimate["rate"])
    ),
    Model = rep(c("Log-normal", "Weibull", "Gamma", "Exponential"), each = length(xs))
  ) |> mutate(Model = factor(Model, levels = c("Log-normal", "Gamma", "Weibull", "Exponential")))
  
  colors    <- c("Log-normal" = "#003366", "Weibull" = "#CC0000",
                 "Gamma" = "#E69F00",      "Exponential" = "#009900")
  linetypes <- c("Log-normal" = "solid",   "Weibull" = "dashed",
                 "Gamma" = "dotdash",      "Exponential" = "dotted")
  
  # Build subtitle from AIC table (ordered as in aic_tbl)
  # Build subtitle from AIC table
  aic_str <- paste(
    sprintf("%s: ΔAIC = %.1f", fits$aic$Distribution, fits$aic$DeltaAIC),
    collapse = "  |  "
  )
  
  ggplot() +
    geom_histogram(aes(x = d, y = after_stat(density)),
                   binwidth = 5, fill = "grey82", color = "grey60",
                   alpha = 0.7, boundary = 0) +
    geom_line(data = lines_df,
              aes(x, Density, color = Model, linetype = Model),
              linewidth = 1.0) +
    scale_color_manual(values = colors) +
    scale_linetype_manual(values = linetypes) +
    coord_cartesian(xlim = c(0, x_max)) +
    labs(
      title    = paste("Fitted distributions —", network_name),
      subtitle = aic_str,
      x = "Degree (k)", y = "Density",
      color = "Model", linetype = "Model"
    ) +
    theme_bw(base_family = "Times New Roman") +
    theme(
      plot.title       = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.subtitle    = element_text(size = 7.5, color = "grey30", hjust = 0.5),
      axis.title       = element_text(size = 11, face = "bold"),
      axis.text        = element_text(size = 9),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_line(color = "grey95"),
      legend.position  = "top",
      legend.title     = element_text(size = 10, face = "bold"),
      legend.text      = element_text(size = 9),
      legend.key.width = unit(1.5, "cm")
    )
}

p_bench <- make_benchmark_panel(fits_A7, "HPV-A7") +
  make_benchmark_panel(fits_A9, "HPV-A9") +
  plot_layout(ncol = 2) +
  plot_annotation(
    title   = "Distribution Benchmark \u2014 CC Co-expression Networks",
    caption = "MLE fitting; AIC model selection; grey bars = empirical density",
    theme   = theme(
      plot.title   = element_text(size = 14, face = "bold", hjust = 0.5,
                                  family = "Times New Roman"),
      plot.caption = element_text(size = 8, color = "grey40",
                                  family = "Times New Roman")
    )
  )

ggsave(
  "~/CESC_Network/5_Network_analisis/5_2_Network_topology/5_2_3_Distribution_benchmark.png",
  p_bench, width = 14, height = 6, dpi = 600
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
  labs(title = "Degree Distribution \u2014 CC Co-expression Networks",
       x = "Degree (k)", y = "Probability p(k)",
       color = "Network", shape = "Network") +
  theme_bw(base_family = "Times New Roman") +
  theme(
    plot.title           = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title           = element_text(size = 12, face = "bold"),
    axis.text            = element_text(size = 10),
    panel.border         = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major     = element_line(color = "grey90"),
    panel.grid.minor     = element_line(color = "grey95"),
    legend.position      = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background    = element_rect(fill = alpha("white", 0.8), color = "black"),
    legend.key           = element_blank(),
    legend.title         = element_text(size = 11, face = "bold"),
    legend.text          = element_text(size = 10)
  )

ggsave(
  "~/CESC_Network/5_Network_analisis/5_2_Network_topology/5_2_2_Degree_distribution_max.png",
  p_pdf, width = 7, height = 6, dpi = 600
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
    "Log-normal fit\n",
    "HPV-A7: μ = %.3f, σ = %.3f\n",
    "HPV-A9: μ = %.3f, σ = %.3f"
  ),
  unname(ks_result$statistic),
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
################################################################################
# 7) CCDF — legend text larger
################################################################################

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
    size = 7,
    label.size = 0.45,
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
    plot.title       = element_text(size = 25, face = "bold", hjust = 0.5),
    axis.title       = element_text(size = 18, face = "bold"),
    axis.text        = element_text(size = 15, color = "black"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_line(color = "grey95", linewidth = 0.2),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1),
    
    # ── Legend: position unchanged, text much larger ────────────────
    legend.position      = c(0.96, 0.96),
    legend.justification = c("right", "top"),
    legend.background    = element_rect(fill = alpha("white", 0.9),
                                        color = "grey30"),
    legend.key           = element_blank(),
    legend.key.width     = unit(1.4, "cm"),   # wider swatch
    legend.key.height    = unit(0.85, "cm"),  # more vertical room
    legend.title         = element_text(face = "bold", size = 18),  # ++ bigger
    legend.text          = element_text(size = 16),                  # ++ bigger
    legend.spacing.y     = unit(0.3, "cm")
  )

p_ccdf


ggsave(
  "~/CESC_Network/5_Network_analisis/5_2_Network_topology/5_1_2_Degree_CCDF_distribution_max.png",
  p_ccdf,
  width = 12, height = 8, dpi = 600, bg = "white"
)


################################################################################
# 8) Plot 4 — Log-normal fit validation
#    QQ plot against theoretical log-normal quantiles
################################################################################

build_lognorm_qq <- function(deg_vec, label, fit_obj) {
  
  d <- sort(as.numeric(deg_vec[deg_vec > 0]))
  
  n <- length(d)
  p <- ppoints(n)
  
  theoretical <- qlnorm(
    p,
    meanlog = fit_obj$estimate["meanlog"],
    sdlog   = fit_obj$estimate["sdlog"]
  )
  
  tibble(
    Theoretical = theoretical,
    Observed    = d,
    Network     = label
  )
}

qq_df <- bind_rows(
  build_lognorm_qq(deg_A7, "HPV_A7", fits_A7$lnorm),
  build_lognorm_qq(deg_A9, "HPV_A9", fits_A9$lnorm)
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
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      size = 10,
      hjust = 0.5,
      color = "grey30"
    ),
    
    axis.title = element_text(
      size = 12,
      face = "bold"
    ),
    
    axis.text = element_text(size = 10),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1
    ),
    
    panel.grid.major = element_line(
      color = "grey90"
    ),
    
    panel.grid.minor = element_line(
      color = "grey95"
    ),
    
    legend.position = c(0.97, 0.03),
    legend.justification = c("right", "bottom"),
    
    legend.background = element_rect(
      fill = alpha("white", 0.9),
      color = "grey40"
    ),
    
    legend.key = element_blank(),
    
    legend.title = element_text(
      face = "bold"
    )
  )

p_lognorm

ggsave(
  "~/CESC_Network/5_Network_analisis/5_2_Network_topology/5_2_4_Lognormal_QQ_fit.png",
  p_lognorm,
  width = 7,
  height = 6,
  dpi = 600,
  bg = "white"
)

message("Saved: Log-normal QQ fit")



