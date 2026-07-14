
################################################################################
# Degree Distribution Analysis for Cervical Cancer Co-expression Networks
# Author: Joel Ruiz Hernandez
#
# This script covers the DESCRIPTIVE degree-distribution outputs and the
# two-network comparison:
#   - Two-sample KS + Wilcoxon comparing the A7 vs A9 degree sequences.
#   - Degree-distribution PDF (log-log).
#   - Complementary CDF / survival (log-log).
#   - Log-normal QQ validation plot.
#
# Formal model selection (power-law / log-normal / stretched-exp / cutoff /
# exponential) following the Clauset, Shalizi & Newman (2009) workflow
# (MLE + KS-bootstrap goodness of fit + likelihood-ratio tests) lives in a
# separate script: 5_2_3_Distribution_fitting_Clauset.R.
################################################################################

suppressPackageStartupMessages({
  library(igraph)
  library(ggplot2)
  library(scales)
  library(dplyr)
  library(tibble)
  library(poweRlaw)
})

################################################################################
# 1) Load graphs
################################################################################

A7_HPV <- readRDS("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_2_4_Subnetwork_A7_max.rds")

# A9 is already MI-sorted and edge-matched to A7 upstream
# (5_1_Build_matched_networks.R). Load its graph directly — no old/individual
# edgelist is used; both graphs come solely from the new matched networks.
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
# 4) Log-normal parameters (discrete MLE, xmin = 1) for the QQ validation plot
#    Full model fitting / benchmark moved to 5_2_3_Distribution_fitting_Clauset.R
################################################################################

lnorm_pars_of <- function(deg_vec) {
  d <- as.integer(deg_vec[deg_vec > 0])
  m <- dislnorm$new(d); m$setXmin(1); m$setPars(estimate_pars(m)$pars)
  m$pars
}
lnorm_pars_A7 <- lnorm_pars_of(deg_A7)
lnorm_pars_A9 <- lnorm_pars_of(deg_A9)

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
    plot.title           = element_text(size = 30, face = "bold", hjust = 0.5),
    axis.title           = element_text(size = 27, face = "bold"),
    axis.text            = element_text(size = 23),
    panel.border         = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major     = element_line(color = "grey90"),
    panel.grid.minor     = element_line(color = "grey95"),
    legend.position      = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background    = element_rect(fill = alpha("white", 0.8), color = "black"),
    legend.key           = element_blank(),
    legend.title         = element_text(size = 25, face = "bold"),
    legend.text          = element_text(size = 23)
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

ann_box <- sprintf(
  "KS test: D = %.3f, p < 2.2×10^-16",
  unname(ks_result$statistic)
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
    size = 8.5,
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
    plot.title       = element_text(size = 30, face = "bold", hjust = 0.5),
    axis.title       = element_text(size = 27, face = "bold"),
    axis.text        = element_text(size = 23, color = "black"),
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
    legend.title         = element_text(face = "bold", size = 25),
    legend.text          = element_text(size = 23),
    legend.spacing.y     = unit(0.3, "cm")
  )

ggsave(
  "~/CESC_Network/5_Network_analisis/5_2_Network_topology/5_1_2_Degree_CCDF_distribution_max.png",
  p_ccdf,
  width = 8, height = 12, dpi = 600, bg = "white"
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
  build_lognorm_qq(deg_A7, "HPV_A7", lnorm_pars_A7),
  build_lognorm_qq(deg_A9, "HPV_A9", lnorm_pars_A9)
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
    plot.title    = element_text(size = 30, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 24, hjust = 0.5, color = "grey30"),
    axis.title    = element_text(size = 27, face = "bold"),
    axis.text     = element_text(size = 23),
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey95"),
    legend.position      = c(0.97, 0.03),
    legend.justification = c("right", "bottom"),
    legend.background    = element_rect(fill = alpha("white", 0.9), color = "grey40"),
    legend.key           = element_blank(),
    legend.title         = element_text(size = 25, face = "bold"),
    legend.text          = element_text(size = 23)
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
