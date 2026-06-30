################################################################################
# Degree Conservation Scatterplot — HPV-A7 vs HPV-A9  (versión mejorada)
# Escala log-log; coherente con distribución log-normal confirmada
# Énfasis visual en correlación baja-moderada
################################################################################

suppressPackageStartupMessages({
  library(igraph)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(dplyr)
})

# ── 1) Cargar grafos ──────────────────────────────────────────────────────────

# A7_HPV <- readRDS("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_2_4_Subnetwork_A7_max.rds")
# A9_HPV <- readRDS("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_2_4_Subnetwork_A9_max.rds")

A7_HPV <- readRDS("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_2_4_Subnetwork_A7_max.rds")
A9_HPV <- readRDS("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_2_4_Subnetwork_A9_max.rds")

# ── 2) Nodos compartidos + grado ─────────────────────────────────────────────

common_nodes <- intersect(V(A7_HPV)$name, V(A9_HPV)$name)

degree_df <- data.frame(
  Gene   = common_nodes,
  deg_A7 = as.numeric(igraph::degree(A7_HPV, v = common_nodes)),
  deg_A9 = as.numeric(igraph::degree(A9_HPV, v = common_nodes)),
  stringsAsFactors = FALSE
) |>
  filter(deg_A7 > 0, deg_A9 > 0)

# ── 3) Color por grupo dominante ──────────────────────────────────────────────

degree_df <- degree_df |>
  mutate(
    delta       = deg_A7 - deg_A9,
    color_group = factor(
      ifelse(delta >= 0, "HPV-A7", "HPV-A9"),
      levels = c("HPV-A7", "HPV-A9")
    )
  )

# ── 4) Spearman ───────────────────────────────────────────────────────────────

cor_test <- cor.test(degree_df$deg_A7, degree_df$deg_A9, method = "spearman")
rho_val  <- round(cor_test$estimate, 3)

rho_label <- dplyr::case_when(
  abs(rho_val) < 0.30 ~ "weak",
  abs(rho_val) < 0.50 ~ "low\u2013moderate",
  abs(rho_val) < 0.70 ~ "moderate",
  TRUE                 ~ "strong"
)

ann_label <- sprintf(
  "Spearman \u03c1 = %.3f  (%s)\np < 2.2 \u00d7 10\u207b\u00b9\u2076\nn = %d",
  rho_val, rho_label, nrow(degree_df)
)

# ── 5) Top 5 hubs por red (solo extremos) ─────────────────────────────────────

top_A7 <- degree_df |>
  filter(color_group == "HPV-A7") |>
  slice_max(deg_A7, n = 5)

top_A9 <- degree_df |>
  filter(color_group == "HPV-A9") |>
  slice_max(deg_A9, n = 5)

# Unión sin duplicados; sin top_discord
label_df <- bind_rows(top_A7, top_A9) |>
  distinct(Gene, .keep_all = TRUE)

# ── 6) Paleta y posición de la anotación ─────────────────────────────────────

pal <- c("HPV-A7" = "#003366",
         "HPV-A9" = "#CC2200")

x_ann <- min(degree_df$deg_A7)^1.05
y_ann <- max(degree_df$deg_A9) * 0.85

# ── 7) Plot ───────────────────────────────────────────────────────────────────

p_scatter <- ggplot(degree_df,
                    aes(x = deg_A7, y = deg_A9, color = color_group)) +
  
  # diagonal de perfecta concordancia
  geom_abline(
    slope     = 1, intercept = 0,
    linetype  = "dashed",
    linewidth = 0.6,
    color     = "grey50"
  ) +
  
  # nube de puntos
  geom_point(alpha = 0.25, size = 1.5,
             key_glyph = draw_key_point) +
  
  # caja de estadístico
  annotate(
    "label",
    x             = x_ann,
    y             = y_ann,
    hjust         = 0, vjust = 1,
    label         = ann_label,
    fontface      = "bold",
    size          = 6.5,
    fill          = alpha("white", 0.92),
    color         = "grey15",
    label.size    = 0.45,
    label.padding = unit(0.5, "lines")
  ) +
  
  # escalas log-log
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x)),
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x)),
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  annotation_logticks(sides = "bl", color = "grey50", linewidth = 0.35) +
  
  scale_color_manual(
    values = pal,
    guide  = guide_legend(override.aes = list(alpha = 0.85, size = 4))
  ) +
  
  labs(
    title    = "Degree Conservation Across HPV Clade Networks",
    subtitle = sprintf(
      "Log-log scale \u2014 diagonal = perfect concordance  |  Spearman \u03c1 = %.3f (%s)",
      rho_val, rho_label
    ),
    x     = "Degree in HPV-A7  (log scale)",
    y     = "Degree in HPV-A9  (log scale)",
    color = NULL
  ) +
  
  theme_bw(base_family = "Times New Roman") +
  theme(
    plot.title       = element_text(size = 24, face = "bold",
                                    hjust = 0.5, margin = margin(b = 4)),
    plot.subtitle    = element_text(size = 18, hjust = 0.5,
                                    color = "grey35", margin = margin(b = 8)),
    axis.title       = element_text(size = 20, face = "bold"),
    axis.text        = element_text(size = 17, color = "black"),
    axis.line        = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white"),
    legend.position      = c(0.87, 0.13),
    legend.justification = c("right", "bottom"),
    legend.background    = element_rect(fill  = alpha("white", 0.88),
                                        color = "grey60", linewidth = 0.4),
    legend.key           = element_blank(),
    legend.key.size      = unit(0.8, "cm"),
    legend.text          = element_text(size = 17),
    legend.spacing.y     = unit(0.3, "cm"),
    plot.margin          = margin(12, 15, 10, 12)
  )

ggsave(
  "~/CESC_Network/5_Network_analisis/5_2_Network_topology/5_2_5_Degree_Scatter.png",
  p_scatter,
  width = 12, height = 9, dpi = 600, bg = "white"
)

message("Saved: Degree conservation scatterplot (log-log)")

