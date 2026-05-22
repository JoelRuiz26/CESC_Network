############################################################
## 6_0_hub_threshold.R
##
## Visualización de la distribución de conectividad (degree)
## de la red para justificar el umbral de hub genes.
##
## Plots:
##   A) Distribución de degree (histograma, escala log)
##      con líneas en percentiles clave
##   B) Curva de degree ordenado (rank plot) — permite ver
##      el "codo" donde emergen los hubs
##   C) Fracción de conexiones acumuladas — qué % del total
##      de enlaces captura cada umbral
##   D) Tabla de umbrales: n genes y degree mínimo por corte
##
## Output: 6_0_hub_threshold.pdf / .png
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
})

out_dir <- "~/CESC_Network/6_DepMap_cacu"
deg_dir <- "~/CESC_Network/5_Network_analisis/5_4_DEG_integration"

# ── 1) Cargar tablas ─────────────────────────────────────────
DE_A7 <- read_tsv(file.path(deg_dir, "5_4_1_DEG_A7_Degree.tsv"),
                  show_col_types = FALSE) %>%
  mutate(network = "HPV-A7")

DE_A9 <- read_tsv(file.path(deg_dir, "5_4_1_DEG_A9_Degree.tsv"),
                  show_col_types = FALSE) %>%
  mutate(network = "HPV-A9")

df_all <- bind_rows(DE_A7, DE_A9) %>%
  mutate(network = factor(network, levels = c("HPV-A7", "HPV-A9")))

# Umbrales a evaluar
CUTS <- c(0.01, 0.02, 0.05, 0.10, 0.20)
CUT_LABELS <- paste0(CUTS * 100, "%")

NET_COL <- c("HPV-A7" = "#457B9D", "HPV-A9" = "#E63946")

# ── 2) Tabla resumen por umbral ───────────────────────────────
cat("\n── Resumen de umbrales ─────────────────────────────────\n")
for (net in c("HPV-A7", "HPV-A9")) {
  d <- df_all %>% filter(network == net) %>% pull(degree)
  n_total <- length(d)
  cat(sprintf("\n%s  (N = %d genes, degree mediana = %.0f, max = %.0f)\n",
              net, n_total, median(d), max(d)))
  cat(sprintf("  %-8s  %-10s  %-12s\n", "Corte", "N genes", "Degree mín"))
  for (i in seq_along(CUTS)) {
    n_hub  <- ceiling(CUTS[i] * n_total)
    d_min  <- sort(d, decreasing = TRUE)[n_hub]
    cat(sprintf("  %-8s  %-10d  %-12.0f\n", CUT_LABELS[i], n_hub, d_min))
  }
}

# ── 3) Plot A — Histograma de degree (log x) ─────────────────
# Percentil 5% por red para línea de referencia
vlines_df <- df_all %>%
  group_by(network) %>%
  summarise(
    d_p05 = quantile(degree, 0.95),  # top 5% = por encima del percentil 95
    d_p10 = quantile(degree, 0.90),
    .groups = "drop"
  )

pA <- ggplot(df_all, aes(x = degree, fill = network, colour = network)) +
  geom_histogram(bins = 60, alpha = 0.5, position = "identity") +
  geom_vline(data = vlines_df,
             aes(xintercept = d_p05, colour = network),
             linetype = "dashed", linewidth = 0.8) +
  geom_vline(data = vlines_df,
             aes(xintercept = d_p10, colour = network),
             linetype = "dotted", linewidth = 0.6) +
  scale_x_log10() +
  scale_fill_manual(values  = NET_COL) +
  scale_colour_manual(values = NET_COL) +
  annotate("text", x = Inf, y = Inf,
           label = "dashed = top 5%\ndotted = top 10%",
           hjust = 1.05, vjust = 1.3, size = 3, colour = "grey30") +
  labs(
    title  = "A) Distribución de degree",
    x      = "Degree (escala log)",
    y      = "Número de genes",
    fill   = NULL, colour = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = c(0.85, 0.80))

# ── 4) Plot B — Rank plot (degree ordenado) ───────────────────
rank_df <- df_all %>%
  group_by(network) %>%
  arrange(desc(degree), .by_group = TRUE) %>%
  mutate(rank_pct = row_number() / n() * 100) %>%
  ungroup()

# Líneas de percentil para el rank plot
hlines_df <- df_all %>%
  group_by(network) %>%
  summarise(
    across(
      .cols = c(),  # placeholder; compute manually below
      .fns  = list()
    )
  )

pct_lines <- df_all %>%
  group_by(network) %>%
  reframe(
    cut_pct   = CUTS * 100,
    cut_label = CUT_LABELS,
    degree_at = {
      d_sorted <- sort(degree, decreasing = TRUE)
      sapply(CUTS, function(p) d_sorted[ceiling(p * length(d_sorted))])
    }
  )

pB <- ggplot(rank_df, aes(x = rank_pct, y = degree, colour = network)) +
  geom_line(linewidth = 0.7, alpha = 0.9) +
  geom_vline(xintercept = CUTS * 100,
             linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  annotate("text",
           x     = CUTS * 100,
           y     = max(rank_df$degree) * 0.95,
           label = CUT_LABELS,
           angle = 90, hjust = 1, vjust = -0.3,
           size  = 3, colour = "grey40") +
  scale_colour_manual(values = NET_COL) +
  scale_y_log10() +
  labs(
    title  = "B) Curva de degree ordenado (rank plot)",
    x      = "Percentil de genes (% más conectados → izquierda)",
    y      = "Degree (escala log)",
    colour = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = c(0.85, 0.80))

# ── 5) Plot C — Fracción de conexiones acumuladas ─────────────
cum_df <- df_all %>%
  group_by(network) %>%
  arrange(desc(degree), .by_group = TRUE) %>%
  mutate(
    rank_pct  = row_number() / n() * 100,
    cum_links = cumsum(degree) / sum(degree) * 100
  ) %>%
  ungroup()

pC <- ggplot(cum_df, aes(x = rank_pct, y = cum_links, colour = network)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = CUTS * 100,
             linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_hline(yintercept = 50, linetype = "dotted",
             colour = "grey50", linewidth = 0.4) +
  scale_colour_manual(values = NET_COL) +
  labs(
    title  = "C) Fracción acumulada de conexiones",
    x      = "Top X% de genes (por degree)",
    y      = "% de conexiones totales capturadas",
    colour = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = c(0.75, 0.30))

# ── 6) Plot D — Tabla visual de umbrales ──────────────────────
tbl_df <- df_all %>%
  group_by(network) %>%
  summarise(
    n_total = n(),
    data    = list(degree),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    rows = list(
      tibble(
        cut_pct   = CUTS * 100,
        cut_label = CUT_LABELS,
        n_hub     = ceiling(CUTS * n_total),
        deg_min   = {
          d_s <- sort(data, decreasing = TRUE)
          sapply(CUTS, function(p) d_s[ceiling(p * length(d_s))])
        },
        pct_links = {
          d_s <- sort(data, decreasing = TRUE)
          sapply(CUTS, function(p) {
            n_h <- ceiling(p * length(d_s))
            sum(d_s[1:n_h]) / sum(d_s) * 100
          })
        }
      )
    )
  ) %>%
  select(network, rows) %>%
  unnest(rows)

pD <- ggplot(tbl_df,
             aes(x = cut_label, y = n_hub, fill = network, colour = network)) +
  geom_col(position = position_dodge(0.7), width = 0.6, alpha = 0.7) +
  geom_text(aes(label = sprintf("%d\n(≥%d)", n_hub, as.integer(deg_min))),
            position = position_dodge(0.7),
            vjust = -0.3, size = 2.8, colour = "grey20") +
  scale_fill_manual(values   = NET_COL) +
  scale_colour_manual(values = NET_COL) +
  labs(
    title    = "D) N hub genes y degree mínimo por umbral",
    subtitle = "Etiquetas: n genes (degree mínimo del grupo)",
    x        = "Umbral (top X%)",
    y        = "N hub genes",
    fill     = NULL, colour = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = c(0.15, 0.80))

# ── 7) Guardar ────────────────────────────────────────────────
fig <- (pA + pB) / (pC + pD) +
  plot_annotation(
    title = "Distribución de conectividad de la red | Selección de umbral Hub genes",
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13))
  )

ggsave(file.path(out_dir, "6_0_hub_threshold.pdf"),
       plot = fig, width = 14, height = 10, device = cairo_pdf)
ggsave(file.path(out_dir, "6_0_hub_threshold.png"),
       plot = fig, width = 14, height = 10, dpi = 300)

cat("\n6_0_hub_threshold.pdf / .png guardados.\n")
