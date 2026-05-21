############################################################
# SCRIPT: 4_2_DGE_plots_A7_A9.R
#
# Purpose:
#   Visualización de resultados de expresión diferencial:
#     1) Volcano plot  — A7 vs A9
#     2) Venn 2 círculos — DEGs A7-vs-SNT y A9-vs-SNT (todos los compartidos,
#        sin filtrar por dirección)
#     3) Venn 3 círculos — DEGs A7-vs-SNT, A9-vs-SNT y A7-vs-A9
#     4) Scatter plot ordenado por logFC con correlación de Spearman
#
# Input:
#   - 3_1_limma_A7_vs_SNT_allGenes.tsv
#   - 3_1_limma_A9_vs_SNT_allGenes.tsv
#   - 3_1_limma_A7_vs_A9_allGenes.tsv
#
# Color palette (homogenizada):
#   Azul  = #003366  (A7 / significativo negativo)
#   Rojo  = #CC0000  (A9 / significativo positivo)
#   Gris  = #AAAAAA  (no significativo)
############################################################

setwd("~/CESC_Network/3_DGE/")

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(ggvenn)
  library(ggrepel)
  library(patchwork)
  library(scales)
})

# ── Paleta global ──────────────────────────────────────────
COL_A7   <- "#003366"   # azul oscuro
COL_A9   <- "#CC0000"   # rojo
COL_GREY <- "#AAAAAA"   # no-significativo
COL_BOTH <- "#800080"   # intersección Venn (opcional; se usa fill_alpha)

# ── Umbrales ──────────────────────────────────────────────
PADJ_THR <- 0.01
LFC_THR  <- 1

############################################################
## 1) Cargar datos
############################################################

DE_A7_full <- readr::read_tsv("3_1_limma_A7_vs_SNT_allGenes.tsv", col_types = cols())
DE_A9_full <- readr::read_tsv("3_1_limma_A9_vs_SNT_allGenes.tsv", col_types = cols())
DE_A7A9    <- readr::read_tsv("3_1_limma_A7_vs_A9_allGenes.tsv",  col_types = cols())

############################################################
## 2) Conjuntos de DEGs significativos
############################################################

# A7 vs SNT
sig_A7 <- DE_A7_full %>%
  filter(adj.P.Val <= PADJ_THR, abs(logFC) >= LFC_THR) %>%
  pull(Gene)

# A9 vs SNT
sig_A9 <- DE_A9_full %>%
  filter(adj.P.Val <= PADJ_THR, abs(logFC) >= LFC_THR) %>%
  pull(Gene)

# A7 vs A9
sig_A7A9 <- DE_A7A9 %>%
  filter(adj.P.Val <= PADJ_THR, abs(logFC) >= LFC_THR) %>%
  pull(Gene)

############################################################
## FIGURA 1 — Volcano plot: A7 vs A9
############################################################

# Añadir columna de significancia
volcano_df <- DE_A7A9 %>%
  mutate(
    significant = case_when(
      adj.P.Val <= PADJ_THR & logFC >=  LFC_THR ~ "Up in A7",
      adj.P.Val <= PADJ_THR & logFC <= -LFC_THR ~ "Up in A9",
      TRUE                                        ~ "NS"
    ),
    neg_log10_p = -log10(P.Value)
  )

# Top genes para etiquetar (top 15 por |logFC| entre los sig)
top_labels <- volcano_df %>%
  filter(significant != "NS") %>%
  slice_max(order_by = abs(logFC), n = 20)

# Conteos para el subtítulo
n_up_A7 <- sum(volcano_df$significant == "Up in A7")
n_up_A9 <- sum(volcano_df$significant == "Up in A9")

p_volcano <- ggplot(volcano_df,
                    aes(x = logFC, y = neg_log10_p, colour = significant)) +
  geom_point(alpha = 0.55, size = 1.4) +
  scale_colour_manual(
    values = c("Up in A7" = COL_A7, "Up in A9" = COL_A9, "NS" = COL_GREY),
    labels = c(
      "Up in A7" = paste0("Up in A7 (n = ", n_up_A7, ")"),
      "Up in A9" = paste0("Up in A9 (n = ", n_up_A9, ")"),
      "NS"       = "Not significant"
    )
  ) +
  geom_hline(yintercept = -log10(PADJ_THR), linetype = "dashed",
             colour = "grey40", linewidth = 0.5) +
  geom_vline(xintercept = c(-LFC_THR, LFC_THR), linetype = "dashed",
             colour = "grey40", linewidth = 0.5) +
  ggrepel::geom_text_repel(
    data          = top_labels,
    aes(label     = Gene),
    size          = 2.8,
    max.overlaps  = 20,
    segment.size  = 0.3,
    segment.alpha = 0.6,
    show.legend   = FALSE
  ) +
  labs(
    title    = "Volcano Plot — HPV A7 vs A9 tumors",
    subtitle = paste0("adj.P.Val ≤ ", PADJ_THR, "  |  |logFC| ≥ ", LFC_THR,
                      "  |  Positive logFC = higher in A7"),
    x        = expression(log[2]~Fold~Change~(A7/A9)),
    y        = expression(-log[10]~(p~value)),
    colour   = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 15),
    plot.subtitle   = element_text(hjust = 0.5, size = 10, colour = "grey35"),
    legend.position = c(0.82, 0.88),
    legend.background = element_rect(fill = alpha("white", 0.7), colour = NA),
    legend.text     = element_text(size = 10)
  )

ggsave(
  filename = "~/CESC_Network/3_DGE/3_1_1_Volcano_A7_vs_A9.png",
  plot     = p_volcano,
  width    = 8, height = 6.5, dpi = 600
)
message("✔ Volcano plot guardado.")

############################################################
## FIGURA 2 — Venn 2 círculos: todos los DEGs compartidos
##            (sin filtrar por dirección)
############################################################

list_2venn <- list(
  HPV_A7_vs_SNT = sig_A7,
  HPV_A9_vs_SNT = sig_A9
)

p_venn2 <- ggvenn(
  list_2venn,
  fill_color      = c(COL_A7, COL_A9),
  fill_alpha      = 0.50,
  stroke_size     = 0.5,
  set_name_size   = 7,
  text_size       = 9,
  show_percentage = TRUE,
  digits          = 1
) +
  labs(title = "Shared DEGs — HPV A7 & A9 vs Normal") +
  theme_void(base_size = 14) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 18,
                                   margin = margin(b = 8)),
    legend.position = "none"
  )

ggsave(
  filename = "~/CESC_Network/3_DGE/3_1_2_Venn2_DEGs_A7_A9.png",
  plot     = p_venn2,
  width    = 8, height = 7, dpi = 600
)
message("✔ Venn 2 círculos guardado.")

############################################################
## FIGURA 3 — Venn 3 círculos: A7-vs-SNT, A9-vs-SNT, A7-vs-A9
############################################################

list_3venn <- list(
  HPV_A7_vs_SNT = sig_A7,
  HPV_A9_vs_SNT = sig_A9,
  A7_vs_A9      = sig_A7A9
)

# ggvenn acepta hasta 4 conjuntos; para 3 usa un color intermedio (morado)
p_venn3 <- ggvenn(
  list_3venn,
  fill_color      = c(COL_A7, COL_A9, "#7B0D8E"),   # azul, rojo, morado
  fill_alpha      = 0.45,
  stroke_size     = 0.5,
  set_name_size   = 5.5,
  text_size       = 7,
  show_percentage = TRUE,
  digits          = 1
) +
  labs(
    title    = "DEGs — A7-vs-SNT · A9-vs-SNT · A7-vs-A9",
    subtitle = paste0("adj.P.Val ≤ ", PADJ_THR, "  |  |logFC| ≥ ", LFC_THR)
  ) +
  theme_void(base_size = 13) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 16,
                                 margin = margin(b = 4)),
    plot.subtitle = element_text(hjust = 0.5, size = 10, colour = "grey40",
                                 margin = margin(b = 8)),
    legend.position = "none"
  )

ggsave(
  filename = "~/CESC_Network/3_DGE/3_1_3_Venn3_DEGs_A7_A9_contrast.png",
  plot     = p_venn3,
  width    = 9, height = 8, dpi = 600
)
message("✔ Venn 3 círculos guardado.")

############################################################
## FIGURA 4 — Scatter plot ordenado por logFC A7-vs-SNT
##            con correlación de Spearman
############################################################

# Genes DE significativos en AMBOS contrastes vs SNT
DE_A7_sig_df <- DE_A7_full %>%
  filter(adj.P.Val <= PADJ_THR, abs(logFC) >= LFC_THR) %>%
  transmute(Gene, logFC_A7 = logFC)

DE_A9_sig_df <- DE_A9_full %>%
  filter(adj.P.Val <= PADJ_THR, abs(logFC) >= LFC_THR) %>%
  transmute(Gene, logFC_A9 = logFC)

scatter_df <- inner_join(DE_A7_sig_df, DE_A9_sig_df, by = "Gene") %>%
  arrange(logFC_A7) %>%
  mutate(rank = row_number())

n_common <- nrow(scatter_df)

# Correlación de Spearman
sp <- cor.test(scatter_df$logFC_A7, scatter_df$logFC_A9, method = "spearman")
rho_label <- paste0(
  "Spearman ρ = ", round(sp$estimate, 3),
  "\np < 2.2×10⁻¹⁶",
  "\nn = ", n_common, " genes"
)

# Colorear según dirección concordante vs discordante
scatter_df <- scatter_df %>%
  mutate(
    direction = case_when(
      logFC_A7 > 0 & logFC_A9 > 0 ~ "Both UP",
      logFC_A7 < 0 & logFC_A9 < 0 ~ "Both DOWN",
      TRUE                          ~ "Discordant"
    )
  )

p_scatter <- ggplot(scatter_df, aes(x = logFC_A7, y = logFC_A9, colour = direction)) +
  geom_point(alpha = 0.65, size = 1.6) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey25",
              linewidth = 0.7, linetype = "dashed") +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.4) +
  scale_colour_manual(
    values = c("Both UP" = COL_A9, "Both DOWN" = COL_A7, "Discordant" = "#888888"),
    name   = "Direction"
  ) +
  annotate(
    "label",
    x     = min(scatter_df$logFC_A7) * 0.85,
    y     = max(scatter_df$logFC_A9) * 0.92,
    label = rho_label,
    hjust = 0, vjust = 1,
    size  = 3.8,
    fill  = "white", colour = "grey20",
    label.size = 0.3
  ) +
  labs(
    title    = "logFC concordance — HPV A7 vs A9 (DEGs shared vs Normal)",
    subtitle = paste0("Genes DE in both A7-vs-SNT and A9-vs-SNT  |  ",
                      "adj.P.Val ≤ ", PADJ_THR, ", |logFC| ≥ ", LFC_THR,
                      "  |  Ordered by logFC A7-vs-SNT"),
    x = expression(log[2]~FC~~(A7~vs~Normal)),
    y = expression(log[2]~FC~~(A9~vs~Normal))
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 9.5, colour = "grey35"),
    legend.position = c(0.88, 0.15),
    legend.background = element_rect(fill = alpha("white", 0.75), colour = NA),
    legend.title  = element_text(face = "bold", size = 10),
    legend.text   = element_text(size = 9)
  )

ggsave(
  filename = "~/CESC_Network/3_DGE/3_1_4_Scatter_logFC_A7_A9.png",
  plot     = p_scatter,
  width    = 8, height = 7, dpi = 600
)

message("✔ Scatter plot guardado.")

############################################################
## Resumen en consola
############################################################

message("\n── Resumen de conjuntos ──────────────────────────────────")
message("DEGs A7-vs-SNT   : ", length(sig_A7))
message("DEGs A9-vs-SNT   : ", length(sig_A9))
message("DEGs A7-vs-A9    : ", length(sig_A7A9))
message("Compartidos A7∩A9 (vs SNT): ", length(intersect(sig_A7, sig_A9)))
message("Scatter (ambos sig): ", n_common)
message("Spearman ρ       : ", round(sp$estimate, 4))
message("──────────────────────────────────────────────────────────\n")
message("Todas las figuras guardadas en ~/CESC_Network/3_DGE/")