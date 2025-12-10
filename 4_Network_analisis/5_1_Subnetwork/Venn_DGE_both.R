rm(list = ls())
gc()

library(tidyverse)
library(ggrepel)
library(readr)

out_dir <- "~/CESC_Network/4_Network_analisis/5_3_DeltaCon/"

# -------------------- CARGA CONTRIBUCIONES -------------------- #
contrib_df <- read_tsv(
  file.path(out_dir, "5_3_2_DeltaCon_gene_contributions_A7_vs_A9.tsv"),
  show_col_types = FALSE
)

# -------------------- CARGA DEGs -------------------- #
DE_A7 <- read_tsv(
  "~/CESC_Network/4_Network_analisis/4_3_limma_A7_vs_SNT_allGenes_annot_PR_DEG.tsv",
  col_types = cols()
)

DE_A9 <- read_tsv(
  "~/CESC_Network/4_Network_analisis/4_3_limma_A9_vs_SNT_allGenes_annot_PR_DEG.tsv",
  col_types = cols()
)

DE_A7_sig <- DE_A7 %>%
  filter(adj.P.Val < 0.01, abs(logFC) > 1) %>%
  transmute(Gene, logFC_A7 = logFC)

DE_A9_sig <- DE_A9 %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 1) %>%
  transmute(Gene, logFC_A9 = logFC)

DE_both <- inner_join(DE_A7_sig, DE_A9_sig, by = "Gene")

# -------------------- FILTRAR |Z| >= 3 -------------------- #
highlight_df <- contrib_df %>%
  filter(abs(Z_robust) >= 3)

label_df <- highlight_df %>%
  slice_max(order_by = abs(Z_robust), n = 20)

# -------------------- MARCAR DGE EN highlight_df -------------------- #
highlight_df <- highlight_df %>%
  mutate(
    DE_status = factor(
      if_else(Gene %in% DE_both$Gene, "DE_both", "Not_DE"),
      levels = c("Not_DE", "DE_both")
    )
  )

label_df <- label_df %>%
  mutate(
    DE_status = factor(
      if_else(Gene %in% DE_both$Gene, "DE_both", "Not_DE"),
      levels = c("Not_DE", "DE_both")
    )
  )

# -------------------- PLOT -------------------- #

p_affinity_Z3 <- ggplot(
  highlight_df,
  aes(x = Affinity_A7, y = Affinity_A9,
      color = Higher_affinity,
      shape = DE_status)
) +
  # puntos: borde según HPV_A7/HPV_A9, relleno solo si DE_both
  geom_point(
    aes(fill = if_else(DE_status == "DE_both",
                       Higher_affinity,
                       NA_character_)),
    alpha  = 0.7,
    size   = 3.5,
    stroke = 0.8
  ) +
  # etiquetas top 20
  geom_text_repel(
    data = label_df,
    aes(label = Gene),
    color = "grey30",
    size = 4.3,
    fontface = "bold",
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = 60,
    segment.alpha = 0.7
  ) +
  # línea identidad
  geom_abline(
    slope = 1, intercept = 0,
    linetype = "dashed", color = "gray50"
  ) +
  # colores de borde (HPV_A7 azul, HPV_A9 rojo)
  scale_color_manual(
    values = c("A7" = "#377EB8",   # HPV_A7 azul
               "A9" = "#E41A1C"),  # HPV_A9 rojo
    breaks = c("A7", "A9"),
    labels = c("HPV_A7", "HPV_A9"),
    name   = "Higher-affinity\nnetwork"
  ) +
  # relleno: mismo color que el borde, solo para DE_both
  scale_fill_manual(
    values   = c("A7" = "#377EB8",
                 "A9" = "#E41A1C"),
    na.value = NA
  ) +
  # forma: ambas categorías son círculos, pero se usan para la leyenda
  scale_shape_manual(
    values = c("Not_DE" = 21, "DE_both" = 21),
    name   = "DE status",
    labels = c("Not_DE (hollow)", "DE_both (filled)"
    )
  ) +
  guides(
    fill  = "none",
    shape = guide_legend(
      override.aes = list(
        color = "black",          # en la leyenda, borde negro
        fill  = c(NA, "black"),   # Not_DE hueco, DE_both lleno
        alpha = 1
      )
    )
  ) +
  labs(
    title    = "Gene-level DeltaCon contribution (HPV_A7 vs HPV_A9)",
    subtitle = "Only genes with |robust Z-score of contribution| ≥ 3",
    x = "Affinity in HPV_A7 network",
    y = "Affinity in HPV_A9 network"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title    = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 13),
    legend.text  = element_text(size = 11)
  )

p_affinity_Z3

# -------------------- GUARDAR -------------------- #
ggsave(file.path(out_dir,
                 "5_3_4_DeltaCon_affinity_Z3_with_DElegend_ahora_siChido.pdf"),
       p_affinity_Z3, width = 9, height = 7)

ggsave(file.path(out_dir,
                 "5_3_4_DeltaCon_affinity_Z3_with_DElegend_ahorasichido.png"),
       p_affinity_Z3, width = 9, height = 7, dpi = 600)

# -------------------- CUÁNTOS DE_both HAY EN ESTE CONJUNTO ------------------ #
n_DE_both_in_Z3 <- highlight_df %>%
  semi_join(DE_both, by = "Gene") %>%
  nrow()

n_DE_both_in_Z3
# [1] 551 

prop_DE_both_in_Z3 <- n_DE_both_in_Z3 / nrow(DE_both)
prop_DE_both_in_Z3
