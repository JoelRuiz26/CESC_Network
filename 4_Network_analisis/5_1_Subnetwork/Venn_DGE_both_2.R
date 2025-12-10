### Libraries ###
pacman::p_load("igraph", "tidyverse", "ggraph", "tidygraph", "svglite", "ggvenn")

DE_A7 <- readr::read_tsv(
  "~/CESC_Network/4_Network_analisis/4_3_limma_A7_vs_SNT_allGenes_annot_PR_DEG.tsv",
  col_types = cols()
)

DE_A9 <- readr::read_tsv(
  "~/CESC_Network/4_Network_analisis/4_3_limma_A9_vs_SNT_allGenes_annot_PR_DEG.tsv",
  col_types = cols()
)

# ===================== GENES CON MISMA DIRECCIÓN + adj.P.Val < 0.01 ===================== #
# Preparamos tablas de DE con columnas renombradas

DE_A7_sig <- DE_A7 %>%
  filter(adj.P.Val < 0.01, abs(logFC) > 1) %>%
  transmute(
    Gene,
    logFC_A7 = logFC,
    adj_A7   = adj.P.Val
  )

DE_A9_sig <- DE_A9 %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 1) %>%
  transmute(
    Gene,
    logFC_A9 = logFC,
    adj_A9   = adj.P.Val
  )


DE_both <- inner_join(DE_A7_sig, DE_A9_sig, by = "Gene")


# Genes con dirección opuesta (uno up y otro down)
genes_opposite_dir <- DE_both %>%
  filter(logFC_A7 * logFC_A9 < 0) %>%
  pull(Gene)

# ----------------------------
# 3) Definir conjuntos para el Venn
#    Quitamos los genes con direcciones opuestas
#    => La intersección serán solo genes con misma dirección
# ----------------------------
genes_A7_clean <- setdiff(DE_A7_sig$Gene, genes_opposite_dir)
genes_A9_clean <- setdiff(DE_A9_sig$Gene, genes_opposite_dir)

list_DGE <- list(
  HPV_A7_vs_SNT = genes_A7_clean,
  HPV_A9_vs_SNT = genes_A9_clean
)

# ----------------------------
# 4) Venn con mismo esquema de colores pero más "bonito"
# ----------------------------
Venn_DGE <- ggvenn(
  list_DGE,
  fill_color      = c("#003366", "#FF0000"),  # mismos colores que tu ejemplo
  fill_alpha      = 0.5,
  stroke_size     = 0.5,
  set_name_size   = 7,
  text_size       = 9,
  show_percentage = TRUE,     # % respecto a la unión
  label_sep       = "\n",     # cuenta en 1a línea, % en 2a
  digits          = 1
) +
  labs(
    title    = "Shared DGE HPV_A7 & HPV_A9"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 20),
    plot.subtitle = element_text(hjust = 0.5, size = 10, lineheight = 1.1),
    plot.caption  = element_text(hjust = 0.5, size = 9),
    legend.position = "none"
  )


Venn_DGE <- Venn_DGE +
  theme(
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    axis.title      = element_blank(),
    panel.grid      = element_blank(),   # quita todas las líneas “cuadros”
    panel.border    = element_blank(),   # quita el borde del panel
    panel.background= element_blank()    # fondo limpio sin cuadro
  )


# ----------------------------
# 5) Guardar figura
# ----------------------------
output_file <- "~/CESC_Network/4_Network_analisis/4_4_Venn_DGE_A7_A9.png"

ggsave(
  filename = output_file,
  plot     = Venn_DGE,
  width    = 8,
  height   = 7,
  dpi      = 600
)

output_file
