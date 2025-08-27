load("~/CESC_Network/6_OCTAD/6_2_4_OCTAD.RData")
setwd("~/CESC_Network/6_OCTAD/")

library(pheatmap)
library(ggplot2)
library(dplyr)
library(tibble)
library(ggrepel)
library(readr)
library(tidyr)
library(purrr)


#SIMILARITY HEATMAP
# Preprocess
rownames(similitud_A7) <- similitud_A7$cell_id
rownames(similitud_A9) <- similitud_A9$cell_id

combined_similitud <- full_join(
  similitud_A7 %>% dplyr::select(cell_id, medcor) %>% dplyr::rename(clado_A7 = medcor),
  similitud_A9 %>% dplyr::select(cell_id, medcor) %>% dplyr::rename(clado_A9 = medcor),
  by = "cell_id"
)

combined_similitud$max_cor <- apply(combined_similitud[, c("clado_A7", "clado_A9")], 1, max, na.rm = TRUE)
combined_similitud <- combined_similitud[order(-combined_similitud$max_cor), ]
rownames(combined_similitud) <- combined_similitud$cell_id

cor_matrix_combined <- as.matrix(combined_similitud[, c("clado_A7", "clado_A9")])
labels <- rownames(cor_matrix_combined)
selected_cells <- union(lineas_similares_A7, lineas_similares_A9)
label_exp <- ifelse(labels %in% selected_cells,
                    paste0("bold('", labels, "')"),
                    paste0("plain('", labels, "')"))

# PDF
pdf("Similarity_Heatmap.pdf", width=10, height=10)
pheatmap(
  cor_matrix_combined,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "Similarity between cell lines and cervical cancer samples",
  na_col = "grey90",
  labels_row = parse(text = label_exp),
  fontsize_row = 10,
  fontsize_col = 14,
  fontface_col = "bold",
  cellwidth = 65,
  angle_col = 0
)
dev.off()

# PNG
png("Similarity_Heatmap.png", width=1000, height=1000, res=150)
pheatmap(
  cor_matrix_combined,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "Similarity between cell lines and cervical cancer samples",
  na_col = "grey90",
  labels_row = parse(text = label_exp),
  fontsize_row = 10,
  fontsize_col = 14,
  fontface_col = "bold",
  cellwidth = 65,
  angle_col = 0
)
dev.off()




# DRUG PREDICTION SCATTER PLOT
# ==========================================================
# SETUP
# ==========================================================
setwd("~/CESC_Network/6_OCTAD/")

library(ggplot2)
library(dplyr)
library(tibble)
library(ggrepel)
library(readr)
library(tidyr)
library(purrr)

# Asegúrate de tener cargado el RData
# load("~/CESC_Network/6_OCTAD/6_3_OCTAD.RData")

# ==========================================================
# DEFINIR CONJUNTOS LÓGICOS DE CATEGORÍAS
# ==========================================================
UNICOS_A7 <- setdiff(result_A7_0.2$pert_iname, result_A9_0.2$pert_iname)
UNICOS_A9 <- setdiff(result_A9_0.2$pert_iname, result_A7_0.2$pert_iname)
SHARED <- intersect(result_A7_0.2$pert_iname, result_A9_0.2$pert_iname)

# ==========================================================
# COMBINAR DATOS Y CATEGORÍAS
# ==========================================================
scatter_data <- full_join(
  result_A7 %>% dplyr::select(pert_iname, A7 = sRGES),
  result_A9 %>% dplyr::select(pert_iname, A9 = sRGES),
  by = "pert_iname"
) %>%
  dplyr::mutate(
    categoria = case_when(
      pert_iname %in% SHARED ~ "Shared",
      pert_iname %in% UNICOS_A7 ~ "clade_A7",
      pert_iname %in% UNICOS_A9 ~ "clade_A9",
      TRUE ~ "Non-significant"
    )
  )

# ==========================================================
# SELECCIÓN DE ETIQUETAS MÁS NEGATIVAS Y SIN BRD
# ==========================================================
seleccionar_etiquetas <- function(df, categoria, var1, var2) {
  sub <- df %>%
    filter(categoria == !!categoria, !is.na({{var1}}), !is.na({{var2}})) %>%
    mutate(score_sum = {{var1}} + {{var2}}) %>%
    arrange(score_sum)
  
  # Primero los que no son BRD
  no_brd <- sub %>% filter(!grepl("^BRD", pert_iname)) %>% slice_head(n = 5)
  
  # Completar con BRD si faltan
  if (nrow(no_brd) < 5) {
    n_extra <- 5 - nrow(no_brd)
    extra_brd <- sub %>% filter(grepl("^BRD", pert_iname)) %>% slice_head(n = n_extra)
    no_brd <- bind_rows(no_brd, extra_brd)
  }
  return(no_brd)
}

set.seed(123)

etiquetas_azul <- seleccionar_etiquetas(scatter_data, "clade_A7", A7, A9)
etiquetas_rojo <- seleccionar_etiquetas(scatter_data, "clade_A9", A7, A9)
etiquetas_verde <- seleccionar_etiquetas(scatter_data, "Shared", A7, A9)

etiquetas_data <- bind_rows(etiquetas_azul, etiquetas_rojo, etiquetas_verde)

# ==========================================================
# LEYENDA CON CUENTAS
# ==========================================================
label_counts <- scatter_data %>%
  group_by(categoria) %>%
  summarise(n = n()) %>%
  mutate(label = paste0(categoria, " (n=", n, ")"))

color_values <- c(
  "Shared" = "darkgreen",
  "clade_A9" = "red",
  "clade_A7" = "steelblue",
  "Non-significant" = "grey70"
)

# ==========================================================
# RANGOS PARA EJE (0 hasta más negativo)
# ==========================================================
min_x <- min(0, scatter_data$A7, na.rm = TRUE)
min_y <- min(0, scatter_data$A9, na.rm = TRUE)

# ==========================================================
# PLOTEO
# ==========================================================
p_final <- ggplot(scatter_data, aes(x = A7, y = A9, color = categoria)) +
  geom_point(alpha = 0.8, size = 4) +
  geom_hline(yintercept = -0.2, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = -0.2, linetype = "dashed", color = "grey40") +
  geom_text_repel(
    data = etiquetas_data,
    aes(label = pert_iname),
    size = 6,
    fontface = "bold",
    color = "black",
    max.overlaps = Inf
  ) +
  scale_color_manual(
    name = "Drug Category",
    values = color_values,
    labels = setNames(label_counts$label, label_counts$categoria)
  ) +
  labs(
    title = "Drugs with strongest perturbation across HPV clades",
    subtitle = "Significant Reverse Gene Expression Score (sRGES ≤ -0.2)",
    x = "sRGES clade A7",
    y = "sRGES clade A9"
  ) +
  coord_cartesian(
    xlim = c(min_x, 0),
    ylim = c(min_y, 0)
  ) +
  theme_minimal(base_size = 20) +
  theme(
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18)
  )

plot(p_final)
# ==========================================================
# GUARDAR PDF Y PNG
# ==========================================================
ggsave("Drug_Prediction_Scatter_Plot_Top5_NoBRD.pdf", plot = p_final, width = 14, height = 12)
ggsave("Drug_Prediction_Scatter_Plot_Top5_NoBRD.png", plot = p_final, width = 14, height = 12, dpi = 300)

# ==========================================================
# MOSTRAR EN VENTANA DE PLOTEO DE R
# ==========================================================
print(p_final)









#. IN SILICO VALIDATION PLOTS
if (exists("p_ic50_A7")) print(p_ic50_A7)
if (exists("p_auc_A7")) print(p_auc_A7)
if (exists("p_ic50_A9")) print(p_ic50_A9)
if (exists("p_auc_A9")) print(p_auc_A9)

if (exists("p_ic50_A7_global")) print(p_ic50_A7_global)
if (exists("p_auc_A7_global")) print(p_auc_A7_global)
if (exists("p_ic50_A9_global")) print(p_ic50_A9_global)
if (exists("p_auc_A9_global")) print(p_auc_A9_global)


#. ENRICHMENT HEATMAP (A7 vs A9)
enriched_A7 <- read_csv("~/CESC_Network/6_OCTAD/6_1_Enriquecimiento_A7/chembl_targets/enriched_chembl_targets.csv") %>%
  filter(padj < 0.05) %>%
  select(target, score) %>%
  rename(score_A7 = score)

enriched_A9 <- read_csv("~/CESC_Network/6_OCTAD/6_1_Enriquecimiento_A9/chembl_targets/enriched_chembl_targets.csv") %>%
  filter(padj < 0.05) %>%
  select(target, score) %>%
  rename(score_A9 = score)

enrichment_matrix <- full_join(enriched_A7, enriched_A9, by = "target") %>%
  replace(is.na(.), 0) %>%
  column_to_rownames("target") %>%
  as.matrix()

pheatmap(
  enrichment_matrix,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  scale = "none",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  main = "CHEMBL Target Enrichment Scores (A7 vs A9)",
  fontsize_row = 8,
  fontsize_col = 12
)







