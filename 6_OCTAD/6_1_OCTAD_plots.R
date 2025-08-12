# ================================================
# 0. Load Environment & Libraries
# ================================================
load("~/CESC_Network/6_OCTAD/6_3_OCTAD.RData")

library(pheatmap)
library(ggplot2)
library(dplyr)
library(tibble)
library(ggrepel)
library(readr)
library(tidyr)
library(purrr)

# ================================================
# 2. SIMILARITY HEATMAP
# ================================================
tryCatch({
  rownames(similitud_A7) <- similitud_A7$cell_id
  rownames(similitud_A9) <- similitud_A9$cell_id
  combined_similitud <- full_join(
    similitud_A7 %>% select(cell_id, medcor) %>% rename(clado_A7 = medcor),
    similitud_A9 %>% select(cell_id, medcor) %>% rename(clado_A9 = medcor),
    by = "cell_id"
  )
  combined_similitud$max_cor <- apply(combined_similitud[, c("clado_A7", "clado_A9")], 1, max, na.rm = TRUE)
  combined_similitud <- combined_similitud[order(-combined_similitud$max_cor), ]
  rownames(combined_similitud) <- combined_similitud$cell_id
  
  cor_matrix_combined <- as.matrix(combined_similitud[, c("clado_A7", "clado_A9")])
  labels <- rownames(cor_matrix_combined)
  selected_cells <- union(lineas_similares_A7, lineas_similares_A9)
  label_exp <- ifelse(labels %in% selected_cells, paste0("bold('", labels, "')"), paste0("plain('", labels, "')"))
  
  # Mostrar el gráfico de Similarity Heatmap
  png("~/CESC_Network/6_OCTAD/Similarity_Heatmap.png", width = 1200, height = 1000)
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
}, error=function(e) message("⚠️ Error en Similarity Heatmap: ", e$message))

# ================================================
# 3. DRUG PREDICTION SCATTER PLOT
# ================================================
# ================================================
# 3. DRUG PREDICTION SCATTER PLOT
# ================================================
# ================================================
# 3. DRUG PREDICTION SCATTER PLOT
# ================================================
# ================================================
# 3. DRUG PREDICTION SCATTER PLOT
# ================================================
# ================================================
# 3. DRUG PREDICTION SCATTER PLOT
# ================================================
# ================================================
# 3. DRUG PREDICTION SCATTER PLOT
# ================================================
tryCatch({
  scatter_data <- full_join(
    result_A7 %>% select(pert_iname, A7 = sRGES),
    result_A9 %>% select(pert_iname, A9 = sRGES),
    by = "pert_iname"
  ) %>%
    mutate(
      categoria = case_when(
        A7 < -0.2 & A9 < -0.2 ~ "Shared",       # verde
        A7 < -0.2 & (is.na(A9) | A9 >= -0.2) ~ "clade_A7",  # azul (Clade A7)
        A9 < -0.2 & (is.na(A7) | A7 >= -0.2) ~ "clade_A9",  # rojo (Clade A9)
        TRUE ~ "Non-significant"
      )
    )
  
  label_counts <- scatter_data %>%
    group_by(categoria) %>%
    summarise(n = n()) %>%
    mutate(label = paste0(categoria, " (n=", n, ")"))
  
  color_values <- c("Shared" = "darkgreen", "clade_A9" = "red", "clade_A7" = "steelblue", "Non-significant" = "grey70")
  
  # Seleccionar los top 5 fármacos de cada categoría sin mezclar entre categorías
  top_shared <- scatter_data %>% filter(categoria == "Shared") %>% top_n(-5, A7) %>% top_n(-5, A9)
  top_clade_A7 <- scatter_data %>% filter(categoria == "clade_A7") %>% top_n(-5, A7)
  top_clade_A9 <- scatter_data %>% filter(categoria == "clade_A9") %>% top_n(-5, A9)
  
  # Unir los top 5 de cada categoría
  top_drugs <- bind_rows(top_shared, top_clade_A7, top_clade_A9)
  
  p <- ggplot(scatter_data, aes(x = A7, y = A9, color = categoria)) +
    geom_point(alpha = 0.8, size = 4) +
    geom_hline(yintercept = -0.2, linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = -0.2, linetype = "dashed", color = "grey40") +
    # Etiquetas para los top 5 fármacos de cada categoría con color negro
    geom_text_repel(
      data = top_drugs,
      aes(label = pert_iname),
      size = 6, 
      fontface = "bold",
      color = "black",  # Etiquetas en negro
      max.overlaps = 10
    ) +
    scale_color_manual(
      name = "Drug category",
      values = color_values,
      labels = setNames(label_counts$label, label_counts$categoria)
    ) +
    labs(
      title = "Drugs with strongest perturbation across HPV clades",
      subtitle = "Significant Reverse Gene Expression Score (sRGES ≤ -0.2)",
      x = "sRGES clade A7",
      y = "sRGES clade A9"
    ) +
    theme_minimal(base_size = 18) +  # Aumentar el tamaño de la letra
    theme(
      legend.title = element_text(size = 20),  # Agrandar el título de la leyenda
      legend.text = element_text(size = 16),   # Agrandar el texto de la leyenda
      axis.title = element_text(size = 16),    # Agrandar los títulos de los ejes
      axis.text = element_text(size = 14)      # Agrandar los números de los ejes
    ) +
    # Limitar los ejes a 0.0 en ambos pero hacia los valores negativos
    xlim(min(c(0, scatter_data$A7), na.rm = TRUE), 0) +
    ylim(min(c(0, scatter_data$A9), na.rm = TRUE), 0)
  
  # Mostrar el gráfico de Drug Prediction Scatter Plot
  ggsave("~/CESC_Network/6_OCTAD/Drug_Prediction_Scatter_Plot.png", plot = p, width = 12, height = 10, dpi = 300)
}, error=function(e) message("⚠️ Error en Drug Scatter Plot: ", e$message))

plot(p)

# ================================================
# 4. IN SILICO VALIDATION PLOTS (IF THEY EXIST)
# ================================================
for (plot_obj in c("p_ic50_A7", "p_auc_A7", "p_ic50_A9", "p_auc_A9",
                   "p_ic50_A7_global", "p_auc_A7_global", "p_ic50_A9_global", "p_auc_A9_global")) {
  tryCatch({
    if (exists(plot_obj)) {
      # Guardar cada gráfico en PNG por separado
      ggsave(paste0("~/CESC_Network/6_OCTAD/", plot_obj, ".png"), plot = get(plot_obj), width = 12, height = 10, dpi = 300)
    } else {
      message("⚠️ Plot not found: ", plot_obj)
    }
  }, error=function(e) message("⚠️ Error printing ", plot_obj, ": ", e$message))
}

# ================================================
# 5. ENRICHMENT HEATMAP (A7 vs A9)
# ================================================
tryCatch({
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
  
  # Mostrar el gráfico de Enrichment Heatmap
  png("~/CESC_Network/6_OCTAD/Enrichment_Heatmap_A7_vs_A9.png", width = 1200, height = 1000)
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
  dev.off()
}, error=function(e) message("⚠️ Error en Enrichment Heatmap: ", e$message))









# ================================================
# 3. DRUG PREDICTION SCATTER PLOT (CORREGIDO PARA AZULES)
# ================================================
scatter_data <- full_join(
  result_A7 %>% select(pert_iname, A7 = sRGES),
  result_A9 %>% select(pert_iname, A9 = sRGES),
  by = "pert_iname"
) %>%
  mutate(
    categoria = case_when(
      A7 < -0.2 & A9 < -0.2 ~ "Shared",
      A7 < -0.2 & (is.na(A9) | A9 >= -0.2) ~ "clade_A7",
      A9 < -0.2 & (is.na(A7) | A7 >= -0.2) ~ "clade_A9",
      TRUE ~ "Non-significant"
    )
  )

# Top 5 bien seleccionados por categoría
top_shared <- scatter_data %>% 
  filter(categoria == "Shared", !is.na(A7), !is.na(A9)) %>% 
  arrange(A7 + A9) %>% 
  slice(1:5)

top_clade_A7 <- scatter_data %>% 
  filter(categoria == "clade_A7", !is.na(A7)) %>% 
  arrange(A7) %>% 
  slice(1:5)
setdiff()
top_clade_A9 <- scatter_data %>% 
  filter(categoria == "clade_A9", !is.na(A9)) %>% 
  arrange(A9) %>% 
  slice(1:5)

top_drugs <- bind_rows(top_shared, top_clade_A7, top_clade_A9)

# Etiquetas de leyenda con (n=)
label_counts <- scatter_data %>%
  group_by(categoria) %>%
  summarise(n = n()) %>%
  mutate(label = paste0(categoria, " (n=", n, ")"))

# Colores
color_values <- c("Shared" = "darkgreen", "clade_A9" = "red", "clade_A7" = "steelblue", "Non-significant" = "grey70")

# Gráfico
p <- ggplot(scatter_data, aes(x = A7, y = A9, color = categoria)) +
  geom_point(alpha = 0.8, size = 4) +
  geom_hline(yintercept = -0.2, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = -0.2, linetype = "dashed", color = "grey40") +
  geom_text_repel(
    data = top_drugs,
    aes(label = pert_iname),
    size = 7,
    fontface = "bold",
    color = "black",
    max.overlaps = 30
  ) +
  scale_color_manual(
    name = "Drug category",
    values = color_values,
    labels = setNames(label_counts$label, label_counts$categoria)
  ) +
  labs(
    title = "Drugs with strongest perturbation across HPV clades",
    subtitle = "Significant Reverse Gene Expression Score (sRGES ≤ -0.2)",
    x = "sRGES clade A7",
    y = "sRGES clade A9"
  ) +
  theme_minimal(base_size = 22) +
  theme(
    legend.title = element_text(size = 26),
    legend.text = element_text(size = 22),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20)
  ) +
  coord_cartesian(
    xlim = c(min(0, scatter_data$A7, na.rm = TRUE), 0),
    ylim = c(min(0, scatter_data$A9, na.rm = TRUE), 0)
  )

# Guardar
ggsave("~/CESC_Network/6_OCTAD/Drug_Prediction_Scatter_Plot_Top5_AZULES_OK.png", plot = p, width = 12, height = 10, dpi = 300)

p







UNICOS <- setdiff(result_A7_0.2$pert_iname, result_A9_0.2$pert_iname)
print(UNICOS)






























# ================================================
# 3. DRUG PREDICTION SCATTER PLOT (CORRECTO CON SETDIFF Y INTERSECT)
# ================================================

# Definir conjuntos
UNICOS_A7 <- setdiff(result_A7_0.2$pert_iname, result_A9_0.2$pert_iname)
UNICOS_A9 <- setdiff(result_A9_0.2$pert_iname, result_A7_0.2$pert_iname)
SHARED <- intersect(result_A7_0.2$pert_iname, result_A9_0.2$pert_iname)

# Construir tabla combinada
scatter_data <- full_join(
  result_A7 %>% select(pert_iname, A7 = sRGES),
  result_A9 %>% select(pert_iname, A9 = sRGES),
  by = "pert_iname"
) %>%
  mutate(
    categoria = case_when(
      pert_iname %in% SHARED ~ "Shared",
      pert_iname %in% UNICOS_A7 ~ "clade_A7",
      pert_iname %in% UNICOS_A9 ~ "clade_A9",
      TRUE ~ "Non-significant"
    )
  )

# Filtrar top 5 de cada categoría real
top_shared <- scatter_data %>%
  filter(categoria == "Shared", !is.na(A7), !is.na(A9)) %>%
  arrange(A7 + A9) %>%
  slice(1:5)

top_clade_A7 <- scatter_data %>%
  filter(categoria == "clade_A7", !is.na(A7)) %>%
  arrange(A7) %>%
  slice(1:5)

top_clade_A9 <- scatter_data %>%
  filter(categoria == "clade_A9", !is.na(A9)) %>%
  arrange(A9) %>%
  slice(1:5)

top_drugs <- bind_rows(top_shared, top_clade_A7, top_clade_A9)

# Etiquetas de leyenda
label_counts <- scatter_data %>%
  group_by(categoria) %>%
  summarise(n = n()) %>%
  mutate(label = paste0(categoria, " (n=", n, ")"))

# Colores
color_values <- c("Shared" = "darkgreen", "clade_A9" = "red", "clade_A7" = "steelblue", "Non-significant" = "grey70")

# Gráfico
p <- ggplot(scatter_data, aes(x = A7, y = A9, color = categoria)) +
  geom_point(alpha = 0.8, size = 4) +
  geom_hline(yintercept = -0.2, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = -0.2, linetype = "dashed", color = "grey40") +
  geom_text_repel(
    data = top_drugs,
    aes(label = pert_iname),
    size = 7,
    fontface = "bold",
    color = "black",
    max.overlaps = 30
  ) +
  scale_color_manual(
    name = "Drug category",
    values = color_values,
    labels = setNames(label_counts$label, label_counts$categoria)
  ) +
  labs(
    title = "Drugs with strongest perturbation across HPV clades",
    subtitle = "Significant Reverse Gene Expression Score (sRGES ≤ -0.2)",
    x = "sRGES clade A7",
    y = "sRGES clade A9"
  ) +
  theme_minimal(base_size = 22) +
  theme(
    legend.title = element_text(size = 26),
    legend.text = element_text(size = 22),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20)
  ) +
  coord_cartesian(
    xlim = c(min(0, scatter_data$A7, na.rm = TRUE), 0),
    ylim = c(min(0, scatter_data$A9, na.rm = TRUE), 0)
  )

# Guardar
ggsave("~/CESC_Network/6_OCTAD/Drug_Prediction_Scatter_Plot_Top5_FINAL_SETDIFF.png", plot = p, width = 12, height = 10, dpi = 300)

p








# ================================================
# 3. DRUG PREDICTION SCATTER PLOT (COMPLETO CON AZULES)
# ================================================

# 1. Definir conjuntos
UNICOS_A7 <- setdiff(result_A7_0.2$pert_iname, result_A9_0.2$pert_iname)
UNICOS_A9 <- setdiff(result_A9_0.2$pert_iname, result_A7_0.2$pert_iname)
SHARED <- intersect(result_A7_0.2$pert_iname, result_A9_0.2$pert_iname)

# 2. Construir tabla combinada
scatter_data <- full_join(
  result_A7 %>% select(pert_iname, A7 = sRGES),
  result_A9 %>% select(pert_iname, A9 = sRGES),
  by = "pert_iname"
) %>%
  mutate(
    categoria = case_when(
      pert_iname %in% SHARED ~ "Shared",
      pert_iname %in% UNICOS_A7 ~ "clade_A7",
      pert_iname %in% UNICOS_A9 ~ "clade_A9",
      TRUE ~ "Non-significant"
    )
  )

# 3. Filtrar top 5 para etiquetar
top_shared <- scatter_data %>%
  filter(categoria == "Shared", !is.na(A7), !is.na(A9)) %>%
  arrange(A7 + A9) %>%
  slice(1:5)

top_clade_A7 <- scatter_data %>%
  filter(pert_iname %in% UNICOS_A7, !is.na(A7)) %>%
  arrange(A7) %>%
  slice(1:5)

top_clade_A9 <- scatter_data %>%
  filter(pert_iname %in% UNICOS_A9, !is.na(A9)) %>%
  arrange(A9) %>%
  slice(1:5)

top_drugs <- bind_rows(top_shared, top_clade_A7, top_clade_A9)

# 4. Leyenda con conteos
label_counts <- scatter_data %>%
  group_by(categoria) %>%
  summarise(n = n()) %>%
  mutate(label = paste0(categoria, " (n=", n, ")"))

# 5. Colores
color_values <- c(
  "Shared" = "darkgreen",
  "clade_A9" = "red",
  "clade_A7" = "steelblue",
  "Non-significant" = "grey70"
)

# 6. Plot
p <- ggplot(scatter_data, aes(x = A7, y = A9, color = categoria)) +
  geom_point(alpha = 0.8, size = 4) +
  geom_hline(yintercept = -0.2, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = -0.2, linetype = "dashed", color = "grey40") +
  geom_text_repel(
    data = top_drugs,
    aes(label = pert_iname),
    size = 7,
    fontface = "bold",
    color = "black",
    max.overlaps = 30
  ) +
  scale_color_manual(
    name = "Drug category",
    values = color_values,
    labels = setNames(label_counts$label, label_counts$categoria)
  ) +
  labs(
    title = "Drugs with strongest perturbation across HPV clades",
    subtitle = "Significant Reverse Gene Expression Score (sRGES ≤ -0.2)",
    x = "sRGES clade A7",
    y = "sRGES clade A9"
  ) +
  theme_minimal(base_size = 22) +
  theme(
    legend.title = element_text(size = 26),
    legend.text = element_text(size = 22),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20)
  ) +
  coord_cartesian(
    xlim = c(min(0, scatter_data$A7, na.rm = TRUE), 0),
    ylim = c(min(0, scatter_data$A9, na.rm = TRUE), 0)
  )

# 7. Guardar
ggsave(
  "~/CESC_Network/6_OCTAD/Drug_Prediction_Scatter_Plot_Final_AZULES_OK.png",
  plot = p,
  width = 12,
  height = 10,
  dpi = 300
)

# 8. Mostrar
p
