# Load saved environment
load("~/CESC_Network/6_OCTAD/6_3_OCTAD.RData")

# Required libraries
library(pheatmap)
library(octad)
library(octad.db)
library(dplyr)
library(tibble)
library(ggplot2)
library(readr)
library(purrr)
library(ggrepel)
library(tidyr)

# === Similarity Heatmap Between Clades ===
# Ensure rownames
rownames(similitud_A7) <- similitud_A7$cell_id
rownames(similitud_A9) <- similitud_A9$cell_id

# Join similarity matrices
combined_similitud <- full_join(
  similitud_A7 %>% dplyr::select(cell_id, medcor) %>% dplyr::rename(clado_A7 = medcor),
  similitud_A9 %>% dplyr::select(cell_id, medcor) %>% dplyr::rename(clado_A9 = medcor),
  by = "cell_id"
)

# Order by max similarity
combined_similitud$max_cor <- apply(combined_similitud[, c("clado_A7", "clado_A9")], 1, max, na.rm = TRUE)
combined_similitud <- combined_similitud[order(-combined_similitud$max_cor), ]
rownames(combined_similitud) <- combined_similitud$cell_id

# Create similarity matrix
cor_matrix_combined <- as.matrix(combined_similitud[, c("clado_A7", "clado_A9")])
labels <- rownames(cor_matrix_combined)
selected_cells <- union(lineas_similares_A7, lineas_similares_A9)
label_exp <- ifelse(labels %in% selected_cells, paste0("bold('", labels, "')"), paste0("plain('", labels, "')"))

pheatmap(cor_matrix_combined,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "Similarity between the expression of cell lines \n and cervical cancer samples",
         na_col = "grey90",
         labels_row = parse(text = label_exp),
         fontsize_row = 10,
         fontsize_col = 14,
         fontface_col = "bold",
         cellwidth = 65,
         angle_col = 0)


# === Drug Prediction Scatter Plot (sRGES) ===
# Preparar los datos
library(dplyr)
library(ggplot2)
library(ggrepel)

# Clasificación de fármacos según sRGES
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

# Conteo por categoría
label_counts <- scatter_data %>%
  group_by(categoria) %>%
  summarise(n = n()) %>%
  mutate(label = paste0(categoria, " (n=", n, ")"))

# Asignar colores consistentes
color_values <- c("Shared" = "darkgreen", "clade_A9" = "red", "clade_A7" = "steelblue", "Non-significant" = "grey70")

ggplot(scatter_data, aes(x = A7, y = A9, color = categoria)) +
  geom_point(alpha = 0.8, size = 4) +
  geom_hline(yintercept = -0.2, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = -0.2, linetype = "dashed", color = "grey40") +
  geom_text_repel(
    data = filter(scatter_data, categoria == "Shared" & A7 < -0.4 & A9 < -0.4),
    aes(label = pert_iname), size = 4, fontface = "bold", max.overlaps = 10
  ) +
  coord_cartesian(xlim = c(min(scatter_data$A7, na.rm = TRUE), 0.03),
                  ylim = c(min(scatter_data$A9, na.rm = TRUE), 0.03)) +
  scale_color_manual(
    name = "Drug category",
    values = color_values,
    labels = setNames(label_counts$label, label_counts$categoria)
  ) +
  labs(
    title = "Drugs with the strongest perturbation in transcriptomic profiles across HPV-related clades",
    subtitle = "Significant Reverse Gene Expression Score (sRGES ≤ -0.2)",
    x = "sRGES clade A7",
    y = "sRGES clade A9"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 20),
    plot.subtitle = element_text(size = 16, face = "italic"),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14)
  )





# =====================
# 1. Load Required Libraries
# =====================
library(octad)
library(octad.db)
library(dplyr)
library(vroom)
library(tibble)
library(purrr)
library(ggplot2)
library(ggrepel)

# =====================
# 12. In Silico Validation: IC50 and AUC
# =====================

# Function to read validation files
read_insilico_file <- function(cell_line, base_dir, type = c("ic50", "auc")) {
  type <- match.arg(type)
  file_path <- file.path(base_dir, cell_line, "CellLineEval", paste0(cell_line, "_", type, "_insilico_data.tsv"))
  if (file.exists(file_path)) {
    df <- read_tsv(file_path, show_col_types = FALSE)
    df$cell_line <- cell_line
    return(df)
  } else {
    return(NULL)
  }
}

# --- A7 ---
base_dir_A7 <- "~/CESC_Network/6_OCTAD/6_2_Validation_A7/"
cell_lines_A7 <- list.dirs(base_dir_A7, full.names = FALSE, recursive = FALSE)

ic50_A7 <- map_dfr(cell_lines_A7, ~read_insilico_file(.x, base_dir_A7, "ic50")) %>%
  filter(!is.na(medIC50), is.finite(medIC50), sRGES < -0.2)
auc_A7 <- map_dfr(cell_lines_A7, ~read_insilico_file(.x, base_dir_A7, "auc")) %>%
  filter(!is.na(medauc), is.finite(medauc), sRGES < -0.2)

p_ic50_A7 <- ggplot(ic50_A7, aes(x = sRGES, y = log10(medIC50), color = StronglyPredicted)) +
  geom_point() +
  geom_text_repel(aes(label = pert_iname), size = 2.3, max.overlaps = 5, segment.color = 'grey70') +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  facet_wrap(~cell_line, scales = "free") +
  labs(title = "A7: sRGES vs log(IC50) by Cell Line", y = "log10(medIC50)", x = "sRGES") +
  theme_bw()

p_auc_A7 <- ggplot(auc_A7, aes(x = sRGES, y = medauc, color = StronglyPredicted)) +
  geom_point() +
  geom_text_repel(aes(label = pert_iname), size = 2.3, max.overlaps = 5, segment.color = 'grey70') +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  facet_wrap(~cell_line, scales = "free") +
  labs(title = "A7: sRGES vs AUC by Cell Line", y = "AUC", x = "sRGES") +
  theme_bw()

print(p_ic50_A7)
print(p_auc_A7)

# --- A9 ---
base_dir_A9 <- "~/CESC_Network/6_OCTAD/6_2_Validation_A9/"
cell_lines_A9 <- list.dirs(base_dir_A9, full.names = FALSE, recursive = FALSE)

ic50_A9 <- map_dfr(cell_lines_A9, ~read_insilico_file(.x, base_dir_A9, "ic50")) %>%
  filter(!is.na(medIC50), is.finite(medIC50), sRGES < -0.2)
auc_A9 <- map_dfr(cell_lines_A9, ~read_insilico_file(.x, base_dir_A9, "auc")) %>%
  filter(!is.na(medauc), is.finite(medauc), sRGES < -0.2)

p_ic50_A9 <- ggplot(ic50_A9, aes(x = sRGES, y = log10(medIC50), color = StronglyPredicted)) +
  geom_point() +
  geom_text_repel(aes(label = pert_iname), size = 2.3, max.overlaps = 5, segment.color = 'grey70') +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  facet_wrap(~cell_line, scales = "free") +
  labs(title = "A9: sRGES vs log(IC50) by Cell Line", y = "log10(medIC50)", x = "sRGES") +
  theme_bw()

p_auc_A9 <- ggplot(auc_A9, aes(x = sRGES, y = medauc, color = StronglyPredicted)) +
  geom_point() +
  geom_text_repel(aes(label = pert_iname), size = 2.3, max.overlaps = 5, segment.color = 'grey70') +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  facet_wrap(~cell_line, scales = "free") +
  labs(title = "A9: sRGES vs AUC by Cell Line", y = "AUC", x = "sRGES") +
  theme_bw()

print(p_ic50_A9)
print(p_auc_A9)

# Optional Save
# ggsave("facet_ic50_A7.png", plot = p_ic50_A7, width = 12, height = 8, dpi = 300)
# ggsave("facet_auc_A7.png", plot = p_auc_A7, width = 12, height = 8, dpi = 300)
# ggsave("facet_ic50_A9.png", plot = p_ic50_A9, width = 12, height = 8, dpi = 300)
# ggsave("facet_auc_A9.png", plot = p_auc_A9, width = 12, height = 8, dpi = 300)






# =====================
# 13. Global In Silico Validation: A7 and A9 (No Facets)
# =====================

# --- A7 Global Plot ---
p_ic50_A7_global <- ggplot(ic50_A7, aes(x = sRGES, y = log10(medIC50), color = StronglyPredicted)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "grey80") +
  labs(title = "A7: Global Correlation sRGES vs log(IC50)",
       y = "log10(medIC50)", x = "sRGES") +
  theme_minimal(base_size = 14)

p_auc_A7_global <- ggplot(auc_A7, aes(x = sRGES, y = medauc, color = StronglyPredicted)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "grey80") +
  labs(title = "A7: Global Correlation sRGES vs AUC",
       y = "AUC", x = "sRGES") +
  theme_minimal(base_size = 14)

# --- A9 Global Plot ---
p_ic50_A9_global <- ggplot(ic50_A9, aes(x = sRGES, y = log10(medIC50), color = StronglyPredicted)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "grey80") +
  labs(title = "A9: Global Correlation sRGES vs log(IC50)",
       y = "log10(medIC50)", x = "sRGES") +
  theme_minimal(base_size = 14)

p_auc_A9_global <- ggplot(auc_A9, aes(x = sRGES, y = medauc, color = StronglyPredicted)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "grey80") +
  labs(title = "A9: Global Correlation sRGES vs AUC",
       y = "AUC", x = "sRGES") +
  theme_minimal(base_size = 14)

# --- Visualize ---
print(p_ic50_A7_global)
print(p_auc_A7_global)
print(p_ic50_A9_global)
print(p_auc_A9_global)

# --- Optional Save ---
# ggsave("A7_global_ic50.png", plot = p_ic50_A7_global, width = 8, height = 6, dpi = 300)
# ggsave("A7_global_auc.png", plot = p_auc_A7_global, width = 8, height = 6, dpi = 300)
# ggsave("A9_global_ic50.png", plot = p_ic50_A9_global, width = 8, height = 6, dpi = 300)
# ggsave("A9_global_auc.png", plot = p_auc_A9_global, width = 8, height = 6, dpi = 300)
