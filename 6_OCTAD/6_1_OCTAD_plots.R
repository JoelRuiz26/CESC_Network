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
# 1. OPEN PDF
# ================================================
pdf("~/CESC_Network/6_OCTAD/All_OCTAD_Plots.pdf", width=12, height=10)

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
}, error=function(e) message("⚠️ Error en Similarity Heatmap: ", e$message))

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
        A7 < -0.2 & A9 < -0.2 ~ "Shared",
        A7 < -0.2 & (is.na(A9) | A9 >= -0.2) ~ "clade_A7",
        A9 < -0.2 & (is.na(A7) | A7 >= -0.2) ~ "clade_A9",
        TRUE ~ "Non-significant"
      )
    )
  
  label_counts <- scatter_data %>%
    group_by(categoria) %>%
    summarise(n = n()) %>%
    mutate(label = paste0(categoria, " (n=", n, ")"))
  
  color_values <- c("Shared" = "darkgreen", "clade_A9" = "red", "clade_A7" = "steelblue", "Non-significant" = "grey70")
  
  p <- ggplot(scatter_data, aes(x = A7, y = A9, color = categoria)) +
    geom_point(alpha = 0.8, size = 4) +
    geom_hline(yintercept = -0.2, linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = -0.2, linetype = "dashed", color = "grey40") +
    geom_text_repel(
      data = filter(scatter_data, categoria == "Shared" & A7 < -0.4 & A9 < -0.4),
      aes(label = pert_iname), size = 4, fontface = "bold", max.overlaps = 10
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
    theme_minimal(base_size = 16)
  print(p)
}, error=function(e) message("⚠️ Error en Drug Scatter Plot: ", e$message))


# ================================================
# 4. IN SILICO VALIDATION PLOTS (IF THEY EXIST)
# ================================================
for (plot_obj in c("p_ic50_A7", "p_auc_A7", "p_ic50_A9", "p_auc_A9",
                   "p_ic50_A7_global", "p_auc_A7_global", "p_ic50_A9_global", "p_auc_A9_global")) {
  tryCatch({
    if (exists(plot_obj)) {
      print(get(plot_obj))
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
}, error=function(e) message("⚠️ Error en Enrichment Heatmap: ", e$message))

# ================================================
# 6. CLOSE PDF
# ================================================
dev.off()
