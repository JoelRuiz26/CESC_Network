############################################################
## SCRIPT: 5_3_2_DeltaCon_gene_contrib_A7_vs_A9.R
############################################################

rm(list = ls())
gc()

library(Matrix)
library(tidyverse)
library(ggrepel)

#====================== PATHS Y CARGA DE MATRICES =============================#

out_dir <- "~/CESC_Network/4_Network_analisis/5_3_DeltaCon/"
S_A7 <- readRDS(file.path(out_dir, "5_3_1_S_A7_full.rds"))
S_A9 <- readRDS(file.path(out_dir, "5_3_1_S_A9_full.rds"))
delta_con_full <- readRDS(file.path(out_dir,
                                    "5_3_1_DeltaCon_distance_A7_A9_full.rds"))

stopifnot(identical(rownames(S_A7), rownames(S_A9)))
stopifnot(identical(colnames(S_A7), colnames(S_A9)))

genes <- rownames(S_A7)
n <- length(genes)
cat("Número de genes:", n, "\n")
cat("DeltaCon FULL cargado:", delta_con_full, "\n\n")

#====================== FUNCIÓN: CONTRIBUCIÓN POR GEN ========================#

compute_gene_stats <- function(S1, S2, block_size = 800) {
  stopifnot(all(dim(S1) == dim(S2)))
  
  n <- nrow(S1)
  blocks <- split(seq_len(n), ceiling(seq_along(seq_len(n)) / block_size))
  
  contrib   <- numeric(n)
  affinity1 <- numeric(n)
  affinity2 <- numeric(n)
  
  cat(">> Calculando afinidades y contribuciones por bloques (",
      length(blocks), " bloques)...\n", sep = "")
  
  for (b in seq_along(blocks)) {
    idx <- blocks[[b]]
    cat("   - Bloque", b, "de", length(blocks),
        "(filas:", min(idx), "-", max(idx), ")\n")
    
    B1 <- as.matrix(S1[idx, , drop = FALSE])
    B2 <- as.matrix(S2[idx, , drop = FALSE])
    
    B1 <- pmax(B1, 0)
    B2 <- pmax(B2, 0)
    
    sqrtB1 <- sqrt(B1)
    sqrtB2 <- sqrt(B2)
    
    affinity1[idx] <- rowSums(sqrtB1)
    affinity2[idx] <- rowSums(sqrtB2)
    
    diff <- sqrtB1 - sqrtB2
    contrib[idx] <- rowSums(diff^2)
    
    rm(B1, B2, sqrtB1, sqrtB2, diff)
    gc()
  }
  
  list(
    contrib   = contrib,
    affinity1 = affinity1,
    affinity2 = affinity2
  )
}

stats_list <- compute_gene_stats(S_A7, S_A9, block_size = 700)

gene_contrib   <- stats_list$contrib
affinity_A7    <- stats_list$affinity1
affinity_A9    <- stats_list$affinity2

rm(stats_list)
gc()

# Comprobación: la suma de contrib debería aproximar DeltaCon^2
delta_from_genes <- sqrt(sum(gene_contrib))
cat("\nDeltaCon reconstruido a partir de contribuciones por gen:",
    delta_from_genes, "\n")
cat("Diferencia con DeltaCon FULL original:",
    delta_from_genes - as.numeric(delta_con_full), "\n\n")

#====================== TABLA: RANKEO Y Z-SCORE ROBUSTO ======================#

total_contrib <- sum(gene_contrib)

contrib_df <- tibble(
  Gene        = genes,
  Affinity_A7 = affinity_A7,
  Affinity_A9 = affinity_A9,
  Contribution = gene_contrib,
  Percent      = 100 * gene_contrib / total_contrib
)

mad_val <- mad(contrib_df$Contribution)
if (mad_val == 0) {
  warning("MAD == 0; usando sd en lugar de MAD para el Z-score.")
  contrib_df <- contrib_df %>%
    mutate(
      Z_robust = (Contribution - mean(Contribution)) / sd(Contribution)
    )
} else {
  contrib_df <- contrib_df %>%
    mutate(
      Z_robust = (Contribution - median(Contribution)) / mad_val
    )
}

contrib_df <- contrib_df %>%
  mutate(
    Higher_affinity = if_else(Affinity_A7 > Affinity_A9, "A7", "A9")
  ) %>%
  arrange(desc(abs(Z_robust)))

write_tsv(
  contrib_df,
  file.path(out_dir, "5_3_2_DeltaCon_gene_contributions_A7_vs_A9.tsv")
)

cat("Top 10 genes por |Z_robust|:\n")
print(contrib_df %>% slice_head(n = 10))

#====================== DATA PARA EL PLOT ====================================#

highlight_df <- contrib_df %>% filter(abs(Z_robust) >= 3)

cat("\nNúmero de genes con |Z_robust| >= 3:", nrow(highlight_df), "\n")

label_df <- highlight_df %>% slice_max(order_by = abs(Z_robust), n = 15)

#====================== PLOT: AFINIDAD A7 vs A9 ==============================#

p_affinity <- ggplot(contrib_df, aes(x = Affinity_A7, y = Affinity_A9)) +
  geom_point(
    aes(color = Higher_affinity),
    alpha = 0.25,
    size  = 1.8
  ) +
  geom_point(
    data = highlight_df,
    aes(color = Higher_affinity),
    alpha = 0.9,
    size  = 3.5
  ) +
  geom_text_repel(
    data = label_df,
    aes(label = Gene),
    size = 3.2,
    fontface = "bold",
    max.overlaps = 50,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "gray55",
    segment.alpha = 0.7
  ) +
  scale_color_manual(
    values = c("A7" = "#E41A1C", "A9" = "#377EB8"),
    name   = "Higher-affinity\nnetwork"
  ) +
  geom_abline(
    slope = 1, intercept = 0,
    linetype = "dashed",
    color = "gray40"
  ) +
  labs(
    title = "Gene-level DeltaCon contribution (A7 vs A9)",
    # 🔴 AQUÍ EL CAMBIO: sin símbolo ≥, solo ASCII
    subtitle = "Highlighted genes: |robust Z-score of contribution| >= 3",
    x = "Affinity in A7 network",
    y = "Affinity in A9 network"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title   = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle= element_text(size = 13, hjust = 0.5),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11)
  )

p_affinity

#====================== GUARDAR PLOTS ========================================#

ggsave(
  filename = file.path(out_dir, "5_3_2_DeltaCon_affinity_A7_vs_A9_chido.pdf"),
  plot     = p_affinity,
  width    = 9,
  height   = 7
)

ggsave(
  filename = file.path(out_dir, "5_3_2_DeltaCon_affinity_A7_vs_A9_chido.png"),
  plot     = p_affinity,
  width    = 9,
  height   = 7,
  dpi      = 300
)

cat("\nPlot y tabla de contribuciones guardados en:\n", out_dir, "\n")
