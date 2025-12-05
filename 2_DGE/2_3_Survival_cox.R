############################################################
## Cox univariado por gen (expresión continua)
## A7_clade y A9_clade
############################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(survival)
  library(vroom) })

## =========================
## 1) Cargar expresión
## =========================

expr_df <- readr::read_tsv(
  "~/CESC_Network/2_Prepro_TCGA_GTEx/2_6_Full_counts_A7_A9_annot.tsv",
  col_types = cols()
)

dim(expr_df)   # 17481 x 281
# Primera columna = Gene (símbolo); el resto son muestras

# Pasar a data.frame con rownames = Gene
expr_mat <- as.data.frame(expr_df)
rownames(expr_mat) <- expr_mat$Gene
expr_mat$Gene <- NULL

## =========================
## 2) Cargar factores
## =========================

factors <- readr::read_tsv(
  "~/CESC_Network/2_Prepro_TCGA_GTEx/2_4_Factors.tsv",
  col_types = cols()
)
factors <- factors %>% filter(!HPV_clade == "Solid Tissue Normal")
# Checar clados / tipos
table(factors$HPV_clade)

#table(factors$sample_type)

## =========================
## 3) Preparar DGE (A7 y A9)
## =========================

DE_A7_vs_SNT <- vroom::vroom(
  "~/CESC_Network/4_Network_analisis/4_3_limma_A7_vs_SNT_allGenes_annot_PR_DEG.tsv",
  show_col_types = FALSE)

DE_A9_vs_SNT <- vroom::vroom(
  "~/CESC_Network/4_Network_analisis/4_3_limma_A9_vs_SNT_allGenes_annot_PR_DEG.tsv",
  show_col_types = FALSE)

# Nos vamos a quedar (por ejemplo) con genes con adj.P.Val < 0.05
# (ajusta el umbral si quieres ser más/menos estricto)
# Filtro + TOP 100 por |logFC|
dge_genes_A7 <- DE_A7_vs_SNT %>%
  filter(logFC < -1 | logFC > 1) %>%      # corte en logFC
  filter(adj.P.Val < 0.05) %>%            # FDR limma
  filter(Pagerank_group == "High_pagerank") %>%  # hubs por PR
  slice_max(order_by = abs(logFC), n = 200, with_ties = TRUE) %>%  # NUEVO: top 100
  pull(Gene) %>%
  unique()

dge_genes_A9 <- DE_A9_vs_SNT %>%
  filter(logFC < -1 | logFC > 1) %>%
  filter(adj.P.Val < 0.05) %>%
  filter(Pagerank_group == "High_pagerank") %>%
  slice_max(order_by = abs(logFC), n = 200, with_ties = TRUE) %>%  # NUEVO: top 100
  pull(Gene) %>%
  unique()

length(dge_genes_A7) #[1] 100
length(dge_genes_A9) #[1] 100

## =========================
## 4) Construir datos de supervivencia para A7
## =========================

# Solo tumores primarios del clado A7
surv_A7 <- factors %>%
  filter(
    HPV_clade   == "A7_clade",
    sample_type == "Primary Tumor"
  ) %>%
  mutate(
    # status: 1 = muere, 0 = censurado
    vital_status = toupper(vital_status),
    status = case_when(
      vital_status %in% c("DECEASED", "DEAD") ~ 1,
      vital_status %in% c("LIVING", "ALIVE")  ~ 0,
      TRUE ~ NA_real_
    ),
    # tiempo de seguimiento en días:
    time_days = dplyr::coalesce(days_to_death, days_to_last_followup)
  ) %>%
  filter(!is.na(time_days), !is.na(status))

nrow(surv_A7)

# Alinear expresión: columnas de expr_mat que correspondan a estos sampleID
common_samples_A7 <- intersect(colnames(expr_mat), surv_A7$sampleID)
length(common_samples_A7)

surv_A7 <- surv_A7 %>%
  filter(sampleID %in% common_samples_A7) %>%
  arrange(match(sampleID, common_samples_A7))

expr_A7 <- expr_mat[, common_samples_A7, drop = FALSE]
stopifnot(all(colnames(expr_A7) == surv_A7$sampleID))

# Intersección de genes DGE significativos con los presentes en la matriz
genes_use_A7 <- intersect(rownames(expr_A7), dge_genes_A7)
length(genes_use_A7)

# Objeto de supervivencia
Surv_A7_obj <- Surv(time = surv_A7$time_days, event = surv_A7$status)

## =========================
## 5) Cox univariado continuo para A7
## =========================

res_list_A7 <- lapply(genes_use_A7, function(g) {
  expr_vec <- as.numeric(expr_A7[g, ])
  
  # Si el gen es constante (sin varianza), saltamos
  if (sd(expr_vec, na.rm = TRUE) == 0) {
    return(NULL)
  }
  
  # Escalamos la expresión (HR será por 1 desviación estándar)
  expr_scaled <- as.numeric(scale(expr_vec))
  
  fit <- tryCatch(
    coxph(Surv_A7_obj ~ expr_scaled),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(NULL)
  }
  
  s <- summary(fit)
  
  data.frame(
    gene       = g,
    HR         = s$coef[1, "exp(coef)"],
    HR_lower95 = s$conf.int[1, "lower .95"],
    HR_upper95 = s$conf.int[1, "upper .95"],
    p_value    = s$coef[1, "Pr(>|z|)"],
    n_patients = s$n,
    stringsAsFactors = FALSE
  )
})

# Unir resultados y calcular FDR
cox_A7 <- do.call(rbind, res_list_A7)
cox_A7$FDR <- p.adjust(cox_A7$p_value, method = "BH")

# Ordenar por FDR
cox_A7 <- cox_A7 %>%
  arrange(FDR)

head(cox_A7)

# Opcional: genes significativos
cox_A7_sig <- cox_A7 %>%
  filter(FDR < 0.05)
head(cox_A7_sig)

## =========================
## 6) Repetir todo para A9
## =========================

surv_A9 <- factors %>%
  filter(
    HPV_clade   == "A9_clade",
    sample_type == "Primary Tumor"
  ) %>%
  mutate(
    vital_status = toupper(vital_status),
    status = case_when(
      vital_status %in% c("DECEASED", "DEAD") ~ 1,
      vital_status %in% c("LIVING", "ALIVE")  ~ 0,
      TRUE ~ NA_real_
    ),
    time_days = dplyr::coalesce(days_to_death, days_to_last_followup)
  ) %>%
  filter(!is.na(time_days), !is.na(status))

nrow(surv_A9)

common_samples_A9 <- intersect(colnames(expr_mat), surv_A9$sampleID)
length(common_samples_A9)

surv_A9 <- surv_A9 %>%
  filter(sampleID %in% common_samples_A9) %>%
  arrange(match(sampleID, common_samples_A9))

expr_A9 <- expr_mat[, common_samples_A9, drop = FALSE]
stopifnot(all(colnames(expr_A9) == surv_A9$sampleID))

genes_use_A9 <- intersect(rownames(expr_A9), dge_genes_A9)
length(genes_use_A9)

Surv_A9_obj <- Surv(time = surv_A9$time_days, event = surv_A9$status)

res_list_A9 <- lapply(genes_use_A9, function(g) {
  expr_vec <- as.numeric(expr_A9[g, ])
  
  if (sd(expr_vec, na.rm = TRUE) == 0) {
    return(NULL)
  }
  
  expr_scaled <- as.numeric(scale(expr_vec))
  
  fit <- tryCatch(
    coxph(Surv_A9_obj ~ expr_scaled),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(NULL)
  }
  
  s <- summary(fit)
  
  data.frame(
    gene       = g,
    HR         = s$coef[1, "exp(coef)"],
    HR_lower95 = s$conf.int[1, "lower .95"],
    HR_upper95 = s$conf.int[1, "upper .95"],
    p_value    = s$coef[1, "Pr(>|z|)"],
    n_patients = s$n,
    stringsAsFactors = FALSE
  )
})

cox_A9 <- do.call(rbind, res_list_A9)
cox_A9$FDR <- p.adjust(cox_A9$p_value, method = "BH")
cox_A9 <- cox_A9 %>%
  arrange(FDR)

head(cox_A9)

cox_A9_sig <- cox_A9 %>%
  filter(FDR < 0.05)

head(cox_A9_sig)
head(cox_A7_sig)

## =========================
## 7) Guardar resultados
## =========================

#readr::write_tsv(cox_A7,     "~/CESC_Network/2_DGE/2_2_Cox_continuous_A7_allDGE.tsv")
#readr::write_tsv(cox_A7_sig, "~/CESC_Network/2_DGE/2_2_Cox_continuous_A7_FDR05.tsv")

#readr::write_tsv(cox_A9,     "~/CESC_Network/2_DGE/2_2_Cox_continuous_A9_allDGE.tsv")
#readr::write_tsv(cox_A9_sig, "~/CESC_Network/2_DGE/2_2_Cox_continuous_A9_FDR05.tsv")
