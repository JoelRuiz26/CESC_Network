# =========================================
# Load Required Libraries
# =========================================
library(octad)
library(octad.db)
library(dplyr)
library(vroom)
setwd("~/OCTAD/")

# =========================================
# Load OCTAD Preprocessed Data
# =========================================
count_matrix <- vroom("~/OCTAD/count_matrix_octad.tsv")
metadata <- vroom("~/OCTAD/metadata_octad.tsv")

# =========================================
# Load LINCS Phenotype Data
# =========================================
phenoDF <- get_ExperimentHub_data("EH7274")

# =========================================
# A7 Clade Processing
# =========================================
muestras_clado_A7 <- metadata %>%
  filter(HPV_clade == "A7") %>%
  pull(cases.submitter_id)

expSet_A7 <- count_matrix[, colnames(count_matrix) %in% muestras_clado_A7]
dim(expSet_A7)  #[1] 11240    66

interseccion_A7 <- intersect(muestras_clado_A7, phenoDF$sample.id)
length(interseccion_A7) #[1] 65
setdiff(muestras_clado_A7, interseccion_A7)  # [1] "TCGA-C5-A7CM-01"


similitud_A7 <- computeCellLine(case_id = interseccion_A7, source = 'octad.small')
similitud_A7 <- similitud_A7 %>% mutate(cell_id = rownames(similitud_A7))
dim(similitud_A7) #[1] 51  2
lineas_similares_A7 <- similitud_A7 %>% filter(medcor > 0.3) %>% pull(cell_id)
length(lineas_similares_A7) #[1] 15

# =========================================
# A9 Clade Processing
# =========================================
muestras_clado_A9 <- metadata %>%
  filter(HPV_clade == "A9") %>%
  pull(cases.submitter_id)

expSet_A9 <- count_matrix[, colnames(count_matrix) %in% muestras_clado_A9]
dim(expSet_A9) #[1] 11240   202

interseccion_A9 <- intersect(muestras_clado_A9, phenoDF$sample.id)
length(interseccion_A9) #[1] 202

similitud_A9 <- computeCellLine(case_id = interseccion_A9, source = 'octad.small')
similitud_A9 <- similitud_A9 %>% mutate(cell_id = rownames(similitud_A9))
lineas_similares_A9 <- similitud_A9 %>% filter(medcor > 0.3) %>% pull(cell_id)
length(lineas_similares_A9) #[1] 11

# =========================================
# Load Differential Expression Signatures
# =========================================
dge_A7 <- read.delim("~/3_DGE/DESeq2/DGE_A7_vs_SNT_GSEA.tsv") %>%
  filter(FDR < 0.05)
firma_A7 <- data.frame(
  Symbol = dge_A7$gene,
  log2FoldChange = dge_A7$logFC,
  regulation = case_when(
    dge_A7$logFC > 1 ~ "up",        
    dge_A7$logFC < -1 ~ "down",     
    TRUE ~ NA_character_)) %>% filter(!is.na(regulation))    
dim(firma_A7) #[1] 2950    3

dge_A9 <- read.delim("~/3_DGE/DESeq2/DGE_A9_vs_SNT_GSEA.tsv") %>%
  filter(FDR < 0.05)
firma_A9 <- data.frame(
  Symbol = dge_A9$gene,
  log2FoldChange = dge_A9$logFC,
  regulation = case_when(
    dge_A9$logFC > 1 ~ "up",        
    dge_A9$logFC < -1 ~ "down",     
    TRUE ~ NA_character_)) %>% filter(!is.na(regulation))    
dim(firma_A9) #[1] 3059    3

# =========================================
# Compute sRGES Scores (LINCS-Based Signature Reversal)
# =========================================
result_A7 <- octad::runsRGES(
  dz_signature = firma_A7,
  cells = lineas_similares_A7,
  weight_cell_line = similitud_A7,
  permutations = 1000,
  choose_fda_drugs = TRUE
)

result_A9 <- runsRGES(
  dz_signature = firma_A9,
  cells = lineas_similares_A9,
  weight_cell_line = similitud_A9,
  permutations = 1000,
  choose_fda_drugs = TRUE
)

# =========================================
# Null Model: Save Inputs for Reuse in Other Scripts
# =========================================
#save.image("~/OCTAD/Input_nullmodel.RData")

# =========================================
# Load Results from Previous Null Model Runs
# =========================================
# result_A7_pre <- readRDS(result_A9_sig, "~/OCTAD/p_empirico_A7.rds")
# result_A9_pre <- readRDS(result_A9_sig, "~/OCTAD/p_empirico_A9.rds")
# result_A7_sig <- readRDS(result_A9_sig, "~/OCTAD/result_A7_sig.rds")
# result_A9_sig <- readRDS(result_A9_sig, "~/OCTAD/result_A9_sig.rds")

# result_A7 <- result_A7_sig
# result_A9 <- result_A9_sig

# =========================================
# Filter Drugs with sRGES < -0.2
# =========================================
result_A7_0.2 <- result_A7 %>% filter(sRGES < -0.15)
dim(result_A7_0.2) #[1] 117   6
result_A9_0.2 <- result_A9 %>% filter(sRGES < -0.15) 
dim(result_A9_0.2) #[1] 140   6

length(intersect(result_A7_0.2$pert_iname,result_A9_0.2$pert_iname)) #[1] 51

# =========================================
# Enrichment Analysis of Drug Targets (CHEMBL)
# =========================================

# A7 enrichment
octadDrugEnrichment(
  sRGES = result_A7_0.2,
  target_type = "chembl_targets",
  enrichFolder = "Enriquecimiento_A7",
  outputFolder = "~/OCTAD"
)
enriched_chembl_A7 <- vroom("~/OCTAD/Enriquecimiento_A7/chembl_targets/enriched_chembl_targets.csv")
dim(enriched_chembl_A7) #[1] 126   4
enriched_chembl_A7_filtered <- enriched_chembl_A7 %>% filter(padj < 0.05)
dim(enriched_chembl_A7_filtered)#[1] 64  4
head(enriched_chembl_A7_filtered)
# A tibble: 6 × 4
#target  score     p  padj
#<chr>   <dbl> <dbl> <dbl>
#1 ABCB1   0.484     0     0
#2 ACHE    0.387     0     0
#3 ADORA2A 0.462     0     0
#4 ADRA1D  0.462     0     0
#5 ADRA2A  0.196     0     0
#6 ADRA2B  0.462     0     0


# A9 enrichment
octadDrugEnrichment(
  sRGES = result_A9_0.2,
  target_type = "chembl_targets",
  enrichFolder = "Enriquecimiento_A9",
  outputFolder = "~/OCTAD"
)
enriched_chembl_A9 <- vroom("~/OCTAD/Enriquecimiento_A9/chembl_targets/enriched_chembl_targets.csv")
enriched_chembl_A9_filtered <- enriched_chembl_A9 %>% filter(padj < 0.05)
dim(enriched_chembl_A9_filtered)
#[1] 71  4
head(enriched_chembl_A9_filtered)
# A tibble: 6 × 4
#target  score     p  padj
#<chr>   <dbl> <dbl> <dbl>
#1 ABCB1   0.291     0     0
#2 ABCG2   0.325     0     0
#3 ACHE    0.287     0     0
#4 ADORA2A 0.297     0     0
#5 ADRA1D  0.287     0     0
#6 ADRA2B  0.297     0     0

# =========================================
# Validation Step (According to OCTAD pipeline)
# =========================================
# OCTAD recommends validation via:
# 1. Literature mining of top-ranked compounds
# 2. Cross-referencing with clinical trial databases (e.g. ClinicalTrials.gov)
# 3. Integration with adverse effect databases (e.g. FAERS)

# You can proceed to validate 'result_A7_0.15' and 'result_A9_0.15' by:
# - Selecting top N perturbagens by sRGES
# - Using functions or scripts for enrichment against side effects or targets
# - Performing network-based clustering of targets for functional coherence

# Example (to be developed in subsequent steps):
# Further validation workflows can include DGIdb, FAERS, or ClinicalTrials integration.

#save.image("~/OCTAD/OCTAD.RData")


# =========================================
# In Silico Validation Using topLineEval
# =========================================
# Definir rutas absolutas
source("topLineEval_custom.R")

base_dir <- "/datos/home/jjruiz/OCTAD"
output_A7 <- file.path(base_dir, "Validation_A7")
output_A9 <- file.path(base_dir, "Validation_A9")

# Crear directorios principales si no existen
if (!dir.exists(output_A7)) dir.create(output_A7, recursive = TRUE)
if (!dir.exists(output_A9)) dir.create(output_A9, recursive = TRUE)

# Validación para líneas celulares del clado A7
for (linea in lineas_similares_A7) {
  output_linea_A7 <- file.path(output_A7, linea)
  if (!dir.exists(output_linea_A7)) dir.create(output_linea_A7, recursive = TRUE)
  
  message(paste("Validando A7:", linea))
  tryCatch({
    topLineEval(
      topline = linea,
      mysRGES = result_A7,
      outputFolder = output_linea_A7
    )
  }, error = function(e) {
    message(paste("⚠️ Error con:", linea, "-", e$message))
  })
}

print(topLineEval)

# Validación para líneas celulares del clado A9
for (linea in lineas_similares_A9) {
  output_linea_A9 <- file.path(output_A9, linea)
  if (!dir.exists(output_linea_A9)) dir.create(output_linea_A9, recursive = TRUE)
  
  message(paste("Validando A9:", linea))
  tryCatch({
    topLineEval(
      topline = linea,
      mysRGES = result_A9,
      outputFolder = output_linea_A9
    )
  }, error = function(e) {
    message(paste("⚠️ Error con:", linea, "-", e$message))
  })
}

#save.image("~/OCTAD/OCTAD.RData")

