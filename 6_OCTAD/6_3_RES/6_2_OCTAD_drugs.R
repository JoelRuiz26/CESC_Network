library(octad)
library(octad.db)
library(dplyr)
library(vroom)
library(tibble)
library(signatureSearch)

setwd("~/CESC_Network/6_OCTAD/")

# =========================================
sample_clado_A7 <- readRDS("~/CESC_Network/6_OCTAD/6_0_A7_Samples.rds")
length(sample_clado_A7) #[1] 66
sample_clado_A9 <- readRDS("~/CESC_Network/6_OCTAD/6_0_A9_Samples.rds")
length(sample_clado_A9) #[1] 202

# =========================================
# Metadata (ExperimentHub: EH7274)
# =========================================
phenoDF <- octad.db::get_ExperimentHub_data("EH7274")  # sample metadata (GTEx/TCGA)  [~19k muestras]
case_A7 <- intersect(sample_clado_A7, phenoDF$sample.id)
case_A9 <- intersect(sample_clado_A9, phenoDF$sample.id)

h5_path <- file.path(getwd(), "octad.counts.and.tpm.h5")
if (!file.exists(h5_path)) {
  url <- "https://chenlab-data-public.s3-us-west-2.amazonaws.com/octad/octad.counts.and.tpm.h5"
  message(">> Descargando HDF5 (~3 GB) a: ", h5_path)
  download.file(url, destfile = h5_path, mode = "wb", quiet = FALSE)
}
stopifnot(file.exists(h5_path))  # ensure

# =========================================
# A7 Clade 
# =========================================
##  usar whole para similitud de líneas (si tienes el HDF5)
similitud_A7 <- computeCellLine(case_id = case_A7, source = "octad.small")#, file = h5_path)
similitud_A7 <- dplyr::mutate(similitud_A7, cell_id = rownames(similitud_A7))
lineas_similares_A7 <- dplyr::filter(similitud_A7, medcor > 0.31) |> dplyr::pull(cell_id)


similitud_A9 <- computeCellLine(case_id = case_A9, source = "octad.small")#, file = h5_path)
similitud_A9 <- dplyr::mutate(similitud_A9, cell_id = rownames(similitud_A9))
lineas_similares_A9 <- dplyr::filter(similitud_A9, medcor > 0.31) |> dplyr::pull(cell_id)


# ========================================= 
# Load Differential Expression Signatures 
# ========================================= 

firma_A7 <- readRDS("6_1_4_DE_A7_signific.rds") 
firma_A9 <- readRDS("6_1_4_DE_A9_signific.rds") 

# ========================================= 
# Compute sRGES Scores (LINCS-Based Signature Reversal) 
# =========================================
## sRGES 
result_A7 <- runsRGES(
  dz_signature     = firma_A7,
  cells            = lineas_similares_A7,
  weight_cell_line = similitud_A7,
  permutations     = 10000,     
  choose_fda_drugs = TRUE,
  output           = FALSE
)

result_A9 <- runsRGES(
  dz_signature     = firma_A9,
  cells            = lineas_similares_A9,
  weight_cell_line = similitud_A9,
  permutations     = 10000,
  choose_fda_drugs = TRUE,
  output           = FALSE
)

# =========================================
# Filter Drugs with sRGES < -0.2
# =========================================
result_A7_0.2 <- result_A7 %>% dplyr::filter(sRGES <= -0.25) %>% dplyr::arrange(sRGES)
result_A9_0.2 <- result_A9 %>% dplyr::filter(sRGES <= -0.25) %>% dplyr::arrange(sRGES)
dim(result_A7_0.2); dim(result_A9_0.2)
###---  small ---###
#[1] 148   6
#[1] 161   6

###---  whole ---###
#### medcor 0.35 #### 
#[1] 91  6
#[1] 73  6

#### medcor 0.33 ####
#[1] 1280    6
#[1] 70  6

### medcor 0.31 ### 
#[1] 1269    6
#[1] 90  6

### medcor 0.3 ### 
#[1] 1267    6
#[1] 94  6

### Save data ###
saveRDS(result_A7_0.2, file = "~/CESC_Network/6_OCTAD/6_2_1_OCTAD_A7_results_0.2.rds")
saveRDS(result_A9_0.2, file = "~/CESC_Network/6_OCTAD/6_2_1_OCTAD_A9_results_0.2.rds")

# =========================================
# Enrichment Analysis of Drug Targets (CHEMBL)
# =========================================
#score= (k/K) / (n/N)
#The score is the fold of the target on your drug list with good sRGES, 
# compared to the universe of all the drugs analyzed.

# A7 enrichment
octadDrugEnrichment(
  sRGES = result_A7_0.2,
  target_type = "chembl_targets",
  enrichFolder = "6_2_2_Enriquecimiento_A7",
  outputFolder = "~/CESC_Network/6_OCTAD/"
)

enriched_chembl_A7 <- vroom("~/CESC_Network/6_OCTAD/6_2_2_Enriquecimiento_A7/chembl_targets/enriched_chembl_targets.csv")
dim(enriched_chembl_A7) #[1] 155   4

enriched_chembl_A7_filtered <- enriched_chembl_A7 %>% filter(padj < 0.05)
dim(enriched_chembl_A7_filtered)#[1] 32  4
head(enriched_chembl_A7_filtered)
# A tibble: 6 × 4
#target  score     p  padj
#<chr>   <dbl> <dbl> <dbl>
#1 CYP2C19 0.434     0     0
#2 CYP2C9  0.157     0     0
#3 CYP3A4  0.189     0     0
#4 EGFR    0.530     0     0
#5 FGFR1   0.530     0     0
#6 FLT1    0.277     0     0

# A9 enrichment
octadDrugEnrichment(
  sRGES = result_A9_0.2,
  target_type = "chembl_targets",
  enrichFolder = "6_2_2_Enriquecimiento_A9",
  outputFolder = "~/CESC_Network/6_OCTAD/"
)
enriched_chembl_A9 <- vroom("~/CESC_Network/6_OCTAD/6_2_2_Enriquecimiento_A9/chembl_targets/enriched_chembl_targets.csv")
enriched_chembl_A9_filtered <- enriched_chembl_A9 %>% filter(padj < 0.05)
dim(enriched_chembl_A9_filtered)
#[1] 26  4
head(enriched_chembl_A9_filtered)
#  A tibble: 6 × 4
#target score     p  padj
#<chr>  <dbl> <dbl> <dbl>
#1 CYP1A2 0.266     0     0
#2 EGFR   0.198     0     0
#3 FGFR1  0.438     0     0
#4 INSR   0.453     0     0
#5 ITK    0.453     0     0
#6 KDR    0.243     0     0

# =========================================
# Validation Step (According to OCTAD pipeline)
# =========================================
# OCTAD recommends validation via:
# 1. Literature mining of top-ranked compounds
# 2. Cross-referencing with clinical trial databases (e.g. ClinicalTrials.gov)
# 3. Integration with adverse effect databases (e.g. FAERS)

# - Selecting top N perturbagens by sRGES
# Further validation workflows can include DGIdb, FAERS, or ClinicalTrials integration.

# =========================================
# In Silico Validation Using topLineEval
# =========================================
# Definir rutas absolutas
source("topLineEval_custom.R")

base_dir <- "~/CESC_Network/6_OCTAD/"
output_A7 <- file.path(base_dir, "6_2_3_Validation_A7")
output_A9 <- file.path(base_dir, "6_2_3_Validation_A9")

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

save.image("~/CESC_Network/6_OCTAD/6_2_4_OCTAD.RData")
#load("~/CESC_Network/6_OCTAD/6_2_4_OCTAD.RData")


#length(unique(result_A7_0.2$pert_iname))
#[1] 91
#sum(grepl("^BRD-", result_A7_0.2$pert_iname))
#[1] 6

#length(unique(result_A9_0.2$pert_iname))
#[1] 73
#sum(grepl("^BRD-", result_A9_0.2$pert_iname))
#[1] 1
