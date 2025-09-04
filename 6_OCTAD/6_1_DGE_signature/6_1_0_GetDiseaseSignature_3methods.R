# =========================================
# Load Required Libraries
# =========================================
# https://github.com/Bin-Chen-Lab/octad
library(octad)
library(octad.db)
library(dplyr)
library(vroom)
library(tibble)
library(ggplot2)


setwd("~/CESC_Network/6_OCTAD/6_1_DGE_signature/")

# =========================================
# Load samples
# =========================================
muestras_clado_A7 <- readRDS("~/CESC_Network/6_OCTAD/6_0_A7_Samples.rds")
length(muestras_clado_A7) #[1] 66
muestras_clado_A9 <- readRDS("~/CESC_Network/6_OCTAD/6_0_A9_Samples.rds")
length(muestras_clado_A9) #[1] 202

# =========================================
# Download HDF5
# =========================================
h5_path <- file.path(getwd(), "octad.counts.and.tpm.h5")
if (!file.exists(h5_path)) {
  url <- "https://chenlab-data-public.s3-us-west-2.amazonaws.com/octad/octad.counts.and.tpm.h5"
  message(">> Descargando HDF5 (~3 GB) a: ", h5_path)
  download.file(url, destfile = h5_path, mode = "wb", quiet = FALSE)
}
stopifnot(file.exists(h5_path))  # ensure

#Get universe of genes (raw matrix before DGE)
#gene_key <- "/meta/transcripts"
#ens_ids  <- rhdf5::h5read(h5_path, gene_key)
#ens_ids  <- sub("\\.\\d+$", "", ens_ids)
# Map HGNC
#map <- AnnotationDbi::select(org.Hs.eg.db,
#                             keys = ens_ids,
#                             keytype = "ENSEMBL",
#                             columns = "SYMBOL")
#universo_hgnc <- unique(na.omit(map$SYMBOL))
#length(universo_hgnc) #[1] 36228
#saveRDS(universo_hgnc, "~/CESC_Network/6_OCTAD/6_1_DGE_signature/6_1_0_Output_rds/6_1_6_Universe_ORA.rds")

# =========================================
# Metadata (ExperimentHub: EH7274)
# =========================================
phenoDF <- octad.db::get_ExperimentHub_data("EH7274")  # sample metadata (GTEx/TCGA)  [~19k muestras]
dim(phenoDF)
#[1] 19127    15

#BiocGenerics::table(phenoDF$sample.type)
#adjacent metastatic     normal    primary  recurrent 
#737        781       7412      10033        164

# Ensure IDs existence available in OCTAD
case_A7 <- intersect(muestras_clado_A7, phenoDF$sample.id)
case_A9 <- intersect(muestras_clado_A9, phenoDF$sample.id)
if (length(case_A7) < length(muestras_clado_A7)) warning("Missing some samples in A7")
if (length(case_A9) < length(muestras_clado_A9)) warning("Missing some samples in A9")

#Missing some samples in A7 
#setdiff(muestras_clado_A7, case_A7)
#[1] "TCGA-C5-A7CM-01"

utils::str(phenoDF, list.len = 20, max.level = 2, vec.len = 8, give.attr = FALSE)


# =========================================
# Controle samples
# =========================================
# === HPV samples ===
cases_union <- unique(c(case_A7, case_A9))
x <- length(cases_union)

# === Add controls: GTEx(normal) + TCGA(adjacent) ===
controls_all <- phenoDF %>%
  dplyr::filter(grepl("CERVIX", biopsy.site, ignore.case = TRUE),
                (data.source == "GTEX" & sample.type == "normal") |
                  (data.source == "TCGA" & sample.type == "adjacent")
  ) %>%
  dplyr::pull(sample.id)


# === Bind without dup ===
control_union <- unique(c(controls_all))
length(control_union) #13

# === Full dataframe of control samples ===
controls_tbl <- phenoDF %>% dplyr::filter(sample.id %in% control_union)

table(controls_tbl$biopsy.site)
#   CERVIX CERVIX - ECTOCERVIX CERVIX - ENDOCERVIX  
#     3                   6                   4                     

table(controls_tbl$data.source)
#GTEX TCGA 
#10    3 

saveRDS(controls_tbl, "~/CESC_Network/6_OCTAD/6_1_DGE_signature/6_1_0_Output_rds/6_1_1_Control_GTEx_TCGAadj.rds")

# =========================================
# DGE with 3 methods
# =========================================
#We decide to use only the 13 OCTAD normal samples because the GTEx portal TPMs (RNA-SeQC v2, GENCODE v26) do not follow 
#the exact computational path as the Toil used by OCTAD; 
#this introduces quantification differences 
#(effective length definitions, expected count vs. TPM, annotation versions, etc.). 
#RUVSeq mitigates batching, but does not fully correct upstream pipeline discrepancies.

# Thresholds for the "significant" signature (tune as needed)
lfc_thr  <- 1.0     # absolute log2 fold change >= 1
padj_thr <- 0.01    # adjusted p-value (FDR) <= 0.01

# --- Small helpers to standardize columns across OCTAD outputs ---
pick_id_col <- function(df){
  # Prefer Ensembl/identifier column when Symbol is missing
  if ("identifier" %in% names(df)) return("identifier")
  if ("ENSEMBL"   %in% names(df)) return("ENSEMBL")
  stop("No gene ID column found (identifier/ENSEMBL).")
}
pick_symbol_col <- function(df){
  # OCTAD usually provides Symbol; fall back to ID if needed
  if ("Symbol" %in% names(df)) return("Symbol")
  if ("Symbol_autho" %in% names(df)) return("Symbol_autho")
  NA_character_
}

# Build the "ALL" (unfiltered) compact signature: Symbol, LFC, padj, up/down
make_signature_all <- function(res_df){
  sym_col <- pick_symbol_col(res_df)
  df <- if (!is.na(sym_col)) {
    res_df %>% dplyr::mutate(Symbol = toupper(.data[[sym_col]]))
  } else {
    res_df %>% dplyr::mutate(Symbol = .data[[pick_id_col(res_df)]])
  }
  df %>%
    dplyr::filter(!is.na(log2FoldChange), !is.na(padj)) %>%
    dplyr::distinct(Symbol, .keep_all = TRUE) %>%
    dplyr::transmute(
      Symbol,
      log2FoldChange,
      regulation = dplyr::if_else(log2FoldChange > 0, "up", "down"),
      padj
    )
}

# Build the "SIGNIFICANT" signature applying |LFC| and FDR thresholds
make_signature_sig <- function(res_df, lfc_thr = 1, padj_thr = 0.01){
  make_signature_all(res_df) %>%
    dplyr::filter(abs(log2FoldChange) >= lfc_thr, padj <= padj_thr)
}

# --- Minimal runner: execute one method and save results for one clade,
run_and_save <- function(clade_label, case_ids, ctrl_ids, method = c("edgeR","DESeq2","limma")){
  method <- match.arg(method)
  method_tag <- switch(method, edgeR="EdgeR", DESeq2="DESeq2", limma="limma")
  
  # In OCTAD, RUVSeq is applied only for edgeR/DESeq2 (not for limma)
  normalize_flag <- method %in% c("edgeR","DESeq2")
  
  # Differential expression via OCTAD::diffExp
  res <- diffExp(
    case_id           = case_ids,
    control_id        = ctrl_ids,
    source            = "octad.whole",
    file              = h5_path,
    normalize_samples = normalize_flag,   # RUVSeq on for edgeR/DESeq2
    k                 = if (normalize_flag) 1 else 0,
    n_topGenes        = 10000,
    DE_method         = method,
    annotate          = TRUE,
    output            = FALSE
  )
  
  # Save full DE table
  saveRDS(res, file.path("~/CESC_Network/6_OCTAD/6_1_DGE_signature/6_1_0_Output_rds/",
                         paste0("6_1_2_DE_full_", clade_label, "_", method_tag, ".rds")))
  
  # Build compact signatures (ALL + significant)
  sig_all <- make_signature_all(res)
  sig_sig <- make_signature_sig(res, lfc_thr = lfc_thr, padj_thr = padj_thr)
  
  saveRDS(sig_all, file.path("~/CESC_Network/6_OCTAD/6_1_DGE_signature/6_1_0_Output_rds/",
                             paste0("6_1_3_DE_", clade_label, "_ALL_", method_tag, ".rds")))
  saveRDS(sig_sig, file.path("~/CESC_Network/6_OCTAD/6_1_DGE_signature/6_1_0_Output_rds/",
                             paste0("6_1_4_DE_", clade_label, "_signific_", method_tag, ".rds")))
  
  # Print UP/DOWN counts
  cnt <- function(x) { tb <- table(x); up <- if ("up" %in% names(tb)) tb[["up"]] else 0L
  down <- if ("down" %in% names(tb)) tb[["down"]] else 0L
  c(up=as.integer(up), down=as.integer(down)) }
  c_all <- cnt(sig_all$regulation)
  c_sig <- cnt(sig_sig$regulation)
  
  cat(clade_label, "|", method_tag,
      "-> TOTAL:", nrow(sig_all), "(UP:", c_all["up"], "| DOWN:", c_all["down"], ")",
      "| SIGNIFIC:", nrow(sig_sig), "(UP:", c_sig["up"], "| DOWN:", c_sig["down"], ")\n")
}


# --- Run for both clades with the three methods ---
for (m in c("edgeR","DESeq2","limma")) run_and_save("A7", case_A7, control_union, m)
for (m in c("edgeR","DESeq2","limma")) run_and_save("A9", case_A9, control_union, m)
cat("Done. Signatures saved ALL and signific per method and clade\n")


### With K=1  ###
#A7
#computing DE via edgeR
#loading whole octad expression data for78samples A7 | EdgeR -> TOTAL: 17900 (UP: 9539 | DOWN: 8361 ) | SIGNIFIC: 4147 (UP: 1883 | DOWN: 2264 )
#computing DE via DESeq2
#fitting model and testingA7 | DESeq2 -> TOTAL: 17900 (UP: 9582 | DOWN: 8318 ) | SIGNIFIC: 7434 (UP: 3956 | DOWN: 3478 )
#computing DE via limma-voom
#computing DE via limmaA7 | limma -> TOTAL: 17900 (UP: 8402 | DOWN: 9498 ) | SIGNIFIC: 2579 (UP: 702 | DOWN: 1877 )

#A9
#computing DE via edgeR
#loading whole octad expression data for215samples A9 | EdgeR -> TOTAL: 19349 (UP: 10640 | DOWN: 8709 ) | SIGNIFIC: 5013 (UP: 2172 | DOWN: 2841 )
#computing DE via DESeq2
#fitting model and testingA9 | DESeq2 -> TOTAL: 19349 (UP: 10277 | DOWN: 9072 ) | SIGNIFIC: 7587 (UP: 3933 | DOWN: 3654 )
#computing DE via limma-voom
#computing DE via limmaA9 | limma -> TOTAL: 19349 (UP: 9338 | DOWN: 10011 ) | SIGNIFIC: 3386 (UP: 1004 | DOWN: 2382 )


save.image("~/CESC_Network/6_OCTAD/6_1_DGE_signature/6_1_0_Output_rds/6_1_5_Image_DiseaseSignature.RData")
#load("~/CESC_Network/6_OCTAD/6_1_DGE_signature/6_1_0_Output_rds/6_1_5_Image_DiseaseSignature.RData")




