
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


setwd("~/CESC_Network/6_OCTAD/")

# =========================================
# Load samples
# =========================================
muestras_clado_A7 <- readRDS("~/CESC_Network/6_OCTAD/6_0_A7_Samples.rds")
length(muestras_clado_A7) #[1] 66
muestras_clado_A9 <- readRDS("~/CESC_Network/6_OCTAD/6_0_A9_Samples.rds")
length(muestras_clado_A9) #[1] 202

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

# =========================================
# Controle samples
# =========================================
# === HPV samples ===
cases_union <- unique(c(case_A7, case_A9))
x <- length(cases_union)
# === Aditional controls: GTEx(normal) + TCGA(adjacent) ===
controls_all <- phenoDF %>%
  dplyr::filter(grepl("CERVIX", biopsy.site, ignore.case = TRUE),
    (data.source == "GTEX" & sample.type == "normal") |
      (data.source == "TCGA" & sample.type == "adjacent")
  ) %>%
  dplyr::pull(sample.id)


# === Top samples autoencoded ===
controls_AE_100 <- computeRefTissue(
  case_id      = cases_union,
  adjacent     = FALSE,      
  source       = "octad",    
  control_size = 100, #Suggestion top 50
  output       = TRUE,
  outputFolder = getwd()
) |> as.character()


# === Filter Famele gender ===
controls_AE <- phenoDF %>% dplyr::filter(sample.id %in% controls_AE_100) %>% dplyr::filter(gender== "Female")
dim(controls_AE) #[1] 40  15

# === Bind without dup ===
control_union <- unique(c(controls_all, controls_AE$sample.id))
length(control_union)

# === Full dataframe of control samples ===
controls_tbl <- phenoDF %>% dplyr::filter(sample.id %in% control_union)

table(controls_tbl$biopsy.site)
#   CERVIX CERVIX - ECTOCERVIX CERVIX - ENDOCERVIX  ESOPHAGUS - MUCOSA       VAGINA 
#     3                   6                   4            32                   7                                  2 

table(controls_tbl$data.source)
#GTEX TCGA 
#49    3 

saveRDS(controls_tbl, "6_1_1_Control_GTEx_TCGAadj_TopAE.rds")

# =========================
# Get tSNE of all samples used
# =========================

tsne <- octad.db::get_ExperimentHub_data("EH7276")

# IDs used (cases + controls)
used_ids <- unique(c(cases_union, controls_tbl$sample.id))

# Top-100 AE not used 
bg_ids <- setdiff(controls_AE_100, used_ids)
tsne_bg <- tsne %>%
  dplyr::filter(sample.id %in% bg_ids) %>%
  dplyr::mutate(is_top100AE = TRUE) 

# Used
ctrl_meta <- controls_tbl %>% dplyr::select(sample.id, data.source, sample.type, biopsy.site)
tsne_used <- tsne %>%
  dplyr::filter(sample.id %in% used_ids) %>%
  dplyr::left_join(ctrl_meta, by = "sample.id") %>%
  dplyr::mutate(
    type = dplyr::case_when(
      sample.id %in% cases_union ~ "Case_TCGA: CERVIX - CANCER",
      data.source == "TCGA" & sample.type == "adjacent" ~ "Control_TCGA: CERVIX - ADJACENT",
      data.source == "GTEX" & sample.type == "normal"   ~ paste0("Control_GTEx: ", biopsy.site),
      TRUE ~ "other_control"
    ),
    is_top100AE = sample.id %in% controls_AE_100
  )

# Plot: first 1 = gray BG (Top100 not used ), second 2 = used
p_tsne <- ggplot() +
  geom_point(
    data = tsne_bg,
    aes(X, Y, shape = is_top100AE),
    color = "grey40", alpha = 0.45, size = 5.2,  
    inherit.aes = FALSE
  ) +
  geom_point(
    data = tsne_used,
    aes(X, Y, color = type, shape = is_top100AE),
    alpha = 0.75, size = 5.2                   
  ) +
  scale_shape_manual(
    values = c(`FALSE` = 16, `TRUE` = 17),   
    breaks = c(TRUE),                        
    labels = c("Control_Top_100_AE"),                   
    name   = "Control_Top_100_AE") +
  guides(
    color = guide_legend(title = "", order = 1),  
    shape = guide_legend(title = "",        order = 2)     
  ) +
  labs(x = "t-SNE 1", y = "t-SNE 2",
       caption = "EH7276 precomputed t-SNE (autoencoder space)") +
  theme_bw() +
  theme(legend.title = element_text(),
        legend.text  = element_text(size = 13),
        legend.position = "right")

print(p_tsne)

ggsave(p_tsne, file="6_1_2_tSNE_samples_used.png")

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

# =========================================
# DGE  (RUVSeq + edgeR)
# =========================================
# Whole-transcriptome DE (A7 vs controls, A9 vs controls)
# RUVSeq normalization recommended when mixing GTEx/TCGA

res_A7 <- diffExp(
  case_id           = case_A7,
  control_id        = control_union,
  source            = "octad.whole",
  file              = h5_path,      
  normalize_samples = TRUE,       # RUVSeq
  k                 = 1,
  n_topGenes        = 10000,
  DE_method         = "edgeR",
  annotate          = TRUE,
  output            = FALSE
)



res_A9 <- diffExp(
  case_id           = case_A9,
  control_id        = control_union,
  source            = "octad.whole",
  file              = h5_path,
  normalize_samples = TRUE,
  k                 = 1,
  n_topGenes        = 10000,
  DE_method         = "edgeR",
  annotate          = TRUE,
  output            = FALSE
)

saveRDS(res_A7, "6_1_3_DE_full_A7.rds")
saveRDS(res_A9, "6_1_3_DE_full_A9.rds")

# =========================================
# Compact signatures (Symbol, LFC, padj, up/down)
# =========================================
pick_id_col <- function(df){
  if ("identifier" %in% names(df)) return("identifier")
  if ("ENSEMBL"   %in% names(df)) return("ENSEMBL")
  stop("No gene ID column found (identifier/ENSEMBL).")
}
pick_symbol_col <- function(df){
  if ("Symbol" %in% names(df)) return("Symbol")
  if ("Symbol_autho" %in% names(df)) return("Symbol_autho")
  NA_character_
}
make_signature <- function(res_df, lfc_thr = 1, padj_thr = 1e-3){
  sym_col <- pick_symbol_col(res_df)
  df <- if (!is.na(sym_col)) {
    res_df %>% dplyr::filter(!is.na(.data[[sym_col]])) %>%
      dplyr::mutate(Symbol = toupper(.data[[sym_col]]))
  } else {
    res_df %>% dplyr::mutate(Symbol = .data[[pick_id_col(res_df)]])  # fallback
  }
  df %>%
    dplyr::filter(abs(log2FoldChange) >= lfc_thr, padj < padj_thr) %>%
    dplyr::distinct(Symbol, .keep_all = TRUE) %>%
    dplyr::transmute(Symbol, log2FoldChange,
                     regulation = dplyr::if_else(log2FoldChange > 0, "up", "down"),
                     padj)
}

sig_A7 <- make_signature(res_A7)
sig_A9 <- make_signature(res_A9)

table(sig_A7$regulation)
#down   up 
#2109 2170 

table(sig_A9$regulation)
#down   up 
#2178 2384 

saveRDS(sig_A7, "6_1_4_DE_A7_signific.rds")
saveRDS(sig_A9, "6_1_4_DE_A9_signific.rds")

cat("Controls:", length(control_union),
    "| A7 sig:", nrow(sig_A7),
    "| A9 sig:", nrow(sig_A9), "\n")
#Controls: 52 | A7 sig: 4279 | A9 sig: 4562 


