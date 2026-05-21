############################################################
# SCRIPT: 4_1_DE_limma_A7_A9_vs_SNT.R
#
# Purpose:
#   Differential gene expression analysis using limma
#   comparing HPV clade–specific tumors vs solid normal tissue,
#   and between HPV clades (A7 vs A9).
#
# Input (preprocessed, batch-corrected data):
#   - 2_6_Full_counts_A7_A9_annot.tsv  (TMM + ARSyN, annotated, genes x samples)
#   - 2_4_Factors.tsv                  (sample metadata with HPV_clade, source, etc.)
#
# Output:
#   - 4_1_limma_A7_vs_SNT_allGenes.tsv
#   - 4_2_limma_A9_vs_SNT_allGenes.tsv
#   - 4_3_limma_A7_vs_A9_allGenes.tsv   
#   - 4_4_limma_DE_results_list.rds
#
# Notes:
#   - Expression matrix contains TMM-normalized, ARSyN-corrected counts
#     in linear scale (not log-transformed).
#   - We apply a log2(x + 0.5) transform and use limma-trend.
############################################################


setwd("~/CESC_Network/2_Prepro_TCGA_GTEx/")

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(vroom)
  library(limma)
})

############################################################
## 1) Load preprocessed data
############################################################

# Expression matrix: rows = genes (symbols), columns = samples
# First column is "Gene" with gene symbols
expr_df <- readr::read_tsv(
  "2_6_Full_counts_A7_A9_annot.tsv",
  col_types = cols())

# Sample metadata: includes at least sampleID and HPV_clade
factors <- readr::read_tsv(
  "2_4_Factors.tsv",
  col_types = cols())

############################################################
## 2) Prepare expression matrix
############################################################

# Sanity check: make sure "Gene" column exists
stopifnot("Gene" %in% colnames(expr_df))

# Extract gene identifiers (gene symbols) and remove column from matrix
gene_ids <- expr_df$Gene

expr_mat <- expr_df %>%
  dplyr::select(-Gene) %>%
  as.matrix()

# Assign gene symbols as rownames
rownames(expr_mat) <- gene_ids

# Check that matrix is numeric
stopifnot(is.numeric(expr_mat[1, 1]))

############################################################
## 3) Align samples between expression matrix and metadata
############################################################

# Keep only samples present in both expr_mat and factors
common_samples <- intersect(colnames(expr_mat), factors$sampleID)
message("Number of samples in common between matrix and metadata: ",
        length(common_samples))

# Subset expression matrix and metadata to common samples
expr_mat   <- expr_mat[, common_samples, drop = FALSE]
factors_al <- factors %>%
  dplyr::filter(sampleID %in% common_samples) %>%
  as.data.frame()

# Reorder metadata rows to match column order of expression matrix
factors_al <- factors_al[match(colnames(expr_mat), factors_al$sampleID), ]

# Final consistency check: columns in expr_mat match rows in factors_al
stopifnot(all(colnames(expr_mat) == factors_al$sampleID))

############################################################
## 4) Define groups for limma design matrix
############################################################

# HPV_clade should contain (at least):
#   - "A7_clade"
#   - "A9_clade"
#   - "Solid Tissue Normal" (controls)
print(table(factors_al$HPV_clade))
#A7_clade            A9_clade Solid Tissue Normal 
#65                 202                  13


# Encode HPV_clade as a factor and set "Solid Tissue Normal" as reference
factors_al$HPV_clade <- factor(factors_al$HPV_clade)
factors_al$HPV_clade <- relevel(factors_al$HPV_clade, ref = "Solid Tissue Normal")

# Build design matrix without intercept (one column per level)
design <- model.matrix(~ 0 + HPV_clade, data = factors_al)

# Remove prefix "HPV_clade" from column names
colnames(design) <- gsub("HPV_clade", "", colnames(design))

# Rename the normal-tissue column to a syntactically valid name ("Normal")
colnames(design)[colnames(design) == "Solid Tissue Normal"] <- "Normal"

# Design matrix should now have columns: "Normal", "A7_clade", "A9_clade"
print(colnames(design))
#[1] "Normal"   "A7_clade" "A9_clade"


############################################################
## 5) Log-transform preprocessed counts for limma-trend
############################################################

# expr_mat contains TMM-normalized, ARSyN-corrected counts in linear scale.
# For limma-trend, we work on a log2 scale to stabilize the variance.
# A small offset (0.5) avoids log(0) for very low counts.
log_expr <- log2(expr_mat + 0.5)

############################################################
## 6) Fit limma linear model and define contrasts
############################################################

# Fit linear model for each gene
fit <- limma::lmFit(log_expr, design)

# Define contrasts:
#   - A7_vs_SNT: A7_clade tumors vs Normal (solid tissue normal)
#   - A9_vs_SNT: A9_clade tumors vs Normal (solid tissue normal)
#   - A7_vs_A9:  A7_clade tumors vs A9_clade tumors (NEW)
#
# Note for A7_vs_A9: a positive logFC means higher expression in A7,
# and a negative logFC means higher expression in A9.
contrast_matrix <- limma::makeContrasts(
  A7_vs_SNT = A7_clade - Normal,
  A9_vs_SNT = A9_clade - Normal,
  A7_vs_A9  = A7_clade - A9_clade,
  levels = design)

# Apply contrasts and empirical Bayes moderation (with trend)
fit2 <- limma::contrasts.fit(fit, contrast_matrix)
fit2 <- limma::eBayes(fit2, trend = TRUE, robust = TRUE)

############################################################
## 7) Extract full differential expression tables
############################################################

# We extract all genes (no thresholds here) and sort by p-value.
# Filtering by adj.P.Val and |logFC| can be done downstream.

DE_A7_vs_SNT <- limma::topTable(
  fit2,
  coef    = "A7_vs_SNT",
  number  = Inf,
  sort.by = "P")

DE_A9_vs_SNT <- limma::topTable(
  fit2,
  coef    = "A9_vs_SNT",
  number  = Inf,
  sort.by = "P")

# NEW: Extract DE results for the A7 vs A9 tumor contrast
DE_A7_vs_A9 <- limma::topTable(
  fit2,
  coef    = "A7_vs_A9",
  number  = Inf,
  sort.by = "P")

# Add explicit gene column (rownames contain gene symbols)
DE_A7_vs_SNT$Gene <- rownames(DE_A7_vs_SNT)
DE_A9_vs_SNT$Gene <- rownames(DE_A9_vs_SNT)
DE_A7_vs_A9$Gene  <- rownames(DE_A7_vs_A9)

# Move Gene column to the first position
DE_A7_vs_SNT <- DE_A7_vs_SNT %>%
  dplyr::relocate(Gene, .before = dplyr::everything())

DE_A9_vs_SNT <- DE_A9_vs_SNT %>%
  dplyr::relocate(Gene, .before = dplyr::everything())

DE_A7_vs_A9 <- DE_A7_vs_A9 %>%
  dplyr::relocate(Gene, .before = dplyr::everything())

############################################################
## 8) Save results (all final annotated genes)
############################################################

# Save full limma tables as TSV files
vroom::vroom_write(
  DE_A7_vs_SNT,
  file  = "~/CESC_Network/3_DGE/3_1_limma_A7_vs_SNT_allGenes.tsv",
  delim = "\t")

vroom::vroom_write(
  DE_A9_vs_SNT,
  file  = "~/CESC_Network/3_DGE/3_1_limma_A9_vs_SNT_allGenes.tsv",
  delim = "\t")

# NEW: Save A7 vs A9 tumor contrast
vroom::vroom_write(
  DE_A7_vs_A9,
  file  = "~/CESC_Network/3_DGE/3_1_limma_A7_vs_A9_allGenes.tsv",
  delim = "\t")

############################################################
## 9) Explore results
############################################################

# Function to count DE genes
count_DE <- function(df, name = "Contrast", padj = 0.01, lfc = 1) {
  
  df_sig  <- df %>% filter(adj.P.Val <= padj, abs(logFC) >= lfc)
  df_up   <- df_sig %>% filter(logFC > 0)
  df_down <- df_sig %>% filter(logFC < 0)
  
  tibble(
    Contrast      = name,
    Significant   = nrow(df_sig),
    Upregulated   = nrow(df_up),
    Downregulated = nrow(df_down)
  )
}

# Apply to A7 vs SNT, A9 vs SNT, and A7 vs A9
res_A7     <- count_DE(DE_A7_vs_SNT, "A7_vs_SNT")
res_A9     <- count_DE(DE_A9_vs_SNT, "A9_vs_SNT")
res_A7_A9  <- count_DE(DE_A7_vs_A9,  "A7_vs_A9")

# Print results
res_A7
res_A9
res_A7_A9


#> res_A7
# A tibble: 1 × 4
#Contrast  Significant Upregulated Downregulated
#<    chr>           <int>       <int>         <int>
#  1 A7_vs_SNT        7033        3506          3527
#> res_A9
# A tibble: 1 × 4
#Contrast  Significant Upregulated Downregulated
#     <chr>           <int>       <int>         <int>
#  1 A9_vs_SNT        7123        3540          3583
#> res_A7_A9
# A tibble: 1 × 4
#Contrast Significant Upregulated Downregulated
#     <chr>          <int>       <int>         <int>
#  1 A7_vs_A9         996         474           522
