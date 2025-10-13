################################################################################
# Running PATHIFIER (Drier et al., 2013)
# Modified from Github: https://github.com/AngelCampos-Miguel Angel Garcia-Campos
################################################################################
# --- Packages ---
# msigdbr: retrieve KEGG pathways; dplyr: wrangling; pathifier: quantification
suppressPackageStartupMessages({
  library(msigdbr); library(dplyr); library(pathifier)
})
setwd("~/Pathifier-cervical-cancer/")
# --- Load expression data (genes in rows; samples in columns) ---
# The first column must contain gene identifiers (named "Gene"), which will be moved to rownames.
exp.matrix <- readRDS("Data/counts_ARSyNSeq.rds") %>% as.data.frame()

# --- Load sample metadata to derive phenotypes (Normal vs Tumor) ---
metadata <- readRDS("Data/metadata.rds")

# --- Retrieve KEGG gene sets from MSigDB ---
# Using the MEDICUS collection (current). Switch to "CP:KEGG_LEGACY" if needed.
kegg_df_all <- msigdbr(
  species     = "Homo sapiens",
  category    = "C2",
  subcategory = "CP:KEGG_MEDICUS"
) %>%
  dplyr::select(gs_name, gs_id, gene_symbol) %>%
  dplyr::distinct()

#unique_types <- unique(sub("^KEGG_MEDICUS_([^_]+).*", "\\1", kegg_df_all$gs_name))
#[1] "ENV"       "PATHOGEN"  "REFERENCE" "VARIANT"  

# Keep only KEGG_MEDICUS_REFERENCE and KEGG_MEDICUS_VARIANT
kegg_df <- kegg_df_all %>%
  dplyr::filter(grepl("KEGG_MEDICUS_REFERENCE|KEGG_MEDICUS_VARIANT", gs_name))

# Check how many remain per type
#table(sub("^KEGG_MEDICUS_([^_]+).*", "\\1", unique(kegg_df_filtered$gs_name)))
#REFERENCE   VARIANT 
#383       158 

# --- Build 'gene_sets' matrix (GMT-like): row = pathway; col1 = name; col2 = id; col3.. = genes ---
gs_list   <- split(kegg_df$gene_symbol, kegg_df$gs_name)         # pathway -> vector of genes
gs_id_map <- tapply(kegg_df$gs_id, kegg_df$gs_name, `[`, 1)      # pathway -> stable ID (first)
max_len   <- max(lengths(gs_list))                                # widest set dictates matrix width

gene_sets <- t(vapply(names(gs_list), function(pw) {
  genes <- unique(gs_list[[pw]])
  c(pw, gs_id_map[[pw]], genes, rep(NA_character_, max_len - length(genes)))
}, FUN.VALUE = character(2 + max_len)))
mode(gene_sets) <- "character"

# --- Build 'gs' list exactly as in the original script (each element = 1-column matrix of genes) ---
gs <- vector("list", length = nrow(gene_sets))
for (i in seq_len(nrow(gene_sets))) {
  a <- as.vector(gene_sets[i, 3:ncol(gene_sets)])    # gene columns only
  a <- a[!is.na(a) & a != ""]                        # drop NA/empty
  gs[[i]] <- matrix(a, ncol = 1)                     # 1-column matrix
}

# --- Pathway names and PATHWAYS container (same fields as original) ---
pathwaynames <- as.list(gene_sets[, 1])
PATHWAYS <- list(); PATHWAYS$gs <- gs; PATHWAYS$pathwaynames <- pathwaynames

# =========================
# Prepare data and parameters
# =========================

# --- Ensure genes are rownames and samples are columns (numeric matrix) ---
# Moves the "Gene" column to rownames and drops it from the data frame.
stopifnot("Gene" %in% colnames(exp.matrix))
rownames(exp.matrix) <- exp.matrix$Gene
exp.matrix <- exp.matrix[, setdiff(colnames(exp.matrix), "Gene"), drop = FALSE]
exp.matrix <- as.matrix(exp.matrix)

# --- Build the 'normals' logical phenotype from metadata ---
# TRUE = "Solid Tissue Normal"; FALSE = all other labels (e.g., "Primary Tumor").
stopifnot(all(c("specimenID","sample_type") %in% colnames(metadata)))
sample_map <- setNames(metadata$sample_type, metadata$specimenID)
common_samples <- intersect(colnames(exp.matrix), names(sample_map))
exp.matrix <- exp.matrix[, common_samples, drop = FALSE]
normals <- as.vector(sample_map[common_samples] == "Solid Tissue Normal")

# --- Calculate min_std on normal samples (25th percentile of gene-wise SD) ---
N.exp.matrix <- exp.matrix[, as.logical(normals), drop = FALSE]
if (ncol(N.exp.matrix) == 0L) stop("No normal samples found in metadata.")
rsd <- apply(N.exp.matrix, 1, sd, na.rm = TRUE)
min_std <- as.numeric(quantile(rsd, 0.25, na.rm = TRUE))

# --- Calculate min_exp globally (10th percentile of all expression values) ---
min_exp <- as.numeric(quantile(as.vector(exp.matrix), 0.10, na.rm = TRUE))

# --- Filter low-value genes and floor at min_exp ---
# Keep genes with >10% of samples above min_exp; then set values < min_exp to min_exp.
over   <- apply(exp.matrix, 1, function(x) x > min_exp)
G.over <- apply(over, 2, mean)
G.over <- names(G.over)[G.over > 0.10]
exp.matrix <- exp.matrix[G.over, , drop = FALSE]
exp.matrix[exp.matrix < min_exp] <- min_exp

# --- Keep up to 5000 most variable genes (heuristic) ---
# Sort by variance and select top 5000 (or all if fewer).
gene_var <- apply(exp.matrix, 1, var)
V <- names(sort(gene_var, decreasing = TRUE))[seq_len(min(5000L, length(gene_var)))]
V <- V[!is.na(V)]
exp.matrix <- exp.matrix[V, , drop = FALSE]
genes <- rownames(exp.matrix)
allgenes <- as.vector(rownames(exp.matrix))

# --- Pack into DATASET (same fields used by quantify_pathways_deregulation) ---
DATASET <- list(); DATASET$allgenes <- allgenes; DATASET$normals <- normals; DATASET$data <- exp.matrix

# =========================
# Light sanity checks (print-only; safe for publication)
# =========================

# Dimensions and phenotype balance
cat("DATASET dimensions (genes x samples): ", nrow(DATASET$data), " x ", ncol(DATASET$data), "\n", sep = "")
#DATASET dimensions (genes x samples): 8000 x 290
cat("Normals / Tumors: ", sum(DATASET$normals), " / ", sum(!DATASET$normals), "\n", sep = "")
#Normals / Tumors: 22 / 268

# Pathway coverage after gene filtering (size of each set post-intersection)
genes_in_sets <- unique(unlist(lapply(PATHWAYS$gs, function(m) as.vector(m[,1]))))
overlap_genes <- length(intersect(DATASET$allgenes, genes_in_sets))
cat("Unique genes in gene sets: ", length(genes_in_sets),
    " | Overlap with DATASET$allgenes: ", overlap_genes, "\n", sep = "")

set_sizes_post <- vapply(PATHWAYS$gs, function(m) sum(DATASET$allgenes %in% as.vector(m[,1])), integer(1))
cat("Pathways (total): ", length(PATHWAYS$gs),
    " | Empty (post-overlap): ", sum(set_sizes_post == 0),
    " | Median size (post-overlap): ", median(set_sizes_post), "\n", sep = "")
#Pathways (total): 658 | Empty (post-overlap): 16 | Median size (post-overlap): 8

# =========================
# Run Pathifier
# =========================

PDS <- quantify_pathways_deregulation(
  DATASET$data,
  DATASET$allgenes,
  PATHWAYS$gs,
  PATHWAYS$pathwaynames,
  DATASET$normals, 
  maximize_stability = TRUE,
  attempts = 10,
  logfile = "logfile.txt",
  min_std = min_std,
  min_exp = min_exp
)

# =========================
# Clean-up and save
# =========================

rm(gene_sets, exp.matrix, allgenes, DATASET, PATHWAYS, rsd, V, over, G.over, 
   N.exp.matrix, gs, genes, min_exp, min_std, pathwaynames)

save.image("~/Pathifier-cervical-cancer/1_1_PDS.RData")
message("DONE")
