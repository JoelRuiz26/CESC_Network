# ===================== SCRIPT OVERVIEW ===================== #
# Title: Prepro TCGA and GTEX data
# Author: Joel Ruiz Hernandez
# Contact: josejoelruizhernandez@gmail.com
# =========================================================== #
setwd("~/CESC_Network/03_preprocessing/")
# ===================== LOAD REQUIRED LIBRARIES ===================== #
library(NOISeq)
library(edgeR)
library(EDASeq)
library(tidyverse)

# ===================== SEED FOR REPRODUCIBILITY ===================== #
set.seed(123)

# ----------------- Load input data ----------------- #
factors <- readRDS("~/CESC_Network/02_cohort_data/1_Output_rds/2_2_CESC_TCGA_GTEx_metadata_clade.rds")

unstranded_counts <- readRDS("~/CESC_Network/02_cohort_data/1_Output_rds/1_1_1_Counts_CERVIX_TCGA_GTEx_raw.rds") %>% as.data.frame()

gene_metadata <- readRDS("~/CESC_Network/02_cohort_data/1_Output_rds/1_2_RowAnnot_CERVIX_GC_Length_hg38.rds")

# ----------------- Clean counts and samples ----------------- #
# Trim spaces in sample IDs
factors$sampleID <- trimws(factors$sampleID)

# Keep only samples present in both metadata and counts
valid_samples <- intersect(factors$sampleID, colnames(unstranded_counts))

factors <- factors %>%
  dplyr::filter(sampleID %in% valid_samples,
                sample_type == "Primary Tumor") %>%   # keep only primary tumor samples (A7/A9)
  as.data.frame()

counts <- unstranded_counts[, factors$sampleID, drop = FALSE]

# Remove genes with all-zero counts
counts <- counts[rowSums(counts) > 0, ]
dim(counts) #[1] 49109   280

# ----------------- Build gene annotation ----------------- #
myannot <- gene_metadata %>%
  dplyr::transmute(
    feature      = identifier,                # Ensembl gene ID (no version)
    length       = as.integer(gene_length_bp),
    gc_percent   = gc_percent,
    hgnc_symbol  = as.character(Symbol),
    gene_biotype = gene_biotype
  ) %>%
  dplyr::mutate(feature = as.character(feature)) %>%
  dplyr::distinct(feature, .keep_all = TRUE)

# If rownames(counts) include versions, uncomment:
# rownames(counts) <- sub("\\..*$", "", rownames(counts))

# ----------------- Sync genes between counts and annotation ----------------- #
common_genes <- rownames(counts)[rownames(counts) %in% myannot$feature]

counts  <- counts[common_genes, , drop = FALSE]
myannot <- myannot[match(common_genes, myannot$feature), ]

stopifnot(identical(rownames(counts), myannot$feature))

# ----------------- Sync sample metadata with counts ----------------- #
factors <- factors %>%
  dplyr::filter(sampleID %in% colnames(counts)) %>%
  as.data.frame()

counts <- counts[, factors$sampleID, drop = FALSE]
rownames(factors) <- factors$sampleID

stopifnot(all(colnames(counts) == rownames(factors)))

dim(counts)
#[1] 46731   280
dim(myannot)
#[1] 46731     5
dim(factors)
#[1] 280  14

# ----------------- Save annotation table ----------------- #
saveRDS(
  myannot,
  file = "~/CESC_Network/03_preprocessing/2_1_Universe_annotation_table.rds"
)

# =====================  CREATE NOISeq OBJECT ===================== #
# Create named vectors for annotations USING ENSG (unique IDs)
mylength <- setNames(myannot$length, myannot$feature)
mygc <- setNames(myannot$gc_percent, myannot$feature)
mybiotype <- setNames(myannot$gene_biotype, myannot$feature)

noiseqData_before <- NOISeq::readData(
  data = counts,
  factors = factors,
  gc = mygc,
  biotype = mybiotype,
  length = mylength)

saveRDS(noiseqData_before, file = "~/CESC_Network/03_preprocessing/2_2_Noiseq_before.rds")

# =====================  IDENTIFY QUALITY CONTROL ===================== #
### 3) Quality control of count data--- ---
# 3.2) Sequencing depth & Expression Quantification
# 3.2.1) Count distribution per sample  
pdf("plots_prepro_data/QC_boxplot_counts_before.pdf", width = 6, height = 6)
mycountsbio_HPV = NOISeq::dat(noiseqData_before, type = "countsbio", factor = "HPV_clade")
explo.plot(mycountsbio_HPV,
           plottype = "boxplot",
           samples = NULL)
dev.off()

# 3.2.2) Sensitivity plot: check for low count genes
pdf("plots_prepro_data/QC_barplot_sensitivity_low_counts_before.pdf", width = 8, height = 6)
explo.plot(mycountsbio_HPV,
           plottype = "barplot",
           samples = NULL)
dev.off()

# Plot global distribution of CPM 
pdf("plots_prepro_data/QC_histogram_CPM_distribution_before.pdf", width = 8, height = 6)
log_cpm_means <- rowMeans(cpm(noiseqData_before@assayData$exprs, log = TRUE))
hist_data <- hist(log_cpm_means, plot = FALSE, breaks = 30)
hist(log_cpm_means,
     ylab = "Number of genes",
     xlab = "log CPM averages",
     col = "lightblue",
     border = "darkblue",
     main = "CPM Distribution",
     breaks = 30,
     xlim = c(min(log_cpm_means), max(log_cpm_means)),
     ylim = c(0, max(hist_data$counts) * 1.1),
     xaxs = "i",
     yaxs = "i",
     cex.axis = 1.5,
     cex.lab = 1.5,
     cex.main = 2)
dev.off()

# 3.3) Sequencing bias detection
# 3.3.1) Length bias
mylengthbias_SampleType = NOISeq::dat(noiseqData_before, factor = "HPV_clade", type = "lengthbias")

pdf("plots_prepro_data/QC_length_bias_global_before.pdf", width = 8, height = 6)
par(cex.axis=1.5, cex.lab=1.5, cex.main=2)
explo.plot(mylengthbias_SampleType, samples = NULL, toplot = "global")
dev.off()

pdf("plots_prepro_data/QC_length_bias_A7_before.pdf", width = 8, height = 6)
par(cex.axis=1.5, cex.lab=1.5, cex.main=2)
explo.plot(mylengthbias_SampleType, samples = 1, toplot = "global")
dev.off()

pdf("plots_prepro_data/QC_length_bias_A9_before.pdf", width = 8, height = 6)
par(cex.axis=1.5, cex.lab=1.5, cex.main=2)
explo.plot(mylengthbias_SampleType, samples = 2, toplot = "global")
dev.off()

# 3.3.2) GC content bias
myGCbias_Sample = dat(noiseqData_before, factor = "HPV_clade", type = "GCbias")

pdf("plots_prepro_data/QC_GC_bias_global_before.pdf", width = 8, height = 6)
par(cex.axis=1.5, cex.lab=1.5, cex.main=2)
explo.plot(myGCbias_Sample, samples = NULL, toplot = "global")
dev.off()

pdf("plots_prepro_data/QC_GC_bias_A7_before.pdf", width = 8, height = 6)
par(cex.axis=1.5, cex.lab=1.5, cex.main=2)
explo.plot(myGCbias_Sample, samples = 1, toplot = "global")
dev.off()

pdf("plots_prepro_data/QC_GC_bias_A9_before.pdf", width = 8, height = 6)
par(cex.axis=1.5, cex.lab=1.5, cex.main=2)
explo.plot(myGCbias_Sample, samples = 2, toplot = "global")
dev.off()

# 3.3.3) RNA composition (slow operation)
mycd <- dat(noiseqData_before, type = "cd", norm = FALSE)
# Verify diagnostic (PASSED/FAILED across all samples)
print(table(mycd@dat$DiagnosticTest[, "Diagnostic Test"]))
pdf("plots_prepro_data/QC_Mvalues_mycd_before.pdf", width = 7, height = 6)
par(mar = c(5, 4, 4, 2) - 1)
par(cex.axis = 1.5, cex.lab = 1.5, cex.main = 2, xaxs = "i")
explo.plot(mycd, samples = 1:12)
dev.off()


# 3.4) PCA
myPCA_sample = dat(noiseqData_before, type = "PCA", norm = F, logtransf = F)
pdf("plots_prepro_data/QC_PCA_scores_HPVclade_before.pdf", width = 8, height = 6)
explo.plot(myPCA_sample, plottype= "scores", factor = "HPV_clade")
dev.off()

# ===================== FILTER FOR EDASeq QC ===================== #
#Filter low counts
#convert as integer
counts <- round(as.matrix(counts))
storage.mode(counts) <- "integer"
stopifnot(all(counts >= 0, na.rm = TRUE))

#counts_filtered <- filtered.data(
#  counts, factor = "HPV_clade", norm = FALSE, depth = NULL, method = 1, cpm = 0, p.adj = "fdr")
#13705 features are to be kept for differential expression analysis with filtering method 1
#length(unique(rownames(counts_filtered))) #[1] 13705

y       <- edgeR::DGEList(counts = round(counts),
                          group  = factors$HPV_clade)
cpm_y   <- edgeR::cpm(y)
keep    <- rowSums(cpm_y > 1) >= min(table(factors$HPV_clade))
counts_filtered  <- counts[keep, , drop = FALSE]
dim(counts_filtered) #[1] 19528   280

# ----------------- Sync gene annotation (order + rownames) ----------------- #
# Take annotation in the SAME order as counts_filtered
myannot_filtered <- myannot[match(rownames(counts_filtered), myannot$feature), ]

# convert to data.frame and set rownames = feature
myannot_filtered_df <- as.data.frame(myannot_filtered)
rownames(myannot_filtered_df) <- myannot_filtered_df$feature

# Make sure everything is aligned
stopifnot(identical(rownames(counts_filtered), rownames(myannot_filtered_df)))

# ----------------- Create EDASeq object ----------------- #
mydataEDA <- EDASeq::newSeqExpressionSet(
  counts      = counts_filtered,
  featureData = Biobase::AnnotatedDataFrame(myannot_filtered_df),
  phenoData   = Biobase::AnnotatedDataFrame(factors)
)

#order for less bias
gcFull <- withinLaneNormalization(mydataEDA, 
                                  "gc_percent", which = "full")#corrects GC bias 
lFull <- withinLaneNormalization(gcFull, "length", which = "full")#corrects length bias 
fullfullTMM <-NOISeq::tmm(normCounts(lFull), long = 1000, lc = 0, k = 0)
#norm.counts <- betweenLaneNormalization(normCounts(lFull),
# which = "median", offset = FALSE)

dim(fullfullTMM) #[1] 19528   280

# ===================== BATCH EFFECT (OMITTED) ===================== #
# ARSyNseq batch-effect correction is intentionally NOT applied:
# analysis uses only primary tumor samples (A7/A9 clades), so the
# GC + length + TMM normalized matrix (fullfullTMM) is the final matrix.

# ===================== FINAL QUALITY CHECK ===================== #
mylength_n <- setNames(myannot_filtered_df$length, myannot_filtered_df$feature)
mygcn_n <- setNames(myannot_filtered_df$gc_percent, myannot_filtered_df$feature)
mybiotype_n <- setNames(myannot_filtered_df$gene_biotype, myannot_filtered_df$feature)

noiseqData_norm_s = NOISeq::readData(data = fullfullTMM,
                                     factors=factors,
                                     gc = mygcn_n,
                                     biotype = mybiotype_n,
                                     length = mylength_n)

# 3.3.3) RNA composition
#mycd_norm_2 <- dat(noiseqData_norm_s, type = "cd", norm = TRUE)
# Verify diagnostic
#table(mycd_norm_2@dat$DiagnosticTest[, "Diagnostic Test"])
#FAILED PASSED bind #b=T source
#   36    243 

# ===================== IDENTIFY QUALITY CONTROL AFTER NORMALIZATION ===================== #
### 3) Quality control of count data--- ---
# 3.2) Sequencing depth & Expression Quantification
# 3.2.1) Count distribution per sample
pdf("plots_prepro_data/QC_boxplot_counts_after.pdf", width = 6, height = 6)
mycountsbio_HPV = NOISeq::dat(noiseqData_norm_s, type = "countsbio", factor = "HPV_clade", norm = TRUE)
explo.plot(mycountsbio_HPV,
           plottype = "boxplot",
           samples = NULL)
dev.off()

# 3.2.2) Sensitivity plot: check for low count genes
pdf("plots_prepro_data/QC_barplot_sensitivity_low_counts_after.pdf", width = 8, height = 6)
explo.plot(mycountsbio_HPV,
           plottype = "barplot",
           samples = NULL)
dev.off()

# Plot global distribution of CPM
pdf("plots_prepro_data/QC_histogram_CPM_distribution_after.pdf", width = 8, height = 6)

log_cpm_means <- rowMeans(cpm(noiseqData_norm_s@assayData$exprs, log = TRUE))
hist_data <- hist(log_cpm_means, plot = FALSE, breaks = 30)
hist(log_cpm_means,
     ylab = "Number of genes",
     xlab = "log CPM averages",
     col = "lightblue",
     border = "darkblue",
     main = "CPM Distribution",
     breaks = 30,
     xlim = c(min(log_cpm_means), max(log_cpm_means)),
     ylim = c(0, max(hist_data$counts) * 1.1),
     xaxs = "i",
     yaxs = "i",
     cex.axis = 1.5,
     cex.lab = 1.5,
     cex.main = 2)
dev.off()

# 3.3) Sequencing bias detection
# 3.3.1) Length bias
mylengthbias_SampleType = NOISeq::dat(noiseqData_norm_s, factor = "HPV_clade", type = "lengthbias", norm = TRUE)
pdf("plots_prepro_data/QC_length_bias_global_after.pdf", width = 8, height = 6)
par(cex.axis=1.5, cex.lab=1.5, cex.main=2)
explo.plot(mylengthbias_SampleType, samples = NULL, toplot = "global")
dev.off()

pdf("plots_prepro_data/QC_length_bias_A7_after.pdf", width = 8, height = 6)
par(cex.axis=1.5, cex.lab=1.5, cex.main=2)
explo.plot(mylengthbias_SampleType, samples = 1, toplot = "global")
dev.off()

pdf("plots_prepro_data/QC_length_bias_A9_after.pdf", width = 8, height = 6)
par(cex.axis=1.5, cex.lab=1.5, cex.main=2)
explo.plot(mylengthbias_SampleType, samples = 2, toplot = "global")
dev.off()

# 3.3.2) GC content bias
myGCbias_Sample = dat(noiseqData_norm_s, factor = "HPV_clade", type = "GCbias", norm = TRUE)

pdf("plots_prepro_data/QC_GC_bias_global_after.pdf", width = 8, height = 6)
par(cex.axis=1.5, cex.lab=1.5, cex.main=2)
explo.plot(myGCbias_Sample, samples = NULL, toplot = "global")
dev.off()

pdf("plots_prepro_data/QC_GC_bias_A7_after.pdf", width = 8, height = 6)
par(cex.axis=1.5, cex.lab=1.5, cex.main=2)
explo.plot(myGCbias_Sample, samples = 1, toplot = "global")
dev.off()

pdf("plots_prepro_data/QC_GC_bias_A9_after.pdf", width = 8, height = 6)
par(cex.axis=1.5, cex.lab=1.5, cex.main=2)
explo.plot(myGCbias_Sample, samples = 2, toplot = "global")
dev.off()

# 3.3.3) RNA composition (slow operation)
mycd_norm <- dat(noiseqData_norm_s, type = "cd", norm = TRUE)
# Verify diagnostic (PASSED/FAILED across all samples)
print(table(mycd_norm@dat$DiagnosticTest[, "Diagnostic Test"]))
pdf("plots_prepro_data/QC_Mvalues_mycd_after.pdf", width = 7, height = 6)
par(mar = c(5, 4, 4, 2) - 1)
par(cex.axis = 1.5, cex.lab = 1.5, cex.main = 2, xaxs = "i")
explo.plot(mycd_norm, samples = 1:12)
dev.off()

# 3.4) PCA
myPCA_sample_after = dat(noiseqData_norm_s, type = "PCA", norm = TRUE, logtransf = FALSE)
pdf("plots_prepro_data/QC_PCA_scores_HPVclade_after.pdf", width = 8, height = 6)
explo.plot(myPCA_sample_after, plottype= "scores", factor = "HPV_clade")
dev.off()

# ===================== SAVE FINAL  ===================== #
Final_SNT_vs_PP <- fullfullTMM %>% as.data.frame()
Final_SNT_vs_PP$Gene <- rownames(Final_SNT_vs_PP)  #Add column for input ARACNE
Final_SNT_vs_PP <- Final_SNT_vs_PP[, c(ncol(Final_SNT_vs_PP), 1:(ncol(Final_SNT_vs_PP)-1))]  
dim(Final_SNT_vs_PP) #[1] 19528   281

# Save full matrix
vroom::vroom_write(Final_SNT_vs_PP, file = "~/CESC_Network/03_preprocessing/2_3_Full_counts_A7_A9_notAnot.tsv", delim = "\t")
# Save metadata
vroom::vroom_write(factors, file = "~/CESC_Network/03_preprocessing/2_4_Factors.tsv", delim = "\t")
#Save annotation
vroom::vroom_write(myannot_filtered_df, file = "~/CESC_Network/03_preprocessing/2_5_Myannot.tsv", delim = "\t")

# Save matrix annoted
myannot_filtered_df <- myannot_filtered_df %>% filter(!str_detect(hgnc_symbol, "\\."))
Final_SNT_vs_PP$Symbol <- myannot_filtered_df$hgnc_symbol[
  match(Final_SNT_vs_PP$Gene, myannot_filtered_df$feature)]

#left rows without symbol
Final_SNT_vs_PP <- Final_SNT_vs_PP[!is.na(Final_SNT_vs_PP$Symbol), ]

# keep highest meadian expression 
expr_cols <- setdiff(colnames(Final_SNT_vs_PP), c("Gene", "Symbol"))

keep_rows_list <- tapply(
  1:nrow(Final_SNT_vs_PP),
  Final_SNT_vs_PP$Symbol,
  function(i) i[ which.max(rowMeans(Final_SNT_vs_PP[i, expr_cols, drop = FALSE])) ]
)

keep_rows <- sort(unlist(keep_rows_list, use.names = FALSE))

Final_SNT_vs_PP_annot <- Final_SNT_vs_PP[keep_rows, ]

# rename
Final_SNT_vs_PP_annot$Gene <- Final_SNT_vs_PP_annot$Symbol
rownames(Final_SNT_vs_PP_annot) <- Final_SNT_vs_PP_annot$Gene
Final_SNT_vs_PP_annot$Symbol <- NULL

vroom::vroom_write(Final_SNT_vs_PP_annot,
                   file = "~/CESC_Network/03_preprocessing/2_6_Full_counts_A7_A9_annot.tsv",
                   delim = "\t")

# Save A7 for network
Factors_A7 <- factors %>% dplyr::filter(HPV_clade == "A7_clade") %>% dplyr::pull(sampleID)
Final_A7_NW <- Final_SNT_vs_PP_annot %>%
  dplyr::select(Gene, dplyr::all_of(Factors_A7))

vroom::vroom_write(Final_A7_NW,
                   file = "~/CESC_Network/03_preprocessing/2_7_Full_counts_A7_annot.tsv",
                   delim = "\t")

# Save A9 for network
Factors_A9 <- factors %>% filter(HPV_clade == "A9_clade") %>% pull(sampleID)
Final_A9_NW <- Final_SNT_vs_PP_annot %>% dplyr::select(Gene, dplyr::all_of(Factors_A9))
vroom::vroom_write(Final_A9_NW,
                   file = "~/CESC_Network/03_preprocessing/2_8_Full_counts_A9_annot.tsv",
                   delim = "\t")

save.image("/STORAGE/csbig/jruiz/cesc_network/2_9_Image_Post_Prepro.RData")

#load("/STORAGE/csbig/jruiz/cesc_network/2_9_Image_Post_Prepro.RData")


# ============================================================================ #
# ===== INDEPENDENT BLOCK: GC-bias "after" with INVERTED within-lane order ==== #
# ============================================================================ #
# Run this block ON ITS OWN in a fresh R session. It is self-contained: it
# reloads the post-prepro image and only re-does the EDASeq within-lane
# normalization with the order INVERTED (length FIRST, then GC), then produces
# the "after" GC-bias plot for that alternative ordering.
# All variables are renamed with the suffix "_inv" and the PDFs use distinct
# names, so NOTHING here overwrites the objects or outputs of the main pipeline.
if (FALSE) {   # guard: skipped when sourcing the full script; run manually

  library(NOISeq)
  library(EDASeq)
  library(Biobase)

  setwd("~/CESC_Network/03_preprocessing/")
  load("/STORAGE/csbig/jruiz/cesc_network/2_9_Image_Post_Prepro.RData")

  # Invert within-lane normalization: length FIRST, then GC (main flow is GC->length)
  lFull_inv       <- withinLaneNormalization(mydataEDA, "length",     which = "full") # length bias
  gcFull_inv      <- withinLaneNormalization(lFull_inv, "gc_percent", which = "full") # GC bias
  fullfullTMM_inv <- NOISeq::tmm(normCounts(gcFull_inv), long = 1000, lc = 0, k = 0)

  # NOISeq object for QC (independent names)
  mylength_inv  <- setNames(myannot_filtered_df$length,       myannot_filtered_df$feature)
  mygc_inv      <- setNames(myannot_filtered_df$gc_percent,   myannot_filtered_df$feature)
  mybiotype_inv <- setNames(myannot_filtered_df$gene_biotype, myannot_filtered_df$feature)

  noiseqData_inv <- NOISeq::readData(data    = fullfullTMM_inv,
                                     factors = factors,
                                     gc      = mygc_inv,
                                     biotype = mybiotype_inv,
                                     length  = mylength_inv)

  # GC-bias plot AFTER only, inverted order -> distinct filenames (no overwrite)
  myGCbias_inv <- dat(noiseqData_inv, factor = "HPV_clade", type = "GCbias", norm = TRUE)

  pdf("plots_prepro_data/QC_GC_bias_global_after_lengthFirst.pdf", width = 8, height = 6)
  par(cex.axis=1.5, cex.lab=1.5, cex.main=2)
  explo.plot(myGCbias_inv, samples = NULL, toplot = "global")
  dev.off()

  pdf("plots_prepro_data/QC_GC_bias_A7_after_lengthFirst.pdf", width = 8, height = 6)
  par(cex.axis=1.5, cex.lab=1.5, cex.main=2)
  explo.plot(myGCbias_inv, samples = 1, toplot = "global")
  dev.off()

  pdf("plots_prepro_data/QC_GC_bias_A9_after_lengthFirst.pdf", width = 8, height = 6)
  par(cex.axis=1.5, cex.lab=1.5, cex.main=2)
  explo.plot(myGCbias_inv, samples = 2, toplot = "global")
  dev.off()

}
