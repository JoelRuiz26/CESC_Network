# ===================== SCRIPT OVERVIEW ===================== #
# Title: Prepro TCGA and GTEX data
# Author: Joel Ruiz Hernandez
# Contact: josejoelruizhernandez@gmail.com
# =========================================================== #
setwd("~/CESC_Network/2_Prepro_TCGA_GTEx/")
# ===================== LOAD REQUIRED LIBRARIES ===================== #
library(biomaRt)
library(NOISeq)
library(edgeR)
library(EDASeq)
library(vroom)
library(tidyverse)
library(ggplot2)
library(limma)

# ===================== SEED FOR REPRODUCIBILITY ===================== #
set.seed(123)
# ===================== LOAD INPUT DATA ===================== #
factors <- readRDS(file = "~/CESC_Network/1_Get_Data/2_Harmoni_metadata_TCGA_GTEx.rds")
dim(factors) #[1] 290   8
colnames(factors) 
#[1] "specimenID"         "source"             "cases.submitter_id"
#[4] "sample_type"        "figo_stage"         "primary_diagnosis" 
#[7] "HPV_type"           "HPV_clade"  

unstranded_counts <- readRDS(file = "~/CESC_Network/1_Get_Data/2_Harmoni_counts_TCGA_GTEx.rds")
dim(unstranded_counts) #[1] 56033   291

# ===================== CLEAN AND ANNOTATE DATA ===================== #
# Clean counts matrix
counts <- as.data.frame(unstranded_counts) %>%
  mutate(gene = str_remove(gene, "\\..*$"))  # remove Ensembl version if present

# Remove genes with all zero counts
counts <- counts[rowSums(counts[,-1]) > 0, ]
dim(counts) #[1] 56012   291

# Get annotation from Ensembl (GENCODE v43 ~ Ensembl v110 default)
mart <- useEnsembl("ensembl", dataset="hsapiens_gene_ensembl")
myannot <- getBM(attributes = c("ensembl_gene_id", 
                                "percentage_gene_gc_content", "gene_biotype",
                                "start_position", "end_position", "hgnc_id", "hgnc_symbol"),
                 filters = "ensembl_gene_id", 
                 values = counts$gene,
                 mart = mart) %>%
  rename(feature = ensembl_gene_id) %>%
  mutate(length = abs(end_position - start_position)) %>%
  filter(hgnc_id != "" & hgnc_symbol != "") %>%
  distinct(feature, .keep_all = TRUE)

# Filter counts to those with valid annotations
counts <- counts %>% filter(gene %in% myannot$feature)
dim(counts) #[1] 39750   291

# Reorder myannot to match counts order
myannot <- myannot %>%
  filter(feature %in% counts$gene) %>%
  arrange(match(feature, counts$gene))
dim(myannot) #[1] 39750     8

#Set rownames as ENSG for analysis
rownames(counts) <- counts$gene
counts <- counts %>% dplyr::select(-gene)
dim(counts) #[1] 39750   290

# ========== SYNC DATA ========== #
factors <- factors %>% filter(specimenID %in% colnames(counts)) %>% as.data.frame()
counts  <- counts[, factors$specimenID, drop = FALSE] 
rownames(factors) <- factors$specimenID
#Verify
stopifnot(all(colnames(counts) == rownames(factors)))
dim(counts) #[1] 39750   290
dim(factors) #[1] 290   8
table(factors$HPV_clade)
#A7_clade            A9_clade Solid Tissue Normal 
#66                 202                  22 

# =====================  CREATE NOISeq OBJECT ===================== #
# Create named vectors for annotations USING ENSG (unique IDs)
mylength <- setNames(myannot$length, myannot$feature)
mygc <- setNames(myannot$percentage_gene_gc_content, myannot$feature)
mybiotype <- setNames(myannot$gene_biotype, myannot$feature)

noiseqData_before <- NOISeq::readData(
  data = counts,
  factors = factors,
  gc = mygc,
  biotype = mybiotype,
  length = mylength)

#saveRDS(noiseqData_before, file = "~/CESC_Network/2_Prepro_TCGA_GTEx/2_2_Noiseq_before.rds")

# =====================  IDENTIFY QUALITY CONTROL ===================== #
### 3) Quality control of count data--- ---
# 3.2) Sequencing depth & Expression Quantification
# 3.2.1) Count distribution per sample  
pdf("plots_prepro_data/QC_boxplot_counts_before.pdf", width = 8, height = 6)
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

pdf("plots_prepro_data/QC_length_bias_SNT.pdf", width = 8, height = 6)
par(cex.axis=1.5, cex.lab=1.5, cex.main=2)
explo.plot(mylengthbias_SampleType, samples = 3, toplot = "global")
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

pdf("plots_prepro_data/QC_GC_bias_SNT_before.pdf", width = 8, height = 6)
par(cex.axis=1.5, cex.lab=1.5, cex.main=2)
explo.plot(myGCbias_Sample, samples = 3, toplot = "global")
dev.off()

# 3.3.3) RNA composition
mycd <- dat(noiseqData_before, type = "cd", norm = FALSE) # Slow operation
pdf("plots_prepro_data/QC_Mvalues_mycd_before.pdf", width = 10, height = 6)
par(mar = c(5, 4, 4, 2) - 1)
par(cex.axis = 1.5, cex.lab = 1.5, cex.main = 2, xaxs = "i")
explo.plot(mycd, samples = c(266:276))
dev.off()
# Verify diagnostic
table(mycd@dat$DiagnosticTest[, "Diagnostic Test"])
#FAILED PASSED 
#288      1 

# 3.4) PCA
myPCA_sample = dat(noiseqData_before, type = "PCA", norm = F, logtransf = F)
pdf("plots_prepro_data/QC_PCA_scores_HPVclade_before.pdf", width = 8, height = 6)
explo.plot(myPCA_sample, plottype= "scores", factor = "HPV_clade")
dev.off()

pdf("plots_prepro_data/QC_PCA_scores_source_before.pdf", width = 8, height = 6)
explo.plot(myPCA_sample, plottype= "scores", factor = "source")
dev.off()


# ===================== FILTER FOR EDASeq QC ===================== #
#Filter low counts
counts_filtered <- filtered.data(
  counts, factor = "HPV_clade", norm = FALSE, depth = NULL, method = 1, cpm = 0, p.adj = "fdr")
#11544 features are to be kept for differential expression analysis with filtering method 1
length(unique(rownames(counts_filtered))) #[1] 11544
#sync
myannot_filtered <- myannot[myannot$feature %in% rownames(counts_filtered), ]
stopifnot(all(rownames(counts_filtered) %in% myannot_filtered$feature))

# Set rownames 
rownames(myannot_filtered) <- myannot_filtered$feature
factors <- factors[colnames(counts_filtered), , drop = FALSE]
rownames(factors) <- factors$specimenID
# Ensure full synchrony
stopifnot(all(rownames(factors) == colnames(counts_filtered)))

# Create EDASeq object
mydataEDA <- newSeqExpressionSet(
  counts = as.matrix(counts_filtered),
  featureData = AnnotatedDataFrame(data = myannot_filtered),
  phenoData = AnnotatedDataFrame(data = factors))

#order for less bias
lFull <- withinLaneNormalization(mydataEDA, "length", which = "full")#corrects length bias 
gcFull <- withinLaneNormalization(lFull, 
                                  "percentage_gene_gc_content", which = "full")#corrects GC bias 
fullfullTMM <-NOISeq::tmm(normCounts(gcFull), long = 1000, lc = 0, k = 0)
#norm.counts <- betweenLaneNormalization(normCounts(lFull),
# which = "median", offset = FALSE)

noiseqData_filtered = NOISeq::readData(data = fullfullTMM, 
                                   factors=factors)

#myPCA_pre = dat(noiseqData_filtered, type = "PCA", norm = T,logtransf = F)
#explo.plot(myPCA_sample, plottype= "scores", factor = "HPV_clade")
#explo.plot(myPCA_sample, plottype= "scores", factor = "source")
#dev.off()

mycd_filtered <- dat(noiseqData_filtered, type = "cd", norm = TRUE ) # Slow operation
# Verify diagnostic
table(myc_filteredd@dat$DiagnosticTest[, "Diagnostic Test"])
#FAILED PASSED 
#34    255

# ===================== SOLVE BATCH EFFECT ===================== #
#Expression counts were preprocessed by GC and length bias correction (EDASeq) 
#and TMM normalization (edgeR). No batch correction was applied at this stage, 
#to retain covariance structure for network inference. 
#Batch (source) was subsequently modeled as a covariate in the 
#differential expression analysis using limma-trend, 
#to adjust for cross-dataset differences without removing biological signal.


#“For variables highly confounded with biological groups,
#batch correction should be done at the model stage (as a covariate), 
#not as data correction beforehand.” (Leek et al., 2010; Law et al., 2014)

#Batch correction with limma::removeBatchEffect was applied to reduce source-related bias 
#while adjusting for sample_type. Given the confounding between batch and biology, 
#biological contrasts may be partially attenuated.

# log2-transform TMM-normalized counts
expr_logCPM <- log2(fullfullTMM + 0.5)

# Correct batch using limma::removeBatchEffect
# We adjust batch (source) while preserving sample_type biological signal
expr_batchcorrected <- removeBatchEffect(
  expr_logCPM,
  batch = factors$source,
  design = model.matrix(~ sample_type, data = factors)
)

#ffTMMARSyn=ARSyNseq(noiseqData_filtered, factor = "HPV_clade", batch = F,
#                    norm = "n",  logtransf = F)

noiseqData_batch = NOISeq::readData(data = expr_batchcorrected, 
                                       factors=factors)

myPCA_after = dat(noiseqData_batch  , type = "PCA", norm = T,logtransf = T)


pdf("plots_prepro_data/QC_PCA_scores_HPVclade_after.pdf", width = 8, height = 6)
explo.plot(myPCA_after, plottype = "scores", factor = "source")
dev.off()

pdf("plots_prepro_data/QC_PCA_scores_source_after.pdf", width = 8, height = 6)
explo.plot(myPCA_after, plottype = "scores", factor = "HPV_clade")
dev.off()

# ===================== FINAL QUALITY CHECK ===================== #
mylength_n <- setNames(myannot_filtered$length, myannot_filtered$feature)
mygcn_n <- setNames(myannot_filtered$percentage_gene_gc_content, myannot_filtered$feature)
mybiotype_n <- setNames(myannot_filtered$gene_biotype, myannot_filtered$feature)

noiseqData_norm = NOISeq::readData(data = expr_batchcorrected, 
                                   factors=factors,
                                   gc = mygcn_n,
                                   biotype = mybiotype_n,
                                   length = mylength_n)

#  _2 = second ARSyNseq (factor = source, batch= T, + factor = HPV_clade, batch = F) 
#ffTMMARSyn_norm=ARSyNseq(noiseqData_norm, factor = "HPV_clade", batch = F,
#                    norm = "n",  logtransf = F)
#noiseqData_norm_2 = NOISeq::readData(data = exprs(ffTMMARSyn_norm), 
#                                   factors=factors,
#                                   gc = mygcn_n,
#                                   biotype = mybiotype_n,
#                                   length = mylength_n)
#myPCA_after_2 = dat(ffTMMARSyn_norm, type = "PCA", norm = T,logtransf = F)
#pdf("plots_prepro_data/QC_PCA_scores_HPVclade_after_2.pdf", width = 8, height = 6)
#explo.plot(myPCA_after_2, plottype = "scores", factor = "source")
#dev.off()
#pdf("plots_prepro_data/QC_PCA_scores_source_after_2.pdf", width = 8, height = 6)
#explo.plot(myPCA_after_2, plottype = "scores", factor = "HPV_clade")
#dev.off()
# 3.3.3) RNA composition
#mycd_norm_2 <- dat(noiseqData_norm_2, type = "cd", norm = TRUE)
# Verify diagnostic
#table(mycd_norm_2@dat$DiagnosticTest[, "Diagnostic Test"])
#FAILED PASSED bind #b=T source and #b=F HPV_clade
#51    238 

# ===================== IDENTIFY QUALITY CONTROL AFTER NORMALIZATION ===================== #
### 3) Quality control of count data--- ---
# 3.2) Sequencing depth & Expression Quantification
# 3.2.1) Count distribution per sample
pdf("plots_prepro_data/QC_boxplot_counts_after.pdf", width = 8, height = 6)
mycountsbio_HPV = NOISeq::dat(noiseqData_norm, type = "countsbio", factor = "HPV_clade", norm = TRUE)
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
log_cpm_means <- rowMeans(cpm(noiseqData_norm@assayData$exprs, log = TRUE))
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
mylengthbias_SampleType = NOISeq::dat(noiseqData_norm, factor = "HPV_clade", type = "lengthbias", norm = TRUE)
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

pdf("plots_prepro_data/QC_length_bias_SNT_after.pdf", width = 8, height = 6)
par(cex.axis=1.5, cex.lab=1.5, cex.main=2)
explo.plot(mylengthbias_SampleType, samples = 3, toplot = "global")
dev.off()

# 3.3.2) GC content bias
myGCbias_Sample = dat(noiseqData_norm, factor = "HPV_clade", type = "GCbias", norm = TRUE)

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

pdf("plots_prepro_data/QC_GC_bias_SNT_after.pdf", width = 8, height = 6)
par(cex.axis=1.5, cex.lab=1.5, cex.main=2)
explo.plot(myGCbias_Sample, samples = 3, toplot = "global")
dev.off()

# 3.3.3) RNA composition
mycd_norm<- dat(noiseqData_norm, type = "cd", norm = TRUE)
# Verify diagnostic
table(mycd_norm@dat$DiagnosticTest[, "Diagnostic Test"])
#FAILED PASSED #source b=T  
#25    264 
#FAILED PASSED #b=F HPV_clade 
#68    221 

#FAILED PASSED    limma corrected (I've used this)
# 67    222 

pdf("plots_prepro_data/QC_Mvalues_mycd_after.pdf", width = 10, height = 6)
par(mar = c(5, 4, 4, 2) - 1)
par(cex.axis = 1.5, cex.lab = 1.5, cex.main = 2, xaxs = "i")
explo.plot(mycd_norm_2, samples = c(266:276))
dev.off()

# 3.4) PCA
myPCA_sample_after = dat(noiseqData_norm, type = "PCA", norm = TRUE, logtransf = FALSE)
pdf("plots_prepro_data/QC_PCA_scores_HPVclade_after.pdf", width = 8, height = 6)
explo.plot(myPCA_sample_after, plottype= "scores", factor = "HPV_clade")
dev.off()

pdf("plots_prepro_data/QC_PCA_scores_source_after.pdf", width = 8, height = 6)
explo.plot(myPCA_sample_after, plottype= "scores", factor = "source")
dev.off()

# ===================== SAVE FINAL  ===================== #
Final_SNT_vs_PP <- expr_batchcorrected %>% as.data.frame()
Final_SNT_vs_PP$Gene <- rownames(Final_SNT_vs_PP)  #Add column for input ARACNE
Final_SNT_vs_PP <- Final_SNT_vs_PP[, c(ncol(Final_SNT_vs_PP), 1:(ncol(Final_SNT_vs_PP)-1))]  
dim(Final_SNT_vs_PP) #[1] 11544   291

# Save full matrix
vroom::vroom_write(Final_SNT_vs_PP, file = "~/CESC_Network/2_Prepro_TCGA_GTEx/2_3_Full_counts_A7_A9_notAnot_11544.tsv", delim = "\t")
# Save metadata
vroom::vroom_write(factors, file = "~/CESC_Network/2_Prepro_TCGA_GTEx/2_4_Factors.tsv", delim = "\t")
#Save annotation
vroom::vroom_write(myannot_filtered, file = "~/CESC_Network/2_Prepro_TCGA_GTEx/2_5_Myannot.tsv", delim = "\t")

# Matrix annoted
Final_SNT_vs_PP$HGNC_symbol <- myannot_filtered$hgnc_symbol[match(Final_SNT_vs_PP$Gene, myannot_filtered$feature)]
Final_SNT_vs_PP_annot <- Final_SNT_vs_PP %>% filter(!is.na(HGNC_symbol) & HGNC_symbol != "")
Final_SNT_vs_PP_annot <- Final_SNT_vs_PP_annot %>%
  dplyr::select(HGNC_symbol, everything(), -Gene) %>%
  dplyr::rename(Gene = HGNC_symbol)
rownames(Final_SNT_vs_PP_annot) <- Final_SNT_vs_PP_annot$Gene

# Save matrix annoted
vroom::vroom_write(Final_SNT_vs_PP_annot,
                   file = "~/CESC_Network/2_Prepro_TCGA_GTEx/2_6_Full_counts_A7_A9_annot.tsv",
                   delim = "\t")

# Save A7 for network
Factors_A7 <- factors %>% filter(HPV_clade == "A7_clade")
Final_A7_NW <- Final_SNT_vs_PP_annot %>%
  dplyr::select(Gene, all_of(Factors_A7$specimenID))
vroom::vroom_write(Final_A7_NW, 
                   file = "~/CESC_Network/2_Prepro_TCGA_GTEx/2_7_Full_counts_A7_annot.tsv", 
                   delim = "\t")

# Save A9 for network
Factors_A9 <- factors %>% filter(HPV_clade == "A9_clade")
Final_A9_NW <- Final_SNT_vs_PP_annot %>%
  dplyr::select(Gene, all_of(Factors_A9$specimenID))
vroom::vroom_write(Final_A9_NW, 
                   file = "~/CESC_Network/2_Prepro_TCGA_GTEx/2_8_Full_counts_A9_annot.tsv", 
                   delim = "\t")

#save.image("~/CESC_Network/2_Prepro_TCGA_GTEx/2_9_Image_Post_Prepro.RData")
#load("~/CESC_Network/2_Prepro_TCGA_GTEx/2_9_Image_Post_Prepro.RData")
