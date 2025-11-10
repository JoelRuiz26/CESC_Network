# ===================== SCRIPT OVERVIEW ===================== #
# Title: Prepro TCGA and GTEX data
# Author: Joel Ruiz Hernandez
# Contact: josejoelruizhernandez@gmail.com
# =========================================================== #
setwd("~/CESC_Network/2_Prepro_TCGA_GTEx/")
# ===================== LOAD REQUIRED LIBRARIES ===================== #
library(NOISeq)
library(edgeR)
library(EDASeq)
library(vroom)
library(tidyverse)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg38)

# ===================== SEED FOR REPRODUCIBILITY ===================== #
set.seed(123)
# ===================== LOAD INPUT DATA ===================== #
factors <- readRDS(file = "~/CESC_Network/1_Get_Data/1_2_1_Metadata_full_TCGA_GTEX.rds")
dim(factors) #[1] 281   9
colnames(factors) 
#[1] "specimenID"        "source"            "case_id_12"       
#[4] "sample_type"       "HPV_type"          "HPV_clade"        
#[7] "figo_stage"        "primary_diagnosis" "race" 

unstranded_counts <- readRDS(file = "~/CESC_Network/1_Get_Data/1_2_2_Cervix_Toil_raw_counts.rds")
dim(unstranded_counts) #[1] 60498   281

# ===================== CLEAN AND ANNOTATE DATA ===================== #
# Clean counts matrix
counts <- unstranded_counts %>%
  dplyr::mutate(gene = str_remove(gene, "\\..*$"))  # remove Ensembl version if present

# Remove genes with all zero counts
counts <- counts[rowSums(counts[,-1]) > 0, ]
dim(counts) #[1] 49109   281

# ===================== ANNOTATION FROM GENCODE v23 (GRCh38) ===================== #
# GTF v23 (alineado a Toil hg38). 
#  dir.create("/STORAGE/csbig/jruiz/cesc_network/refs/", showWarnings = FALSE, recursive = TRUE)
#download.file(
#  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_23/gencode.v23.annotation.gtf.gz",
#  destfile = "/STORAGE/csbig/jruiz/cesc_network/refs/gencode.v23.annotation.gtf.gz", mode = "wb"
#)

gtf <- "/STORAGE/csbig/jruiz/cesc_network/refs/gencode.v23.annotation.gtf.gz"

# --- TxDb y exones por gen ---
txdb <- makeTxDbFromGFF(gtf, format = "gtf")
# UCSC para compatibilidad con BSgenome
try({ GenomeInfoDb::seqlevelsStyle(txdb) <- "UCSC" }, silent = TRUE)

ex_by_gene_all <- exonsBy(txdb, by = "gene")

# --- ENSG sin versión y filtrado al universo de 'counts' ---
gene_ids_orig   <- names(ex_by_gene_all)                 # e.g. ENSG000001.12
gene_ids_novers <- sub("\\..*$", "", gene_ids_orig)      # ENSG000001
keep_idx        <- gene_ids_novers %in% counts$gene

ex_by_gene  <- ex_by_gene_all[keep_idx]
feature_ids <- gene_ids_novers[keep_idx]

# --- Genoma de referencia ---
hg38 <- BSgenome.Hsapiens.UCSC.hg38

# --- Función: longitud exónica no redundante + %GC ---
calc_len_gc <- function(feature_id, gr) {
  # Por seguridad, estilo UCSC
  try({ GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC" }, silent = TRUE)
  rr <- GenomicRanges::reduce(gr, ignore.strand = TRUE)
  if (length(rr) == 0L) {
    return(tibble(feature = feature_id, length = 0L, gc_percent = NA_real_))
  }
  seqs <- getSeq(hg38, rr)
  total_len <- sum(width(rr))
  lf <- letterFrequency(seqs, letters = c("G", "C"), as.prob = FALSE)
  gc_count <- sum(lf[, "G"] + lf[, "C"])
  tibble(
    feature    = feature_id,
    length     = as.integer(total_len),
    gc_percent = ifelse(total_len > 0, (gc_count / total_len) * 100, NA_real_)
  )
}

# --- Cálculo en bloque ---
ex_by_gene_list <- as.list(ex_by_gene)
annot_len_gc <- purrr::map2_dfr(feature_ids, ex_by_gene_list, calc_len_gc)

# ===================== SYMBOL & BIOTYPE FROM THE SAME GTF ===================== #
gtf_gr <- import(gtf)
try({ GenomeInfoDb::seqlevelsStyle(gtf_gr) <- "UCSC" }, silent = TRUE)

gene_rows <- gtf_gr[gtf_gr$type == "gene"]
gene_meta <- tibble(
  feature_with_ver = mcols(gene_rows)$gene_id,   # ENSG con versión
  hgnc_symbol      = mcols(gene_rows)$gene_name,
  gene_biotype     = mcols(gene_rows)$gene_type
) %>%
  mutate(feature = sub("\\..*$", "", feature_with_ver)) %>%  # ENSG sin versión
  distinct(feature, .keep_all = TRUE) %>%
  dplyr::select(feature, hgnc_symbol, gene_biotype)

# ===================== JOIN: LENGTH/GC + SYMBOL/BIOTYPE ===================== #
myannot <- annot_len_gc %>%
  dplyr::left_join(gene_meta, by = "feature") %>%
  dplyr::distinct(feature, .keep_all = TRUE)
dim(myannot)
#[1] 49109     5
length(unique(myannot$hgnc_symbol))
#[1] 48505

# ===================== SYNC WITH COUNTS & FINALIZE ===================== #
# Filtra counts a genes con anotación válida
counts <- counts %>% dplyr::filter(gene %in% myannot$feature)
dim(counts)
#[1] 49109   281

# Reordena anotación para que coincida con el orden de 'counts'
myannot <- myannot %>%
  dplyr::filter(feature %in% counts$gene) %>%
  dplyr::arrange(match(feature, counts$gene))

# Rownames = ENSG y elimina columna 'gene' en matriz
rownames(counts) <- counts$gene
counts <- counts %>% dplyr::select(-gene)

# Dimensiones finales (aprox. ~35k genes)
dim(counts)
#[1] 49109   280
dim(myannot)
#[1] 49109     5

# Opcional: guardar para ORA/otros análisis
#saveRDS(myannot, file = "~/CESC_Network/2_Prepro_TCGA_GTEx/2_1_Universe_annotation_table.rds")

# ========== SYNC DATA ========== #
factors$specimenID <- sub("^(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}).*$", "\\1", factors$specimenID)
factors <- factors %>% filter(specimenID %in% colnames(counts)) %>% as.data.frame()
counts  <- counts[, factors$specimenID, drop = FALSE] 
rownames(factors) <- factors$specimenID
#Verify
stopifnot(all(colnames(counts) == rownames(factors)))
dim(counts) #[1] 34577   280
dim(factors) #[1] 280   9
table(factors$HPV_clade)
#A7_clade            A9_clade Solid Tissue Normal 
#65                 202                  13

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
explo.plot(mycd, samples = c(268:279))
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
#convert as integer
counts <- round(as.matrix(counts))
storage.mode(counts) <- "integer"
stopifnot(all(counts >= 0, na.rm = TRUE))

counts_filtered <- filtered.data(
  counts, factor = "HPV_clade", norm = FALSE, depth = NULL, method = 1, cpm = 0, p.adj = "fdr")
#13705 features are to be kept for differential expression analysis with filtering method 1
length(unique(rownames(counts_filtered))) #[1] 13705

#sync
myannot_filtered <- myannot[myannot$feature %in% rownames(counts_filtered), ]
stopifnot(all(rownames(counts_filtered) %in% myannot_filtered$feature))

# Set rownames 
rownames(myannot_filtered) <- myannot_filtered$feature
factors <- factors[colnames(counts_filtered), , drop = FALSE]
rownames(factors) <- factors$specimenID
# Ensure full synchrony
stopifnot(all(rownames(factors) == colnames(counts_filtered)))
counts_filtered <- counts_filtered[rowSums(is.na(counts_filtered))==0 & rowSums(counts_filtered)>0, colSums(is.na(counts_filtered))==0 & colSums(counts_filtered)>0, drop=FALSE]
#
myannot_filtered_df <- myannot_filtered %>%
  as.data.frame() %>%
  mutate(
    length     = as.integer(length),
    gc_percent = as.numeric(gc_percent)
  )
rownames(myannot_filtered_df) <- myannot_filtered_df$feature

stopifnot(all(rownames(counts_filtered) %in% rownames(myannot_filtered_df)))
myannot_filtered_df <- myannot_filtered_df[rownames(counts_filtered), , drop = FALSE]

# Create EDASeq object
mydataEDA <- newSeqExpressionSet(
  counts = as.matrix(counts_filtered),
  featureData = AnnotatedDataFrame(data = myannot_filtered_df),
  phenoData = AnnotatedDataFrame(data = factors))

#order for less bias
gcFull <- withinLaneNormalization(mydataEDA, 
                                  "gc_percent", which = "full")#corrects GC bias 
lFull <- withinLaneNormalization(gcFull, "length", which = "full")#corrects length bias 
fullfullTMM <-NOISeq::tmm(normCounts(lFull), long = 1000, lc = 0, k = 0)
#norm.counts <- betweenLaneNormalization(normCounts(lFull),
# which = "median", offset = FALSE)

noiseqData_filtered = NOISeq::readData(data = fullfullTMM, 
                                   factors=factors)

dim(fullfullTMM) #[1] 13705   280
#myPCA_pre = dat(noiseqData_filtered, type = "PCA", norm = T,logtransf = F)
#explo.plot(myPCA_sample, plottype= "scores", factor = "HPV_clade")
#explo.plot(myPCA_sample, plottype= "scores", factor = "source")
#dev.off()
mycd_filtered <- dat(noiseqData_filtered, type = "cd", norm = TRUE ) # Slow operation

# Verify diagnostic
table(mycd_filtered@dat$DiagnosticTest[, "Diagnostic Test"])
#FAILED PASSED 
#29    260 

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
#expr_logCPM <- log2(fullfullTMM + 1)

# Correct batch using limma::removeBatchEffect
# We adjust batch (source) while preserving sample_type biological signal

#expr_batchcorrected <- removeBatchEffect(
#  expr_logCPM,
#  batch = factors$source,
#  design = model.matrix(~ sample_type, data = factors))

###
#ffTMMARSyn=ARSyNseq(noiseqData_filtered, factor = "HPV_type", batch = F,
#                    norm = "n",  logtransf = F)

ffTMMARSyn_s=ARSyNseq(noiseqData_filtered, factor = "source", batch = T,
                    norm = "n",  logtransf = F)

###
#noiseqData_batch = NOISeq::readData(data = expr_batchcorrected, 
#                                       factors=factors)

#myPCA_after = dat(ffTMMARSyn  , type = "PCA", norm = T,logtransf = T)

#pdf("plots_prepro_data/QC_PCA_scores_HPVclade_after.pdf", width = 8, height = 6)
#explo.plot(myPCA_after, plottype = "scores", factor = "source")
#dev.off()

#pdf("plots_prepro_data/QC_PCA_scores_source_after.pdf", width = 8, height = 6)
#explo.plot(myPCA_after, plottype = "scores", factor = "HPV_clade")
#dev.off()


myPCA_after_s = dat(ffTMMARSyn_s  , type = "PCA", norm = T,logtransf = T)

pdf("plots_prepro_data/QC_PCA_scores_HPVclade_after_s.pdf", width = 8, height = 6)
explo.plot(myPCA_after_s, plottype = "scores", factor = "source")
dev.off()

pdf("plots_prepro_data/QC_PCA_scores_source_after_s.pdf", width = 8, height = 6)
explo.plot(myPCA_after_s, plottype = "scores", factor = "HPV_clade")
dev.off()

# ===================== FINAL QUALITY CHECK ===================== #
mylength_n <- setNames(myannot_filtered_df$length, myannot_filtered_df$feature)
mygcn_n <- setNames(myannot_filtered_df$gc_percent, myannot_filtered_df$feature)
mybiotype_n <- setNames(myannot_filtered_df$gene_biotype, myannot_filtered_df$feature)

noiseqData_norm_s = NOISeq::readData(data = ffTMMARSyn_s, 
                                   factors=factors,
                                   gc = mygcn_n,
                                   biotype = mybiotype_n,
                                   length = mylength_n)

# 3.3.3) RNA composition
mycd_norm_2 <- dat(noiseqData_norm_s, type = "cd", norm = TRUE)
# Verify diagnostic
table(mycd_norm_2@dat$DiagnosticTest[, "Diagnostic Test"])
#FAILED PASSED bind #b=T source
#   36    243 

# ===================== IDENTIFY QUALITY CONTROL AFTER NORMALIZATION ===================== #
### 3) Quality control of count data--- ---
# 3.2) Sequencing depth & Expression Quantification
# 3.2.1) Count distribution per sample
pdf("plots_prepro_data/QC_boxplot_counts_after.pdf", width = 8, height = 6)
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

pdf("plots_prepro_data/QC_length_bias_SNT_after.pdf", width = 8, height = 6)
par(cex.axis=1.5, cex.lab=1.5, cex.main=2)
explo.plot(mylengthbias_SampleType, samples = 3, toplot = "global")
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

pdf("plots_prepro_data/QC_GC_bias_SNT_after.pdf", width = 8, height = 6)
par(cex.axis=1.5, cex.lab=1.5, cex.main=2)
explo.plot(myGCbias_Sample, samples = 3, toplot = "global")
dev.off()

# 3.3.3) RNA composition
mycd_norm<- dat(noiseqData_norm_s, type = "cd", norm = TRUE)
# Verify diagnostic
table(mycd_norm@dat$DiagnosticTest[, "Diagnostic Test"])

#FAILED PASSED #b=F HPV_clade 
# 62    217 

pdf("plots_prepro_data/QC_Mvalues_mycd_after.pdf", width = 10, height = 6)
par(mar = c(5, 4, 4, 2) - 1)
par(cex.axis = 1.5, cex.lab = 1.5, cex.main = 2, xaxs = "i")
explo.plot(mycd_norm, samples = c(266:276))
dev.off()

# 3.4) PCA
myPCA_sample_after = dat(noiseqData_norm_s, type = "PCA", norm = TRUE, logtransf = FALSE)
pdf("plots_prepro_data/QC_PCA_scores_HPVclade_after.pdf", width = 8, height = 6)
explo.plot(myPCA_sample_after, plottype= "scores", factor = "HPV_clade")
dev.off()

pdf("plots_prepro_data/QC_PCA_scores_source_after.pdf", width = 8, height = 6)
explo.plot(myPCA_sample_after, plottype= "scores", factor = "source")
dev.off()

# ===================== SAVE FINAL  ===================== #
Final_SNT_vs_PP <- exprs(ffTMMARSyn_s) %>% as.data.frame()
Final_SNT_vs_PP$Gene <- rownames(Final_SNT_vs_PP)  #Add column for input ARACNE
Final_SNT_vs_PP <- Final_SNT_vs_PP[, c(ncol(Final_SNT_vs_PP), 1:(ncol(Final_SNT_vs_PP)-1))]  
dim(Final_SNT_vs_PP) #[1] 13705   281

# Save full matrix
vroom::vroom_write(Final_SNT_vs_PP, file = "~/CESC_Network/2_Prepro_TCGA_GTEx/2_3_Full_counts_A7_A9_notAnot_13705.tsv", delim = "\t")
# Save metadata
vroom::vroom_write(factors, file = "~/CESC_Network/2_Prepro_TCGA_GTEx/2_4_Factors.tsv", delim = "\t")
#Save annotation
vroom::vroom_write(myannot_filtered, file = "~/CESC_Network/2_Prepro_TCGA_GTEx/2_5_Myannot.tsv", delim = "\t")

# Matrix annoted
Final_SNT_vs_PP$HGNC_symbol <- myannot_filtered_df$hgnc_symbol[match(Final_SNT_vs_PP$Gene, myannot_filtered_df$feature)]
Final_SNT_vs_PP_annot <- Final_SNT_vs_PP %>% filter(!is.na(HGNC_symbol) & HGNC_symbol != "")
Final_SNT_vs_PP_annot <- Final_SNT_vs_PP_annot %>%
  dplyr::select(HGNC_symbol, everything(), -Gene) %>%
  dplyr::rename(Gene = HGNC_symbol)
Final_SNT_vs_PP_annot <- Final_SNT_vs_PP_annot %>% dplyr::group_by(Gene) %>% dplyr::summarise(dplyr::across(where(is.numeric), ~mean(., na.rm=TRUE)), .groups="drop"); rownames(Final_SNT_vs_PP_annot) <- Final_SNT_vs_PP_annot$Gene
dim(Final_SNT_vs_PP_annot)
#[1] 13677   281

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

save.image("/STORAGE/csbig/jruiz/cesc_network/2_9_Image_Post_Prepro.RData")

#load("/STORAGE/csbig/jruiz/cesc_network/2_9_Image_Post_Prepro.RData")

