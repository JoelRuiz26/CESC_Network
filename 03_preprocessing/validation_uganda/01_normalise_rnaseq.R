### Joel Ruiz
### Preprocess the Uganda (HTMCP-CC) validation cohort exactly as
### 03_preprocessing/01_normalise_rnaseq.R does for TCGA-CESC: NOISeq QC,
### low-count filtering, EDASeq GC + length correction, TMM. No batch
### correction is applied here either -- the Fisher's exact check
### (02_cohort_data/validation_uganda/02_fisher_clinical_vs_clade.R) found
### HIV status independent of HPV clade (p = 0.373), so there is no
### confounding to remove, matching the reasoning already used for CESC
### (03_preprocessing/01_normalise_rnaseq.R, "BATCH EFFECT (OMITTED)").
###
### Reuses CESC's own gene annotation (03_preprocessing/2_1_Universe_
### annotation_table.rds: Ensembl ID -> length, GC%, symbol, biotype) rather
### than re-deriving it, since these are static genome properties, not
### cohort-specific, and using the same values keeps the two cohorts'
### normalisation directly comparable.
setwd("~/CESC_Network/03_preprocessing/validation_uganda")
suppressPackageStartupMessages({
  library(NOISeq); library(edgeR); library(EDASeq); library(tidyverse)
})
set.seed(123)
dir.create("plots_prepro_data", showWarnings = FALSE)

## ----------------- Load input data ----------------- #
raw <- readRDS("~/CESC_Network/02_cohort_data/validation_uganda/1_Output/3_1_raw_gene_counts_A7A9.rds")
meta <- vroom::vroom("~/CESC_Network/02_cohort_data/validation_uganda/1_Output/1_3_HTMCP_clinical_HPV_HIV_with_rnaseq_flag.tsv",
                     show_col_types = FALSE) |>
  filter(has_rnaseq, HPV.clade %in% c("A7", "A9")) |>
  transmute(sampleID = Patient, HPV_clade = HPV.clade)

gene_metadata <- readRDS("~/CESC_Network/03_preprocessing/2_1_Universe_annotation_table.rds")

## ----------------- Clean counts and samples ----------------- #
counts <- raw$counts
## PAR_Y genes (44, pseudoautosomal duplicates on X/Y) share their stripped
## Ensembl ID with their X-chromosome counterpart; dropped rather than
## arbitrarily merged
counts <- counts[!grepl("PAR_Y", rownames(counts)), , drop = FALSE]
rownames(counts) <- sub("\\..*$", "", rownames(counts))

stopifnot(all(meta$sampleID %in% colnames(counts)))
factors <- meta[match(colnames(counts), meta$sampleID), ] |> as.data.frame()
rownames(factors) <- factors$sampleID
stopifnot(identical(factors$sampleID, colnames(counts)))

counts <- counts[rowSums(counts) > 0, ]
dim(counts) # genes x 109 patients

## ----------------- Build gene annotation ----------------- #
myannot <- gene_metadata |>
  transmute(feature = as.character(feature), length = as.integer(length),
           gc_percent = gc_percent, hgnc_symbol = as.character(hgnc_symbol),
           gene_biotype = gene_biotype) |>
  distinct(feature, .keep_all = TRUE)

common_genes <- rownames(counts)[rownames(counts) %in% myannot$feature]
counts  <- counts[common_genes, , drop = FALSE]
myannot <- myannot[match(common_genes, myannot$feature), ]
stopifnot(identical(rownames(counts), myannot$feature))
dim(counts); dim(myannot); dim(factors)

saveRDS(myannot, "2_1_Universe_annotation_table_uganda.rds")

## ===================== NOISeq object (before) ===================== #
mylength   <- setNames(myannot$length, myannot$feature)
mygc       <- setNames(myannot$gc_percent, myannot$feature)
mybiotype  <- setNames(myannot$gene_biotype, myannot$feature)

noiseqData_before <- NOISeq::readData(data = counts, factors = factors,
                                      gc = mygc, biotype = mybiotype, length = mylength)
saveRDS(noiseqData_before, "2_2_Noiseq_before.rds")

## ===================== QC (before) ===================== #
pdf("plots_prepro_data/QC_boxplot_counts_before.pdf", width = 6, height = 6)
mycountsbio <- NOISeq::dat(noiseqData_before, type = "countsbio", factor = "HPV_clade")
explo.plot(mycountsbio, plottype = "boxplot", samples = NULL); dev.off()

pdf("plots_prepro_data/QC_barplot_sensitivity_low_counts_before.pdf", width = 8, height = 6)
explo.plot(mycountsbio, plottype = "barplot", samples = NULL); dev.off()

pdf("plots_prepro_data/QC_histogram_CPM_distribution_before.pdf", width = 8, height = 6)
log_cpm <- rowMeans(cpm(noiseqData_before@assayData$exprs, log = TRUE))
hist(log_cpm, ylab = "Number of genes", xlab = "log CPM averages", col = "lightblue",
    border = "darkblue", main = "CPM Distribution", breaks = 30,
    cex.axis = 1.5, cex.lab = 1.5, cex.main = 2); dev.off()

mylengthbias <- NOISeq::dat(noiseqData_before, factor = "HPV_clade", type = "lengthbias")
pdf("plots_prepro_data/QC_length_bias_global_before.pdf", width = 8, height = 6)
par(cex.axis=1.5, cex.lab=1.5, cex.main=2); explo.plot(mylengthbias, samples = NULL, toplot = "global"); dev.off()
pdf("plots_prepro_data/QC_length_bias_A7_before.pdf", width = 8, height = 6)
par(cex.axis=1.5, cex.lab=1.5, cex.main=2); explo.plot(mylengthbias, samples = 1, toplot = "global"); dev.off()
pdf("plots_prepro_data/QC_length_bias_A9_before.pdf", width = 8, height = 6)
par(cex.axis=1.5, cex.lab=1.5, cex.main=2); explo.plot(mylengthbias, samples = 2, toplot = "global"); dev.off()

myGCbias <- dat(noiseqData_before, factor = "HPV_clade", type = "GCbias")
pdf("plots_prepro_data/QC_GC_bias_global_before.pdf", width = 8, height = 6)
par(cex.axis=1.5, cex.lab=1.5, cex.main=2); explo.plot(myGCbias, samples = NULL, toplot = "global"); dev.off()
pdf("plots_prepro_data/QC_GC_bias_A7_before.pdf", width = 8, height = 6)
par(cex.axis=1.5, cex.lab=1.5, cex.main=2); explo.plot(myGCbias, samples = 1, toplot = "global"); dev.off()
pdf("plots_prepro_data/QC_GC_bias_A9_before.pdf", width = 8, height = 6)
par(cex.axis=1.5, cex.lab=1.5, cex.main=2); explo.plot(myGCbias, samples = 2, toplot = "global"); dev.off()

mycd <- dat(noiseqData_before, type = "cd", norm = FALSE)
print(table(mycd@dat$DiagnosticTest[, "Diagnostic Test"]))
pdf("plots_prepro_data/QC_Mvalues_mycd_before.pdf", width = 7, height = 6)
par(mar = c(5,4,4,2)-1, cex.axis=1.5, cex.lab=1.5, cex.main=2, xaxs="i")
explo.plot(mycd, samples = 1:12); dev.off()

myPCA <- dat(noiseqData_before, type = "PCA", norm = FALSE, logtransf = FALSE)
pdf("plots_prepro_data/QC_PCA_scores_HPVclade_before.pdf", width = 8, height = 6)
explo.plot(myPCA, plottype = "scores", factor = "HPV_clade"); dev.off()

## ===================== Filter low counts ===================== #
counts <- round(as.matrix(counts)); storage.mode(counts) <- "integer"
stopifnot(all(counts >= 0, na.rm = TRUE))

y     <- edgeR::DGEList(counts = counts, group = factors$HPV_clade)
cpm_y <- edgeR::cpm(y)
keep  <- rowSums(cpm_y > 1) >= min(table(factors$HPV_clade))
counts_filtered <- counts[keep, , drop = FALSE]
dim(counts_filtered)

myannot_filtered <- myannot[match(rownames(counts_filtered), myannot$feature), ]
myannot_filtered_df <- as.data.frame(myannot_filtered)
rownames(myannot_filtered_df) <- myannot_filtered_df$feature
stopifnot(identical(rownames(counts_filtered), rownames(myannot_filtered_df)))

## ===================== EDASeq: GC, length, TMM ===================== #
mydataEDA <- EDASeq::newSeqExpressionSet(
  counts = counts_filtered,
  featureData = Biobase::AnnotatedDataFrame(myannot_filtered_df),
  phenoData   = Biobase::AnnotatedDataFrame(factors))

gcFull   <- withinLaneNormalization(mydataEDA, "gc_percent", which = "full")
lFull    <- withinLaneNormalization(gcFull, "length", which = "full")
fullfullTMM <- NOISeq::tmm(normCounts(lFull), long = 1000, lc = 0, k = 0)
dim(fullfullTMM)

## no batch-effect correction (see header): only Primary Tumor samples, and
## HIV status is not confounded with clade

## ===================== QC (after) ===================== #
mylength_n  <- setNames(myannot_filtered_df$length, myannot_filtered_df$feature)
mygc_n      <- setNames(myannot_filtered_df$gc_percent, myannot_filtered_df$feature)
mybiotype_n <- setNames(myannot_filtered_df$gene_biotype, myannot_filtered_df$feature)
noiseqData_norm <- NOISeq::readData(data = fullfullTMM, factors = factors,
                                    gc = mygc_n, biotype = mybiotype_n, length = mylength_n)

pdf("plots_prepro_data/QC_boxplot_counts_after.pdf", width = 6, height = 6)
mycountsbio2 <- NOISeq::dat(noiseqData_norm, type = "countsbio", factor = "HPV_clade", norm = TRUE)
explo.plot(mycountsbio2, plottype = "boxplot", samples = NULL); dev.off()

pdf("plots_prepro_data/QC_barplot_sensitivity_low_counts_after.pdf", width = 8, height = 6)
explo.plot(mycountsbio2, plottype = "barplot", samples = NULL); dev.off()

pdf("plots_prepro_data/QC_histogram_CPM_distribution_after.pdf", width = 8, height = 6)
log_cpm2 <- rowMeans(cpm(noiseqData_norm@assayData$exprs, log = TRUE))
hist(log_cpm2, ylab = "Number of genes", xlab = "log CPM averages", col = "lightblue",
    border = "darkblue", main = "CPM Distribution", breaks = 30,
    cex.axis = 1.5, cex.lab = 1.5, cex.main = 2); dev.off()

mylengthbias2 <- NOISeq::dat(noiseqData_norm, factor = "HPV_clade", type = "lengthbias", norm = TRUE)
pdf("plots_prepro_data/QC_length_bias_global_after.pdf", width = 8, height = 6)
par(cex.axis=1.5, cex.lab=1.5, cex.main=2); explo.plot(mylengthbias2, samples = NULL, toplot = "global"); dev.off()
pdf("plots_prepro_data/QC_length_bias_A7_after.pdf", width = 8, height = 6)
par(cex.axis=1.5, cex.lab=1.5, cex.main=2); explo.plot(mylengthbias2, samples = 1, toplot = "global"); dev.off()
pdf("plots_prepro_data/QC_length_bias_A9_after.pdf", width = 8, height = 6)
par(cex.axis=1.5, cex.lab=1.5, cex.main=2); explo.plot(mylengthbias2, samples = 2, toplot = "global"); dev.off()

myGCbias2 <- dat(noiseqData_norm, factor = "HPV_clade", type = "GCbias", norm = TRUE)
pdf("plots_prepro_data/QC_GC_bias_global_after.pdf", width = 8, height = 6)
par(cex.axis=1.5, cex.lab=1.5, cex.main=2); explo.plot(myGCbias2, samples = NULL, toplot = "global"); dev.off()
pdf("plots_prepro_data/QC_GC_bias_A7_after.pdf", width = 8, height = 6)
par(cex.axis=1.5, cex.lab=1.5, cex.main=2); explo.plot(myGCbias2, samples = 1, toplot = "global"); dev.off()
pdf("plots_prepro_data/QC_GC_bias_A9_after.pdf", width = 8, height = 6)
par(cex.axis=1.5, cex.lab=1.5, cex.main=2); explo.plot(myGCbias2, samples = 2, toplot = "global"); dev.off()

mycd_norm <- dat(noiseqData_norm, type = "cd", norm = TRUE)
print(table(mycd_norm@dat$DiagnosticTest[, "Diagnostic Test"]))
pdf("plots_prepro_data/QC_Mvalues_mycd_after.pdf", width = 7, height = 6)
par(mar = c(5,4,4,2)-1, cex.axis=1.5, cex.lab=1.5, cex.main=2, xaxs="i")
explo.plot(mycd_norm, samples = 1:12); dev.off()

myPCA_after <- dat(noiseqData_norm, type = "PCA", norm = TRUE, logtransf = FALSE)
pdf("plots_prepro_data/QC_PCA_scores_HPVclade_after.pdf", width = 8, height = 6)
explo.plot(myPCA_after, plottype = "scores", factor = "HPV_clade"); dev.off()

## ===================== SAVE FINAL ===================== #
Final <- fullfullTMM |> as.data.frame()
Final$Gene <- rownames(Final)
Final <- Final[, c(ncol(Final), 1:(ncol(Final) - 1))]
dim(Final)

vroom::vroom_write(Final, "2_3_Full_counts_A7_A9_notAnot.tsv", delim = "\t")
vroom::vroom_write(factors, "2_4_Factors.tsv", delim = "\t")
vroom::vroom_write(myannot_filtered_df, "2_5_Myannot.tsv", delim = "\t")

Final$Symbol <- myannot_filtered_df$hgnc_symbol[match(Final$Gene, myannot_filtered_df$feature)]
Final <- Final[!is.na(Final$Symbol) & !grepl("\\.", Final$Symbol), ]

expr_cols <- setdiff(colnames(Final), c("Gene", "Symbol"))
keep_rows <- sort(unlist(tapply(seq_len(nrow(Final)), Final$Symbol,
                                function(i) i[which.max(rowMeans(Final[i, expr_cols, drop = FALSE]))]),
                        use.names = FALSE))
Final_annot <- Final[keep_rows, ]
Final_annot$Gene <- Final_annot$Symbol
rownames(Final_annot) <- Final_annot$Gene
Final_annot$Symbol <- NULL
vroom::vroom_write(Final_annot, "2_6_Full_counts_A7_A9_annot.tsv", delim = "\t")

Samples_A7 <- factors$sampleID[factors$HPV_clade == "A7"]
Samples_A9 <- factors$sampleID[factors$HPV_clade == "A9"]
vroom::vroom_write(Final_annot[, c("Gene", Samples_A7)], "2_7_Full_counts_A7_annot.tsv", delim = "\t")
vroom::vroom_write(Final_annot[, c("Gene", Samples_A9)], "2_8_Full_counts_A9_annot.tsv", delim = "\t")

message(sprintf("[OK] A7: %d genes x %d samples | A9: %d genes x %d samples",
                nrow(Final_annot), length(Samples_A7), nrow(Final_annot), length(Samples_A9)))
