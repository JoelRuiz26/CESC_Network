# ===================== SCRIPT OVERVIEW ===================== #
# Title: Prepro TCGA and GTEX data
# Author: Joel Ruiz Hernandez
# Contact: josejoelruizhernandez@gmail.com
# =========================================================== #

# ===================== LOAD REQUIRED LIBRARIES ===================== #
library(biomaRt)
library(NOISeq)
library(edgeR)
library(EDASeq)
library(vroom)
library(tidyverse)
library(ggplot2)

# ===================== SEED FOR REPRODUCIBILITY ===================== #
set.seed(123)

# ===================== LOAD INPUT DATA ===================== #
factors <- readRDS(file = "~/1_Get_Data_TCGA/2_Harmoni_metadata_TCGA_GTEx.rds")
unstranded_counts <- readRDS(file = "~/1_Get_Data_TCGA/2_Harmoni_counts_TCGA_GTEx.rds")

# ===================== STEP 1: CLEAN AND ANNOTATE DATA ===================== #
# Clean counts matrix
counts <- as.data.frame(unstranded_counts) %>%
  mutate(gene = str_remove(gene, "\\..*$"))  # remove Ensembl version if present

# Remove genes with all zero counts
counts <- counts[rowSums(counts[,-1]) > 0, ]

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

# Reorder myannot to match counts order
myannot <- myannot %>%
  filter(feature %in% counts$gene) %>%
  arrange(match(feature, counts$gene))

# Añadir anotación HGNC a la matriz de conteos
annot_hgnc <- myannot %>% dplyr::select(feature, hgnc_symbol)
counts_annotated <- counts %>%
  dplyr::rename(feature = gene) %>%
  inner_join(annot_hgnc, by = "feature") %>%
  dplyr::select(hgnc_symbol, everything(), -feature) %>%
  dplyr::rename(Gene = hgnc_symbol) %>%
  distinct(Gene, .keep_all = TRUE)
rownames(counts_annotated) <- counts_annotated$Gene
counts_annotated <- counts_annotated %>% dplyr::select(-Gene)

# Filter and sync metadata
matched_samples <- factors$specimenID[factors$specimenID %in% colnames(counts_annotated)]
factors <- factors %>% filter(specimenID %in% matched_samples)
# Reorder columns of counts to match order of factors
counts <- counts_annotated %>% dplyr::select(all_of(factors$specimenID))
# Convert factors to data.frame, drop extra columns, and align rownames
factors <- droplevels(as.data.frame(factors))
rownames(factors) <- factors$specimenID
factors$specimenID <- rownames(factors)
# Ensure synchrony
stopifnot(all(colnames(counts) == rownames(factors)))
# Prepare final matrix
counts <- as.data.frame(counts)

# ===================== STEP 2: CREATE NOISeq OBJECT ===================== #
# Create named vectors for annotations
mylength <- setNames(myannot$length, myannot$hgnc_symbol)
mygc <- setNames(myannot$percentage_gene_gc_content, myannot$hgnc_symbol)
mybiotype <- setNames(myannot$gene_biotype, myannot$hgnc_symbol)

# Create NOISeq object
noiseqData_before <- NOISeq::readData(
  data = counts,
  factors = factors,
  gc = mygc,
  biotype = mybiotype,
  length = mylength
)
#save.image("~/2_Prepro_Data_TCGA/2_Image_Pre_Prepro.RData")
#load("~/2_Prepro_Data_TCGA/2_Image_Pre_Prepro.RData")

### 3) Quality control of count data--- ---
# 3.2) Sequencing depth & Expression Quantification
# 3.2.1) Count distribution per sample  
mycountsbio_HPV = NOISeq::dat(noiseqData_beforeNormal, type = "countsbio", factor = "HPV_clade")# counts per million with if norm=true
explo.plot(mycountsbio_HPV,
           plottype = "boxplot", #type of plot
           samples = NULL)

# 3.2.2) Sensitivity ploT: check for low count genes
explo.plot(mycountsbio_HPV,
           plottype = "barplot", # tipo de gráfico
           samples = NULL) 
dev.off()

# Plot global distribution of CPM 
log_cpm_means <- rowMeans(cpm(noiseqData_beforeNormal@assayData$exprs, log = TRUE))
# Calcula el histograma sin graficarlo
hist_data <- hist(log_cpm_means, plot = FALSE, breaks = 30)
hist(log_cpm_means,
     ylab = "Number of genes",
     xlab = "log CPM averages",
     col = "lightblue",
     border = "darkblue",
     main = "CPM Distribution",
     breaks = 30,
     xlim = c(min(log_cpm_means), max(log_cpm_means)),
     ylim = c(0, max(hist_data$counts) * 1.1),  # +10% por estética
     xaxs = "i",
     yaxs = "i",
     cex.axis = 1.5,
     cex.lab = 1.5,
     cex.main = 2)

# 3.3) Sequencing bias detection
# 3.3.1) Length bias
#Relationship between the gene length and the expression values
#If the model p-value is significant and R2 value is high (more than 70%), the
#expression depends on the feature length and the curve shows the type of dependence.
# Verifica cuántos nombres coinciden
#sum(rownames(counts) %in% names(mylength)) == nrow(counts)  # debería dar TRUE

mylengthbias_SampleType = NOISeq::dat(noiseqData_before, factor = "HPV_clade", type = "lengthbias")
par(cex.axis=1.5,   
    cex.lab=1.5,     
    cex.main=2)
explo.plot(mylengthbias_SampleType, samples = NULL, toplot = "global")
explo.plot(mylengthbias_SampleType, samples = c(1), toplot = "global")
explo.plot(mylengthbias_SampleType, samples = c(2), toplot = "global")
explo.plot(mylengthbias_SampleType, samples = c(3), toplot = "global")

dev.off()


# 3.3.2) GC content bias
# Relationship between the gene GC content and the expression values
# If the model p-value is signifficant and R2 value is high (more than 70%), the expression will depend on
#the feature GC content and the curve will show the type of dependence.

myGCbias_Sample = dat(noiseqData_before, factor = "HPV_clade", type = "GCbias")
par(cex.axis=1.5,   
    cex.lab=1.5,     
    cex.main=2)
explo.plot(myGCbias_Sample, samples = NULL, toplot = "global")
explo.plot(myGCbias_Sample, samples = 1, toplot = "global")
explo.plot(myGCbias_Sample, samples = 2, toplot = "global")
explo.plot(myGCbias_Sample, samples = 3, toplot = "global")

dev.off() 

# 3.3.3) RNA composition
#Diagnostic of data

#each sample "s" is compared to a reference "r" (which can be arbitrarily chosen).
#by computing M values=log2(counts = countsr). 
#Confidence intervals (CI) for the M median is computed by bootstrapping.
#If the median of M values for each comparison is not in the CI, the deviation
# of the sample is significant, therefore, normalization is needed 
#"cd" means "Cumulative Distribution.

mycd <- dat(noiseqData_before, type = "cd", norm = FALSE) # Slow operation

#[1] "Warning: 4201 features with 0 counts in all samples are to be removed for this analysis."
#[1] "Reference sample is: TCGA-C5-A1M5-01A-11R-A13Y-07"
#[1] "Confidence intervals for median of M:"
#0.01%                   99.99%  
#[1] "Diagnostic test: FAILED. Normalization is required to correct this bias."

# Verify diagnostic
table(mycd@dat$DiagnosticTest[, "Diagnostic Test"])
#FAILED PASSED 
#264     25 

#Plot cpm vs counts
#png("7_Mvalues_mycd.png", width=1000, height=600)
par(mar = c(5, 4, 4, 2) - 1)  # reduce márgenes, especialmente izquierda/derecha
par(cex.axis = 1.5, cex.lab = 1.5, cex.main = 2, xaxs = "i")
explo.plot(mycd, samples = c(266:277)) # Between both kind of samples
#dev.off()

# 3.4) PCA
#used to visualize if the experimental samples are clustered according to the experimental desig
myPCA_sample = dat(noiseqData_before, type = "PCA", norm = F,logtransf = F)
explo.plot(myPCA_sample, plottype= "scores", factor = "HPV_clade")
explo.plot(myPCA_sample, plottype= "scores", factor = "source")
dev.off()


# ===================== STEP 3: FILTER FOR EDASeq QC ===================== #
counts_filtered <- filtered.data(
  counts, factor = "HPV_clade", norm = FALSE, depth = NULL, method = 1, cpm = 0, p.adj = "fdr")
#11531 features are to be kept for differential expression analysis with filtering method 1

myannot_filtered <- myannot[myannot$hgnc_symbol %in% rownames(counts_filtered), ]
stopifnot(all(rownames(counts_filtered) %in% myannot_filtered$hgnc_symbol))

# Set rownames for annotation and metadata
rownames(myannot_filtered) <- myannot_filtered$hgnc_symbol
factors <- factors[colnames(counts_filtered), , drop = FALSE]
rownames(factors) <- factors$specimenID

# Ensure full synchrony
stopifnot(all(rownames(factors) == colnames(counts_filtered)))

# Create EDASeq object
mydataEDA <- newSeqExpressionSet(
  counts = as.matrix(counts_filtered),
  featureData = AnnotatedDataFrame(data = myannot_filtered),
  phenoData = AnnotatedDataFrame(data = factors)
)
#order for less bias
lFull <- withinLaneNormalization(mydataEDA, "length", which = "full")#corrects length bias 
gcFull <- withinLaneNormalization(lFull, 
                                  "percentage_gene_gc_content", which = "full")#corrects GC bias 
fullfullTMM <-NOISeq::tmm(normCounts(gcFull), long = 1000, lc = 0, k = 0)
#norm.counts <- betweenLaneNormalization(normCounts(lFull),
# which = "median", offset = FALSE)

noiseqData_filtered = NOISeq::readData(data = fullfullTMM, 
                                   factors=factors)

myPCA_pre = dat(noiseqData_filtered, type = "PCA", norm = T,logtransf = F)
explo.plot(myPCA_sample, plottype= "scores", factor = "HPV_clade")
explo.plot(myPCA_sample, plottype= "scores", factor = "source")
dev.off()


mycd <- dat(noiseqData_filtered, type = "cd", norm = TRUE ) # Slow operation

#[1] "Warning: 4201 features with 0 counts in all samples are to be removed for this analysis."
#[1] "Reference sample is: TCGA-C5-A1M5-01A-11R-A13Y-07"
#[1] "Confidence intervals for median of M:"
#0.01%                   99.99%  
#[1] "Diagnostic test: FAILED. Normalization is required to correct this bias."

# Verify diagnostic
table(mycd@dat$DiagnosticTest[, "Diagnostic Test"])
#FAILED PASSED 
#28    261



#############################SOLVE BATCH EFFECT#######################################################


ffTMMARSyn=ARSyNseq(noiseqData_filtered, factor = "source", batch = T,
                    norm = "n",  logtransf = T)

myPCA = dat(ffTMMARSyn, type = "PCA", norm = T,logtransf = T)
explo.plot(myPCA, plottype = "scores", factor = "source")
explo.plot(myPCA, plottype = "scores", factor = "HPV_clade")


#############################FINAL QUALITY CHECK#######################################################


mylength_n <- setNames(myannot_filtered$length, myannot_filtered$hgnc_symbol)
mygcn_n <- setNames(myannot_filtered$percentage_gene_gc_content, myannot_filtered$hgnc_symbol)
mybiotype_n <- setNames(myannot_filtered$gene_biotype, myannot_filtered$hgnc_symbol)

noiseqData_norm = NOISeq::readData(data = exprs(ffTMMARSyn), 
                                   factors=factors,
                                   gc = mygcn_n,
                                   biotype = mybiotype_n,
                                   length = mylength_n)


################3Plots
### 3) Quality control of count data--- ---
# 3.2) Sequencing depth & Expression Quantification
# 3.2.1) Count distribution per sample  
mycountsbio_HPV = NOISeq::dat(noiseqData_norm, type = "countsbio", factor = "HPV_clade", norm = T)# counts per million with if norm=true
explo.plot(mycountsbio_HPV,
           plottype = "boxplot", #type of plot
           samples = NULL)

# 3.2.2) Sensitivity ploT: check for low count genes
explo.plot(mycountsbio_HPV,
           plottype = "barplot", # tipo de gráfico
           samples = NULL) 
dev.off()

# Plot global distribution of CPM 
log_cpm_means <- rowMeans(cpm(noiseqData_norm@assayData$exprs, log = TRUE))
# Calcula el histograma sin graficarlo
hist_data <- hist(log_cpm_means, plot = FALSE, breaks = 30)
hist(log_cpm_means,
     ylab = "Number of genes",
     xlab = "log CPM averages",
     col = "lightblue",
     border = "darkblue",
     main = "CPM Distribution",
     breaks = 30,
     xlim = c(min(log_cpm_means), max(log_cpm_means)),
     ylim = c(0, max(hist_data$counts) * 1.1),  # +10% por estética
     xaxs = "i",
     yaxs = "i",
     cex.axis = 1.5,
     cex.lab = 1.5,
     cex.main = 2)

# 3.3) Sequencing bias detection
# 3.3.1) Length bias
#Relationship between the gene length and the expression values
#If the model p-value is significant and R2 value is high (more than 70%), the
#expression depends on the feature length and the curve shows the type of dependence.
# Verifica cuántos nombres coinciden
#sum(rownames(counts) %in% names(mylength)) == nrow(counts)  # debería dar TRUE

mylengthbias_SampleType = NOISeq::dat(noiseqData_norm, factor = "HPV_clade", type = "lengthbias",norm = T)
par(cex.axis=1.5,   
    cex.lab=1.5,     
    cex.main=2)
explo.plot(mylengthbias_SampleType, samples = NULL, toplot = "global")
explo.plot(mylengthbias_SampleType, samples = c(1), toplot = "global")
explo.plot(mylengthbias_SampleType, samples = c(2), toplot = "global")
explo.plot(mylengthbias_SampleType, samples = c(3), toplot = "global")
dev.off()


# 3.3.2) GC content bias
# Relationship between the gene GC content and the expression values
# If the model p-value is signifficant and R2 value is high (more than 70%), the expression will depend on
#the feature GC content and the curve will show the type of dependence.

myGCbias_Sample = dat(noiseqData_norm, factor = "HPV_clade", type = "GCbias",norm = T)
par(cex.axis=1.5,   
    cex.lab=1.5,     
    cex.main=2)
explo.plot(myGCbias_Sample, samples = NULL, toplot = "global")
explo.plot(myGCbias_Sample, samples = 1, toplot = "global")
explo.plot(myGCbias_Sample, samples = 2, toplot = "global")
explo.plot(myGCbias_Sample, samples = 3, toplot = "global")

dev.off() 

# 3.3.3) RNA composition
#Diagnostic of data

#each sample "s" is compared to a reference "r" (which can be arbitrarily chosen).
#by computing M values=log2(counts = countsr). 
#Confidence intervals (CI) for the M median is computed by bootstrapping.
#If the median of M values for each comparison is not in the CI, the deviation
# of the sample is significant, therefore, normalization is needed 
#"cd" means "Cumulative Distribution.

mycd <- dat(noiseqData_norm, type = "cd", norm = TRUE ) # Slow operation

#[1] "Warning: 4201 features with 0 counts in all samples are to be removed for this analysis."
#[1] "Reference sample is: TCGA-C5-A1M5-01A-11R-A13Y-07"
#[1] "Confidence intervals for median of M:"
#0.01%                   99.99%  
#[1] "Diagnostic test: FAILED. Normalization is required to correct this bias."

# Verify diagnostic
table(mycd@dat$DiagnosticTest[, "Diagnostic Test"])
#FAILED PASSED 
#27    262

#Plot cpm vs counts
#png("7_Mvalues_mycd.png", width=1000, height=600)
par(mar = c(5, 4, 4, 2) - 1)  # reduce márgenes, especialmente izquierda/derecha
par(cex.axis = 1.5, cex.lab = 1.5, cex.main = 2, xaxs = "i")
explo.plot(mycd, samples = c(266:277)) # Between both kind of samples
#dev.off()


# ===================== STEP 4: SAVE FINAL  ===================== #

Final_SNT_vs_PP <- exprs(ffTMMARSyn) %>% as.data.frame()
Final_SNT_vs_PP$Gene <- rownames(Final_SNT_vs_PP)  # copiar rownames a columna Gene
Final_SNT_vs_PP <- Final_SNT_vs_PP[, c(ncol(Final_SNT_vs_PP), 1:(ncol(Final_SNT_vs_PP)-1))]  

# Guardar matriz completa
vroom::vroom_write(Final_SNT_vs_PP, file = "~/2_Prepro_Data_TCGA/1_Full_counts_A7_A9_t.tsv", delim = "\t")

# Guardar metadatos sincronizados
vroom::vroom_write(factors, file = "~/2_Prepro_Data_TCGA/1_Factors.tsv", delim = "\t")
#Save annotation
vroom::vroom_write(myannot_filtered, file = "~/2_Prepro_Data_TCGA/1_Myannot.tsv", delim = "\t")

# Subconjunto A7
Factors_A7 <- factors %>% filter(HPV_clade == "A7_clade")
Final_A7_NW <- Final_SNT_vs_PP %>%
  dplyr::select(Gene, all_of(Factors_A7$specimenID))
vroom::vroom_write(Final_A7_NW, file = "~/2_Prepro_Data_TCGA/1_Final_A7_NW_t.tsv", delim = "\t")

# Subconjunto A9
Factors_A9 <- factors %>% filter(HPV_clade == "A9_clade")
Final_A9_NW <- Final_SNT_vs_PP %>%
  dplyr::select(Gene, all_of(Factors_A9$specimenID))
vroom::vroom_write(Final_A9_NW, file = "~/2_Prepro_Data_TCGA/1_Final_A9_NW_t.tsv", delim = "\t")

#save.image("~/2_Prepro_Data_TCGA/2_Image_Post_Prepro.RData")


