# ===================== SCRIPT OVERVIEW ===================== #
# Title: Differential Expression Analysis by HPV_clade with Limma
# Author: Joel Ruiz Hernandez (adapted)
#
# This script assumes that the expression input matrix is already:
# - Filtered (low-count genes removed)
# - Normalized (GC, length, TMM)
# - Log2-transformed
# - Batch-corrected for source while preserving sample_type

# =========================================================== #

# ===================== LOAD REQUIRED LIBRARIES ===================== #

library(limma)
library(vroom)
library(ggplot2)
library(dplyr)

# ===================== SET WORKING DIRECTORY ===================== #

setwd("~/CESC_Network/3_DGE/")

# ===================== LOAD DATA ===================== #
expr <- vroom::vroom("~/CESC_Network/2_Prepro_TCGA_GTEx/2_6_Full_counts_A7_A9_annot.tsv")
expr_mat <- as.data.frame(expr)
rownames(expr_mat) <- expr_mat$Gene
expr_mat$Gene <- NULL

factors <- vroom::vroom("~/CESC_Network/2_Prepro_TCGA_GTEx/2_4_Factors.tsv")
factors <- as.data.frame(factors)
rownames(factors) <- factors$specimenID
expr_mat <- expr_mat[, rownames(factors)]

stopifnot(all(colnames(expr_mat) == rownames(factors)))

# ===================== DESIGN MATRIX ===================== #
design <- model.matrix(~ 0 + HPV_clade, data=factors)
colnames(design) <- make.names(colnames(design))
print(colnames(design))

# ===================== CONTRASTS ===================== #
contrast.matrix <- makeContrasts(
  A7vsSNT = HPV_cladeA7_clade - HPV_cladeSolid.Tissue.Normal,
  A9vsSNT = HPV_cladeA9_clade - HPV_cladeSolid.Tissue.Normal,
  levels=design
)
print(contrast.matrix)

# ===================== LIMMA ANALYSIS ===================== #
fit <- lmFit(expr_mat, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# ===================== THRESHOLDS ===================== #
fdr_threshold <- 0.01
logfc_threshold <- 1

# ===================== MAKE OUTPUT FOLDER FOR PLOTS ===================== #
if (!dir.exists("plots_DGE")) dir.create("plots_DGE")

# ===================== ANALYSIS FOR A7 vs SNT ===================== #
# Extract results
top_A7 <- topTable(fit2, coef="A7vsSNT", number=Inf, sort.by="P")
top_A7$Gene <- rownames(top_A7)
dim(top_A7) #[1] 11544     7

# Save ALL results
write.table(top_A7, "3_1_DGE_A7_vs_SNT_AllResults.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# Filter significant DEGs
A7_DEGs <- top_A7 %>%
  filter(adj.P.Val < fdr_threshold & abs(logFC) >= logfc_threshold)
dim(A7_DEGs) #[1] 3149    7

write.table(A7_DEGs, "3_2_DGE_A7_vs_SNT_Significant_logFC1.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# Split into Up and Down
A7_Up <- A7_DEGs %>% filter(logFC >= logfc_threshold)
A7_Down <- A7_DEGs %>% filter(logFC <= -logfc_threshold)
nrow(A7_Up) #[1] 1703
nrow(A7_Down) #[1] 1446

# Save Up/Down tables
write.table(A7_Up, "3_2_1_DGE_A7_vs_SNT_Upregulated.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(A7_Down, "3_2_2DGE_A7_vs_SNT_Downregulated.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# Volcano Plot
A7_DEGs$Regulation <- ifelse(A7_DEGs$logFC >= logfc_threshold, "Upregulated", "Downregulated")
volcano_colors <- c("Upregulated" = "red", "Downregulated" = "blue")

pdf("plots_DGE/VolcanoPlot_A7_vs_SNT.pdf", width=8, height=6)
ggplot(top_A7, aes(x=logFC, y=-log10(adj.P.Val))) +
  geom_point(alpha=0.4, color="gray80") +
  geom_point(data=A7_DEGs, aes(color=Regulation), alpha=0.8, size=2) +
  scale_color_manual(values=volcano_colors) +
  theme_bw(base_size=18) +
  labs(title="Volcano Plot: A7 vs SNT",
       x="log2 Fold Change",
       y="-log10 Adjusted P-Value (FDR)",
       color="Regulation") +
  geom_vline(xintercept=c(-logfc_threshold, logfc_threshold), linetype="dashed") +
  geom_hline(yintercept=-log10(fdr_threshold), linetype="dashed")
dev.off()

# ===================== ANALYSIS FOR A9 vs SNT ===================== #
# Extract results
top_A9 <- topTable(fit2, coef="A9vsSNT", number=Inf, sort.by="P")
top_A9$Gene <- rownames(top_A9)
dim(top_A9) #[1] 11544     7

# Save ALL results
write.table(top_A9, "3_3_DGE_A9_vs_SNT_AllResults.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# Filter significant DEGs
A9_DEGs <- top_A9 %>%
  filter(adj.P.Val < fdr_threshold & abs(logFC) >= logfc_threshold)
dim(A9_DEGs) #[1] 3279    7
write.table(A9_DEGs, "3_4_DGE_A9_vs_SNT_Significant_logFC1.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# Split into Up and Down
A9_Up <- A9_DEGs %>% filter(logFC >= logfc_threshold)
A9_Down <- A9_DEGs %>% filter(logFC <= -logfc_threshold)
nrow(A9_Up) #[1] 1785
nrow(A9_Down) #[1] 1494

# Save Up/Down tables
write.table(A9_Up, "3_4_1_DGE_A9_vs_SNT_Upregulated.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(A9_Down, "3_4_1_DGE_A9_vs_SNT_Downregulated.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# Volcano Plot
A9_DEGs$Regulation <- ifelse(A9_DEGs$logFC >= logfc_threshold, "Upregulated", "Downregulated")

pdf("plots_DGE/VolcanoPlot_A9_vs_SNT.pdf", width=8, height=6)
ggplot(top_A9, aes(x=logFC, y=-log10(adj.P.Val))) +
  geom_point(alpha=0.4, color="gray80") +
  geom_point(data=A9_DEGs, aes(color=Regulation), alpha=0.8, size=2) +
  scale_color_manual(values=volcano_colors) +
  theme_bw(base_size=18) +
  labs(title="Volcano Plot: A9 vs SNT",
       x="log2 Fold Change",
       y="-log10 Adjusted P-Value (FDR)",
       color="Regulation") +
  geom_vline(xintercept=c(-logfc_threshold, logfc_threshold), linetype="dashed") +
  geom_hline(yintercept=-log10(fdr_threshold), linetype="dashed")
dev.off()


#Smyth 2015
#https://pmc.ncbi.nlm.nih.gov/articles/PMC4402510/#:~:text=overcome%20the%20problem%20of%20small,Recently
#Es el método recomendado en Bioconductor para datos ya transformados/log2CPM.
