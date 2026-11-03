### ============================================================
### 06_latent_blocks.R
###
### Recover SBM-style BLOCKS directly from the ASE embeddings (the second
### half of Sussman et al. 2012's own method, doi:10.1080/01621459.2012.
### 699795: embed, then cluster the embedded points to recover block
### membership), via a Gaussian mixture model (mclust), G chosen by BIC.
###
### FIRST ATTEMPT (discarded): fit mclust SEPARATELY on X7 and on
### X9_aligned. Two independently-fit GMMs place their block boundaries
### independently, so near a boundary a gene can get a different block
### label between the two fits for no biological reason at all -- this
### gave an implausible 86% of genes "changing block", which is mostly
### boundary-placement noise, not signal.
###
### USED INSTEAD: fit ONE mclust model jointly on rbind(X7, X9_aligned),
### i.e. one shared partition of the latent space, then read off each
### network's block membership from that SAME model -- a gene can only
### look like it "changed block" if it moved to a genuinely different
### region of one shared space, not because two separate models drew their
### boundaries differently.
###
### Note: BIC does not plateau here (it keeps improving up to the search
### ceiling) -- i.e. this latent space behaves like a continuous manifold
### (the degree-driven "arch" already visible in the PCA projection) rather
### than a small number of well-separated Gaussian blocks. The ceiling
### (G_MAX) is therefore a pragmatic computational cap, not a discovered
### "true" number of blocks; this is itself informative and reported as
### such.
###
### -> results/06_latent_blocks.tsv  (gene, block_A7, block_A9)
### ============================================================
suppressPackageStartupMessages({ library(mclust); library(dplyr); library(readr) })

set.seed(123)
res_dir <- "~/CESC_Network/11_latent_embedding/results"
G_MAX   <- 15

X7         <- readRDS(file.path(res_dir, "01_ase_X7.rds"))
X9_aligned <- readRDS(file.path(res_dir, "01_ase_X9_aligned.rds"))
genes      <- read_tsv(file.path(res_dir, "01_shared_gene_order.tsv"), show_col_types = FALSE)$gene

## -- ONE joint GMM on both networks' embeddings stacked together --
fit <- Mclust(rbind(X7, X9_aligned), G = 2:G_MAX, modelNames = c("EII", "VII", "EEI", "VVI"))
cat(sprintf("Joint model: BIC-selected G = %d (model %s)%s\n", fit$G, fit$modelName,
            if (fit$G == G_MAX) " -- hit the search ceiling, BIC still rising (no sharp natural block count)" else ""))

relabel <- setNames(paste0("B", seq_along(sort(table(fit$classification), decreasing = TRUE))),
                     names(sort(table(fit$classification), decreasing = TRUE)))
block_all <- unname(relabel[as.character(fit$classification)])
n <- length(genes)
block_A7 <- block_all[1:n]
block_A9 <- block_all[(n + 1):(2 * n)]

blocks <- tibble(gene = genes, block_A7 = block_A7, block_A9 = block_A9)
write_tsv(blocks, file.path(res_dir, "06_latent_blocks.tsv"))

cat("\nBlock sizes (A7):\n");  print(sort(table(block_A7), decreasing = TRUE))
cat("\nBlock sizes (A9):\n");  print(sort(table(block_A9), decreasing = TRUE))
cat(sprintf("\nGenes that changed latent block (A7 -> A9): %d / %d (%.1f%%)\n",
            sum(block_A7 != block_A9), nrow(blocks), 100 * mean(block_A7 != block_A9)))
message("\n[OK] Saved: 06_latent_blocks.tsv")
