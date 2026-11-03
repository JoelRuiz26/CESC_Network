### ============================================================
### 02_permutation_significance.R
###
### Is the observed A7-vs-A9 latent distance smaller than you'd get from
### two UNRELATED gene embeddings? A full parametric graph-resampling
### bootstrap (the test proposed in Tang et al. 2017, doi:10.1080/
### 10618600.2016.1193505) needs to re-simulate whole graphs from a dense
### n x n Bernoulli probability matrix many times; at n ~ 15,000 genes that
### is not computationally tractable (that test was validated on
### connectomes of a few hundred nodes). Used here instead: a cheaper,
### honest vertex-label permutation test on the ALREADY-COMPUTED embeddings
### from 01_ase_procrustes_distance.R. It answers a related but narrower
### question: "do genes keep their relative latent position across clades
### more than if gene identity were scrambled?" -- not "are the two whole
### graphs significantly different" in the full Tang et al. sense.
###
### Null: permute WHICH gene's A9 embedding is matched to which gene's A7
### embedding, redo the Procrustes alignment for that scrambled pairing,
### and recompute the global distance. Repeat B times (embarrassingly
### parallel -> parallel::mclapply, mc.cores = 9, leaving 3 of 12 free).
### One-sided p-value: P(null mean distance <= observed mean distance).
###
### -> results/02_permutation_null.tsv   (B null statistics)
### -> results/02_permutation_pvalue.tsv (observed stat, p-value)
### ============================================================
suppressPackageStartupMessages({ library(parallel); library(dplyr); library(readr) })

set.seed(123)
out_dir <- "~/CESC_Network/11_latent_embedding/results"
B <- 2000
N_CORES <- 9   # 12 detected on this machine; leaves 3 free

X7  <- readRDS(file.path(out_dir, "01_ase_X7.rds"))
X9  <- readRDS(file.path(out_dir, "01_ase_X9_raw.rds"))
obs <- read_tsv(file.path(out_dir, "01_global_distance_summary.tsv"), show_col_types = FALSE)

procrustes_align <- function(X_source, X_target) {
  sv <- svd(crossprod(X_source, X_target))
  X_source %*% (sv$u %*% t(sv$v))
}

one_perm <- function(i) {
  perm <- sample.int(nrow(X9))
  X9_perm_aligned <- procrustes_align(X9[perm, ], X7)
  mean(sqrt(rowSums((X7 - X9_perm_aligned)^2)))
}

cat(sprintf("Running %d permutations on %d cores...\n", B, N_CORES))
null_stats <- unlist(mclapply(seq_len(B), one_perm, mc.cores = N_CORES))

p_value <- (1 + sum(null_stats <= obs$mean_distance)) / (B + 1)

write_tsv(tibble(permutation = seq_len(B), mean_distance = null_stats),
          file.path(out_dir, "02_permutation_null.tsv"))
result <- tibble(observed_mean_distance = obs$mean_distance,
                  null_mean = mean(null_stats), null_sd = sd(null_stats),
                  B = B, p_value = p_value)
write_tsv(result, file.path(out_dir, "02_permutation_pvalue.tsv"))

cat("\n== Permutation test: observed vs scrambled gene-identity null ==\n")
print(as.data.frame(result))
message("\n[OK] Saved: 02_permutation_null.tsv, 02_permutation_pvalue.tsv")
