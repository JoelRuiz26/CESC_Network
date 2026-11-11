#!/usr/bin/env Rscript
############################################################
# 03_spearman_networks.R
#
# The 200 leave-x-out networks rebuilt as Pearson relevance networks
# (Butte & Kohane 2000), from the SAME subsample matrices ARACNe used.
#
# Edge criterion: Bonferroni-corrected significance of the correlation.
#   t = r*sqrt(n-2)/sqrt(1-r^2),  df = n-2 = 58
#   alpha = 0.05 / 116,762,121 pairs  ->  |r| >= 0.7015
# The threshold depends only on n and the number of pairs, so it is identical
# for every network and for both clades; nothing about it is tuned to the
# data. Each network's edge count is therefore an output, not a setting.
#
# Spearman variant: genes are rank-transformed before correlating, so this is
# invariant to monotone transforms (the property MI has and Pearson lacks).
# Same Bonferroni threshold -- it depends only on n and the pair count.
#
# Phase 1 builds and saves the networks, phase 2 computes all 200x200 Jaccard
# values with one sparse crossprod, phase 3 summarises.

suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(parallel); library(RhpcBLASctl)
})
options(bitmapType = "cairo")

REPO <- path.expand("~/CESC_Network")
BASE <- file.path(REPO, "04_network_inference", "correlation_networks")
OUT  <- file.path(BASE, "networks_spearman")
dir.create(file.path(OUT, "networks", "A7"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "networks", "A9"), recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "")

NWORKER <- 6L; BLAS_PER <- 8L            # 6 x 8 = 48 cores, under the 70 cap
NREP <- 100L

sub_file <- function(cl, i)
  sprintf("%s/%s/submatrices/2_%s_Full_counts_%s_r%02d.tsv", BASE, cl,
          if (cl == "A7") "7" else "8", cl, i)

genes <- fread(sub_file("A7", 1), sep = "\t", header = TRUE, select = 1L)[[1]]
NG <- length(genes); NP <- NG * (NG - 1) / 2
NSAMP <- ncol(fread(sub_file("A7", 1), sep = "\t", header = TRUE, nrows = 1)) - 1L
THR <- { tq <- qt(1 - (0.05 / NP) / 2, NSAMP - 2); tq / sqrt(NSAMP - 2 + tq^2) }
log_msg("genes=%d  n=%d  pares=%s  ->  umbral Bonferroni |r| >= %.4f",
        NG, NSAMP, format(NP, big.mark = ","), THR)

## ---- phase 1: build ------------------------------------------------------

build_one <- function(job) {
  blas_set_num_threads(BLAS_PER); omp_set_num_threads(BLAS_PER)
  cl <- job$cl; i <- job$i
  dt <- fread(sub_file(cl, i), sep = "\t", header = TRUE)
  X  <- log2(as.matrix(dt[, -1, with = FALSE]) + 1); rm(dt); X <- t(apply(X, 1, rank))
  R  <- stats::cor(t(X)); rm(X); gc(FALSE)
  R[!is.finite(R)] <- 0
  R[lower.tri(R, diag = TRUE)] <- 0
  idx <- which(R >= THR | R <= -THR)
  ij  <- arrayInd(idx, dim(R)); r <- R[idx]; rm(R); gc(FALSE)
  o <- order(abs(r), decreasing = TRUE)
  ij <- ij[o, , drop = FALSE]; r <- r[o]
  fwrite(data.table(Regulator = genes[ij[, 1]], Target = genes[ij[, 2]], r = r),
         sprintf("%s/networks/%s/spearman_r%02d.txt", OUT, cl, i), sep = "\t")
  keys <- as.numeric(ij[, 1]) * (NG + 1) + ij[, 2]
  saveRDS(keys, sprintf("%s/networks/%s/keys_r%02d.rds", OUT, cl, i))
  data.table(clade = cl, rep = i, edges = length(keys),
             nodes = length(unique(c(ij[, 1], ij[, 2]))),
             min_abs_r = min(abs(r)), max_abs_r = max(abs(r)))
}

jobs <- c(lapply(seq_len(NREP), function(i) list(cl = "A7", i = i)),
          lapply(seq_len(NREP), function(i) list(cl = "A9", i = i)))
log_msg("construyendo %d redes con %d workers ...", length(jobs), NWORKER)
t0 <- Sys.time()
stats_dt <- rbindlist(mclapply(jobs, build_one, mc.cores = NWORKER, mc.preschedule = FALSE))
log_msg("listo en %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
fwrite(stats_dt, file.path(OUT, "network_stats.csv"))
print(stats_dt[, .(redes = .N, aristas_media = round(mean(edges)), aristas_min = min(edges),
                   aristas_max = max(edges), nodos_media = round(mean(nodes)),
                   pct_nodos = round(100 * mean(nodes) / NG, 1)), by = clade])

## ---- phase 2: all 200x200 Jaccard in one sparse crossprod ---------------

log_msg("ensamblando matriz dispersa de patrones ...")
labs <- c(sprintf("A7_r%02d", 1:NREP), sprintf("A9_r%02d", 1:NREP))
K <- lapply(jobs, function(j) readRDS(sprintf("%s/networks/%s/keys_r%02d.rds", OUT, j$cl, j$i)))
uni <- sort(unique(unlist(K)))
log_msg("universo de aristas: %s", format(length(uni), big.mark = ","))
M <- sparseMatrix(i = unlist(lapply(K, function(k) match(k, uni))),
                  j = rep(seq_along(K), lengths(K)),
                  x = 1, dims = c(length(uni), length(K)), dimnames = list(NULL, labs))
rm(K); gc(FALSE)
cnt <- Matrix::colSums(M)
I <- as.matrix(Matrix::crossprod(M)); rm(M); gc(FALSE)
U <- outer(cnt, cnt, "+") - I
J <- I / U; J[U == 0] <- NA_real_
dimnames(J) <- list(labs, labs)
saveRDS(J, file.path(OUT, "jaccard_matrix.rds"))

## ---- phase 3: summary ----------------------------------------------------

up <- upper.tri(J)
grp <- outer(substr(labs, 1, 2), substr(labs, 1, 2), paste, sep = "-")
summ <- data.table(
  comparacion = c("intra-A7", "intra-A9", "inter-clado"),
  n_pares = c(sum(up & grp == "A7-A7"), sum(up & grp == "A9-A9"), NREP * NREP),
  jaccard_medio = c(mean(J[up & grp == "A7-A7"]), mean(J[up & grp == "A9-A9"]),
                    mean(J[1:NREP, (NREP + 1):(2 * NREP)])),
  sd = c(sd(J[up & grp == "A7-A7"]), sd(J[up & grp == "A9-A9"]),
         sd(as.vector(J[1:NREP, (NREP + 1):(2 * NREP)]))))
log_msg("=== Jaccard, red Spearman/Bonferroni ===")
print(summ)
fwrite(summ, file.path(OUT, "jaccard_summary.csv"))

log_msg("Done -> %s", OUT)
