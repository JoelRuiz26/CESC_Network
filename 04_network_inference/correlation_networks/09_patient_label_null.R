#!/usr/bin/env Rscript
############################################################
# 46_patient_label_null.R
#
# The null the design actually needs, replacing the QAP of 44_qap_contrasts.R.
#
# WHY THE QAP IS THE WRONG TEST HERE. It permutes the clade label of the 200
# finished networks, which assumes the 200 are exchangeable under H0. They are
# not. A7 has 65 patients, so two A7 subsamples of 60 share 55-58 of them; two
# A9 subsamples share 8-29 of 202; cross-clade pairs share none, because the
# clades are disjoint patient sets. Even if A7 and A9 were biologically
# identical, intra-A7 Jaccard would sit near 0.70 and inter near 0.30 from
# patient overlap alone. Permuting network labels destroys that structure, so
# the permutation null is far too flat and the P value hits the 1/(B+1) floor
# whatever the data say. It tests the subsampling design, not the clade.
#
# THE CORRECT NULL. Permute the clade label of the PATIENTS. Partition the 267
# tumours at random into a pseudo-A7 of 65 and a pseudo-A9 of 202 -- the real
# sizes -- then redo the leave-x-out subsampling, rebuild the networks and
# recompute the four metrics. Overlap structure is preserved exactly (60 of 65
# against 60 of 202, disjoint pseudo-clades); only clade identity is destroyed.
# This cannot reuse the existing networks: a network over a mixed set of 60
# patients has to be computed, the similarity between two networks depends on
# which patients went into each.
#
# Cost is why B is 100 and not 20,000: each network is a 15,282 x 15,282
# rank correlation, 23.7 s with 7 BLAS threads. The floor is therefore
# P = 1/101 = 0.0099. A P at that floor under an honest null says more than
# 5e-5 under a null that was never valid.
#
# Everything else is held identical to the real networks: same 15,282 genes,
# same n = 60, same Bonferroni threshold |rho| >= 0.701528 (it depends only on
# n and the pair count, so it does not drift between permutations), same
# Spearman ranks, same four metrics, same two contrasts.
#
# The observed side reuses 20 of the existing A7 and 20 of the A9 networks:
# the null is built at 20 + 20, and the observed statistic has to come from the
# same design or the null variance is wrong.

suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(igraph)
  library(matrixStats); library(parallel); library(RhpcBLASctl)
})

BASE <- "/STORAGE/csbig/jruiz/cesc_network/ARACNE_AP/Leave_x_out"
OUT  <- file.path(BASE, "46_patient_null"); dir.create(OUT, showWarnings = FALSE)
## args: <perm_from> <perm_to> <workers>.  The permutations are seeded
## individually (set.seed(10000 + b)), so a range can be run in a separate
## process and appended without changing any result.
.a <- commandArgs(TRUE)
NPERM   <- 100L                          # P floor = 1/101 = 0.0099
PFROM   <- if (length(.a) >= 1) as.integer(.a[1]) else 1L
PTO     <- if (length(.a) >= 2) as.integer(.a[2]) else NPERM
NWORKER <- if (length(.a) >= 3) as.integer(.a[3]) else 10L
BLAS_PER <- max(1L, as.integer(floor(70 / NWORKER)))
MCORE   <- 8L                            # for the pairwise metric loops
NREP    <- 20L                           # networks per pseudo-clade
N7      <- 65L; N9 <- 202L; NSAMP <- 60L
set.seed(2026)
log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "")

## ---- data ---------------------------------------------------------------

A <- fread(file.path(BASE, "2_7_Full_counts_A7_annot.tsv"))
B9 <- fread(file.path(BASE, "2_8_Full_counts_A9_annot.tsv"))
stopifnot(identical(A[[1]], B9[[1]]))
GENES <- A[[1]]; NG <- length(GENES)
X <- cbind(as.matrix(A[, -1, with = FALSE]), as.matrix(B9[, -1, with = FALSE]))
rm(A, B9); gc(FALSE)
stopifnot(ncol(X) == N7 + N9)
NP  <- NG * (NG - 1) / 2
THR <- { tq <- qt(1 - (0.05 / NP) / 2, NSAMP - 2); tq / sqrt(NSAMP - 2 + tq^2) }
log_msg("genes=%d  pacientes=%d  n=%d  umbral |rho| >= %.6f", NG, ncol(X), NSAMP, THR)
log_msg("permutaciones %d..%d  |  %d workers x %d hilos BLAS", PFROM, PTO, NWORKER, BLAS_PER)

## ---- one network from a set of patient columns --------------------------

net_keys <- function(cols) {
  blas_set_num_threads(BLAS_PER); omp_set_num_threads(BLAS_PER)
  Y <- matrixStats::rowRanks(X[, cols, drop = FALSE], ties.method = "average")
  R <- stats::cor(t(Y)); rm(Y)
  R[!is.finite(R)] <- 0
  R[lower.tri(R, diag = TRUE)] <- 0
  idx <- which(R >= THR | R <= -THR); rm(R); gc(FALSE)
  ij <- arrayInd(idx, c(NG, NG))
  as.numeric(ij[, 1]) * (NG + 1) + ij[, 2]
}

## ---- the four metrics and the two contrasts, from a list of key vectors --

ari2 <- function(x, y) {
  tab <- table(x, y); n <- sum(tab)
  sij <- sum(choose(tab, 2)); sa <- sum(choose(rowSums(tab), 2)); sb <- sum(choose(colSums(tab), 2))
  expct <- sa * sb / choose(n, 2); mx <- (sa + sb) / 2
  if (mx == expct) return(1); (sij - expct) / (mx - expct)
}

contrasts_from_keys <- function(K, grp) {
  NN <- length(K); IJ <- t(combn(NN, 2))
  ## Jaccard
  uni <- sort(unique(unlist(K)))
  Inc <- sparseMatrix(i = rep(seq_along(K), lengths(K)),
                      j = unlist(lapply(K, function(k) match(k, uni))), x = 1,
                      dims = c(NN, length(uni)))
  Int <- as.matrix(Inc %*% Matrix::t(Inc)); de <- diag(Int)
  JAC <- Int / (outer(de, de, "+") - Int); diag(JAC) <- NA
  rm(Inc, Int, uni); gc(FALSE)
  ## degrees over the fixed gene universe, zeros kept
  DEG <- matrix(0L, NG, NN)
  EL <- lapply(K, function(k) { b <- k %% (NG + 1); a <- (k - b) / (NG + 1)
                                list(a = as.integer(a), b = as.integer(b)) })
  for (i in seq_len(NN)) DEG[, i] <- tabulate(c(EL[[i]]$a, EL[[i]]$b), nbins = NG)
  SPD <- cor(DEG, method = "spearman"); diag(SPD) <- NA
  ## KS D
  KS <- matrix(NA_real_, NN, NN)
  v <- unlist(mclapply(seq_len(nrow(IJ)), function(r)
    suppressWarnings(ks.test(DEG[, IJ[r, 1]], DEG[, IJ[r, 2]])$statistic),
    mc.cores = MCORE))
  KS[IJ] <- v; KS[IJ[, c(2, 1)]] <- v
  ## Louvain + ARI.  character edgelist so the vertices carry names
  MEM <- lapply(EL, function(e) {
    g <- simplify(graph_from_edgelist(cbind(as.character(e$a), as.character(e$b)),
                                      directed = FALSE))
    mm <- cluster_louvain(g)$membership; names(mm) <- V(g)$name; mm
  })
  v <- unlist(mclapply(seq_len(nrow(IJ)), function(r) {
    mi <- MEM[[IJ[r, 1]]]; mj <- MEM[[IJ[r, 2]]]
    cm <- intersect(names(mi), names(mj))
    if (length(cm) < 10) return(NA_real_)
    ari2(mi[cm], mj[cm])
  }, mc.cores = MCORE))
  ARIm <- matrix(NA_real_, NN, NN); ARIm[IJ] <- v; ARIm[IJ[, c(2, 1)]] <- v
  rm(DEG, EL, MEM); gc(FALSE)
  ## the two contrasts: inter - intra, exactly as in 44_qap_contrasts.R
  up <- upper.tri(matrix(0, NN, NN))
  s7 <- up & outer(grp, grp, function(a, b) a == "g7" & b == "g7")
  s9 <- up & outer(grp, grp, function(a, b) a == "g9" & b == "g9")
  si <- up & outer(grp, grp, `!=`)
  out <- list()
  for (nm in c("Jaccard", "Spearman", "ARI", "KS")) {
    M <- switch(nm, Jaccard = JAC, Spearman = SPD, ARI = ARIm, KS = KS)
    out[[nm]] <- c(g7 = mean(M[si], na.rm = TRUE) - mean(M[s7], na.rm = TRUE),
                   g9 = mean(M[si], na.rm = TRUE) - mean(M[s9], na.rm = TRUE),
                   m7 = mean(M[s7], na.rm = TRUE), m9 = mean(M[s9], na.rm = TRUE),
                   mi = mean(M[si], na.rm = TRUE))
  }
  out
}

flat <- function(o, tag, b) rbindlist(lapply(names(o), function(nm)
  data.table(source = tag, perm = b, metric = nm,
             T_g7 = o[[nm]][["g7"]], T_g9 = o[[nm]][["g9"]],
             mean_g7 = o[[nm]][["m7"]], mean_g9 = o[[nm]][["m9"]],
             mean_inter = o[[nm]][["mi"]])))

## ---- observed, from the existing networks, same 20 + 20 design ----------

grp <- rep(c("g7", "g9"), each = NREP)
if (PFROM > 1L && file.exists(file.path(OUT, "observed.csv"))) {
  OBS <- fread(file.path(OUT, "observed.csv"))
  log_msg("observado leido de observed.csv")
} else {
K_obs <- c(lapply(seq_len(NREP), function(i)
             readRDS(sprintf("%s/27_spearman_networks/networks/A7/keys_r%02d.rds", BASE, i))),
           lapply(seq_len(NREP), function(i)
             readRDS(sprintf("%s/27_spearman_networks/networks/A9/keys_r%02d.rds", BASE, i))))
log_msg("observado: %d redes, aristas media %s", length(K_obs),
        format(round(mean(lengths(K_obs))), big.mark = ","))
OBS <- flat(contrasts_from_keys(K_obs, grp), "observed", 0L)
fwrite(OBS, file.path(OUT, "observed.csv")); print(OBS)
rm(K_obs); gc(FALSE)
}

## ---- the null -----------------------------------------------------------

nullfile <- file.path(OUT, "null_stats.csv")
for (b in seq(PFROM, PTO)) {
  t0 <- Sys.time()
  set.seed(10000 + b)
  pp <- sample(ncol(X))
  p7 <- pp[1:N7]; p9 <- pp[(N7 + 1):(N7 + N9)]
  jobs <- c(lapply(seq_len(NREP), function(i) sample(p7, NSAMP)),
            lapply(seq_len(NREP), function(i) sample(p9, NSAMP)))
  K <- mclapply(jobs, net_keys, mc.cores = NWORKER, mc.preschedule = FALSE)
  bad <- !vapply(K, is.numeric, logical(1))
  if (any(bad)) { log_msg("perm %d: %d redes fallaron, se omite", b, sum(bad)); next }
  st <- flat(contrasts_from_keys(K, grp), "null", b)
  fwrite(st, nullfile, append = file.exists(nullfile))
  log_msg("perm %3d/%d  (%.1f min)  aristas media %s | T_g7 Jaccard %+.4f  T_g9 %+.4f",
          b, NPERM, as.numeric(difftime(Sys.time(), t0, units = "mins")),
          format(round(mean(lengths(K))), big.mark = ","),
          st[metric == "Jaccard", T_g7], st[metric == "Jaccard", T_g9])
  rm(K); gc(FALSE)
}

## ---- P values -----------------------------------------------------------

NU <- fread(nullfile)
res <- rbindlist(lapply(c("Jaccard", "Spearman", "ARI", "KS"), function(nm) {
  o <- OBS[metric == nm]; n <- NU[metric == nm]
  data.table(metric = nm, B = nrow(n),
             obs_g7 = o$T_g7, null_mean_g7 = mean(n$T_g7), null_sd_g7 = sd(n$T_g7),
             n_ge_g7 = sum(abs(n$T_g7) >= abs(o$T_g7)),
             p_g7 = (sum(abs(n$T_g7) >= abs(o$T_g7)) + 1) / (nrow(n) + 1),
             ## signed: (obs - null mean)/null sd.  Taking abs(obs) against a
             ## signed null mean mixes the signs and inflates z wildly for the
             ## metrics whose statistic is negative.
             z_g7 = (o$T_g7 - mean(n$T_g7)) / sd(n$T_g7),
             obs_g9 = o$T_g9, null_mean_g9 = mean(n$T_g9), null_sd_g9 = sd(n$T_g9),
             n_ge_g9 = sum(abs(n$T_g9) >= abs(o$T_g9)),
             p_g9 = (sum(abs(n$T_g9) >= abs(o$T_g9)) + 1) / (nrow(n) + 1),
             z_g9 = (o$T_g9 - mean(n$T_g9)) / sd(n$T_g9))
}))
fwrite(res, file.path(OUT, "patient_null_pvalues.csv"))
print(res)
log_msg("Done -> %s", OUT)
