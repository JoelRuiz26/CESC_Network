#!/usr/bin/env Rscript
############################################################
# 04_pairwise_metrics.R
#
# The complete Figure-2 validation panel for the correlation networks:
# the four pairwise metrics already used for ARACNe, recomputed identically
# for Pearson and Spearman on the same leave-x-out subsamples.
#
#   Jaccard      edge overlap                       (05_bootstrap_null_validation.R)
#   Spearman rho degree-sequence similarity, over the FIXED 15,282-gene
#                universe with zeros included so absent genes still count
#   ARI          Louvain module-partition agreement on the shared nodes
#   KS D         two-sample KS statistic between the two degree vectors
#                (the D statistic as effect size, not its p-value -- degrees
#                 are not independent draws, so the p-value is meaningless at
#                 this many nodes; this is the reasoning in 05c)
#
# Note on direction: Jaccard, rho and ARI are SIMILARITIES (higher = more
# alike) while KS D is a DISTANCE (lower = more alike). The QAP statistic is
# always computed as inter - intra, so a clade effect gives a negative value
# for the three similarities and a positive one for KS D. The plot labels say
# which way each panel reads.
#
# ARACNe is recomputed here from its full consolidated networks rather than
# read from the old RDS, so all three methods get byte-identical treatment.

options(bitmapType = "cairo")
suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(igraph); library(parallel)
})

REPO <- path.expand("~/CESC_Network")
BASE <- file.path(REPO, "04_network_inference", "correlation_networks")
OUT  <- file.path(BASE, "pairwise_metrics"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
set.seed(42); NREP <- 100L; NCORE <- 12L
log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "")

GENES <- fread(file.path(BASE, "A7", "submatrices", "2_7_Full_counts_A7_r01.tsv"),
               sep = "\t", header = TRUE, select = 1L)[[1]]
NG <- length(GENES); gidx <- setNames(seq_len(NG), GENES)

edges_of <- function(method, cl, i) {
  if (method == "aracne") {
    e <- fread(sprintf("%s/networks_consolidated/%s/network_r%02d.txt", BASE, cl, i),
               sep = "\t", header = TRUE, select = c("Regulator", "Target"))
    a <- gidx[e$Regulator]; b <- gidx[e$Target]
  } else {
    dir <- if (method == "pearson") "networks_pearson" else "networks_spearman"
    k <- readRDS(sprintf("%s/%s/networks/%s/keys_r%02d.rds", BASE, dir, cl, i))
    b <- k %% (NG + 1); a <- (k - b) / (NG + 1)
  }
  ok <- !is.na(a) & !is.na(b) & a != b
  list(a = as.integer(a[ok]), b = as.integer(b[ok]))
}

ari <- function(x, y) {
  tab <- table(x, y); n <- sum(tab)
  sij <- sum(choose(tab, 2)); sa <- sum(choose(rowSums(tab), 2)); sb <- sum(choose(colSums(tab), 2))
  expct <- sa * sb / choose(n, 2); mx <- (sa + sb) / 2
  if (mx == expct) return(1); (sij - expct) / (mx - expct)
}

qap <- function(mat, labels, B = 5000) {
  same <- outer(labels, labels, `==`); diag(same) <- NA
  keep <- !is.na(mat)
  T_obs <- mean(mat[keep & !same], na.rm = TRUE) - mean(mat[keep & same], na.rm = TRUE)
  T_perm <- vapply(seq_len(B), function(b) {
    lp <- sample(labels); sm <- outer(lp, lp, `==`); diag(sm) <- NA
    mean(mat[keep & !sm], na.rm = TRUE) - mean(mat[keep & sm], na.rm = TRUE)
  }, numeric(1))
  up <- upper.tri(mat)
  iv <- mat[same & up]; ev <- mat[!same & up]
  pooled <- sqrt((var(iv, na.rm = TRUE) + var(ev, na.rm = TRUE)) / 2)
  list(p = (sum(abs(T_perm) >= abs(T_obs)) + 1) / (B + 1),
       d = (mean(ev, na.rm = TRUE) - mean(iv, na.rm = TRUE)) / pooled)
}

labs  <- c(sprintf("A7_r%02d", 1:NREP), sprintf("A9_r%02d", 1:NREP))
clade <- rep(c("A7", "A9"), each = NREP)
jobs  <- c(lapply(1:NREP, function(i) list(cl = "A7", i = i)),
           lapply(1:NREP, function(i) list(cl = "A9", i = i)))
NN <- length(jobs); IJ <- t(combn(NN, 2))

methods <- c(aracne = "ARACNe  (MI + DPI)", pearson = "Pearson  (Bonferroni)", spearman = "Spearman  (Bonferroni)")
## ARACNe is included only when its consolidated networks are on disk; the
## correlation panels do not depend on it.
if (!dir.exists(file.path(BASE, "networks_consolidated", "A9"))) {
  log_msg("networks_consolidated/ ausente -> solo Pearson y Spearman")
  methods <- methods[c("pearson", "spearman")]
}
all_pairs <- list(); all_stats <- list()

for (m in names(methods)) {
  log_msg("================ %s ================", m)
  ED <- lapply(jobs, function(j) edges_of(m, j$cl, j$i))
  log_msg("aristas/red: media %s", format(round(mean(sapply(ED, function(e) length(e$a)))), big.mark = ","))

  ## --- Jaccard (sparse incidence product) ---
  K <- lapply(ED, function(e) as.numeric(pmin(e$a, e$b)) * (NG + 1) + pmax(e$a, e$b))
  uni <- sort(unique(unlist(K)))
  Inc <- sparseMatrix(i = rep(seq_along(K), lengths(K)),
                      j = unlist(lapply(K, function(k) match(k, uni))), x = 1,
                      dims = c(NN, length(uni)))
  Inter <- as.matrix(Inc %*% Matrix::t(Inc)); de <- diag(Inter)
  JAC <- Inter / (outer(de, de, "+") - Inter); diag(JAC) <- NA
  rm(Inc, Inter, K, uni); gc(FALSE)
  log_msg("  Jaccard listo")

  ## --- degree matrix over the fixed gene universe ---
  DEG <- matrix(0L, NG, NN)
  for (i in seq_len(NN)) {
    tb <- tabulate(c(ED[[i]]$a, ED[[i]]$b), nbins = NG)
    DEG[, i] <- tb
  }
  SPD <- cor(DEG, method = "spearman"); diag(SPD) <- NA
  log_msg("  Spearman de grado listo")

  ## --- KS D between degree vectors ---
  KS <- matrix(NA_real_, NN, NN)
  ksv <- unlist(mclapply(seq_len(nrow(IJ)), function(r)
    suppressWarnings(ks.test(DEG[, IJ[r, 1]], DEG[, IJ[r, 2]])$statistic),
    mc.cores = NCORE))
  KS[IJ] <- ksv; KS[IJ[, c(2, 1)]] <- ksv
  log_msg("  KS D listo")

  ## --- Louvain + pairwise ARI ---
  ## The edge list MUST be character. With an integer matrix igraph creates
  ## vertices with no $name, intersect(names(mi), names(mj)) then returns
  ## nothing, every pair falls through the "<10 common nodes" guard and the
  ## whole ARI column comes out NaN.
  MEM <- mclapply(ED, function(e) {
    el <- cbind(as.character(e$a), as.character(e$b))
    g  <- simplify(graph_from_edgelist(el, directed = FALSE))
    mm <- cluster_louvain(g)$membership; names(mm) <- V(g)$name; mm
  }, mc.cores = NCORE)
  stopifnot(!is.null(names(MEM[[1]])))
  arv <- unlist(mclapply(seq_len(nrow(IJ)), function(r) {
    mi <- MEM[[IJ[r, 1]]]; mj <- MEM[[IJ[r, 2]]]
    cm <- intersect(names(mi), names(mj))
    if (length(cm) < 10) return(NA_real_)
    ari(mi[cm], mj[cm])
  }, mc.cores = NCORE))
  ARIm <- matrix(NA_real_, NN, NN); ARIm[IJ] <- arv; ARIm[IJ[, c(2, 1)]] <- arv
  stopifnot(is.numeric(arv), sum(!is.na(arv)) == nrow(IJ))
  log_msg("  ARI listo")
  rm(ED, DEG, MEM); gc(FALSE)

  ## --- long table + QAP per metric ---
  grp <- ifelse(clade[IJ[, 1]] == clade[IJ[, 2]],
                paste0("intra-", clade[IJ[, 1]]), "inter A7-A9")
  pd <- data.table(method = methods[[m]], method_key = m, group = grp,
                   Jaccard = JAC[IJ], Spearman = SPD[IJ], ARI = ARIm[IJ], KS = KS[IJ])
  all_pairs[[m]] <- pd

  for (met in c("Jaccard", "Spearman", "ARI", "KS")) {
    M <- switch(met, Jaccard = JAC, Spearman = SPD, ARI = ARIm, KS = KS)
    q <- qap(M, clade, 5000)
    all_stats[[paste(m, met)]] <- data.table(
      method = methods[[m]], method_key = m, metric = met,
      intra_A7 = mean(pd[group == "intra-A7"][[met]], na.rm = TRUE),
      intra_A9 = mean(pd[group == "intra-A9"][[met]], na.rm = TRUE),
      inter    = mean(pd[group == "inter A7-A9"][[met]], na.rm = TRUE),
      QAP_p = q$p, cohen_d = q$d)
    log_msg("  %-9s intraA7=%.4f intraA9=%.4f inter=%.4f | QAP p=%.4g d=%.2f",
            met, all_stats[[paste(m, met)]]$intra_A7, all_stats[[paste(m, met)]]$intra_A9,
            all_stats[[paste(m, met)]]$inter, q$p, q$d)
  }
  saveRDS(list(JAC = JAC, SPD = SPD, ARI = ARIm, KS = KS), file.path(OUT, sprintf("matrices_%s.rds", m)))
  rm(JAC, SPD, ARIm, KS); gc(FALSE)
}

PD <- rbindlist(all_pairs); ST <- rbindlist(all_stats)
fwrite(PD, file.path(OUT, "pairwise_metrics_long.csv.gz"), compress = "gzip")
fwrite(ST, file.path(OUT, "pairwise_metrics_summary.csv"))
log_msg("Done -> %s", OUT)
print(ST)
