#!/usr/bin/env Rscript
############################################################
# 06_qap_contrasts.R
#
# The QAP test restated as the TWO contrasts the design actually asks about:
#
#     C1   A7:A7  vs  A7:A9
#     C2   A9:A9  vs  A7:A9
#
# and NOT A7:A7 vs A9:A9. Those two intra groups are not on equal footing:
# two A7 subsamples share at least 55 of 60 patients (92%) because A7 has only
# 65, while two A9 subsamples share 12-29 of 60. Any difference between them is
# dominated by that, so the earlier single statistic -- inter minus the pooled
# intra -- invited exactly the misreading a reviewer would flag.
#
# Test. The clade label of each of the 200 networks is permuted B = 5000 times
# keeping the 100/100 split and each contrast is recomputed from the same fixed
# pairwise matrix. The permutation unit is the NETWORK, not the pair: the
# 19,900 pairwise values are not independent.
#
# Effect size is Cohen's d between the two groups of that contrast only, with
# the pooled s.d. of those two groups.

suppressPackageStartupMessages({ library(data.table) })
REPO <- path.expand("~/CESC_Network")
BASE <- file.path(REPO, "04_network_inference", "correlation_networks")
OUT  <- file.path(BASE, "figure"); dir.create(OUT, showWarnings = FALSE)
set.seed(42); NREP <- 100L; B <- 5000L
clade <- rep(c("A7", "A9"), each = NREP)

contrasts_qap <- function(M, labels, B) {
  up <- upper.tri(M)
  sel <- function(lb, what) {
    o7 <- outer(lb, lb, function(a, b) a == "A7" & b == "A7")
    o9 <- outer(lb, lb, function(a, b) a == "A9" & b == "A9")
    oi <- outer(lb, lb, `!=`)
    switch(what, A7 = up & o7, A9 = up & o9, inter = up & oi)
  }
  stat <- function(lb, g) {
    mean(M[sel(lb, "inter")], na.rm = TRUE) - mean(M[sel(lb, g)], na.rm = TRUE)
  }
  obs <- c(A7 = stat(labels, "A7"), A9 = stat(labels, "A9"))
  perm <- vapply(seq_len(B), function(b) {
    lp <- sample(labels); c(stat(lp, "A7"), stat(lp, "A9"))
  }, numeric(2))
  ev <- M[sel(labels, "inter")]
  out <- rbindlist(lapply(c("A7", "A9"), function(g) {
    iv <- M[sel(labels, g)]
    pooled <- sqrt((var(iv, na.rm = TRUE) + var(ev, na.rm = TRUE)) / 2)
    i <- if (g == "A7") 1L else 2L
    data.table(contrast = sprintf("%s:%s vs A7:A9", g, g),
               n_intra = sum(!is.na(iv)), n_inter = sum(!is.na(ev)),
               mean_intra = mean(iv, na.rm = TRUE), mean_inter = mean(ev, na.rm = TRUE),
               diff = obs[[i]],
               QAP_p = (sum(abs(perm[i, ]) >= abs(obs[[i]])) + 1) / (B + 1),
               cohen_d = (mean(ev, na.rm = TRUE) - mean(iv, na.rm = TRUE)) / pooled)
  }))
  out
}

res <- list()
for (meth in c("spearman", "pearson")) {
  MM <- readRDS(file.path(BASE, "pairwise_metrics", sprintf("matrices_%s.rds", meth)))
  for (met in c("Jaccard", "Spearman", "ARI", "KS")) {
    M <- MM[[switch(met, Jaccard = "JAC", Spearman = "SPD", ARI = "ARI", KS = "KS")]]
    r <- contrasts_qap(M, clade, B)
    r[, `:=`(method = meth, metric = met)]
    res[[paste(meth, met)]] <- r
    cat(sprintf("[%s] %-9s %s\n", meth, met,
        paste(sprintf("%s: p=%.4g d=%+.2f", r$contrast, r$QAP_p, r$cohen_d), collapse = " | ")))
  }
}
R <- rbindlist(res)[, .(method, metric, contrast, n_intra, n_inter,
                        mean_intra, mean_inter, diff, QAP_p, cohen_d)]
fwrite(R, file.path(OUT, "qap_contrasts.csv"))
cat("\n-> ", file.path(OUT, "qap_contrasts.csv"), "\n")
