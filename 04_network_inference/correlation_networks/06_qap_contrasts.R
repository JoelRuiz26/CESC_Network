#!/usr/bin/env Rscript
############################################################
# 44_qap_contrasts.R
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
# Test. The clade label of each of the 200 networks is permuted B = 20000 times
# keeping the 100/100 split and each contrast is recomputed from the same fixed
# pairwise matrix. The permutation unit is the NETWORK, not the pair: the
# 19,900 pairwise values are not independent.
#
# Effect size is Cohen's d between the two groups of that contrast only, with
# the pooled s.d. of those two groups.
#
# On the P value. A permutation test cannot return anything below 1/(B+1): if
# no permutation reaches the observed difference the answer is the floor, and
# raising B only lowers the floor. So the table reports the exact empirical
# value AND, for the cases that sit at the floor, how far out the observation
# is: n_ge (permutations at least as extreme), perm_mean / perm_sd of the null,
# and perm_z = (|obs| - mean) / sd. The z is the informative number when the
# P value is resolution-limited; p_gauss is its normal-tail equivalent, given
# only as a reference and never as the primary statistic.

suppressPackageStartupMessages({ library(data.table) })
BASE <- "/STORAGE/csbig/jruiz/cesc_network/ARACNE_AP/Leave_x_out"
OUT  <- file.path(BASE, "43_figure_paper"); dir.create(OUT, showWarnings = FALSE)
set.seed(42); NREP <- 100L; B <- 20000L
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
               n_ge = sum(abs(perm[i, ]) >= abs(obs[[i]])),
               QAP_p = (sum(abs(perm[i, ]) >= abs(obs[[i]])) + 1) / (B + 1),
               perm_mean = mean(perm[i, ]), perm_sd = sd(perm[i, ]),
               perm_z = (abs(obs[[i]]) - mean(perm[i, ])) / sd(perm[i, ]),
               p_gauss = 2 * pnorm((abs(obs[[i]]) - mean(perm[i, ])) / sd(perm[i, ]),
                                   lower.tail = FALSE),
               cohen_d = (mean(ev, na.rm = TRUE) - mean(iv, na.rm = TRUE)) / pooled)
  }))
  out
}

res <- list()
for (meth in c("spearman", "pearson")) {
  MM <- readRDS(file.path(BASE, "29_full_panel", sprintf("matrices_%s.rds", meth)))
  for (met in c("Jaccard", "Spearman", "ARI", "KS")) {
    M <- MM[[switch(met, Jaccard = "JAC", Spearman = "SPD", ARI = "ARI", KS = "KS")]]
    r <- contrasts_qap(M, clade, B)
    r[, `:=`(method = meth, metric = met)]
    res[[paste(meth, met)]] <- r
    cat(sprintf("[%s] %-9s %s\n", meth, met,
        paste(sprintf("%s: p=%.4g (%d/%d) z=%.1f d=%+.2f",
                      r$contrast, r$QAP_p, r$n_ge, B, r$perm_z, r$cohen_d),
              collapse = " | ")))
  }
}
R <- rbindlist(res)[, .(method, metric, contrast, n_intra, n_inter,
                        mean_intra, mean_inter, diff, n_ge, B_perm = B, QAP_p,
                        perm_mean, perm_sd, perm_z, p_gauss, cohen_d)]
fwrite(R, file.path(OUT, "qap_contrasts.csv"))
cat("\n-> ", file.path(OUT, "qap_contrasts.csv"), "\n")
