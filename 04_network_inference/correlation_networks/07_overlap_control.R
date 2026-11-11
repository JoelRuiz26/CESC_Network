#!/usr/bin/env Rscript
############################################################
# 07_overlap_control.R
#
# Does the intra-vs-inter gap survive the patient-overlap confound?
#
# Inter-clade pairs share ZERO patients by construction: A7 and A9 are disjoint
# patient sets. Intra-clade pairs always share some. So "intra more similar
# than inter" is partly guaranteed by design, and the permutation test does not
# remove it -- permuting labels destroys clade structure but not the fact that
# same-clade networks were drawn from a common patient pool.
#
# The control: within A9, the number of shared patients varies across pairs
# (8-29 of 60). Regress the metric on shared count and extrapolate to 0 shared.
# If the extrapolated intra-A9 value still exceeds inter, the clade effect is
# not an overlap artefact.
#
# A7 cannot be controlled this way and this script does not pretend otherwise:
# its pairs share 55-58 of 60, so extrapolating to 0 is a 55-patient
# extrapolation outside the observed range. A7 is reported for completeness
# with its slope, never with an extrapolated value.

suppressPackageStartupMessages({ library(data.table) })
REPO <- path.expand("~/CESC_Network")
BASE <- file.path(REPO, "04_network_inference", "correlation_networks")
OUT  <- file.path(BASE, "figure"); dir.create(OUT, showWarnings = FALSE)
NREP <- 100L
clade <- rep(c("A7", "A9"), each = NREP)

sub_file <- function(cl, i)
  sprintf("%s/%s/submatrices/2_%s_Full_counts_%s_r%02d.tsv", BASE, cl,
          if (cl == "A7") "7" else "8", cl, i)
samples_of <- function(cl, i) strsplit(readLines(sub_file(cl, i), n = 1), "\t")[[1]][-1]
S <- list(A7 = lapply(1:NREP, function(i) samples_of("A7", i)),
          A9 = lapply(1:NREP, function(i) samples_of("A9", i)))

res <- list()
for (meth in c("spearman", "pearson")) {
  MM <- readRDS(file.path(BASE, "pairwise_metrics", sprintf("matrices_%s.rds", meth)))
  for (met in c("Jaccard", "Spearman", "ARI", "KS")) {
    M <- MM[[switch(met, Jaccard = "JAC", Spearman = "SPD", ARI = "ARI", KS = "KS")]]
    inter <- mean(M[1:NREP, (NREP + 1):(2 * NREP)], na.rm = TRUE)
    for (cl in c("A7", "A9")) {
      ij <- t(combn(NREP, 2)); off <- if (cl == "A7") 0L else NREP
      ov <- apply(ij, 1, function(p) length(intersect(S[[cl]][[p[1]]], S[[cl]][[p[2]]])))
      vv <- M[cbind(off + ij[, 1], off + ij[, 2])]
      ok <- !is.na(vv); ov <- ov[ok]; vv <- vv[ok]
      fit <- lm(vv ~ ov); ci <- suppressWarnings(predict(fit, data.frame(ov = 0), interval = "confidence"))
      extrapolable <- min(ov) <= 30
      res[[paste(meth, met, cl)]] <- data.table(
        method = meth, metric = met, clade = cl,
        shared_min = min(ov), shared_median = median(ov), shared_max = max(ov),
        intra_observed = mean(vv), slope_per_patient = coef(fit)[2],
        intra_at_zero = if (extrapolable) ci[1] else NA_real_,
        ci_low = if (extrapolable) ci[2] else NA_real_,
        ci_high = if (extrapolable) ci[3] else NA_real_,
        inter = inter,
        survives = if (extrapolable) {
          if (met == "KS") ci[3] < inter else ci[2] > inter
        } else NA)
    }
  }
}
R <- rbindlist(res)
fwrite(R, file.path(OUT, "overlap_control.csv"))
cat("\n=== A9: la unica donde el control es posible (solape 8-29 de 60) ===\n")
print(R[clade == "A9", .(method, metric, intra_observed = round(intra_observed, 4),
                         intra_at_zero = round(intra_at_zero, 4),
                         ci = sprintf("[%.4f, %.4f]", ci_low, ci_high),
                         inter = round(inter, 4), survives)])
cat("\n=== A7: solape 55-58 de 60, no extrapolable ===\n")
print(R[clade == "A7", .(method, metric, shared_min, shared_max,
                         intra_observed = round(intra_observed, 4),
                         slope = signif(slope_per_patient, 3), inter = round(inter, 4))])
