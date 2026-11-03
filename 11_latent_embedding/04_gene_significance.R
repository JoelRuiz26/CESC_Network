### ============================================================
### 04_gene_significance.R
###
### Per-gene significance for the A7-vs-A9 latent displacement
### (X7 - X9_aligned), beyond the global permutation test in 02_ (which
### already shows the identity-matched alignment fits far better than
### chance, p=0.0005 -- i.e. MOST genes show no true displacement, so the
### bulk of genes can serve as an empirical null, Efron-style: Efron B
### (2004) "Large-scale simultaneous hypothesis testing: the choice of a
### null hypothesis." JASA 99(465):96-104. doi:10.1198/016214504000000089).
###
### FIRST ATTEMPT (discarded, kept here as a documented negative result):
### fit one global robust covariance (MCD) to the displacement vectors and
### use a Mahalanobis-distance / chi-squared(df=d) test. That assumes
### homoskedastic, roughly Gaussian noise -- checking it showed it is
### badly wrong here: distance vs. gene degree has Spearman rho ~ 0.62-0.66
### (ASE embedding noise is well known to scale with a node's own degree --
### ties to the CLT covariance in Athreya et al. 2016 "A Limit Theorem for
### Scaled Eigenvectors of Random Dot Product Graphs", Sankhya A,
### doi:10.1007/s13171-015-0071-x), so a single pooled covariance massively
### over-flagged high-degree hub genes as "significant" just for being
### hubs (42% of all genes at FDR<0.05 -- not a believable result).
###
### USED INSTEAD: a degree-stratified empirical null. Genes are binned into
### 20 degree-matched strata (by min(degree_A7, degree_A9), the noisier of
### the two embeddings), and each gene's p-value is its empirical rank
### among ONLY the genes in its own stratum -- no distributional
### assumption, and it controls for the degree confound directly instead of
### assuming it away. BH-adjusted across all genes.
###
### Also cross-checked against the Louvain modules already computed in
### 06_modules: does a significant gene actually change module (not just
### continuous latent position) between A7 and A9?
###
### -> results/04_gene_significance.tsv
### ============================================================
suppressPackageStartupMessages({ library(igraph); library(dplyr); library(readr) })

res_dir  <- "~/CESC_Network/11_latent_embedding/results"
topo_dir <- "~/CESC_Network/05_topology/results"
mod_dir  <- "~/CESC_Network/06_modules/results"
ALPHA    <- 0.01
N_BINS   <- 20

X7         <- readRDS(file.path(res_dir, "01_ase_X7.rds"))
X9_aligned <- readRDS(file.path(res_dir, "01_ase_X9_aligned.rds"))
genes      <- read_tsv(file.path(res_dir, "01_shared_gene_order.tsv"), show_col_types = FALSE)$gene
distance   <- sqrt(rowSums((X7 - X9_aligned)^2))

## -- Degree-matched strata (the noisier/lower-degree side of each gene) --
g7 <- readRDS(file.path(topo_dir, "5_2_4_Subnetwork_A7_max.rds"))
g9 <- readRDS(file.path(topo_dir, "5_2_4_Subnetwork_A9_max.rds"))
min_degree <- pmin(degree(g7)[genes], degree(g9)[genes])
cat(sprintf("Spearman(distance, min degree) = %.3f -- confirms ASE noise scales with degree, so a single\n", cor(distance, min_degree, method = "spearman")))
cat("pooled covariance/chi2 test would just flag hubs; stratifying by degree instead.\n")

strata <- cut(min_degree, breaks = quantile(min_degree, probs = seq(0, 1, length.out = N_BINS + 1)),
              include.lowest = TRUE, labels = FALSE)

## -- Degree-standardised residual, then a parametric upper-tail p-value --
## SECOND failed attempt, also discarded as a documented negative result:
## ranking each gene's z-score among the SAME 14,998 genes being tested and
## calling that rank/(n+1) a "p-value". That self-referential construction
## has a hard floor of 1/(n+1), so even the single most extreme gene gets
## BH-adjusted back to ~1 (p_(1) x m/1 = [1/(n+1)] x n is ~1 by
## construction, whatever the real effect size) -- it needs an external
## reference distribution, not a rank among the tested items themselves.
##
## USED INSTEAD: the per-stratum median/MAD z-score is, by the classical
## definition of MAD (R's default constant=1.4826 makes it a consistent SD
## estimator under normality), already calibrated so the BULK of each
## stratum approximates z ~ N(0,1). A one-sided upper-tail p-value from
## that reference (pnorm), THEN BH-adjusted across all 14,998 genes, fixes
## the resolution problem (pnorm has no floor) while still controlling for
## the degree confound via the per-stratum standardisation.
sig <- tibble(gene = genes, distance = distance, min_degree = min_degree, stratum = strata) |>
  group_by(stratum) |>
  mutate(stratum_median = median(distance), stratum_mad = mad(distance),
         z = (distance - stratum_median) / stratum_mad, n_in_stratum = n()) |>
  ungroup() |>
  mutate(p_value = pnorm(z, lower.tail = FALSE),   # one-sided: large positive z = unusually far
         p_adj = p.adjust(p_value, method = "BH")) |>
  arrange(p_adj, desc(distance))

## -- Cross-check against Louvain module membership (06_modules) --
gm7 <- readRDS(file.path(mod_dir, "5_1_1_Louvain_A7_graph_max.rds"))
gm9 <- readRDS(file.path(mod_dir, "5_1_1_Louvain_A9_graph_max.rds"))
mod7 <- setNames(V(gm7)$Community_name, V(gm7)$name)
mod9 <- setNames(V(gm9)$Community_name, V(gm9)$name)
sig <- sig |> mutate(module_A7 = mod7[gene], module_A9 = mod9[gene],
                      module_changed = !is.na(module_A7) & !is.na(module_A9) & module_A7 != module_A9)

write_tsv(sig, file.path(res_dir, "04_gene_significance.tsv"))

n_sig <- sum(sig$p_adj < ALPHA)
cat(sprintf("\nDegree-stratified empirical test (%d strata), BH-adjusted.\n", N_BINS))
cat(sprintf("Significant genes at FDR < %.2f: %d / %d (%.1f%%)\n", ALPHA, n_sig, nrow(sig), 100 * n_sig / nrow(sig)))
cat(sprintf("Of those, module membership changed (Louvain, A7 vs A9): %d / %d\n",
            sum(sig$module_changed[sig$p_adj < ALPHA], na.rm = TRUE), n_sig))
cat("\nTop 15 significant genes:\n")
print(as.data.frame(head(sig, 15) |> dplyr::select(gene, distance, min_degree, z, p_adj, module_A7, module_A9, module_changed)))
message("\n[OK] Saved: 04_gene_significance.tsv")
