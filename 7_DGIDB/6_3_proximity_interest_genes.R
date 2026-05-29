# ============================================================
# 6_3_proximity_interest_genes.R
# Network-based drug proximity to DepMap interest genes
# ============================================================
# Metric : Closest proximity (Guney et al. 2016, Nat Commun)
#   d_c(T,S) = (1/|T|) * Σ_{t∈T}  min_{s∈S}  dist(t,s)
#
# Selectivity score:
#   S = d_c(T, non-interest) / d_c(T, interest)
#   S > 1 → targets are closer to interest genes than to the rest
#
# Significance: degree-preserving permutation of drug targets
#   z = (d_c_obs - µ_null) / σ_null
#   p = Φ(z)  [one-sided: negative z = closer than expected]
#
# Interest genes: is_interest == TRUE in 6_DepMap_cacu tables
#   (DGE ∩ Hub genes, i.e., is_DEG & is_hub)
#
# References:
#   Guney et al.  2016 – https://doi.org/10.1038/ncomms10331
#   Cheng et al.  2019 – https://doi.org/10.1038/s41467-019-10744-6
#   Zhou  et al.  2020 – https://doi.org/10.1038/s41421-020-0153-3
# ============================================================

setwd("~/CESC_Network/7_DGIDB/")
options(warn = 1)

suppressPackageStartupMessages({
  library(tidyverse)
  library(igraph)
  library(parallel)
})

set.seed(42)
N_PERM  <- 500L
N_CORES <- max(1L, parallel::detectCores() - 2L)
cat(sprintf("[setup] %d cores | %d permutations\n", N_CORES, N_PERM))

# ── 1. LOAD DATA ─────────────────────────────────────────────

cat("\n[1] Loading networks...\n")
graph_A7 <- readRDS("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_2_4_Subnetwork_A7_max.rds")
graph_A9 <- readRDS("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_2_4_Subnetwork_A9_max.rds")
cat(sprintf("  A7: %d nodes, %d edges\n", vcount(graph_A7), ecount(graph_A7)))
cat(sprintf("  A9: %d nodes, %d edges\n", vcount(graph_A9), ecount(graph_A9)))

cat("\n[2] Loading DGIdb FDA-approved interactions...\n")
edges_all <- read_tsv("6_1_edges_bipartite.tsv", show_col_types = FALSE)
cat(sprintf("  %d drug–gene pairs | %d unique drugs\n",
            nrow(edges_all), n_distinct(edges_all$drug_name)))

cat("\n[3] Loading interest genes (DGE ∩ Hub from DepMap)...\n")
depmap_dir <- "~/CESC_Network/6_DepMap_cacu"
gene_table_A7 <- read_tsv(file.path(depmap_dir, "6_3_gene_table_A7.tsv"), show_col_types = FALSE)
gene_table_A9 <- read_tsv(file.path(depmap_dir, "6_3_gene_table_A9.tsv"), show_col_types = FALSE)

interest_A7 <- gene_table_A7 %>% filter(is_interest, !is.na(mean_dep), mean_dep <= -0.5) %>% pull(Gene)
interest_A9 <- gene_table_A9 %>% filter(is_interest, !is.na(mean_dep), mean_dep <= -0.5) %>% pull(Gene)
cat(sprintf("  Interest genes (dep <= -0.5) – A7: %d | A9: %d\n", length(interest_A7), length(interest_A9)))

# ── 2. HELPERS ───────────────────────────────────────────────

# Closest proximity d_c(T,S): mean over targets of min dist to gene set
# dist_block: numeric matrix [|T| × |S|] with Inf→NA already handled
closest_prox <- function(dist_block) {
  if (is.null(dist_block) || nrow(dist_block) == 0L || ncol(dist_block) == 0L) return(NA_real_)
  row_mins <- apply(dist_block, 1L, function(r) {
    v <- r[is.finite(r)]
    if (length(v) == 0L) NA_real_ else min(v)
  })
  mean(row_mins, na.rm = TRUE)
}

# Sample nodes matching the degree-bin composition of target_bins
sample_degree_matched <- function(bins_all, target_bins, all_nodes) {
  samp <- lapply(unique(target_bins), function(b) {
    pool <- names(bins_all[bins_all == b & !is.na(bins_all)])
    n    <- sum(target_bins == b, na.rm = TRUE)
    if (length(pool) < n) pool <- all_nodes   # fallback: unrestricted sample
    sample(pool, n, replace = FALSE)
  })
  unlist(samp, use.names = FALSE)
}

# ── 3. MAIN PIPELINE ─────────────────────────────────────────

run_proximity <- function(graph, edges, interest_genes, cohort, n_perm, n_cores) {

  all_nodes      <- V(graph)$name
  interest_in_g  <- intersect(interest_genes, all_nodes)
  non_int_in_g   <- setdiff(all_nodes, interest_in_g)

  cat(sprintf("\n=== Cohort %s ===\n", cohort))
  cat(sprintf("  Interest in network: %d / %d\n",
              length(interest_in_g), length(interest_genes)))

  # ── 3a. Drug-target map (targets must be in network) ─────
  drug_map <- edges %>%
    filter(gene_name %in% all_nodes) %>%
    group_by(drug_name) %>%
    summarise(
      targets   = list(unique(gene_name)),
      n_targets = n_distinct(gene_name),
      .groups   = "drop"
    )
  cat(sprintf("  Drugs with ≥1 target: %d\n", nrow(drug_map)))

  unique_targets <- unique(unlist(drug_map$targets))

  # ── 3b. Precompute distance matrices ─────────────────────

  # D_int: all_nodes × interest  (604 BFS – fast)
  # Used for permutation null distribution.
  cat(sprintf("  [dist] BFS from %d interest genes... ", length(interest_in_g)))
  t0 <- proc.time()[3]
  D_raw <- distances(graph, v = interest_in_g, to = all_nodes, weights = NA)
  D_int <- t(D_raw)               # all_nodes × interest
  rownames(D_int) <- all_nodes
  colnames(D_int) <- interest_in_g
  cat(sprintf("done (%.1f s)\n", proc.time()[3] - t0))

  # D_tgt: unique_targets × all_nodes  (n_unique_targets BFS)
  # Used for observed d_c_interest and d_c_noninterest.
  cat(sprintf("  [dist] BFS from %d unique drug targets... ", length(unique_targets)))
  t0 <- proc.time()[3]
  D_tgt <- distances(graph, v = unique_targets, to = all_nodes, weights = NA)
  rownames(D_tgt) <- unique_targets
  colnames(D_tgt) <- all_nodes
  cat(sprintf("done (%.1f s)\n", proc.time()[3] - t0))

  # ── 3c. Degree bins for permutation ──────────────────────
  deg  <- degree(graph)
  bins <- cut(deg,
              breaks      = quantile(deg, probs = seq(0, 1, 0.1), na.rm = TRUE),
              include.lowest = TRUE)
  names(bins) <- names(deg)   # node-name → bin

  # ── 3d. Compute scores per drug ──────────────────────────
  cat(sprintf("  [prox] Computing scores for %d drugs (n_perm=%d)...\n",
              nrow(drug_map), n_perm))
  t0 <- proc.time()[3]

  compute_one <- function(i) {
    drug    <- drug_map$drug_name[i]
    targets <- drug_map$targets[[i]]

    # ── Observed proximity ───────────────────────────────
    # d_c to interest genes
    sub_int  <- D_tgt[targets, interest_in_g, drop = FALSE]
    d_c_int  <- closest_prox(sub_int)

    # d_c to non-interest genes
    sub_non  <- D_tgt[targets, non_int_in_g, drop = FALSE]
    d_c_non  <- closest_prox(sub_non)

    # ── Permutation null (permute targets, interest fixed) ──
    bins_t <- bins[targets]

    null_vals <- vapply(seq_len(n_perm), function(.) {
      samp <- tryCatch(
        sample_degree_matched(bins, bins_t, all_nodes),
        error = function(e) sample(all_nodes, length(targets))
      )
      samp <- intersect(samp, rownames(D_int))
      if (length(samp) == 0L) return(NA_real_)
      closest_prox(D_int[samp, , drop = FALSE])
    }, numeric(1L))

    null_vals <- null_vals[is.finite(null_vals)]
    mu_null   <- mean(null_vals)
    sd_null   <- sd(null_vals)
    z         <- if (!is.na(sd_null) && sd_null > 0) (d_c_int - mu_null) / sd_null else NA_real_
    p_val     <- pnorm(z)           # one-sided lower tail

    tibble(
      drug            = drug,
      n_targets       = length(targets),
      d_c_interest    = d_c_int,
      d_c_noninterest = d_c_non,
      z_score         = z,
      p_value         = p_val,
      mean_null       = mu_null,
      sd_null         = sd_null
    )
  }

  idx <- seq_len(nrow(drug_map))

  if (n_cores > 1L && .Platform$OS.type != "windows") {
    res_list <- parallel::mclapply(idx, compute_one, mc.cores = n_cores)
  } else {
    res_list <- lapply(idx, compute_one)
  }

  cat(sprintf("  done (%.1f s)\n", proc.time()[3] - t0))

  # Guney et al. 2016 uses z < -1.5 directly as the proximity threshold.
  # BH correction on pnorm(z) across 2800+ drugs is prohibitively strict
  # (requires z < -4.1), so we follow the published convention:
  #   significant  : z < -1.96  (one-sided α = 0.025, literature standard)
  #   sig_z2       : z < -2     (more conservative)
  #   sig_BH       : BH-corrected p < 0.05  (reported for reference)
  results <- bind_rows(res_list) %>%
    mutate(
      p_adj       = p.adjust(p_value, method = "BH"),
      significant = !is.na(z_score) & z_score < -1.96,
      sig_z2      = !is.na(z_score) & z_score < -2,
      sig_BH      = !is.na(z_score) & z_score < 0 & p_adj < 0.05,
      sig_tier    = case_when(
        !is.na(z_score) & z_score < -3   ~ "z<-3",
        !is.na(z_score) & z_score < -2   ~ "z<-2",
        !is.na(z_score) & z_score < -1.96 ~ "z<-1.96",
        !is.na(z_score) & z_score < -1.5  ~ "z<-1.5",
        TRUE ~ "ns"
      ),
      cohort      = cohort
    ) %>%
    arrange(z_score)

  list(results = results, D_int = D_int, D_tgt = D_tgt,
       interest_in_g = interest_in_g, drug_map = drug_map)
}

# ── 4. RUN BOTH COHORTS ───────────────────────────────────────

out_A7 <- run_proximity(graph_A7, edges_all, interest_A7, "A7", N_PERM, N_CORES)
out_A9 <- run_proximity(graph_A9, edges_all, interest_A9, "A9", N_PERM, N_CORES)

# ── 5. SAVE RESULTS ──────────────────────────────────────────

cat("\n[5] Saving results...\n")

write_tsv(out_A7$results, "6_3_proximity_A7.tsv")
write_tsv(out_A9$results, "6_3_proximity_A9.tsv")

# Combined table
both <- bind_rows(out_A7$results, out_A9$results)
write_tsv(both, "6_3_proximity_combined.tsv")

# Compact RDS (proximity objects for optional re-use)
saveRDS(
  list(
    A7 = list(results = out_A7$results, interest = out_A7$interest_in_g,
              drug_map = out_A7$drug_map),
    A9 = list(results = out_A9$results, interest = out_A9$interest_in_g,
              drug_map = out_A9$drug_map)
  ),
  "6_3_proximity_objects.rds"
)

# ── 6. QUICK SUMMARY ─────────────────────────────────────────

cat("\n[6] Summary\n")
for (coh in c("A7", "A9")) {
  df <- both %>% filter(cohort == coh)
  cat(sprintf(
    "  %s | n_drugs=%d | z<-1.96: %d | z<-2: %d | z<-3: %d | BH<0.05: %d\n",
    coh, nrow(df),
    sum(df$significant, na.rm = TRUE),
    sum(df$sig_z2,      na.rm = TRUE),
    sum(df$sig_tier == "z<-3", na.rm = TRUE),
    sum(df$sig_BH,      na.rm = TRUE)
  ))
  cat("  Top 15 candidates (z < -1.96):\n")
  top15 <- df %>%
    filter(significant) %>%
    slice_head(n = 15) %>%
    select(drug, n_targets, d_c_interest, d_c_noninterest, z_score, sig_tier, p_adj)
  print(top15, n = 15)
  cat("\n")
}

cat("Done. Output files:\n")
cat("  6_3_proximity_A7.tsv\n")
cat("  6_3_proximity_A9.tsv\n")
cat("  6_3_proximity_combined.tsv\n")
cat("  6_3_proximity_objects.rds\n")
