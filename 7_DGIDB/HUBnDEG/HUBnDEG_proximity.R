# ============================================================
# HUBnDEG_proximity.R
# Network-based drug proximity to HUBnDEG genes
# ============================================================
# Gene set  : intersect(DEGs, top_N_connected)
#             where N = number of DEGs per cohort
#             A7: top 7086 connected nodes ∩ 7086 DEGs
#             A9: top 7160 connected nodes ∩ 7160 DEGs
#
# Rationale : Uses the same gene-set size as the full DEG run
#             but restricts to genes that are also highly
#             connected in the network — a biologically
#             motivated hub-aware DEG intersection.
#
# Metric    : Closest proximity (Guney et al. 2016, Nat Commun)
#   d_c(T,S) = (1/|T|) * Σ_{t∈T}  min_{s∈S}  dist(t,s)
#
# Significance: degree-preserving permutation of drug targets
#   z = (d_c_obs - µ_null) / σ_null
#   p = Φ(z)  [one-sided: negative z = closer than expected]
#
# References:
#   Guney et al.  2016 – https://doi.org/10.1038/ncomms10331
#   Cheng et al.  2019 – https://doi.org/10.1038/s41467-019-10744-6
#   Zhou  et al.  2020 – https://doi.org/10.1038/s41421-020-0153-3
# ============================================================

setwd("~/CESC_Network/7_DGIDB/HUBnDEG/")
options(warn = 1)

suppressPackageStartupMessages({
  library(tidyverse)
  library(igraph)
  library(parallel)
})

set.seed(42)
N_PERM  <- 500L
N_CORES <- 9L
BATCH   <- N_CORES * 3L  # 27 drugs per batch — limits peak RAM per wave
cat(sprintf("[setup] %d cores | %d permutations | batch=%d\n", N_CORES, N_PERM, BATCH))

# ── 1. LOAD DATA ─────────────────────────────────────────────

cat("\n[1] Loading networks...\n")
graph_A7 <- readRDS("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_2_4_Subnetwork_A7_max.rds")
graph_A9 <- readRDS("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_2_4_Subnetwork_A9_max.rds")
cat(sprintf("  A7: %d nodes, %d edges\n", vcount(graph_A7), ecount(graph_A7)))
cat(sprintf("  A9: %d nodes, %d edges\n", vcount(graph_A9), ecount(graph_A9)))

cat("\n[2] Loading DGIdb FDA-approved interactions...\n")
edges_all <- read_tsv("../6_1_edges_bipartite.tsv", show_col_types = FALSE)
cat(sprintf("  %d drug–gene pairs | %d unique drugs\n",
            nrow(edges_all), n_distinct(edges_all$drug_name)))

cat("\n[3] Loading DEG gene sets and computing HUBnDEG intersection...\n")
depmap_dir <- "~/CESC_Network/6_DepMap_cacu"
gene_table_A7 <- read_tsv(file.path(depmap_dir, "6_3_gene_table_A7.tsv"), show_col_types = FALSE)
gene_table_A9 <- read_tsv(file.path(depmap_dir, "6_3_gene_table_A9.tsv"), show_col_types = FALSE)

deg_A7 <- gene_table_A7 %>% filter(is_DEG) %>% pull(Gene)
deg_A9 <- gene_table_A9 %>% filter(is_DEG) %>% pull(Gene)

# Top-N most connected nodes (N = number of DEGs per cohort)
top_connected <- function(graph, n) {
  deg_order <- sort(degree(graph), decreasing = TRUE)
  names(deg_order)[seq_len(n)]
}

top_A7 <- top_connected(graph_A7, length(deg_A7))
top_A9 <- top_connected(graph_A9, length(deg_A9))

hubndeg_A7 <- intersect(deg_A7, top_A7)
hubndeg_A9 <- intersect(deg_A9, top_A9)

cat(sprintf("  DEGs           – A7: %d | A9: %d\n", length(deg_A7),      length(deg_A9)))
cat(sprintf("  Top-N connected– A7: %d | A9: %d\n", length(top_A7),      length(top_A9)))
cat(sprintf("  HUBnDEG        – A7: %d | A9: %d\n", length(hubndeg_A7),  length(hubndeg_A9)))

# ── 2. HELPERS ───────────────────────────────────────────────

closest_prox <- function(dist_block) {
  if (is.null(dist_block) || nrow(dist_block) == 0L || ncol(dist_block) == 0L) return(NA_real_)
  row_mins <- apply(dist_block, 1L, function(r) {
    v <- r[is.finite(r)]
    if (length(v) == 0L) NA_real_ else min(v)
  })
  mean(row_mins, na.rm = TRUE)
}

sample_degree_matched <- function(bins_all, target_bins, all_nodes) {
  samp <- lapply(unique(target_bins), function(b) {
    pool <- names(bins_all[bins_all == b & !is.na(bins_all)])
    n    <- sum(target_bins == b, na.rm = TRUE)
    if (length(pool) < n) pool <- all_nodes
    sample(pool, n, replace = FALSE)
  })
  unlist(samp, use.names = FALSE)
}

# ── 3. MAIN PIPELINE ─────────────────────────────────────────

run_proximity_hubndeg <- function(graph, edges, hubndeg_genes, cohort, n_perm, n_cores, batch) {

  all_nodes      <- V(graph)$name
  hubndeg_in_g   <- intersect(hubndeg_genes, all_nodes)
  non_hubndeg_in_g <- setdiff(all_nodes, hubndeg_in_g)

  cat(sprintf("\n=== Cohort %s ===\n", cohort))
  cat(sprintf("  HUBnDEG in network : %d / %d\n", length(hubndeg_in_g), length(hubndeg_genes)))
  cat(sprintf("  Non-HUBnDEG nodes  : %d\n",       length(non_hubndeg_in_g)))

  # ── 3a. Drug-target map ──────────────────────────────────
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

  # D_hnd: all_nodes × hubndeg_genes  (used for permutation null)
  cat(sprintf("  [dist] BFS from %d HUBnDEG genes... ", length(hubndeg_in_g)))
  t0 <- proc.time()[3]
  D_raw <- distances(graph, v = hubndeg_in_g, to = all_nodes, weights = NA)
  D_hnd <- t(D_raw)               # all_nodes × hubndeg
  rownames(D_hnd) <- all_nodes
  colnames(D_hnd) <- hubndeg_in_g
  cat(sprintf("done (%.1f s)\n", proc.time()[3] - t0))

  # D_tgt: unique_targets × all_nodes
  cat(sprintf("  [dist] BFS from %d unique drug targets... ", length(unique_targets)))
  t0 <- proc.time()[3]
  D_tgt <- distances(graph, v = unique_targets, to = all_nodes, weights = NA)
  rownames(D_tgt) <- unique_targets
  colnames(D_tgt) <- all_nodes
  cat(sprintf("done (%.1f s)\n", proc.time()[3] - t0))

  # ── 3c. Degree bins for permutation ──────────────────────
  deg  <- degree(graph)
  bins <- cut(deg,
              breaks         = quantile(deg, probs = seq(0, 1, 0.1), na.rm = TRUE),
              include.lowest = TRUE)
  names(bins) <- names(deg)

  # ── 3d. Compute scores per drug ──────────────────────────
  cat(sprintf("  [prox] %d drugs | %d perms | batch=%d\n",
              nrow(drug_map), n_perm, batch))
  t0 <- proc.time()[3]

  compute_one <- function(i) {
    drug    <- drug_map$drug_name[i]
    targets <- drug_map$targets[[i]]

    # Observed d_c to HUBnDEG genes
    sub_hnd  <- D_tgt[targets, hubndeg_in_g,     drop = FALSE]
    d_c_hnd  <- closest_prox(sub_hnd)

    # Observed d_c to non-HUBnDEG genes
    sub_non  <- D_tgt[targets, non_hubndeg_in_g, drop = FALSE]
    d_c_non  <- closest_prox(sub_non)

    # Permutation null (permute targets, HUBnDEG set fixed)
    bins_t <- bins[targets]

    null_vals <- vapply(seq_len(n_perm), function(.) {
      samp <- tryCatch(
        sample_degree_matched(bins, bins_t, all_nodes),
        error = function(e) sample(all_nodes, length(targets))
      )
      samp <- intersect(samp, rownames(D_hnd))
      if (length(samp) == 0L) return(NA_real_)
      closest_prox(D_hnd[samp, , drop = FALSE])
    }, numeric(1L))

    null_vals <- null_vals[is.finite(null_vals)]
    mu_null   <- mean(null_vals)
    sd_null   <- sd(null_vals)
    z         <- if (!is.na(sd_null) && sd_null > 0) (d_c_hnd - mu_null) / sd_null else NA_real_
    p_val     <- pnorm(z)

    tibble(
      drug             = drug,
      n_targets        = length(targets),
      d_c_HUBnDEG      = d_c_hnd,
      d_c_nonHUBnDEG   = d_c_non,
      z_score          = z,
      p_value          = p_val,
      mean_null        = mu_null,
      sd_null          = sd_null
    )
  }

  idx      <- seq_len(nrow(drug_map))
  res_list <- vector("list", length(idx))

  if (n_cores > 1L && .Platform$OS.type != "windows") {
    n_batches <- ceiling(length(idx) / batch)
    for (b in seq_len(n_batches)) {
      i1        <- (b - 1L) * batch + 1L
      i2        <- min(b * batch, length(idx))
      batch_idx <- idx[i1:i2]
      tmp       <- parallel::mclapply(batch_idx, compute_one,
                                      mc.cores = n_cores, mc.cleanup = TRUE)
      res_list[batch_idx] <- tmp
      rm(tmp); gc(verbose = FALSE)
      cat(sprintf("  batch %d/%d | drugs %d-%d / %d (%.0f%%) | %.1f s\n",
                  b, n_batches, i1, i2, length(idx),
                  i2 / length(idx) * 100, proc.time()[3] - t0))
    }
  } else {
    res_list <- lapply(idx, compute_one)
  }

  cat(sprintf("  done (%.1f s)\n", proc.time()[3] - t0))

  results <- bind_rows(res_list) %>%
    mutate(
      p_adj       = p.adjust(p_value, method = "BH"),
      significant = !is.na(z_score) & z_score < -1.96,
      sig_z2      = !is.na(z_score) & z_score < -2,
      sig_BH      = !is.na(z_score) & z_score < 0 & p_adj < 0.05,
      sig_tier    = case_when(
        !is.na(z_score) & z_score < -3    ~ "z<-3",
        !is.na(z_score) & z_score < -2    ~ "z<-2",
        !is.na(z_score) & z_score < -1.96 ~ "z<-1.96",
        !is.na(z_score) & z_score < -1.5  ~ "z<-1.5",
        TRUE ~ "ns"
      ),
      cohort = cohort
    ) %>%
    arrange(z_score)

  list(results = results, D_hnd = D_hnd, D_tgt = D_tgt,
       hubndeg_in_g = hubndeg_in_g, drug_map = drug_map)
}

# ── 4. RUN BOTH COHORTS ───────────────────────────────────────

out_A7 <- run_proximity_hubndeg(graph_A7, edges_all, hubndeg_A7, "A7", N_PERM, N_CORES, BATCH)
out_A9 <- run_proximity_hubndeg(graph_A9, edges_all, hubndeg_A9, "A9", N_PERM, N_CORES, BATCH)

# ── 5. SAVE RESULTS ──────────────────────────────────────────

cat("\n[5] Saving results...\n")

write_tsv(out_A7$results, "proximity_HUBnDEG_A7.tsv")
write_tsv(out_A9$results, "proximity_HUBnDEG_A9.tsv")

both <- bind_rows(out_A7$results, out_A9$results)
write_tsv(both, "proximity_HUBnDEG_combined.tsv")

saveRDS(
  list(
    A7 = list(results = out_A7$results, hubndeg = out_A7$hubndeg_in_g,
              drug_map = out_A7$drug_map),
    A9 = list(results = out_A9$results, hubndeg = out_A9$hubndeg_in_g,
              drug_map = out_A9$drug_map)
  ),
  "proximity_HUBnDEG_objects.rds"
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
    select(drug, n_targets, d_c_HUBnDEG, d_c_nonHUBnDEG, z_score, sig_tier, p_adj)
  print(top15, n = 15)
  cat("\n")
}

cat("Done. Output files:\n")
cat("  proximity_HUBnDEG_A7.tsv\n")
cat("  proximity_HUBnDEG_A9.tsv\n")
cat("  proximity_HUBnDEG_combined.tsv\n")
cat("  proximity_HUBnDEG_objects.rds\n")
