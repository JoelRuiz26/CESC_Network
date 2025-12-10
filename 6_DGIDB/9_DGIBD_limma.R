############################################################
## SCRIPT: 8_5_DGIDB_ClosestProximity_Zscore_ModuleHeatmap_OPT.R
##
## Goal:
##   Apply "closest proximity" (dc) and degree-controlled
##   z-score (zc) between drug targets and key DGE genes,
##   then visualize module-level proximities in A7 and A9
##   Louvain networks.
##
## Key changes vs versión anterior:
##   - SIN parallel::future: compute_drug_proximity usa lapply.
##   - dist_to_keys se mantiene local (no se exporta a workers).
##   - Se guardan matrices mat_A7 y mat_A9 (TSV + RDS).
############################################################

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(igraph)
  library(dplyr)
  library(readr)
  library(vroom)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(ComplexHeatmap)
  library(circlize)
  library(future.apply)  # ya no se usa future_lapply, pero no estorba
  library(grid)
  library(tibble)
})

set.seed(123)

# ----------------------------
# Paths
# ----------------------------
outdir <- "~/CESC_Network/6_DGIDB/"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 1) Load Louvain community graphs
# ----------------------------
graph_A7 <- readRDS("~/CESC_Network/5_Modularity/5_1_Modules/5_2_1_Louvain_A7_graph_max.rds")
graph_A9 <- readRDS("~/CESC_Network/5_Modularity/5_1_Modules/5_2_1_Louvain_A9_graph_max.rds")

A7_genes <- V(graph_A7)$name
A9_genes <- V(graph_A9)$name
genes_universe <- unique(c(A7_genes, A9_genes))

# ----------------------------
# 2) Load limma DGE tables and define key genes
# ----------------------------
DGE_A7 <- read_tsv(
  "~/CESC_Network/4_Network_analisis/4_3_limma_A7_vs_SNT_allGenes_annot_PR_DEG.tsv",
  col_types = cols()
)

DGE_A9 <- read_tsv(
  "~/CESC_Network/4_Network_analisis/4_3_limma_A9_vs_SNT_allGenes_annot_PR_DEG.tsv",
  col_types = cols()
)

padj_cutoff <- 0.01

key_genes_A7 <- DGE_A7 %>%
  filter(!is.na(adj.P.Val), adj.P.Val < padj_cutoff) %>%
  filter(abs(logFC) > 1) %>%
  pull(Gene) %>%
  unique() %>%
  intersect(V(graph_A7)$name)

key_genes_A9 <- DGE_A9 %>%
  filter(!is.na(adj.P.Val), adj.P.Val < padj_cutoff) %>%
  filter(abs(logFC) > 1) %>%
  pull(Gene) %>%
  unique() %>%
  intersect(V(graph_A9)$name)

# ----------------------------
# 3) Download DGIdb interactions (latest)
# ----------------------------
url_interactions <- "https://dgidb.org/data/latest/interactions.tsv"
dgidb_interactions <- read_tsv(url_interactions, show_col_types = FALSE)

dgidb_fda <- dgidb_interactions %>%
  filter(approved == "TRUE") %>%
  mutate(
    interaction_score = ifelse(interaction_score == "NULL", "0", interaction_score),
    interaction_score = as.numeric(interaction_score),
    drug_name = as.character(drug_name),
    gene_name = as.character(gene_name)
  )

edges_all <- dgidb_fda %>%
  filter(gene_name %in% genes_universe) %>%
  distinct(drug_name, gene_name)

write_tsv(edges_all, file.path(outdir, "8_5_edges_all_FDA_in_universe.tsv"))

# ============================================================
# 4) Utility functions
# ============================================================

closest_dc_graph <- function(graph, targets, disease_genes) {
  targets <- intersect(targets, V(graph)$name)
  disease_genes <- intersect(disease_genes, V(graph)$name)
  if (length(targets) == 0 || length(disease_genes) == 0) return(NA_real_)
  dmat <- distances(graph, v = targets, to = disease_genes)
  mins <- apply(dmat, 1, min, na.rm = TRUE)
  mean(mins, na.rm = TRUE)
}

make_degree_bins <- function(graph, n_bins = 10) {
  deg <- igraph::degree(graph)
  deg <- as.numeric(deg)
  names(deg) <- V(graph)$name
  
  qs <- stats::quantile(
    deg,
    probs = seq(0, 1, length.out = n_bins + 1),
    na.rm = TRUE
  )
  qs <- unique(qs)
  if (length(qs) < 2) {
    stop("Not enough distinct degree values to build bins.")
  }
  bins <- cut(deg, breaks = qs, include.lowest = TRUE, right = TRUE)
  names(bins) <- V(graph)$name
  bins
}

sample_degree_matched <- function(nodes, bins) {
  nodes <- intersect(nodes, names(bins))
  if (length(nodes) == 0) return(character(0))
  
  bins_nodes <- bins[nodes]
  all_names  <- names(bins)
  
  sampled <- unlist(lapply(levels(bins), function(b) {
    n_b <- sum(bins_nodes == b, na.rm = TRUE)
    if (n_b == 0) return(character(0))
    pool <- all_names[bins == b]
    if (length(pool) == 0) {
      return(sample(all_names, n_b))
    }
    if (length(pool) < n_b) {
      return(sample(pool, n_b, replace = TRUE))
    }
    sample(pool, n_b, replace = FALSE)
  }))
  
  sampled
}

# ============================================================
# 5) Cached dc and zc using precomputed distances (genes -> key genes)
# ============================================================

precompute_dist_to_keys <- function(graph, key_genes) {
  key_genes <- intersect(key_genes, V(graph)$name)
  if (length(key_genes) == 0) stop("No key genes in graph.")
  dmat <- distances(graph, v = V(graph), to = key_genes)
  rownames(dmat) <- V(graph)$name
  colnames(dmat) <- key_genes
  dmat
}

closest_dc_cached <- function(targets, dist_to_keys) {
  targets <- intersect(targets, rownames(dist_to_keys))
  if (length(targets) == 0) return(NA_real_)
  sub <- dist_to_keys[targets, , drop = FALSE]
  mins <- apply(sub, 1, min, na.rm = TRUE)
  mean(mins, na.rm = TRUE)
}

closest_zc_cached <- function(targets, dist_to_keys, bins, nperm = 100) {
  targets <- intersect(targets, rownames(dist_to_keys))
  if (length(targets) == 0) return(NA_real_)
  
  obs <- closest_dc_cached(targets, dist_to_keys)
  
  rand_vals <- replicate(nperm, {
    rt <- sample_degree_matched(targets, bins)
    rt <- intersect(rt, rownames(dist_to_keys))
    if (length(rt) == 0) return(NA_real_)
    sub <- dist_to_keys[rt, , drop = FALSE]
    mins <- apply(sub, 1, min, na.rm = TRUE)
    mean(mins, na.rm = TRUE)
  })
  
  mu  <- mean(rand_vals, na.rm = TRUE)
  sdv <- sd(rand_vals, na.rm = TRUE)
  if (is.na(sdv) || sdv == 0) return(NA_real_)
  
  (obs - mu) / sdv
}

# ============================================================
# 6) Compute per-drug dc and zc vs ALL key genes (SIN parallel)
# ============================================================

compute_drug_proximity <- function(graph, edges, key_genes, label,
                                   nperm = 100, n_bins = 10) {
  
  bins <- make_degree_bins(graph, n_bins = n_bins)
  dist_to_keys <- precompute_dist_to_keys(graph, key_genes)
  
  drug_targets <- edges %>%
    dplyr::filter(gene_name %in% V(graph)$name) %>%
    dplyr::group_by(drug_name) %>%
    dplyr::summarise(
      targets = list(unique(gene_name)),
      .groups = "drop"
    )
  
  # ---- AQUÍ EL CAMBIO IMPORTANTE: lapply en lugar de future_lapply ----
  res_list <- lapply(seq_len(nrow(drug_targets)), function(i) {
    drug    <- drug_targets$drug_name[i]
    targets <- drug_targets$targets[[i]]
    
    dc_val <- closest_dc_cached(targets, dist_to_keys)
    zc_val <- closest_zc_cached(targets, dist_to_keys, bins = bins, nperm = nperm)
    
    tibble(
      Cohort    = label,
      drug      = drug,
      n_targets = length(intersect(targets, rownames(dist_to_keys))),
      dc        = dc_val,
      zc        = zc_val
    )
  })
  
  res <- bind_rows(res_list) %>%
    arrange(zc, dc)
  
  res
}

# ============================================================
# 7) Run for A7 and A9
# ============================================================

message("Computing drug proximity for A7...")
results_A7 <- compute_drug_proximity(
  graph = graph_A7,
  edges = edges_all,
  key_genes = key_genes_A7,
  label = "A7",
  nperm = 100
)

message("Computing drug proximity for A9...")
results_A9 <- compute_drug_proximity(
  graph = graph_A9,
  edges = edges_all,
  key_genes = key_genes_A9,
  label = "A9",
  nperm = 100
)

write_tsv(results_A7, file.path(outdir, "8_5_ClosestProximity_A7_FULL.tsv"))
write_tsv(results_A9, file.path(outdir, "8_5_ClosestProximity_A9_FULL.tsv"))

# ============================================================
# 8) Select Top 20 drugs by best zc
# ============================================================

top_n <- 20

top20_A7 <- results_A7 %>%
  filter(!is.na(zc)) %>%
  slice_head(n = top_n) %>%
  pull(drug)

top20_A9 <- results_A9 %>%
  filter(!is.na(zc)) %>%
  slice_head(n = top_n) %>%
  pull(drug)

write_tsv(tibble(drug = top20_A7),
          file.path(outdir, "8_5_A7_Top20_by_zc.tsv"))
write_tsv(tibble(drug = top20_A9),
          file.path(outdir, "8_5_A9_Top20_by_zc.tsv"))

# ============================================================
# 9) Build module × drug dc matrix
# ============================================================

module_dc_matrix <- function(graph, edges, key_genes, drugs,
                             min_key_genes_per_module = 5) {
  
  comm_name <- V(graph)$Community_name
  if (is.null(comm_name)) {
    comm_name <- as.character(V(graph)$Community)
  }
  gene2comm <- setNames(comm_name, V(graph)$name)
  
  key_genes <- intersect(key_genes, V(graph)$name)
  key_by_comm <- split(key_genes, gene2comm[key_genes])
  
  key_by_comm <- key_by_comm[sapply(key_by_comm, length) >= min_key_genes_per_module]
  comms <- names(key_by_comm)
  if (length(comms) == 0) {
    warning("No modules with enough key genes.")
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }
  
  edges_g <- edges %>%
    filter(gene_name %in% V(graph)$name, drug_name %in% drugs)
  
  drug_targets <- edges_g %>%
    group_by(drug_name) %>%
    summarise(targets = list(unique(gene_name)), .groups = "drop")
  
  M <- sapply(drug_targets$drug_name, function(drug) {
    targets <- drug_targets$targets[[which(drug_targets$drug_name == drug)]]
    sapply(comms, function(cn) {
      S_mod <- key_by_comm[[cn]]
      closest_dc_graph(graph, targets, S_mod)
    })
  })
  
  M <- as.matrix(M)
  rownames(M) <- comms
  colnames(M) <- drug_targets$drug_name
  M
}

mat_A7 <- module_dc_matrix(graph_A7, edges_all, key_genes_A7, top20_A7)
mat_A9 <- module_dc_matrix(graph_A9, edges_all, key_genes_A9, top20_A9)

# ---- Guardar matrices para reusar plots sin recalcular proximidad ----
saveRDS(list(
  mat_A7 = mat_A7,
  mat_A9 = mat_A9,
  results_A7 = results_A7,
  results_A9 = results_A9
), file.path(outdir, "8_5_ModuleDrug_Proximity_Objects.rds"))

mat_A7_df <- as.data.frame(mat_A7) %>%
  rownames_to_column("Module")
mat_A9_df <- as.data.frame(mat_A9) %>%
  rownames_to_column("Module")

write_tsv(mat_A7_df, file.path(outdir, "8_5_A7_ModuleDrug_dc_matrix.tsv"))
write_tsv(mat_A9_df, file.path(outdir, "8_5_A9_ModuleDrug_dc_matrix.tsv"))

# ============================================================
# 10) Heatmaps (bicolor + monocromático estilo ejemplo)
# ============================================================

plot_and_save_heatmap <- function(mat, results_df, cohort_label, outfile_prefix) {
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    warning("Empty matrix for ", cohort_label)
    return(invisible(NULL))
  }
  
  z_map <- setNames(results_df$zc, results_df$drug)
  col_z <- z_map[colnames(mat)]
  
  ha <- HeatmapAnnotation(
    zc = col_z,
    annotation_name_side = "left"
  )
  
  rng <- range(mat, na.rm = TRUE)
  col_fun <- colorRamp2(
    c(rng[1], mean(rng), rng[2]),
    c("#08306B", "#F7F7F7", "#B30000")
  )
  
  ht <- Heatmap(
    mat,
    name = "dc (closest)",
    col = col_fun,
    top_annotation = ha,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 9),
    column_title = paste0("HPV ", cohort_label,
                          " | Module-level drug proximity (Top 20)"),
    column_title_gp = gpar(fontsize = 14, fontface = "bold")
  )
  
  pdf(file.path(outdir, paste0(outfile_prefix, ".pdf")), width = 12, height = 8)
  draw(ht)
  dev.off()
  
  png(file.path(outdir, paste0(outfile_prefix, ".png")),
      width = 1400, height = 900, res = 160)
  draw(ht)
  dev.off()
  
  invisible(ht)
}

annotate_heatmap <- function(mat, results_df, cohort_label, outfile_prefix) {
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    warning("Empty matrix for ", cohort_label)
    return(invisible(NULL))
  }
  
  z_map <- setNames(results_df$zc, results_df$drug)
  
  new_colnames <- ifelse(
    !is.na(z_map[colnames(mat)]) & z_map[colnames(mat)] < -1.5,
    paste0(colnames(mat), " *"),
    colnames(mat)
  )
  colnames(mat) <- new_colnames
  
  proximity_palette <- colorRamp2(
    c(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE)),
    c("#08306B", "#FFFAEC")
  )
  
  ht <- Heatmap(
    mat,
    name = "Mean Proximity",
    col = proximity_palette,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 14, fontface = "bold"),
    column_names_gp = gpar(fontsize = 12, fontface = "bold"),
    heatmap_legend_param = list(
      title = "Mean Proximity Distance",
      legend_height = unit(4, "cm")
    ),
    column_title = paste0(
      "HPV Clade ", cohort_label,
      ": Drug–Gene Community Proximity\n"
    ),
    column_title_gp = gpar(fontsize = 16, fontface = "bold")
  )
  
  pdf(file.path(outdir, paste0(outfile_prefix, "_mono.pdf")), width = 12, height = 8)
  draw(ht)
  dev.off()
  
  png(file.path(outdir, paste0(outfile_prefix, "_mono.png")),
      width = 1400, height = 900, res = 160)
  draw(ht)
  dev.off()
  
  invisible(ht)
}

ht_A7 <- plot_and_save_heatmap(
  mat_A7,
  results_A7,
  "A7",
  "8_5_A7_ClosestProximity_ModuleHeatmap_Top20"
)

ht_A9 <- plot_and_save_heatmap(
  mat_A9,
  results_A9,
  "A9",
  "8_5_A9_ClosestProximity_ModuleHeatmap_Top20"
)

ht_A7_mono <- annotate_heatmap(
  mat_A7,
  results_A7 %>% filter(drug %in% colnames(mat_A7)),
  "A7",
  "8_5_A7_ClosestProximity_ModuleHeatmap_Top20"
)

ht_A9_mono <- annotate_heatmap(
  mat_A9,
  results_A9 %>% filter(drug %in% colnames(mat_A9)),
  "A9",
  "8_5_A9_ClosestProximity_ModuleHeatmap_Top20"
)

message("Done!")
message("Saved outputs in: ", outdir)
