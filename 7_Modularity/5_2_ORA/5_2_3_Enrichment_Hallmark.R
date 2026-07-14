### Over-Representation Analysis with MSigDB Hallmark (H, 50 gene sets)
### Louvain modules, HPV-A7 and HPV-A9
### Outputs: per-network Hallmark ORA TSVs + per-module term-count summary
start_time <- Sys.time()

suppressPackageStartupMessages({
  library(igraph)
  library(dplyr)
  library(clusterProfiler)
  library(msigdbr)
  library(readr)
})

set.seed(123)
P_CUT   <- 0.05   # BH-adjusted p-value cutoff
Q_CUT   <- 0.20
out_dir <- "~/CESC_Network/7_Modularity/5_2_ORA/"

### Load module graphs (communities named M1, M2, ... by size) ###
graph_A7 <- readRDS("~/CESC_Network/7_Modularity/5_1_Modules/5_1_1_Louvain_A7_graph_max.rds")
graph_A9 <- readRDS("~/CESC_Network/7_Modularity/5_1_Modules/5_1_1_Louvain_A9_graph_max.rds")

### Universe (same background as the GO ORA) ###
universo <- readr::read_rds("~/CESC_Network/2_Prepro_TCGA_GTEx/2_1_Universe_annotation_table.rds") |>
  dplyr::pull(hgnc_symbol) |> unique()
cat(sprintf("[1] Universe: %d genes\n", length(universo)))

### Hallmark gene sets -> TERM2GENE (50 sets) ###
t2g <- msigdbr(species = "Homo sapiens", collection = "H") |>
  dplyr::select(gs_name, gene_symbol) |>
  dplyr::distinct()
cat(sprintf("[2] Hallmark gene sets: %d (%d gene-set memberships)\n",
            dplyr::n_distinct(t2g$gs_name), nrow(t2g)))

### ORA (enricher) for every module of one network ###
run_hallmark_ora <- function(graph, label) {
  comms <- sort(unique(V(graph)$Community))
  cat(sprintf("[ORA] %s — %d modules\n", label, length(comms)))
  res <- lapply(comms, function(cm) {
    genes <- V(graph)[V(graph)$Community == cm]$name
    cname <- unique(V(graph)[V(graph)$Community == cm]$Community_name)
    csize <- length(genes)
    if (length(genes) < 5) return(NULL)
    eg <- tryCatch(
      enricher(gene = genes, universe = universo, TERM2GENE = t2g,
               pvalueCutoff = P_CUT, qvalueCutoff = Q_CUT, pAdjustMethod = "BH"),
      error = function(e) NULL)
    if (is.null(eg) || nrow(as.data.frame(eg)) == 0) return(NULL)
    as.data.frame(eg) |>
      arrange(p.adjust, desc(Count)) |>
      mutate(Community = cm, Community_Name = cname, Community_size = csize, Network = label)
  })
  bind_rows(res)
}

ora7 <- run_hallmark_ora(graph_A7, "A7")
ora9 <- run_hallmark_ora(graph_A9, "A9")

write.table(ora7, file.path(out_dir, "5_2_3_Hallmark7_louv_max.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(ora9, file.path(out_dir, "5_2_3_Hallmark9_louv_max.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

### Per-module term-count summary (how many Hallmark terms per module) ###
term_counts <- bind_rows(ora7, ora9) |>
  count(Network, Community_Name, Community_size, name = "n_hallmark_terms") |>
  arrange(Network, desc(n_hallmark_terms))
write.table(term_counts, file.path(out_dir, "5_2_3_Hallmark_module_term_counts.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

### Report ###
report <- function(ora, total_mods, label) {
  enriched <- ora |> distinct(Community_Name) |> nrow()
  per_mod  <- ora |> count(Community_Name, Community_size, name = "n_terms") |>
    arrange(desc(n_terms), desc(Community_size))
  cat(sprintf("\n══ %s ══\n", label))
  cat(sprintf("Modules with >= 1 enriched Hallmark: %d / %d\n", enriched, total_mods))
  cat("Modules ranked by number of enriched Hallmark terms:\n")
  print(as.data.frame(per_mod), row.names = FALSE)
}
report(ora7, length(unique(V(graph_A7)$Community)), "HPV-A7")
report(ora9, length(unique(V(graph_A9)$Community)), "HPV-A9")

cat(sprintf("\n[OK] Saved Hallmark ORA TSVs + term-count summary. Time: %.1f min\n",
            as.numeric(difftime(Sys.time(), start_time, units = "mins"))))
