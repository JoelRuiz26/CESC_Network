### ============================================================
### 04_gsea_hallmark_generic.R
###
### Generic rank-based GSEA (Gene Set Enrichment Analysis) against the
### MSigDB Hallmark gene sets. Takes ANY differential-expression table as
### input -- it only needs a gene-symbol column and a ranking-statistic
### column (e.g. limma's t, or logFC) -- so it is not tied to one specific
### comparison; call run_gsea_hallmark() on whatever DEG table you have.
###
### Gene sets: MSigDB Hallmark (collection "H", 50 gene sets). This is the
### SAME collection already used elsewhere in this repository
### (02_hallmark_enrichment.R here, and
### 09_functional_similarity_degs/05_gsea_hallmark.R) -- the curated,
### non-redundant "Hallmark Molecular Signatures" built by Liberzon et al.
### 2015 (Cell Systems, doi:10.1016/j.cels.2015.12.004) by computationally
### merging thousands of overlapping "founder" gene sets from the other
### MSigDB collections (C1-C8) into 50 representative sets, one per
### well-defined biological process/state. This is the collection the
### manuscript was finally written around, after an earlier side-by-side
### comparison against a larger, less curated candidate set (~600 "core"
### cancer genes, and a ~3,000-gene list compiled directly from a
### Hallmarks-of-Cancer review): MSigDB Hallmark was preferred because it
### is versioned, citable, already split into 50 coherent, non-redundant
### sets (rather than one flat gene list), and directly comparable across
### every other enrichment already run in this project.
###
### Gene sets are never downloaded by hand: the msigdbr package fetches
### them programmatically and ships them version-pinned with the package
### itself (msigdbr v26.1.0 here -> MSigDB db_version "2026.1.Hs"), so the
### exact gene-set definitions are fully reproducible from the R package
### version alone (Dolgalev I, msigdbr: MSigDB Gene Sets for Multiple
### Organisms in a Tidy Data Format, CRAN).
###
### Usage: source this file, then e.g.
###   deg <- readr::read_tsv("~/path/to/your_deg_table.tsv")
###   run_gsea_hallmark(deg, gene_col = "gene", stat_col = "t", label = "my_comparison")
### -> results/04_gsea_hallmark_<label>.tsv
### ============================================================
suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(readr); library(clusterProfiler); library(msigdbr)
})

set.seed(123)
out_dir <- "~/CESC_Network/06_modules/results"

## -- Hallmark gene sets (H collection, 50 sets), fetched via msigdbr and
## reshaped into the two-column TERM2GENE format clusterProfiler expects
t2g <- msigdbr(species = "Homo sapiens", collection = "H") |>
  dplyr::select(gs_name, gene_symbol) |>
  dplyr::distinct()
cat(sprintf("[Hallmark] %d gene sets (%d gene-set memberships)\n",
            dplyr::n_distinct(t2g$gs_name), nrow(t2g)))

## -- Generic GSEA runner. `deg` needs a gene-symbol column (`gene_col`)
## and a numeric ranking statistic (`stat_col`); ties on gene symbol keep
## the first row. Returns (and saves) the full GSEA result table.
run_gsea_hallmark <- function(deg, gene_col = "gene", stat_col = "t", label = "run") {
  rnk <- deg |>
    dplyr::filter(!is.na(.data[[stat_col]])) |>
    dplyr::distinct(.data[[gene_col]], .keep_all = TRUE) |>
    dplyr::arrange(dplyr::desc(.data[[stat_col]])) |>
    (\(d) setNames(d[[stat_col]], d[[gene_col]]))()

  set.seed(123)
  gs <- GSEA(geneList = rnk, TERM2GENE = t2g, pvalueCutoff = 1, eps = 0, seed = TRUE,
             pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500)

  res <- as.data.frame(gs) |>
    tibble::as_tibble() |>
    dplyr::mutate(HM = gsub("_", " ", sub("^HALLMARK_", "", ID)), run = label) |>
    dplyr::arrange(p.adjust)

  readr::write_tsv(res, file.path(out_dir, sprintf("04_gsea_hallmark_%s.tsv", label)))
  cat(sprintf("[%s] %d/%d Hallmark terms significant (adj.P<0.05)\n",
              label, sum(res$p.adjust < 0.05), nrow(res)))
  res
}

## ---- example (edit and uncomment to run) ----
# deg <- readr::read_tsv("~/path/to/your_deg_table.tsv")
# run_gsea_hallmark(deg, gene_col = "gene", stat_col = "t", label = "my_comparison")
