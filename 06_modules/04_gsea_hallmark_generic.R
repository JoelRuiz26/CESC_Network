### ============================================================
### 04_gsea_hallmark_generic.R
###
### Generic rank-based GSEA (Gene Set Enrichment Analysis) against the
### literature-curated "Cancer Hallmark" gene sets -- NOT MSigDB's Hallmark
### collection. Takes ANY differential-expression table as input -- it only
### needs a gene-symbol column and a ranking-statistic column (e.g. limma's
### t, or logFC) -- so it is not tied to one specific comparison.
###
### Gene sets: Menyhart O, Jayasekara Kothalawala W, Gyorffy B. "A gene set
### enrichment analysis for cancer hallmarks." J Pharm Anal. 2025;15(5):
### 101065. doi:10.1016/j.jpha.2024.101065 (companion web tool:
### cancerhallmarks.com). The authors reconciled 7 independent hallmark-gene
### resources against Hanahan & Weinberg's 10 hallmarks of cancer and
### release two gene sets per hallmark:
###   "core"       genes backed by >=2 of the 7 source resources
###                (n = 1,574 genes total) -- the higher-confidence,
###                cross-database CONSENSUS set, used here by default.
###   "integrated" genes backed by >=1 resource, i.e. the union of all 7
###                (n = 6,763 genes total) -- a looser, larger set.
### Both are plain-text GMT-style files (one line per hallmark: hallmark
### name, then its member gene symbols, tab-separated) served directly by
### cancerhallmarks.com; fetched here on the fly with readLines() -- no
### manual download and no account/API key needed.
###
### Usage: source this file, then e.g.
###   deg <- readr::read_tsv("~/path/to/your_deg_table.tsv")
###   run_gsea_hallmarks(deg, gene_col = "gene", stat_col = "t", label = "my_comparison")
### -> results/04_gsea_hallmarks_<label>.tsv
### ============================================================
suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(readr); library(clusterProfiler)
})

set.seed(123)
out_dir <- "~/CESC_Network/06_modules/results"

HALLMARK_URLS <- c(
  core       = "https://cancerhallmarks.com/download_file/Menyhart_JPA_CancerHallmarks_core.txt",
  integrated = "https://cancerhallmarks.com/download_file/Menyhart_JPA_CancerHallmarks_integrated.txt"
)

## -- Fetch one of the two hallmark gene sets and reshape it into the
## two-column TERM2GENE format clusterProfiler expects (gs_name, gene_symbol)
get_cancer_hallmarks <- function(geneset = c("core", "integrated")) {
  geneset <- match.arg(geneset)
  lines <- readLines(HALLMARK_URLS[[geneset]], warn = FALSE)
  bind_rows(lapply(lines, function(ln) {
    tok <- strsplit(ln, "\t")[[1]]
    tok <- tok[nzchar(tok)]
    tibble(gs_name = tok[1], gene_symbol = tok[-1])
  })) |> distinct()
}

t2g <- get_cancer_hallmarks("core")
cat(sprintf("[Cancer hallmarks: core] %d hallmark categories, %d unique genes\n",
            n_distinct(t2g$gs_name), n_distinct(t2g$gene_symbol)))

## -- Generic GSEA runner. `deg` needs a gene-symbol column (`gene_col`) and
## a numeric ranking statistic (`stat_col`); ties on gene symbol keep the
## first row. maxGSSize is raised well above clusterProfiler's 500 default
## because these hallmark categories are large (up to ~3,574 genes for the
## "integrated" set). Returns (and saves) the full GSEA result table.
run_gsea_hallmarks <- function(deg, gene_col = "gene", stat_col = "t", label = "run") {
  rnk <- deg |>
    filter(!is.na(.data[[stat_col]])) |>
    distinct(.data[[gene_col]], .keep_all = TRUE) |>
    arrange(desc(.data[[stat_col]])) |>
    (\(d) setNames(d[[stat_col]], d[[gene_col]]))()

  set.seed(123)
  gs <- GSEA(geneList = rnk, TERM2GENE = t2g, pvalueCutoff = 1, eps = 0, seed = TRUE,
             pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 5000)

  res <- as.data.frame(gs) |>
    as_tibble() |>
    mutate(run = label) |>
    arrange(p.adjust)

  write_tsv(res, file.path(out_dir, sprintf("04_gsea_hallmarks_%s.tsv", label)))
  cat(sprintf("[%s] %d/%d hallmark categories significant (adj.P<0.05)\n",
              label, sum(res$p.adjust < 0.05), nrow(res)))
  res
}

## ---- example (edit and uncomment to run) ----
# deg <- readr::read_tsv("~/path/to/your_deg_table.tsv")
# run_gsea_hallmarks(deg, gene_col = "gene", stat_col = "t", label = "my_comparison")
