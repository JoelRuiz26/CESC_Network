### ============================================================
### 01_deg_a7_vs_a9.R
###
### Differential expression HPV-A7 vs HPV-A9, computed SEPARATELY within
### each cohort (CESC and HTMCP-CC) -- never comparing cohorts directly
### against each other, to avoid batch effects. Uses limma on log2(x+1) of
### the already GC/length-bias-normalised counts from 03_preprocessing
### (03_preprocessing/2_7_/2_8_Full_counts_A{7,9}_annot.tsv for CESC, and
### the validation_uganda/ equivalents for HTMCP-CC).
###
### HTMCP-CC only: HIV status is a known, imbalanced confounder in this
### cohort (A7 = 34/52 HIV+, A9 = 33/57 HIV+) and is corrected for with
### ComBat (ComBat(..., mod = model.matrix(~group))) BEFORE limma, so the
### A7-vs-A9 biological signal is explicitly protected from removal. The
### single sample with HIV status "Unknown" is dropped (ComBat needs a
### clean 2-level batch). CESC has no HIV (or equivalent) covariate
### available in its clinical table, so no analogous correction applies
### there.
###
### DEG definition: |logFC| > 1 & adj.P.Val < 0.05 (BH), A7 vs A9 (positive
### logFC = higher in HPV-A7).
###
### -> results/01_deg_full_{cesc,htmcp}.tsv   (full topTable, all genes)
### -> results/01_deg_significant_{cesc,htmcp}.tsv (DEGs only)
### ============================================================
suppressPackageStartupMessages({
  library(tidyverse); library(vroom); library(limma); library(sva)
})

out_dir <- "~/CESC_Network/09_functional_similarity_degs"
LFC_CUT <- 1; P_CUT <- 0.05

read_expr <- function(path) {
  m <- vroom(path, show_col_types = FALSE)
  mat <- as.matrix(m[, -1]); rownames(mat) <- m[[1]]
  mat
}

run_deg <- function(logexpr, group, label) {
  design <- model.matrix(~group)
  fit <- eBayes(lmFit(logexpr, design))
  tt <- topTable(fit, coef = "groupA7", number = Inf, sort.by = "none") %>%
    as_tibble(rownames = "gene") %>%
    mutate(cohort = label)
  cat(sprintf("[%s] n_A7=%d n_A9=%d genes=%d\n", label, sum(group == "A7"), sum(group == "A9"), nrow(tt)))
  tt
}

## -- CESC: no covariate correction (no HIV/equivalent field available)
exprA7_cesc <- read_expr("~/CESC_Network/03_preprocessing/2_7_Full_counts_A7_annot.tsv")
exprA9_cesc <- read_expr("~/CESC_Network/03_preprocessing/2_8_Full_counts_A9_annot.tsv")
stopifnot(identical(rownames(exprA7_cesc), rownames(exprA9_cesc)))
logexpr_cesc <- log2(cbind(exprA7_cesc, exprA9_cesc) + 1)
group_cesc <- factor(c(rep("A7", ncol(exprA7_cesc)), rep("A9", ncol(exprA9_cesc))), levels = c("A9", "A7"))
deg_cesc <- run_deg(logexpr_cesc, group_cesc, "CESC")

## -- HTMCP-CC: ComBat-adjust for HIV status first, protecting the A7/A9
## signal via the `mod` argument
exprA7_htmcp <- read_expr("~/CESC_Network/03_preprocessing/validation_uganda/2_7_Full_counts_A7_annot.tsv")
exprA9_htmcp <- read_expr("~/CESC_Network/03_preprocessing/validation_uganda/2_8_Full_counts_A9_annot.tsv")
stopifnot(identical(rownames(exprA7_htmcp), rownames(exprA9_htmcp)))

clin_htmcp <- vroom("~/CESC_Network/02_cohort_data/validation_uganda/1_Output/1_1_HTMCP_clinical_HPV_HIV_table.tsv", show_col_types = FALSE)
all_samples <- c(colnames(exprA7_htmcp), colnames(exprA9_htmcp))
group_all <- factor(c(rep("A7", ncol(exprA7_htmcp)), rep("A9", ncol(exprA9_htmcp))), levels = c("A9", "A7"))
hiv <- clin_htmcp$HIV.status[match(all_samples, clin_htmcp$Patient)]
keep <- !is.na(hiv) & hiv %in% c("Positive", "Negative")
cat(sprintf("[HTMCP-CC] dropping %d sample(s) with missing/Unknown HIV status before ComBat\n", sum(!keep)))

logexpr_htmcp_raw <- log2(cbind(exprA7_htmcp, exprA9_htmcp)[, keep] + 1)
group_htmcp <- droplevels(group_all[keep])
hiv_htmcp <- factor(hiv[keep])
mod_htmcp <- model.matrix(~group_htmcp)
logexpr_htmcp <- ComBat(dat = logexpr_htmcp_raw, batch = hiv_htmcp, mod = mod_htmcp, par.prior = TRUE)
deg_htmcp <- run_deg(logexpr_htmcp, group_htmcp, "HTMCP-CC")

write_tsv(deg_cesc,  file.path(out_dir, "results/01_deg_full_cesc.tsv"))
write_tsv(deg_htmcp, file.path(out_dir, "results/01_deg_full_htmcp.tsv"))

sig_cesc  <- deg_cesc  %>% filter(adj.P.Val < P_CUT, abs(logFC) > LFC_CUT) %>% arrange(adj.P.Val)
sig_htmcp <- deg_htmcp %>% filter(adj.P.Val < P_CUT, abs(logFC) > LFC_CUT) %>% arrange(adj.P.Val)
write_tsv(sig_cesc,  file.path(out_dir, "results/01_deg_significant_cesc.tsv"))
write_tsv(sig_htmcp, file.path(out_dir, "results/01_deg_significant_htmcp.tsv"))

cat(sprintf("\nDEGs (|logFC|>%.0f, adj.P<%.2f): CESC=%d (up-A7=%d, up-A9=%d) | HTMCP-CC=%d (up-A7=%d, up-A9=%d)\n",
            LFC_CUT, P_CUT, nrow(sig_cesc), sum(sig_cesc$logFC > 0), sum(sig_cesc$logFC < 0),
            nrow(sig_htmcp), sum(sig_htmcp$logFC > 0), sum(sig_htmcp$logFC < 0)))
message("Saved: results/01_deg_full_{cesc,htmcp}.tsv, results/01_deg_significant_{cesc,htmcp}.tsv")
