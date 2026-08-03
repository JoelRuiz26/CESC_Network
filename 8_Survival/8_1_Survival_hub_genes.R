### ============================================================
### 8_1_Survival_hub_genes.R
###
### Prognostic value of the SELECTIVE essential hub genes for overall survival
### in the full cervical-cancer cohort (both HPV clades pooled = "el total").
###
### Genes: the 90 unique selective essential hubs (A7:47, A9:54, shared:11),
### defined exactly as in 6_14 / 6_15.
###
### Endpoint: overall survival. Time and vital status come from the clinical
### table the user manually audited against TCGA + Almorox (8_3_Clinical_mixed),
### NOT from any inferred status.
###
### Design (classical survival analysis, no machine learning):
###   A. UNIVARIATE screen: Cox of each gene alone. Reports any survival role,
###      reported whether or not it survives multiple testing.
###   B. CLINICAL multivariable base model built from the established clinical
###      determinants of cervical-cancer survival: FIGO stage, age, histology.
###      HPV clade is deliberately NOT included as a covariate.
###   C. INCREMENTAL value of each gene, individually: clinical base + gene,
###      tested against the clinical base by likelihood-ratio test (nested
###      models). Does the gene add prognostic information beyond the clinic?
###
### Expression = log2(TPM+1), z-scored per gene, so every HR is per 1 SD.
### Multiple testing: Benjamini-Hochberg FDR across the 90-gene family.
###
### Outputs (8_Survival/results/):
###   8_1_univariate_cox.tsv                per-gene univariate Cox + FDR
###   8_1_forest_univariate_nominal.png/pdf forest of nominal (p<0.05) genes
###   8_1_clinical_base_model.txt           clinical multivariable model summary
###   8_1_gene_incremental_over_clinical.tsv per-gene adjusted HR + LRT + FDR
###   8_1_forest_incremental_nominal.png/pdf forest of genes adding value
###   8_1_cohort_and_audit.txt              cohort flow + data-audit log
### ============================================================
suppressPackageStartupMessages({
  library(tidyverse); library(vroom); library(igraph)
  library(survival); library(survminer)
})

base    <- "~/CESC_Network"
prep    <- file.path(base, "2_Prepro_TCGA_GTEx")
res_dir <- file.path(base, "6_DepMap_cacu/results")
out_dir <- file.path(base, "8_Survival/results")
clin_f  <- "/home/joelr/backup_sefi/8_Survival/8_3_Clinical_mixed.tsv"   # user-audited
cmx_f   <- file.path(base, "1_Get_Data/1_Output_rds/2_0_TCGA_CESC_clinicalMatrix_full.rds")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

col_A7 <- "#4E79A7"; col_A9 <- "#CB6A5C"; col_shared <- "#8C79B5"
clade_col <- c(`A7-only` = col_A7, `A9-only` = col_A9, shared = col_shared)
audit <- character(0)
log <- function(...) { m <- sprintf(...); audit <<- c(audit, m); cat(m, "\n") }

### ---- 1. selective essential hub gene list (same rule as 6_14 / 6_15) ------
elbow <- function(v) { n <- length(v); x <- (seq_len(n) - 1) / (n - 1)
  y <- (v - v[1]) / (v[n] - v[1]); v[which.max(abs(y - x))] }
tab <- vroom(file.path(res_dir, "6_2_gene_dep_degree_table.tsv"), show_col_types = FALSE) %>%
  filter(!is.na(degree), !is.na(mean_dep)) %>% distinct(network, gene, .keep_all = TRUE)
th <- tab %>% group_by(network) %>%
  summarise(k = elbow(sort(degree, decreasing = TRUE)),
            d = elbow(sort(mean_dep)), .groups = "drop")
ess <- tab %>% left_join(th, by = "network") %>%
  filter(degree >= k, mean_dep <= d) %>% mutate(gene = toupper(gene))
pan <- toupper(vroom(file.path(res_dir, "6_1_common_essential_genes.tsv"),
                     delim = "\t", show_col_types = FALSE)$gene_u)
sel  <- ess %>% filter(!gene %in% pan)
sel7 <- sort(unique(sel$gene[sel$network == "HPV-A7"]))
sel9 <- sort(unique(sel$gene[sel$network == "HPV-A9"]))
shared <- intersect(sel7, sel9); genes <- union(sel7, sel9)
clade_of <- setNames(ifelse(genes %in% shared, "shared",
                     ifelse(genes %in% sel7, "A7-only", "A9-only")), genes)
log("Selective essential hubs -> A7:%d A9:%d shared:%d union tested:%d",
    length(sel7), length(sel9), length(shared), length(genes))

### ---- 2. clinical / survival table (VALIDATED endpoint) --------------------
cm  <- vroom(clin_f, show_col_types = FALSE)
fac <- vroom(file.path(prep, "2_4_Factors.tsv"), show_col_types = FALSE)
cmx <- readRDS(cmx_f)
age_tbl <- cmx %>%
  transmute(patient = bcr_patient_barcode,
            age = suppressWarnings(as.numeric(age_at_initial_pathologic_diagnosis))) %>%
  distinct(patient, .keep_all = TRUE)

recon <- fac %>% transmute(patient = bcr_patient_barcode,
                           fac_status = ifelse(vital_status == "DECEASED", "Dead", "Alive")) %>%
  left_join(cm %>% select(case_submitter_id, vital_status),
            by = c("patient" = "case_submitter_id"))
disc <- recon %>% filter(fac_status != vital_status)
log("Vital-status agreement Factors vs validated clinical: %d/%d; %d discrepant (validated used).",
    sum(recon$fac_status == recon$vital_status), nrow(recon), nrow(disc))

surv <- fac %>%
  transmute(sampleID, patient = bcr_patient_barcode,
            clade = recode(HPV_clade, A7_clade = "A7", A9_clade = "A9"),
            stage = na_if(Stage_group, "(NOT REPORTED)"),
            histology = histological_group) %>%
  left_join(cm %>% select(patient = case_submitter_id, vital_status,
                          days_to_death, days_to_last_follow_up), by = "patient") %>%
  left_join(age_tbl, by = "patient") %>%
  mutate(event   = vital_status == "Dead",
         os_time = if_else(event, days_to_death, days_to_last_follow_up))

n0 <- nrow(surv)
surv <- surv %>% filter(!is.na(os_time), os_time > 0, !is.na(event))
surv$stage     <- factor(surv$stage, levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))
surv$histology <- factor(surv$histology, levels = c("Squamous", "Adenocarcinoma", "Adenosquamous"))
surv$age10     <- surv$age / 10                       # HR per 10 years
log("Cohort: %d tumours -> %d with valid OS(>0). Deaths: %d. Age NA:%d Stage NA:%d.",
    n0, nrow(surv), sum(surv$event), sum(is.na(surv$age)), sum(is.na(surv$stage)))

### ---- 3. expression: log2(TPM+1), z-scored per gene ------------------------
expr <- vroom(file.path(prep, "2_6_Full_counts_A7_A9_annot.tsv"), show_col_types = FALSE)
stopifnot(all(genes %in% expr$Gene))
em <- expr %>% filter(Gene %in% genes) %>% column_to_rownames("Gene") %>% as.matrix()
em <- em[genes, surv$sampleID, drop = FALSE]
stopifnot(identical(colnames(em), surv$sampleID))
logtpm <- log2(em + 1)
zmat   <- t(scale(t(logtpm)))                          # genes x samples, z per gene
Sy <- Surv(surv$os_time, surv$event)

### ============================================================
### A. UNIVARIATE Cox screen (every gene alone)
### ============================================================
hr_ci <- function(m, term = "z") {
  s <- summary(m); co <- s$coefficients
  if (!term %in% rownames(co)) return(c(hr = NA, lo = NA, hi = NA, p = NA))
  ci <- s$conf.int[term, c("lower .95", "upper .95")]
  c(hr = unname(co[term, "exp(coef)"]), lo = ci[[1]], hi = ci[[2]],
    p = unname(co[term, "Pr(>|z|)"]))
}
uni <- map_dfr(genes, function(g) {
  d <- surv; d$z <- as.numeric(zmat[g, ]); m <- coxph(Sy ~ z, data = d)
  a <- hr_ci(m)
  tibble(gene = g, clade = clade_of[g], n = m$n, nevent = m$nevent,
         HR = a["hr"], lo = a["lo"], hi = a["hi"], p = a["p"],
         ph_p = tryCatch(cox.zph(m)$table["z", "p"], error = function(e) NA_real_))
}) %>% mutate(FDR = p.adjust(p, "BH")) %>% arrange(p)
write_tsv(uni, file.path(out_dir, "8_1_univariate_cox.tsv"))
uni_nom <- uni %>% filter(p < 0.05)
log("A. Univariate: %d/%d genes p<0.05; FDR<0.05: %d (min FDR=%.2f).",
    nrow(uni_nom), nrow(uni), sum(uni$FDR < 0.05), min(uni$FDR))

forest <- function(df, file, title, sub) {
  df <- df %>% arrange(HR) %>% mutate(gene = factor(gene, levels = gene))
  g <- ggplot(df, aes(HR, gene, colour = clade)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey55") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, linewidth = 0.7) +
    geom_point(size = 2.8) + scale_colour_manual(values = clade_col, name = "hub type") +
    scale_x_log10() +
    labs(x = "Hazard ratio per SD of log2 expression (95% CI)", y = NULL,
         title = title, subtitle = sub) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
  ggsave(paste0(file, ".png"), g, width = 8, height = 0.36 * nrow(df) + 1.9, dpi = 400)
  ggsave(paste0(file, ".pdf"), g, width = 8, height = 0.36 * nrow(df) + 1.9)
}
if (nrow(uni_nom))
  forest(uni_nom, file.path(out_dir, "8_1_forest_univariate_nominal"),
         "Univariate survival association of selective essential hubs",
         sprintf("%d/%d genes with uncorrected p<0.05; min BH-FDR = %.2f (none <0.05).",
                 nrow(uni_nom), nrow(uni), min(uni$FDR)))

### ============================================================
### B. CLINICAL multivariable base model (stage + age + histology; no clade)
### ============================================================
cc <- surv %>% filter(!is.na(stage), !is.na(age), !is.na(histology))
Sc <- Surv(cc$os_time, cc$event)
clin_base <- coxph(Sc ~ stage + age10 + histology, data = cc)
sc <- summary(clin_base)
log("B. Clinical base model: n=%d, events=%d. Concordance=%.3f. LR p=%.2g.",
    clin_base$n, clin_base$nevent, sc$concordance[1], sc$logtest["pvalue"])
cb_txt <- capture.output({
  cat("Clinical multivariable base model for overall survival\n")
  cat("Covariates: FIGO stage (ref Stage I), age per 10y, histology (ref Squamous).\n")
  cat("HPV clade intentionally excluded.\n\n"); print(sc)
  cat("\nProportional hazards (cox.zph):\n"); print(cox.zph(clin_base)$table)
})
writeLines(cb_txt, file.path(out_dir, "8_1_clinical_base_model.txt"))

### ============================================================
### C. INCREMENTAL value of each gene over the clinical base (individually)
###    nested LRT: clinical  vs  clinical + gene
### ============================================================
zc <- zmat[, cc$sampleID, drop = FALSE]
inc <- map_dfr(genes, function(g) {
  d <- cc; d$z <- as.numeric(zc[g, ])
  m0 <- coxph(Sc ~ stage + age10 + histology, data = d)
  m1 <- coxph(Sc ~ stage + age10 + histology + z, data = d)
  a  <- hr_ci(m1, "z")
  lrt <- anova(m0, m1, test = "LRT")
  tibble(gene = g, clade = clade_of[g],
         HR = a["hr"], lo = a["lo"], hi = a["hi"], p_adj = a["p"],
         lrt_p = lrt$`Pr(>|Chi|)`[2],
         dConcordance = summary(m1)$concordance[1] - summary(m0)$concordance[1])
}) %>% mutate(fdr_adj = p.adjust(p_adj, "BH"),
             fdr_lrt = p.adjust(lrt_p, "BH")) %>% arrange(p_adj)
write_tsv(inc, file.path(out_dir, "8_1_gene_incremental_over_clinical.tsv"))
inc_nom <- inc %>% filter(p_adj < 0.05)
log("C. Incremental over clinical: %d/%d genes adjusted p<0.05; FDR<0.05: %d (min FDR=%.2f).",
    nrow(inc_nom), nrow(inc), sum(inc$fdr_adj < 0.05), min(inc$fdr_adj))
if (nrow(inc_nom))
  forest(inc_nom %>% rename(p = p_adj), file.path(out_dir, "8_1_forest_incremental_nominal"),
         "Selective hubs adding prognostic value over clinical factors",
         sprintf("Cox adjusted for FIGO stage + age + histology; %d/%d genes p<0.05, min BH-FDR=%.2f.",
                 nrow(inc_nom), nrow(inc), min(inc$fdr_adj)))

### ---- audit log ------------------------------------------------------------
hdr <- c("Prognostic analysis of selective essential hub genes (overall survival)",
  "Cohort: full A7/A9 cervical-cancer tumour set, both clades pooled.",
  sprintf("Endpoint from user-audited file: %s", clin_f),
  "Covariates from 2_4_Factors.tsv (stage/histology) + clinicalMatrix (age).",
  "Expression: batch-corrected TPM, log2(TPM+1), z per gene; HR per 1 SD.",
  "FDR = Benjamini-Hochberg across the 90-gene family.",
  strrep("-", 70))
writeLines(c(hdr, audit), file.path(out_dir, "8_1_cohort_and_audit.txt"))
message("\nDone. Results in ", out_dir)
