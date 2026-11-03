### ============================================================
### 02_hallmark_enrichment.R
###
### Over-representation analysis (MSigDB Hallmark, 50 gene sets) on the
### A7-vs-A9 DEGs of each cohort separately (from 01_deg_a7_vs_a9.R),
### using the same universe/enricher() pattern as 06_modules/02_hallmark_
### enrichment.R for consistency. DEGs tested as one combined set per
### cohort (both up-in-A7 and up-in-A9 together), so the comparison speaks
### to "which biological programs distinguish A7 from A9" in each cohort.
###
### -> results/02_hallmark_ora_{cesc,htmcp}.tsv
### -> figures/02_hallmark_shared_dotplot.png/.pdf
### ============================================================
suppressPackageStartupMessages({
  library(tidyverse); library(vroom); library(clusterProfiler); library(msigdbr)
})

out_dir <- "~/CESC_Network/09_functional_similarity_degs"
P_CUT <- 0.05; Q_CUT <- 0.20
FF <- "Times New Roman"

t2g <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  select(gs_name, gene_symbol) %>% distinct()

sig_cesc  <- vroom(file.path(out_dir, "results/01_deg_significant_cesc.tsv"), show_col_types = FALSE)
sig_htmcp <- vroom(file.path(out_dir, "results/01_deg_significant_htmcp.tsv"), show_col_types = FALSE)

universo_cesc  <- read_rds("~/CESC_Network/03_preprocessing/2_1_Universe_annotation_table.rds") %>% pull(hgnc_symbol) %>% unique()
universo_htmcp <- read_rds("~/CESC_Network/03_preprocessing/validation_uganda/2_1_Universe_annotation_table_uganda.rds") %>% pull(hgnc_symbol) %>% unique()

run_ora <- function(sig_tbl, universo, label) {
  eg <- tryCatch(enricher(gene = sig_tbl$gene, universe = universo, TERM2GENE = t2g,
                          pvalueCutoff = P_CUT, qvalueCutoff = Q_CUT, pAdjustMethod = "BH"),
                 error = function(e) NULL)
  if (is.null(eg) || nrow(as.data.frame(eg)) == 0) return(tibble())
  as.data.frame(eg) %>% as_tibble() %>%
    mutate(HM = gsub("_", " ", sub("^HALLMARK_", "", ID)), cohort = label) %>%
    arrange(p.adjust)
}

ora_cesc  <- run_ora(sig_cesc,  universo_cesc,  "CESC")
ora_htmcp <- run_ora(sig_htmcp, universo_htmcp, "HTMCP-CC")
write_tsv(ora_cesc,  file.path(out_dir, "results/02_hallmark_ora_cesc.tsv"))
write_tsv(ora_htmcp, file.path(out_dir, "results/02_hallmark_ora_htmcp.tsv"))

cat(sprintf("Hallmark terms enriched among A7-vs-A9 DEGs: CESC=%d, HTMCP-CC=%d\n", nrow(ora_cesc), nrow(ora_htmcp)))
shared_terms <- intersect(ora_cesc$HM, ora_htmcp$HM)
cat(sprintf("Shared across both cohorts: %d\n", length(shared_terms)))
cat(paste(shared_terms, collapse = ", "), "\n")

## ── comparison dot-plot: all terms enriched in EITHER cohort, side by side
gr_num <- function(x) vapply(strsplit(x, "/"), function(p) as.numeric(p[1]) / as.numeric(p[2]), numeric(1))
both <- bind_rows(ora_cesc, ora_htmcp) %>%
  mutate(GeneRatio = gr_num(GeneRatio), neglogp = -log10(p.adjust))
term_order <- both %>% group_by(HM) %>% summarise(mx = max(neglogp), .groups = "drop") %>% arrange(desc(mx)) %>% pull(HM)
both <- both %>% mutate(HM = factor(HM, rev(term_order)), cohort = factor(cohort, c("CESC", "HTMCP-CC")),
                        shared = HM %in% shared_terms)

p <- ggplot(both, aes(cohort, HM)) +
  geom_point(aes(size = GeneRatio, colour = neglogp)) +
  scale_colour_viridis_c(option = "magma", begin = 0.1, end = 0.9, direction = -1, name = expression(-log[10]~italic(p)[adj])) +
  scale_size_continuous(range = c(2, 9), name = "Gene ratio") +
  labs(title = "Programas Hallmark enriquecidos entre los DEGs HPV-A7 vs HPV-A9",
       subtitle = sprintf("%d de %d términos mostrados son compartidos entre ambas cohortes", length(shared_terms), n_distinct(both$HM)),
       x = NULL, y = NULL) +
  theme_bw(base_family = FF) +
  theme(axis.text.y = element_text(size = 11, face = ifelse(levels(both$HM) %in% shared_terms, "bold", "plain")),
        axis.text.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"), plot.subtitle = element_text(size = 12),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 11))

ggsave(file.path(out_dir, "figures/02_hallmark_shared_dotplot.png"), p, width = 11, height = max(6, 0.28 * n_distinct(both$HM)), dpi = 300, bg = "white")
ggsave(file.path(out_dir, "figures/02_hallmark_shared_dotplot.pdf"), p, width = 11, height = max(6, 0.28 * n_distinct(both$HM)), device = cairo_pdf, bg = "white")
message("Saved: figures/02_hallmark_shared_dotplot.png / .pdf")
