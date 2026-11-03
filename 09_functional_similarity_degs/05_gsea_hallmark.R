### ============================================================
### 05_gsea_hallmark.R
###
### GSEA (rank-based, not ORA) against MSigDB Hallmark (H, 50 gene sets),
### run SEPARATELY on the full A7-vs-A9 DEG ranking of each cohort (CESC,
### HTMCP-CC) -- genes ranked by limma's moderated t-statistic (standard
### GSEA input), using the same Hallmark collection as
### 06_modules/02_hallmark_enrichment.R but here testing the whole ranked
### gene list rather than per-module gene sets.
###
### -> results/05_gsea_hallmark_{cesc,htmcp}.tsv
### -> figures/05_gsea_hallmark_barplot.png / .pdf
### ============================================================
suppressPackageStartupMessages({
  library(tidyverse); library(vroom); library(clusterProfiler); library(msigdbr); library(patchwork)
})

out_dir <- "~/CESC_Network/09_functional_similarity_degs"
FF <- "Times New Roman"
set.seed(123)

t2g <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  select(gs_name, gene_symbol) %>% distinct()

full_cesc  <- vroom(file.path(out_dir, "results/01_deg_full_cesc.tsv"),  show_col_types = FALSE)
full_htmcp <- vroom(file.path(out_dir, "results/01_deg_full_htmcp.tsv"), show_col_types = FALSE)

run_gsea <- function(tt, label) {
  rnk <- tt %>% filter(!is.na(t)) %>% distinct(gene, .keep_all = TRUE) %>%
    arrange(desc(t)) %>% { setNames(.$t, .$gene) }
  set.seed(123)
  gs <- GSEA(geneList = rnk, TERM2GENE = t2g, pvalueCutoff = 1, eps = 0, seed = TRUE,
            pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500)
  as.data.frame(gs) %>% as_tibble() %>%
    mutate(HM = gsub("_", " ", sub("^HALLMARK_", "", ID)), cohort = label) %>%
    arrange(p.adjust)
}

gsea_cesc  <- run_gsea(full_cesc,  "CESC")
gsea_htmcp <- run_gsea(full_htmcp, "HTMCP-CC")
write_tsv(gsea_cesc,  file.path(out_dir, "results/05_gsea_hallmark_cesc.tsv"))
write_tsv(gsea_htmcp, file.path(out_dir, "results/05_gsea_hallmark_htmcp.tsv"))

cat(sprintf("Significant Hallmark terms (adj.P<0.05): CESC=%d, HTMCP-CC=%d\n",
            sum(gsea_cesc$p.adjust < 0.05), sum(gsea_htmcp$p.adjust < 0.05)))
shared_terms <- intersect(gsea_cesc$HM[gsea_cesc$p.adjust < 0.05], gsea_htmcp$HM[gsea_htmcp$p.adjust < 0.05])
cat(sprintf("Shared significant terms: %d\n", length(shared_terms))); print(shared_terms)

## ── barplot: significant terms (adj.P<0.05) in EITHER cohort, NES on x,
## coloured by direction (NES>0 = enriched in A7, NES<0 = enriched in A9),
## faceted by cohort, ordered by CESC's NES for visual alignment
COL_A7 <- "#2166AC"; COL_A9 <- "#B2182B"

sig_terms <- union(gsea_cesc$HM[gsea_cesc$p.adjust < 0.05], gsea_htmcp$HM[gsea_htmcp$p.adjust < 0.05])
both <- bind_rows(gsea_cesc, gsea_htmcp) %>% filter(HM %in% sig_terms) %>%
  mutate(direction = if_else(NES > 0, "Enriched in A7", "Enriched in A9"),
         cohort = factor(cohort, c("CESC", "HTMCP-CC")))
term_order <- gsea_cesc %>% filter(HM %in% sig_terms) %>% arrange(NES) %>% pull(HM)
missing_in_cesc <- setdiff(sig_terms, term_order)
term_order <- c(missing_in_cesc, term_order)
both <- both %>% mutate(HM = factor(HM, term_order))

p <- ggplot(both, aes(NES, HM, fill = direction)) +
  geom_col(width = 0.72, colour = "grey25", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.4) +
  facet_wrap(~cohort, nrow = 1) +
  scale_fill_manual(values = c("Enriched in A7" = COL_A7, "Enriched in A9" = COL_A9), name = NULL) +
  labs(title = "GSEA (Hallmark) — HPV-A7 vs HPV-A9, CESC and HTMCP-CC",
       x = "Normalized Enrichment Score (NES)", y = NULL) +
  theme_bw(base_family = FF) +
  theme(strip.background = element_rect(fill = "#16213E", colour = NA),
        strip.text = element_text(colour = "white", size = 15, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold", colour = "grey10"),
        axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
        panel.border = element_rect(colour = "grey30", linewidth = 0.7),
        panel.spacing = unit(1.2, "lines"),
        legend.position = "bottom", legend.text = element_text(size = 13))

ggsave(file.path(out_dir, "figures/05_gsea_hallmark_barplot.png"), p, width = 13,
       height = max(6, 0.4 * length(sig_terms)), dpi = 300, bg = "white")
ggsave(file.path(out_dir, "figures/05_gsea_hallmark_barplot.pdf"), p, width = 13,
       height = max(6, 0.4 * length(sig_terms)), device = cairo_pdf, bg = "white")
message("Saved: figures/05_gsea_hallmark_barplot.png / .pdf")
