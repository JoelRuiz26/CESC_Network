### ============================================================
### 03_shared_degs_and_lfc_correlation.R
###
### Identity of A7-vs-A9 DEGs shared between CESC and HTMCP-CC (from
### 01_deg_a7_vs_a9.R), and the correlation of their logFC between cohorts
### -- the key evidence for "both cohorts look functionally similar":
### if the same genes come up as A7-vs-A9 DEGs in both, AND their fold
### changes point the same direction and scale together, that is a
### concordant transcriptional signature independent of cohort.
###
### -> results/03_shared_degs.tsv
### -> results/03_lfc_correlation_stats.tsv
### -> figures/03_lfc_correlation_scatter.png/.pdf
### ============================================================
suppressPackageStartupMessages({
  library(tidyverse); library(vroom); library(ggrepel)
})

out_dir <- "~/CESC_Network/09_functional_similarity_degs"
FF <- "Times New Roman"

full_cesc  <- vroom(file.path(out_dir, "results/01_deg_full_cesc.tsv"),  show_col_types = FALSE)
full_htmcp <- vroom(file.path(out_dir, "results/01_deg_full_htmcp.tsv"), show_col_types = FALSE)
sig_cesc   <- vroom(file.path(out_dir, "results/01_deg_significant_cesc.tsv"),  show_col_types = FALSE)
sig_htmcp  <- vroom(file.path(out_dir, "results/01_deg_significant_htmcp.tsv"), show_col_types = FALSE)

shared_genes <- intersect(sig_cesc$gene, sig_htmcp$gene)
shared <- inner_join(sig_cesc %>% select(gene, logFC_cesc = logFC, padj_cesc = adj.P.Val),
                     sig_htmcp %>% select(gene, logFC_htmcp = logFC, padj_htmcp = adj.P.Val),
                     by = "gene") %>%
  mutate(concordant = sign(logFC_cesc) == sign(logFC_htmcp)) %>%
  arrange(desc(concordant), padj_cesc)

write_tsv(shared, file.path(out_dir, "results/03_shared_degs.tsv"))
cat(sprintf("DEGs compartidos entre CESC (n=%d) y HTMCP-CC (n=%d): %d\n", nrow(sig_cesc), nrow(sig_htmcp), nrow(shared)))
cat(sprintf("De estos, %d (%.1f%%) tienen el MISMO sentido (misma dirección del cambio) en ambas cohortes\n",
            sum(shared$concordant), 100 * mean(shared$concordant)))

## correlation of logFC for the shared DEGs, using each cohort's full
## topTable (not just the significant subset) so the correlation isn't
## artificially inflated by the significance filter itself
lfc_all <- inner_join(full_cesc %>% select(gene, logFC_cesc = logFC),
                      full_htmcp %>% select(gene, logFC_htmcp = logFC), by = "gene")
lfc_shared <- lfc_all %>% filter(gene %in% shared_genes)

ct_all    <- cor.test(lfc_all$logFC_cesc, lfc_all$logFC_htmcp, method = "spearman")
ct_shared <- cor.test(lfc_shared$logFC_cesc, lfc_shared$logFC_htmcp, method = "spearman")

stats_tbl <- tibble(
  set = c("Todos los genes en común", "DEGs compartidos (A7 vs A9, ambas cohortes)"),
  n = c(nrow(lfc_all), nrow(lfc_shared)),
  spearman_rho = c(unname(ct_all$estimate), unname(ct_shared$estimate)),
  p_value = c(ct_all$p.value, ct_shared$p.value)
)
write_tsv(stats_tbl, file.path(out_dir, "results/03_lfc_correlation_stats.tsv"))
print(stats_tbl)

## ── scatter: logFC CESC vs logFC HTMCP-CC for the shared DEGs
lab <- shared %>% mutate(mag = abs(logFC_cesc) + abs(logFC_htmcp)) %>% slice_max(mag, n = 15)
rho_lab <- sprintf("Spearman ρ = %.3f\np = %.1e\nn = %d genes", ct_shared$estimate, ct_shared$p.value, nrow(lfc_shared))

p <- ggplot(shared, aes(logFC_cesc, logFC_htmcp)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.5) +
  geom_point(aes(colour = concordant), size = 2.6, alpha = 0.85) +
  geom_text_repel(data = lab, aes(label = gene), family = FF, size = 4, fontface = "bold",
                  max.overlaps = Inf, seed = 1, box.padding = 0.4, segment.size = 0.25, segment.colour = "grey50") +
  annotate("label", x = min(shared$logFC_cesc), y = max(shared$logFC_htmcp), hjust = 0, vjust = 1,
           label = rho_lab, family = FF, fontface = "bold", size = 4.4, fill = alpha("white", 0.85),
           label.padding = unit(0.35, "lines")) +
  scale_colour_manual(values = c(`TRUE` = "#2C7BB6", `FALSE` = "#D7191C"),
                      labels = c(`TRUE` = "Mismo sentido", `FALSE` = "Sentido opuesto"), name = "Concordancia") +
  labs(title = "DEGs HPV-A7 vs HPV-A9 compartidos entre cohortes: correlación de logFC",
       subtitle = "CESC (TCGA) vs HTMCP-CC (Uganda), cada una con su propio análisis A7-vs-A9",
       x = "logFC (A7 vs A9) — CESC", y = "logFC (A7 vs A9) — HTMCP-CC") +
  theme_bw(base_family = FF) +
  theme(plot.title = element_text(size = 16, face = "bold"), plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 13, face = "bold"), axis.text = element_text(size = 11),
        legend.position = "bottom", legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11), panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "figures/03_lfc_correlation_scatter.png"), p, width = 9, height = 9, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "figures/03_lfc_correlation_scatter.pdf"), p, width = 9, height = 9, device = cairo_pdf, bg = "white")
message("Saved: figures/03_lfc_correlation_scatter.png / .pdf")
