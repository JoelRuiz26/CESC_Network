############################################################
## 6_1_DepMap_download.R
##
## Descarga metadata y scores CRISPR de DepMap para
## líneas celulares de cáncer cervicouterino HPV A7/A9.
##
## Filtro 1: primary_disease == "Cervical Cancer"
## Filtro 2: HPV genotipo verificado A7 o A9 (Cellosaurus /
##           literatura primaria; ver tabla hpv_annot abajo)
##
## Excluidas:
##   - Líneas HPV-negativas: C-33A, HT-3, HCA-1
##   - HPV Alpha-6 (HPV-56): HOKUG
##   - HPV desconocido: OMC-1
##   - Origen no cervical: BGC823, GLC82, SF767, SW954
##
## Outputs (en este mismo directorio):
##   6_1_cervix_metadata.tsv      — metadata de las líneas
##   6_1_cervix_dep_mean.tsv      — score medio global (todas las líneas)
##   6_1_cervix_dep_mean_A7.tsv   — score medio por gen, líneas A7
##   6_1_cervix_dep_mean_A9.tsv   — score medio por gen, líneas A9
##   6_1_crispr_lines.tsv         — líneas con datos CRISPR y su grupo
##   6_1_common_essential_genes.tsv — set pan-esencial oficial DepMap
##                                    (flag common_essential, Chronos_Combined)
############################################################

suppressPackageStartupMessages({
  library(depmap)
  library(ExperimentHub)
  library(dplyr)
  library(readr)
})

out_dir <- "~/CESC_Network/07_dependency"
res_dir <- file.path(out_dir, "results")           # outputs kept apart from scripts
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

# ── 1) Tabla de anotación HPV (verificada por Cellosaurus / literatura) ──
# Fuentes: Cellosaurus, Skok et al. 2021 (PMC8051715),
#          Fernández-Martínez et al. 2023 (PMC10417116)
hpv_annot <- tribble(
  ~stripped_cell_line_name, ~hpv_genotype, ~hpv_group,
  "HELA",        "HPV-18", "A7",
  "ME180",       "HPV-68", "A7",
  "MS751",       "HPV-45", "A7",
  "SW756",       "HPV-18", "A7",
  "C4I",         "HPV-18", "A7",
  "C4II",        "HPV-18", "A7",
  "SISO",        "HPV-18", "A7",
  "HCS2",        "HPV-18", "A7",
  "HCSC1",       "HPV-18", "A7",
  "SIHA",        "HPV-16", "A9",
  "CASKI",       "HPV-16", "A9",
  "DOTC24510",   "HPV-16", "A9",
  "SKGIIIA",     "HPV-16", "A9",
  "BOKU",        "HPV-16", "A9",
  "TCYIK",       "HPV-16", "A9"
)

cat("Líneas HPV A7/A9 a retener:", nrow(hpv_annot), "\n")

# ── 2) Metadata ──────────────────────────────────────────────
meta <- depmap_metadata()

cervix_meta <- meta %>%
  filter(grepl("cervical", primary_disease, ignore.case = TRUE)) %>%
  inner_join(hpv_annot, by = "stripped_cell_line_name")

cat("\nLíneas cervicales HPV A7/A9 encontradas en DepMap:", nrow(cervix_meta), "\n")
print(cervix_meta[, c("depmap_id", "cell_line", "hpv_genotype", "hpv_group")])

cervix_ids <- cervix_meta$depmap_id
ids_A7     <- cervix_meta %>% filter(hpv_group == "A7") %>% pull(depmap_id)
ids_A9     <- cervix_meta %>% filter(hpv_group == "A9") %>% pull(depmap_id)

# ── 3) CRISPR dependency scores ──────────────────────────────
crispr_raw <- depmap_crispr() %>%
  filter(depmap_id %in% cervix_ids) %>%
  mutate(gene_u = toupper(gene_name)) %>%
  rename(dep_score = dependency)

cat("\nFilas CRISPR:", nrow(crispr_raw), "\n")
cat("Líneas con datos CRISPR:", n_distinct(crispr_raw$cell_line), "\n")

# Líneas con datos CRISPR, por grupo
crispr_lines <- crispr_raw %>%
  distinct(depmap_id, cell_line) %>%
  left_join(cervix_meta %>% select(depmap_id, hpv_group), by = "depmap_id")

lines_A7 <- crispr_lines %>% filter(hpv_group == "A7") %>% pull(cell_line)
lines_A9 <- crispr_lines %>% filter(hpv_group == "A9") %>% pull(cell_line)

cat("\nA7 lines con CRISPR (n =", length(lines_A7), "):", paste(lines_A7, collapse = ", "), "\n")
cat("A9 lines con CRISPR (n =", length(lines_A9), "):", paste(lines_A9, collapse = ", "), "\n")

# Líneas en metadata sin datos CRISPR
missing_crispr <- cervix_meta %>% filter(!depmap_id %in% crispr_raw$depmap_id)
if (nrow(missing_crispr) > 0) {
  cat("\nLíneas sin datos CRISPR en este release:\n")
  print(missing_crispr[, c("cell_line", "hpv_genotype", "hpv_group")])
}

# ── 4) Score promedio por gen ─────────────────────────────────
# Global
dep_mean <- crispr_raw %>%
  group_by(gene_u) %>%
  summarise(
    mean_dep   = mean(dep_score, na.rm = TRUE),
    median_dep = median(dep_score, na.rm = TRUE),
    n_lines    = sum(!is.na(dep_score)),
    .groups    = "drop"
  )

# Por grupo alpha
dep_mean_A7 <- crispr_raw %>%
  filter(depmap_id %in% ids_A7) %>%
  group_by(gene_u) %>%
  summarise(mean_dep = mean(dep_score, na.rm = TRUE), .groups = "drop")

dep_mean_A9 <- crispr_raw %>%
  filter(depmap_id %in% ids_A9) %>%
  group_by(gene_u) %>%
  summarise(mean_dep = mean(dep_score, na.rm = TRUE), .groups = "drop")

cat("\nGenes únicos con score (global):", nrow(dep_mean), "\n")
cat("Genes únicos con score (A7):", nrow(dep_mean_A7), "\n")
cat("Genes únicos con score (A9):", nrow(dep_mean_A9), "\n")

# ── 4b) Set pan-esencial oficial de DepMap ───────────────────
# Flag "common_essential" del resumen por gen, dataset Chronos_Combined.
# IMPORTANTE: depmap_crispr() carga el CRISPR gene effect INTEGRADO de
# DepMap (CRISPR_gene_effect.csv, ~1086 lineas), que en gene_summary
# corresponde a "Chronos_Combined", NO a "Chronos_Achilles" (~941 lineas).
# Se toma el common_essential del MISMO dataset de los scores para que el
# corte pan-esencial sea consistente con la dependencia usada arriba.
# NO se re-deriva por nosotros. Requerido por 6_5 y 6_6.
common_ess <- depmap_gene_summary() %>%
  filter(dataset == "Chronos_Combined", common_essential) %>%
  transmute(gene_u = toupper(gene_name)) %>%
  distinct()

cat("\nGenes pan-esenciales (Chronos_Combined):", nrow(common_ess), "\n")

# ── 5) Guardar outputs ───────────────────────────────────────
write_tsv(cervix_meta,  file.path(res_dir, "6_1_cervix_metadata.tsv"))
write_tsv(dep_mean,     file.path(res_dir, "6_1_cervix_dep_mean.tsv"))
write_tsv(dep_mean_A7,  file.path(res_dir, "6_1_cervix_dep_mean_A7.tsv"))
write_tsv(dep_mean_A9,  file.path(res_dir, "6_1_cervix_dep_mean_A9.tsv"))
write_tsv(crispr_lines, file.path(res_dir, "6_1_crispr_lines.tsv"))
write_tsv(common_ess,   file.path(res_dir, "6_1_common_essential_genes.tsv"))

cat("\nOutputs guardados en:", out_dir, "\n")
cat("  6_1_cervix_metadata.tsv\n")
cat("  6_1_cervix_dep_mean.tsv\n")
cat("  6_1_cervix_dep_mean_A7.tsv\n")
cat("  6_1_cervix_dep_mean_A9.tsv\n")
cat("  6_1_crispr_lines.tsv\n")
cat("  6_1_common_essential_genes.tsv\n")
cat("✓ Descarga completa.\n")
