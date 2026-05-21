############################################################
## 6_1_DepMap_download.R
##
## Descarga metadata y scores CRISPR de DepMap para
## líneas celulares de cáncer cervicouterino.
##
## Filtro: primary_disease == "Cervical Cancer" (DepMap)
## Outputs (en este mismo directorio):
##   6_1_cervix_metadata.tsv   — metadata de las líneas
##   6_1_cervix_dep_mean.tsv   — score medio por gen
############################################################

suppressPackageStartupMessages({
  library(depmap)
  library(ExperimentHub)
  library(dplyr)
  library(readr)
})

out_dir <- "~/CESC_Network/6_DepMap_cacu"

# ── 1) Metadata ─────────────────────────────────────────────
meta <- depmap_metadata()

cat("Valores únicos en primary_disease (muestra):\n")
print(sort(unique(meta$primary_disease))[1:30])

cervix_meta <- meta %>%
  filter(grepl("cervical", primary_disease, ignore.case = TRUE))

cat("\nLíneas celulares de cáncer cervicouterino encontradas:", nrow(cervix_meta), "\n")
print(cervix_meta[, c("depmap_id", "cell_line", "primary_disease",
                      "sample_collection_site")])

cervix_ids <- cervix_meta$depmap_id

# ── 2) CRISPR dependency scores ──────────────────────────────
crispr_raw <- depmap_crispr() %>%
  filter(depmap_id %in% cervix_ids) %>%
  mutate(gene_u = toupper(gene_name)) %>%
  rename(dep_score = dependency)

cat("\nFilas CRISPR (líneas cervicales):", nrow(crispr_raw), "\n")
cat("Líneas celulares con datos CRISPR:", n_distinct(crispr_raw$cell_line), "\n")

# Score promedio por gen (media sobre líneas celulares disponibles)
dep_mean <- crispr_raw %>%
  group_by(gene_u) %>%
  summarise(
    mean_dep   = mean(dep_score, na.rm = TRUE),
    median_dep = median(dep_score, na.rm = TRUE),
    n_lines    = sum(!is.na(dep_score)),
    .groups    = "drop"
  )

cat("Genes únicos con score:", nrow(dep_mean), "\n")

# ── 3) Guardar outputs ───────────────────────────────────────
write_tsv(cervix_meta, file.path(out_dir, "6_1_cervix_metadata.tsv"))
write_tsv(dep_mean,    file.path(out_dir, "6_1_cervix_dep_mean.tsv"))

cat("\nOutputs guardados en:", out_dir, "\n")
cat("  6_1_cervix_metadata.tsv\n")
cat("  6_1_cervix_dep_mean.tsv\n")
cat("✓ Descarga completa.\n")
