### Joel Ruiz
### Is FIGO stage or histological group associated with HPV clade (A7 vs A9)?
### Chi-squared, or Fisher's exact if any expected cell < 5. Same check and
### same stage/histology grouping as the OCTAD drug-repositioning paper.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(vroom) })

meta <- readRDS("~/CESC_Network/02_cohort_data/1_Output_rds/2_2_CESC_TCGA_GTEx_metadata_clade.rds") %>%
  filter(HPV_clade %in% c("A7_clade", "A9_clade"), Stage_group != "(NOT REPORTED)") %>%
  mutate(HPV_clade = sub("_clade", "", HPV_clade))

test_feature <- function(var) {
  tab <- table(meta$HPV_clade, meta[[var]])
  chi <- suppressWarnings(chisq.test(tab))
  use_fisher <- any(chi$expected < 5)
  p <- if (use_fisher) fisher.test(tab)$p.value else chi$p.value
  v <- sqrt(unname(chi$statistic) / (sum(tab) * (min(dim(tab)) - 1)))
  as.data.frame(tab) |>
    rename(HPV_clade = Var1, Level = Var2, n = Freq) |>
    group_by(HPV_clade) |> mutate(pct = 100 * n / sum(n)) |> ungroup() |>
    mutate(cell = sprintf("%d (%.1f%%)", n, pct)) |>
    pivot_wider(id_cols = Level, names_from = HPV_clade, values_from = c(n, cell)) |>
    transmute(Feature = var, Level, A7 = cell_A7, A9 = cell_A9,
             Total = n_A7 + n_A9,
             method = if (use_fisher) "Fisher's exact" else "Chi-squared",
             p_value = signif(p, 3), cramers_v = signif(v, 3))
}

result <- bind_rows(test_feature("Stage_group"), test_feature("histological_group"))
vroom_write(result, "~/CESC_Network/02_cohort_data/1_Output_rds/3_1_fisher_clinical_vs_clade.tsv", delim = "\t")
print(as.data.frame(result), row.names = FALSE)
