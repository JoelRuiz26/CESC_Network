### Joel Ruiz
### Is FIGO stage, histological group or HIV status associated with HPV
### clade (A7 vs A9) in the Uganda (HTMCP-CC) validation cohort? Same check,
### same method as 02_cohort_data/03_fisher_clinical_vs_clade.R for CESC:
### chi-squared, or Fisher's exact if any expected cell < 5.
### Population: patients with an open-access RNA-seq file and a clade call
### of A7 or A9 (the population the clade-specific networks will use).
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(vroom) })

meta <- vroom("~/CESC_Network/02_cohort_data/validation_uganda/1_Output/1_3_HTMCP_clinical_HPV_HIV_with_rnaseq_flag.tsv",
              show_col_types = FALSE) %>%
  filter(has_rnaseq, HPV.clade %in% c("A7", "A9")) %>%
  transmute(HPV_clade = HPV.clade,
           Stage_group = gsub("Stage ([IV]+).*", "Stage \\1", Stage),
           Stage_group = ifelse(grepl("^Stage [IV]+$", Stage_group), Stage_group, NA_character_),
           histological_group = Final.histology,
           HIV_status = HIV.status)

test_feature <- function(var) {
  d <- meta %>% filter(!is.na(.data[[var]]))
  tab <- table(d$HPV_clade, d[[var]])
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

result <- bind_rows(test_feature("Stage_group"), test_feature("histological_group"),
                    test_feature("HIV_status"))
vroom_write(result, "~/CESC_Network/02_cohort_data/validation_uganda/1_Output/2_1_fisher_clinical_vs_clade.tsv", delim = "\t")
print(as.data.frame(result), row.names = FALSE)
