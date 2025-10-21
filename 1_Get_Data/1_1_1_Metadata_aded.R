library(ggplot2)
library(dplyr)
library(vroom)

factors_full <- vroom("~/CESC_Network/2_Prepro_TCGA_GTEx/2_4_Factors.tsv")

# Filtrar y limpiar los datos
factors <- factors_full %>%
  # Reemplazar NAs por "(NOT REPORTED)" antes de filtrar
  mutate(
    HPV_clade  = ifelse(is.na(HPV_clade), "(NOT REPORTED)", HPV_clade),
    figo_stage = ifelse(is.na(figo_stage) | trimws(figo_stage) == "", "(NOT REPORTED)", figo_stage)
  ) %>%
  filter(HPV_clade %in% c("A7_clade", "A9_clade", "(NOT REPORTED)"),
         figo_stage != "Solid Tissue Normal") %>%  # Eliminar SNT
  
  # Extraer el número romano principal del estadio usando regex
  mutate(Stage_group = gsub("Stage ([IV]+).*", "Stage \\1", figo_stage),
         Stage_group = ifelse(grepl("^Stage [IV]+$", Stage_group), 
                              Stage_group, 
                              "(NOT REPORTED)"))  # Conservar también los no coincidentes como "(NOT REPORTED)"


histological_type <- vroom("~/OCTAD_Cervical_Cancer/histological_type.tsv")

# 1) Asegura un mapeo único muestra -> histological_type
hist_map <- histological_type %>%
  mutate(sample = trimws(sample)) %>%
  select(sample, histological_type) %>%
  filter(!is.na(sample), sample != "") %>%
  distinct(sample, .keep_all = TRUE)

# 2) Limpia claves y hace LEFT JOIN (conserva exactamente las filas de factors)
factors_joined <- factors %>%
  mutate(specimenID = trimws(specimenID)) %>%
  left_join(hist_map, by = c("specimenID" = "sample"))



BiocGenerics::table(factors_joined$histological_type,factors_joined$HPV_clade)
#                                                  A7_clade A9_clade
#Adenosquamous                                        2        1
#Cervical Squamous Cell Carcinoma                    53      175
#Endocervical Adenocarcinoma of the Usual Type        1        4
#Endocervical Type of Adenocarcinoma                  4       13
#Endometrioid Adenocarcinoma of Endocervix            1        0
#Mucinous Adenocarcinoma of Endocervical Type         4        8


#glimpse(factors_joined)
#Rows: 252
#Columns: 13
#$ specimenID         <chr> "TCGA-C5-A1M5-01", "TCGA-DG-A2KK-01", "TCGA-C5-A1ML-01", "TCGA-R2-A69V-01…
#$ source             <chr> "TCGA", "TCGA", "TCGA", "TCGA", "TCGA", "TCGA", "TCGA", "TCGA", "TCGA", "…
#$ case_id_12         <chr> "TCGA-C5-A1M5", "TCGA-DG-A2KK", "TCGA-C5-A1ML", "TCGA-R2-A69V", "TCGA-C5-…
#$ sample_type        <chr> "Primary Tumor", "Primary Tumor", "Primary Tumor", "Primary Tumor", "Prim…
#$ HPV_type           <chr> "HPV33", "HPV16", "HPV16", "HPV39", "HPV16", "HPV18", "HPV16", "HPV16", "…
#$ HPV_clade          <chr> "A9_clade", "A9_clade", "A9_clade", "A7_clade", "A9_clade", "A7_clade", "…
#$ figo_stage         <chr> "Stage IB", "Stage IIIB", "Stage IB2", "Stage IB", "Stage IB", "Stage IB1…
#$ primary_diagnosis  <chr> "Squamous cell carcinoma, NOS", "Adenocarcinoma, NOS", "Squamous cell car…
#$ race               <chr> "white", "not reported", "white", "white", "black or african american", "…
#$ specimenID_portion <chr> "TCGA-C5-A1M5-01", "TCGA-DG-A2KK-01", "TCG_

BiocGenerics::table(factors_joined$Stage_group ,factors_joined$HPV_clade)

BiocGenerics::table(factors_joined$HPV_type,factors_joined$HPV_clade)

library(dplyr)
# Supongamos que tu data frame se llama factors_joined

factors_joined2 <- factors_joined %>%
  mutate(
    histological_group = case_when(
      histological_type == "Cervical Squamous Cell Carcinoma" ~ "Squamous",
      histological_type == "Adenosquamous" ~ "Adenosquamous",
      # Para los que sean adenocarcinomas, agrúpalos todos como "Adenocarcinoma"
      histological_type %in% c(
        "Endocervical Adenocarcinoma of the Usual Type",
        "Endocervical Type of Adenocarcinoma",
        "Mucinous Adenocarcinoma of Endocervical Type",
        "Endometrioid Adenocarcinoma of Endocervix"
      ) ~ "Adenocarcinoma",
      TRUE ~ NA_character_   # en caso de valores inesperados
    )
  )

# Ver tabla cruzada nueva
table(factors_joined2$histological_group, factors_joined2$HPV_clade)



#####################
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(gt)
  library(rlang)
  library(xml2)
})

# ---------- HELPERS (definir antes de usarlos) ----------
make_block <- function(df, var, feature_name, level_order = NULL,
                       test = c("auto","none")) {
  test <- match.arg(test)
  stopifnot(all(c("HPV_clade", var) %in% names(df)))
  
  var_sym <- rlang::sym(var)
  
  # Orden sugerido de niveles (si aplica)
  if (!is.null(level_order)) {
    df <- dplyr::mutate(df, !!var_sym := base::factor(!!var_sym, levels = level_order))
  }
  
  # El test NO incluye NAs
  df_nona <- df %>% dplyr::filter(!is.na(.data$HPV_clade), !is.na(!!var_sym))
  
  # Conteos y % por clado
  long <- df_nona %>%
    dplyr::count(HPV_clade, !!var_sym, name = "n") %>%
    dplyr::group_by(HPV_clade) %>%
    dplyr::mutate(pct = 100 * n / base::sum(n)) %>%
    dplyr::ungroup()
  
  # Determinar todos los niveles a mostrar
  all_levels <- if (!is.null(level_order)) level_order else
    long %>% dplyr::distinct(!!var_sym) %>% dplyr::pull(!!var_sym) %>% as.character()
  
  # Completar niveles/clados ausentes con 0
  long <- long %>%
    tidyr::complete(HPV_clade, !!var_sym := all_levels, fill = list(n = 0, pct = 0)) %>%
    dplyr::mutate(cell = sprintf("%d (%.1f%%)", n, pct)) %>%
    dplyr::select(-pct)
  
  # Ancho
  wide <- long %>%
    dplyr::select(!!var_sym, HPV_clade, cell, n) %>%
    tidyr::pivot_wider(
      names_from  = HPV_clade,
      values_from = c(cell, n),
      values_fill = list(cell = "0 (0.0%)", n = 0)
    )
  
  # Total por nivel
  n_cols <- grep("^n_", names(wide), value = TRUE)
  if (length(n_cols) == 0) {
    wide$Total <- 0L
  } else {
    wide <- dplyr::mutate(wide, Total = base::rowSums(dplyr::across(dplyr::all_of(n_cols))))
  }
  
  out <- wide %>%
    dplyr::transmute(
      Feature  = feature_name,
      Level    = !!var_sym,
      `A7_clade` = dplyr::coalesce(.data[["cell_A7_clade"]], "0 (0.0%)"),
      `A9_clade` = dplyr::coalesce(.data[["cell_A9_clade"]], "0 (0.0%)"),
      Total
    )
  
  # Prueba (si procede)
  if (test == "auto") {
    tab <- table(df_nona$HPV_clade, df_nona[[var]], useNA = "no")
    # Si solo hay una columna o fila no se puede testear: devuelve NA
    if (nrow(tab) < 2 || ncol(tab) < 2) {
      pval <- NA_real_; method <- "Not applicable"
    } else {
      chi <- suppressWarnings(stats::chisq.test(tab))
      use_fisher <- any(chi$expected < 5) || any(chi$expected == 0)
      test_obj <- if (use_fisher) stats::fisher.test(tab) else chi
      pval   <- signif(test_obj$p.value, 3)
      method <- if (use_fisher) "Fisher's exact" else "Chi-squared"
    }
  } else {
    pval <- NA_real_; method <- "Descriptive only (no test)"
  }
  
  attr(out, "p.value") <- pval
  attr(out, "method")  <- method
  out
}

save_table_xml <- function(df, file = "Table_clade_features.xml",
                           root_tag = "table", row_tag = "row", cell_tag = "cell") {
  df_chr <- df
  for (cn in names(df_chr)) df_chr[[cn]] <- as.character(df_chr[[cn]])
  doc <- xml2::xml_new_root(root_tag)
  header <- xml2::xml_add_child(doc, "header")
  for (cn in names(df_chr)) {
    col_node <- xml2::xml_add_child(header, "colname")
    xml2::xml_set_text(col_node, cn)
  }
  for (i in seq_len(nrow(df_chr))) {
    rnode <- xml2::xml_add_child(doc, row_tag)
    for (cn in names(df_chr)) {
      cnode <- xml2::xml_add_child(rnode, cell_tag, cn)
      xml2::xml_set_text(cnode, df_chr[[cn]][i])
    }
  }
  xml2::write_xml(doc, file); invisible(file)
}

# ---------- 1) Crear "histological_group" a partir de tus etiquetas ----------
factors <- factors_joined %>%
  mutate(
    histological_group = case_when(
      histological_type == "Cervical Squamous Cell Carcinoma" ~ "Squamous",
      histological_type == "Adenosquamous" ~ "Adenosquamous",
      histological_type %in% c(
        "Endocervical Adenocarcinoma of the Usual Type",
        "Endocervical Type of Adenocarcinoma",
        "Mucinous Adenocarcinoma of Endocervical Type",
        "Endometrioid Adenocarcinoma of Endocervix"
      ) ~ "Adenocarcinoma",
      TRUE ~ NA_character_
    ),
    # Ordenes sugeridos
    Stage_group        = factor(Stage_group, levels = c("Stage I","Stage II","Stage III","Stage IV")),
    histological_group = factor(histological_group, levels = c("Squamous","Adenocarcinoma","Adenosquamous"))
  )

# ---------- 2) Construir bloques y correr tests (solo FIGO y Histology) ----------
blk_stage <- make_block(
  df          = factors,
  var         = "Stage_group",
  feature_name= "FIGO stage",
  level_order = c("Stage I","Stage II","Stage III","Stage IV"),
  test        = "auto"
)

blk_hist <- make_block(
  df          = factors,
  var         = "histological_group",
  feature_name= "Histological group",
  level_order = c("Squamous","Adenocarcinoma","Adenosquamous"),
  test        = "auto"
)

# ---------- 3) Combinar, anotar p-valores y mostrar ----------
clin_tbl <- dplyr::bind_rows(blk_stage, blk_hist) %>%
  mutate(
    Level   = as.character(Level),
    Feature = factor(Feature, levels = c("FIGO stage","Histological group"))
  )

pv <- tibble::tibble(
  Feature = c("FIGO stage","Histological group"),
  method  = c(attr(blk_stage,"method"), attr(blk_hist,"method")),
  p_value = c(attr(blk_stage,"p.value"), attr(blk_hist,"p.value"))
)

note_txt <- sprintf(
  "**A7 vs A9** — FIGO stage: *%s*, p = %s; Histological group: *%s*, p = %s.",
  pv$method[pv$Feature=="FIGO stage"], pv$p_value[pv$Feature=="FIGO stage"],
  pv$method[pv$Feature=="Histological group"], pv$p_value[pv$Feature=="Histological group"]
)

gt_clin <- clin_tbl |>
  gt::gt(rowname_col = "Level", groupname_col = "Feature") |>
  gt::cols_label(
    `A7_clade` = gt::md("A7 clade"),
    `A9_clade` = gt::md("A9 clade"),
    Total      = gt::md("Total (n)")
  ) |>
  gt::fmt_number(columns = "Total", decimals = 0) |>
  gt::tab_options(table.font.size = gt::px(12),
                  data_row.padding = gt::px(4),
                  row_group.font.weight = "bold") |>
  gt::tab_source_note(source_note = gt::md(note_txt))

gt_clin

# ---------- 4) (Opcional) guardar XML ----------
# save_table_xml(clin_tbl, file = "Table_clade_stage_histology.xml")

# ---------- 5) (Opcional) imprimir p-valores en consola ----------
message("FIGO stage: ", pv$method[pv$Feature=="FIGO stage"], " p = ", pv$p_value[pv$Feature=="FIGO stage"])
message("Histological group: ", pv$method[pv$Feature=="Histological group"], " p = ", pv$p_value[pv$Feature=="Histological group"])


# Tamaño de efecto (Cramer's V) para las dos tablas
# FIGO
tab_figo <- table(factors$HPV_clade, factors$Stage_group, useNA="no")
cramer_figo <- DescTools::CramerV(tab_figo)   # requiere DescTools

# Histology
tab_hist <- table(factors$HPV_clade, factors$histological_group, useNA="no")
cramer_hist <- DescTools::CramerV(tab_hist)

cramer_figo; cramer_hist


