# ============================================================
# Script: Safe mapping of human-readable drug names onto OCTAD sRGES hits
# Goal  : Add 'display_name' ONLY when identity is unambiguous
#         (single InChIKey across all pert_id for a pert_iname).
#         Preserve original columns; add QC flags and metadata.
# Author: Joel
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(octad.db)
  library(signatureSearch)   # lincs_pert_info2 (pert_id, inchi_key, pref_name, aliases)
})

# -----------------------------
# 1) Load OCTAD sRGES results
# -----------------------------
res_A7 <- readRDS("~/CESC_Network/6_OCTAD/6_2_1_OCTAD_A7_results_0.2.rds")
res_A9 <- readRDS("~/CESC_Network/6_OCTAD/6_2_1_OCTAD_A9_results_0.2.rds")

# -----------------------------
# 2) LINCS metadata
#    - EH7270: relaciones pert_iname <-> pert_id (instancia-level)
#    - lincs_pert_info2: info por pert_id (InChIKey, pref_name, aliases)
# -----------------------------
eh <- octad.db::get_ExperimentHub_data("EH7270") %>%
  filter(pert_type == "trt_cp") %>%
  distinct(pert_iname, pert_id)           # NO colapsar aún

data("lincs_pert_info2", package = "signatureSearch")
info <- lincs_pert_info2 %>%
  distinct(pert_id, .keep_all = TRUE) %>%
  dplyr::select(pert_id, inchi_key, pref_name, compound_aliases)

# -----------------------------
# 3) Construir un RESUMEN 1:1 por pert_iname
#    Criterio de mapeo SEGURO:
#    - Todos los pert_id ligados a 'pert_iname' comparten UN solo InChIKey no-NA
#    => entonces podemos tomar un nombre preferido para ESE InChIKey.
# -----------------------------
collapse_unique <- function(x) {
  x <- unique(na.omit(x))
  if (length(x) == 0) NA_character_ else paste(x, collapse = "; ")
}

map_safe <- eh %>%
  left_join(info, by = "pert_id") %>%
  group_by(pert_iname) %>%
  summarise(
    n_pert_ids           = n_distinct(pert_id),
    pert_ids             = collapse_unique(pert_id),
    inchi_keys           = collapse_unique(inchi_key),
    n_distinct_inchi     = n_distinct(na.omit(inchi_key)),
    # candidatos de nombre (no se usan si hay ambigüedad)
    cand_pref_name       = collapse_unique(pref_name),
    cand_aliases         = collapse_unique(compound_aliases),
    .groups = "drop"
  ) %>%
  mutate(
    # condición estricta: un solo InChIKey y no NA
    identity_unambiguous = (!is.na(inchi_keys)) & (n_distinct_inchi == 1),
    # si identidad es inequívoca, proponer nombre preferido; si no, NA
    proposed_name        = if_else(identity_unambiguous & !is.na(cand_pref_name),
                                   cand_pref_name,
                                   if_else(identity_unambiguous & is.na(cand_pref_name) & !is.na(cand_aliases),
                                           cand_aliases, NA_character_))
  )

# -----------------------------
# 4) Función para NORMALIZAR texto (para detectar cambios triviales)
# -----------------------------
norm_name <- function(x){
  x <- tolower(trimws(x))
  if (requireNamespace("stringi", quietly = TRUE)) {
    x <- stringi::stri_trans_general(x, "Latin-ASCII")
  }
  x <- gsub("[[:punct:]]+", " ", x)
  x <- gsub("\\s+", " ", x)
  x
}

# -----------------------------
# 5) Aplicar mapeo SEGURO a una tabla de sRGES
# -----------------------------
apply_safe_map <- function(res_tbl){
  res_tbl %>%
    left_join(map_safe, by = "pert_iname") %>%
    mutate(
      # nombre a mostrar: sólo si identidad es inequívoca y hay nombre propuesto
      display_name = if_else(identity_unambiguous & !is.na(proposed_name),
                             proposed_name, pert_iname),
      changed      = display_name != pert_iname,
      trivial_fmt  = norm_name(display_name) == norm_name(pert_iname),
      real_change  = changed & !trivial_fmt
    )
}

resA7_final <- apply_safe_map(res_A7)
resA9_final <- apply_safe_map(res_A9)

# -----------------------------
# 6) QC mínimo (conteos)
# -----------------------------
qc_counts <- function(df){
  df %>%
    summarise(
      total = n(),
      unambiguous = sum(identity_unambiguous, na.rm = TRUE),
      renamed = sum(real_change, na.rm = TRUE),
      trivial_changes = sum(changed & !real_change, na.rm = TRUE),
      unchanged = sum(!changed, na.rm = TRUE)
    )
}
cat("A7 QC:\n"); print(qc_counts(resA7_final))
cat("A9 QC:\n"); print(qc_counts(resA9_final))

# (Opcional) listar los cambios REALES aplicados
A7_changes <- resA7_final %>%
  filter(real_change) %>%
  dplyr::select(pert_iname, display_name, pert_ids, inchi_keys, sRGES, n) %>%
  arrange(sRGES)
A9_changes <- resA9_final %>%
  filter(real_change) %>%
  dplyr::select(pert_iname, display_name, pert_ids, inchi_keys, sRGES, n) %>%
  arrange(sRGES)

# print(A7_changes, n = Inf); print(A9_changes, n = Inf)

# -----------------------------
# 7) Guardar resultados
# -----------------------------
saveRDS(resA7_final, "~/CESC_Network/6_OCTAD/6_2_1_1_OCTAD_A7_results_0.2_named_safe.rds")
saveRDS(resA9_final, "~/CESC_Network/6_OCTAD/6_2_1_1_OCTAD_A9_results_0.2_named_safe.rds")
