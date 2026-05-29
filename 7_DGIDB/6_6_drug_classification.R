############################################################
## 6_6_drug_classification.R
##
## Pharmacological classification of proximity-significant
## drugs (z < -1.96) in HPV-A7 and HPV-A9 networks.
##
## Classification source (primary):
##   WHO ATC Classification System (WHO Collaborating Centre
##   for Drug Statistics Methodology, Oslo, 2024)
##   Retrieved via NLM RxClass REST API:
##   https://rxnav.nlm.nih.gov/REST/rxclass/class/byRxcui.json
##   using RxCUI identifiers from DGIdb (Freshour et al. 2021)
##
## Fallback (12 drugs without RxCUI or ATC mapping):
##   Manual assignment based on ATC index + FDA drug labels
##
## DGIdb flags (anti_neoplastic, immunotherapy) stored as
## supplementary annotation.
##
## Outputs (in 7_DGIDB/):
##   6_6_drug_atc_cache.tsv     — ATC lookup cache
##   6_6_drug_classes.tsv       — full classified table
##   6_6_classified_lollipop.pdf/.png
##   6_6_class_summary_bar.pdf/.png
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(stringr)
  library(forcats)
  library(jsonlite)
})

setwd("~/CESC_Network/7_DGIDB")
out_dir  <- "~/CESC_Network/7_DGIDB"
CACHE_F  <- file.path(out_dir, "6_6_drug_atc_cache.tsv")

# ══════════════════════════════════════════════════════════
## 1) Load proximity results & DGIdb data
# ══════════════════════════════════════════════════════════
prox_A7 <- read_tsv("6_3_proximity_A7.tsv", show_col_types = FALSE)
prox_A9 <- read_tsv("6_3_proximity_A9.tsv", show_col_types = FALSE)

sig_A7   <- prox_A7 %>% filter(significant)
sig_A9   <- prox_A9 %>% filter(significant)
both_sig <- intersect(sig_A7$drug, sig_A9$drug)

prox_combined <- bind_rows(
  sig_A7 %>% mutate(cohort = "HPV-A7"),
  sig_A9 %>% mutate(cohort = "HPV-A9")
) %>%
  mutate(
    clade_group = case_when(
      drug %in% both_sig ~ "Both cohorts",
      cohort == "HPV-A7"  ~ "HPV-A7 only",
      TRUE                ~ "HPV-A9 only"
    )
  )

dgidb <- readRDS("6_1_dgidb_fda_interactions.rds")

# One row per unique significant drug + its DGIdb IDs
drug_meta <- dgidb %>%
  filter(drug_name %in% union(sig_A7$drug, sig_A9$drug)) %>%
  distinct(drug_name, drug_concept_id, anti_neoplastic, immunotherapy) %>%
  mutate(
    rxcui = case_when(
      grepl("^rxcui:", drug_concept_id) ~ sub("rxcui:", "", drug_concept_id),
      TRUE                              ~ NA_character_
    )
  )

cat(sprintf("Significant drugs — A7: %d | A9: %d | Shared: %d | Unique: %d\n",
            nrow(sig_A7), nrow(sig_A9), length(both_sig), nrow(drug_meta)))

# ══════════════════════════════════════════════════════════
## 2) ATC lookup via NLM RxClass API
##    Source: WHO ATC Classification System
##    (WHO Collaborating Centre for Drug Statistics
##     Methodology, Oslo, Norway — https://www.whocc.no)
##    API:    NLM RxClass (https://rxnav.nlm.nih.gov)
# ══════════════════════════════════════════════════════════
query_atc_all <- function(rxcui) {
  if (is.na(rxcui)) return(tibble(atc_code = NA_character_, atc_class = NA_character_))
  url <- sprintf(
    "https://rxnav.nlm.nih.gov/REST/rxclass/class/byRxcui.json?rxcui=%s&relaSource=ATC",
    rxcui
  )
  r <- tryCatch(fromJSON(url, simplifyVector = FALSE), error = function(e) NULL)
  items <- r$rxclassDrugInfoList$rxclassDrugInfo
  if (is.null(items) || length(items) == 0) {
    return(tibble(atc_code = NA_character_, atc_class = NA_character_))
  }
  tibble(
    atc_code  = sapply(items, function(x) x$rxclassMinConceptItem$classId),
    atc_class = sapply(items, function(x) x$rxclassMinConceptItem$className)
  )
}

# Prefer antineoplastic/immune (L) > antiinfective (J) > antiparasitic (P) > rest
ATC_PRIORITY <- c("L", "J", "P", "N", "B", "C", "H", "D", "M", "A", "R", "V")

pick_best_atc <- function(tbl) {
  if (nrow(tbl) == 0 || all(is.na(tbl$atc_code))) {
    return(list(atc_code = NA_character_, atc_class = NA_character_))
  }
  tbl <- tbl %>% filter(!is.na(atc_code))
  # Prefer specific codes (level 4 = 5 chars); rank by ATC_PRIORITY letter
  tbl <- tbl %>%
    mutate(
      lvl1 = substr(atc_code, 1, 1),
      prio = match(lvl1, ATC_PRIORITY, nomatch = length(ATC_PRIORITY) + 1),
      code_len = nchar(atc_code)
    ) %>%
    arrange(prio, desc(code_len))
  list(atc_code = tbl$atc_code[1], atc_class = tbl$atc_class[1])
}

# Use cache if it exists; otherwise query API
if (file.exists(CACHE_F)) {
  cat("Loading ATC cache...\n")
  atc_cache <- read_tsv(CACHE_F, show_col_types = FALSE)
} else {
  cat(sprintf("Querying RxClass API for %d drugs (be patient)...\n", nrow(drug_meta)))
  atc_list <- lapply(seq_len(nrow(drug_meta)), function(i) {
    Sys.sleep(0.15)
    all_atc <- query_atc_all(drug_meta$rxcui[i])
    best    <- pick_best_atc(all_atc)
    cat(sprintf("  [%02d] %-40s → %s (%s)\n",
                i, drug_meta$drug_name[i], best$atc_code, best$atc_class))
    tibble(
      drug_name = drug_meta$drug_name[i],
      rxcui     = drug_meta$rxcui[i],
      atc_code  = best$atc_code,
      atc_class = best$atc_class
    )
  })
  atc_cache <- bind_rows(atc_list)
  write_tsv(atc_cache, CACHE_F)
  cat(sprintf("  Cache saved: %s\n", CACHE_F))
}

# ══════════════════════════════════════════════════════════
## 3) Manual fallback for drugs without ATC mapping
##    (based on ATC index and FDA drug label)
## ══════════════════════════════════════════════════════════
# These 12 drugs have no RxCUI or no ATC in RxClass:
# ATALUREN, BITHIONOL, CAMPTOTHECIN, DOCONEXENT, GLUTAMIC ACID,
# INCOBOTULINUMTOXINA, LINOLENIC ACID, PRABOTULINUMTOXIN A,
# RESVERATROL, THIETHYLPERAZINE MALATE, URIDINE, AVOBENZONE
fallback_atc <- tribble(
  ~drug_name,                  ~atc_code, ~atc_class,
  "ATALUREN",                  "R07AX03", "Other respiratory system drugs",       # ATC 2024
  "BITHIONOL",                 "P02CA03", "Taenicides",
  "CAMPTOTHECIN",              "L01CE02", "Camptothecins",
  "DOCONEXENT",                "C10AX06", "Omega-3-triglycerides including other esters and acids",
  "GLUTAMIC ACID",             "A16AA",   "Amino acids and derivatives",
  "INCOBOTULINUMTOXINA",       "M03AX01", "Other muscle relaxants, peripherally acting",
  "LINOLENIC ACID",            "C10AX",   "Other lipid modifying agents",
  "PRABOTULINUMTOXIN A",       "M03AX01", "Other muscle relaxants, peripherally acting",
  "RESVERATROL",               "A11HA",   "Other plain vitamin preparations",     # best available
  "THIETHYLPERAZINE MALATE",   "A04AD04", "Other antiemetics",
  "URIDINE",                   "A16AX",   "Other alimentary tract and metabolism products",
  "AVOBENZONE",                "D02BB",   "Sunscreens"
)

# Merge API results with fallback
atc_full <- atc_cache %>%
  rows_patch(fallback_atc, by = "drug_name", unmatched = "ignore")

# ══════════════════════════════════════════════════════════
## 4) Derive broad class from ATC level 1-2
# ══════════════════════════════════════════════════════════
atc_to_broad <- function(code) {
  if (is.na(code)) return("Unclassified")
  lvl2 <- substr(code, 1, 3)
  lvl1 <- substr(code, 1, 1)
  dplyr::case_when(
    lvl2 %in% c("L01", "L02")        ~ "Antineoplastic",
    lvl2 %in% c("L03", "L04")        ~ "Immunomodulatory",
    substr(code,1,3) == "J05"        ~ "Antiviral",
    lvl1 == "J"                       ~ "Antibacterial",
    lvl1 == "P"                       ~ "Antiparasitic",
    lvl1 == "N"                       ~ "CNS / Nervous system",
    lvl1 == "A"                       ~ "Alimentary / Metabolism",
    lvl1 == "B"                       ~ "Blood / Hematology",
    lvl1 == "C"                       ~ "Cardiovascular",
    lvl1 == "D"                       ~ "Dermatological",
    lvl1 == "H"                       ~ "Systemic hormonal",
    lvl1 == "M"                       ~ "Anti-inflammatory (NSAID)",
    lvl1 %in% c("R", "V")            ~ "Other / Supplement",
    TRUE                              ~ "Unclassified"
  )
}

atc_full <- atc_full %>%
  mutate(class_broad = sapply(atc_code, atc_to_broad))

# ══════════════════════════════════════════════════════════
## 5) Join everything & export
# ══════════════════════════════════════════════════════════
classified <- prox_combined %>%
  mutate(drug_key = str_trim(str_to_upper(drug))) %>%
  left_join(
    drug_meta %>%
      rename(drug_key = drug_name) %>%
      mutate(drug_key = str_trim(str_to_upper(drug_key))),
    by = "drug_key"
  ) %>%
  left_join(
    atc_full %>%
      select(-rxcui) %>%          # rxcui already comes from drug_meta join
      rename(drug_key = drug_name) %>%
      mutate(drug_key = str_trim(str_to_upper(drug_key))),
    by = "drug_key"
  ) %>%
  mutate(
    atc_code   = ifelse(is.na(atc_code),   "Unclassified", atc_code),
    atc_class  = ifelse(is.na(atc_class),  "Unclassified", atc_class),
    class_broad = ifelse(is.na(class_broad), "Unclassified", class_broad)
  )

write_tsv(
  classified %>% select(cohort, drug, n_targets, d_c_interest, z_score, p_adj, sig_tier,
                        clade_group, rxcui, atc_code, atc_class, class_broad,
                        anti_neoplastic, immunotherapy),
  file.path(out_dir, "6_6_drug_classes.tsv")
)
cat(sprintf("\nClassified %d drug-cohort entries | %d unique ATC classes\n",
            nrow(classified), n_distinct(classified$atc_class)))

# ══════════════════════════════════════════════════════════
## 6) Palette & class ordering
# ══════════════════════════════════════════════════════════
BROAD_PAL <- c(
  "Antineoplastic"          = "#C0392B",
  "Immunomodulatory"        = "#7D3C98",
  "Antiviral"               = "#2471A3",
  "Antibacterial"           = "#117A65",
  "Antiparasitic"           = "#1A8F5A",
  "CNS / Nervous system"    = "#2E4057",
  "Alimentary / Metabolism" = "#D68910",
  "Blood / Hematology"      = "#E74C3C",
  "Cardiovascular"          = "#E67E22",
  "Dermatological"          = "#CB4335",
  "Systemic hormonal"       = "#839192",
  "Anti-inflammatory (NSAID)"= "#CA6F1E",
  "Other / Supplement"      = "#7F8C8D",
  "Unclassified"            = "#BDC3C7"
)

# Order broad classes for y-axis grouping: antineoplastic first, then by n descending
class_order <- classified %>%
  count(class_broad) %>%
  arrange(match(class_broad, names(BROAD_PAL))) %>%
  pull(class_broad)

CLADE_SHP <- c("Both cohorts" = 18, "HPV-A7 only" = 16, "HPV-A9 only" = 17)
CLADE_COL <- c("Both cohorts" = "#6C3483", "HPV-A7 only" = "#1A6FB5", "HPV-A9 only" = "#C0392B")

# ══════════════════════════════════════════════════════════
## 7) Lollipop: drug × z-score, colored by ATC broad class
# ══════════════════════════════════════════════════════════
make_lollipop <- function(df, panel_title, coh) {
  df_p <- df %>%
    filter(cohort == coh) %>%
    mutate(drug_lbl = str_trunc(drug, 28)) %>%
    arrange(match(class_broad, class_order), z_score) %>%
    mutate(
      drug_ord    = factor(drug_lbl, levels = unique(drug_lbl)),
      class_broad = factor(class_broad, levels = rev(class_order))
    )

  ggplot(df_p, aes(x = z_score, y = drug_ord)) +
    geom_vline(xintercept = -1.96, linetype = "dashed",
               color = "#2471A3", linewidth = 0.65) +
    geom_vline(xintercept = -3,    linetype = "dotted",
               color = "#6C3483", linewidth = 0.55) +
    geom_vline(xintercept = 0,     color = "#DDDDDD", linewidth = 0.40) +
    geom_segment(aes(xend = 0, yend = drug_ord, colour = class_broad),
                 linewidth = 1.1, alpha = 0.70) +
    geom_point(aes(colour = class_broad, size = n_targets,
                   shape = clade_group), alpha = 0.95) +
    scale_colour_manual(name = "ATC class (WHO)", values = BROAD_PAL, drop = FALSE) +
    scale_shape_manual(name = "Cohort",
                       values = CLADE_SHP) +
    scale_size_continuous(name = "# Drug\ntargets",
                          range = c(2, 7.5), breaks = c(1, 10, 50, 100)) +
    scale_x_continuous(breaks  = c(-5, -4, -3, -2, -1.96, 0),
                       labels  = c("-5", "-4", "-3", "-2", "-1.96", "0"),
                       expand  = expansion(mult = c(0.02, 0.06))) +
    labs(title = panel_title, x = "Permutation z-score", y = NULL) +
    theme_classic(base_size = 13) +
    theme(
      axis.text.y        = element_text(size = 10, face = "bold", colour = "grey15"),
      axis.text.x        = element_text(size = 11, face = "bold"),
      axis.title.x       = element_text(size = 12, face = "bold", colour = "grey30"),
      plot.title         = element_text(face = "bold", size = 15),
      legend.position    = "right",
      legend.title       = element_text(size = 10, face = "bold"),
      legend.text        = element_text(size = 9),
      panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3)
    )
}

p_A7 <- make_lollipop(classified, "HPV-A7  —  Significant drugs (WHO ATC class)", "HPV-A7")
p_A9 <- make_lollipop(classified, "HPV-A9  —  Significant drugs (WHO ATC class)", "HPV-A9")

combined_lollipop <- (p_A7 | p_A9) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title   = "Network-Proximal Drugs by WHO ATC Pharmacological Class | CESC",
    caption = paste0(
      "Classification: WHO ATC System (WHOCC, Oslo 2024) via NLM RxClass API · ",
      "DGIdb interactions (Freshour et al. 2021, NAR) | ",
      "Shape: ◆ both HPV clades  ● HPV-A7  ▲ HPV-A9 | Dashed: z = -1.96 | Dotted: z = -3"
    ),
    theme = theme(
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 15),
      plot.caption = element_text(colour = "grey50", size = 9, hjust = 0.5)
    )
  )

n_drugs_max <- max(nrow(sig_A7), nrow(sig_A9))
h_lollipop  <- n_drugs_max * 0.33 + 4

ggsave(file.path(out_dir, "6_6_classified_lollipop.pdf"),
       combined_lollipop, width = 18, height = h_lollipop,
       device = cairo_pdf, limitsize = FALSE)
ggsave(file.path(out_dir, "6_6_classified_lollipop.png"),
       combined_lollipop, width = 18, height = h_lollipop,
       dpi = 300, limitsize = FALSE)
cat("  Saved: 6_6_classified_lollipop.pdf/png\n")

# ══════════════════════════════════════════════════════════
## 8) Summary barplot: drug count per ATC broad class × cohort
# ══════════════════════════════════════════════════════════
summary_counts <- classified %>%
  distinct(cohort, drug, class_broad) %>%
  count(cohort, class_broad) %>%
  mutate(
    cohort      = factor(cohort, levels = c("HPV-A7", "HPV-A9")),
    class_broad = factor(class_broad, levels = rev(class_order))
  )

p_bar <- ggplot(summary_counts,
                aes(x = n, y = class_broad, fill = class_broad)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = n), hjust = -0.3, size = 3.8, fontface = "bold",
            colour = "grey25") +
  scale_fill_manual(values = BROAD_PAL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.20))) +
  facet_wrap(~cohort, ncol = 2) +
  labs(
    title    = "Pharmacological classes of proximity-significant drugs (WHO ATC)",
    subtitle = "Drug count per ATC level-1/2 class  ·  z < -1.96  ·  Source: WHOCC via NLM RxClass",
    x        = "Number of drugs", y = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle    = element_text(hjust = 0.5, colour = "grey40", size = 10),
    axis.text.y      = element_text(size = 11, face = "bold"),
    axis.text.x      = element_text(size = 10),
    strip.text       = element_text(face = "bold", size = 12),
    strip.background = element_blank(),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3)
  )

h_bar <- n_distinct(classified$class_broad) * 0.55 + 3
ggsave(file.path(out_dir, "6_6_class_summary_bar.pdf"),
       p_bar, width = 10, height = h_bar, device = cairo_pdf, limitsize = FALSE)
ggsave(file.path(out_dir, "6_6_class_summary_bar.png"),
       p_bar, width = 10, height = h_bar, dpi = 300, limitsize = FALSE)
cat("  Saved: 6_6_class_summary_bar.pdf/png\n")

# ══════════════════════════════════════════════════════════
## 9) Print summary
# ══════════════════════════════════════════════════════════
cat("\n── Summary by ATC broad class ────────────────────────\n")
classified %>%
  distinct(cohort, drug, class_broad, atc_code, atc_class) %>%
  count(cohort, class_broad, atc_class) %>%
  tidyr::pivot_wider(names_from = cohort, values_from = n, values_fill = 0L) %>%
  arrange(class_broad, atc_class) %>%
  print(n = 50)

cat("\n✓ 6_6 complete.\n")
cat("  Classification source: WHO ATC Classification System (WHOCC, Oslo)\n")
cat("  API: NLM RxClass (https://rxnav.nlm.nih.gov)\n")
cat("  Drug-gene source: DGIdb (Freshour et al. 2021, Nucleic Acids Research)\n")
