## ─────────────────────────────────────────────────────────
## 0) Paquetes y carga del cache completo
## ─────────────────────────────────────────────────────────
# install.packages(c("ctrialsgov","dplyr","readr","stringr","tidyr","purrr"))
library(ctrialsgov)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)

# Carga el cache completo (~150 MB). Si ya lo cargaste en esta sesión, no hace falta repetir.
ctgov_load_cache()

## (Opcional) mira los términos categóricos válidos
terms <- ctgov_query_terms()
str(terms, max.level = 1)

#> str(terms, max.level = 1)
#List of 11
#$ study_type         : chr [1:4] "Interventional" "Observational" "Observational [Patient Registry]" "Expanded Access"
#$ allocation         : chr [1:3] "Randomized" "N/A" "Non-Randomized"
#$ intervention_model : chr [1:5] "Parallel Assignment" "Single Group Assignment" "Crossover Assignment" "Sequential Assignment" ...
#$ observational_model: chr [1:8] "Cohort" "Case-Control" "Case-Only" "Other" ...
#$ primary_purpose    : chr [1:10] "Treatment" "Prevention" "Basic Science" "Supportive Care" ...
#$ time_perspective   : chr [1:4] "Prospective" "Retrospective" "Cross-Sectional" "Other"
#$ masking_description: chr [1:5] "None (Open Label)" "Single" "Double" "Triple" ...
#$ sampling_method    : chr [1:3] "" "Non-Probability Sample" "Probability Sample"
#$ phase              : chr [1:8] "N/A" "Early Phase 1" "Phase 1" "Phase 1/Phase 2" ...
#$ gender             : chr [1:3] "All" "Female" "Male"
#$ sponsor_type       : chr [1:4] "Industry" "NIH" "U.S. Fed" "Other"



# 1) Trae el universo completo vía un filtro categórico “total”.
#    (Usamos todos los study_type válidos para abarcar todo el dataset.)
df_all <- ctgov_query(
  study_type = terms$study_type,   # => Interventional, Observational, etc. (todos)
  ignore_case = TRUE, match_all = FALSE
)

dim(df_all) #[1] 389454     32
# 2) Verifica columnas disponibles
print(names(df_all))
#[1] "nct_id"                         "start_date"                    
#[3] "phase"                          "enrollment"                    
#[5] "brief_title"                    "official_title"                
#[7] "primary_completion_date"        "study_type"                    
#[9] "rec_status"                     "completion_date"               
#[11] "last_update"                    "description"                   
#[13] "eudract_num"                    "other_id"                      
#[15] "allocation"                     "intervention_model"            
#[17] "observational_model"            "primary_purpose"               
#[19] "time_perspective"               "masking_description"           
#[21] "intervention_model_description" "sampling_method"               
#[23] "gender"                         "minimum_age"                   
#[25] "maximum_age"                    "population"                    
#[27] "criteria"                       "sponsor"                       
#[29] "sponsor_type"                   "conditions"                    
#[31] "interventions"                  "outcomes" 


library(dplyr)
library(stringr)
library(readr)

# --- Regex de inclusión (núcleo y opcional CIN) y exclusión (falsos positivos) ---
include_core_rx <- "(?i)\\b(cervic|uterine\\s+cervix|cervix\\s+uteri)\\b.*\\b(cancer|carcinoma|neoplasm|neoplasia|tumou?r|adenocarcinoma|squamous)\\b"
include_cin_rx  <- "(?i)\\bCIN\\b|cervical\\s+intraepithelial\\s+neoplas"  # opcional (precáncer)
exclude_rx      <- "(?i)cervical\\s+(spine|spondyl|radiculopathy|dystonia)|neck\\s+pain"

# --- Filtrado sobre 'conditions' (solo esa columna) ---
cc_from_conditions <- df_all %>%
  filter(!is.na(conditions) & nzchar(conditions)) %>%
  filter(
    str_detect(conditions, include_core_rx) |    # núcleo de cáncer cervicouterino
      str_detect(conditions, include_cin_rx)       # opcional: CIN (ajustable)
  ) %>%
  filter(!str_detect(conditions, exclude_rx))     # quitar no-ginecológico

head(cc_from_conditions)

# A tibble: 6 × 32
#nct_id      start_date phase enrollment brief_title official_title primary_completion_d…¹ study_type
#<chr>       <date>     <chr>      <int> <chr>       <chr>          <date>                 <chr>     
#1 NCT05022511 2021-09-30 N/A         6250 Three Bird… Three Birds W… 2022-09-30             Intervent…
#2 NCT04712851 2021-06-30 Phas…         25 Pembrolizu… A Phase II Op… 2023-01-31             Intervent…
#3 NCT04895020 2021-05-28 Phas…       1200 Immunobrid… Evaluating th… 2021-12-31             Intervent…
#4 NCT04910802 2021-05-03 Phas…     150000 Concomitan… Concomitant H… 2025-12-31             Intervent…
#5 NCT04831437 2021-04-30 Phas…         60 Clinical R… Comparison of… 2023-03-31             Intervent…
#6 NCT05045755 2021-04-13 NA          1339 The Durabi… The Durabilit… 2023-08-10             Observati…
# ℹ abbreviated name: ¹​primary_completion_date
# ℹ 24 more variables: rec_status <chr>, completion_date <date>, last_update <date>,
#   description <chr>, eudract_num <chr>, other_id <chr>, allocation <chr>, intervention_model <chr>,
#   observational_model <chr>, primary_purpose <chr>, time_perspective <chr>,
#   masking_description <chr>, intervention_model_description <chr>, sampling_method <chr>,
#   gender <chr>, minimum_age <dbl>, maximum_age <dbl>, population <chr>, criteria <chr>,
#   sponsor <chr>, sponsor_type <chr>, conditions <chr>, interventions <list<tibble[,3]>>, …


# Guarda el universo completo
saveRDS(df_all, "~/CESC_Network/9_ClinicalTrials/ctgov_all_studies.rds")

# Guarda solo los filtrados a cáncer cervicouterino
saveRDS(cc_from_conditions, "~/CESC_Network/9_ClinicalTrials/ctgov_cc_from_conditions.rds")



#############################

# --- 1) Inspecciona la estructura interna de 'interventions'
cc_from_conditions$interventions[[1]] %>% glimpse()  # mira nombres reales (p.ej. type, name, description)




library(dplyr)
library(tidyr)
library(stringr)
library(readr)

# 1) Desanidar y limpiar nombres
iv_long <- cc_from_conditions %>%
  dplyr::select(nct_id, brief_title, phase, study_type,rec_status,primary_purpose, interventions) %>%
  unnest(interventions) %>%  # columnas: intervention_type, name, description
  mutate(
    iv_type = intervention_type,
    iv_name_raw = name,
    # quita el prefijo "Tipo: " si viene incluido en 'name'
    iv_name = str_replace(iv_name_raw, "^[^:]+:\\s*", "") %>% str_squish()
  )

BiocGenerics::table(iv_long$iv_type)
#Behavioral          Biological Combination Product              Device     Diagnostic Test  Dietary Supplement 
#5                  55                   3                  44                   7                   1 
#Drug             Genetic               Other           Procedure           Radiation 
#74                   1                  35                  64                   4 

BiocGenerics::table(iv_long$primary_purpose)
#Basic Science       Diagnostic Health Services Research                    Other               Prevention 
#2                       18                       17                        1                       49 
#Screening          Supportive Care                Treatment 
#33                        2                      133 

iv_core <- iv_long %>%
  filter(!primary_purpose %in% c("Prevention", "Diagnostic"))

drug_types <- c("Drug", "Biological", "Combination Product", "Dietary Supplement", "Genetic")

iv_drugs <- iv_long %>%
  filter(!is.na(iv_type), iv_type %in% drug_types) %>%
  filter(!is.na(iv_name), nzchar(iv_name)) %>%
  group_by(nct_id) %>%
  summarise(drugs = paste(sort(unique(iv_name)), collapse = " | "),
            .groups = "drop")


unique(iv_drugs$drugs)
#[1] "motexafin lutetium | photodynamic therapy"                                                                               
#[2] "HspE7"                                                                                                                   
#[3] "Celecoxib"                                                                                                               
#[4] "EM-1421"                                                                                                                 
#[5] "Pemetrexed"                                                                                                              
#[6] "progesterone"                                                                                                            
#[7] "A007 | placebo"                                                                                                          
#[8] "Iopamidol 370 mgI/mL"                                                                                                    
#[9] "defined green tea catechin extract | placebo"                                                                            
#[10] "GSK Biologicals' HPV-16/18 Vaccine (Cervarix™) | Havrix™"                                                                
#[11] "Cervarix™"                                                                                                               
#[12] "lapatinib (GW572016) | pazopanib (GW786034)"                                                                             
#[13] "HspE7 and Poly-ICLC"                                                                                                     
#[14] "A-007"                                                                                                                   
#[15] "Hexaminolevulinate (HAL) | Placebo"                                                                                      
#[16] "Topical imiquimod therapy"                                                                                               
#[17] "DNA vaccination | imiquimod | intra-lesional vaccine administration | intramuscular vaccination"                         
#[18] "Placebo | RO5217790"                                                                                                     
#[19] "ADXS11-001 (Lm-LLO-E7) | Placebo Control"                                                                                
#[20] "cisplatin | paclitaxel | topotecan"                                                                                      
#[21] "HAL 0.2% with illumination | HAL 1% with illumination | HAL 5% with illumination | Placebo ointment without illumination"
#[22] "Topical Imiquimod"                                                                                                       
#[23] "Cisplatin | Nimotuzumab"                                                                                                 
#[49] "9-valent HPV vaccine | Imiquimod"                                                                                        
#[50] "CERVARIX | GARDASIL"                                                                                                     
#[51] "Vvax001 therapeutic cancer vaccine"                                                                                      
#[52] "SOR007 (Uncoated Nanoparticulate Paclitaxel) Ointment"                                                                   
#[53] "Imiquimod | Topical Fluorouracil"                                                                                        
#[54] "GX-188E, GX-I7 | GX-188E, Imiquimod"                                                                                     
#[55] "ABI-1968"                                                                                                                
#[56] "BLS-ILB-E710c | Placebo"                                                                                                 
#[57] "ICG | TMTP1"                                                                                                             
#[58] "FRD"                                                                                                                     
#[59] "HPV (Types 16, 18) Vaccine, Adsorbed"                                                                                    
#[60] "9-valent HPV Recombinant Vaccine | Placebo"                                                                              
#[61] "Topical ABI-1968"                                                                                                        
#[62] "Composite Gel Containing Black Raspberry | placebo"                                                                      
#[63] "Gardasil-9 | Placebo"                                                                                                    
#[64] "Vaccinated group"                                                                                                        
#[65] "Artesunate vaginal insert | Placebo vaginal insert"                                                                      
#[66] "pNGVL4aCRTE6E7L2"                                                                                                        
#[67] "E7 TCR"                                                                                                                  
#[68] "Cevira® | Placebo"                                                                                                       
#[69] "Bivalent HPV vaccine | Nonavalent HPV vaccine"                                                                           
#[70] "Recombinant Human Papillomavirus Bivalent (Types 16, 18) Vaccine (Escherichia coli)"                                     
#[71] "Pembrolizumab"                                                                                                           
#[72] "9-valent HPV vaccine"                                                                                                    
#[73] "Gardasil9"    

saveRDS(iv_long,   "ctgov_cc_interventions_long.rds")
saveRDS(iv_drugs,  "ctgov_cc_drugs_by_nct.rds")



A7_drugs <- readRDS("~/CESC_Network/6_OCTAD/6_3_RES/A7/RES_A7_common3_collapsed_FDA_Launched_0.20.rds")
A9_drugs <- readRDS("~/CESC_Network/6_OCTAD/6_3_RES/A9/RES_A9_common3_collapsed_FDA_Launched_0.20.rds")
#> head(A7_drugs)
# A tibble: 6 × 2
#pert_iname       sRGES
#<chr>            <dbl>
#1 pentamidine     -0.330
#2 maraviroc       -0.327
#3 raloxifene      -0.315
#4 triflupromazine -0.288
#5 perphenazine    -0.283
#6 fulvestrant     -0.275








library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(readr)

# Normalizador genérico (sin alias predefinidos)
normalize_drug <- function(x) {
  x %>%
    tolower() %>%
    str_replace_all("[®™©]", "") %>%
    str_replace_all("\\b\\d+(\\.\\d+)?\\s*(mg|mcg|µg|g|ml|mL|%)\\b", "") %>%
    str_replace_all("\\b(iv|po|sc|im|topical)\\b", "") %>%
    str_replace_all("(ointment|gel|cream|solution|suspension|capsule|tablet|injection|vaginal( insert)?)", "") %>%
    str_replace_all("\\([^)]*\\)", "") %>%   # quita paréntesis
    str_replace_all("\\bplacebo\\b", "") %>%
    str_replace_all("[-–_]", " ") %>%
    str_squish()
}

# (opcional) tabla de alias que TÚ puedes ir llenando; empieza vacía
alias_tbl <- tibble(from = character(), to = character())

apply_alias <- function(x, alias_tbl) {
  if (nrow(alias_tbl) == 0) return(x)
  out <- x
  for (i in seq_len(nrow(alias_tbl))) {
    out <- if_else(out == alias_tbl$from[i], alias_tbl$to[i], out)
  }
  out
}




# Usa tus tipos reales; si ya definiste 'drug_types', úsalo; si no, infiere de iv_long
# drug_types <- c("Drug","Biological","Combination Product","Dietary Supplement","Genetic")

iv_drugs_norm <- iv_long %>%
  filter(!is.na(iv_type), iv_type %in% unique(iv_long$iv_type)) %>%   # conserva tu set real
  filter(!is.na(iv_name), nzchar(iv_name)) %>%
  mutate(
    drug_norm = normalize_drug(iv_name),
    drug_norm = apply_alias(drug_norm, alias_tbl)
  ) %>%
  filter(nzchar(drug_norm)) %>%
  distinct(nct_id, drug_norm)




# Une A7 y A9 y normaliza pert_iname igual que arriba
lincs_all <- bind_rows(
  A7_drugs %>% mutate(source_set = "A7"),
  A9_drugs %>% mutate(source_set = "A9")
) %>%
  transmute(pert_iname, sRGES, source_set) %>%
  mutate(
    pert_norm = normalize_drug(pert_iname),
    pert_norm = apply_alias(pert_norm, alias_tbl)
  ) %>%
  filter(nzchar(pert_norm))

# Coincidencias exactas por nombre normalizado
matches <- iv_drugs_norm %>%
  inner_join(lincs_all, by = c("drug_norm" = "pert_norm")) %>%
  arrange(nct_id, drug_norm, desc(abs(sRGES)))

# Fármacos de tus ensayos que no aparecen en tu lista LINCS
no_match <- iv_drugs_norm %>%
  anti_join(lincs_all, by = c("drug_norm" = "pert_norm"))

# Resumen por estudio (agrega fase/estatus desde tu cc_from_conditions)
cc_match_summary <- matches %>%
  left_join(cc_from_conditions %>% select(nct_id, phase, rec_status, brief_title), by = "nct_id") %>%
  group_by(nct_id, brief_title, phase, rec_status) %>%
  summarise(
    lincs_hits = paste(unique(drug_norm), collapse = " | "),
    sets = paste(unique(source_set), collapse = " | "),
    top_min_sRGES = suppressWarnings(min(sRGES, na.rm = TRUE)),
    .groups = "drop"
  )




# Guarda resultados
write_csv(matches,         "cc_lincs_matches_by_nct.csv")
write_csv(no_match,        "cc_drugs_without_lincs_match.csv")
write_csv(cc_match_summary,"cc_lincs_matches_by_nct_summary.csv")

# Vistas rápidas
head(matches, 20)
head(no_match, 20)
head(cc_match_summary, 20)


