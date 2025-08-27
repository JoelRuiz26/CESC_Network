# 6_0 
#Get the input signatures cleaned for octad tool:

# =========================================
# Load Required Libraries
# =========================================
library(dplyr)
library(vroom)
library(tibble)

setwd("~/CESC_Network/6_OCTAD/")

# =========================================
# Load OCTAD Preprocessed Data
# =========================================
metadata_full <- vroom("~/CESC_Network/2_Prepro_TCGA_GTEx/2_4_Factors.tsv")

#Clean and filter metadata (remove normal tissue samples)
metadata_full$sample_type <- trimws(metadata_full$sample_type)
metadata <- metadata_full %>%
  filter(sample_type != "Solid Tissue Normal")

# Also update metadata with "-01" appended
metadata <- metadata %>%
  mutate(
    cases.submitter_id = ifelse(
      grepl("-01$", cases.submitter_id),
      cases.submitter_id,
      paste0(cases.submitter_id, "-01")
    )
  )


dim(metadata)   #[1] 268   8

# =========================================
# Sample selection
# =========================================
muestras_clado_A7 <- metadata %>%
  filter(HPV_clade == "A7_clade") %>%
  pull(cases.submitter_id)

muestras_clado_A9 <- metadata %>%
  filter(HPV_clade == "A9_clade") %>%
  pull(cases.submitter_id)


# =========================================
# Save samples
# =========================================

saveRDS(muestras_clado_A7,"~/CESC_Network/6_OCTAD/6_0_A7_Samples.rds")
saveRDS(muestras_clado_A9,"~/CESC_Network/6_OCTAD/6_0_A9_Samples.rds")



