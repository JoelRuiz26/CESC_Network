#Script for deepmap exploration
#https://pmc.ncbi.nlm.nih.gov/articles/PMC5709193/
#Meyers RM, Bryan JG, McFarland JM, Weir BA, Sizemore AE, Xu H, Dharia NV, Montgomery PG, Cowley GS, Pantel S, Goodale A, Lee Y, Ali LD, Jiang G, 
#Lubonja R, Harrington WF, Strickland M, Wu T, Hawes DC, Zhivich VA, Wyatt MR, Kalani Z, Chang JJ, Okamoto M, Stegmaier K, Golub TR, Boehm JS, Vazquez F, Root DE, Hahn WC, Tsherniak A. 
#Computational correction of copy number effect improves specificity of CRISPR-Cas9 essentiality screens in cancer cells. 
#Nat Genet. 2017 Dec;49(12):1779-1784. doi: 10.1038/ng.3984. Epub 2017 Oct 30. PMID: 29083409; PMCID: PMC5709193.


library(depmap)
library(dplyr)

#Get metadata
cell_metadata <- depmap_metadata()
glimpse(cell_metadata)
#Rows: 1,840
#Columns: 29
#$ depmap_id                  <chr> "ACH-000016", "ACH-000032", "ACH-000033", "ACH-000043", …
#$ cell_line_name             <chr> "SLR 21", "MHH-CALL-3", "NCI-H1819", "Hs 895.T", "HEK TE…
#$ stripped_cell_line_name    <chr> "SLR21", "MHHCALL3", "NCIH1819", "HS895T", "HEKTE", "TE6…
#$ cell_line                  <chr> "SLR21_KIDNEY", "MHHCALL3_HAEMATOPOIETIC_AND_LYMPHOID_TI…
#$ aliases                    <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "MUTZ3",…
#$ cosmic_id                  <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, …
#$ sex                        <chr> NA, "Female", "Female", "Female", NA, "Female", "Male", …
#$ source                     <chr> "Academic lab", "DSMZ", "Academic lab", "ATCC", "Academi…
#$ RRID                       <chr> "CVCL_V607", "CVCL_0089", "CVCL_1497", "CVCL_0993", "CVC…
#$ WTSI_master_cell_ID        <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, …
#$ sample_collection_site     <chr> "kidney", "bone_marrow", "lymph_node", "fibroblast", "ki…
#$ primary_or_metastasis      <chr> "Metastasis", NA, "Metastasis", "Metastasis", NA, "Prima…
#$ primary_disease            <chr> "Kidney Cancer", "Leukemia", "Lung Cancer", "Non-Cancero…
#$ subtype_disease            <chr> "Renal Cell Carcinoma", "Acute Lymphoblastic Leukemia (A…
#$ age                        <chr> NA, "11", "55", "48", NA, "1", NA, "61", "16", "35", "50…
#$ sanger_id                  <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, …
#$ additional_info            <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, …
#$ lineage                    <chr> "kidney", "blood", "lung", "fibroblast", "kidney", "soft…
#$ lineage_subtype            <chr> "renal_cell_carcinoma", "ALL", "NSCLC", "fibroblast_skin…
#$ lineage_sub_subtype        <chr> NA, "b_cell", "NSCLC_adenocarcinoma", NA, NA, NA, NA, "b…
#$ lineage_molecular_subtype  <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, …
#$ default_growth_pattern     <chr> NA, NA, NA, "2D: adherent", NA, "2D: adherent", NA, NA, …
#$ model_manipulation         <chr> NA, NA, NA, NA, "immortalized", NA, "immortalized", NA, …
#$ model_manipulation_details <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, …
##$ patient_id                 <chr> "PT-JnARLB", "PT-p2KOyI", "PT-9p1WQv", "PT-rTUVZQ", "PT-…
#$ parent_depmap_id           <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, …
#$ Cellosaurus_NCIt_disease   <chr> "Clear cell renal cell carcinoma", "Childhood B acute ly…
#$ Cellosaurus_NCIt_id        <chr> "C4033", "C9140", "C3512", "C3224", NA, "C3359", NA, "C4…
#$ Cellosaurus_issues         <chr> NA, NA, NA, NA, "No information is available about this …


length(unique(cell_metadata$cell_line_name))
#[1] 1748

############################################################################################
#Drugs per cell line
drug_data <- depmap_drug_sensitivity() 
#glimpse(drug_data) #[1] 2708508      14
#Rows: 2,708,508
#Columns: 14
#$ depmap_id    <chr> "ACH-000001", "ACH-000007", "ACH-000008", "ACH-000010", "ACH-000011", …
#$ cell_line    <chr> "NIHOVCAR3_OVARY", "LS513_LARGE_INTESTINE", "A101D_SKIN", "NCIH2077_LU…
#$ compound     <chr> "BRD-A00077618-236-07-6::2.5::HTS", "BRD-A00077618-236-07-6::2.5::HTS"…
#$ dependency   <dbl> -0.01557664, -0.09573033, 0.37948042, 0.11889032, 0.14534607, 0.103348…
#$ broad_id     <chr> "BRD-A00077618-236-07-6", "BRD-A00077618-236-07-6", "BRD-A00077618-236…
#$ name         <chr> "8-bromo-cGMP", "8-bromo-cGMP", "8-bromo-cGMP", "8-bromo-cGMP", "8-bro…
#$ dose         <dbl> 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, …
#$ screen_id    <chr> "HTS", "HTS", "HTS", "HTS", "HTS", "HTS", "HTS", "HTS", "HTS", "HTS", …
#$ moa          <chr> "PKA activator", "PKA activator", "PKA activator", "PKA activator", "P…
#$ target       <chr> "PRKG1", "PRKG1", "PRKG1", "PRKG1", "PRKG1", "PRKG1", "PRKG1", "PRKG1"…
#$ disease_area <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
#$ indication   <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
#$ smiles       <chr> "Nc1nc(O)c2nc(Br)n([C@@H]3O[C@@H]4COP(O)(=O)O[C@H]4[C@H]3O)c2n1", "Nc1…
#$ phase        <chr> "Preclinical", "Preclinical", "Preclinical", "Preclinical", "Preclinic…

length(unique(drug_data$cell_line))
#[1] 578

ncol(BiocGenerics::table(drug_data$cell_line, drug_data$compound))/
  nrow(BiocGenerics::table(drug_data$cell_line, drug_data$compound)) 
#[1] 8.107266 compounds per cell



############################################################################################
#Transcriptomic data (RNA-seq) 
rna_data <- depmap::depmap_TPM() %>% 
  filter(depmap_id %in% cell_metadata$depmap_id)

#Cell lines RNA
length(unique(rna_data$cell_line))  #[1] 1405
length(unique(rna_data$gene))
#[1] 19221 genes by cell line


############################################################################################
#In a CRISPR-Cas9 dropout screen:
#  A pooled library of sgRNAs targets nearly all genes.
#At baseline (t₀), guides are present at similar abundances.
#After several cell doublings (tfinal), sgRNAs targeting essential genes drop out because those cells grow slower or die.
#sgRNA abundance is measured by sequencing, and the depletion is summarized as log fold change (LFC).

# CERES model formulation for CRISPR dependency scores
#
# For each sgRNA g in cell line c:
#
#   LFC[g,c] = D[gene,c] * E[g] + C[gene,c] + ε
#
# where:
#   - LFC[g,c] : observed log fold change of sgRNA g in cell line c
#   - D[gene,c]: dependency score of the targeted gene in that cell line (the main parameter to estimate)
#   - E[g]     : estimated efficiency of sgRNA g (learned from the model)
#   - C[gene,c]: copy-number effect correction (adjusts for DNA amplifications leading to extra breaks)
#   - ε        : residual error (noise not explained by the model)
#
# The CERES algorithm fits this equation jointly across all sgRNAs, genes, 
# and cell lines, solving a regression problem that simultaneously estimates:
#   1. Gene dependency scores (D)
#   2. sgRNA efficiencies (E)
#   3. Copy-number corrections (C)
#
# After estimation, dependency scores are normalized such that:
#   - Pan-essential genes have an average score ≈ -1
#   - Non-essential genes have an average score ≈ 0
#
# Interpretation:
#   - 0   ≈ knockout has no effect
#   - -1  ≈ knockout is lethal, similar to core essential genes
#   - < -1 indicates even stronger essentiality

# CRISPR 
crispr_data <- depmap_crispr() %>% dplyr::filter(depmap_id %in% cell_metadata$depmap_id)
glimpse(crispr_data)
#Rows: 18,881,196
#Columns: 6
#$ depmap_id  <chr> "ACH-000001", "ACH-000004", "ACH-000005", "ACH-000007", "ACH-000009", "A…
#$ gene       <chr> "A1BG (1)", "A1BG (1)", "A1BG (1)", "A1BG (1)", "A1BG (1)", "A1BG (1)", …
#$ dependency <dbl> -0.1348083358, 0.0818526684, -0.0941960326, -0.0115440496, -0.0507823258…
#$ entrez_id  <int> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, …
#$ gene_name  <chr> "A1BG", "A1BG", "A1BG", "A1BG", "A1BG", "A1BG", "A1BG", "A1BG", "A1BG", …
#$ cell_line  <chr> "NIHOVCAR3_OVARY", "HEL_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "HEL9217_HA…

length(unique(crispr_data$gene_name))
#[1] 17386

#All genes are in all cells? 
crispr_data %>% 
  dplyr::group_by(cell_line) %>% 
  dplyr::summarise(n_genes = dplyr::n_distinct(gene_name)) %>% 
  arrange(desc(n_genes))

genes_por_cellline <- crispr_data %>%
  group_by(cell_line) %>%
  summarise(n_genes = n_distinct(gene_name), .groups = "drop")
head(genes_por_cellline)
# A tibble: 6 × 2
#cell_line                      n_genes
#<chr>                            <int>
#1 143B_BONE                        17386
#2 170MGBA_CENTRAL_NERVOUS_SYSTEM   17386
#3 21NT_BREAST                      17386
#4 22RV1_PROSTATE                   17386


summary(crispr_A7_A9$dependency)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-2.91687 -0.15701 -0.03954 -0.15429  0.04015  0.79146 


