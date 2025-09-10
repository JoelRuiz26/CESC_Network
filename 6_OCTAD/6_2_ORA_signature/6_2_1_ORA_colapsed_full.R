# ================================
# ORA for DGE
# - Input: consensus_tbl (clade, regulation, gene), universe
# - Output:
#     * FULL (all ora_terms):        6_2_2_Ora_tbl.rds
#     * COLLAPSED(Grouped ORA_terms) 6_2_2_Ora_tbl_0.95.rds
# ================================

suppressPackageStartupMessages({
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(tibble)   
})

# ---- Inputs ----
consensus_tbl <- readRDS("~/CESC_Network/6_OCTAD/6_1_DGE_signature/6_1_1_Output_plots/6_10_Consensus_genes_3methods.rds")
universe      <- readRDS("~/CESC_Network/6_OCTAD/6_1_DGE_signature/6_1_0_Output_rds/6_1_6_Universe_ORA.rds")

# ---- Settings ----
ontologies       <- c("BP","CC","MF")
min_gene_list_n  <- 10
p_cut            <- 0.05
q_cut            <- 0.20
key_type         <- "SYMBOL"

out_dir <- "~/CESC_Network/6_OCTAD/6_2_ORA_signature/"

# --------------------------------------------------------------------
# Function: run_ora()
# - Runs enrichGO across clade × regulation × ontology
# - If simplify = TRUE, applies semantic collapsing with cutoff
# - Returns a single data.frame with all results (or an empty tibble)
# --------------------------------------------------------------------
run_ora <- function(consensus_tbl, universe,
                    ontologies = c("BP","CC","MF"),
                    min_gene_list_n = 10,
                    p_cut = 0.05, q_cut = 0.20,
                    key_type = "SYMBOL",
                    simplify = FALSE, simplify_cutoff = 0.95,
                    simplify_by = "p.adjust", simplify_measure = "Wang",
                    simplify_select = min) {
  
  out_rows <- list()
  
  for (cl in unique(consensus_tbl$clade)) {
    for (dr in c("up","down")) {
      
      genes_sym <- consensus_tbl %>%
        filter(clade == cl, regulation == dr) %>%
        pull(gene) %>%
        unique() %>%
        { .[!is.na(.) & nzchar(.)] }
      
      if (length(genes_sym) < min_gene_list_n) next
      
      for (ont in ontologies) {
        eg <- tryCatch(
          enrichGO(
            gene          = genes_sym,
            universe      = universe,
            OrgDb         = org.Hs.eg.db,
            keyType       = key_type,
            ont           = ont,
            pAdjustMethod = "BH",
            pvalueCutoff  = p_cut,
            qvalueCutoff  = q_cut,
            readable      = FALSE
          ),
          error = function(e) NULL
        )
        if (is.null(eg)) next
        
        if (isTRUE(simplify)) {
          eg <- tryCatch(
            simplify(eg, cutoff = simplify_cutoff, by = simplify_by,
                     select_fun = simplify_select, measure = simplify_measure),
            error = function(e) eg # fallback to non-simplified if simplify fails
          )
        }
        
        df <- as.data.frame(eg) %>%
          arrange(p.adjust, desc(Count)) %>%
          mutate(clade = cl, regulation = dr, ontology = ont, .before = 1)
        
        if (nrow(df) > 0) out_rows[[length(out_rows) + 1]] <- df
      }
    }
  }
  
  if (length(out_rows)) dplyr::bind_rows(out_rows) else
    tibble(clade=character(), regulation=character(), ontology=character(),
           ID=character(), Description=character(), GeneRatio=character(),
           BgRatio=character(), pvalue=double(), p.adjust=double(), qvalue=double(),
           geneID=character(), Count=integer())
}

# ---- Run FULL ----
ora_tbl_full <- run_ora(
  consensus_tbl, universe,
  ontologies = ontologies,
  min_gene_list_n = min_gene_list_n,
  p_cut = p_cut, q_cut = q_cut,
  key_type = key_type,
  simplify = FALSE
)
saveRDS(ora_tbl_full, file.path(out_dir, "6_2_2_Ora_tbl.rds"))

# ---- Run COLLAPSED (simplify = TRUE) ----
ora_tbl_grouped <- run_ora(
  consensus_tbl, universe,
  ontologies = ontologies,
  min_gene_list_n = min_gene_list_n,
  p_cut = p_cut, q_cut = q_cut,
  key_type = key_type,
  simplify = TRUE, simplify_cutoff = 0.25 #Change for semantic terms grouping
)

saveRDS(ora_tbl_grouped, file.path(out_dir, "6_2_3_Ora_tbl_0_25.rds"))

save.image("~/CESC_Network/6_OCTAD/6_2_ORA_signature/6_2_4_Image_ORA.RData")
cat("FULL   unique terms:", length(unique(ora_tbl_full$Description)),    "\n")
cat("GROUP  unique terms:", length(unique(ora_tbl_grouped$Description)), "\n")

### ###
#load("~/CESC_Network/6_OCTAD/6_2_ORA_signature/6_2_4_Image_ORA.RData")

#### Check different descriptions

a7_only <- setdiff(unique(ora_tbl_grouped$Description[ora_tbl_grouped$clade == "A7"]),
                   unique(ora_tbl_grouped$Description[ora_tbl_grouped$clade == "A9"]))

a9_only <- setdiff(unique(ora_tbl_grouped$Description[ora_tbl_grouped$clade == "A9"]),
                   unique(ora_tbl_grouped$Description[ora_tbl_grouped$clade == "A7"]))

# A7 
a7_only_df <- ora_tbl_grouped %>%
  filter(clade == "A7", Description %in% a7_only) %>%
  dplyr::select(Description, regulation)

#Description regulation
#GO:0008017 microtubule binding         up
#GO:0140666  annealing activity         up

# A9 
a9_only_df <- ora_tbl_grouped %>%
  filter(clade == "A9", Description %in% a9_only) %>%
  dplyr::select(Description, regulation)
#Description regulation
#GO:0140888          interferon-mediated signaling pathway         up
#GO:0003678                          DNA helicase activity         up
#GO:0140801               histone H2AXY142 kinase activity       down
#GO:0141147 intracellularly calcium-gated channel activity       down


