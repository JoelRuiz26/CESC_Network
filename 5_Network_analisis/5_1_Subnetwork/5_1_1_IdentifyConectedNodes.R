# --------------------------------------------------------
# Determine the minimal number of edges required to 
# cover all nodes in Co-expression Networks (ARACNE sorted)
#
# Author: Joel Ruiz Hernandez
#
# Description:
#   - For each network: identifies the row (edge count) 
#     at which all unique genes are included at least once.
# Input: 
#   -Edge list (.sort) with highest mutual information 
# Output:
#   - Reports summary statistics to the console.
#   - Saves the maximum cutoff value between both networks
#     as a single .rds object.
# --------------------------------------------------------

# --------------------------------------------------------
# Load required libraries
# --------------------------------------------------------
library(dplyr)
library(vroom)
library(igraph)
library(tools)

# --------------------------------------------------------
# Function to analyze a single network file
# Assumes the input file is pre-sorted by MI
# --------------------------------------------------------
#                               Set a huge slice of edgelist that could have all nodes
analyze_network <- function(file_path, max_rows = 10000000) {
  
  # Load the edge list (up to max_rows)
  edge_list <- vroom(
    file_path, 
    col_names = c("GenA", "GenB", "MI"),
    show_col_types = FALSE
  ) %>%
    na.omit() %>%
    slice(1:max_rows)
  
  # Extract all unique nodes (genes)
  all_nodes <- unique(c(edge_list$GenA, edge_list$GenB))
  total_nodes <- length(all_nodes)
  
  # Determine the cutoff row where all nodes are covered
  covered_nodes <- character()
  cutoff_row <- NA
  
  for (i in 1:nrow(edge_list)) {
    covered_nodes <- unique(c(covered_nodes, edge_list$GenA[i], edge_list$GenB[i]))
    
    if (length(covered_nodes) == total_nodes) {
      cutoff_row <- i
      break
    }
  }
  
  # Return structured summary
  return(list(
    network_name = file_path_sans_ext(basename(file_path)),
    total_nodes = total_nodes,
    cutoff_row = cutoff_row,
    covered_nodes = length(covered_nodes),
    coverage_percent = round(length(covered_nodes) / total_nodes * 100, 2),
    sample_size = nrow(edge_list)
  ))
}


# --------------------------------------------------------
# Analyze both networks
# --------------------------------------------------------
results_A7 <- analyze_network(
  "~/CESC_Network/4_Network/ARACNE-multicore/launch/2_7_Full_counts_A7_annot.sort"
)

results_A9 <- analyze_network(
  "~/CESC_Network/4_Network/ARACNE-multicore/launch/2_8_Full_counts_A9_annot.sort"
)

# --------------------------------------------------------
# Function to generate and print a detailed report
# --------------------------------------------------------
generate_report <- function(results) {
  report <- paste0(
    "\n════════════════════════════════════\n",
    " NETWORK ANALYSIS REPORT - ", results$network_name, "\n",
    "════════════════════════════════════\n",
    "• Total unique nodes: ", results$total_nodes, "\n",
    "• Rows analyzed: ", results$sample_size, "\n",
    "• Cutoff row found: ", 
    ifelse(is.na(results$cutoff_row), "NOT COMPLETE", results$cutoff_row), "\n",
    "• Nodes covered: ", results$covered_nodes, "/", results$total_nodes, "\n",
    "• Coverage percent: ", results$coverage_percent, "%\n",
    "════════════════════════════════════\n"
  )
  
  cat(report)
  invisible(report)
}

# --------------------------------------------------------
# Print summary reports to console
# --------------------------------------------------------
generate_report(results_A7)
generate_report(results_A9)

# --------------------------------------------------------
# Determine the maximum cutoff value across both networks
# This ensures consistent edge density in comparative analysis
# --------------------------------------------------------
cutoff_A7 <- results_A7$cutoff_row
cutoff_A9 <- results_A9$cutoff_row

max_cutoff_links <- max(c(cutoff_A7, cutoff_A9), na.rm = TRUE)

cat("\n════════════════════════════════════\n")
cat(" FINAL COMMON THRESHOLD SELECTED\n")
cat("════════════════════════════════════\n")
cat("Maximum cutoff row between networks: ", max_cutoff_links, "\n")
cat("════════════════════════════════════\n")

# --------------------------------------------------------
# Save the single cutoff value as .rds
# This can be easily loaded in downstream analyses:
#   readRDS(".../max_cutoff_links.rds")
# --------------------------------------------------------

# --------------------------------------------------------
# Save both edgelists sliced to the common maximum cutoff
# --------------------------------------------------------
A7_common <- vroom(
  "~/CESC_Network/4_Network/ARACNE-multicore/launch/2_7_Full_counts_A7_annot.sort",
  col_names = c("GenA", "GenB", "MI"),
  show_col_types = FALSE
) %>% na.omit() %>% dplyr::slice(1:max_cutoff_links)

A9_common <- vroom(
  "~/CESC_Network/4_Network/ARACNE-multicore/launch/2_8_Full_counts_A9_annot.sort",
  col_names = c("GenA", "GenB", "MI"),
  show_col_types = FALSE
) %>% na.omit() %>% dplyr::slice(1:max_cutoff_links)

saveRDS(A7_common, "~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_1_2_A7_edgelist_allnodes.rds")
saveRDS(A9_common, "~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_1_2_A9_edgelist_allnodes.rds")
cat("[INFO] Saved both A7 and A9 edgelists with common cutoff:", max_cutoff_links, "\n")

saveRDS(
  cutoff_A7, 
  file = "~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_1_1_A7_cutoff_links.rds"
)

saveRDS(
  cutoff_A9, 
  file = "~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_1_1_A9_cutoff_links.rds"
)


cat("Saved object 'max_cutoff_links' as .rds file for future use.\n")


#════════════════════════════════════
#NETWORK ANALYSIS REPORT - 2_7_Full_counts_A7_annot
#════════════════════════════════════
#• Total unique nodes: 11544
#• Rows analyzed: 10000000
#• Cutoff row found: 2690573
#• Nodes covered: 11544/11544
#• Coverage percent: 100%
#════════════════════════════════════

#════════════════════════════════════
#NETWORK ANALYSIS REPORT - 2_8_Full_counts_A9_annot
#════════════════════════════════════
#• Total unique nodes: 11544
#• Rows analyzed: 10000000
#• Cutoff row found: 4967805
#• Nodes covered: 11544/11544
#• Coverage percent: 100%
#════════════════════════════════════

#════════════════════════════════════
#FINAL COMMON THRESHOLD SELECTED
#════════════════════════════════════
#Maximum cutoff row between networks:  4967805 
#════════════════════════════════════


