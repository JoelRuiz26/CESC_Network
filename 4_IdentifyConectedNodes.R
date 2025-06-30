# Load required libraries
library(dplyr)
library(vroom)
library(igraph)

# Function to analyze the network without reordering (assuming input is already sorted)
analyze_network <- function(file_path, max_rows = 300000) {
  # Load the edge list (up to max_rows rows)
  edge_list <- vroom(file_path, 
                     col_names = c("GenA", "GenB", "MI"),
                     show_col_types = FALSE) %>%
    na.omit() %>%
    slice(1:max_rows)  # Limit number of rows if necessary
  
  # Get all unique nodes in the network
  all_nodes <- unique(c(edge_list$GenA, edge_list$GenB))
  total_nodes <- length(all_nodes)
  
  # Determine cutoff row where all nodes are covered
  covered_nodes <- character()
  cutoff_row <- NA
  
  for (i in 1:nrow(edge_list)) {
    covered_nodes <- unique(c(covered_nodes, edge_list$GenA[i], edge_list$GenB[i]))
    
    if (length(covered_nodes) == total_nodes) {
      cutoff_row <- i
      break
    }
  }
  
  # Return a structured summary
  return(list(
    network_name = tools::file_path_sans_ext(basename(file_path)),
    total_nodes = total_nodes,
    cutoff_row = cutoff_row,
    covered_nodes = length(covered_nodes),
    coverage_percent = round(length(covered_nodes) / total_nodes * 100, 2),
    sample_size = nrow(edge_list)
  ))
}

# Analyze both networks
results_A7 <- analyze_network("~/CESC_Network/4_Network/ARACNE-multicore/launch/2_7_Full_counts_A7_annot.sort")
results_A9 <- analyze_network("~/CESC_Network/4_Network/ARACNE-multicore/launch/2_8_Full_counts_A9_annot.sort")

# ==================================
# PRINT REPORTS TO CONSOLE
# ==================================

# Function to generate and print a detailed report
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
  
  # Print to console only
  cat(report)
  
  return(invisible(report))
}

# Print the reports
generate_report(results_A7)
generate_report(results_A9)


#════════════════════════════════════
#NETWORK ANALYSIS REPORT - Final_A7_Annot
#════════════════════════════════════
#• Total unique nodes: 11240
#• Rows analyzed: 1000000
#• Cutoff row found: 341747
#• Nodes covered: 11240/11240
#• Coverage percent: 100%
#════════════════════════════════════

#════════════════════════════════════
#NETWORK ANALYSIS REPORT - Final_A9_Annot
#════════════════════════════════════
#• Total unique nodes: 11240
#• Rows analyzed: 1000000
#• Cutoff row found: 549320
#• Nodes covered: 11240/11240R  
#• Coverage percent: 100%
#════════════════════════════════════

# ==================================
# SAVE R OBJECTS (OPTIONAL)
# ==================================
# Uncomment the following line if you want to save the R environment
save.image(file = "~/CESC_Network/5_Network_analisis/1_Subnetworks/Results_cut_A7A9.RData")
#load(file = "~/5_Global_Analisis/1_Subnetworks/Results_cut_A7A9.RData")
