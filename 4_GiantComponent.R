# Load required libraries
library(dplyr)
library(vroom)
library(igraph)

# Function to analyze the network without reordering (assuming input is already sorted)
analyze_network <- function(file_path, max_rows = 2000000) {
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
results_A7 <- analyze_network("~/4_MI/ARACNE-multicore/launch/Final_A7_Annot.sort")
results_A9 <- analyze_network("~/4_MI/ARACNE-multicore/launch/Final_A9_Annot.sort")

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

# ==================================
# SAVE R OBJECTS (OPTIONAL)
# ==================================
# Uncomment the following line if you want to save the R environment
save.image(file = "~/5_Global_Analisis/1_Subnetworks/ARACNE_fullnetworks/Results_cut_A7A9.RData")

# ==================================
# VERIFYING NODE COVERAGE THRESHOLD
# ==================================

# Load top 600,000 rows for A7
A7 <- vroom(file = "~/4_MI/ARACNE-multicore/launch/Final_A7_Annot.sort",
            col_names = c("GenA", "GenB", "MI"), show_col_types = FALSE) %>%
  slice(1:600000)

# Extract edges up to threshold (cutoff row)
A7_th <- A7 %>% slice(1:341747)
all_nodes_A7 <- unique(c(A7_th$GenA, A7_th$GenB))
length(all_nodes_A7)  # Should be 11240

# One row before threshold
A7_under_th <- A7 %>% slice(1:341746)
length(unique(c(A7_under_th$GenA, A7_under_th$GenB)))  # Should be 11239

# Same verification for A9
A9 <- vroom(file = "~/4_MI/ARACNE-multicore/launch/Final_A9_Annot.sort",
            col_names = c("GenA", "GenB", "MI"), show_col_types = FALSE) %>%
  slice(1:600000)

A9_th <- A9 %>% slice(1:549320)
length(unique(c(A9_th$GenA, A9_th$GenB)))  # Should be 11240

A9_under_th <- A9 %>% slice(1:549319)
length(unique(c(A9_under_th$GenA, A9_under_th$GenB)))  # Should be 11239

# ==================================
# ANALYZE CONNECTED COMPONENTS
# ==================================

# Function to calculate connected components
get_connected_components <- function(edge_list_df) {
  # Create an undirected graph
  g <- graph_from_data_frame(edge_list_df, directed = FALSE)
  
  # Compute connected components
  comps <- components(g)
  
  return(list(
    number_components = comps$no,
    membership = comps$membership,
    component_sizes = comps$csize
  ))
}

# Apply to A7 network
A7_components <- get_connected_components(A7_th)
cat("Number of connected components in A7:", A7_components$number_components, "\n")
#Number of connected components in A7: 1 

# Apply to A9 network
A9_components <- get_connected_components(A9_th)
cat("Number of connected components in A9:", A9_components$number_components, "\n")
#Number of connected components in A9: 1 
