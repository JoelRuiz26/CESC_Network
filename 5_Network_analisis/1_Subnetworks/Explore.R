
# ==================================
# VERIFYING NODE COVERAGE THRESHOLD
# ==================================
path_A7 <- "/home/jruiz/CESC_Network/4_Network/ARACNE-multicore/launch/2_7_Full_counts_A7_annot.sort"
path_A9 <- "/home/jruiz/CESC_Network/4_Network/ARACNE-multicore/launch/2_8_Full_counts_A9_annot.sort"

# Load top 600,000 rows for A7
A7 <- vroom(file = path_A7,
            col_names = c("GenA", "GenB", "MI"), show_col_types = FALSE) %>%
  dplyr::slice(1:600000)

all_nodes_A7 <- unique(c(A7$GenA, A7$GenB))
length(all_nodes_A7)  # [1] 9699


# Extract edges up to threshold (cutoff row)
A7_th <- A7 %>% slice(1:500000)
all_nodes_A7 <- unique(c(A7_th$GenA, A7_th$GenB))
length(all_nodes_A7)  # [1] 8328


# Same verification for A9
A9 <- vroom(file = path_A9,
            col_names = c("GenA", "GenB", "MI"), show_col_types = FALSE) %>%
  slice(1:600000)
all_nodes_A9 <- unique(c(A9$GenA, A9$GenB))
length(all_nodes_A9)  # [1] 9258



A9_th <- A9 %>% slice(1:500000)

length(unique(c(A9_th$GenA, A9_th$GenB))) # [1] 7689



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

A7_components <- get_connected_components(A7)
cat("Number of connected components in A7:", A7_components$number_components, "\n")
#Number of connected components in A7: 52 

# Apply to A9 network
A9_components <- get_connected_components(A9)
cat("Number of connected components in A9:", A9_components$number_components, "\n")
#Number of connected components in A9: 66 


# Apply to A7 network
A7_components <- get_connected_components(A7_th)
cat("Number of connected components in A7:", A7_components$number_components, "\n")
#Number of connected components in A7: 119 

# Apply to A9 network
A9_components <- get_connected_components(A9_th)
cat("Number of connected components in A9:", A9_components$number_components, "\n")
#Number of connected components in A9: 149 

# ==================================
# SAVE_SUBGRAPH
# ==================================

A9_subgraph <- A9
A7_subgraph <- A7

write_graph(graph_from_data_frame(A7_subgraph, directed = FALSE), "~/5_Global_Analisis/1_Subnetworks/A7_subgraph.graphml", format = "graphml")
write_graph(graph_from_data_frame(A9_subgraph, directed = FALSE), "~/5_Global_Analisis/1_Subnetworks/A9_subgraph.graphml", format = "graphml")
