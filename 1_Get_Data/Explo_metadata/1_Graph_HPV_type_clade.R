#Graph for create HPV-Clade distribution plot in RNA samples
#Joel Ruiz Hernandez
#Input: metadata (1_2_Metadata.tsv) from 1_Get_data_RNA.R script
#Output: graph "Distribution of HPV Types Present in RNAseq Samples of CC"

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(vroom)

# Set working directory
setwd("~/1_Get_Data_TCGA/Explo_metadata/")

metadata <- vroom(file = "~/1_Get_Data_TCGA/1_2_Metadata.tsv")
# Filter data, group by clade and HPV type, and count samples
combined_counts <- metadata %>%
  dplyr::filter(HPV_clade != "Solid Tissue Normal") %>%
  dplyr::count(HPV_clade, HPV_type)

# Get the total samples by clade and reorder the clades
clade_order <- combined_counts %>%
  group_by(HPV_clade) %>%
  summarise(total_count = sum(n)) %>%
  arrange(desc(total_count)) %>%
  pull(HPV_clade)

# Reorder clades and HPV types based on the number of samples
combined_counts <- combined_counts %>%
  mutate(HPV_clade = factor(HPV_clade, levels = clade_order)) %>%
  arrange(HPV_clade, desc(n)) %>%
  mutate(HPV_label = factor(paste(HPV_type, "(n =", n, ")"), levels = unique(paste(HPV_type, "(n =", n, ")"))))

combined_counts <- combined_counts %>% dplyr::slice(1:10) #Se seleccionan solo los HPV de alto riesgo, se excluye clado "Otro"
# A tibble: 10 × 4
#HPV_clade HPV_type     n HPV_label       
#<fct>     <chr>    <int> <fct>           
#  1 A9        HPV16      166 HPV16 (n = 166 )
#2 A9        HPV33        8 HPV33 (n = 8 )  
#3 A9        HPV52        8 HPV52 (n = 8 )  
#4 A9        HPV31        7 HPV31 (n = 7 )  
#5 A9        HPV58        7 HPV58 (n = 7 )  
#6 A9        HPV35        6 HPV35 (n = 6 )  
#7 A7        HPV18       38 HPV18 (n = 38 ) 
#8 A7        HPV45       20 HPV45 (n = 20 ) 
#9 A7        HPV39        5 HPV39 (n = 5 )  
#10 A7        HPV59        3 HPV59 (n = 3 )  

# Add the total samples by clade
total_by_clade <- combined_counts %>%
  group_by(HPV_clade) %>%
  summarise(total_count = sum(n))
# Plot
plot <- ggplot(combined_counts, aes(x = HPV_clade, y = n, fill = HPV_label)) +
  geom_bar(stat = "identity", position = "stack", width = 0.9) +  # Thinner bars
  scale_fill_viridis_d(option = "C", end = 0.9) +  # Distinct color palette
  geom_text(data = total_by_clade, aes(x = HPV_clade, y = total_count, label = paste("n =", total_count)), 
            inherit.aes = FALSE, vjust = -0.5, size = 5) +  # Total labels by clade with "n="
  scale_x_discrete(expand = c(0.5, 0)) +  # Reduce space between bars
  labs(
    title = "Distribución de tipos de HPV presentes en muestras de RNAseq de CC",
    x = "Clado de VPH",
    y = "Frecuencia",
    fill = "Tipo de VPH"
  ) +
  theme_minimal(base_size = 18) +  # Larger text
  theme(
    plot.title = element_text(face = "bold", hjust = 0.1, size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels if needed
    legend.position = "right"
  )

# Show plot
print(plot)

# Save plot
ggsave("1_1_Graph_HPV_type_clade.png", plot, width = 9, height = 7)
