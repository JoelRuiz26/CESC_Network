library(ggplot2)
library(dplyr)

factors <- readRDS("~/2_Prepro_Data_TCGA/SNT_Full_A7A9/factors.rds")

# Filtrar y limpiar los datos
filtered_data <- factors %>%
  filter(HPV_clade %in% c("A7", "A9"),  # Mantener solo A7 y A9
         !is.na(figo_stage),            # Eliminar NA
         figo_stage != "Solid Tissue Normal") %>%  # Eliminar SNT
  
  # Extraer el número romano principal del estadio usando regex
  mutate(Stage_group = gsub("Stage ([IV]+).*", "Stage \\1", figo_stage),
         Stage_group = ifelse(grepl("^Stage [IV]+$", Stage_group), 
                              Stage_group, 
                              NA)) %>%  # Conservar solo los grupos principales
  filter(!is.na(Stage_group))  # Eliminar estadios que no coincidan

# Verificar los datos agrupados
table(filtered_data$HPV_clade, filtered_data$Stage_group)
#     Stage I Stage II Stage III Stage IV
#A7      36       17        10        2
#A9     106       43        32       15
table(filtered_data$Stage_group)

#Stage I  Stage II Stage III  Stage IV 
#142        60        42        17

ggplot(filtered_data, aes(x = Stage_group, fill = HPV_clade)) +
  geom_bar(position = "dodge") +
  labs(title = "Frecuencia de estadios FIGO agrupados por clado HPV",
       x = "Estadio",
       y = "Frecuencia",
       fill = "Clado filogenético HPV") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),  # Tamaño texto eje X
    axis.text.y = element_text(size = 20),                         # Tamaño texto eje Y
    axis.title = element_text(size = 22),                          # Tamaño títulos ejes
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5), # Título
    legend.text = element_text(size = 20),                         # Texto leyenda
    legend.title = element_text(size = 22)                          # Título leyenda
  ) +
  scale_fill_manual(values = c("A7" = "blue", "A9" = "red")) +
  scale_x_discrete(limits = c("Stage I", "Stage II", "Stage III", "Stage IV"))



# Calcular porcentajes por clado
percentage_data <- filtered_data %>%
  group_by(HPV_clade, Stage_group) %>%
  summarise(Count = n()) %>%
  group_by(HPV_clade) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()

# Mostrar la tabla de porcentajes
percentage_data

# A tibble: 8 × 4
# HPV_clade Stage_group Count Percentage
#<chr>     <chr>       <int>      <dbl>
#1 A7        Stage I        36      55.4 
#2 A7        Stage II       17      26.2 
#3 A7        Stage III      10      15.4 
#4 A7        Stage IV        2       3.08

#5 A9        Stage I       106      54.1 
#6 A9        Stage II       43      21.9 
#7 A9        Stage III      32      16.3 
#8 A9        Stage IV       15       7.65



ggplot(percentage_data, aes(x = HPV_clade, y = Percentage, fill = Stage_group)) +
  geom_col() +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)), 
            position = position_stack(vjust = 0.5), 
            color = "black", size = 8) +  # Aumentado a size = 8 (equivalente al aumento proporcional)
  labs(title = "Distribución porcentual de estadios FIGO por clado HPV",
       x = "Clado filogenético HPV",
       y = "Porcentaje",
       fill = "Estadio") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 20, face = "bold"),  # Tamaño texto eje X
    axis.text.y = element_text(size = 20),  # Tamaño texto eje Y
    axis.title.x = element_text(size = 22), # Título eje X
    axis.title.y = element_text(size = 22), # Título eje Y
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5), # Título
    legend.text = element_text(size = 20),   # Texto leyenda
    legend.title = element_text(size = 22),  # Título leyenda
    legend.position = "bottom"
  ) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) # Pequeño ajuste para mejor visualización
