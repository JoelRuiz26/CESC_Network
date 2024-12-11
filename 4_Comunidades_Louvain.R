# Cargar librerías necesarias
library(dplyr)
library(igraph)
library(vroom)
library(ggvenn)

set.seed(123)
# Cargar las redes
setwd("~/5_Network_features/")
graph_A7_10k <- read_graph("2_4_Graph_10k_A7.graphml", format = "graphml")
graph_A9_10k <- read_graph("2_4_Graph_10k_A9.graphml", format = "graphml")

# Función para detectar comunidades y extraer nodos por comunidad
procesar_comunidades <- function(graph) {
  # Detectar comunidades usando Louvain
  comunidades <- cluster_louvain(graph)
  
  # Agrupar nodos por comunidad
  tibble(
    comunidad = as.factor(membership(comunidades)),
    nodo = V(graph)$name
  ) %>%
    group_by(comunidad) %>%
    summarize(nodos = list(sort(nodo)), .groups = "drop")
}

# Procesar redes
comunidades_A7 <- procesar_comunidades(graph_A7_10k)
comunidades_A9 <- procesar_comunidades(graph_A9_10k)

# Identificar comunidades únicas y compartidas
set_comunidades_A7 <- comunidades_A7$nodos %>% lapply(sort) %>% unique()
set_comunidades_A9 <- comunidades_A9$nodos %>% lapply(sort) %>% unique()

comunidades_compartidas <- intersect(set_comunidades_A7, set_comunidades_A9)
comunidades_unicas_A7 <- setdiff(set_comunidades_A7, comunidades_compartidas)
comunidades_unicas_A9 <- setdiff(set_comunidades_A9, comunidades_compartidas)

# Crear tabla informativa
tabla_informativa <- tibble(
  Red = c("VPH-A7", "VPH-A9"),
  `Total de comunidades` = c(length(set_comunidades_A7), length(set_comunidades_A9)),
  `Comunidades únicas` = c(length(comunidades_unicas_A7), length(comunidades_unicas_A9)),
  `Comunidades compartidas` = c(length(comunidades_compartidas), length(comunidades_compartidas))
)

# Guardar la tabla informativa
vroom_write(tabla_informativa, "Tabla_Informativa_Comunidades.tsv")

# Crear diagrama de Venn
venn_data <- list(
  `VPH-A7` = lapply(set_comunidades_A7, paste, collapse = ","),
  `VPH-A9` = lapply(set_comunidades_A9, paste, collapse = ",")
)

Venn_identidadGenes <- ggvenn(
  venn_data,
  fill_color = c("#59A5F9", "red"),
  text_size = 15  # Aumentar el tamaño de texto dentro del diagrama
) +
  ggtitle("Diagrama de Venn de comunidades únicas y compartidas") +
  theme(
    plot.title = element_text(size = 30, face = "bold"), # Título más grande y en negrita
    text = element_text(size = 30)                      # Texto general más grande
  )

# Mostrar el gráfico
plot(Venn_identidadGenes)

# Guardar el diagrama en un archivo PNG
ggsave(
  filename = "~/5_Network_features/7_Chid_VenFinal.png",
  plot = Venn_identidadGenes,
  width = 12,
  height = 9,
  dpi = 300
)
