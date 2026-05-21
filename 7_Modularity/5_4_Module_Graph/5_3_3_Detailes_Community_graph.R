library(igraph)

# Directorios
input_dir  <- "/home/joelr/CESC_Network/5_Modularity/5_4_Module_Graph"
output_dir <- "/STORAGE/csbig/jruiz/cesc_network"

# Leer los objetos RDS
g_A9 <- readRDS(file.path(input_dir,
                          "5_3_2_1_Community_Graph_A9_louv_with_GO.rds"))
g_A7 <- readRDS(file.path(input_dir,
                          "5_3_2_1_Community_Graph_A7_louv_with_GO.rds"))

# (Opcional) verifica que sean igraph
print(g_A9)
print(g_A7)

# -------------------------------------------------
# Función helper para aplicar todas las transformaciones
# -------------------------------------------------
process_graph <- function(g) {
  # 1) Asignar size = 1 a nodos GO (BP, MF, CC)
  if (!is.null(V(g)$Ontology) && !is.null(V(g)$size)) {
    V(g)$size <- ifelse(
      V(g)$Ontology %in% c("BP", "MF", "CC"),
      1,
      V(g)$size
    )
  }
  
  # 2) Reemplazar NA por "Community" en Ontology
  if (!is.null(V(g)$Ontology)) {
    V(g)$Ontology <- ifelse(
      is.na(V(g)$Ontology),
      "Community",
      V(g)$Ontology
    )
  }
  
  # 3) Rescalar avg_pagerank al rango [0.1, 1] en nueva columna
  if (!is.null(V(g)$avg_pagerank)) {
    ap  <- V(g)$avg_pagerank
    rng <- range(ap, na.rm = TRUE)
    
    # Evitar división por cero si todos los valores son iguales
    if (rng[1] == rng[2]) {
      ap_rescaled <- rep(1, length(ap))  # todo a 1 (o 0.1, si prefieres)
      ap_rescaled[is.na(ap)] <- NA_real_
    } else {
      ap_rescaled <- (ap - rng[1]) / (rng[2] - rng[1]) * (1 - 0.1) + 0.1
      ap_rescaled[is.na(ap)] <- NA_real_
    }
    
    V(g)$avg_pagerank_rescaled <- ap_rescaled
  }
  
  # 4) Crear NES_range: NEGATIVE / POSITIVE / ZERO / NA
  if (!is.null(V(g)$NES)) {
    nes <- V(g)$NES
    nes_range <- rep(NA_character_, length(nes))
    
    nes_range[!is.na(nes) & nes < 0] <- "NEGATIVE"
    nes_range[!is.na(nes) & nes > 0] <- "POSITIVE"
    nes_range[!is.na(nes) & nes == 0] <- "ZERO"
    
    V(g)$NES_range <- nes_range
  }
  
  return(g)
}

# Aplicar a ambas redes
g_A7 <- process_graph(g_A7)
g_A9 <- process_graph(g_A9)

# (Opcional) checar rápido
 table(V(g_A7)$Ontology)
 table(V(g_A9)$Ontology)
 summary(V(g_A7)$avg_pagerank_rescaled)
 summary(V(g_A9)$avg_pagerank_rescaled)
 table(V(g_A7)$NES_range, useNA = "ifany")
 table(V(g_A9)$NES_range, useNA = "ifany")

# -------------------------------------------------
# Exportar a GraphML (ya con todos los cambios)
# -------------------------------------------------
write_graph(
  g_A9,
  file.path(output_dir, "Community_Graph_A9_louv_with_GO.graphml"),
  format = "graphml"
)

write_graph(
  g_A7,
  file.path(output_dir, "Community_Graph_A7_louv_with_GO.graphml"),
  format = "graphml"
)

# (Opcional) guardar también las versiones RDS actualizadas
saveRDS(g_A7, file.path(output_dir, "Community_Graph_A7_louv_with_GO_updated.rds"))
saveRDS(g_A9, file.path(output_dir, "Community_Graph_A9_louv_with_GO_updated.rds"))
