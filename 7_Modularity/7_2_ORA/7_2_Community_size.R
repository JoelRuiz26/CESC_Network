#Explore_enrichment

### Cargar librerías ###
library(igraph)
library(tidyverse)
library(vroom)
library(ggplot2)
library(patchwork)

###Cargar resultados ###
ora_A7 <- vroom(file = "~/CESC_Network/7_Modularity/7_2_ORA/7_2_1_Ora7.tsv")
ora_A9 <- vroom(file = "~/CESC_Network/7_Modularity/7_2_ORA/7_2_1_Ora9.tsv")
# ====== PREPARAR DATA ======
# A7
df_A7 <- ora_A7 %>%
  distinct(Community_Name, Community_size) %>%
  arrange(desc(Community_size)) %>%
  mutate(Community_Name = factor(Community_Name, levels = Community_Name))

# A9
df_A9 <- ora_A9 %>%
  distinct(Community_Name, Community_size) %>%
  arrange(desc(Community_size)) %>%
  mutate(Community_Name = factor(Community_Name, levels = Community_Name))

# ====== PLOTS INDIVIDUALES ======
# ====== PLOTS INDIVIDUALES ======
plot_A7 <- ggplot(df_A7, aes(x = Community_Name, y = Community_size)) +
  geom_col(fill = "steelblue") +
  theme_minimal(base_size = 14) +
  labs(
    title = "HPV_A7",
    x = "Community Name",
    y = "Number of Genes"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold")
  )

plot_A9 <- ggplot(df_A9, aes(x = Community_Name, y = Community_size)) +
  geom_col(fill = "firebrick") +
  theme_minimal(base_size = 14) +
  labs(
    title = "HPV_A9",
    x = "Community Name",
    y = "Number of Genes"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold")
  )

# ====== COMBINAR CON TITULO GENERAL ======
combined_plot <- (plot_A7 + plot_A9) +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "Community Size",
    theme = theme(
      plot.title = element_text(size = 20, face = "bold")
    )
  )

# ====== MOSTRAR EN PANTALLA ======
combined_plot

# ====== GUARDAR EN PDF ======
ggsave("~/CESC_Network/7_Modularity/7_2_ORA/7_2_1_Community_Size_Barplots.png", combined_plot, width = 14, height = 7)
