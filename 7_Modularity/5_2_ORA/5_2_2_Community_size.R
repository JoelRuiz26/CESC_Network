### Community size plot — only modules with GO enrichment (SIMPLIFIED ORA)
suppressPackageStartupMessages({
  library(tidyverse)
  library(vroom)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

ora_A7 <- vroom("~/CESC_Network/7_Modularity/5_2_ORA/5_2_2_Ora7_louv_max_SIMPLIFIED_0_2.tsv",
                show_col_types = FALSE)
ora_A9 <- vroom("~/CESC_Network/7_Modularity/5_2_ORA/5_2_2_Ora9_louv_max_SIMPLIFIED_0_2.tsv",
                show_col_types = FALSE)

df_A7 <- ora_A7 %>%
  distinct(Community_Name, Community_size) %>%
  arrange(desc(Community_size)) %>%
  mutate(Community_Name = factor(Community_Name, levels = Community_Name))

df_A9 <- ora_A9 %>%
  distinct(Community_Name, Community_size) %>%
  arrange(desc(Community_size)) %>%
  mutate(Community_Name = factor(Community_Name, levels = Community_Name))

base_theme <- theme_bw(base_family = "Times New Roman") +
  theme(
    axis.text.x        = element_text(size = 18, face = "bold", angle = 45, hjust = 1),
    axis.text.y        = element_text(size = 20, face = "bold"),
    axis.title.x       = element_text(size = 22, face = "bold", margin = margin(t = 8)),
    axis.title.y       = element_text(size = 22, face = "bold", margin = margin(r = 8)),
    plot.title         = element_text(size = 26, face = "bold", hjust = 0.5,
                                      margin = margin(b = 8)),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border       = element_rect(color = "black", linewidth = 0.8)
  )

pA7 <- ggplot(df_A7, aes(x = Community_Name, y = Community_size)) +
  geom_col(fill = "#003366", alpha = 0.88, width = 0.8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), labels = comma) +
  labs(title = sprintf("HPV-A7  (%d modules with GO terms)", nrow(df_A7)),
       x = NULL, y = "Number of genes") +
  base_theme

pA9 <- ggplot(df_A9, aes(x = Community_Name, y = Community_size)) +
  geom_col(fill = "#CC2200", alpha = 0.88, width = 0.8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), labels = comma) +
  labs(title = sprintf("HPV-A9  (%d modules with GO terms)", nrow(df_A9)),
       x = NULL, y = "Number of genes") +
  base_theme

combined <- pA7 + pA9 + plot_layout(ncol = 2) +
  plot_annotation(
    title = "Size of Enriched Louvain Modules — HPV Co-expression Networks",
    theme = theme(plot.title = element_text(size = 30, face = "bold", hjust = 0.5,
                                             family = "Times New Roman",
                                             margin = margin(b = 12)))
  )

ggsave("~/CESC_Network/7_Modularity/5_2_ORA/5_2_2_1_Community_Size_louvain_max_simpl.png",
       combined, width = 22, height = 10, dpi = 300, bg = "white")
message("[OK] Saved: 5_2_2_1_Community_Size_louvain_max_simpl.png")
