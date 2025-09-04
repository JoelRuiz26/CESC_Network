# ================================
# Plots de ORA (top-10 por panel)
# - Entrada: ora_tbl (RDS)
# - Salidas: 
#    1) 6_2_2_Ora_dotplot_ALL_top10.png  (facet)
#    2) PNGs por combinación clade_reg_ont (p.ej., A7_up_BP_top10.png)
# ================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(forcats)
})

# ---- Cargar resultados ORA ----
ora_tbl <- readRDS("~/CESC_Network/6_OCTAD/6_2_ORA_signature/6_2_2_Ora_tbl.rds")

# Carpeta de salida
out_dir <- "~/CESC_Network/6_OCTAD/6_2_ORA_signature/plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Preparar columnas numéricas útiles ----
# GeneRatio viene como "k/N": lo convertimos a numérico
ratio_num <- function(x) {
  sapply(strsplit(x, "/"), function(z) as.numeric(z[1]) / as.numeric(z[2]))
}

ora_clean <- ora_tbl %>%
  mutate(
    GeneRatio_num = ratio_num(GeneRatio),
    minusLog10Padj = -log10(p.adjust),
    # Acortar descripciones para que quepan bien sin perder legibilidad
    Description_wrap = str_wrap(Description, width = 50),
    clade = factor(clade, levels = sort(unique(clade))),
    regulation = factor(regulation, levels = c("up","down")),
    ontology = factor(ontology, levels = c("BP","CC","MF"))
  )

# ---- Top-10 por clade × regulation × ontology (ordenado por p.adjust) ----
top10 <- ora_clean %>%
  group_by(clade, regulation, ontology) %>%
  arrange(p.adjust, desc(Count), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

# ---- DOTPLOT FACETEADO (todos los paneles con top-10) ----
p_all <- ggplot(
  top10,
  aes(x = GeneRatio_num,
      y = fct_reorder(Description_wrap, GeneRatio_num),
      size = Count,
      color = minusLog10Padj)
) +
  geom_point() +
  facet_grid(ontology ~ interaction(clade, regulation, sep = " / "),
             scales = "free_y", space = "free_y") +
  labs(
    x = "GeneRatio",
    y = NULL,
    color = "-log10(p.adjust)",
    size = "Genes",
    title = "GO ORA (Top-10 por panel)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.spacing.y = unit(6, "pt"),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(out_dir, "6_2_2_Ora_dotplot_ALL_top10.png"),
       p_all, width = 12, height = 10, dpi = 300)

# ---- PNGs individuales por combinación clade × regulation × ontology ----
#     Mismo estilo de dotplot, 10 términos por combinación
make_safe <- function(x) gsub("[^A-Za-z0-9_]+", "_", x)

combos <- top10 %>%
  distinct(clade, regulation, ontology) %>%
  arrange(clade, regulation, ontology)

for (i in seq_len(nrow(combos))) {
  cl  <- combos$clade[i]
  dr  <- combos$regulation[i]
  ont <- combos$ontology[i]
  
  df <- top10 %>%
    filter(clade == cl, regulation == dr, ontology == ont) %>%
    # reordenar ejes por GeneRatio_num para este panel
    mutate(Description_wrap = fct_reorder(Description_wrap, GeneRatio_num))
  
  if (nrow(df) == 0) next
  
  p <- ggplot(
    df,
    aes(x = GeneRatio_num,
        y = Description_wrap,
        size = Count,
        color = minusLog10Padj)
  ) +
    geom_point() +
    labs(
      x = "GeneRatio",
      y = NULL,
      color = "-log10(p.adjust)",
      size = "Genes",
      title = paste0("GO ORA Top-10 — ", cl, " / ", dr, " / ", as.character(ont))
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y = element_text(size = 9),
      plot.title = element_text(face = "bold")
    )
  
  fname <- sprintf("%s_%s_%s_top10.png",
                   make_safe(as.character(cl)),
                   make_safe(as.character(dr)),
                   make_safe(as.character(ont)))
  ggsave(file.path(out_dir, fname), p, width = 9, height = 6, dpi = 300)
}

cat("Listo.\n- Facet: ", file.path(out_dir, "6_2_2_Ora_dotplot_ALL_top10.png"),
    "\n- PNGs individuales por combinación en: ", out_dir, "\n")
