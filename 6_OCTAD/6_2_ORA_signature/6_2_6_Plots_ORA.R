# ================================
# Libraries
# ================================
suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(forcats)
  library(ggplot2); library(scales)
})

# ================================
# Helper: save PNG + PDF (vector)
# - Ensures same width/height for both
# - Uses Cairo PDF to keep text metrics consistent
# ================================
save_png_pdf <- function(path_png, plot, width, height,
                         units = "in", dpi = 1200, limitsize = FALSE,
                         use_cairo = TRUE) {
  
  # Forzar fondo blanco en el objeto
  plot_white <- plot +
    theme(
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
  
  # 1) PNG (con bg blanco explícito)
  ggplot2::ggsave(
    filename = path_png, plot = plot_white,
    width = width, height = height, units = units,
    dpi = dpi, limitsize = limitsize, bg = "white"
  )
  
  # 2) PDF (vector; Cairo para métricas consistentes)
  path_pdf <- sub("\\.png$", ".pdf", path_png)
  if (use_cairo) {
    ggplot2::ggsave(
      filename = path_pdf, plot = plot_white,
      width = width, height = height, units = units,
      device = grDevices::cairo_pdf, limitsize = limitsize
    )
  } else {
    ggplot2::ggsave(
      filename = path_pdf, plot = plot_white,
      width = width, height = height, units = units,
      limitsize = limitsize
    )
  }
  invisible(path_pdf)
}


# ================================
# I/O
# ================================
ora_tbl <- readRDS("~/CESC_Network/6_OCTAD/6_2_ORA_signature/6_2_3_Ora_tbl_0_25.rds")
out_dir <- "~/CESC_Network/6_OCTAD/6_2_ORA_signature/6_2_7_Output_plots/"

# ================================
# Utils
# ================================
ratio_num <- function(x) sapply(strsplit(x, "/"), \(z) as.numeric(z[1]) / as.numeric(z[2]))

# ================================
# Clean & transform
# ================================
ora_clean <- ora_tbl %>%
  mutate(
    GeneRatio_num    = ratio_num(GeneRatio),
    minusLog10Padj   = -log10(p.adjust),
    Description_wrap = stringr::str_wrap(Description, width = 55),
    clade            = factor(clade, levels = sort(unique(clade))),
    regulation       = factor(regulation, levels = c("up", "down")),
    ontology         = factor(ontology, levels = c("BP","CC","MF")),
    clade_lbl        = dplyr::recode(clade, A7 = "HPV-A7", A9 = "HPV-A9")
  )

# Top terms per clade/regulation/ontology
top40 <- ora_clean %>%
  group_by(clade_lbl, regulation, ontology) %>%
  arrange(p.adjust, desc(Count), .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup()

# Nice breaks for x-axis
xb <- pretty(c(0, max(top40$GeneRatio_num, na.rm = TRUE)), n = 4)

# ================================
# Plot
# ================================
p_all <- ggplot(top40,
                aes(x = GeneRatio_num,
                    y = forcats::fct_reorder(Description_wrap, GeneRatio_num),
                    size = Count, color = minusLog10Padj)) +
  geom_point(alpha = 0.95) +
  facet_grid(ontology ~ interaction(clade_lbl, regulation, sep = " / "),
             scales = "free_y", space = "free_y") +
  labs(x = "GeneRatio", y = NULL, color = "-log10(p.adjust)", size = "Genes") +
  scale_color_viridis_c(direction = -1) +
  scale_size_area(max_size = 12) +
  scale_x_continuous(
    breaks = xb,
    labels = scales::label_number(accuracy = 0.02),
    expand = expansion(mult = c(0.01, 0.03)),
    guide  = guide_axis(check.overlap = TRUE)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text   = element_text(size = 14, face = "bold"),
    strip.text.y = element_text(angle = 0, face = "bold"),
    axis.text.y  = element_text(size = 13, face = "bold", lineheight = 1.3),
    panel.spacing.y    = grid::unit(12, "pt"),
    panel.spacing.x    = grid::unit(15, "pt"),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor   = element_blank(),
    plot.margin        = margin(10, 20, 10, 40)
  ) +
  coord_cartesian(clip = "off")

# ================================
# Save: PNG + PDF (same geometry)
# ================================
save_png_pdf(file.path(out_dir, "6_2_8_Ora_dotplot_clean.png"),
             p_all, width = 16, height = 9, dpi = 1200, limitsize = FALSE)

p_all

# (optional) Save workspace
save.image("~/CESC_Network/6_OCTAD/6_2_ORA_signature/6_2_7_Output_plots/6_2_9_Image_plots.RData")
