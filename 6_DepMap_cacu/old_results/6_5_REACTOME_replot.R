############################################################
## 6_5_REACTOME_replot.R
##
## Re-genera SOLO el plot de 6_5_MSigDB_REACTOME_barplot
## a partir del cache 6_5_REACTOME_oradata.rds, sin volver
## a correr enrichGO/enricher (que es lo lento de 6_5_MSigDB_ORA.R).
##
## Uso: ajustar los parámetros de la sección "Parámetros del
## plot" y correr este script para iterar rápido sobre estilo.
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(patchwork)
})

out_dir <- "~/CESC_Network/6_DepMap_cacu"

reorder_within <- function(x, by, within, ...) {
  new_x <- paste(x, within, sep = "___")
  stats::reorder(new_x, by, ...)
}
scale_y_reordered <- function(..., sep = "___") {
  ggplot2::scale_y_discrete(
    labels = function(x) sub(paste0(sep, ".*"), "", x), ...)
}

GROUP_COL <- c("HPV-A7" = "#1A6FB5", "HPV-A9" = "#C0392B")

make_panel <- function(grp, accent_col, df_grp, col_lo, col_hi, font_mult = 1) {
  if (nrow(df_grp) == 0) return(NULL)

  ggplot(df_grp, aes(x = GeneRatio_num, y = term_ord, fill = log_padj)) +
    geom_col(width = 0.68, colour = NA) +
    scale_y_reordered() +
    scale_fill_gradient(
      low  = col_lo, high = col_hi,
      name = expression(bold(-log[10](p[adj])))
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.06)),
                       n.breaks = 4,
                       guide = guide_axis(check.overlap = TRUE)) +
    labs(title = grp, x = "Gene Ratio", y = NULL) +
    theme_classic(base_size = 16 * font_mult) +
    theme(
      plot.title         = element_text(face = "bold", size = 19 * font_mult,
                                        colour = accent_col),
      axis.text.y        = element_text(size = 14 * font_mult, colour = "grey10",
                                        face = "bold"),
      axis.text.x        = element_text(size = 13 * font_mult, face = "bold",
                                        colour = "grey20"),
      axis.title.x       = element_text(size = 14 * font_mult, face = "bold",
                                        colour = "grey30"),
      axis.line          = element_line(colour = "grey70"),
      legend.title       = element_text(size = 13 * font_mult, face = "bold"),
      legend.text        = element_text(size = 12 * font_mult, face = "bold"),
      panel.grid.major.x = element_line(colour = "grey93", linewidth = 0.35)
    )
}

# ══════════════════════════════════════════════════════════════
## Cache
# ══════════════════════════════════════════════════════════════
cache <- readRDS(file.path(out_dir, "6_5_REACTOME_oradata.rds"))
ora_comb_raw <- cache$ora_comb
meta         <- cache$meta
universe_n   <- cache$universe_n

# ══════════════════════════════════════════════════════════════
## Parámetros del plot (ajustar aquí para iterar)
# ══════════════════════════════════════════════════════════════
top_n      <- meta$top_n
font_mult  <- 2.2
wrap_width <- 24
ncol_panels <- meta$ncol
plot_w     <- 26
row_h      <- 1.25

# ══════════════════════════════════════════════════════════════
## Plot
# ══════════════════════════════════════════════════════════════
ora_comb <- ora_comb_raw %>%
  mutate(Desc_wrap = str_wrap(Description, width = wrap_width))

plot_df <- ora_comb %>%
  group_by(group) %>%
  slice_max(GeneRatio_num, n = top_n, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(term_ord = reorder_within(Desc_wrap, GeneRatio_num, group))

p_A7 <- make_panel("HPV-A7", GROUP_COL["HPV-A7"],
                   plot_df %>% filter(group == "HPV-A7"),
                   meta$lo, meta$hi, font_mult)
p_A9 <- make_panel("HPV-A9", GROUP_COL["HPV-A9"],
                   plot_df %>% filter(group == "HPV-A9"),
                   meta$lo, meta$hi, font_mult)
panels <- Filter(Negate(is.null), list(p_A7, p_A9))

title_mult  <- min(font_mult, 1.4)
title_width <- if (ncol_panels == 2) 70 else 42

combined <- wrap_plots(panels, ncol = ncol_panels) +
  plot_annotation(
    title   = str_wrap(sprintf("%s — Hub ∩ DE genes | Cervical Cancer (CESC)", meta$label), width = title_width),
    caption = sprintf(
      "MSigDB · Universe: %d genes · BH-adjusted p < 0.05 · Top %d terms by Gene Ratio",
      universe_n, top_n
    ),
    theme = theme(
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 17 * title_mult),
      plot.caption = element_text(colour = "grey50", size = 11 * title_mult, hjust = 0.5)
    )
  )

n_grp  <- plot_df %>% count(group) %>% pull(n) %>% max(na.rm = TRUE)
plot_h <- n_grp * row_h + 4

ggsave(file.path(out_dir, paste0(meta$file, ".pdf")),
       plot = combined, width = plot_w, height = plot_h, device = cairo_pdf)
ggsave(file.path(out_dir, paste0(meta$file, ".png")),
       plot = combined, width = plot_w, height = plot_h, dpi = 150)

cat(sprintf("Guardado: %s.pdf/png (%.1f x %.1f in)\n", meta$file, plot_w, plot_h))
