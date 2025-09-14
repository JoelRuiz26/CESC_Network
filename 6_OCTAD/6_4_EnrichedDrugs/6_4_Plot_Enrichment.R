suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(ggthemes)
  library(viridis)
  library(grid)  # for unit()
})


# --- Paths ---
out_base <- "~/CESC_Network/6_OCTAD/6_4_EnrichedDrugs"
rds_in   <- file.path(out_base, "6_4_1_enrichment_sig_with_n_present.rds")

# --- Load table (must have: target, score, padj, n_present, signature, target_type) ---
enrich_sig_n <- readRDS(rds_in)
enrich_sig_n_score05 <- enrich_sig_n %>% filter(score >= 0.5) %>% filter(n_present > 3)

# A tibble: 6 × 7
#target score     p  padj signature target_type    n_present
#<chr>  <dbl> <dbl> <dbl> <chr>     <chr>              <int>
#1 ABCB11 0.641     0     0 HPV-A7    chembl_targets         2
#2 ABCC1  0.673     0     0 HPV-A7    chembl_targets         2
#3 ABCC4  0.589     0     0 HPV-A7    chembl_targets         1
#4 ABCC8  0.791     0     0 HPV-A7    chembl_targets         1


# --- Replace signature values with nicer labels ---
enrich_sig_n <- enrich_sig_n %>%
  dplyr::mutate(
    signature = dplyr::case_when(
      signature == "A7" ~ "HPV-A7",
      signature == "A9" ~ "HPV-A9",
      TRUE ~ signature  # keep other values unchanged
    )
  )


# Helper: shorten long labels for Y axis



# Helper: shorten long labels for Y axis
shorten <- function(x, width = 60) stringr::str_trunc(x, width = width, side = "right", ellipsis = "…")

plot_bubble_horizontal_pretty <- function(df,
                                          signatures      = c("HPV-A7","HPV-A9"),
                                          target_types    = c("mesh"),
                                          alpha           = 0.05,
                                          min_present     = 3,
                                          pdf_out         = file.path(out_base, "6_4_plot_bubble.pdf"),
                                          png_out         = file.path(out_base, "6_4_plot_bubble.png"),
                                          show            = TRUE,
                                          base_size       = 16,
                                          # slimmer & taller figure
                                          width_in        = 8.0,
                                          height_min      = 8.5,
                                          height_per_cat  = 0.34,
                                          dpi_png         = 900,
                                          legend_breaks_n = 3,
                                          legend_title    = "n") {
  
  # ---- Filter, order and prettify labels
  d <- df %>%
    dplyr::filter(signature %in% signatures,
                  target_type %in% target_types,
                  padj <= alpha,
                  n_present > min_present) %>%
    dplyr::group_by(signature, target_type) %>%
    dplyr::arrange(padj, dplyr::desc(score), .by_group = TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(target_short = shorten(target, 60)) %>%
    dplyr::group_by(signature) %>%
    dplyr::mutate(target_short = factor(target_short,
                                        levels = rev(unique(target_short[order(-score)])))) %>%
    dplyr::ungroup()
  
  if (nrow(d) == 0L) {
    warning("No rows to plot with current filters (padj <= alpha, n_present > min_present).")
    return(invisible(NULL))
  }
  
  # ---- Dynamic figure height by number of categories per facet
  n_cat_per_sig <- d %>%
    dplyr::group_by(signature) %>%
    dplyr::summarise(n = dplyr::n_distinct(target_short), .groups = "drop") %>%
    dplyr::pull(n)
  height_in <- max(height_min, max(n_cat_per_sig, na.rm = TRUE) * height_per_cat)
  
  # ---- Size legend: 3 breaks between min & max n_present
  np_min <- min(d$n_present, na.rm = TRUE)
  np_max <- max(d$n_present, na.rm = TRUE)
  brks <- unique(round(seq(np_min, np_max, length.out = legend_breaks_n)))
  brks <- brks[brks >= np_min & brks <= np_max]
  
  # ---- Compact X with extra left/right breathing room
  xr  <- range(d$score, na.rm = TRUE)
  pad <- diff(xr) * 0.02
  xlim_use <- c(xr[1] - pad, xr[2] + pad)
  
  p <- ggplot(d, aes(x = score, y = target_short)) +
    geom_point(aes(size = n_present, fill = signature),
               shape = 21, color = "grey20", stroke = 0.25, alpha = 0.95) +
    scale_size_continuous(name = legend_title, breaks = brks, range = c(3.2, 10)) +
    scale_fill_viridis_d(option = "D", end = 0.85, name = "Signature") +
    labs(
      title    = NULL,
      subtitle = NULL,
      caption  = "All points are significant (FDR \u2264 0.05)",
      x        = "Enrichment score",
      y        = NULL
    ) +
    facet_grid(signature ~ ., scales = "free_y", space = "free_y") +
    ggthemes::theme_clean(base_size = base_size) +
    theme(
      axis.text.x       = element_text(size = base_size - 1),
      # slightly smaller y labels + right margin so they don't touch the panel
      axis.text.y       = element_text(size = base_size - 3, lineheight = 0.95,
                                       margin = margin(r = 6)),
      plot.caption      = element_text(size = base_size - 2, hjust = 0),
      legend.position   = "right",
      panel.grid.minor  = element_blank(),
      panel.grid.major.y= element_blank(),
      # more separation between HPV-A7 and HPV-A9 panels
      panel.spacing.y   = unit(1.2, "lines"),
      # larger left/bottom margin so labels/caption never overlap
      plot.margin       = margin(t = 14, r = 20, b = 26, l = 30)
    ) +
    # compact X + extra space; generous Y padding so bubbles never clip
    scale_x_continuous(limits = xlim_use, expand = expansion(mult = c(0.06, 0.03))) +
    scale_y_discrete(expand = expansion(add = 1.0)) +
    coord_cartesian(clip = "off")
  
  if (isTRUE(show)) print(p)
  
  # Save (PDF vector + high-DPI PNG)
  ggsave(pdf_out, p, width = width_in, height = height_in, dpi = 300, useDingbats = FALSE)
  ggsave(png_out, p, width = width_in, height = height_in, dpi = dpi_png)
  
  message("Saved: ", pdf_out)
  message("Saved: ", png_out)
  invisible(p)
}

# ---- Call
plot_bubble_horizontal_pretty(
  df           = enrich_sig_n,
  signatures   = c("HPV-A7","HPV-A9"),
  target_types = c("mesh"),
  alpha        = 0.01,
  min_present  = 3
)



# === CHEMBL targets  ===
plot_bubble_horizontal_pretty(
  df           = enrich_sig_n,
  signatures   = c("HPV-A7","HPV-A9"),
  target_types = c("chembl_targets"),   # <-- aquí cambiamos a CHEMBL
  alpha        = 0.01,                  # mismo umbral que usaste arriba
  min_present  = 3,
  pdf_out      = file.path(out_base, "6_4_plot_bubble_chembl.pdf"),
  png_out      = file.path(out_base, "6_4_plot_bubble_chembl.png"),
  show         = TRUE
)
