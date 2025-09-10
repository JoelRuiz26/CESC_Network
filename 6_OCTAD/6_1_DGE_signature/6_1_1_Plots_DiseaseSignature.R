# =========================================
# Load libraries
# =========================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(patchwork)
  library(ggVennDiagram)
})

# =========================================
# Parameters
# =========================================
lfc_thr  <- 1
padj_thr <- 0.01
base_dir <- "~/CESC_Network/6_OCTAD/6_1_DGE_signature/6_1_0_Output_rds/"
out_dir  <- "~/CESC_Network/6_OCTAD/6_1_DGE_signature/6_1_1_Output_plots/"

# High-contrast palette (consistent across volcano & Venn)
# EdgeR = deep red, DESeq2 = deep blue, limma = teal-green
method_cols <- c(
  EdgeR  = "#B40426",
  DESeq2 = "#3B4CC0",
  limma  = "#1FA187"
)

# =========================================
# Helper to save PNG + PDF (vector)
# =========================================
save_png_pdf <- function(path_png, plot, width, height,
                         units = "in", dpi = 1200, limitsize = FALSE,
                         use_cairo = TRUE) {
  # 1) Save PNG (hi-DPI)
  ggplot2::ggsave(filename = path_png, plot = plot,
                  width = width, height = height, units = units,
                  dpi = dpi, limitsize = limitsize)
  # 2) Save vector PDF (good for journals)
  path_pdf <- sub("\\.png$", ".pdf", path_png)
  if (use_cairo) {
    ggplot2::ggsave(filename = path_pdf, plot = plot,
                    width = width, height = height, units = units,
                    device = grDevices::cairo_pdf, limitsize = limitsize)
  } else {
    ggplot2::ggsave(filename = path_pdf, plot = plot,
                    width = width, height = height, units = units,
                    limitsize = limitsize)
  }
  invisible(path_pdf)
}

# =========================================
# I/O helpers (unchanged logic)
# =========================================

# Read the raw DE table saved previously by your pipeline
read_de_full <- function(clade, method_tag){
  f <- file.path(base_dir, sprintf("6_1_2_DE_full_%s_%s.rds", clade, method_tag))
  if (!file.exists(f)) stop("File not found: ", f)
  readRDS(f)
}

# Choose a stable gene identifier for Venns
pick_id_col <- function(df){
  if ("identifier" %in% names(df)) return("identifier")
  if ("ENSEMBL"   %in% names(df)) return("ENSEMBL")
  NA_character_
}
pick_symbol_col <- function(df){
  if ("Symbol" %in% names(df)) return("Symbol")
  if ("Symbol_autho" %in% names(df)) return("Symbol_autho")
  NA_character_
}

# Minimal frame for Venn sets
read_de_for_venn <- function(clade, method_tag){
  x <- read_de_full(clade, method_tag)
  id_col  <- pick_id_col(x)
  sym_col <- pick_symbol_col(x)
  if (!is.na(sym_col)) {
    x <- x %>% dplyr::mutate(gene_id = toupper(.data[[sym_col]]))
  } else if (!is.na(id_col)) {
    x <- x %>% dplyr::mutate(gene_id = .data[[id_col]])
  } else {
    stop("No gene identifier column found in DE result.")
  }
  stopifnot(all(c("gene_id","log2FoldChange","padj") %in% names(x)))
  x %>%
    dplyr::select(gene_id, log2FoldChange, padj) %>%
    dplyr::filter(is.finite(log2FoldChange), is.finite(padj)) %>%
    dplyr::mutate(method = method_tag)
}

# Minimal frame for volcano plotting
read_de_for_volcano <- function(clade, method_tag){
  read_de_for_venn(clade, method_tag) %>%     
    dplyr::distinct(gene_id, .keep_all = TRUE) %>% # keep one row per gene
    dplyr::select(log2FoldChange, padj) %>%
    dplyr::mutate(method = method_tag)
}

# =========================================
# Volcano prep (unchanged logic)
# =========================================
prep_volcano_df <- function(clade){
  methods <- c("EdgeR","DESeq2","limma")
  df <- dplyr::bind_rows(lapply(methods, \(m) read_de_for_volcano(clade, m))) %>%
    dplyr::mutate(
      method   = factor(method, levels = names(method_cols)),
      mlog10p  = -log10(padj),
      sig      = dplyr::case_when(
        log2FoldChange >=  lfc_thr & padj < padj_thr ~ "up",
        log2FoldChange <= -lfc_thr & padj < padj_thr ~ "down",
        TRUE ~ "ns"
      ),
      # Zones to control transparency explicitly
      zone = dplyr::case_when(
        abs(log2FoldChange) <  lfc_thr & padj >= padj_thr ~ "central_ns", 
        abs(log2FoldChange) <  lfc_thr & padj <  padj_thr ~ "lfc_only_ns", 
        abs(log2FoldChange) >= lfc_thr & padj >= padj_thr ~ "p_only_ns",   
        TRUE ~ "sig"
      ),
      alpha_plot = dplyr::case_when(
        zone == "central_ns" ~ 0.02,  
        zone == "lfc_only_ns" ~ 0.06,
        zone == "p_only_ns"   ~ 0.06,
        TRUE                  ~ 0.92  
      )
    )
  df
}

# Count up/down by method to annotate inside the plot
annot_counts_df <- function(df){
  counts <- df %>%
    dplyr::filter(sig != "ns") %>%
    dplyr::count(method, sig, name = "n") %>%
    tidyr::complete(method, sig = c("up","down"), fill = list(n = 0))
  
  ymax <- max(df$mlog10p, na.rm = TRUE)
  xmax <- max(abs(df$log2FoldChange), na.rm = TRUE)
  
  methods <- levels(df$method)
  vstep <- (ymax * 0.05)
  base_y <- ymax * 0.98
  idx <- setNames(seq_along(methods)-1, methods)
  
  counts %>%
    dplyr::mutate(
      x = dplyr::if_else(sig == "up",  xmax*0.985, -xmax*0.985),
      y = base_y - idx[as.character(method)]*vstep,
      hjust = dplyr::if_else(sig == "up", 1, 0),
      # Formal label
      label = paste0(method, " ", ifelse(sig=="up","up","down"), ": ", n)
    )
}

# =========================================
# Volcano overlay (per clade) — identical look
# =========================================
make_overlay_volcano <- function(clade, out_png = NULL){
  df   <- prep_volcano_df(clade)
  ann  <- annot_counts_df(df)
  ythr <- -log10(padj_thr)
  
  xmax   <- max(abs(df$log2FoldChange), na.rm = TRUE)
  x_lims <- c(-xmax, xmax) * 1.03
  y_max  <- max(df$mlog10p, na.rm = TRUE) * 1.02
  
  df_ns  <- dplyr::filter(df, zone != "sig")  # draw beneath
  df_sig <- dplyr::filter(df, zone == "sig")  # draw on top
  
  p <- ggplot() +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0,    ymax = ythr,
             fill = "grey50", alpha = 0.05) +
    annotate("rect", xmin = -lfc_thr, xmax = lfc_thr, ymin = 0, ymax = Inf,
             fill = "grey50", alpha = 0.05) +
    geom_point(data = df_ns,
               aes(x = log2FoldChange, y = mlog10p, color = method, alpha = alpha_plot),
               size = 0.8, stroke = 0, show.legend = FALSE) +
    geom_point(data = df_sig,
               aes(x = log2FoldChange, y = mlog10p, color = method),
               size = 1.1, alpha = 0.92, stroke = 0, show.legend = TRUE) +
    scale_color_manual(values = method_cols, name = "Method") +
    scale_alpha_identity() +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed", linewidth = 0.5) +
    geom_hline(yintercept = ythr,                    linetype = "dashed", linewidth = 0.5) +
    geom_text(data = ann,
              aes(x = x, y = y, label = label, color = method, hjust = hjust),
              inherit.aes = FALSE, size = 4.4, fontface = "bold") +
    scale_y_continuous(limits = c(0, y_max), expand = expansion(mult = c(0, 0))) +
    scale_x_continuous(limits = x_lims,     expand = expansion(mult = c(0, 0.01))) +
    labs(
      title = paste0("DGE (HPV-", clade, " clade)"),
      x = expression(log[2]~Fold~Change),
      y = expression(-log[10](adjusted~p))
    ) +
    theme_bw(base_size = 20) +
    theme(
      legend.position   = "right",
      legend.title      = element_text(size = 20),
      legend.text       = element_text(size = 18),
      axis.title        = element_text(size = 19),
      axis.text         = element_text(size = 18),
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(linewidth = 0.25, color = "white"),
      plot.title        = element_text(face = "bold", size = 18)
    ) +
    guides(color = guide_legend(override.aes = list(size = 3.5, alpha = 1)))
  
  if (is.null(out_png)) {
    out_png <- file.path(out_dir, paste0("6_1_6_Volcano_HPV_", clade, "_clade_overlay.png"))
  }
  # >>> Save PNG + PDF (instead of PNG only) <<<
  save_png_pdf(out_png, p, width = 7.6, height = 4.9, dpi = 1200)
  p
}

# =========================================
# Build volcanos (unchanged)
# =========================================
pA7_ov <- make_overlay_volcano("A7")
pA9_ov <- make_overlay_volcano("A9")

# Grid with shared legend and larger text
g_volcano <- (pA7_ov + pA9_ov) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.title    = element_text(size = 20),
    legend.text     = element_text(size = 18)
  )

# >>> Save PNG + PDF for the grid <<<
save_png_pdf(file.path(out_dir, "6_1_7_Volcano_HPV_A7_A9_clades_overlay_grid.png"),
             g_volcano, width = 13.5, height = 5.6, dpi = 1200)

# =========================================
# Venns (Up/Down by clade) — identical logic
# =========================================
get_sig_sets <- function(clade){
  methods <- c("EdgeR","DESeq2","limma")
  dfs <- lapply(methods, \(m) read_de_for_venn(clade, m))
  names(dfs) <- methods
  up   <- lapply(dfs, \(d) unique(d$gene_id[d$log2FoldChange >=  lfc_thr & d$padj < padj_thr]))
  down <- lapply(dfs, \(d) unique(d$gene_id[d$log2FoldChange <= -lfc_thr & d$padj < padj_thr]))
  list(up = up, down = down)
}

plot_venn <- function(sets_list, title, out_png,
                      count_size = 4.5, setlabel_size = 4.2,
                      width = 8, height = 8, dpi = 1200,
                      fill_tone = "#F7FAFF") {
  # Ensure order matches your palette
  stopifnot(exists("method_cols"))
  sets_list <- sets_list[names(method_cols)]
  
  # Base Venn
  p <- ggVennDiagram(
    sets_list,
    label = "count",
    set_color      = unname(method_cols[names(sets_list)]),
    set_edge_color = unname(method_cols[names(sets_list)]),
    set_edge_size  = 1.2
  ) +
    # Flat very-light fill; no legend
    scale_fill_gradient(low = fill_tone, high = fill_tone, guide = "none") +
    labs(title = title) +
    theme_void(base_size = 16) +
    theme(
      plot.title      = element_text(face = "bold", hjust = 0.5, size = 16, margin = margin(b = 6)),
      text            = element_text(size = 16),
      legend.position = "none",
      plot.margin     = margin(t = 6, r = 30, b = 6, l = 30)
    ) +
    coord_cartesian(clip = "off") +
    theme(aspect.ratio = 1)
  
  # Tweak text layers:
  # - make counts bigger & bold
  # - make set labels bigger & bold
  # - REMOVE label rectangles (no fill, no border)
  layers <- p$layers
  for (i in seq_along(layers)) {
    dat <- layers[[i]]$data
    g   <- layers[[i]]$geom
    if (inherits(g, "GeomText") || inherits(g, "GeomLabel") ||
        inherits(g, "GeomSfText") || inherits(g, "GeomSfLabel")) {
      if (is.null(layers[[i]]$aes_params)) layers[[i]]$aes_params <- list()
      
      # Interior counts have a "count" column in the layer data
      if (is.data.frame(dat) && ("count" %in% names(dat))) {
        layers[[i]]$aes_params$size     <- count_size
        layers[[i]]$aes_params$fontface <- "bold"
      } else {
        layers[[i]]$aes_params$size     <- setlabel_size
        layers[[i]]$aes_params$fontface <- "bold"
      }
      
      # Kill the label box if this layer is a label-type geom
      if (inherits(g, "GeomLabel") || inherits(g, "GeomSfLabel")) {
        layers[[i]]$aes_params$fill       <- NA  # transparent background
        layers[[i]]$aes_params$label.size <- 0   # no rectangle border
        # Optional: tighten padding
        # layers[[i]]$aes_params$label.padding <- grid::unit(0.05, "lines")
      }
    }
  }
  p$layers <- layers
  
  # Save PNG + vector PDF
  ggplot2::ggsave(out_png, p, width = width, height = height,
                  dpi = dpi, units = "in", limitsize = FALSE)
  ggplot2::ggsave(sub("\\.png$", ".pdf", out_png), p,
                  width = width, height = height, units = "in",
                  device = grDevices::cairo_pdf)
  p
}


# A7
sets_A7 <- get_sig_sets("A7")
pA7_venn_down <- plot_venn(
  sets_A7$down,
  "HPV-A7 clade \
  Down-regulated genes",
  file.path(out_dir, "6_1_6_Venn_HPV_A7_clade_DOWN.png")
)

pA7_venn_up <- plot_venn(
  sets_A7$up,
  "HPV-A7 clade \
  Up-regulated genes",
  file.path(out_dir, "6_1_6_Venn_HPV_A7_clade_UP.png")
)

# A9
sets_A9 <- get_sig_sets("A9")
pA9_venn_down <- plot_venn(
  sets_A9$down,
  "HPV-A9 clade \
  Down-regulated genes",
  file.path(out_dir, "6_1_6_Venn_HPV_A9_clade_DOWN.png")
)

pA9_venn_up <- plot_venn(
  sets_A9$up,
  "HPV-A9 clade \
  Up-regulated genes",
  file.path(out_dir, "6_1_6_Venn_HPV_A9_clade_UP.png")
)

# Optional grids (per clade) with larger canvas
gA7_venn <- pA7_venn_down + pA7_venn_up + plot_layout(ncol = 2)
save_png_pdf(file.path(out_dir, "6_1_8_Venn_HPV_A7_clade_grid.png"),
             gA7_venn, width = 12.0, height = 5.2, dpi = 1200)

gA9_venn <- pA9_venn_down + pA9_venn_up + plot_layout(ncol = 2)
save_png_pdf(file.path(out_dir, "6_1_8_Venn_HPV_A9_clade_grid.png"),
             gA9_venn, width = 12.0, height = 5.2, dpi = 1200)

gA7_venn
gA9_venn

# =========================================
# Combo figure: volcanos (top) + (Venn A7 | Venn A9) (bottom)
# =========================================
combo_fig <- g_volcano / (gA7_venn | gA9_venn) +
  plot_layout(heights = c(1, 1)) &
  theme(plot.title = element_text(margin = margin(b = 10)))

combo_fig

save_png_pdf(file.path(out_dir, "6_1_9_Volcano_Venns_combo.png"),
             combo_fig, width = 12, height = 10, dpi = 1200, limitsize = FALSE)

# =========================================
# Save intersection genes (unchanged logic)
# =========================================
sets_A7 <- get_sig_sets("A7")
sets_A9 <- get_sig_sets("A9")

cons_up_A7   <- Reduce(intersect, sets_A7$up)
cons_down_A7 <- Reduce(intersect, sets_A7$down)
cons_up_A9   <- Reduce(intersect, sets_A9$up)
cons_down_A9 <- Reduce(intersect, sets_A9$down)

consensus_tbl <- dplyr::bind_rows(
  tibble::tibble(clade = "A7", regulation = "up",   gene = sort(cons_up_A7)),
  tibble::tibble(clade = "A7", regulation = "down", gene = sort(cons_down_A7)),
  tibble::tibble(clade = "A9", regulation = "up",   gene = sort(cons_up_A9)),
  tibble::tibble(clade = "A9", regulation = "down", gene = sort(cons_down_A9))
)

saveRDS(consensus_tbl, file.path(out_dir, "6_10_Consensus_genes_3methods.rds"))
table(consensus_tbl$clade, consensus_tbl$regulation)

# =========================================
# A7 vs A9 Venns for consensus sets (UP/DOWN)
# =========================================

# Helpers
cons_sets_by <- function(tbl, reg){
  a7 <- tbl %>% dplyr::filter(clade=="A7", regulation==reg) %>% dplyr::pull(gene) %>% unique()
  a9 <- tbl %>% dplyr::filter(clade=="A9", regulation==reg) %>% dplyr::pull(gene) %>% unique()
  list(A7 = a7, A9 = a9)
}

shared_stats <- function(a, b){
  inter <- intersect(a, b)
  uni   <- union(a, b)
  c(nA   = length(a),
    nB   = length(b),
    nI   = length(inter),
    pctA = 100 * length(inter)/max(1, length(a)),
    pctB = 100 * length(inter)/max(1, length(b)),
    jacc = 100 * length(inter)/max(1, length(uni)))
}

plot_clade_venn <- function(sets, title_prefix, out_png,
                            count_size = 7, setlabel_size = 5,
                            width = 6.5, height = 6.5, dpi = 1200,
                            fill_tone = "#F7FAFF") {
  stats <- shared_stats(sets$A7, sets$A9)
  subtitle <- sprintf("Shared = %d | %.1f%% of A7, %.1f%% of A9 | Jaccard = %.1f%%",
                      stats["nI"], stats["pctA"], stats["pctB"], stats["jacc"])
  
  sets_plot <- list(`HPV-A7` = sets$A7, `HPV-A9` = sets$A9)
  
  p <- ggVennDiagram(
    sets_plot, label = "count",
    set_color      = c("#6A3D9A", "#FF7F00"),
    set_edge_color = c("#6A3D9A", "#FF7F00"),
    set_edge_size  = 1.2
  ) +
    scale_fill_gradient(low = fill_tone, high = fill_tone, guide = "none") +
    labs(title = paste0(title_prefix, " (consensus 3 methods)"),
         subtitle = subtitle) +
    theme_void(base_size = 16) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 16, margin = margin(b = 2)),
      plot.subtitle = element_text(hjust = 0.5, size = 13),
      legend.position = "none",
      plot.margin   = margin(t = 6, r = 26, b = 6, l = 26)
    ) +
    coord_cartesian(clip = "off") +
    theme(aspect.ratio = 1) 
  
  # Enlarge counts and labels
  layers <- p$layers
  for (i in seq_along(layers)) {
    dat <- layers[[i]]$data
    cls <- class(layers[[i]]$geom)
    if (any(cls %in% c("GeomText", "GeomLabel", "GeomSfText", "GeomSfLabel"))) {
      if (is.null(layers[[i]]$aes_params)) layers[[i]]$aes_params <- list()
      if (is.data.frame(dat) && ("count" %in% names(dat))) {
        layers[[i]]$aes_params$size <- count_size
        layers[[i]]$aes_params$fontface <- "bold"
      } else {
        layers[[i]]$aes_params$size <- setlabel_size
        layers[[i]]$aes_params$fontface <- "bold"
      }
    }
  }
  p$layers <- layers
  
  # Save PNG + PDF (vector)
  ggsave(out_png, p, width = width, height = height, dpi = dpi, units = "in", limitsize = FALSE)
  ggsave(sub("\\.png$", ".pdf", out_png), p, width = width, height = height, units = "in",
         device = grDevices::cairo_pdf)
  p
}

# Build A7 vs A9 consensus Venns (UP and DOWN)
sets_cons_up   <- cons_sets_by(consensus_tbl, "up")
sets_cons_down <- cons_sets_by(consensus_tbl, "down")

p_clades_up <- plot_clade_venn(
  sets_cons_up,
  "UP-regulated genes (HPV-A7 vs HPV-A9)",
  file.path(out_dir, "6_1_10_Venn_A7_vs_A9_UP_consensus.png")
)

p_clades_down <- plot_clade_venn(
  sets_cons_down,
  "DOWN-regulated genes (HPV-A7 vs HPV-A9)",
  file.path(out_dir, "6_1_10_Venn_A7_vs_A9_DOWN_consensus.png")
)

# Grid for these two consensus Venns
g_clade_venns <- p_clades_up + p_clades_down + plot_layout(ncol = 2)
save_png_pdf(file.path(out_dir, "6_1_10_Venn_A7_vs_A9_consensus_grid.png"),
             g_clade_venns, width = 13.0, height = 6.8, dpi = 1200)

g_clade_venns

# Optional: save workspace (unchanged)
save.image("~/CESC_Network/6_OCTAD/6_1_DGE_signature/6_1_1_Output_plots/6_1_9_Image_plots.RData")
#load("~/CESC_Network/6_OCTAD/6_1_DGE_signature/6_1_1_Output_plots/6_1_9_Image_plots.RData")
