# =========================================
# Collapse to a single RGES-like row per drug and draw a heatmap (A7 vs A9)
# =========================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# --------------------------
# Paths & config
# --------------------------
base_res   <- "~/CESC_Network/6_OCTAD/6_3_RES"
cutfile    <- "cut-0.20.rds"
signatures <- c("A7","A9")
methods    <- c("DESeq2","EdgeR","limma")

# Set topN to Inf to include ALL drugs (otherwise a number like 80)
topN <- Inf

# --------------------------
# 1) Read & combine per-signature/per-method RDS
# --------------------------
read_one <- function(sig, mth) {
  fp <- file.path(base_res, sig, mth, paste0("RES_", sig, "_", mth, "_", cutfile))
  if (!file.exists(fp)) return(NULL)
  x <- readRDS(fp)
  if (!"signature" %in% names(x)) x$signature <- sig
  if (!"method"    %in% names(x)) x$method    <- mth
  x
}

sig_dfs <- lapply(signatures, function(sig) {
  lst <- lapply(methods, function(m) read_one(sig, m))
  lst <- Filter(Negate(is.null), lst)
  if (length(lst) == 0) return(NULL)
  dplyr::bind_rows(lst)
})
sig_dfs <- Filter(Negate(is.null), sig_dfs)

all_res <- bind_rows(sig_dfs) %>%
  arrange(signature, method, sRGES)

saveRDS(all_res, file.path(base_res, paste0("RES_ALL_", cutfile)))
message("Combined rows: ", nrow(all_res))

# --------------------------
# 2) Keep drugs present in >= 3 methods within each signature
# --------------------------
# Keep drugs present in >= k methods within a given signature
common_by_signature <- function(df, sig, k = 3) {
  df_sig <- dplyr::filter(df, signature == sig)
  
  # cuáles fármacos tienen al menos k métodos distintos
  keep <- df_sig %>%
    dplyr::group_by(pert_iname) %>%
    dplyr::summarise(n_methods = dplyr::n_distinct(method), .groups = "drop") %>%
    dplyr::filter(n_methods >= k) %>%
    dplyr::pull(pert_iname)
  
  # devuelve solamente esos fármacos (ordenados)
  df_sig %>%
    dplyr::filter(pert_iname %in% keep) %>%
    dplyr::arrange(pert_iname, method, sRGES)
}


common_A7 <- common_by_signature(all_res, "A7", k = 3)
common_A9 <- common_by_signature(all_res, "A9", k = 3)

saveRDS(common_A7, file.path(base_res, "A7", paste0("RES_A7_common_3methods_", cutfile)))
saveRDS(common_A9, file.path(base_res, "A9", paste0("RES_A9_common_3methods_", cutfile)))

# --------------------------
# 3) Collapse to one RGES-like row per drug (mean, NA-aware)
# --------------------------
collapse_num <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

collapse_rges <- function(df_sig_common) {
  df_sig_common %>%
    group_by(pert_iname) %>%
    summarise(
      mean   = collapse_num(mean),
      n      = collapse_num(n),        # use round(n) if you prefer integer
      median = collapse_num(median),
      sd     = collapse_num(sd),
      sRGES  = collapse_num(sRGES),
      .groups = "drop"
    ) %>%
    arrange(sRGES)
}

collapsed_A7 <- collapse_rges(common_A7)
collapsed_A9 <- collapse_rges(common_A9)

#saveRDS(collapsed_A7, file.path(base_res, "A7", paste0("RES_A7_common3_collapsed_", cutfile)))
#saveRDS(collapsed_A9, file.path(base_res, "A9", paste0("RES_A9_common3_collapsed_", cutfile)))

# --------------------------
# 4) Build heatmap input (A7 vs A9), membership, ranking
# --------------------------
wide <- full_join(
  collapsed_A7 %>% dplyr::select(pert_iname, A7 = sRGES),
  collapsed_A9 %>% dplyr::select(pert_iname, A9 = sRGES),
  by = "pert_iname"
) %>%
  mutate(
    membership = case_when(
      !is.na(A7) & !is.na(A9) ~ "Common",
      !is.na(A7) &  is.na(A9) ~ "A7-only",
      is.na(A7) & !is.na(A9) ~ "A9-only",
      TRUE ~ NA_character_
    ),
    # Row-wise minimum sRGES across A7/A9 (keep NA if both are NA)
    min_sRGES = ifelse(is.na(A7) & is.na(A9), NA_real_, pmin(A7, A9, na.rm = TRUE))
  )

top_tbl <- wide %>%
  filter(!is.na(min_sRGES)) %>%
  mutate(
    # Show full labels in the legend
    membership = factor(membership,
                        levels = c("Common", "A7-only",   "A9-only"),
                        labels = c("Common", "HPV-A7-only", "HPV-A9-only"))
  ) %>%
  arrange(membership, min_sRGES)

if (is.finite(topN)) {
  top_tbl <- slice(top_tbl, 1:min(topN, nrow(top_tbl)))
}

mat <- top_tbl %>%
  dplyr::select(pert_iname, A7, A9) %>%
  column_to_rownames("pert_iname") %>%
  as.matrix()

# Column labels as requested
colnames(mat) <- c("  HPV-A7", "  HPV-A9")

# --------------------------
# 5) Heatmap settings
#    Color scale anchored at 0: blue (more negative) -> light-yellow (~0) -> red (positive)
# --------------------------
vals    <- as.numeric(mat)
neg_min <- min(vals, na.rm = TRUE)
pos_max <- max(vals, na.rm = TRUE)

neg_col  <- "#2166AC"  # blue for strong negatives
zero_col <- "#FFF3B0"  # light yellow ~0
pos_col  <- "#D7301F"  # red/orange for positives

col_fun <- if (pos_max <= 0 || !is.finite(pos_max)) {
  circlize::colorRamp2(c(neg_min, 0),        c(neg_col, zero_col))
} else {
  circlize::colorRamp2(c(neg_min, 0, pos_max), c(neg_col, zero_col, pos_col))
}

row_ha <- rowAnnotation(
  membership = top_tbl$membership,
  col = list(membership = c("Common"="#16A085", "HPV-A7-only"="#F39C12", "HPV-A9-only"="#8E44AD")),
  annotation_legend_param = list(title = "membership")
)

# añade row_names_max_width para que el layout deje espacio suficiente
ht <- Heatmap(
  mat,
  name = "sRGES",
  col = col_fun,
  na_col = "grey90",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_title = if (is.finite(topN)) paste0("Top-", topN, " drugs by lowest sRGES") else "",
  left_annotation = row_ha,
  border = NA,
  row_names_gp = gpar(fontsize = 12, fontface = "bold"),
  row_names_side = "left",  # <— mueve etiquetas de filas a la IZQUIERDA
  column_names_gp = gpar(fontsize = 16, fontface = "bold"),
  column_names_rot = 0,
  row_names_max_width = max_text_width(rownames(mat), gp = gpar(fontsize = 12, fontface = "bold")) * 1.15
)

draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")


# --------------------------
# 6) Save heatmap as vector graphics (max quality, no size specified)
# ---------
suffix  <- if (is.finite(topN)) paste0("top", topN) else "ALL"
pdf_out <- file.path(base_res, paste0("6_3_3_RES_heatmap_A7_A9_", suffix, ".pdf"))
svg_out <- file.path(base_res, paste0("6_3_3_RES_heatmap_A7_A9_", suffix, ".svg"))
png_out <- file.path(base_res, paste0("6_3_3_RES_heatmap_A7_A9_", suffix, ".png"))

n_rows <- nrow(mat)

# (ajustes suaves) -> un poco más ancho y menos alto
lab_w_in   <- grid::convertWidth(
  max_text_width(rownames(mat), gp = gpar(fontsize = 12, fontface = "bold")) * 1.2,
  "in", valueOnly = TRUE
)
# (ajustes) -> un poco más ancho y menos alto
body_w_in   <- 3.8     # antes 3.6
legend_w_in <- 2.7     # antes 2.4   (más espacio para la leyenda)
pdf_w       <- body_w_in + lab_w_in + legend_w_in

# menos alto aún por fila (y mínimo global más bajo)
pdf_h <- max(5.5, n_rows * 0.165)  # antes: max(6, n_rows * 0.18)

# PDF (vector)
pdf(file = pdf_out, width = pdf_w, height = pdf_h, onefile = TRUE)
draw(ht,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     merge_legend = TRUE,
     padding = unit(c(2, 12, 2, 2), "mm"))  # ↑ más separación a la derecha
dev.off()

# SVG (vector editable)
svg(filename = svg_out, width = pdf_w, height = pdf_h)
draw(ht,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     merge_legend = TRUE,
     padding = unit(c(2, 12, 2, 2), "mm"))
dev.off()

# PNG (1200 dpi)
if (requireNamespace("ragg", quietly = TRUE)) {
  ragg::agg_png(png_out, width = pdf_w, height = pdf_h, units = "in", res = 1200)
} else {
  png(filename = png_out, width = pdf_w, height = pdf_h, units = "in",
      res = 1200, type = "cairo-png")
}
draw(ht,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     merge_legend = TRUE,
     padding = unit(c(2, 12, 2, 2), "mm"))
dev.off()

