###############################################################################
# Heatmap: Top 10 pathways by tumor median PDS z-score (paper-style)
# Logic: row-wise z over ALL samples; rank by tumor median z; show top-10
###############################################################################

options(stringsAsFactors = FALSE)
setwd("~/Pathifier-cervical-cancer/")

suppressPackageStartupMessages({
  if (!requireNamespace("gplots", quietly = TRUE))       install.packages("gplots")
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer")
  library(gplots); library(RColorBrewer)
})

K     <- 10L   # top-k pathways
CLIP  <- 3.0   # clip z to [-CLIP, +CLIP] for color scale
OUT   <- "Output_plots/heatmap_top10_tumorMedianZ_allZ.png"

# --- Load PDS and normals -----------------------------------------------------
load("1_1_PDS.RData")
stopifnot(exists("PDS"))
PDSmatrix <- t(mapply(FUN = c, PDS$scores))
if (is.null(rownames(PDSmatrix))) rownames(PDSmatrix) <- names(PDS$scores)

if (!exists("normals") || length(normals) != ncol(PDSmatrix)) {
  stopifnot(exists("metadata"))
  smap <- setNames(metadata$sample_type, metadata$specimenID)
  common <- intersect(colnames(PDSmatrix), names(smap))
  PDSmatrix <- PDSmatrix[, common, drop = FALSE]
  normals <- smap[common] == "Solid Tissue Normal"
}
stopifnot(any(normals), any(!normals))

# --- 1) z-score per pathway ACROSS ALL samples --------------------------------
zAll <- t(scale(t(PDSmatrix)))             # mean=0, sd=1 per row over ALL columns
zAll[!is.finite(zAll)] <- NA

# --- 2) rank by tumor median z (directional) ----------------------------------
tum_med <- apply(zAll[, !normals, drop = FALSE], 1, median, na.rm = TRUE)
ord     <- order(tum_med, decreasing = TRUE, na.last = NA)
keep    <- ord[seq_len(min(K, length(ord)))]
zTop    <- zAll[keep, , drop = FALSE]

# --- 3) reorder samples: normals first (like the paper visual) ----------------
o_cols  <- c(which(normals), which(!normals))
zTop    <- zTop[, o_cols, drop = FALSE]
normals <- normals[o_cols]

# --- 4) clip for color contrast & cluster ------------------------------------
zTop[zTop >  CLIP] <-  CLIP
zTop[zTop < -CLIP] <- -CLIP

row.distance <- if (nrow(zTop) > 1) dist(zTop, method = "euclidean") else NULL
row.cluster  <- if (!is.null(row.distance)) hclust(row.distance, method = "ward.D2") else NULL
col.distance <- if (ncol(zTop) > 1) dist(t(zTop), method = "euclidean") else NULL
col.cluster  <- if (!is.null(col.distance)) hclust(col.distance, method = "ward.D2") else NULL

# --- palette, breaks, column colour key ---------------------------------------
my_palette <- rev(colorRampPalette(brewer.pal(11, "Spectral"))(1000))
my_breaks  <- seq(-CLIP, CLIP, length.out = 1001)
colLabels  <- ifelse(normals, "#377EB8", "#E41A1C")   # blue=Normal, red=Tumor

# --- 5) draw ------------------------------------------------------------------
png(OUT, width = 3400, height = 2500, units = "px", res = 300, pointsize = 12)
heatmap.2(zTop,
          main = "Top 10 deregulated pathways by tumor median PDS z-score",
          density.info = "none", trace = "none",
          margins = c(6, 22),      # más espacio en eje Y para etiquetas
          col = my_palette, breaks = my_breaks,
          Rowv = if (!is.null(row.cluster)) as.dendrogram(row.cluster) else FALSE,
          Colv = if (!is.null(col.cluster)) as.dendrogram(col.cluster) else FALSE,
          keysize = 1.0,
          ColSideColors = colLabels,
          na.color = "grey95",
          cexRow = 1.0,            # aumentar letra de las vías
          cexCol = 0.6             # mantener columnas pequeñas
)

par(lend = 1, xpd = NA)
legend("topright", legend = c("Normals","Tumors"),
       col = c("#377EB8","#E41A1C"), lty = 1, lwd = 5, bty = "n", inset = 0.02)
dev.off()
message("Saved: ", OUT)
