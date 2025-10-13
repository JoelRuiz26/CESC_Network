###############################################################################
## Heatmap for Pathifier results in R
### Author: Joel Ruiz. Taked and modified from: 
#           Angel Garcia-Campos https://github.com/AngelCampos
###############################################################################

###############################################################################
### Installing and/or loading required packages
###############################################################################
setwd("~/Pathifier-cervical-cancer/")
suppressPackageStartupMessages({
  if (!require("gplots"))         { install.packages("gplots", dependencies = TRUE);         library(gplots) }
  if (!require("RColorBrewer"))   { install.packages("RColorBrewer", dependencies = TRUE);   library(RColorBrewer) }
})

###############################################################################
## Load Pathifier results and turn into a matrix
###############################################################################
load("1_1_PDS.RData")
PDSmatrix <- t(mapply(FUN = c, PDS$scores))

###############################################################################
## Creating Custom Palette
###############################################################################

# creates a own color palette passing from blue, green yellow to dark red
my_palette <- rev(colorRampPalette(brewer.pal(11, "Spectral"))(n = 1000))

###############################################################################
## Clustering Methods
###############################################################################

# If you want to change the default clustering method (complete linkage method
# with Euclidean distance measure), this can be done as follows: For non-square
# matrix, we can define the distance and cluster based on our matrix data by

row.distance = dist(PDSmatrix, method = "euclidean")
row.cluster = hclust(row.distance, method = "ward.D2")

col.distance = dist(t(PDSmatrix), method = "euclidean")
col.cluster = hclust(col.distance, method = "ward.D2")

# Arguments for the dist() function are: euclidean (default), maximum, canberra,
# binary, minkowski, manhattan
# And arguments for hclust(): complete (default), single, average, mcquitty,
# median, centroid, ward.D2

###############################################################################
## Assign Column labels (Optional)
###############################################################################

colLabels <- as.character(normals)
colLabels[colLabels == "TRUE"] <- "#377EB8"
colLabels[colLabels == "FALSE"] <- "#E41A1C"

###############################################################################
## Plotting the Heatmap!! (where all colorful things happen...)
###############################################################################
if (is.null(rownames(PDSmatrix))) rownames(PDSmatrix) <- names(PDS$scores)
rownames(PDSmatrix) <- sub("^KEGG_MEDICUS_", "", rownames(PDSmatrix))


png("Output_plots/1_3_heatmap.png",
    width = 3400, height = 2500, units = "px",
    res = 300,            # 300 ppi (más estable que 1200 para márgenes/texto)
    pointsize = 12)       # un poco mayor para legibilidad

hm <- heatmap.2(PDSmatrix,
                main = "PDS-Heatmap",
                density.info = "none",
                trace = "none",
                margins = c(6, 30),      # evita "figure margins too large"
                col = my_palette,
                Rowv = as.dendrogram(row.cluster),
                Colv = as.dendrogram(col.cluster),
                keysize = 1.0,
                ColSideColors = colLabels,
                na.color = "grey95")     # por si hay NA

## Leyenda: debe ir ANTES de dev.off()
par(lend = 1, xpd = NA)            # xpd=NA deja que la leyenda sobresalga si hace falta
legend("topright",
       legend = c("Normals", "Tumors"),
       col = c("dodgerblue", "firebrick1"),
       lty = 1, lwd = 5, bty = "n", inset = 0.02)

dev.off()

#graphics.off()
#heatmap.2(PDSmatrix,
#          main="PDS-Heatmap", density.info="none", trace="none",
#          margins=c(6,10), col=my_palette,
#          Rowv=as.dendrogram(row.cluster),
#          Colv=as.dendrogram(col.cluster),
#          keysize=1.0, ColSideColors=colLabels, na.color="grey95")





