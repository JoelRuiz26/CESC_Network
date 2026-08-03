### ============================================================
### 6_9_C_Selective_hub_heatmap.R
###
### Per-cell-line CRISPR gene effect of the non pan-essential essential
### hubs, i.e. the same 90 genes coloured in 6_9_B.
###   columns = the 10 cervical lines, split manually by HPV clade
###   rows    = the 90 genes, split by the network(s) in which the gene
###             is an essential hub (HPV-A7 / both / HPV-A9)
###
### Thresholds are the data-driven ones of 6_10 (elbow of each ranked
### curve): A7 degree >= 35 and gene effect <= -0.317; A9 degree >= 42
### and gene effect <= -0.316.
###
### Columns are ordered and split by clade, never clustered: with 10
### lines an unsupervised dendrogram carries no evidential weight, and
### the genes were selected using clade mean gene effect, so clustering
### the same lines on the same genes would be circular.
###
###   6_9_C_selective_hub_heatmap.png/.pdf
### ============================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(vroom)
  library(depmap)
  library(ComplexHeatmap)
  library(circlize)
})

res_dir <- "~/CESC_Network/6_DepMap_cacu/results"
FF  <- "Times New Roman"
A7c <- "#4E79A7"
A9c <- "#CB6A5C"

## ── Cell lines, pan-essential reference ───────────────────────────────────────
lines <- vroom(file.path(res_dir, "6_1_crispr_lines.tsv"), show_col_types = FALSE) %>%
  mutate(label = str_remove(cell_line, "_CERVIX$"))

pan <- vroom(file.path(res_dir, "6_1_common_essential_genes.tsv"),
             delim = "\t", show_col_types = FALSE)$gene_u %>% toupper()

## ── The selective set, rebuilt with the 6_10 thresholds ──────────────────────
elbow <- function(v) {
  n <- length(v)
  x <- (seq_len(n) - 1) / (n - 1)
  y <- (v - v[1]) / (v[n] - v[1])
  v[which.max(abs(y - x))]
}

net <- vroom(file.path(res_dir, "6_2_gene_dep_degree_table.tsv"), show_col_types = FALSE) %>%
  filter(!is.na(degree), !is.na(mean_dep)) %>%
  distinct(network, gene, .keep_all = TRUE)

th <- net %>%
  group_by(network) %>%
  summarise(k_cut = elbow(sort(degree, decreasing = TRUE)),
            d_cut = elbow(sort(mean_dep)), .groups = "drop")

sel <- net %>%
  left_join(th, by = "network") %>%
  mutate(ess = degree >= k_cut & mean_dep <= d_cut,
         net = str_remove(network, "HPV-")) %>%
  filter(ess) %>%
  select(gene, net) %>%
  distinct() %>%
  pivot_wider(names_from = net, values_from = net,
              values_fn = ~ TRUE, values_fill = FALSE) %>%
  mutate(gene = toupper(gene),
         block = case_when(A7 & A9 ~ "Both", A7 ~ "HPV-A7", TRUE ~ "HPV-A9")) %>%
  filter(!(gene %in% pan)) %>%                      # non pan-essential only
  select(gene, block)

message("selective essential hubs: ", nrow(sel))
print(table(sel$block))

## ── Per-line gene effect matrix ──────────────────────────────────────────────
cr <- depmap_crispr() %>%
  filter(depmap_id %in% lines$depmap_id, toupper(gene_name) %in% sel$gene) %>%
  transmute(gene = toupper(gene_name), depmap_id, effect = dependency)

m <- cr %>%
  left_join(lines %>% select(depmap_id, label), by = "depmap_id") %>%
  select(gene, label, effect) %>%
  pivot_wider(names_from = label, values_from = effect) %>%
  column_to_rownames("gene") %>%
  as.matrix()

## ── Splits and ordering ───────────────────────────────────────────────────────
blk <- setNames(sel$block, sel$gene)[rownames(m)]
row_split <- factor(blk, levels = c("HPV-A7", "Both", "HPV-A9"))

clade <- setNames(lines$hpv_group, lines$label)[colnames(m)]
col_split <- factor(clade, levels = c("A7", "A9"))

## strongest first inside each block and each clade
row_ord <- order(row_split, rowMeans(m, na.rm = TRUE))
col_ord <- order(col_split, colMeans(m, na.rm = TRUE))
m         <- m[row_ord, col_ord]
row_split <- row_split[row_ord]
col_split <- col_split[col_ord]

## ── Colours: a single-hue green ramp, deliberately sharing no hue with the
##    clade annotation so the two legends cannot be confused ───────────────────
col_fun <- colorRamp2(c(0.2, 0, -0.5, -1.0, -1.8),
                      c("#F5F2E9", "#FAFAF7", "#C7E6C4", "#4EA96B", "#0B4C2C"))

top_ann <- HeatmapAnnotation(
  Clade = clade[colnames(m)],
  col   = list(Clade = c(A7 = A7c, A9 = A9c)),
  simple_anno_size = unit(4, "mm"),
  annotation_name_gp = gpar(fontfamily = FF, fontsize = 11, fontface = "bold"),
  annotation_legend_param = list(
    Clade = list(title = "HPV clade",
                 title_gp  = gpar(fontfamily = FF, fontsize = 12, fontface = "bold"),
                 labels_gp = gpar(fontfamily = FF, fontsize = 11))),
  border = TRUE
)

ht <- Heatmap(
  m,
  name = "effect",
  col  = col_fun,
  na_col = "grey92",

  cluster_rows = FALSE, cluster_columns = FALSE,
  row_split = row_split, column_split = col_split,
  row_gap = unit(3, "mm"), column_gap = unit(4, "mm"),
  row_title_rot = 0, row_title_side = "right",
  row_title_gp    = gpar(fontfamily = FF, fontsize = 13, fontface = "bold"),
  column_title_gp = gpar(fontfamily = FF, fontsize = 15, fontface = "bold"),

  top_annotation = top_ann,
  rect_gp = gpar(col = "white", lwd = 0.6),
  border = TRUE,

  row_names_side = "left",
  row_names_gp = gpar(fontfamily = FF, fontsize = 8.5, fontface = "bold"),
  column_names_gp = gpar(fontfamily = FF, fontsize = 12, fontface = "bold"),
  column_names_rot = 45,

  heatmap_legend_param = list(
    title = "CRISPR\ngene effect",
    at = c(0, -0.5, -1, -1.5),
    labels = c("0", "-0.5", "-1.0", "-1.5"),
    legend_height = unit(38, "mm"),
    title_gp  = gpar(fontfamily = FF, fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontfamily = FF, fontsize = 11),
    border = "grey40")
)

draw_it <- function() {
  draw(ht,
       column_title = "Per-line gene effect of the non pan-essential essential hubs",
       column_title_gp = gpar(fontfamily = FF, fontsize = 17, fontface = "bold"),
       merge_legend = TRUE,
       padding = unit(c(4, 4, 4, 4), "mm"))
}

H <- 3.2 + nrow(m) * 0.155                      # scale height with the gene count
png(file.path(res_dir, "6_9_C_selective_hub_heatmap.png"),
    width = 9.5, height = H, units = "in", res = 300, type = "cairo", bg = "white")
draw_it(); dev.off()

cairo_pdf(file.path(res_dir, "6_9_C_selective_hub_heatmap.pdf"),
          width = 9.5, height = H)
draw_it(); dev.off()

message("✅ genes = ", nrow(m), " | lines = ", ncol(m), " | saved to: ", res_dir)
print(table(row_split))
print(table(col_split))
