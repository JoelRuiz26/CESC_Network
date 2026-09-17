### ============================================================
### 5_2_10_Figure2_module_overlap_hallmark.R
###
### Figure 2 of the manuscript, one single composite figure:
###
###   a) Module-membership overlap between the HPV-A7 and HPV-A9 Louvain
###      partitions (Jaccard index). The full 16 x 28 matrix is computed and
###      clustered, then the modules that show no enrichment in panel b are
###      dropped: they carry no biology in this figure and only make the panel
###      sparse. A module belonging to a pair the text discusses is kept even
###      when it has no enrichment, so nothing the manuscript names disappears.
###      The near-zero background that carries the divergence message is still
###      visible, and every cell above J_SHOW carries its value, whether or not
###      the text names that pair. The panel is deliberately free of any biological coding: at this
###      point of the Results only membership overlap has been introduced.
###      Every cell has the same size, as in 5_1_3B, and nothing is outlined:
###      the printed value is what separates an overlap from the background.
###
###      The matrix is the one of 5_1_Modules/5_1_3B_jaccard_heatmap (verified
###      identical cell by cell, max difference 0) and it supersedes it, with
###      the module-size bars of that figure kept. Rows and columns follow its
###      dendrogram order, computed on the complete matrix before trimming,
###      reproduced here from the matrix itself with pheatmap's defaults
###      (euclidean distance, complete linkage, modules entered largest first),
###      so both are recognisable as the same matrix. The trees themselves are
###      not drawn: the order is what makes them comparable, not the trees.
###
###   b) Hallmark over-representation per module (A7 | A9). Every significant
###      enrichment is drawn. Rows are blocked by the MSigDB process category of
###      each gene set (Liberzon et al. 2015, Cell Syst 1:417-425, Table 1 - see
###      5_2_10_hallmark_process_categories.tsv), NOT by any grouping of ours.
###      Every gene set and every module is drawn at the same size: the panel
###      shows the enrichment landscape, the text says which parts of it matter.
###
### Nothing about the data is typed in: gene-set names are the MSigDB IDs from
### the ORA output (HALLMARK_ dropped, "_" -> " "), Jaccard values are read from
### the matrix computed here and the row blocks come from the reference table.
### The only hand-written elements are the editorial selection of which module
### pairs the manuscript spotlights, the axis titles and the legend titles.
###
### Inputs:
###   06_modules/5_1_1_Louvain_A{7,9}_graph_max.rds
###   06_modules/5_2_3_Hallmark{7,9}_louv_max.tsv
###   06_modules/5_2_10_hallmark_process_categories.tsv
###   -> 5_2_10_Figure2_module_overlap_hallmark.png / .pdf
### ============================================================
suppressPackageStartupMessages({
  library(tidyverse); library(vroom); library(igraph)
  library(ggh4x);    library(patchwork); library(viridis); library(ggnewscale)
})

mod_dir <- "~/CESC_Network/06_modules/results"
ora_dir <- "~/CESC_Network/06_modules/results"
FAM     <- "Nimbus Roman"          # Times-metric face actually present here
J_SHOW  <- 0.1                     # a cell is worth a printed value above this

## ---- editorial selection (manuscript text) ---------------------------------
## Only identifiers are written here; every value and every grouping below is
## derived from the data or from the MSigDB reference table.
focus_pair <- tribble(
  ~A7,   ~A9,
  "M6",  "M3",
  "M1",  "M1",
  "M2",  "M2",
  "M8",  "M13",
  "M10", "M19",
  "M12", "M17",
  "M15", "M23")

num <- function(x) as.integer(sub("M", "", x))   # module number, used by both panels

## ====================== enrichment table (both panels) ======================
## read before panel a: which modules carry an enrichment decides which rows
## and columns of the Jaccard matrix are worth drawing
gr_num <- function(x) vapply(strsplit(x, "/"),
                             function(p) as.numeric(p[1]) / as.numeric(p[2]), numeric(1))

d <- bind_rows(
  vroom(file.path(ora_dir, "5_2_3_Hallmark7_louv_max.tsv"), show_col_types = FALSE) |>
    mutate(Network = "HPV-A7"),
  vroom(file.path(ora_dir, "5_2_3_Hallmark9_louv_max.tsv"), show_col_types = FALSE) |>
    mutate(Network = "HPV-A9")) |>
  mutate(HM        = sub("^HALLMARK_", "", ID),
         GeneRatio = gr_num(GeneRatio),
         neglogp   = -log10(p.adjust),
         Network   = factor(Network, c("HPV-A7", "HPV-A9")),
         col       = paste(Network, Community_Name))
## columns run M1, M2, ... inside each network, as in panel a
d <- d |> mutate(col = factor(col, d |> distinct(Network, Community_Name, col) |>
                                     arrange(Network, num(Community_Name)) |> pull(col)))

## =============================== panel a =====================================
g7 <- readRDS(file.path(mod_dir, "5_1_1_Louvain_A7_graph_max.rds"))
g9 <- readRDS(file.path(mod_dir, "5_1_1_Louvain_A9_graph_max.rds"))
s7 <- split(V(g7)$name, V(g7)$Community_name)
s9 <- split(V(g9)$name, V(g9)$Community_name)
## modules enter largest first, as in 5_1_3B, so the clustering ties break the same way
sz <- function(sp) names(sp)[order(lengths(sp), decreasing = TRUE)]
s7 <- s7[sz(s7)]; s9 <- s9[sz(s9)]

J <- outer(seq_along(s7), seq_along(s9), Vectorize(function(i, j)
  length(intersect(s7[[i]], s9[[j]])) / length(union(s7[[i]], s9[[j]]))))
dimnames(J) <- list(names(s7), names(s9))

## same axis order as 5_1_3B: pheatmap's defaults, euclidean + complete linkage.
## The tree is built on the complete matrix, so the order is the one of that
## figure; only afterwards are the modules that panel b never shows dropped.
ord_r <- rownames(J)[hclust(dist(J),    method = "complete")$order]
ord_c <- colnames(J)[hclust(dist(t(J)), method = "complete")$order]

## a module with no enrichment underneath adds nothing here, so it is removed;
## the members of a marked pair are kept even so, because the text discusses them
keep7 <- union(d$Community_Name[d$Network == "HPV-A7"], focus_pair$A7)
keep9 <- union(d$Community_Name[d$Network == "HPV-A9"], focus_pair$A9)
message("[check] panel a trimmed to the modules panel b shows: ",
        length(ord_r), " -> ", sum(ord_r %in% keep7), " HPV-A7 rows (dropped ",
        paste(setdiff(ord_r, keep7), collapse = ", "), "), ",
        length(ord_c), " -> ", sum(ord_c %in% keep9), " HPV-A9 columns (dropped ",
        paste(setdiff(ord_c, keep9), collapse = ", "), ")")
ord_r <- ord_r[ord_r %in% keep7]
ord_c <- ord_c[ord_c %in% keep9]
J     <- J[ord_r, ord_c]

## numeric positions, so the module-size bars of 5_1_3B can be drawn next to the
## matrix; the cells themselves keep a single uniform size
pos7 <- setNames(seq(length(ord_r), 1), ord_r)   # first row of the tree on top
pos9 <- setNames(seq_along(ord_c), ord_c)
XL   <- length(ord_c); YL <- length(ord_r)

## The size bars are drawn in data units, and a cell is wider than it is tall,
## so the same number in x and in y would give two bars of different thickness.
## CELL_AR is the cell aspect measured on the rendered panel (129 / 220 px at
## the W x H set at the end of this script); it makes both bars come out equal.
## Re-measure it if the canvas or the number of modules changes a lot.
BAR     <- 0.62                      # bar thickness, in cell heights
GAP     <- 0.19                      # matrix-to-bar clearance, in cell heights
CELL_AR <- 129 / 220
## the module names sit at the bottom and on the right, as in 5_1_3B, so the
## size bars take the two free sides
bar_y   <- YL + 0.5 + GAP + BAR / 2               # A9 bar, above the matrix
bar_x   <- 0.5 - (GAP + BAR / 2) * CELL_AR        # A7 bar, left of the matrix

jac <- as.data.frame(as.table(J)) |> as_tibble() |>
  rlang::set_names(c("A7", "A9", "J")) |>
  mutate(A7 = as.character(A7), A9 = as.character(A9),
         x  = pos9[A9], y = pos7[A7])

## module size = number of genes in the community, taken from the graphs
sz7 <- tibble(A7 = ord_r, y = pos7[ord_r], n = as.integer(lengths(s7)[ord_r]))
sz9 <- tibble(A9 = ord_c, x = pos9[ord_c], n = as.integer(lengths(s9)[ord_c]))

## every spotlighted pair must exist, and its value comes from the matrix
stopifnot(all(focus_pair$A7 %in% ord_r), all(focus_pair$A9 %in% ord_c))
fp <- focus_pair |> left_join(jac, by = c("A7", "A9"))
message("[check] pairs the text discusses (values read from the matrix):\n",
        paste(sprintf("        %s-%s  J = %.2f", fp$A7, fp$A9, fp$J), collapse = "\n"))

## every overlap worth reading is printed, not only the ones the text names
num <- jac |> filter(J > J_SHOW) |> arrange(desc(J))
stopifnot(all(fp$J > J_SHOW))        # nothing the text names may go unlabelled
message("[check] ", nrow(num), " cells with J > ", J_SHOW, " carry a printed value; ",
        nrow(num) - nrow(fp), " of them are not named in the text:\n",
        paste(sprintf("        %s-%s  J = %.2f", num$A7, num$A9, num$J), collapse = "\n"))

pa <- ggplot() +
  geom_tile(data = jac, aes(x, y, fill = J), colour = "grey92", linewidth = 0.15) +
  geom_text(data = num, aes(x, y, label = sprintf("%.2f", J)),
            family = FAM, fontface = "bold", size = 5.0,
            colour = ifelse(num$J > 0.45, "white", "grey10")) +
  scale_fill_gradientn(colours = c("white", "#FFF3CD", "#E8602D", "#8B0000"),
                       limits = c(0, max(jac$J)), name = "Jaccard\nindex",
                       guide = guide_colourbar(order = 1)) +
  annotate("rect", xmin = 0.5, xmax = XL + 0.5, ymin = 0.5, ymax = YL + 0.5,
           fill = NA, colour = "grey30", linewidth = 0.7) +
  ## the two module-size bars of 5_1_3B, each on the side its axis is not using
  new_scale_fill() +
  geom_tile(data = sz9, aes(x, bar_y, fill = n), height = BAR, width = 1,
            colour = "grey55", linewidth = 0.2) +
  scale_fill_distiller(palette = "BuPu", direction = 1,
                       name = "Module size\n(HPV-A9)", guide = guide_colourbar(order = 2)) +
  new_scale_fill() +
  geom_tile(data = sz7, aes(bar_x, y, fill = n), width = BAR * CELL_AR, height = 1,
            colour = "grey55", linewidth = 0.2) +
  scale_fill_distiller(palette = "Greens", direction = 1,
                       name = "Module size\n(HPV-A7)", guide = guide_colourbar(order = 3)) +
  scale_x_continuous(breaks = pos9[ord_c], labels = ord_c,
                     limits = c(bar_x - BAR * CELL_AR / 2 - 0.06, XL + 0.5), expand = c(0, 0)) +
  scale_y_continuous(breaks = pos7[ord_r], labels = ord_r, position = "right",
                     limits = c(0.5, bar_y + BAR / 2 + 0.06), expand = c(0, 0)) +
  labs(x = "HPV-A9 modules", y = "HPV-A7 modules") +
  coord_cartesian(clip = "off") +
  theme_minimal(base_family = FAM) +
  theme(
    axis.text        = element_text(size = 16, face = "bold", colour = "grey15"),
    axis.ticks       = element_blank(),
    axis.title       = element_text(size = 24, face = "bold"),
    axis.title.x     = element_text(margin = margin(t = 8)),
    axis.title.y.right = element_text(margin = margin(l = 10)),
    panel.grid       = element_blank(),
    legend.title     = element_text(size = 17, face = "bold"),
    legend.text      = element_text(size = 16),
    legend.key.height = unit(0.62, "cm"),
    plot.margin      = margin(30, 48, 14, 48))


## =============================== panel b =====================================
## ---- row blocks: the MSigDB process categories, read from the reference -----
cat_tbl <- vroom(file.path(ora_dir, "5_2_10_hallmark_process_categories.tsv"),
                 comment = "#", show_col_types = FALSE)
stopifnot(nrow(cat_tbl) == 50, !anyDuplicated(cat_tbl$HM), all(d$HM %in% cat_tbl$HM))
## blocks are ranked by the enrichment they carry, not alphabetically: strongest
## significance first, gene ratio breaking the ties. The same rule orders the
## gene sets inside each block, so the panel reads top-down by strength.
hm_stat <- d |> group_by(HM) |>
  summarise(mx = max(neglogp), gr = max(GeneRatio), .groups = "drop") |>
  left_join(cat_tbl, by = "HM")
cat_lv  <- hm_stat |> group_by(category) |>
  summarise(mx = max(mx), gr = max(gr), .groups = "drop") |>
  arrange(desc(mx), desc(gr)) |> pull(category)
## strip fills: one mechanical hue per category, darkened so white text reads
cat_col <- setNames(colorRampPalette(c("grey25", "white"))(2)[1], NULL)
cat_col <- setNames(mapply(function(h) {
  m <- col2rgb(h) / 255
  rgb(0.55 * m[1] + 0.45 * 0.18, 0.55 * m[2] + 0.45 * 0.18, 0.55 * m[3] + 0.45 * 0.18)
}, hcl.colors(length(cat_lv), palette = "Dark 3")), cat_lv)

message("[check] ", n_distinct(d$HM), " enriched gene sets, blocked into ", length(cat_lv),
        " MSigDB process categories (Liberzon et al. 2015, Table 1)\n",
        "[check] block order, strongest first: ", paste(cat_lv, collapse = " > "))

key <- hm_stat |> mutate(category = factor(category, cat_lv)) |>
  arrange(category, mx, gr)         # ascending: y levels run bottom-to-top

d <- d |> left_join(key |> select(HM, category), by = "HM") |>
  mutate(HM = factor(HM, key$HM))

## panel b carries the same modules as panel a, in the same dendrogram order:
## a module with no enrichment keeps its column, empty, so the two panels can be
## read column by column instead of only sharing a sequence
lv    <- c(paste("HPV-A7", ord_r), paste("HPV-A9", ord_c))
d     <- d |> mutate(col = factor(as.character(col), lv))
## every column of every panel, so the axis is trained from the level order and
## not from the order in which the layers happen to be drawn
pad   <- tibble(col = factor(lv, lv)) |>
  mutate(Network = factor(sub(" .*$", "", col), c("HPV-A7", "HPV-A9")))
empty <- setdiff(lv, as.character(d$col))
x_lab <- setNames(sub("^HPV-A[79] ", "", lv), lv)
message("[check] module order taken from the panel a clustering:\n",
        "        HPV-A7  ", paste(ord_r, collapse = " "), "\n",
        "        HPV-A9  ", paste(ord_c, collapse = " "), "\n",
        "[check] columns kept empty (no enrichment): ",
        paste(sub("^HPV-", "", empty), collapse = ", "))

pb <- ggplot(d, aes(col, HM)) +
  geom_blank(data = pad, aes(x = col), inherit.aes = FALSE) +
  geom_point(aes(size = GeneRatio, colour = neglogp)) +
  facet_grid2(category ~ Network, scales = "free", space = "free",
              strip = strip_themed(
                background_x = elem_list_rect(fill = "grey25"),
                background_y = elem_list_rect(fill = cat_col[cat_lv]),
                text_x = elem_list_text(colour = "white", size = 20, face = "bold"),
                text_y = elem_list_text(colour = "white", size = 18, face = "bold"))) +
  scale_colour_viridis(option = "magma", begin = 0.12, end = 0.95, direction = -1,
                       name = expression(-log[10]~italic(p)[adj])) +
  scale_size_continuous(range = c(2.5, 11), name = "Gene ratio") +
  scale_x_discrete(labels = function(z) x_lab[z]) +
  scale_y_discrete(labels = function(z) gsub("_", " ", z)) +
  labs(x = "Module", y = NULL) +
  theme_bw(base_family = FAM) +
  theme(
    axis.text.x      = element_text(size = 15, face = "bold", colour = "grey15",
                                    angle = 45, hjust = 1, vjust = 1),
    axis.text.y      = element_text(size = 15, colour = "grey15"),
    axis.title.x     = element_text(size = 19, face = "bold", margin = margin(t = 8)),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.spacing.x  = unit(0.8, "lines"),
    panel.spacing.y  = unit(0.28, "lines"),
    panel.border     = element_rect(colour = "grey30", fill = NA, linewidth = 0.7),
    ## horizontal category labels: a two-row block would clip a rotated one
    strip.text.y.right = element_text(angle = 0),
    legend.title     = element_text(size = 17, face = "bold"),
    legend.text      = element_text(size = 16),
    legend.key.height = unit(0.85, "cm"),
    plot.margin      = margin(14, 6, 4, 6))

## =============================== compose =====================================
## the heatmap is taken out of the alignment altogether so that it spans the
## whole page width instead of starting where the Hallmark names end
fig <- wrap_elements(full = pa) / pb +
  plot_layout(heights = c(7.4, 13.1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = FAM, size = 30, face = "bold"))

W <- 16; H <- 20.5
fig_dir <- "~/CESC_Network/figures"
out <- file.path(fig_dir, "Figure_3_modules_hallmark")
ggsave(paste0(out, ".png"), fig, width = W, height = H, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(paste0(out, ".pdf"), fig, width = W, height = H, device = cairo_pdf, bg = "white", limitsize = FALSE)
message("[OK] saved ", out, ".png / .pdf")
