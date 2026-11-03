### ============================================================
### figure_3_combined.R
###
### CESC + HTMCP-CC combined view of Figure 3. Does NOT touch the original
### CESC-only figure. Grouped by ANALYSIS TYPE: each row is one panel type
### (Hallmark dot-plot, module-overlap Jaccard heatmap), CESC and HTMCP-CC
### side by side under the same tag letter.
###
###   A) Hallmark over-representation per module   [CESC | HTMCP-CC]
###   B) Module-overlap Jaccard heatmap            [CESC | HTMCP-CC]
###
### CESC keeps its historical focus_pair editorial selection and purple/
### blue/red palette; HTMCP-CC has no manuscript text yet to cite specific
### module pairs, so that editorial block is skipped for it, and it uses
### one consistent orange/green pair with an olive "shared" ramp.
###
### -> figures/validation_uganda/Figure_3_modules_hallmark_uganda.png/.pdf
### ============================================================
suppressPackageStartupMessages({
  library(tidyverse); library(vroom); library(igraph)
  library(ggh4x);    library(patchwork); library(viridis); library(ggnewscale)
})

ref_dir <- "~/CESC_Network/06_modules/results"
fig_dir <- "~/CESC_Network/figures/validation_uganda"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

FAM     <- "Nimbus Roman"
J_SHOW  <- 0.1
num <- function(x) as.integer(sub("M", "", x))

build_cohort <- function(mod_dir, ora_dir, col_A7, col_A9, COHORT, cohort_label) {

focus_pair <- if (COHORT != "uganda") tribble(
  ~A7,   ~A9,
  "M6",  "M3", "M1",  "M1", "M2",  "M2", "M8",  "M13",
  "M10", "M19", "M12", "M17", "M15", "M23") else tibble(A7 = character(), A9 = character())

gr_num <- function(x) vapply(strsplit(x, "/"), function(p) as.numeric(p[1]) / as.numeric(p[2]), numeric(1))
d <- bind_rows(
  vroom(file.path(ora_dir, "5_2_3_Hallmark7_louv_max.tsv"), show_col_types = FALSE) |> mutate(Network = "HPV-A7"),
  vroom(file.path(ora_dir, "5_2_3_Hallmark9_louv_max.tsv"), show_col_types = FALSE) |> mutate(Network = "HPV-A9")) |>
  mutate(HM = sub("^HALLMARK_", "", ID), GeneRatio = gr_num(GeneRatio), neglogp = -log10(p.adjust),
         Network = factor(Network, c("HPV-A7", "HPV-A9")), col = paste(Network, Community_Name))
d <- d |> mutate(col = factor(col, d |> distinct(Network, Community_Name, col) |> arrange(Network, num(Community_Name)) |> pull(col)))

g7 <- readRDS(file.path(mod_dir, "5_1_1_Louvain_A7_graph_max.rds"))
g9 <- readRDS(file.path(mod_dir, "5_1_1_Louvain_A9_graph_max.rds"))
s7 <- split(V(g7)$name, V(g7)$Community_name); s9 <- split(V(g9)$name, V(g9)$Community_name)
sz <- function(sp) names(sp)[order(lengths(sp), decreasing = TRUE)]
s7 <- s7[sz(s7)]; s9 <- s9[sz(s9)]
J <- outer(seq_along(s7), seq_along(s9), Vectorize(function(i, j) length(intersect(s7[[i]], s9[[j]])) / length(union(s7[[i]], s9[[j]]))))
dimnames(J) <- list(names(s7), names(s9))

order_pairs <- function(Jg, Jt) {
  Tz <- ifelse(is.na(Jt), 0, Jt); S <- (Jg > J_SHOW) + (Tz > J_SHOW) + (Jg + Tz) / 4
  dimnames(S) <- dimnames(Jg); r <- rownames(S); c <- colnames(S); R <- character(); C <- character()
  while (length(r) && length(c)) {
    k <- which(S[r, c, drop = FALSE] == max(S[r, c, drop = FALSE]), arr.ind = TRUE)[1, ]
    R <- c(R, r[k[1]]); C <- c(C, c[k[2]]); r <- r[-k[1]]; c <- c[-k[2]]
  }
  list(r = c(R, if (length(r)) r[order(apply(S[r, , drop = FALSE], 1, max), decreasing = TRUE)]),
       c = c(C, if (length(c)) c[order(apply(S[, c, drop = FALSE], 2, max), decreasing = TRUE)]))
}
ord_r <- rownames(J); ord_c <- colnames(J)
keep7 <- d$Community_Name[d$Network == "HPV-A7"]; keep9 <- d$Community_Name[d$Network == "HPV-A9"]
ord_r <- ord_r[ord_r %in% keep7]; ord_c <- ord_c[ord_c %in% keep9]; J <- J[ord_r, ord_c]

tt7 <- with(filter(d, Network == "HPV-A7"), split(as.character(HM), Community_Name))
tt9 <- with(filter(d, Network == "HPV-A9"), split(as.character(HM), Community_Name))
gset <- function(L, k) if (is.null(L[[k]])) character(0) else L[[k]]
jj <- function(a, b) { u <- length(union(a, b)); if (!u) NA_real_ else length(intersect(a, b)) / u }
JT <- outer(ord_r, ord_c, Vectorize(function(a, b) jj(gset(tt7, a), gset(tt9, b)))); dimnames(JT) <- list(ord_r, ord_c)
o <- order_pairs(J, JT); ord_r <- o$r; ord_c <- o$c; J <- J[ord_r, ord_c]

pos7 <- setNames(seq(length(ord_r), 1), ord_r); pos9 <- setNames(seq_along(ord_c), ord_c)
XL <- length(ord_c); YL <- length(ord_r)
BAR <- 0.62; GAP <- 0.19; CELL_AR <- 129 / 220
bar_y <- YL + 0.5 + GAP + BAR / 2; bar_x <- 0.5 - (GAP + BAR / 2) * CELL_AR
jac <- as.data.frame(as.table(J)) |> as_tibble() |> rlang::set_names(c("A7", "A9", "J")) |>
  mutate(A7 = as.character(A7), A9 = as.character(A9), x = pos9[A9], y = pos7[A7])
sz7 <- tibble(A7 = ord_r, y = pos7[ord_r], n = as.integer(lengths(s7)[ord_r]))
sz9 <- tibble(A9 = ord_c, x = pos9[ord_c], n = as.integer(lengths(s9)[ord_c]))

if (COHORT == "uganda") {
  mid_rgb <- round((col2rgb(col_A7) + col2rgb(col_A9)) / 2)
  col_mid <- rgb(mid_rgb[1], mid_rgb[2], mid_rgb[3], maxColorValue = 255)
  JAC_RAMP <- colorRampPalette(c("white", "#EAEAD0", "#B8C08A", col_mid, "#3F4A1E"))(64)
  A7_RAMP <- colorRampPalette(c("white", col_A7))(64); A9_RAMP <- colorRampPalette(c("white", col_A9))(64)
} else {
  JAC_RAMP <- colorRampPalette(c("white", "#EDE8F4", "#C9BCDF", "#8C79B5", "#5E3C99", "#33215A"))(64)
  A7_RAMP <- colorRampPalette(c("white", "#4E79A7", "#003366"))(64)
  A9_RAMP <- colorRampPalette(c("white", "#CB6A5C", "#CC2200"))(64)
}
jac_lum2 <- function(v) { m <- grDevices::colorRamp(JAC_RAMP)(pmin(pmax(v, 0), 1)) / 255; 0.299*m[,1]+0.587*m[,2]+0.114*m[,3] }

t7 <- with(filter(d, Network == "HPV-A7"), split(as.character(HM), Community_Name))
t9 <- with(filter(d, Network == "HPV-A9"), split(as.character(HM), Community_Name))
term_of <- function(L, k) if (is.null(L[[k]])) character(0) else L[[k]]
jac_of <- function(a, b) { u <- length(union(a, b)); if (!u) NA_real_ else length(intersect(a, b)) / u }
both <- jac |> rowwise() |> mutate(J_term = { a <- term_of(t7, A7); b <- term_of(t9, A9); if (!length(a) || !length(b)) NA_real_ else jac_of(a, b) }) |>
  ungroup() |> rename(J_gene = J)
tri <- function(dat, value, half) {
  dx <- list(gene = c(-.5,.5,-.5), term = c(-.5,.5,.5))[[half]]; dy <- list(gene = c(-.5,.5,.5), term = c(-.5,-.5,.5))[[half]]
  dat |> mutate(id = paste(A7, A9, half)) |> slice(rep(seq_len(n()), each = 3)) |>
    mutate(vx = x + rep(dx, times = nrow(dat)), vy = y + rep(dy, times = nrow(dat)), V = .data[[value]])
}
poly <- bind_rows(tri(both, "J_gene", "gene"), tri(both, "J_term", "term"))
lab_g <- both |> filter(J_gene > J_SHOW); lab_t <- both |> filter(!is.na(J_term), J_term > J_SHOW)

pa <- ggplot() +
  geom_polygon(data = poly, aes(vx, vy, group = id, fill = V), colour = "grey92", linewidth = 0.15) +
  geom_text(data = lab_g, aes(x - 0.19, y + 0.20, label = sprintf("%.2f", J_gene)), family = FAM, fontface = "bold",
            size = 4.2, colour = ifelse(jac_lum2(lab_g$J_gene) < 0.55, "white", "grey10")) +
  geom_text(data = lab_t, aes(x + 0.19, y - 0.20, label = sprintf("%.2f", J_term)), family = FAM, fontface = "bold",
            size = 4.2, colour = ifelse(jac_lum2(lab_t$J_term) < 0.55, "white", "grey10")) +
  scale_fill_gradientn(colours = JAC_RAMP, limits = c(0, 1), na.value = "grey85", name = "Jaccard index", guide = guide_colourbar(order = 1)) +
  geom_point(data = tibble(x = rep(mean(range(jac$x)), 2), y = rep(1, 2),
                          half = factor(c("Gene membership", "Enriched gene sets"), c("Gene membership", "Enriched gene sets"))),
             aes(x, y, shape = half), alpha = 0, show.legend = TRUE, fill = JAC_RAMP[34], colour = "grey35", stroke = 0.5) +
  scale_shape_manual(values = c(24, 25), name = "Cell half", guide = guide_legend(order = 4, override.aes = list(alpha = 1, size = 6))) +
  annotate("rect", xmin = 0.5, xmax = XL + 0.5, ymin = 0.5, ymax = YL + 0.5, fill = NA, colour = "grey30", linewidth = 0.7) +
  new_scale_fill() +
  geom_tile(data = sz9, aes(x, bar_y, fill = n), height = BAR, width = 1, colour = "grey55", linewidth = 0.2) +
  scale_fill_gradientn(colours = A9_RAMP, name = "Module size\n(HPV-A9)", guide = guide_colourbar(order = 2)) +
  new_scale_fill() +
  geom_tile(data = sz7, aes(bar_x, y, fill = n), width = BAR * CELL_AR, height = 1, colour = "grey55", linewidth = 0.2) +
  scale_fill_gradientn(colours = A7_RAMP, name = "Module size\n(HPV-A7)", guide = guide_colourbar(order = 3)) +
  scale_x_continuous(breaks = pos9[ord_c], labels = ord_c, limits = c(bar_x - BAR * CELL_AR / 2 - 0.06, XL + 0.5), expand = c(0, 0)) +
  scale_y_continuous(breaks = pos7[ord_r], labels = ord_r, position = "right", limits = c(0.5, bar_y + BAR / 2 + 0.06), expand = c(0, 0)) +
  labs(x = "HPV-A9 modules", y = "HPV-A7 modules", title = cohort_label) +
  coord_cartesian(clip = "off") + theme_minimal(base_family = FAM) +
  theme(axis.text = element_text(size = 17, face = "bold", colour = "grey15"), axis.ticks = element_blank(),
        axis.title = element_text(size = 21, face = "bold"), axis.title.x = element_text(margin = margin(t = 8)),
        axis.title.y.right = element_text(margin = margin(l = 10)), panel.grid = element_blank(),
        legend.title = element_text(size = 17, face = "bold"), legend.text = element_text(size = 16),
        legend.key.height = unit(0.42, "cm"), legend.spacing.y = unit(0.16, "cm"), legend.margin = margin(2,2,2,2),
        legend.position = "left", plot.margin = margin(20, 20, 8, 130),
        plot.title = element_text(size = 23, face = "bold", hjust = 0.5))

d_raw <- d   # snapshot before panel-b's per-cohort category/HM re-leveling, for the merged dot-plot

cat_tbl <- vroom(file.path(ref_dir, "5_2_10_hallmark_process_categories.tsv"), comment = "#", show_col_types = FALSE)
hm_stat <- d |> group_by(HM) |> summarise(mx = max(neglogp), gr = max(GeneRatio), .groups = "drop") |> left_join(cat_tbl, by = "HM")
cat_lv <- hm_stat |> group_by(category) |> summarise(mx = max(mx), gr = max(gr), .groups = "drop") |> arrange(desc(mx), desc(gr)) |> pull(category)
cat_col <- setNames(mapply(function(h) { m <- col2rgb(h)/255; rgb(0.55*m[1]+0.45*0.18, 0.55*m[2]+0.45*0.18, 0.55*m[3]+0.45*0.18) },
                           hcl.colors(length(cat_lv), palette = "Dark 3")), cat_lv)
key <- hm_stat |> mutate(category = factor(category, cat_lv)) |> arrange(category, mx, gr)
d <- d |> left_join(key |> select(HM, category), by = "HM") |> mutate(HM = factor(HM, key$HM))
lv <- c(paste("HPV-A7", ord_r), paste("HPV-A9", ord_c)); d <- d |> mutate(col = factor(as.character(col), lv))
pad <- tibble(col = factor(lv, lv)) |> mutate(Network = factor(sub(" .*$", "", col), c("HPV-A7", "HPV-A9")))
x_lab <- setNames(sub("^HPV-A[79] ", "", lv), lv)
tint <- function(h, k = 0.5) { m <- col2rgb(h)/255; rgb(m[1]+(1-m[1])*k, m[2]+(1-m[2])*k, m[3]+(1-m[3])*k) }
STRIP_X <- c(tint(col_A7), tint(col_A9))

pb <- ggplot(d, aes(col, HM)) +
  geom_blank(data = pad, aes(x = col), inherit.aes = FALSE) +
  geom_point(aes(size = GeneRatio, colour = neglogp)) +
  facet_grid2(category ~ Network, scales = "free", space = "free",
              strip = strip_themed(background_x = elem_list_rect(fill = STRIP_X), background_y = elem_list_rect(fill = cat_col[cat_lv]),
                                   text_x = elem_list_text(colour = "grey10", size = 15, face = "bold"),
                                   text_y = elem_list_text(colour = "white", size = 13, face = "bold"))) +
  scale_colour_viridis(option = "magma", begin = 0.12, end = 0.95, direction = -1, name = expression(-log[10]~italic(p)[adj])) +
  scale_size_continuous(range = c(1.8, 8), name = "Gene ratio") +
  scale_x_discrete(labels = function(z) x_lab[z]) + scale_y_discrete(labels = function(z) gsub("_", " ", z)) +
  labs(title = cohort_label, x = "Module", y = NULL) + theme_bw(base_family = FAM) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 5)),
        axis.text.x = element_text(size = 10, face = "bold", colour = "grey15", angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 10, colour = "grey15"), axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 6)),
        panel.grid.major = element_line(colour = "grey90", linewidth = 0.3), panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0.5, "lines"), panel.spacing.y = unit(0.18, "lines"),
        panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.7),
        strip.text.y.right = element_text(angle = 0), legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11), legend.key.height = unit(0.6, "cm"), plot.margin = margin(10, 4, 4, 4))

list(pa = pa, pb = pb, d = d_raw, ord_r = ord_r, ord_c = ord_c)
}

cesc   <- build_cohort("~/CESC_Network/06_modules/results", "~/CESC_Network/06_modules/results", "#4E79A7", "#CB6A5C", "cesc", "CESC")
uganda <- build_cohort("~/CESC_Network/06_modules/validation_uganda/results", "~/CESC_Network/06_modules/validation_uganda/results", "#D55E00", "#009E73", "uganda", "HTMCP-CC")

## ── merged Hallmark dot-plot: ONE shared category column (right), instead of
## one per cohort -- the 50 MSigDB Hallmark gene sets and their process
## category are the same reference table for both cohorts, so showing that
## column twice was redundant. The four module-group columns (CESC-A7,
## CESC-A9, Uganda-A7, Uganda-A9) keep each network's own colour and module
## order (from that cohort's own Jaccard-based clustering in its heatmap,
## panel B/rows below, which stay separate since the two cohorts' modules
## are not comparable one-to-one).
build_merged_dotplot <- function(cesc, uganda) {
  tag <- function(d, cohort) d |> mutate(Cohort = cohort,
    grp = factor(paste(cohort, Network), levels = paste(cohort, c("HPV-A7", "HPV-A9"))))
  d_all <- bind_rows(tag(cesc$d, "CESC"), tag(uganda$d, "HTMCP-CC"))

  ## column order: CESC-A7 (its own clustering order), CESC-A9, HTMCP-CC-A7, HTMCP-CC-A9
  mk_lv <- function(cohort, net, ord) paste(cohort, net, ord)
  lv <- c(mk_lv("CESC", "HPV-A7", cesc$ord_r),   mk_lv("CESC", "HPV-A9", cesc$ord_c),
         mk_lv("HTMCP-CC", "HPV-A7", uganda$ord_r), mk_lv("HTMCP-CC", "HPV-A9", uganda$ord_c))
  d_all <- d_all |> mutate(col = factor(paste(Cohort, Network, Community_Name), levels = lv))
  pad <- tibble(col_chr = lv) |>
    separate(col_chr, into = c("Cohort", "Network", "Mod"), sep = " ", remove = FALSE, extra = "merge") |>
    mutate(grp = factor(paste(Cohort, Network), levels = levels(d_all$grp)), col = factor(col_chr, lv))

  ## category blocks and Hallmark order: recomputed ONCE from the UNION of both
  ## cohorts' enriched gene sets, so the row blocking is identical for all four columns
  cat_tbl <- vroom(file.path(ref_dir, "5_2_10_hallmark_process_categories.tsv"), comment = "#", show_col_types = FALSE)
  hm_stat <- d_all |> group_by(HM) |> summarise(mx = max(neglogp), gr = max(GeneRatio), .groups = "drop") |> left_join(cat_tbl, by = "HM")
  cat_lv <- hm_stat |> group_by(category) |> summarise(mx = max(mx), gr = max(gr), .groups = "drop") |> arrange(desc(mx), desc(gr)) |> pull(category)
  cat_col <- setNames(mapply(function(h) { m <- col2rgb(h)/255; rgb(0.55*m[1]+0.45*0.18, 0.55*m[2]+0.45*0.18, 0.55*m[3]+0.45*0.18) },
                             hcl.colors(length(cat_lv), palette = "Dark 3")), cat_lv)
  key <- hm_stat |> mutate(category = factor(category, cat_lv)) |> arrange(category, mx, gr)
  d_all <- d_all |> left_join(key |> select(HM, category), by = "HM") |> mutate(HM = factor(HM, key$HM))

  x_lab <- setNames(sub("^\\S+ HPV-A[79] ", "", lv), lv)
  tint <- function(h, k = 0.5) { m <- col2rgb(h)/255; rgb(m[1]+(1-m[1])*k, m[2]+(1-m[2])*k, m[3]+(1-m[3])*k) }
  STRIP_X <- c(tint("#4E79A7"), tint("#CB6A5C"), tint("#D55E00"), tint("#009E73"))
  names(STRIP_X) <- levels(d_all$grp)

  ggplot(d_all, aes(col, HM)) +
    geom_blank(data = pad, aes(x = col), inherit.aes = FALSE) +
    geom_point(aes(size = GeneRatio, colour = neglogp)) +
    facet_grid2(category ~ grp, scales = "free", space = "free",
                strip = strip_themed(background_x = elem_list_rect(fill = STRIP_X), background_y = elem_list_rect(fill = cat_col[cat_lv]),
                                     text_x = elem_list_text(colour = "grey10", size = 18, face = "bold"),
                                     text_y = elem_list_text(colour = "white", size = 20, face = "bold"))) +
    scale_colour_viridis(option = "magma", begin = 0.12, end = 0.95, direction = -1, name = expression(-log[10]~italic(p)[adj])) +
    scale_size_continuous(range = c(1.8, 8.5), name = "Gene ratio") +
    scale_x_discrete(labels = function(z) x_lab[z]) + scale_y_discrete(labels = function(z) gsub("_", " ", z)) +
    labs(x = NULL, y = NULL) + theme_bw(base_family = FAM) +
    theme(axis.text.x = element_text(size = 15, face = "bold", colour = "grey15", angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(size = 19, face = "bold", colour = "grey15"),
          panel.grid.major = element_line(colour = "grey90", linewidth = 0.3), panel.grid.minor = element_blank(),
          panel.spacing.x = unit(0.4, "lines"), panel.spacing.y = unit(0.16, "lines"),
          panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.7),
          strip.text.y.right = element_text(angle = 0), legend.title = element_text(size = 19, face = "bold"),
          legend.text = element_text(size = 17), legend.key.height = unit(0.85, "cm"), plot.margin = margin(10, 4, 4, 4))
}

rowA <- wrap_elements(full = build_merged_dotplot(cesc, uganda))
rowB <- wrap_elements(full = (wrap_elements(full = cesc$pa) | wrap_elements(full = uganda$pa)))

fig <- (rowA / rowB) +
  plot_layout(heights = c(2.1, 0.85)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = FAM, size = 28, face = "bold"))

ggsave(file.path(fig_dir, "Figure_3_modules_hallmark_uganda.png"), fig, width = 24, height = 23, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(file.path(fig_dir, "Figure_3_modules_hallmark_uganda.pdf"), fig, width = 24, height = 23, device = cairo_pdf, bg = "white", limitsize = FALSE)
message("[OK] saved Figure_3_modules_hallmark_uganda.png / .pdf  (A=Hallmark dotplot, B=Jaccard heatmap; each CESC | HTMCP-CC)")
