### ============================================================
### 5_2_3_2_Hallmark_chord_A7_A9.R
###
### Same information as 5_2_3_1_Hallmark_dotplot_A7_A9_topo.png, drawn as an
### attractive tripartite chord ("circos-like") diagram, DECLUTTERED so the
### signal reads clearly:
###   * LEFT column  = HPV-A7 network modules
###   * CENTRE       = Hallmark processes (shared axis)
###   * RIGHT column = HPV-A9 network modules
### Each curved ribbon = one enriched (module x Hallmark) pair.
###   ribbon COLOUR = -log10(adjusted p)   (magma, same scale as the dot-plot)
###   ribbon WIDTH  = GeneRatio            (same encoding as dot size)
###
### Only the strongest associations are drawn (p.adjust < P_CUT). Every hallmark
### in the ORA tables is already significant at 0.05, so a stricter cut is what
### actually removes the visual noise; hallmarks / modules left without any
### surviving ribbon are dropped from the layout.
###
###   5_2_3_2_Hallmark_chord_A7_A9.png / .pdf
### ============================================================
suppressPackageStartupMessages({
  library(tidyverse); library(vroom); library(ggforce); library(viridis)
})

ora_dir <- "~/CESC_Network/7_Modularity/5_2_ORA"
col_A7 <- "#4E79A7"; col_A9 <- "#CB6A5C"
P_CUT  <- 1e-3                         # significance threshold for drawn ribbons

gr_num <- function(x) vapply(strsplit(x, "/"),
                             function(p) as.numeric(p[1]) / as.numeric(p[2]), numeric(1))

h7 <- vroom(file.path(ora_dir, "5_2_3_Hallmark7_louv_max.tsv"), show_col_types = FALSE) |>
  mutate(Network = "HPV-A7")
h9 <- vroom(file.path(ora_dir, "5_2_3_Hallmark9_louv_max.tsv"), show_col_types = FALSE) |>
  mutate(Network = "HPV-A9")

d <- bind_rows(h7, h9) |>
  mutate(Hallmark  = gsub("_", " ", sub("^HALLMARK_", "", Description)),
         GeneRatio = gr_num(GeneRatio),
         neglogp   = -log10(p.adjust),
         mod_num   = as.integer(sub("M", "", Community_Name))) |>
  filter(p.adjust < P_CUT)             # <-- declutter: keep only the strongest pairs

## ---- node ordering & vertical positions -----------------------------------
## Hallmarks (centre): most broadly enriched at the TOP, matching the dot-plot.
hall_order <- d |> count(Hallmark) |> arrange(desc(n), Hallmark) |> pull(Hallmark)
nH <- length(hall_order)
y_hall <- setNames(seq(1, 0, length.out = nH), hall_order)

## Modules on each side: keep only those that still carry a ribbon, ordered by number.
mods7 <- d |> filter(Network == "HPV-A7") |> distinct(Community_Name, mod_num) |> arrange(mod_num)
mods9 <- d |> filter(Network == "HPV-A9") |> distinct(Community_Name, mod_num) |> arrange(mod_num)
y_m7 <- setNames(seq(0.97, 0.03, length.out = nrow(mods7)), mods7$Community_Name)
y_m9 <- setNames(seq(0.97, 0.03, length.out = nrow(mods9)), mods9$Community_Name)

## x geometry: modules at +/-1, ribbons stop short of the centre label band.
XM <- 1.0; XI <- 0.17
edges <- d |>
  mutate(
    x    = if_else(Network == "HPV-A7", -XM,  XI),
    y    = if_else(Network == "HPV-A7", y_m7[Community_Name], y_hall[Hallmark]),
    xend = if_else(Network == "HPV-A7", -XI,  XM),
    yend = if_else(Network == "HPV-A7", y_hall[Hallmark], y_m9[Community_Name]),
    id   = row_number()) |>
  arrange(neglogp)                                   # strongest ribbons drawn last (on top)

## ---- node data frames -----------------------------------------------------
nodes_mod <- bind_rows(
  tibble(mod = mods7$Community_Name, x = -XM, y = y_m7[mods7$Community_Name], net = "HPV-A7"),
  tibble(mod = mods9$Community_Name, x =  XM, y = y_m9[mods9$Community_Name], net = "HPV-A9"))
nodes_hall <- tibble(Hallmark = hall_order, x = 0, y = y_hall[hall_order])

## ---- plot -----------------------------------------------------------------
p <- ggplot() +
  ggforce::geom_diagonal(
    data = edges,
    aes(x = x, y = y, xend = xend, yend = yend, group = id,
        colour = neglogp, linewidth = GeneRatio, alpha = neglogp),
    strength = 0.55, lineend = "round", n = 120) +
  scale_colour_viridis(option = "magma", begin = 0.12, end = 0.95, direction = -1,
                       name = expression(-log[10]~italic(p)[adj])) +
  scale_linewidth_continuous(range = c(0.35, 3.0), name = "Gene ratio") +
  scale_alpha_continuous(range = c(0.55, 0.98), guide = "none") +
  ## module nodes + outer labels
  geom_point(data = nodes_mod, aes(x, y, fill = net), shape = 21, size = 5.5,
             colour = "white", stroke = 0.8, show.legend = FALSE) +
  geom_text(data = filter(nodes_mod, net == "HPV-A7"), aes(x, y, label = mod),
            hjust = 1, nudge_x = -0.035, size = 4.6, fontface = "bold", colour = col_A7) +
  geom_text(data = filter(nodes_mod, net == "HPV-A9"), aes(x, y, label = mod),
            hjust = 0, nudge_x = 0.035, size = 4.6, fontface = "bold", colour = col_A9) +
  scale_fill_manual(values = c("HPV-A7" = col_A7, "HPV-A9" = col_A9)) +
  ## centre hallmark labels
  geom_label(data = nodes_hall, aes(x, y, label = Hallmark),
             size = 3.3, fontface = "bold", label.size = 0, label.padding = unit(0.07, "lines"),
             fill = alpha("white", 0.88), colour = "grey15") +
  ## side headers
  annotate("text", x = -XM, y = 1.05, label = "HPV-A7 modules", fontface = "bold",
           size = 6, colour = col_A7, hjust = 0.5) +
  annotate("text", x =  XM, y = 1.05, label = "HPV-A9 modules", fontface = "bold",
           size = 6, colour = col_A9, hjust = 0.5) +
  annotate("text", x = 0, y = 1.05, label = "Cancer Hallmarks", fontface = "bold",
           size = 6, colour = "grey20", hjust = 0.5) +
  annotate("text", x = 0, y = -0.05,
           label = sprintf("Ribbons shown: enriched module x hallmark pairs with p[adj] < %s", format(P_CUT, scientific = FALSE)),
           size = 3.6, colour = "grey40", hjust = 0.5, fontface = "italic") +
  coord_cartesian(xlim = c(-1.5, 1.5), ylim = c(-0.06, 1.08), clip = "off") +
  guides(colour = guide_colourbar(order = 1, barheight = 9),
         linewidth = guide_legend(order = 2, override.aes = list(colour = "grey30"))) +
  theme_void(base_size = 14) +
  theme(legend.position = "right",
        legend.title = element_text(size = 13, face = "bold"),
        legend.text  = element_text(size = 11),
        plot.margin  = margin(10, 20, 10, 20))

out_png <- file.path(ora_dir, "5_2_3_2_Hallmark_chord_A7_A9.png")
out_pdf <- file.path(ora_dir, "5_2_3_2_Hallmark_chord_A7_A9.pdf")
ggsave(out_png, p, width = 13, height = 13, dpi = 350, bg = "white", limitsize = FALSE)
ggsave(out_pdf, p, width = 13, height = 13, device = cairo_pdf, bg = "white", limitsize = FALSE)
message("[OK] Saved chord diagram: ", out_png,
        sprintf("  (p.adj < %g -> %d hallmarks; A7 %d mods / %d edges; A9 %d mods / %d edges)",
                P_CUT, nH, nrow(mods7), sum(d$Network == "HPV-A7"),
                nrow(mods9), sum(d$Network == "HPV-A9")))
