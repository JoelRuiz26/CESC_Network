### Hallmark enrichment dot-plot (module x Hallmark) for HPV-A7 and HPV-A9
### Y = Hallmark processes (shared axis), X = modules, size = GeneRatio,
### colour = significance (-log10 adjusted p). Facets: A7 (left) | A9 (right).
suppressPackageStartupMessages({
  library(tidyverse)
  library(vroom)
  library(viridis)
})

ora_dir <- "~/CESC_Network/7_Modularity/5_2_ORA/"

gr_num <- function(x) vapply(strsplit(x, "/"),
                             function(p) as.numeric(p[1]) / as.numeric(p[2]), numeric(1))

h7 <- vroom(file.path(ora_dir, "5_2_3_Hallmark7_louv_max.tsv"), show_col_types = FALSE) |>
  mutate(Network = "HPV-A7")
h9 <- vroom(file.path(ora_dir, "5_2_3_Hallmark9_louv_max.tsv"), show_col_types = FALSE) |>
  mutate(Network = "HPV-A9")

d <- bind_rows(h7, h9) |>
  mutate(
    Hallmark  = gsub("_", " ", sub("^HALLMARK_", "", Description)),
    GeneRatio = gr_num(GeneRatio),
    neglogp   = -log10(p.adjust),
    mod_num   = as.integer(sub("M", "", Community_Name))
  )

## Axis ordering ---------------------------------------------------------------
# X: modules by number (M1, M2, ... within each facet)
mod_levels <- paste0("M", sort(unique(d$mod_num)))
d$Community_Name <- factor(d$Community_Name, levels = mod_levels)
# Y: Hallmarks ordered by how many (module x network) pairs they appear in,
#    so the most broadly enriched processes sit at the top
hall_order <- d |> count(Hallmark) |> arrange(n) |> pull(Hallmark)
d$Hallmark <- factor(d$Hallmark, levels = hall_order)
d$Network  <- factor(d$Network, levels = c("HPV-A7", "HPV-A9"))

## Plot ------------------------------------------------------------------------
n_hall <- nlevels(d$Hallmark)
n_cols <- nrow(dplyr::distinct(d, Network, Community_Name))  # columns actually drawn (A7 + A9)

p <- ggplot(d, aes(x = Community_Name, y = Hallmark)) +
  geom_point(aes(size = GeneRatio, colour = neglogp)) +
  facet_grid(. ~ Network, scales = "free_x", space = "free_x") +
  scale_colour_gradientn(
    colours = c("#F0C9A8", "#E98D6B", "#D14A61", "#8F3371", "#4A1E66"),
    name = expression(-log[10]~italic(p)[adj])) +
  scale_size_continuous(range = c(3, 12), name = "Gene ratio") +
  labs(x = "Module", y = NULL,
       title = "Hallmark enrichment across co-expression modules") +
  theme_bw(base_family = "Times New Roman") +
  theme(
    plot.title      = element_text(size = 30, face = "bold", hjust = 0.5, margin = margin(b = 12)),
    axis.title.x    = element_text(size = 26, face = "bold", margin = margin(t = 10)),
    axis.text.x     = element_text(size = 20, face = "bold", angle = 0, colour = "black"),
    axis.text.y     = element_text(size = 19, colour = "black"),
    strip.text      = element_text(size = 26, face = "bold", colour = "white"),
    strip.background = element_rect(fill = "grey25"),
    panel.grid.major = element_line(colour = "grey88", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.spacing   = unit(1.1, "lines"),
    panel.border    = element_rect(colour = "black", fill = NA, linewidth = 1),
    legend.title    = element_text(size = 22, face = "bold"),
    legend.text     = element_text(size = 19),
    legend.key.size = unit(1.1, "cm"),
    legend.position = "right"
  )

outfile <- file.path(ora_dir, "5_2_3_1_Hallmark_dotplot_A7_A9.png")
ggsave(outfile, p,
       width  = 8 + n_cols * 0.62,
       height = 3 + n_hall * 0.42,
       dpi = 300, bg = "white", limitsize = FALSE)
message("[OK] Saved: ", outfile, sprintf("  (%d hallmarks x %d columns)", n_hall, n_cols))
