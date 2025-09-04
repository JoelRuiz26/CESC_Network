# ================================
# 4-set Venn (A7_up, A7_down, A9_up, A9_down)
# - Input: ora_tbl (RDS with clade, regulation, ontology, Description)
# - Unifies ALL ontologies; uses UNIQUE terms (Description)
# - Output: single high-res PNG + preview in RStudio
# ================================

suppressPackageStartupMessages({
  library(dplyr)
  library(VennDiagram)  # supports up to 5-set Venns
  library(grid)         # for grid.draw() to preview in RStudio
})

# ---- Load ORA table ----
ora_tbl <- readRDS("~/CESC_Network/6_OCTAD/6_2_ORA_signature/6_2_2_Ora_tbl.rds")

# ---- Output folder ----
out_dir <- "~/CESC_Network/6_OCTAD/6_2_ORA_signature/6_2_7_Output_plots/"

# ---- Helper: get UNIQUE terms per clade × regulation (all ontologies) ----
get_terms_all <- function(df, cl, dr) {
  x <- df %>%
    filter(clade == cl, regulation == dr) %>%
    pull(Description) %>%
    unique()
  x[!is.na(x) & nzchar(x)]
}

# ---- Build the 4 sets (unique terms) ----
set_A7_UP   <- get_terms_all(ora_tbl, "A7", "up")
set_A7_DOWN <- get_terms_all(ora_tbl, "A7", "down")
set_A9_UP   <- get_terms_all(ora_tbl, "A9", "up")
set_A9_DOWN <- get_terms_all(ora_tbl, "A9", "down")

sets_list <- list(
  "A7_UP"   = set_A7_UP,
  "A7_DOWN" = set_A7_DOWN,
  "A9_UP"   = set_A9_UP,
  "A9_DOWN" = set_A9_DOWN
)

# ---- Build grob in memory (preview) and save high-res PNG ----
plot_title <- "GO terms: HPV-A7 / HPV-A9"
outfile    <- file.path(out_dir, "6_2_6_Venn_4sets_A7A9_UPDOWN.png")

# Create grob (no file written yet)
g <- VennDiagram::venn.diagram(
  x = sets_list,
  filename    = NULL,           # return a grob so we can preview and save ourselves
  fill        = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"),
  alpha       = c(0.45, 0.45, 0.45, 0.45),
  lwd         = 2,
  col         = "grey25",
  cex         = 1.2,             
  fontface    = "bold",
  cat.cex     = 1.2,            
  cat.fontface= "bold",
  main        = plot_title,
  main.cex    = 1.3
)

# Preview in RStudio Plots
grid::grid.newpage()
grid::grid.draw(g)

# Save as very high-resolution PNG (Cairo anti-aliased)
png(filename = outfile, width = 3200, height = 3200, res = 1200, type = "cairo-png")
grid::grid.draw(g)
dev.off()

cat("Ready. PNG saved at:\n", outfile, "\n")
