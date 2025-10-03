# =========================================
# A7 vs A9 (FLIPPED AXES) — cells ordered ascending by mean medcor
# =========================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(ggplot2)
  library(octad); library(octad.db);library(forcats)  # optional if you prefer fct_* helpers
})

setwd("~/CESC_Network/6_OCTAD/")

# ---- Load sample IDs ----
sample_clado_A7 <- readRDS("6_0_A7_Samples.rds")
sample_clado_A9 <- readRDS("6_0_A9_Samples.rds")

# ---- Metadata and case IDs ----
phenoDF <- octad.db::get_ExperimentHub_data("EH7274")
case_A7 <- intersect(sample_clado_A7, phenoDF$sample.id)
case_A9 <- intersect(sample_clado_A9, phenoDF$sample.id)

# ---- Cell-line similarity (median correlation) ----
src <- "octad.small"  # or "octad.whole" if available
sim_A7 <- computeCellLine(case_id = case_A7, source = src)
sim_A9 <- computeCellLine(case_id = case_A9, source = src)

sim_A7_df <- sim_A7 %>%
  tibble::rownames_to_column("cell_line") %>%
  transmute(cell_line, A7 = medcor)

sim_A9_df <- sim_A9 %>%
  tibble::rownames_to_column("cell_line") %>%
  transmute(cell_line, A9 = medcor)
# ---- Join & order ASCENDING by mean medcor ----
sim_both <- full_join(sim_A7_df, sim_A9_df, by = "cell_line") %>%
  rowwise() %>%
  mutate(mean_medcor = mean(c_across(c(A7, A9)), na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(mean_medcor) %>%   # ascending order
  mutate(
    # Keep ascending order directly (no rev)
    cell_line = factor(cell_line, levels = cell_line),
    rank = row_number()
  )

sim_long <- sim_both %>%
  pivot_longer(c(A7, A9), names_to = "clade", values_to = "medcor") %>%
  mutate(clade = factor(clade, levels = c("A7", "A9")))

# ---- POINTS (flipped axes) ----
p_pts_flip <- ggplot(sim_long, aes(x = cell_line, y = medcor, color = clade)) +
  geom_point(position = position_dodge(width = 0.6), size = 2.2, na.rm = TRUE) +
  coord_flip() +
  labs(
    title = "Cell-line similarity (medcor)",
    x = "Cell line",
    y = "Median correlation (medcor)",
    color = "Clade"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 10)
  )

print(p_pts_flip)



# Optional save
out_dir <- "6_OCTAD_QCplots"; dir.create(out_dir, showWarnings = FALSE)
ggsave(file.path(out_dir, "A7_A9_medcor_points_Xcell_Ymedcor.png"),
       plot = p_pts, width = 14, height = 6.5, dpi = 1200, bg = "white")
ggsave(file.path(out_dir, "A7_A9_medcor_points_Xcell_Ymedcor.pdf"),
       plot = p_pts, width = 14, height = 6.5, device = cairo_pdf, bg = "white")





