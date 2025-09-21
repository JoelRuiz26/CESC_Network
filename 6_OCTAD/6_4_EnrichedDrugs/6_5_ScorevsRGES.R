# ============================================================
# Build a long table (category × drug) and plot sRGES per category
# - Uses the same filters as your bubble plot (padj, n_present, score_min)
# - Table columns: signature, target_type, target, score, pert_iname, sRGES
# - Plot: Y = category (with score), X = sRGES (one point per drug)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(viridis)
  library(octad.db)   # for get_ExperimentHub_data()
})

# -------------------- Paths & inputs --------------------
out_base    <- "~/CESC_Network/6_OCTAD/6_4_EnrichedDrugs"
enrich_all  <- readRDS(file.path(out_base, "6_4_1_enrichment_all_with_n_present.rds"))  # all + n_present
sRES_A7     <- readRDS("~/CESC_Network/6_OCTAD/6_3_RES/A7/RES_A7_common3_collapsed_FDA_Launched.rds")
sRES_A9     <- readRDS("~/CESC_Network/6_OCTAD/6_3_RES/A9/RES_A9_common3_collapsed_FDA_Launched.rds")

# -------------------- Filters (match your bubble plot) --------------------
alpha_cut   <- 0.01
min_present <- 3
score_min   <- 0.35      # set to 0 if you want *no* score cutoff
sRGES_thr   <- -0.20     # visual threshold line in the plot

# -------------------- Helper: load sets --------------------
load_cmpd_sets <- function(target_type = c("chembl_targets","mesh","ChemCluster")){
  target_type <- match.arg(target_type)
  eh_map <- c(ChemCluster="EH7266", chembl_targets="EH7267", mesh="EH7268")
  x <- octad.db::get_ExperimentHub_data(eh_map[[target_type]])
  sets <- x$cmpd.sets
  names(sets) <- x$cmpd.set.names
  sets
}

# -------------------- Normalize drug names --------------------
norm_names <- function(df) df %>% mutate(pert_iname = toupper(trimws(pert_iname)))
sA7 <- norm_names(sRES_A7)
sA9 <- norm_names(sRES_A9)

# -------------------- Expand membership: (signature, target_type, target) × pert_iname + sRGES --------------------
members_tbl_one <- function(signature = c("A7","A9"), target_type = c("mesh","chembl_targets")){
  signature   <- match.arg(signature)
  target_type <- match.arg(target_type)
  srg <- if (signature == "A7") sA7 else sA9
  sets <- load_cmpd_sets(target_type)
  sets <- lapply(sets, function(v) toupper(trimws(v)))
  purrr::imap_dfr(sets, ~ tibble(target = .y, pert_iname = .x)) %>%
    inner_join(srg, by = "pert_iname") %>%
    mutate(signature = signature, target_type = target_type)
}

members_all <- bind_rows(
  members_tbl_one("A7", "mesh"),
  members_tbl_one("A7", "chembl_targets"),
  members_tbl_one("A9", "mesh"),
  members_tbl_one("A9", "chembl_targets")
)

# -------------------- Keep exactly the categories used in your bubble plot --------------------
# (padj <= alpha_cut, n_present >= min_present, score >= score_min; MeSH & ChEMBL only)
bubble_sets <- enrich_all %>%
  filter(target_type %in% c("mesh","chembl_targets"),
         signature %in% c("A7","A9"),
         padj <= alpha_cut,
         n_present >= min_present,
         score >= score_min) %>%
  dplyr::select(signature, target_type, target, score) %>%
  distinct()

# -------------------- Final long table: add sRGES for each drug in each kept category --------------------
final_tbl <- members_all %>%
  inner_join(bubble_sets, by = c("signature","target_type","target")) %>%
  # category label for plotting (includes target_type, signature, and score)
  mutate(tt_short = dplyr::recode(target_type, mesh = "MeSH", chembl_targets = "ChEMBL"),
         cat_lab  = sprintf("[%s] %s | %s (score=%.3f)", tt_short, target, signature, score))

sRGES_thr <- -0.20
out_base  <- "~/CESC_Network/6_OCTAD/6_4_EnrichedDrugs"

# Split the long table by signature
dA7 <- final_tbl %>% filter(signature == "A7")
dA9 <- final_tbl %>% filter(signature == "A9")

# Order categories by score (descending) within each signature
lvl_A7 <- dA7 %>% distinct(cat_lab, score) %>% arrange(desc(score)) %>% pull(cat_lab)
lvl_A9 <- dA9 %>% distinct(cat_lab, score) %>% arrange(desc(score)) %>% pull(cat_lab)

# ---------- Plot for A7 ----------
p_A7 <- ggplot(dA7, aes(x = sRGES, y = factor(cat_lab, levels = lvl_A7))) +
  geom_vline(xintercept = sRGES_thr, linetype = "dashed", color = "grey40") +
  geom_point(aes(color = signature, shape = target_type),
             size = 1.9,
             alpha = ifelse(dA7$sRGES <= sRGES_thr, 1, 0.25),
             stroke = 0) +
  scale_color_viridis_d(option = "D", end = 0.85, name = "Signature") +
  scale_shape_discrete(name = "Target type") +
  labs(
    title = "Per-category drug sRGES — A7",
    x = "Drug sRGES", y = NULL,
    caption = sprintf("Dashed line: sRGES ≤ %.2f", sRGES_thr)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "right",
    plot.margin = margin(12, 20, 18, 20)
  )

# Dynamic height (about 0.28 in per category, min 6 in)
h_A7 <- max(6, length(lvl_A7) * 0.28)

#ggsave(file.path(out_base, "per_category_drug_sRGES_A7.pdf"),
#       p_A7, width = 12, height = h_A7, units = "in", dpi = 1200, useDingbats = FALSE)
#ggsave(file.path(out_base, "per_category_drug_sRGES_A7.png"),
#       p_A7, width = 12, height = h_A7, units = "in", dpi = 1200)

# ---------- Plot for A9 ----------
p_A9 <- ggplot(dA9, aes(x = sRGES, y = factor(cat_lab, levels = lvl_A9))) +
  geom_vline(xintercept = sRGES_thr, linetype = "dashed", color = "grey40") +
  geom_point(aes(color = signature, shape = target_type),
             size = 1.9,
             alpha = ifelse(dA9$sRGES <= sRGES_thr, 1, 0.25),
             stroke = 0) +
  scale_color_viridis_d(option = "D", end = 0.85, name = "Signature") +
  scale_shape_discrete(name = "Target type") +
  labs(
    title = "Per-category drug sRGES — A9",
    x = "Drug sRGES", y = NULL,
    caption = sprintf("Dashed line: sRGES ≤ %.2f", sRGES_thr)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "right",
    plot.margin = margin(12, 20, 18, 20)
  )

h_A9 <- max(6, length(lvl_A9) * 0.28)

#ggsave(file.path(out_base, "per_category_drug_sRGES_A9.pdf"),
#       p_A9, width = 12, height = h_A9, units = "in", dpi = 1200, useDingbats = FALSE)
#ggsave(file.path(out_base, "per_category_drug_sRGES_A9.png"),
#       p_A9, width = 12, height = h_A9, units = "in", dpi = 1200)

# (Optional) Print to screen
print(p_A7); print(p_A9)
