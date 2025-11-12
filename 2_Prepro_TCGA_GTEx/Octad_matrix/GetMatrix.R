## =========================
## EXTRACT + ANNOTATE OCTAD (counts) + GC% + GENE LENGTH (robust)
## =========================

## --- Required libraries ---
suppressPackageStartupMessages({
  library(octad)
  library(octad.db)
  library(rhdf5)
  library(dplyr)
  library(tibble)
  library(biomaRt)
  ## [ADDED] for GC%/length from sequence
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(Biostrings)
})

options(stringsAsFactors = FALSE)

## --- Paths and output folder ---
setwd("~/CESC_Network/2_Prepro_TCGA_GTEx/Octad_matrix/")
out_dir <- "Output_rds"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## --- Load sample lists (A7/A9) ---
muestras_clado_A7 <- readRDS("~/CESC_Network/6_OCTAD/6_0_A7_Samples.rds")
muestras_clado_A9 <- readRDS("~/CESC_Network/6_OCTAD/6_0_A9_Samples.rds")

## --- HDF5 (already present in your path) ---
h5_path <- "/STORAGE/csbig/jruiz/Octad/octad.counts.and.tpm.h5"
stopifnot(file.exists(h5_path))

## --- Sample metadata (EH7274) ---
phenoDF <- octad.db::get_ExperimentHub_data("EH7274")

## --- Safe intersections against EH7274 ---
case_A7 <- intersect(muestras_clado_A7, phenoDF$sample.id)
case_A9 <- intersect(muestras_clado_A9, phenoDF$sample.id)
if (length(case_A7) < length(muestras_clado_A7)) warning("Some A7 samples are not present in OCTAD/EH7274.")
if (length(case_A9) < length(muestras_clado_A9)) warning("Some A9 samples are not present in OCTAD/EH7274.")

## --- Cases (HPV A7 + A9) ---
cases_union <- unique(c(case_A7, case_A9))

## --- Controls: GTEx normal + TCGA adjacent (cervix) ---
controls_all <- phenoDF %>%
  filter(grepl("CERVIX", biopsy.site, ignore.case = TRUE),
         (data.source == "GTEX" & sample.type == "normal") |
           (data.source == "TCGA" & sample.type == "adjacent")) %>%
  pull(sample.id)
control_union <- unique(controls_all)

## --- Save a controls table (optional) ---
controls_tbl <- phenoDF %>% filter(sample.id %in% control_union)
saveRDS(controls_tbl, file.path(out_dir, "Control_GTEx_TCGAadj.rds"))

## =========================
## Unified COUNTS matrix (cases + controls)
## =========================
sample_vector_all <- unique(c(cases_union, control_union))

counts_all <- octad::loadOctadCounts(
  sample_vector = sample_vector_all,
  type = "counts",
  file = h5_path
)

## Optional subsets
counts_cases    <- counts_all[, colnames(counts_all) %in% cases_union,    drop = FALSE]
counts_controls <- counts_all[, colnames(counts_all) %in% control_union,  drop = FALSE]

## Save matrices
saveRDS(counts_all,      file.path(out_dir, "Counts_ALL_cases_controls.rds"))
saveRDS(counts_cases,    file.path(out_dir, "Counts_ONLY_cases.rds"))
saveRDS(counts_controls, file.path(out_dir, "Counts_ONLY_controls.rds"))

## =========================
## ROW ANNOTATION (EH7272) + [ADDED] GC% & gene length (from hg38 sequence)
## =========================

## 0) Row IDs (prefer matrix rownames in case HDF5 differs)
tr_ids_h5 <- as.character(rhdf5::h5read(h5_path, "meta/transcripts"))
row_ids   <- rownames(counts_all)
if (!identical(tr_ids_h5, row_ids)) {
  message("Note: HDF5 transcript IDs differ from matrix rownames; using matrix rownames.")
}
tr_ids <- row_ids

## 1) OCTAD official mapping (EH7272): Ensembl -> Symbol/biotype/...
merged_gene_info <- octad.db::get_ExperimentHub_data("EH7272")
merged_gene_info$ensembl <- as.character(merged_gene_info$ensembl)
merged_gene_info$Symbol  <- merged_gene_info$gene
merged_gene_info$gene    <- NULL

## Join by Ensembl (exact or stripped)
row_annot <- tibble(identifier = tr_ids) |>
  left_join(merged_gene_info, by = dplyr::join_by(identifier == ensembl))

if (sum(is.na(row_annot$Symbol)) > 0) {
  map_stripped <- merged_gene_info |>
    mutate(ensembl_stripped = sub("\\.\\d+$", "", ensembl))
  row_annot <- tibble(identifier = tr_ids,
                      id_stripped = sub("\\.\\d+$", "", tr_ids)) |>
    left_join(map_stripped, by = dplyr::join_by(id_stripped == ensembl_stripped)) |>
    select(identifier, Symbol, everything(), -id_stripped, -ensembl_stripped)
}

stopifnot(identical(rownames(counts_all), tr_ids))
save.image("/STORAGE/csbig/jruiz/MatrixOctad.RData")

## -------------------------
## [ADDED / FIXED] GC% and gene length without relying on BioMart GC%
## -------------------------
## Strategy:
## - Query Ensembl (current) just for coordinates (no GC% attribute).
## - Compute gc_percent and gene_length_bp from hg38 sequence.

## Prepare Ensembl IDs without version
ensembl_gene_novers <- sub("\\.\\d+$", "", tr_ids)

## Connect to Ensembl BioMart (current; no archive host)
mart <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

## Retrieve coordinates and minimal gene info
bm_attrs <- c("ensembl_gene_id",
              "external_gene_name",
              "gene_biotype",
              "chromosome_name",
              "start_position",
              "end_position",
              "strand")

bm_tbl <- biomaRt::getBM(
  attributes = bm_attrs,
  filters    = "ensembl_gene_id",
  values     = unique(ensembl_gene_novers),
  mart       = mart
)

## Keep standard chromosomes; coerce types
std_chrs <- c(as.character(1:22), "X", "Y", "MT")
bm_tbl <- bm_tbl %>%
  filter(chromosome_name %in% std_chrs) %>%
  mutate(
    start_position = as.integer(start_position),
    end_position   = as.integer(end_position),
    strand         = as.integer(strand)
  )

## Compute gene_length_bp from coordinates
bm_tbl <- bm_tbl %>%
  mutate(gene_length_bp = pmax(1L, end_position - start_position + 1L))

## Build GRanges (convert to UCSC style: "1"->"chr1", "MT"->"chrM")
chr_ucsc <- ifelse(bm_tbl$chromosome_name == "MT", "chrM", paste0("chr", bm_tbl$chromosome_name))
gr_genes <- GenomicRanges::GRanges(
  seqnames = chr_ucsc,
  ranges   = IRanges::IRanges(start = bm_tbl$start_position, end = bm_tbl$end_position),
  strand   = ifelse(bm_tbl$strand == 1, "+", "-")
)

## Extract sequences and compute GC%
seqs <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, gr_genes)
gc_counts  <- rowSums(Biostrings::letterFrequency(seqs, c("G","C")))
gc_percent <- 100 * (gc_counts / pmax(1, width(seqs)))

## Assemble tidy table keyed by Ensembl (no version)
bm_tbl2 <- bm_tbl %>%
  mutate(
    ensembl_gene_id_nover = bm_tbl$ensembl_gene_id,
    gc_percent            = as.numeric(gc_percent)
  ) %>%
  select(ensembl_gene_id_nover,
         external_gene_name, gene_biotype,
         chromosome_name, start_position, end_position, strand,
         gene_length_bp, gc_percent)

## Merge with row_annot keeping matrix row order
row_annot <- row_annot %>%
  mutate(ensembl_gene_id_nover = sub("\\.\\d+$", "", identifier)) %>%
  left_join(bm_tbl2, by = "ensembl_gene_id_nover") %>%
  ## Prefer OCTAD Symbol; keep BioMart as fallback
  mutate(Symbol = dplyr::coalesce(Symbol, external_gene_name)) %>%
  select(identifier, Symbol, everything(), -external_gene_name)

## QC messages
coverage_coords <- mean(!is.na(row_annot$gene_length_bp))
coverage_gc     <- mean(!is.na(row_annot$gc_percent))
message(sprintf("Coordinate coverage: %.1f%% | GC%% coverage: %.1f%%",
                100*coverage_coords, 100*coverage_gc))

## Save final row annotation
saveRDS(row_annot, file.path(out_dir, "RowAnnot_EH7272_PLUS_GC_Length_hg38_from_coords.rds"))

## =========================
## COLUMN ANNOTATION (sample metadata)
## =========================
h5_samples <- as.character(rhdf5::h5read(h5_path, "meta/samples"))
col_ids    <- intersect(colnames(counts_all), h5_samples)

col_annot <- phenoDF %>%
  filter(sample.id %in% col_ids) %>%
  select(sample.id, data.source, sample.type, biopsy.site, disease, tissue)

## Exact order as in the matrix
col_annot <- col_annot[match(colnames(counts_all), col_annot$sample.id), ]

## Final alignment check
stopifnot(identical(colnames(counts_all), col_annot$sample.id))

saveRDS(col_annot, file.path(out_dir, "ColAnnot_from_EH7274.rds"))

message("✓ Done: counts_all + row_annot (with GC% & length) + col_annot saved in: ", normalizePath(out_dir))
