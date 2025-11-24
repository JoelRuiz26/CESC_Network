## =========================
## CERVIX (TCGA + GTEx):
##   - Count matrix
##   - Sample metadata (cervix)
##   - Row annotation (GC%, gene length, Ensembl/Symbol/biotype)
## =========================

suppressPackageStartupMessages({
  library(octad)
  library(octad.db)
  library(dplyr)
  library(tibble)
  library(biomaRt)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(Biostrings)
})

options(stringsAsFactors = FALSE)

setwd("~/CESC_Network/1_Get_Data/Octad_matrix/")
out_dir <- "1_Output_rds"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

h5_path <- "/STORAGE/csbig/jruiz/Octad/octad.counts.and.tpm.h5"
stopifnot(file.exists(h5_path))

## -------------------------
## 1) Metadatos y filtro CÉRVIX
## -------------------------
phenoDF <- octad.db::get_ExperimentHub_data("EH7274")

cervix <- phenoDF %>%
  filter(
    grepl("cervix", biopsy.site, ignore.case = TRUE),
    data.source %in% c("TCGA", "GTEX")
  )

sample_vector_cervix <- cervix$sample.id

## -------------------------
## 2) Matriz de conteos CÉRVIX
## -------------------------
counts_cervix <- octad::loadOctadCounts(
  sample_vector = sample_vector_cervix,
  type          = "counts",
  file          = h5_path
)

## Alinear pheno a columnas de la matriz
cervix <- cervix[match(colnames(counts_cervix), cervix$sample.id), , drop = FALSE]
stopifnot(identical(colnames(counts_cervix), cervix$sample.id))

## =========================
## DES–TRANSFORMAR: log2(x+1) → x
#    doi: 10.1038/s41525-020-00167-4. PMID: 33495453; PMCID: PMC7835362.
## =========================
# counts_cervix está en log2(x+1)
expr_raw <- 2^counts_cervix - 1
storage.mode(expr_raw) <- "double"

# Asegurar nombres correctos
expr_raw <- expr_raw[rownames(counts_cervix), colnames(counts_cervix)]

saveRDS(cervix,        file.path(out_dir, "1_0_Cervix_pheno_TCGA_GTEx.rds"))
saveRDS(counts_cervix, file.path(out_dir, "1_1_Counts_CERVIX_TCGA_GTEx.rds"))
saveRDS(expr_raw,      file.path(out_dir, "1_1_1_Counts_CERVIX_TCGA_GTEx_raw.rds"))

## -------------------------
## 3) ROW ANNOTATION: EH7272 + GC% + length
## -------------------------
row_ids <- rownames(counts_cervix)

## 3.1 Info OCTAD (EH7272): Ensembl -> Symbol / biotipo / etc.
merged_gene_info <- octad.db::get_ExperimentHub_data("EH7272") %>%
  mutate(
    ensembl = as.character(ensembl),
    Symbol  = gene
  ) %>%
  dplyr::select(-gene)

row_annot <- tibble(identifier = row_ids) %>%
  left_join(merged_gene_info, dplyr::join_by(identifier == ensembl))

## Intento 2: empatar sin versión si faltan símbolos
if (any(is.na(row_annot$Symbol))) {
  map_stripped <- merged_gene_info %>%
    mutate(ensembl_stripped = sub("\\.\\d+$", "", ensembl))
  
  row_annot <- tibble(
    identifier  = row_ids,
    id_stripped = sub("\\.\\d+$", "", row_ids)
  ) %>%
    left_join(map_stripped, dplyr::join_by(id_stripped == ensembl_stripped)) %>%
    dplyr::select(identifier, Symbol, everything(), -id_stripped, -ensembl_stripped)
}

## 3.2 Coordenadas + GC% + longitud desde hg38 (SOLO COMO COVARIABLES)
ensembl_novers <- sub("\\.\\d+$", "", row_ids)

mart <- biomaRt::useEnsembl(biomart = "genes",
                            dataset = "hsapiens_gene_ensembl")

bm_tbl <- biomaRt::getBM(
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "gene_biotype",
    "chromosome_name",
    "start_position",
    "end_position",
    "strand"
  ),
  filters    = "ensembl_gene_id",
  values     = unique(ensembl_novers),
  mart       = mart
)

bm_tbl <- bm_tbl %>%
  dplyr::filter(chromosome_name %in% c(as.character(1:22), "X", "Y", "MT")) %>%
  dplyr::mutate(
    start_position = as.integer(start_position),
    end_position   = as.integer(end_position),
    strand         = as.integer(strand),
    gene_length_bp = pmax(1L, end_position - start_position + 1L)
  )

chr_ucsc <- ifelse(bm_tbl$chromosome_name == "MT", "chrM",
                   paste0("chr", bm_tbl$chromosome_name))

gr_genes <- GenomicRanges::GRanges(
  seqnames = chr_ucsc,
  ranges   = IRanges::IRanges(start = bm_tbl$start_position,
                              end   = bm_tbl$end_position),
  strand   = ifelse(bm_tbl$strand == 1, "+", "-")
)

genome_hg38 <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
seqs        <- BSgenome::getSeq(genome_hg38, gr_genes)

gc_counts  <- rowSums(Biostrings::letterFrequency(seqs, c("G", "C")))
gc_percent <- 100 * (gc_counts / pmax(1, Biostrings::width(seqs)))

bm_tbl2 <- bm_tbl %>%
  dplyr::mutate(
    ensembl_gene_id_nover = ensembl_gene_id,
    gc_percent            = as.numeric(gc_percent)
  ) %>%
  dplyr::select(
    ensembl_gene_id_nover,
    external_gene_name,      # lo guardamos aparte, NO sustituye a Symbol
    gene_biotype,
    chromosome_name,
    start_position,
    end_position,
    strand,
    gene_length_bp,
    gc_percent
  )

## IMPORTANTE: aquí NO tocamos Symbol de OCTAD
row_annot <- row_annot %>%
  dplyr::mutate(
    ensembl_gene_id_nover = sub("\\.\\d+$", "", identifier)
  ) %>%
  dplyr::left_join(bm_tbl2, by = "ensembl_gene_id_nover") %>%
  # Symbol se queda tal cual viene de EH7272 / OCTAD
  # external_gene_name se queda como columna aparte, si la quieres
  dplyr::select(identifier, Symbol, everything(), -ensembl_gene_id_nover)

#Ginecological cancer, exclude
row_annot <- row_annot %>% filter(chromosome_name != "Y")

saveRDS(row_annot, file.path(out_dir, "1_2_RowAnnot_CERVIX_GC_Length_hg38.rds"))

message("Listo: counts_cervix + cervix (pheno) + row_annot guardados en ", normalizePath(out_dir))
