# ========= Configuración de rutas =========
base_dir <- "/home/jjruiz/OCTAD_Cervical_Cancer/1_DGE_signature/1_0_Output_rds"

p_A7_deseq2 <- file.path(base_dir, "1_4_DE_A7_signific_DESeq2.rds")
p_A7_edger  <- file.path(base_dir, "1_4_DE_A7_signific_EdgeR.rds")
p_A7_limma  <- file.path(base_dir, "1_4_DE_A7_signific_limma.rds")

p_A9_deseq2 <- file.path(base_dir, "1_4_DE_A9_signific_DESeq2.rds")
p_A9_edger  <- file.path(base_dir, "1_4_DE_A9_signific_EdgeR.rds")
p_A9_limma  <- file.path(base_dir, "1_4_DE_A9_signific_limma.rds")

out_dir <- file.path(base_dir, "1_5_CoreGenes_min")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ========= Helper mínimo: leer símbolos =========
read_syms <- function(path) {
  x <- readRDS(path)
  # Asegurar carácter y únicos; quitar NA si hubiera
  unique(as.character(x$Symbol[!is.na(x$Symbol)]))
}

# ========= Intersecciones =========
core_A7 <- Reduce(intersect, list(
  read_syms(p_A7_deseq2),
  read_syms(p_A7_edger),
  read_syms(p_A7_limma)
))
core_A7 <- sort(core_A7)

core_A9 <- Reduce(intersect, list(
  read_syms(p_A9_deseq2),
  read_syms(p_A9_edger),
  read_syms(p_A9_limma)
))
core_A9 <- sort(core_A9)

# ========= Guardar "pull" (vectores) =========
writeLines(core_A7, file.path(out_dir, "core_genes_A7.txt"))
writeLines(core_A9, file.path(out_dir, "core_genes_A9.txt"))
saveRDS(core_A7, file.path(out_dir, "core_genes_A7.rds"))
saveRDS(core_A9, file.path(out_dir, "core_genes_A9.rds"))

# (Opcional) mensaje en consola
cat("# core A7:", length(core_A7), "\n")
cat("# core A9:", length(core_A9), "\n")
cat("Guardado en:\n", out_dir, "\n")
