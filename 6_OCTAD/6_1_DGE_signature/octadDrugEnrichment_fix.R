# Requiere: dplyr, readr, GSVA, limma, octad.db (cárgalos antes)
octadDrugEnrichment_fix <- function (sRGES = NULL, target_type = "chembl_targets",
                                     enrichFolder = "enrichFolder",
                                     outputFolder = NULL, outputRank = FALSE) {
  if (missing(sRGES)) stop("sRGES input not found")
  if (is.null(sRGES$sRGES) | is.null(sRGES$pert_iname))
    stop("Either sRGES or pert_iname column in Disease signature is missing")
  
  if (is.null(outputFolder)) {
    outputFolder <- tempdir()
    message("outputFolder is NULL, writing output to tempdir()")
  }
  if (!dir.exists(outputFolder)) {
    stop("Looks like output path ", outputFolder,
         " is either non-existent or obstructed. Check outputFolder option.")
  }
  if (!dir.exists(file.path(outputFolder, enrichFolder))) {
    dir.create(file.path(outputFolder, enrichFolder), recursive = TRUE)
  }
  
  # Mapa ExperimentHub como en OCTAD
  eh_dataframe <- data.frame(title = c("cmpd_sets_ChemCluster",
                                       "cmpd_sets_chembl_targets", "cmpd_sets_mesh"))
  row.names(eh_dataframe) <- c("EH7266", "EH7267", "EH7268")
  random_gsea_score <- octad.db::get_ExperimentHub_data("EH7275")
  
  for (target_type_selected in target_type) {
    message("Running enrichment for ", target_type_selected, "\n")
    enrichFolder.n <- file.path(outputFolder, enrichFolder, target_type_selected)
    if (!dir.exists(enrichFolder.n)) dir.create(enrichFolder.n, recursive = TRUE)
    
    eh_dataframe$object <- row.names(eh_dataframe)
    cmpd_sets <- octad.db::get_ExperimentHub_data(
      (eh_dataframe[eh_dataframe$title == paste0("cmpd_sets_", target_type_selected), "object"])
    )
    cmpdSets <- cmpd_sets$cmpd.sets
    names(cmpdSets) <- cmpd_sets$cmpd.set.names
    
    drug_pred <- sRGES
    rgess <- matrix(-1 * drug_pred$sRGES, ncol = 1)
    if (is.null(dim(rgess))) {
      warning("rgess have zero rows, recompute it and try again")
      next
    }
    
    rownames(rgess) <- drug_pred$pert_iname
    rgess <- cbind(rgess, rgess)  # igual que OCTAD
    gsva_params <- GSVA::ssgseaParam(rgess, cmpdSets, normalize = TRUE)
    gsea_results <- GSVA::gsva(gsva_params, verbose = FALSE)
    gsea_results <- gsea_results[-1, , drop = FALSE]  # descarta la 1a fila/col según OCTAD
    
    # Empareja con la distribución nula precomputada
    gsea_results <- merge(random_gsea_score[[target_type_selected]], gsea_results, by = "row.names")
    if (is.null(dim(gsea_results)) || nrow(gsea_results) == 0L) {
      warning("gsea_results have zero rows, recompute it and try again")
      next
    }
    
    rownames(gsea_results) <- gsea_results$Row.names
    gsea_results$Row.names <- NULL
    
    # ------- p empírico con corrección de Monte Carlo: (r+1)/(N+1) -------
    Nnull <- ncol(random_gsea_score[[target_type_selected]])
    obs_col <- Nnull + 1L
    
    # Convierte a matriz numérica solo las columnas relevantes (nulos + observado)
    mat_num <- as.matrix(gsea_results[, seq_len(obs_col), drop = FALSE])
    storage.mode(mat_num) <- "numeric"
    
    # r = # de nulos > observado, por fila
    r <- rowSums(mat_num[, seq_len(Nnull), drop = FALSE] > mat_num[, obs_col])
    p_emp <- (r + 1) / (Nnull + 1)  # nunca 0
    score_obs <- as.numeric(mat_num[, obs_col])
    
    gsea_p <- dplyr::tibble(
      target = rownames(gsea_results),
      score  = score_obs,
      p      = as.numeric(p_emp),
      padj   = p.adjust(p_emp, method = "fdr"),
      Nnull  = Nnull,
      p_min  = 1/(Nnull + 1)
    )
    
    # Ordena empujando NAs al final (robusto)
    gsea_p <- gsea_p[order(gsea_p$padj, na.last = TRUE), , drop = FALSE]
    
    # Escribe CSV (alta precisión, sin redondeos silenciosos)
    readr::write_csv(gsea_p, file.path(enrichFolder.n,
                                       paste0("enriched_", target_type_selected, ".csv")))
    
    # ----------------- Export de gráficas y ranks (robusto a NA) -----------------
    if (nrow(gsea_p) > 0) {
      # Máscara de significancia sin NAs
      sig_mask <- is.finite(gsea_p$padj) & (gsea_p$padj <= 0.05)
      
      # Índices a graficar: si no hay significativos, grafica el primero; tope 50
      idx <- which(sig_mask)
      if (length(idx) == 0L) idx <- 1L
      idx <- head(idx, 50L)
      
      for (i in idx) {
        top_target <- as.character(gsea_p$target[i])
        sRGES$rank <- rank(sRGES$sRGES)
        target_drugs_score <- sRGES$rank[sRGES$pert_iname %in% cmpdSets[[top_target]]]
        if (length(target_drugs_score) < 3) next
        
        pdf(file.path(enrichFolder.n, paste0("/top_enriched_", top_target, "_",
                                             target_type_selected, ".pdf")))
        limma::barcodeplot(sRGES$sRGES, target_drugs_score, main = top_target, xlab = "sRGES")
        dev.off()
        
        if (outputRank) {
          readr::write_csv(
            sRGES[sRGES$pert_iname %in% cmpdSets[[top_target]], c("pert_iname","rank")],
            file.path(enrichFolder.n, paste0("/top_enriched_", top_target, "_",
                                             target_type_selected, ".csv"))
          )
        }
      }
      
      # Export de clusters para ChemCluster usando la misma máscara segura
      if (target_type_selected == "ChemCluster") {
        clusternames <- as.character(gsea_p$target[which(sig_mask)])
        if (length(clusternames) != 0) {
          topclusterlist <- cmpdSets[clusternames]
          clusterdf <- as.data.frame(sapply(topclusterlist, toString))
          clusterdf$cluster <- clusternames
          clusterdf$pval <- gsea_p$padj[match(clusternames, gsea_p$target)]
          colnames(clusterdf)[1] <- "drugs.in.cluster"
          readr::write_csv(clusterdf, file = file.path(enrichFolder.n, "drugstructureclusters.csv"))
        }
      }
      
      message("Done for ", target_type_selected,
              " for ", sum(sig_mask, na.rm = TRUE), " genesets")
    } else {
      message("No significant enrichment found for ", target_type_selected)
    }
  }
}
