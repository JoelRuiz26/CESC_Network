topLineEval <- function(topline = NULL, mysRGES = NULL, outputFolder = NULL) {
  if (missing(mysRGES)) stop("sRGES signature input not found")
  if (length(topline) == 0) stop("No cell lines found in the input")
  if (is.null(outputFolder)) {
    outputFolder <- tempdir()
    message("outputFolder is NULL, writing output to tempdir()")
  }
  toplineName <- paste(topline, collapse = "_")
  cell.line.folder <- file.path(outputFolder, "CellLineEval")
  if (!dir.exists(cell.line.folder)) dir.create(cell.line.folder)
  
  mysRGES$pert_iname <- toupper(mysRGES$pert_iname)
  
  # IC50
  CTRPv2.ic50 <- data.table::dcast.data.table(
    octad.db::get_ExperimentHub_data("EH7264"),
    drugid ~ cellid, value.var = "ic50_recomputed", fun.aggregate = median
  )
  colnames(CTRPv2.ic50) <- gsub("[^0-9A-Za-z///' ]", "", colnames(CTRPv2.ic50))
  colnames(CTRPv2.ic50) <- toupper(colnames(CTRPv2.ic50))
  colnames(CTRPv2.ic50) <- gsub(" ", "", colnames(CTRPv2.ic50))
  CTRP.IC50 <- subset(as.data.frame(CTRPv2.ic50), select = c("DRUGID", topline))
  CTRP.IC50.m <- reshape2::melt(CTRP.IC50, id.vars = "DRUGID")
  CTRP.IC50.medianIC50 <- aggregate(CTRP.IC50.m[3], by = list(CTRP.IC50.m$DRUGID), FUN = median)
  colnames(CTRP.IC50.medianIC50) <- c("DRUGID", "medIC50")
  
  # Merge with mysRGES
  testdf <- merge(mysRGES, CTRP.IC50.medianIC50, by.x = "pert_iname", by.y = "DRUGID")
  IC50.cortest <- cor.test(testdf$sRGES, log10(testdf$medIC50))
  testdf$StronglyPredicted <- ifelse(testdf$sRGES < -0.2, "Yes", "No")
  
  # Save IC50 table
  write.table(testdf, file = file.path(cell.line.folder, paste0(topline, "_ic50_insilico_data.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Plot IC50
  p <- ggplot2::ggplot(testdf, aes(x = sRGES, y = log10(medIC50))) +
    ggplot2::geom_point(aes(color = StronglyPredicted, text = paste("Drug: ", pert_iname, "<br>sRGES: ", sRGES))) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5) +
    ggplot2::labs(x = "sRGES", y = "log(IC50)", title = "Top Line recomputed log(IC50) vs sRGES") +
    ggplot2::scale_color_discrete(name = "Strongly\nPredicted") +
    ggplot2::theme_minimal()
  
  ic50graph <- plotly::ggplotly(p, tooltip = "text")
  htmlwidgets::saveWidget(ic50graph,
                          file.path(cell.line.folder, paste0(topline, "_ic50_insilico_validation.html")),
                          selfcontained = FALSE
  )
  
  # AUC
  CTRPv2.auc <- data.table::dcast.data.table(
    octad.db::get_ExperimentHub_data("EH7264"),
    drugid ~ cellid, value.var = "auc_recomputed", fun.aggregate = median
  )
  colnames(CTRPv2.auc) <- gsub("[^0-9A-Za-z///' ]", "", colnames(CTRPv2.auc))
  colnames(CTRPv2.auc) <- toupper(colnames(CTRPv2.auc))
  colnames(CTRPv2.auc) <- gsub(" ", "", colnames(CTRPv2.auc))
  CTRP.auc <- subset(as.data.frame(CTRPv2.auc), select = c("DRUGID", topline))
  CTRP.auc.m <- reshape2::melt(CTRP.auc, id.vars = "DRUGID")
  CTRP.auc.medianauc <- aggregate(CTRP.auc.m[3], by = list(CTRP.auc.m$DRUGID), FUN = median)
  colnames(CTRP.auc.medianauc) <- c("DRUGID", "medauc")
  
  testdf2 <- merge(mysRGES, CTRP.auc.medianauc, by.x = "pert_iname", by.y = "DRUGID")
  testdf2$StronglyPredicted <- ifelse(testdf2$sRGES < -0.2, "Yes", "No")
  
  # Save AUC table
  write.table(testdf2, file = file.path(cell.line.folder, paste0(topline, "_auc_insilico_data.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Plot AUC
  p2 <- ggplot2::ggplot(testdf2, aes(x = sRGES, y = medauc)) +
    ggplot2::geom_point(aes(color = StronglyPredicted, text = paste("Drug: ", pert_iname, "<br>sRGES: ", sRGES))) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5) +
    ggplot2::labs(x = "sRGES", y = "AUC", title = "Top Line recomputed AUC vs sRGES") +
    ggplot2::scale_color_discrete(name = "Strongly\nPredicted") +
    ggplot2::theme_minimal()
  
  aucgraph <- plotly::ggplotly(p2, tooltip = "text")
  htmlwidgets::saveWidget(aucgraph,
                          file.path(cell.line.folder, paste0(topline, "_auc_insilico_validation.html")),
                          selfcontained = FALSE
  )
}
