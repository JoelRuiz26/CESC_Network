# ======================================================
# Simple in-silico validation 
# ======================================================

# ---------- 0) Base folder & file layout ----------
base_dir  <- "~/CESC_Network/6_OCTAD/6_5_Validation"   

# Output folders
output_A7 <- file.path(base_dir, "Validation_A7")
output_A9 <- file.path(base_dir, "Validation_A9")

# ---------- 1) Source function ----------

setwd(base_dir)
source("topLineEval_custom.R")

# Optional safety check: confirm the function exists
if (!exists("topLineEval")) {
  stop("topLineEval() was not found after sourcing 'topLineEval_custom.R'. ",
       "Make sure that file is located in: ", base_dir)
}

# ---------- 2) Read inputs ----------
RGEs_A7 <- readRDS("~/CESC_Network/6_OCTAD/6_3_RES/A7/RES_A7_common3_collapsed_FULL.rds")
result_A7 <- RGEs_A7[, c("pert_iname", "sRGES")]

RGEs_A9 <- readRDS("~/CESC_Network/6_OCTAD/6_3_RES/A9/RES_A9_common3_collapsed_FULL.rds")
result_A9 <- RGEs_A9[, c("pert_iname", "sRGES")]

# Similar cell lines (coerce directly to unique character vector)
lineas_similares_A7_p <- unique(as.character(unlist(
  readRDS("~/CESC_Network/6_OCTAD/6_3_RES/A7/SimilarCellLines_A7_medcor_gt0.30.rds"),
  use.names = FALSE
)))

lineas_similares_A9_p <- unique(as.character(unlist(
  readRDS("~/CESC_Network/6_OCTAD/6_3_RES/A9/SimilarCellLines_A9_medcor_gt0.30.rds"),
  use.names = FALSE
)))

clado_a7 <- c("HELA","SW756","MS751","ME-180","ME180","C4-I","C4-1","C4-II","C4-2")
clado_a9 <- c("SIHA", "CASKI", "SISO","HELA") #"HELA" is HPV18, but i added here too

lineas_similares_A7 <- unique(c(lineas_similares_A7_p, clado_a7))
lineas_similares_A9 <- unique(c(lineas_similares_A9_p, clado_a9))


# ---------- 4) Minimal runner (no aliases, no filtering) ----------
# For each line, call topLineEval and catch errors so the batch continues.
run_validation <- function(clade, mysRGES, lines_vec, out_dir) {
  # Returns a simple log data.frame with a status per line.
  log <- data.frame(
    clade   = character(),
    line    = character(),
    status  = character(),
    message = character(),
    stringsAsFactors = FALSE
  )
  if (length(lines_vec) == 0) {
    message("No similar lines provided for clade ", clade, ".")
    return(log)
  }
  
  for (ln in lines_vec) {
    message(sprintf("[%s] Validating %s ...", clade, ln))
    out_line_dir <- file.path(out_dir, ln)
    dir.create(out_line_dir, recursive = TRUE, showWarnings = FALSE)
    
    st  <- "success"
    msg <- "OK"
    
    # Call your custom evaluator. If it returns a data.frame, we save it.
    tryCatch({
      ret <- topLineEval(
        topline      = ln,
        mysRGES      = mysRGES,
        outputFolder = out_line_dir
      )
      if (is.data.frame(ret)) {
        utils::write.csv(
          ret,
          file = file.path(out_line_dir, paste0("summary_", ln, ".csv")),
          row.names = FALSE
        )
      }
    }, error = function(e) {
      st  <<- "error"
      msg <<- conditionMessage(e)
      message(sprintf("  -> ERROR with %s: %s", ln, msg))
    })
    
    log <- rbind(
      log,
      data.frame(clade = clade, line = ln, status = st, message = msg,
                 stringsAsFactors = FALSE)
    )
  }
  log
}

# ---------- 5) Run by clade ----------
log_a7  <- run_validation("A7", result_A7, lineas_similares_A7, output_A7)
log_a9  <- run_validation("A9", result_A9, lineas_similares_A9, output_A9)
log_all <- rbind(log_a7, log_a9)

save.image(file.path(base_dir, "6_5_2_Image_Validation.RData"))

