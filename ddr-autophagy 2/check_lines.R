#!/usr/bin/env Rscript
## check_lines.R -------------------------------------------------------------
## Look up your own cell lines in DepMap before going to the bench.
##
##   Rscript check_lines.R                      # uses ~/depmap
##   Rscript check_lines.R /path/to/depmap      # or pass the folder
##
## Prints, for MDA-MB-231 and MDA-MB-468:
##   - autophagy and DDR gene dependencies (Chronos: negative = dependent)
##   - each line's rank among breast lines for the core ATG genes
##   - the 231 vs 468 difference, which is your bench prediction
## ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
DATA_DIR <- path.expand(if (length(args)) args[1] else "~/depmap")

## locate modules whether they are in R/ or flat alongside this script
here <- tryCatch({
  sf <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(sf)) dirname(normalizePath(sub("^--file=", "", sf[1]))) else getwd()
}, error = function(e) getwd())
.src <- function(f) {
  for (p in c(file.path(here, "R", f), file.path(here, f))) {
    if (file.exists(p)) { source(p); return(invisible(p)) }
  }
  stop("Cannot find ", f)
}
.src("01_genesets.R"); .src("02_load_depmap.R")

if (!dir.exists(DATA_DIR))
  stop("Data directory not found: ", DATA_DIR,
       "\nPass it explicitly: Rscript check_lines.R /full/path/to/depmap")

dm <- load_depmap(DATA_DIR, load_expression = FALSE)

## ---- find the lines -------------------------------------------------------
## DepMap strips punctuation in StrippedCellLineName, so MDA-MB-231 becomes
## MDAMB231. Match on both fields to be safe.
find_line <- function(pattern, label) {
  hit <- dm$model[grepl(pattern, StrippedCellLineName, ignore.case = TRUE) |
                  grepl(pattern, CellLineName, ignore.case = TRUE)]
  if (!nrow(hit)) { message("[!] not found in Model.csv: ", label); return(NULL) }
  if (nrow(hit) > 1) message("[!] multiple matches for ", label, "; using first")
  hit[1]
}
l231 <- find_line("MDAMB231|MDA-MB-231", "MDA-MB-231")
l468 <- find_line("MDAMB468|MDA-MB-468", "MDA-MB-468")

rows <- rbindlist(list(l231, l468), fill = TRUE)
if (!nrow(rows)) stop("Neither line found. Check Model.csv.")
cat("\n==================== LINES FOUND ====================\n")
print(rows[, intersect(c("ModelID","CellLineName","StrippedCellLineName",
                         "OncotreeLineage","OncotreeSubtype",
                         "PrimaryOrMetastasis"), names(rows)), with = FALSE])

ids <- rows$ModelID
ids <- ids[ids %in% rownames(dm$gene_effect)]
if (!length(ids)) stop("Lines are in Model.csv but absent from CRISPRGeneEffect.csv (not screened).")

## ---- dependency profile ---------------------------------------------------
genes <- intersect(c(AUTOPHAGY_PRIMARY, "PIK3C3", "BECN1", "RB1CC1", "ULK1",
                     "PARP1", "ATR", "CHEK1", "WEE1", "POLQ", "PRKDC",
                     "PTEN", "MTOR"),
                   colnames(dm$gene_effect))

prof <- dm$gene_effect[ids, genes, drop = FALSE]
rownames(prof) <- rows$StrippedCellLineName[match(ids, rows$ModelID)]

cat("\n============ CHRONOS GENE EFFECT (negative = dependent) ============\n")
print(round(t(prof), 3))

## ---- rank among breast lines ---------------------------------------------
br <- intersect(dm$model[OncotreeLineage == "Breast", ModelID], rownames(dm$gene_effect))
cat(sprintf("\n============ RANK AMONG %d BREAST LINES ============\n", length(br)))
cat("(rank 1 = most dependent of all breast lines)\n\n")

rk <- rbindlist(lapply(intersect(AUTOPHAGY_PRIMARY, colnames(dm$gene_effect)), function(g) {
  v <- sort(dm$gene_effect[br, g], na.last = NA)
  data.table(gene = g,
             n_breast = length(v),
             rank_231 = if (length(l231)) which(names(v) == l231$ModelID)[1] else NA_integer_,
             rank_468 = if (length(l468)) which(names(v) == l468$ModelID)[1] else NA_integer_,
             ge_231   = if (length(l231)) unname(dm$gene_effect[l231$ModelID, g]) else NA_real_,
             ge_468   = if (length(l468)) unname(dm$gene_effect[l468$ModelID, g]) else NA_real_)
}))
rk[, pct_231 := round(100 * rank_231 / n_breast, 1)]
rk[, pct_468 := round(100 * rank_468 / n_breast, 1)]
print(rk, digits = 3)

## ---- the prediction -------------------------------------------------------
if (length(l231) && length(l468)) {
  d <- mean(rk$ge_468 - rk$ge_231, na.rm = TRUE)
  cat("\n==================== PREDICTION ====================\n")
  cat(sprintf("mean ATG gene effect, MDA-MB-231 : %+.3f\n", mean(rk$ge_231, na.rm=TRUE)))
  cat(sprintf("mean ATG gene effect, MDA-MB-468 : %+.3f\n", mean(rk$ge_468, na.rm=TRUE)))
  cat(sprintf("difference (468 - 231)           : %+.3f\n\n", d))
  if (is.na(d)) {
    cat("Cannot compare -- missing data.\n")
  } else if (d > 0.05) {
    cat("468 is LESS autophagy-dependent than 231.\n")
    cat("Consistent with the PTEN/mTOR prediction: PTEN-null 468 should show\n")
    cat("weaker chloroquine-olaparib synergy. That is a falsifiable bench test.\n")
  } else if (d < -0.05) {
    cat("468 is MORE autophagy-dependent than 231.\n")
    cat("OPPOSITE to the PTEN/mTOR prediction. If your bench data then shows\n")
    cat("weaker synergy in 468, baseline dependency is not what drives synergy\n")
    cat("-- which argues for induced rather than basal autophagy. Also useful.\n")
  } else {
    cat("The two lines are essentially equivalent in baseline ATG dependency.\n")
    cat("So any synergy difference you measure cannot be attributed to baseline\n")
    cat("autophagy dependence, and points to induced flux instead.\n")
  }
  cat("\nNote: CRISPR knockout is chronic; chloroquine is acute. These are\n")
  cat("related but not interchangeable readouts.\n")
  cat("====================================================\n\n")
}
