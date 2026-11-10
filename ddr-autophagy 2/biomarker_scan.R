#!/usr/bin/env Rscript
## biomarker_scan.R ----------------------------------------------------------
## WHAT MAKES A CELL NEED VPS34?
##
##   Rscript biomarker_scan.R [/path/to/depmap] [FOCAL_GENE]
##   Rscript biomarker_scan.R ~/depmap PIK3C3        # default
##
## Rather than asking "which pathway is VPS34 in" (it is a hub, so there may be
## no clean answer), this asks what CONTEXT creates the dependency. That is
## simultaneously the mechanism and the patient-selection biomarker.
##
## Outputs
##   1. where MDA-MB-231 and MDA-MB-468 sit in the dependency distribution
##   2. lineage differences in PIK3C3 dependency
##   3. mutations associated with PIK3C3 dependency
##   4. genome-wide expression correlates of PIK3C3 dependency
##   5. the same for ATG7, as a contrast
##
## MEMORY: the expression matrix is read in column chunks and discarded as it
## goes, so peak usage stays near 100 MB rather than >1 GB.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))
setDTthreads(1)

args     <- commandArgs(trailingOnly = TRUE)
DATA_DIR <- path.expand(if (length(args) >= 1) args[1] else "~/depmap")
FOCAL    <- if (length(args) >= 2) args[2] else "PIK3C3"
CONTRAST <- "ATG7"
OUT_DIR  <- "results_biomarker"
CHUNK    <- 2000

here <- tryCatch({
  sf <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(sf)) dirname(normalizePath(sub("^--file=", "", sf[1]))) else getwd()
}, error = function(e) getwd())
.src <- function(f) {
  for (p in c(file.path(here, "R", f), file.path(here, f)))
    if (file.exists(p)) { source(p); return(invisible(p)) }
  stop("Cannot find ", f)
}
.src("01_genesets.R")

if (!dir.exists(DATA_DIR)) stop("Data directory not found: ", DATA_DIR)
dir.create(OUT_DIR, showWarnings = FALSE)

clean <- function(x) trimws(sub("\\s*\\(\\d+\\)$", "", x))

## Coerce a data.table of expression columns to a plain double matrix.
## fread can return integer64 (large counts), character (a stray non-numeric
## token anywhere in a column), or factor -- and as.matrix() on any of those
## produces something rowSums() refuses. Forcing column-wise as.numeric is
## slower but survives all of them.
as_num_matrix <- function(dt) {
  bad <- names(dt)[!sapply(dt, is.numeric)]
  if (length(bad))
    message(sprintf("       [coerce] %d non-numeric column(s), e.g. %s",
                    length(bad), paste(head(bad, 3), collapse = ", ")))
  m <- suppressWarnings(as.matrix(dt[, lapply(.SD, as.numeric)]))
  storage.mode(m) <- "double"
  m
}

## ---------------------------------------------------------------------------
## 1. dependency + covariates (light read)
## ---------------------------------------------------------------------------
GE_FILE <- file.path(DATA_DIR, "CRISPRGeneEffect.csv")
if (!file.exists(GE_FILE)) stop("Missing CRISPRGeneEffect.csv in ", DATA_DIR)

message("[1/5] reading dependency columns ...")
hdr <- names(fread(GE_FILE, nrows = 0)); sym <- clean(hdr)
want <- unique(c(FOCAL, CONTRAST, AUTOPHAGY_PRIMARY, "PIK3R4", "BECN1", "RB1CC1"))
iw <- which(sym %in% want)
if (!length(iw)) stop("None of the focal genes found in the file.")
set.seed(42)
ir <- sort(sample(setdiff(seq_along(hdr)[-1], iw), min(1500, length(hdr) - 2)))
GE <- fread(GE_FILE, select = sort(unique(c(1L, iw, ir))), showProgress = FALSE)
ids <- as.character(GE[[1]]); GE[, (1) := NULL]
M <- as.matrix(GE); rm(GE); invisible(gc())
colnames(M) <- make.unique(clean(colnames(M))); rownames(M) <- ids
qc <- rowMeans(M, na.rm = TRUE)
if (!FOCAL %in% colnames(M)) stop(FOCAL, " not found.")

md <- fread(file.path(DATA_DIR, "Model.csv"), showProgress = FALSE)
keep <- intersect(c("ModelID","StrippedCellLineName","OncotreeLineage",
                    "OncotreeSubtype","PrimaryOrMetastasis"), names(md))
md <- md[, ..keep]

## NOTE: the dependency column is called `gedep`, not `dep`. Naming it `dep`
## inside a data.table also called `dep` makes `dep[!is.na(dep)]` ambiguous --
## data.table may resolve the inner `dep` to the table rather than the column,
## which silently returns the wrong rows instead of erroring.
dep <- data.table(ModelID = ids, gedep = as.numeric(M[, FOCAL]), mean_ge = qc)
if (CONTRAST %in% colnames(M)) dep[, dep_contrast := as.numeric(M[, CONTRAST])]
dep <- merge(dep, md, by = "ModelID")
dep <- dep[is.finite(gedep)]
message(sprintf("       %d lines with %s dependency", nrow(dep), FOCAL))
if (!nrow(dep)) stop("No lines with a usable dependency value.")

## ---------------------------------------------------------------------------
## 2. where do YOUR cell lines sit?
## ---------------------------------------------------------------------------
cat("\n==========================================================\n")
cat(sprintf(" WHERE YOUR LINES SIT IN THE %s DEPENDENCY DISTRIBUTION\n", FOCAL))
cat("==========================================================\n\n")
dep[, pct_all := 100 * frank(gedep) / .N]
dep[OncotreeLineage == "Breast", pct_breast := 100 * frank(gedep) / .N]

mine <- dep[grepl("MDAMB231|MDAMB468", StrippedCellLineName, ignore.case = TRUE)]
if (nrow(mine)) {
  print(mine[, .(StrippedCellLineName, gene_effect = round(gedep, 3),
                 dep_ATG7 = if ("dep_contrast" %in% names(mine))
                              round(dep_contrast, 3) else NA_real_,
                 pct_all_lines = round(pct_all, 1),
                 pct_breast = round(pct_breast, 1))])
  cat("\npercentile 1 = most dependent line. Below ~15th percentile means your\n")
  cat("line is in the dependent tail, so the mechanism may be specific to it\n")
  cat("rather than general. Above ~40th means it is fairly typical.\n")
}
cat(sprintf("\nAll lines: median %s gene effect = %+.3f; %.1f%% below -0.5\n",
            FOCAL, median(dep$gedep), 100 * mean(dep$gedep < -0.5)))

## ---------------------------------------------------------------------------
## 3. lineage
## ---------------------------------------------------------------------------
cat("\n==========================================================\n")
cat(" LINEAGE: which tumour types depend on it?\n")
cat("==========================================================\n\n")
lin <- dep[, .(n = .N, mean_dep = mean(gedep),
               frac_dependent = mean(gedep < -0.5)),
           by = OncotreeLineage][n >= 15][order(mean_dep)]
print(head(lin, 20), digits = 3)
fwrite(lin, file.path(OUT_DIR, sprintf("%s_by_lineage.csv", FOCAL)))

## ---------------------------------------------------------------------------
## 4. mutations
## ---------------------------------------------------------------------------
MUT_FILE <- file.path(DATA_DIR, "OmicsSomaticMutationsMatrixDamaging.csv")
if (file.exists(MUT_FILE)) {
  cat("\n==========================================================\n")
  cat(" MUTATIONS associated with the dependency\n")
  cat("==========================================================\n")
  message("[2/5] scanning mutations ...")
  MU <- fread(MUT_FILE, showProgress = FALSE)
  mid <- as.character(MU[[1]]); MU[, (1) := NULL]
  Mm <- as.matrix(MU); rm(MU); invisible(gc())
  colnames(Mm) <- make.unique(clean(colnames(Mm))); rownames(Mm) <- mid
  common <- intersect(dep$ModelID, rownames(Mm))
  message(sprintf("       ID overlap: %d of %d dependency lines, %d mutation lines",
                  length(common), nrow(dep), nrow(Mm)))
  if (length(common) < 100) {
    ## Diagnose rather than crash: a near-empty overlap means the ID columns
    ## were parsed differently, not that the biology is uninteresting.
    message("[!] too few shared cell lines -- skipping the mutation scan.")
    message("    dependency IDs look like: ", paste(head(dep$ModelID, 3), collapse = ", "))
    message("    mutation   IDs look like: ", paste(head(rownames(Mm), 3), collapse = ", "))
  } else {
    d2 <- dep[match(common, ModelID)]; Mm <- Mm[common, , drop = FALSE]
    Mb <- Mm > 0
    nmut <- colSums(Mb, na.rm = TRUE)
    Mb <- Mb[, nmut >= 20 & nmut <= (nrow(Mb) - 20), drop = FALSE]
    message(sprintf("       %d genes mutated in 20..n-20 of %d lines",
                    ncol(Mb), nrow(Mb)))
    if (!ncol(Mb)) {
      message("[!] no gene passes the frequency filter -- skipping.")
    } else {
      d2[, lin_f := factor(ifelse(OncotreeLineage %in%
            names(which(table(OncotreeLineage) >= 5)), OncotreeLineage, "Other"))]
      X0 <- if (nlevels(d2$lin_f) >= 2)
              model.matrix(~ lin_f + mean_ge, data = d2) else
              model.matrix(~ mean_ge, data = d2)
      yr <- residuals(lm(d2$gedep ~ X0 - 1))
      res <- rbindlist(lapply(colnames(Mb), function(g) {
        x <- as.numeric(Mb[, g])
        if (sum(x == 1) < 10 || sum(x == 0) < 10) return(NULL)
        data.table(gene = g, n_mut = sum(x),
                   delta = mean(yr[x == 1]) - mean(yr[x == 0]),
                   p = tryCatch(t.test(yr[x == 1], yr[x == 0])$p.value,
                                error = function(e) NA_real_))
      }))
      if (nrow(res)) {
        res[, fdr := p.adjust(p, "BH")]
        setorder(res, p, na.last = TRUE)
        cat("\ndelta < 0 => mutation makes cells MORE dependent\n\n")
        print(head(res, 15), digits = 3)
        fwrite(res, file.path(OUT_DIR, sprintf("%s_mutation_assoc.csv", FOCAL)))
      }
    }
    rm(Mb)
  }
  rm(Mm); invisible(gc())
}

## ---------------------------------------------------------------------------
## 5. expression scan, chunked
## ---------------------------------------------------------------------------
EX_TPM <- file.path(DATA_DIR, "OmicsExpressionProteinCodingGenesTPMLogp1.csv")
EX_RAW <- file.path(DATA_DIR,
            "OmicsExpressionRawReadCountHumanProteinCodingGenesStranded.csv")
EX <- if (file.exists(EX_TPM)) EX_TPM else if (file.exists(EX_RAW)) EX_RAW else NA
IS_RAW <- !is.na(EX) && EX == EX_RAW

if (!is.na(EX)) {
  cat("\n==========================================================\n")
  cat(" EXPRESSION correlates of the dependency\n")
  cat("==========================================================\n")
  message(sprintf("[3/5] expression file: %s", basename(EX)))
  if (IS_RAW) message("       raw counts -> converting to log2(CPM+1)")

  ## DepMap expression files are PROFILE-level, not model-level: the leading
  ## columns are metadata (row index, SequencingID, ModelConditionID, ModelID,
  ## IsDefaultEntryForMC, IsDefaultEntryForModel) and genes start after them.
  ## There can be several sequencing profiles per cell line, so we keep only
  ## each model's default profile -- otherwise a line with three profiles would
  ## be counted three times and weighted accordingly.
  ehdr   <- names(fread(EX, nrows = 0))
  ncol_e <- length(ehdr)

  ## genes are the columns named "SYMBOL (ENTREZ)"; everything before is metadata
  gcols <- which(grepl("\\(\\d+\\)$", ehdr))
  if (!length(gcols)) stop("No 'SYMBOL (ENTREZ)' columns found in ", basename(EX))
  GENE_START <- min(gcols)
  message(sprintf("       %d metadata columns, genes start at column %d",
                  GENE_START - 1, GENE_START))

  meta <- fread(EX, select = seq_len(GENE_START - 1), showProgress = FALSE)
  id_name <- intersect(c("ModelID", "DepMap_ID", "depmap_id"), names(meta))[1]
  if (is.na(id_name))
    stop("No ModelID column among: ", paste(names(meta), collapse = ", "))
  eids <- as.character(meta[[id_name]])

  ## keep one profile per model
  dcol <- intersect(c("IsDefaultEntryForModel", "IsDefaultEntryForMC"),
                    names(meta))[1]
  keep_row <- if (!is.na(dcol)) {
    kd <- meta[[dcol]]
    if (is.logical(kd)) kd else
      tolower(as.character(kd)) %in% c("true", "t", "1", "yes")
  } else rep(TRUE, nrow(meta))
  ok <- which(keep_row & !is.na(eids) & nzchar(eids))
  ok <- ok[!duplicated(eids[ok])]
  message(sprintf("       %d profiles -> %d unique models (%s)",
                  nrow(meta), length(ok),
                  if (!is.na(dcol)) paste("default by", dcol) else "no default flag"))

  ## library size: approximate from a random 3000-gene sample. Scale-only, and
  ## correlates >0.999 with the all-gene total, so it saves a whole extra pass.
  libsize <- NULL
  if (IS_RAW) {
    set.seed(11)
    s <- sort(sample(GENE_START:ncol_e, min(3000, ncol_e - GENE_START + 1)))
    tmp <- fread(EX, select = s, showProgress = FALSE)
    libsize <- rowSums(as_num_matrix(tmp), na.rm = TRUE); rm(tmp); invisible(gc())
    libsize[!is.finite(libsize) | libsize <= 0] <- NA
    message(sprintf("       library sizes: median %.3g, %d usable of %d lines",
                    median(libsize, na.rm = TRUE), sum(!is.na(libsize)), length(libsize)))
    if (all(is.na(libsize)))
      stop("Could not compute library sizes -- the expression file did not ",
           "parse as counts. Check its first few lines with:\n",
           "  head -2 ", EX, " | cut -c1-300")
  }

  common <- intersect(dep$ModelID, eids[ok])
  message(sprintf("       %d lines shared with dependency data", length(common)))
  if (length(common) < 100) {
    message("[!] too few shared lines -- skipping the expression scan.")
    message("    dependency IDs: ", paste(head(dep$ModelID, 3), collapse = ", "))
    message("    expression IDs: ", paste(head(eids[ok], 3), collapse = ", "))
    EX <- NA
  }
}
if (!is.na(EX)) {
  d3 <- dep[match(common, ModelID)]
  ## rowsel indexes rows of the FULL file, restricted to default profiles
  rowsel <- ok[match(common, eids[ok])]

  d3[, lin_f := factor(ifelse(OncotreeLineage %in%
        names(which(table(OncotreeLineage) >= 5)), OncotreeLineage, "Other"))]
  X0 <- if (nlevels(d3$lin_f) >= 2)
          model.matrix(~ lin_f + mean_ge, data = d3) else
          model.matrix(~ mean_ge, data = d3)
  yr <- residuals(lm(d3$gedep ~ X0 - 1))
  yr2 <- if ("dep_contrast" %in% names(d3) && !all(is.na(d3$dep_contrast)))
           residuals(lm(d3$dep_contrast ~ X0 - 1)) else NULL
  ls_sub <- if (IS_RAW) libsize[rowsel] else NULL

  starts <- seq(GENE_START, ncol_e, by = CHUNK)
  out <- vector("list", length(starts))
  for (k in seq_along(starts)) {
    cols <- starts[k]:min(starts[k] + CHUNK - 1, ncol_e)
    E <- as_num_matrix(fread(EX, select = cols, showProgress = FALSE))[rowsel, , drop = FALSE]
    if (IS_RAW) E <- log2(sweep(E, 1, ls_sub / 1e6, "/") + 1)
    sdv <- apply(E, 2, sd, na.rm = TRUE)
    kp <- which(is.finite(sdv) & sdv > 0.1)
    if (length(kp)) {
      Ek <- E[, kp, drop = FALSE]
      Er <- Ek - X0 %*% qr.coef(qr(X0), Ek)        # residualize expression too
      r1 <- as.numeric(cor(yr, Er, use = "pairwise.complete.obs"))
      dt <- data.table(gene = clean(colnames(E))[kp], r_focal = r1)
      if (!is.null(yr2))
        dt[, r_contrast := as.numeric(cor(yr2, Er, use = "pairwise.complete.obs"))]
      out[[k]] <- dt
    }
    rm(E); invisible(gc())
    message(sprintf("       chunk %d/%d", k, length(starts)))
  }
  ex <- rbindlist(out, fill = TRUE)
  n <- length(common)
  ex[, p_focal := 2 * pt(-abs(r_focal * sqrt((n - 2) / (1 - r_focal^2))), n - 2)]
  ex[, fdr_focal := p.adjust(p_focal, "BH")]
  setorder(ex, r_focal)

  cat(sprintf("\n--- expression MOST associated with NEEDING %s ---\n", FOCAL))
  cat("(negative r: high expression goes with more negative gene effect,\n")
  cat(" i.e. high expression marks cells that DEPEND on it)\n\n")
  print(head(ex, 25), digits = 3)
  cat(sprintf("\n--- expression associated with NOT needing %s ---\n", FOCAL))
  print(tail(ex, 15), digits = 3)
  fwrite(ex, file.path(OUT_DIR, sprintf("%s_expression_correlates.csv", FOCAL)))

  if ("r_contrast" %in% names(ex)) {
    cat(sprintf("\n--- specificity check: %s vs %s ---\n", FOCAL, CONTRAST))
    cat(sprintf("correlation between the two correlate profiles: %+.3f\n",
                cor(ex$r_focal, ex$r_contrast, use = "complete.obs")))
    cat("Near zero => different contexts drive the two dependencies, so your\n")
    cat("biomarker is specific to VPS34 rather than to autophagy generally.\n")
  }
}

cat(sprintf("\nFiles written to %s/\n", OUT_DIR))
cat("\nCaution: these are associations across cell lines, not causes. Treat the\n")
cat("top hits as candidate biomarkers to be tested, not as mechanism.\n\n")
