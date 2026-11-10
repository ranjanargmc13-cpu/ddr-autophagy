## 02_load_depmap.R ----------------------------------------------------------
## Robust loaders for DepMap flat files.
##
## DepMap renames files and reshuffles columns between quarterly releases, so
## every loader validates its input and fails LOUDLY rather than silently
## returning something mis-shaped. A silent schema change is the single most
## common way this kind of analysis produces a confident wrong answer.
##
## SIGN CONVENTION, stated once and relied on everywhere downstream:
##   Chronos gene effect: 0 = no effect, NEGATIVE = loss of fitness.
##   So MORE NEGATIVE = MORE DEPENDENT.
##   We therefore define  vulnerability = -1 * gene_effect  wherever a
##   "higher means more dependent" orientation is easier to read, and we always
##   label which orientation a column is in.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

## ---------------------------------------------------------------------------
## Expected files. Download from:
##   https://depmap.org/portal/data_page/?tab=allData
## Filenames below match DepMap 24Q4 / 25Q2. If your release differs, change
## them here -- this is the only place they are named.
## ---------------------------------------------------------------------------
DEPMAP_FILES <- list(
  gene_effect = "CRISPRGeneEffect.csv",
  model       = "Model.csv",
  mut_damaging= "OmicsSomaticMutationsMatrixDamaging.csv",
  cn_gene     = "OmicsCNGene.csv",
  expression  = "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
)

#' Strip DepMap's "SYMBOL (ENTREZID)" column names down to "SYMBOL".
#' Duplicated symbols are disambiguated rather than dropped, so you always
#' know if a gene appeared twice.
.clean_gene_names <- function(x) {
  sym <- sub("\\s*\\(\\d+\\)$", "", x)
  sym <- trimws(sym)
  dup <- duplicated(sym)
  if (any(dup)) {
    warning(sprintf("[load] %d duplicated gene symbols after cleaning; suffixing with .dupN",
                    sum(dup)))
    sym <- make.unique(sym, sep = ".dup")
  }
  sym
}

#' Read a DepMap wide matrix file (rows = ModelID, cols = genes).
#' @return numeric matrix with ModelID rownames and clean gene colnames.
read_depmap_matrix <- function(path, value_name = "value", max_cols = NULL) {
  if (!file.exists(path)) {
    stop(sprintf(paste0(
      "[load] Missing file: %s\n",
      "       Download it from https://depmap.org/portal/data_page/?tab=allData\n",
      "       and place it in your data directory."), path), call. = FALSE)
  }
  message(sprintf("[load] reading %s ...", basename(path)))
  dt <- data.table::fread(path, check.names = FALSE, showProgress = FALSE,
                          nrows = if (is.null(max_cols)) Inf else Inf)

  ## First column holds ModelID. DepMap sometimes ships it unnamed, sometimes
  ## as "V1", sometimes as "ModelID" / "depmap_id". Handle all of them.
  first <- names(dt)[1]
  if (!first %in% c("ModelID", "depmap_id", "DepMap_ID")) {
    if (first == "" || grepl("^V1$", first)) {
      data.table::setnames(dt, 1L, "ModelID")
    } else if (is.character(dt[[1]]) && all(grepl("^ACH-", head(dt[[1]], 20)))) {
      ## Column is clearly model IDs even though the header says otherwise.
      data.table::setnames(dt, 1L, "ModelID")
    } else {
      stop(sprintf("[load] %s: cannot identify the ModelID column (first col = '%s')",
                   basename(path), first), call. = FALSE)
    }
  } else {
    data.table::setnames(dt, 1L, "ModelID")
  }

  ids <- as.character(dt[["ModelID"]])
  if (!all(grepl("^ACH-", ids))) {
    warning(sprintf("[load] %s: %d/%d row IDs do not look like ACH-###### model IDs",
                    basename(path), sum(!grepl("^ACH-", ids)), length(ids)))
  }

  dt[, ModelID := NULL]
  gnames <- .clean_gene_names(names(dt))

  m <- as.matrix(dt)
  storage.mode(m) <- "double"
  dimnames(m) <- list(ids, gnames)

  message(sprintf("[load]   -> %d models x %d genes (%s)",
                  nrow(m), ncol(m), value_name))
  m
}

#' Read Model.csv (the cell line annotation table).
read_model_table <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("[load] Missing file: %s", path), call. = FALSE)
  }
  message("[load] reading Model.csv ...")
  mt <- data.table::fread(path, check.names = FALSE, showProgress = FALSE)

  if (!"ModelID" %in% names(mt)) {
    stop("[load] Model.csv has no ModelID column -- unexpected schema.", call. = FALSE)
  }
  ## OncotreeLineage is the lineage field in recent releases; older ones used
  ## "lineage" or "primary_disease".
  lin_col <- intersect(c("OncotreeLineage", "lineage", "primary_disease"), names(mt))[1]
  if (is.na(lin_col)) {
    stop("[load] Model.csv has no recognizable lineage column.", call. = FALSE)
  }
  if (lin_col != "OncotreeLineage") {
    message(sprintf("[load]   using '%s' as lineage (OncotreeLineage absent)", lin_col))
    mt[, OncotreeLineage := get(lin_col)]
  }

  keep <- intersect(c("ModelID", "CellLineName", "StrippedCellLineName",
                      "OncotreeLineage", "OncotreePrimaryDisease",
                      "OncotreeSubtype", "PrimaryOrMetastasis",
                      "GrowthPattern", "Sex", "Age"), names(mt))
  mt <- mt[, ..keep]
  message(sprintf("[load]   -> %d models, %d lineages",
                  nrow(mt), length(unique(mt$OncotreeLineage))))
  mt
}

#' Load everything the analysis needs.
#' @param data_dir directory containing the DepMap CSVs
#' @param load_expression expression is the largest file (~1-2 GB). Set FALSE
#'        to skip it; the replication-stress signature will then be unavailable.
load_depmap <- function(data_dir, load_expression = TRUE) {
  stopifnot(dir.exists(data_dir))
  p <- function(k) file.path(data_dir, DEPMAP_FILES[[k]])

  out <- list()
  out$gene_effect <- read_depmap_matrix(p("gene_effect"), "Chronos gene effect")
  out$model       <- read_model_table(p("model"))

  out$mut <- tryCatch(
    read_depmap_matrix(p("mut_damaging"), "damaging mutation count"),
    error = function(e) { warning("[load] damaging-mutation matrix unavailable: ",
                                  conditionMessage(e)); NULL })

  out$cn <- tryCatch(
    read_depmap_matrix(p("cn_gene"), "relative copy number"),
    error = function(e) { warning("[load] copy-number matrix unavailable: ",
                                  conditionMessage(e)); NULL })

  out$expr <- NULL
  if (load_expression) {
    out$expr <- tryCatch(
      read_depmap_matrix(p("expression"), "log2(TPM+1)"),
      error = function(e) { warning("[load] expression matrix unavailable: ",
                                    conditionMessage(e)); NULL })
  }

  ## ---- Screen-quality / global-essentiality covariate ---------------------
  ## Per-line mean gene effect across ALL genes. This absorbs differences in
  ## screen efficiency, Cas9 activity and overall fitness sensitivity. Omitting
  ## it is the most common source of spurious co-dependency: badly behaved
  ## screens make every gene look like a dependency at once, which manufactures
  ## correlation between any two gene pairs.
  out$line_qc <- data.table(
    ModelID = rownames(out$gene_effect),
    mean_gene_effect = rowMeans(out$gene_effect, na.rm = TRUE),
    n_measured = rowSums(!is.na(out$gene_effect))
  )

  out
}

#' Subset a gene set to genes actually present, reporting what is missing.
#' Never silently drops genes -- a quietly shrinking gene set changes the
#' meaning of a score without changing its name.
match_genes <- function(gene_set, available, set_name = "set") {
  hit  <- intersect(gene_set, available)
  miss <- setdiff(gene_set, available)
  if (length(miss)) {
    message(sprintf("[match] %s: %d/%d present. Missing: %s",
                    set_name, length(hit), length(gene_set),
                    paste(miss, collapse = ", ")))
  } else {
    message(sprintf("[match] %s: all %d genes present.", set_name, length(hit)))
  }
  if (!length(hit)) stop(sprintf("[match] %s: no genes found at all.", set_name),
                         call. = FALSE)
  hit
}
