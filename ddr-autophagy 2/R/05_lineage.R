## 05_lineage.R --------------------------------------------------------------
## Translate per-line results into a ranked list of TUMOUR TYPES -- the actual
## deliverable ("which tumours would benefit").
##
## The composite priority score is reported alongside all of its components and
## never on its own. A single number that nobody can decompose is not evidence,
## and the components frequently disagree in informative ways.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

#' Per-lineage summary of autophagy vulnerability, DDR burden and their coupling.
#'
#' @param autophagy_score data.table(ModelID, score) vulnerability-oriented
#' @param features        DDR feature table
#' @param covars          data.table(ModelID, lineage, mean_gene_effect)
#' @param gene_effect     Chronos matrix (for within-lineage coupling)
#' @param coupling_genes  DDR genes whose dependency defines "coupling".
#'                        Default is the replication-stress checkpoint axis,
#'                        which is where the mechanism predicts the interaction.
#' @param min_n           minimum cell lines per lineage
rank_lineages <- function(autophagy_score, features, covars, gene_effect,
                          coupling_genes = c("ATR", "CHEK1", "WEE1", "PARP1", "POLQ"),
                          min_n = 10) {

  ## Carry raw_mean_ge through as well: `auto` is in z-score units (good for
  ## ranking) but "fraction of lines that are dependent" must be judged in
  ## Chronos units against the conventional -0.5 cutoff. Mixing the two is an
  ## easy and completely silent error.
  d <- merge(autophagy_score[, .(ModelID, auto = score, raw_ge = raw_mean_ge)],
             covars, by = "ModelID")
  d <- merge(d, features, by = "ModelID", all.x = TRUE)

  cg <- intersect(coupling_genes, colnames(gene_effect))
  if (length(cg)) {
    ## Vulnerability orientation for the DDR side too, so that a POSITIVE
    ## coupling correlation = hypothesis-consistent. This is the one place we
    ## flip the DDR sign; it makes the lineage table readable.
    ddr_vuln <- -rowMeans(gene_effect[, cg, drop = FALSE], na.rm = TRUE)
    d[, ddr_vuln := ddr_vuln[match(ModelID, names(ddr_vuln))]]
  } else {
    d[, ddr_vuln := NA_real_]
  }

  keep <- d[, .N, by = lineage][N >= min_n, lineage]
  message(sprintf("[lineage] %d lineages with n >= %d (of %d total)",
                  length(keep), min_n, uniqueN(d$lineage)))
  d <- d[lineage %in% keep]

  ## Guarantee every column the aggregation references actually exists.
  ## Relying on dead `if` branches inside a data.table j-expression works but
  ## is fragile across data.table versions; materializing NA columns is not.
  for (col in c("hrr_any_mut", "mut_TP53", "mut_ATM", "amp_CCNE1",
                "repstress_score")) {
    if (!col %in% names(d)) d[[col]] <- NA_real_
  }

  summ <- d[, {
    n <- .N
    ## within-lineage coupling: do the autophagy-dependent lines in THIS
    ## lineage also depend on the DDR axis?
    cpl <- if (sum(!is.na(auto) & !is.na(ddr_vuln)) >= 8) {
      suppressWarnings(cor(auto, ddr_vuln, method = "spearman",
                           use = "pairwise.complete.obs"))
    } else NA_real_
    .(
      n_lines            = n,
      autophagy_vuln     = mean(auto, na.rm = TRUE),
      autophagy_vuln_sd  = sd(auto, na.rm = TRUE),
      ## judged on Chronos units: gene effect < -0.5 = conventional "dependent"
      frac_autophagy_dep = mean(raw_ge < CHRONOS_DEP_THRESHOLD, na.rm = TRUE),
      mean_chronos       = mean(raw_ge, na.rm = TRUE),
      hrr_mut_frac       = mean(hrr_any_mut,     na.rm = TRUE),
      tp53_mut_frac      = mean(mut_TP53,        na.rm = TRUE),
      atm_mut_frac       = mean(mut_ATM,         na.rm = TRUE),
      ccne1_amp_frac     = mean(amp_CCNE1,       na.rm = TRUE),
      repstress          = mean(repstress_score, na.rm = TRUE),
      ddr_vuln_mean      = mean(ddr_vuln, na.rm = TRUE),
      coupling_rho       = cpl
    )
  }, by = lineage]

  ## ---- composite priority ------------------------------------------------
  ## Equal-weighted sum of z-scored components. Equal weighting is a choice,
  ## not a result: there is no principled basis for weighting these against
  ## each other, so we keep it flat and transparent rather than tuning it.
  z <- function(x) {
    if (all(is.na(x))) return(rep(NA_real_, length(x)))
    as.numeric(scale(x))
  }
  summ[, z_autophagy := z(autophagy_vuln)]
  summ[, z_ddr       := z(ddr_vuln_mean)]
  summ[, z_repstress := z(repstress)]
  summ[, z_coupling  := z(coupling_rho)]

  comp <- c("z_autophagy", "z_ddr", "z_repstress", "z_coupling")
  comp <- comp[sapply(comp, function(c) !all(is.na(summ[[c]])))]
  message(sprintf("[lineage] priority score built from: %s", paste(comp, collapse = ", ")))

  summ[, priority := rowMeans(as.matrix(.SD), na.rm = TRUE), .SDcols = comp]
  summ[, n_components := rowSums(!is.na(as.matrix(.SD))), .SDcols = comp]
  setorder(summ, -priority)
  summ[]
}

## Conventional DepMap cutoff for calling a line "dependent" on a gene, in
## Chronos units. Named once so it cannot drift between modules.
CHRONOS_DEP_THRESHOLD <- -0.5

#' Write a compact, human-readable report of the lineage ranking.
report_lineages <- function(summ, top_n = 15, file = NULL) {
  cols <- intersect(c("lineage", "n_lines", "priority", "mean_chronos",
                      "frac_autophagy_dep", "ddr_vuln_mean", "coupling_rho",
                      "hrr_mut_frac", "ccne1_amp_frac", "repstress",
                      "n_components"),
                    names(summ))
  tbl <- head(summ[, ..cols], top_n)
  cat("\n================ LINEAGE PRIORITY RANKING ================\n")
  cat("priority = mean of z-scored components (higher = stronger\n")
  cat("rationale for DDR-agent + autophagy-blockade combination)\n\n")
  print(tbl, digits = 3)
  cat("\nReminder: this ranks HYPOTHESIS STRENGTH from cell-line data.\n")
  cat("It is not evidence of clinical benefit.\n")
  cat("==========================================================\n\n")
  if (!is.null(file)) {
    data.table::fwrite(summ, file)
    message(sprintf("[lineage] full table written to %s", file))
  }
  invisible(tbl)
}
