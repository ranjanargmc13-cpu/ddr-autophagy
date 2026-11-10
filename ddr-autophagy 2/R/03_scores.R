## 03_scores.R ---------------------------------------------------------------
## Build the two sides of the hypothesis:
##   (a) a per-cell-line AUTOPHAGY DEPENDENCY score, and
##   (b) a per-cell-line DDR FEATURE matrix.
##
## The coherence diagnostics in build_dependency_score() are not optional
## decoration. Autophagy is largely dispensable for proliferation in nutrient-
## rich culture, so most ATG genes show gene effect near zero in most lines.
## If that is true in your release, a "mean ATG z-score" is averaging noise and
## every downstream p-value is meaningless. The diagnostics tell you whether a
## coherent dependency axis exists BEFORE you interpret anything.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

#' Build a per-line dependency score from a gene set, with coherence diagnostics.
#'
#' @param gene_effect numeric matrix, models x genes, Chronos (negative = dependent)
#' @param gene_set    character vector of symbols
#' @param set_name    label for messages
#' @param min_sd      genes with across-line SD below this carry no selective
#'                    signal; they are reported and excluded from the z-score.
#' @param orient      "vulnerability" (higher = more dependent, sign-flipped) or
#'                    "chronos" (raw, negative = more dependent)
#' @return list(score, per_gene, diagnostics)
build_dependency_score <- function(gene_effect, gene_set, set_name,
                                   min_sd = 0.05,
                                   orient = c("vulnerability", "chronos")) {
  orient <- match.arg(orient)
  genes  <- match_genes(gene_set, colnames(gene_effect), set_name)
  M      <- gene_effect[, genes, drop = FALSE]

  ## ---- per-gene summary -------------------------------------------------
  per_gene <- data.table(
    gene      = genes,
    mean_ge   = colMeans(M, na.rm = TRUE),
    sd_ge     = apply(M, 2, sd, na.rm = TRUE),
    ## fraction of lines that are meaningfully dependent. -0.5 is the
    ## conventional DepMap threshold for "dependent"; -1.0 approximates the
    ## median common-essential gene.
    frac_dep_0.5 = colMeans(M < -0.5, na.rm = TRUE),
    frac_dep_1.0 = colMeans(M < -1.0, na.rm = TRUE),
    n_measured   = colSums(!is.na(M))
  )

  flat <- per_gene$gene[per_gene$sd_ge < min_sd]
  if (length(flat)) {
    message(sprintf("[score] %s: %d gene(s) below SD %.3f, excluded from z-score: %s",
                    set_name, length(flat), min_sd, paste(flat, collapse = ", ")))
    genes <- setdiff(genes, flat)
    if (!length(genes))
      stop(sprintf("[score] %s: every gene is flat. No dependency axis exists; stop here.",
                   set_name), call. = FALSE)
    M <- gene_effect[, genes, drop = FALSE]
  }

  ## ---- z-score each gene across lines, then average ----------------------
  Z <- scale(M)                      # centers and scales per column (gene)
  score <- rowMeans(Z, na.rm = TRUE)

  ## ---- COHERENCE DIAGNOSTICS (the part that decides credibility) ---------
  cmat <- suppressWarnings(cor(M, use = "pairwise.complete.obs"))
  offdiag <- cmat[upper.tri(cmat)]
  mean_r  <- mean(offdiag, na.rm = TRUE)

  ## Cronbach's alpha on the gene set: does it behave as one construct?
  k <- length(genes)
  alpha <- if (k > 1 && !is.na(mean_r) && mean_r > -1/(k-1)) {
    (k * mean_r) / (1 + (k - 1) * mean_r)
  } else NA_real_

  ## Variance explained by PC1: a single dominant axis supports a single score.
  pc1_var <- tryCatch({
    Mc <- M[complete.cases(M), , drop = FALSE]
    if (nrow(Mc) > k) {
      pr <- prcomp(Mc, center = TRUE, scale. = TRUE)
      (pr$sdev[1]^2) / sum(pr$sdev^2)
    } else NA_real_
  }, error = function(e) NA_real_)

  diagnostics <- list(
    set_name          = set_name,
    n_genes_used      = k,
    genes_used        = genes,
    genes_excluded_flat = flat,
    mean_pairwise_r   = mean_r,
    cronbach_alpha    = alpha,
    pc1_var_explained = pc1_var,
    median_mean_ge    = median(per_gene$mean_ge, na.rm = TRUE),
    max_frac_dep_0.5  = max(per_gene$frac_dep_0.5, na.rm = TRUE)
  )

  message(sprintf(paste0(
    "[score] %s diagnostics:\n",
    "        genes used        : %d\n",
    "        mean pairwise r   : %+.3f\n",
    "        Cronbach alpha    : %s\n",
    "        PC1 var explained : %s\n",
    "        median mean GE    : %+.3f  (Chronos units)\n",
    "        max %% lines dep.  : %.1f%% (at GE < -0.5)"),
    set_name, k, mean_r,
    ifelse(is.na(alpha), "NA", sprintf("%+.3f", alpha)),
    ifelse(is.na(pc1_var), "NA", sprintf("%.1f%%", 100 * pc1_var)),
    diagnostics$median_mean_ge, 100 * diagnostics$max_frac_dep_0.5))

  if (!is.na(mean_r) && mean_r < 0.1) {
    warning(sprintf(paste0(
      "[score] %s: mean pairwise r = %.3f (< 0.10). The genes do NOT co-vary, ",
      "so there is no coherent dependency axis. Treat the composite score as ",
      "uninterpretable and analyse genes individually instead."),
      set_name, mean_r), call. = FALSE)
  }

  if (orient == "vulnerability") score <- -score  # higher = more dependent

  list(
    score = data.table(ModelID = rownames(M),
                       score = as.numeric(score),
                       raw_mean_ge = rowMeans(M, na.rm = TRUE)),
    per_gene = per_gene,
    diagnostics = diagnostics,
    orientation = orient
  )
}

#' Build the DDR feature matrix (genotype + transcriptional).
#'
#' @param dm output of load_depmap()
#' @param amp_quantile copy-number amplification is called within-dataset at
#'        this upper quantile rather than at a fixed log2 threshold, because
#'        DepMap rescales OmicsCNGene between releases and a hard-coded cutoff
#'        silently changes meaning.
build_ddr_features <- function(dm, amp_quantile = 0.90, del_quantile = 0.10) {
  ids <- rownames(dm$gene_effect)
  ft  <- data.table(ModelID = ids)

  ## ---- HRR damaging-mutation burden (HRD PROXY, not a true HRD score) ----
  if (!is.null(dm$mut)) {
    hrr <- intersect(HRR_GENES, colnames(dm$mut))
    message(sprintf("[ddr] HRR mutation genes present: %d/%d", length(hrr), length(HRR_GENES)))
    Mh <- dm$mut[, hrr, drop = FALSE]
    Mh <- Mh[match(ids, rownames(Mh)), , drop = FALSE]
    ft[, hrr_n_mut   := rowSums(Mh > 0, na.rm = TRUE)]
    ft[, hrr_any_mut := as.integer(hrr_n_mut > 0)]

    for (g in intersect(GENOTYPE_SNVS, colnames(dm$mut))) {
      v <- dm$mut[match(ids, rownames(dm$mut)), g]
      ft[[paste0("mut_", g)]] <- as.integer(v > 0)
    }
  } else {
    message("[ddr] no mutation matrix; skipping genotype features")
  }

  ## ---- copy number: amplifications and deletions -------------------------
  if (!is.null(dm$cn)) {
    cn <- dm$cn[match(ids, rownames(dm$cn)), , drop = FALSE]
    for (g in intersect(GENOTYPE_AMPS, colnames(cn))) {
      v  <- cn[, g]
      th <- quantile(v, amp_quantile, na.rm = TRUE)
      ft[[paste0("amp_", g)]] <- as.integer(v >= th)
      message(sprintf("[ddr] amp_%s threshold = %.3f (q%.2f), %d lines called",
                      g, th, amp_quantile, sum(v >= th, na.rm = TRUE)))
    }
    for (g in intersect(GENOTYPE_DELS, colnames(cn))) {
      v  <- cn[, g]
      th <- quantile(v, del_quantile, na.rm = TRUE)
      ft[[paste0("del_", g)]] <- as.integer(v <= th)
    }
    ft[, cn_gene_median := apply(cn, 1, median, na.rm = TRUE)]
  } else {
    message("[ddr] no copy-number matrix; skipping CN features")
  }

  ## ---- replication-stress transcriptional score --------------------------
  if (!is.null(dm$expr)) {
    rs <- match_genes(REPSTRESS_SIG, colnames(dm$expr), "replication-stress signature")
    E  <- dm$expr[match(ids, rownames(dm$expr)), rs, drop = FALSE]
    ft[, repstress_score := rowMeans(scale(E), na.rm = TRUE)]
    message(sprintf("[ddr] replication-stress score from %d genes", length(rs)))
  } else {
    message("[ddr] no expression matrix; skipping replication-stress score")
  }

  ft
}
