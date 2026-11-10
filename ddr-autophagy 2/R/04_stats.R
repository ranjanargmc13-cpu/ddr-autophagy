## 04_stats.R ----------------------------------------------------------------
## The statistical core.
##
## Two tests, deliberately kept separate because they answer different questions:
##
##   TEST A (co-dependency): across cell lines, does reliance on autophagy track
##     reliance on a given DDR gene? This is the functional-genetic version of
##     the hypothesis and needs no HRD score. It is the stronger design.
##
##   TEST B (genotype association): do lines carrying a DDR lesion depend more
##     on autophagy? This is what a trial would stratify on, but it leans on a
##     genotype proxy for a phenotype, so it is the weaker design.
##
## Both control for lineage and screen quality. Both are calibrated against an
## empirical null built from random gene sets of matched size, because the
## parametric p-value from a composite score is optimistic -- the score's own
## genes are correlated and the effective degrees of freedom are lower than n.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

#' Collapse lineages with too few lines into "Other".
#' A lineage represented by 1-2 lines cannot support a lineage fixed effect and
#' will otherwise consume a degree of freedom to perfectly fit its own points.
collapse_rare_lineages <- function(lineage, min_n = 10) {
  tab  <- table(lineage)
  rare <- names(tab)[tab < min_n]
  out  <- as.character(lineage)
  out[out %in% rare] <- "Other"
  if (length(rare))
    message(sprintf("[stats] collapsed %d lineage(s) with n < %d into 'Other': %s",
                    length(rare), min_n, paste(rare, collapse = ", ")))
  factor(out)
}

#' Residualize y on a covariate data.frame. Returns NA where covariates are NA.
residualize <- function(y, covars) {
  df <- cbind(data.frame(.y = as.numeric(y)), covars)
  ok <- stats::complete.cases(df)
  out <- rep(NA_real_, length(y))
  if (sum(ok) < (ncol(df) + 5L)) return(out)   # not enough data to fit
  fit <- stats::lm(.y ~ ., data = df[ok, , drop = FALSE])
  out[ok] <- stats::residuals(fit)
  out
}

## ---------------------------------------------------------------------------
## TEST A: autophagy dependency vs DDR-gene dependency (co-dependency)
## ---------------------------------------------------------------------------

#' @param autophagy_score data.table(ModelID, score)  [vulnerability-oriented]
#' @param gene_effect     models x genes Chronos matrix
#' @param covars          data.table(ModelID, lineage, mean_gene_effect)
#' @param targets         DDR gene symbols to test
#' @return data.table ranked by FDR
test_codependency <- function(autophagy_score, gene_effect, covars,
                              targets = DDR_TARGETS, method = "spearman") {
  tg <- intersect(targets, colnames(gene_effect))
  message(sprintf("[testA] testing %d/%d DDR target genes", length(tg), length(targets)))

  d <- merge(autophagy_score[, .(ModelID, auto = score)],
             covars, by = "ModelID")
  d <- d[ModelID %in% rownames(gene_effect)]
  setorder(d, ModelID)
  GE <- gene_effect[d$ModelID, tg, drop = FALSE]

  cv <- data.frame(lineage = d$lineage, mean_ge = d$mean_gene_effect)
  auto_r <- residualize(d$auto, cv)

  res <- rbindlist(lapply(tg, function(g) {
    ## Note: the DDR gene's own effect is residualized on the SAME covariates,
    ## so the reported correlation is a partial correlation given lineage and
    ## screen quality.
    gr <- residualize(GE[, g], cv)
    ok <- !is.na(auto_r) & !is.na(gr)
    n  <- sum(ok)
    if (n < 30) return(data.table(gene = g, n = n, rho = NA_real_, p = NA_real_,
                                  ddr_mean_ge = NA_real_, ddr_sd_ge = NA_real_))
    ct <- suppressWarnings(stats::cor.test(auto_r[ok], gr[ok], method = method))
    data.table(gene = g, n = n,
               rho = unname(ct$estimate), p = ct$p.value,
               ddr_mean_ge = mean(GE[, g], na.rm = TRUE),
               ddr_sd_ge   = sd(GE[, g], na.rm = TRUE))
  }))

  ## IMPORTANT sign note for interpretation:
  ## auto is vulnerability-oriented (higher = more autophagy-dependent);
  ## gr is raw Chronos (more negative = more DDR-gene-dependent).
  ## So a NEGATIVE rho means: more autophagy-dependent lines are ALSO more
  ## DDR-gene-dependent. That is the direction the hypothesis predicts.
  res[, hypothesis_consistent := !is.na(rho) & rho < 0]
  res[, fdr := stats::p.adjust(p, method = "BH")]
  setorder(res, p, na.last = TRUE)
  res[]
}

## ---------------------------------------------------------------------------
## TEST B: DDR genotype/phenotype feature vs autophagy dependency
## ---------------------------------------------------------------------------

#' @param features data.table(ModelID, <feature columns>)
#' @param feature_cols which columns to test; default = all numeric non-ID
test_ddr_features <- function(autophagy_score, features, covars,
                              feature_cols = NULL, min_group = 15) {
  d <- merge(autophagy_score[, .(ModelID, auto = score)], covars, by = "ModelID")
  d <- merge(d, features, by = "ModelID")

  if (is.null(feature_cols)) {
    feature_cols <- setdiff(names(features), "ModelID")
  }
  feature_cols <- intersect(feature_cols, names(d))
  message(sprintf("[testB] testing %d DDR features", length(feature_cols)))

  res <- rbindlist(lapply(feature_cols, function(f) {
    x <- d[[f]]
    if (all(is.na(x))) return(NULL)
    is_bin <- all(x %in% c(0, 1, NA))

    ## Binary features need both groups adequately populated, otherwise the
    ## coefficient is driven by a handful of lines.
    if (is_bin) {
      n1 <- sum(x == 1, na.rm = TRUE); n0 <- sum(x == 0, na.rm = TRUE)
      if (min(n1, n0) < min_group)
        return(data.table(feature = f, type = "binary", n = sum(!is.na(x)),
                          n_pos = n1, beta = NA_real_, se = NA_real_,
                          t = NA_real_, p = NA_real_,
                          note = sprintf("skipped: smaller group n=%d < %d",
                                         min(n1, n0), min_group)))
    }

    df <- data.frame(auto = d$auto, x = x,
                     lineage = d$lineage, mean_ge = d$mean_gene_effect)
    ok <- stats::complete.cases(df)
    if (sum(ok) < 50) return(NULL)
    fit <- try(stats::lm(auto ~ x + mean_ge + lineage, data = df[ok, ]), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    cf <- summary(fit)$coefficients
    if (!"x" %in% rownames(cf)) return(NULL)

    data.table(feature = f,
               type  = if (is_bin) "binary" else "continuous",
               n     = sum(ok),
               n_pos = if (is_bin) sum(df$x[ok] == 1) else NA_integer_,
               beta  = cf["x", 1], se = cf["x", 2],
               t     = cf["x", 3], p = cf["x", 4],
               note  = NA_character_)
  }), fill = TRUE)

  if (!nrow(res)) return(res)
  ## beta > 0 means the DDR lesion raises autophagy vulnerability -- the
  ## hypothesis-consistent direction.
  res[, hypothesis_consistent := !is.na(beta) & beta > 0]
  res[, fdr := stats::p.adjust(p, method = "BH")]
  setorder(res, p, na.last = TRUE)
  res[]
}

## ---------------------------------------------------------------------------
## Empirical null from random gene sets
## ---------------------------------------------------------------------------

#' Rebuild the score from B random gene sets of matched size and re-run TEST A.
#' Gives an empirical p-value per DDR target that accounts for the fact that
#' ANY correlated composite score will show some structure.
#'
#' @param eligible_genes pool to sample from. Default: genes with comparable
#'        across-line variance to the real set, so the null is matched on
#'        variance rather than being trivially flat.
empirical_null_codependency <- function(gene_effect, covars, real_result,
                                        set_size, targets = DDR_TARGETS,
                                        B = 200, seed = 1,
                                        eligible_genes = NULL,
                                        match_variance = TRUE,
                                        real_sd = NULL) {
  set.seed(seed)
  tg <- intersect(targets, colnames(gene_effect))

  if (is.null(eligible_genes)) {
    sds <- apply(gene_effect, 2, sd, na.rm = TRUE)
    pool <- names(sds)[is.finite(sds) & sds > 0.05]
    if (match_variance && !is.null(real_sd)) {
      lo <- quantile(real_sd, 0.10, na.rm = TRUE)
      hi <- quantile(real_sd, 0.90, na.rm = TRUE)
      pool2 <- names(sds)[is.finite(sds) & sds >= lo & sds <= hi]
      if (length(pool2) > 10 * set_size) pool <- pool2
      message(sprintf("[null] variance-matched pool: %d genes (SD %.3f-%.3f)",
                      length(pool), lo, hi))
    }
    eligible_genes <- setdiff(pool, c(AUTOPHAGY_ALL, LYSOSOME_SET, tg))
  }
  message(sprintf("[null] drawing %d random sets of %d genes from a pool of %d",
                  B, set_size, length(eligible_genes)))

  nulls <- matrix(NA_real_, nrow = B, ncol = length(tg),
                  dimnames = list(NULL, tg))

  for (b in seq_len(B)) {
    gs <- sample(eligible_genes, set_size)
    sc <- suppressWarnings(suppressMessages(
      build_dependency_score(gene_effect, gs, sprintf("null_%d", b),
                             min_sd = 0, orient = "vulnerability")$score))
    rb <- suppressWarnings(suppressMessages(
      test_codependency(sc, gene_effect, covars, targets = tg)))
    nulls[b, rb$gene] <- rb$rho
    if (b %% 25 == 0) message(sprintf("[null]   %d/%d", b, B))
  }

  out <- rbindlist(lapply(tg, function(g) {
    nv <- nulls[, g]; nv <- nv[is.finite(nv)]
    obs <- real_result[gene == g, rho]
    if (!length(obs) || is.na(obs) || !length(nv))
      return(data.table(gene = g, rho_obs = NA_real_, null_mean = NA_real_,
                        null_sd = NA_real_, z_vs_null = NA_real_,
                        p_empirical = NA_real_))
    ## Two-sided empirical p with the standard +1 correction.
    ## CEILING ON RESOLUTION: the smallest value this can return is
    ## 1/(B+1). With B=50 that is 0.0196, which after BH across ~47 targets
    ## cannot reach FDR<0.05 no matter how strong the effect is. Choose B
    ## with the multiple-testing burden in mind: to clear BH at rank k out of
    ## m targets you need 1/(B+1) <= (k/m)*alpha.
    p_emp <- (sum(abs(nv) >= abs(obs)) + 1) / (length(nv) + 1)
    zz    <- (obs - mean(nv)) / sd(nv)
    ## Normal-approximation p from the null's own mean/SD. Not floored by B,
    ## because it uses the null's location and spread rather than counting
    ## exceedances. The null rho distribution is near-normal here (it is a
    ## correlation of residualized composites), so this is well behaved -- but
    ## it is an approximation, so report it ALONGSIDE p_empirical, never
    ## instead of it. Where they disagree, trust p_empirical and raise B.
    p_z <- 2 * stats::pnorm(-abs(zz))
    data.table(gene = g, rho_obs = obs,
               null_mean = mean(nv), null_sd = sd(nv),
               z_vs_null = zz,
               p_empirical = p_emp,
               p_empirical_floor = 1 / (length(nv) + 1),
               p_z = p_z)
  }))
  out[, fdr_empirical := stats::p.adjust(p_empirical, method = "BH")]
  out[, fdr_z         := stats::p.adjust(p_z,         method = "BH")]
  setorder(out, p_z, na.last = TRUE)

  ## Warn when B is too small for the test to be able to succeed.
  floor_p <- 1 / (B + 1)
  need_p  <- (1 / length(tg)) * 0.05
  if (floor_p > need_p) {
    warning(sprintf(paste0(
      "[null] B=%d gives a minimum empirical p of %.4f, but clearing BH at ",
      "FDR<0.05 across %d targets needs p <= %.4f for the top hit. No target ",
      "CAN reach empirical significance at this B. Raise B to >= %d, or read ",
      "fdr_z / z_vs_null instead."),
      B, floor_p, length(tg), need_p, ceiling(1 / need_p)), call. = FALSE)
  }
  out[]
}
