## selftest_interactive.R ----------------------------------------------------
## Run the self-test from inside an interactive R / RStudio session.
##
## Use this INSTEAD of source("run_all.R"): run_all.R is a command-line script
## and calls quit() when it gets no --data_dir, which would end your session.
##
##   setwd("path/to/ddr-autophagy")
##   source("selftest_interactive.R")
##
## Everything stays in your workspace afterwards, so you can poke at `A`,
## `auto`, `lin` etc. directly.
## ---------------------------------------------------------------------------

library(data.table)

## ---- locate and source the modules ----------------------------------------
## Works whether the files are in R/ + tests/ or all flat in one folder.
.src <- function(f, sub = NULL) {
  cands <- c(if (!is.null(sub)) file.path(sub, f), f)
  for (p in cands) if (file.exists(p)) { source(p); return(invisible(p)) }
  stop("Cannot find ", f, " -- is your working directory right? getwd() = ", getwd())
}
.src("01_genesets.R",   "R")
.src("02_load_depmap.R","R")
.src("03_scores.R",     "R")
.src("04_stats.R",      "R")
.src("05_lineage.R",    "R")
.src("make_fixtures.R", "tests")

cat("\n==================== SELF-TEST ====================\n")

## ---- 1. synthetic data with a known answer --------------------------------
dm <- make_synthetic_depmap(n_lines = 900, planted_rho = -0.35, seed = 42)
planted <- attr(dm, "planted_targets")

covars <- merge(dm$line_qc[, .(ModelID, mean_gene_effect)],
                dm$model[, .(ModelID, lineage_raw = OncotreeLineage)],
                by = "ModelID")
covars[, lineage := collapse_rare_lineages(lineage_raw, min_n = 10)]

## ---- 2. scores -------------------------------------------------------------
cat("\n--- autophagy score (primary) ---\n")
auto <- build_dependency_score(dm$gene_effect, AUTOPHAGY_PRIMARY,
                               "autophagy_conjugation")

cat("\n--- negative control score ---\n")
negc <- build_dependency_score(dm$gene_effect, NEGATIVE_CONTROL_SET,
                               "negative_control")

## ---- 3. TEST A -------------------------------------------------------------
cat("\n--- TEST A: co-dependency (negative rho = hypothesis-consistent) ---\n")
A <- test_codependency(auto$score, dm$gene_effect, covars)
print(head(A[!is.na(rho)], 15), digits = 3)

Aneg <- test_codependency(negc$score, dm$gene_effect, covars)

## ---- 4. DDR features + TEST B ---------------------------------------------
feat <- build_ddr_features(dm)
B <- test_ddr_features(auto$score, feat, covars)
if (nrow(B)) {
  cat("\n--- TEST B: DDR features (positive beta = hypothesis-consistent) ---\n")
  print(head(B, 10), digits = 3)
}

## ---- 5. lineage ranking ----------------------------------------------------
lin <- rank_lineages(auto$score, feat, covars, dm$gene_effect, min_n = 10)
report_lineages(lin, top_n = 10)

## ---- 6. VERDICT ------------------------------------------------------------
## Three criteria. All must pass, or the code is wrong -- not the biology.
top10     <- head(A[!is.na(rho)]$gene, 10)
recovered <- intersect(planted, top10)
right_dir <- A[gene %in% planted & !is.na(rho), all(rho < 0)]
n_neg_sig <- sum(Aneg$fdr < 0.05, na.rm = TRUE)
n_pos_sig <- sum(A$fdr    < 0.05, na.rm = TRUE)

c1 <- length(recovered) >= 3
c2 <- isTRUE(right_dir)
c3 <- n_neg_sig <= 2

cat("\n==================== VERDICT ====================\n")
cat(sprintf("planted targets        : %s\n", paste(planted, collapse = ", ")))
cat(sprintf("recovered in top 10    : %s\n",
            if (length(recovered)) paste(recovered, collapse = ", ") else "(none)"))
cat(sprintf("\n[%s] 1. recovers >=3 planted targets in top 10   (got %d)\n",
            ifelse(c1, "PASS", "FAIL"), length(recovered)))
cat(sprintf("[%s] 2. planted targets all have negative rho    \n",
            ifelse(c2, "PASS", "FAIL")))
cat(sprintf("[%s] 3. negative control finds ~nothing          (got %d at FDR<0.05)\n",
            ifelse(c3, "PASS", "FAIL"), n_neg_sig))
cat(sprintf("\n    (real autophagy set found %d at FDR<0.05, for contrast)\n", n_pos_sig))

if (c1 && c2 && c3) {
  cat("\n>>> ALL PASS. Pipeline recovers a known effect and controls the null.\n")
  cat(">>> Safe to run on real DepMap data:\n")
  cat(">>>   Rscript run_all.R --data_dir=/path/to/depmap\n")
} else {
  cat("\n>>> FAILURE. Do NOT run on real data yet.\n")
  cat(">>> Inspect: A, Aneg, auto$diagnostics, covars\n")
  cat(">>> Most likely causes: data.table version behaviour in rank_lineages(),\n")
  cat(">>> or a sign flip in build_dependency_score().\n")
}
cat("=================================================\n\n")

cat("Objects left in your workspace: dm, covars, auto, negc, A, Aneg, feat, B, lin\n\n")
