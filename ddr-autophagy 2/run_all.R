#!/usr/bin/env Rscript
## run_all.R -----------------------------------------------------------------
## DDR x autophagy co-dependency analysis -- full pipeline.
##
## USAGE
##   Rscript run_all.R --data_dir=/path/to/depmap --out_dir=results
##   Rscript run_all.R --selftest                  # synthetic data, no download
##
## OPTIONS
##   --data_dir=PATH   directory holding the DepMap CSVs       (required)
##   --out_dir=PATH    where results go                        (default: results)
##   --no_expression   skip the big expression file
##   --null_B=N        random gene sets for the empirical null (default 200, 0 = skip)
##   --min_lineage_n=N minimum cell lines per lineage          (default 10)
##   --selftest        run on synthetic data with a planted effect
##
## REQUIRED DEPMAP FILES  (https://depmap.org/portal/data_page/?tab=allData)
##   CRISPRGeneEffect.csv
##   Model.csv
##   OmicsSomaticMutationsMatrixDamaging.csv
##   OmicsCNGene.csv
##   OmicsExpressionProteinCodingGenesTPMLogp1.csv   (optional)
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

## ---- argument parsing -----------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
getopt <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[1])
}
hasflag <- function(name) any(args == paste0("--", name))

SELFTEST      <- hasflag("selftest")
DATA_DIR      <- getopt("data_dir")
OUT_DIR       <- getopt("out_dir", "results")
LOAD_EXPR     <- !hasflag("no_expression")
NULL_B        <- as.integer(getopt("null_B", "200"))
MIN_LINEAGE_N <- as.integer(getopt("min_lineage_n", "10"))

if (!SELFTEST && is.null(DATA_DIR)) {
  cat("ERROR: --data_dir is required (or use --selftest).\n\n")
  cat("  Rscript run_all.R --data_dir=/path/to/depmap\n")
  cat("  Rscript run_all.R --selftest\n\n")
  quit(status = 1)
}

## ---- source modules -------------------------------------------------------
here <- tryCatch({
  sf <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(sf)) dirname(normalizePath(sub("^--file=", "", sf[1]))) else getwd()
}, error = function(e) getwd())

for (f in c("01_genesets.R", "02_load_depmap.R", "03_scores.R",
            "04_stats.R", "05_lineage.R")) {
  path <- file.path(here, "R", f)
  if (!file.exists(path)) path <- file.path(here, f)
  if (!file.exists(path)) stop("Cannot find module: ", f)
  source(path)
}

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("\n")
cat("===========================================================\n")
cat(" DDR x AUTOPHAGY CO-DEPENDENCY ANALYSIS\n")
cat(sprintf(" mode: %s\n", if (SELFTEST) "SELF-TEST (synthetic)" else "real DepMap data"))
cat(sprintf(" out : %s\n", normalizePath(OUT_DIR, mustWork = FALSE)))
cat("===========================================================\n\n")

## ---------------------------------------------------------------------------
## 1. LOAD
## ---------------------------------------------------------------------------
if (SELFTEST) {
  src <- file.path(here, "tests", "make_fixtures.R")
  if (!file.exists(src)) src <- file.path(here, "make_fixtures.R")
  source(src)
  dm <- make_synthetic_depmap(n_lines = 900, planted_rho = -0.35, seed = 42)
  cat(sprintf("\n[selftest] planted autophagy<->ATR coupling: rho = %.2f\n",
              attr(dm, "planted_rho")))
} else {
  dm <- load_depmap(DATA_DIR, load_expression = LOAD_EXPR)
}

## ---- covariates -----------------------------------------------------------
covars <- merge(
  dm$line_qc[, .(ModelID, mean_gene_effect)],
  dm$model[, .(ModelID, lineage_raw = OncotreeLineage)],
  by = "ModelID"
)
covars[, lineage := collapse_rare_lineages(lineage_raw, min_n = MIN_LINEAGE_N)]
cat(sprintf("\n[main] %d models with lineage + QC covariates\n", nrow(covars)))

## ---------------------------------------------------------------------------
## 2. SCORES
## ---------------------------------------------------------------------------
cat("\n--- autophagy dependency score (PRIMARY: conjugation cascade) ---\n")
auto_primary <- build_dependency_score(dm$gene_effect, AUTOPHAGY_PRIMARY,
                                       "autophagy_conjugation")

cat("\n--- autophagy dependency score (SENSITIVITY: all modules) ---\n")
auto_all <- build_dependency_score(dm$gene_effect, AUTOPHAGY_ALL,
                                   "autophagy_all_modules")

cat("\n--- lysosome dependency score (what CQ/HCQ actually hit) ---\n")
lyso <- tryCatch(build_dependency_score(dm$gene_effect, LYSOSOME_SET, "lysosome"),
                 error = function(e) { message("[main] lysosome score failed: ",
                                               conditionMessage(e)); NULL })

cat("\n--- NEGATIVE CONTROL score (should show nothing) ---\n")
negctl <- tryCatch(build_dependency_score(dm$gene_effect, NEGATIVE_CONTROL_SET,
                                          "negative_control"),
                   error = function(e) { message("[main] neg control failed: ",
                                                 conditionMessage(e)); NULL })

fwrite(auto_primary$per_gene, file.path(OUT_DIR, "autophagy_per_gene_stats.csv"))
fwrite(auto_primary$score,    file.path(OUT_DIR, "autophagy_score_per_line.csv"))

## ---------------------------------------------------------------------------
## 3. DDR FEATURES
## ---------------------------------------------------------------------------
cat("\n--- DDR feature matrix ---\n")
feat <- build_ddr_features(dm)
fwrite(feat, file.path(OUT_DIR, "ddr_features_per_line.csv"))

## ---------------------------------------------------------------------------
## 4. TEST A -- co-dependency
## ---------------------------------------------------------------------------
cat("\n--- TEST A: autophagy vs DDR-gene co-dependency ---\n")
cat("    (negative rho = hypothesis-consistent; see 04_stats.R)\n\n")
A_primary <- test_codependency(auto_primary$score, dm$gene_effect, covars)
print(head(A_primary[!is.na(rho)], 20), digits = 3)
fwrite(A_primary, file.path(OUT_DIR, "testA_codependency_primary.csv"))

A_all <- test_codependency(auto_all$score, dm$gene_effect, covars)
fwrite(A_all, file.path(OUT_DIR, "testA_codependency_allmodules.csv"))

if (!is.null(negctl)) {
  A_neg <- test_codependency(negctl$score, dm$gene_effect, covars)
  fwrite(A_neg, file.path(OUT_DIR, "testA_codependency_negcontrol.csv"))
  cat(sprintf("\n[control] negative-control set: %d/%d targets at FDR<0.05 (expect ~0)\n",
              sum(A_neg$fdr < 0.05, na.rm = TRUE), nrow(A_neg)))
  cat(sprintf("[control] primary autophagy set: %d/%d targets at FDR<0.05\n",
              sum(A_primary$fdr < 0.05, na.rm = TRUE), nrow(A_primary)))
}

if (!is.null(lyso)) {
  A_lyso <- test_codependency(lyso$score, dm$gene_effect, covars)
  fwrite(A_lyso, file.path(OUT_DIR, "testA_codependency_lysosome.csv"))
  cat(sprintf("[contrast] lysosome set: %d/%d targets at FDR<0.05\n",
              sum(A_lyso$fdr < 0.05, na.rm = TRUE), nrow(A_lyso)))
  cat("[contrast] If lysosome >> autophagy here, CQ/HCQ effects are likely\n")
  cat("           lysosomal rather than autophagic. Report this either way.\n")
}

## ---------------------------------------------------------------------------
## 5. TEST B -- DDR genotype features
## ---------------------------------------------------------------------------
cat("\n--- TEST B: DDR features vs autophagy dependency ---\n")
cat("    (positive beta = hypothesis-consistent)\n\n")
B_primary <- test_ddr_features(auto_primary$score, feat, covars)
if (nrow(B_primary)) print(head(B_primary, 20), digits = 3)
fwrite(B_primary, file.path(OUT_DIR, "testB_ddr_features.csv"))

## ---------------------------------------------------------------------------
## 6. EMPIRICAL NULL
## ---------------------------------------------------------------------------
if (NULL_B > 0) {
  cat(sprintf("\n--- Empirical null (%d random gene sets) ---\n", NULL_B))
  cat("    Parametric p-values from composite scores are optimistic.\n")
  cat("    Trust p_empirical over p.\n\n")
  nullres <- empirical_null_codependency(
    dm$gene_effect, covars, A_primary,
    set_size = auto_primary$diagnostics$n_genes_used,
    B = NULL_B, real_sd = auto_primary$per_gene$sd_ge)
  print(head(nullres, 15), digits = 3)
  fwrite(nullres, file.path(OUT_DIR, "testA_empirical_null.csv"))

  merged <- merge(A_primary, nullres[, .(gene, null_mean, z_vs_null,
                                         p_empirical, p_empirical_floor,
                                         fdr_empirical, p_z, fdr_z)],
                  by = "gene", all.x = TRUE)
  setorder(merged, p_z, na.last = TRUE)
  fwrite(merged, file.path(OUT_DIR, "testA_with_empirical_null.csv"))
  cat(sprintf("\n[null] surviving empirical FDR<0.05 : %d\n",
              sum(merged$fdr_empirical < 0.05, na.rm = TRUE)))
  cat(sprintf("[null] surviving z-based  FDR<0.05 : %d\n",
              sum(merged$fdr_z < 0.05, na.rm = TRUE)))
  cat(sprintf("[null] empirical p floor at B=%d   : %.4f\n", NULL_B, 1/(NULL_B+1)))
  if (1/(NULL_B+1) > 0.05/length(intersect(DDR_TARGETS, colnames(dm$gene_effect)))) {
    cat("[null] NOTE: B is too small for ANY target to reach empirical\n")
    cat("[null]       FDR<0.05. Read z_vs_null / fdr_z, or raise --null_B.\n")
  }
} else {
  cat("\n[main] empirical null skipped (--null_B=0)\n")
}

## ---------------------------------------------------------------------------
## 7. LINEAGE RANKING
## ---------------------------------------------------------------------------
cat("\n--- Lineage priority ranking ---\n")
lin <- rank_lineages(auto_primary$score, feat, covars, dm$gene_effect,
                     min_n = MIN_LINEAGE_N)
report_lineages(lin, top_n = 15, file = file.path(OUT_DIR, "lineage_ranking.csv"))

## ---------------------------------------------------------------------------
## 8. SESSION PROVENANCE
## ---------------------------------------------------------------------------
sink(file.path(OUT_DIR, "session_info.txt"))
cat("DDR x autophagy analysis\n")
cat("run at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("mode  : ", if (SELFTEST) "selftest" else "real", "\n", sep = "")
cat("data  : ", if (SELFTEST) "synthetic" else normalizePath(DATA_DIR), "\n\n", sep = "")
cat("primary autophagy genes used: ",
    paste(auto_primary$diagnostics$genes_used, collapse = ", "), "\n\n", sep = "")
cat("coherence of primary score:\n")
cat("  mean pairwise r  : ", auto_primary$diagnostics$mean_pairwise_r, "\n", sep = "")
cat("  Cronbach alpha   : ", auto_primary$diagnostics$cronbach_alpha, "\n", sep = "")
cat("  PC1 var explained: ", auto_primary$diagnostics$pc1_var_explained, "\n\n", sep = "")
print(sessionInfo())
sink()

cat("\n===========================================================\n")
cat(" DONE. Files written to: ", OUT_DIR, "\n", sep = "")
cat("\n READ THESE FIRST, IN THIS ORDER:\n")
cat("  1. session_info.txt            -- is the score coherent at all?\n")
cat("  2. autophagy_per_gene_stats.csv-- do ATG genes show ANY dependency?\n")
cat("  3. testA_codependency_negcontrol.csv -- is the null behaving?\n")
cat("  4. testA_with_empirical_null.csv     -- the actual result\n")
cat("  5. lineage_ranking.csv               -- the deliverable\n")
cat("===========================================================\n\n")
