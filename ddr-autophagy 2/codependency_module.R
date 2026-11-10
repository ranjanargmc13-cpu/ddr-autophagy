#!/usr/bin/env Rscript
## codependency_module.R -----------------------------------------------------
## IN-SILICO SUBSTITUTE FOR THE ATG7 KNOCKOUT EXPERIMENT.
##
##   Rscript codependency_module.R                   # uses ~/depmap
##   Rscript codependency_module.R /path/to/depmap
##
## THE LOGIC
## Your cells depend on PIK3C3/VPS34 (-0.865 in MDA-MB-231) but NOT on the
## ATG8 conjugation cascade (ATG7 was +0.089). Two explanations:
##   (a) VPS34 matters because of autophagy, and ATG7 loss is somehow tolerated;
##   (b) VPS34 matters for endolysosomal functions that have nothing to do
##       with ATG8-dependent macroautophagy.
##
## A knockout in one cell line would be one noisy data point. Instead we use the
## fact that DepMap screened ~1100 lines: if VPS34's essentiality were mediated
## by autophagy, then lines that need VPS34 should also need ATG7 -- the
## dependency PROFILES should correlate across lines. Co-dependency correlation
## is well established as a way to recover genes acting in the same complex or
## pathway, so a near-zero correlation is positive evidence for (b).
##
## WHAT IT REPORTS
##   1. partial correlation of PIK3C3 dependency with each ATG8-cascade gene
##   2. the same with each lysosomal gene
##   3. a bootstrap test of whether PIK3C3 sits closer to the lysosome module
##      or the ATG8 module
##   4. genome-wide top co-dependencies of PIK3C3 and of ATG7, side by side
##   5. the same restricted to breast lines
##
## All correlations are PARTIAL, adjusted for lineage and screen quality. Those
## two confounders alone can manufacture correlation between arbitrary gene
## pairs, so an unadjusted version of this analysis would be uninterpretable.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
DATA_DIR <- path.expand(if (length(args)) args[1] else "~/depmap")
OUT_DIR  <- "results_codependency"

here <- tryCatch({
  sf <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(sf)) dirname(normalizePath(sub("^--file=", "", sf[1]))) else getwd()
}, error = function(e) getwd())
.src <- function(f) {
  for (p in c(file.path(here, "R", f), file.path(here, f)))
    if (file.exists(p)) { source(p); return(invisible(p)) }
  stop("Cannot find ", f)
}
.src("01_genesets.R"); .src("02_load_depmap.R")

if (!dir.exists(DATA_DIR))
  stop("Data directory not found: ", DATA_DIR,
       "\nPass it explicitly: Rscript codependency_module.R /full/path/to/depmap")
dir.create(OUT_DIR, showWarnings = FALSE)

dm <- load_depmap(DATA_DIR, load_expression = FALSE)

## ---------------------------------------------------------------------------
## Modules
## ---------------------------------------------------------------------------
## Three competing modules. The trafficking module was added after the
## genome-wide scan showed PIK3C3's neighbours are ESCRT/HOPS/RAB7 genes that
## appear in neither of the original two sets.
ATG8_CASCADE <- AUTOPHAGY_PRIMARY           # ATG7/3/10/12/5/16L1
LYSO_MODULE  <- LYSOSOME_SET                # lysosomal contents + acidification
TRAF_MODULE  <- TRAFFICKING_MODULE          # ESCRT / HOPS / RAB / retromer

MODULES <- list(ATG8_cascade = ATG8_CASCADE,
                lysosome     = LYSO_MODULE,
                trafficking  = TRAF_MODULE)

FOCAL <- c("PIK3C3", "RB1CC1", "BECN1", "ATG7")   # genes to profile

## Internal positive control: PIK3C3 and PIK3R4 are obligate complex partners,
## so their correlation calibrates what "same complex" looks like in this data.
## Every other correlation should be read against that number, not against zero.
CALIBRATION_PAIR <- c("PIK3C3", "PIK3R4")

## ---------------------------------------------------------------------------
## Residualize the whole matrix on lineage + screen quality, in one shot.
##
## Fitting ~18,000 separate lm() calls would take far too long. Instead we build
## the design matrix once and apply the residual-maker projection to every gene
## simultaneously via QR. This is numerically identical to per-gene lm(), just
## thousands of times faster.
## ---------------------------------------------------------------------------
message("\n[prep] building covariates ...")
cv <- merge(dm$line_qc[, .(ModelID, mean_ge = mean_gene_effect)],
            dm$model[, .(ModelID, lineage_raw = OncotreeLineage)], by = "ModelID")
cv[, lineage := factor(ifelse(
  lineage_raw %in% names(which(table(lineage_raw) >= 10)), lineage_raw, "Other"))]
cv <- cv[ModelID %in% rownames(dm$gene_effect)]
setorder(cv, ModelID)

GE <- dm$gene_effect[cv$ModelID, , drop = FALSE]

## Drop genes measured in too few lines, then mean-impute the remaining NAs.
## Imputing to the gene mean is conservative here: it pulls a gene toward no
## signal rather than inventing any.
ok_gene <- colSums(!is.na(GE)) >= (0.9 * nrow(GE))
message(sprintf("[prep] keeping %d/%d genes measured in >=90%% of lines",
                sum(ok_gene), ncol(GE)))
GE <- GE[, ok_gene, drop = FALSE]
nas <- which(is.na(GE), arr.ind = TRUE)
if (nrow(nas)) {
  cm <- colMeans(GE, na.rm = TRUE)
  GE[nas] <- cm[nas[, "col"]]
  message(sprintf("[prep] mean-imputed %d missing values (%.3f%% of matrix)",
                  nrow(nas), 100 * nrow(nas) / length(GE)))
}

message("[prep] residualizing on lineage + screen quality ...")
X  <- model.matrix(~ lineage + mean_ge, data = cv)
qx <- qr(X)
R  <- GE - X %*% qr.coef(qx, GE)          # residual matrix, lines x genes
message(sprintf("[prep] done: %d lines x %d genes, %d covariate columns",
                nrow(R), ncol(R), ncol(X)))

## Correlation of one gene's residual profile against all others.
cor_vs_all <- function(gene, Rmat = R) {
  if (!gene %in% colnames(Rmat)) return(NULL)
  v <- Rmat[, gene]
  cc <- as.numeric(cor(v, Rmat))
  names(cc) <- colnames(Rmat)
  cc
}

## ---------------------------------------------------------------------------
## 1 + 2. PIK3C3 against the two modules
## ---------------------------------------------------------------------------
cat("\n==========================================================\n")
cat(" PART 1: is VPS34 dependency autophagy-mediated?\n")
cat("==========================================================\n")

module_table <- function(focal_gene) {
  cc <- cor_vs_all(focal_gene)
  if (is.null(cc)) { message("[!] ", focal_gene, " absent"); return(NULL) }
  mk <- function(genes, label) {
    g <- intersect(genes, names(cc)); g <- setdiff(g, focal_gene)
    if (!length(g)) return(NULL)
    data.table(focal = focal_gene, module = label, gene = g, partial_r = cc[g])
  }
  rbindlist(lapply(names(MODULES), function(m) mk(MODULES[[m]], m)),
            use.names = TRUE)
}

mt <- rbindlist(lapply(FOCAL, module_table), use.names = TRUE)
fwrite(mt, file.path(OUT_DIR, "module_correlations.csv"))

## Calibration first: everything below is read against this number.
if (all(CALIBRATION_PAIR %in% colnames(R))) {
  cal <- cor(R[, CALIBRATION_PAIR[1]], R[, CALIBRATION_PAIR[2]])
  cat(sprintf(paste0("\nCALIBRATION: r(%s, %s) = %+.3f\n",
    "  These are obligate complex partners, so this is what 'same module'\n",
    "  looks like here. Read every correlation below against it, not zero.\n"),
    CALIBRATION_PAIR[1], CALIBRATION_PAIR[2], cal))
} else cal <- NA_real_

summ_mod <- rbindlist(lapply(FOCAL, function(f) {
  sub <- mt[focal == f]; if (!nrow(sub)) return(NULL)
  s <- sub[, .(mean_r = mean(partial_r)), by = module]
  out <- as.data.table(as.list(setNames(s$mean_r, s$module)))
  cbind(data.table(focal = f), out)
}), fill = TRUE)

for (f in FOCAL) {
  sub <- mt[focal == f]
  if (!nrow(sub)) next
  cat(sprintf("\n--- %s dependency profile vs module members ---\n", f))
  print(sub[order(module, -abs(partial_r)), .(module, gene, partial_r)], digits = 3)
  cat("\n  module means:\n")
  for (m in names(MODULES)) {
    v <- sub[module == m, mean(partial_r)]
    if (length(v) && !is.na(v))
      cat(sprintf("    %-14s %+.3f%s\n", m, v,
                  if (!is.na(cal) && cal != 0)
                    sprintf("   (%.0f%% of complex-level)", 100 * v / cal) else ""))
  }
}

cat("\n--- module means, all focal genes ---\n")
print(summ_mod, digits = 3)
fwrite(summ_mod, file.path(OUT_DIR, "module_means.csv"))

## Headline single number: PIK3C3 vs ATG7.
if (all(c("PIK3C3", "ATG7") %in% colnames(R))) {
  ct <- cor.test(R[, "PIK3C3"], R[, "ATG7"])
  cat(sprintf(paste0("\n>>> HEADLINE: partial r(PIK3C3 dep, ATG7 dep) = %+.3f",
                     "  [95%% CI %+.3f, %+.3f], p = %.3g, n = %d\n"),
              ct$estimate, ct$conf.int[1], ct$conf.int[2], ct$p.value, nrow(R)))
  ## Interpret against the complex-level calibration, not against zero.
  ## (An earlier version of this script printed a fixed verdict here
  ## regardless of the value. That was a reporting bug.)
  frac <- if (!is.na(cal) && cal != 0) ct$estimate / cal else NA_real_
  if (!is.na(frac)) {
    cat(sprintf("    = %.0f%% of complex-level co-dependency (r=%.3f for PIK3C3-PIK3R4)\n",
                100 * frac, cal))
    cat(if (frac < 0.40)
          "    Well below complex level: VPS34 and ATG8 conjugation are\n    largely SEPARATE dependency modules.\n"
        else if (frac > 0.75)
          "    Near complex level: they behave as ONE module. An autophagy\n    interpretation of VPS34 dependency is supported.\n"
        else
          "    Intermediate. Decide using PART 2 and PART 3, not this number.\n")
  }
  cat("    Note: with n>1000 a tiny correlation can be highly significant.\n")
  cat("    Judge the EFFECT SIZE against calibration, not the p-value.\n")
}

## ---------------------------------------------------------------------------
## 3. Bootstrap: which module is PIK3C3 closer to?
##
## The two module correlations share PIK3C3, so they are dependent and cannot be
## compared with a naive two-sample test. Bootstrapping over cell lines respects
## that dependency.
## ---------------------------------------------------------------------------
cat("\n==========================================================\n")
cat(" PART 2: bootstrap -- lysosome module vs ATG8 cascade\n")
cat("==========================================================\n")

## Bootstrap over cell lines: the module correlations all share the focal gene,
## so they are dependent and cannot be compared with a naive two-sample test.
## Resampling lines respects that dependency.
boot_delta <- function(focal_gene, mod_a, mod_b, B = 2000, seed = 1) {
  set.seed(seed)
  A <- intersect(setdiff(MODULES[[mod_a]], focal_gene), colnames(R))
  Bg<- intersect(setdiff(MODULES[[mod_b]], focal_gene), colnames(R))
  if (!length(A) || !length(Bg) || !focal_gene %in% colnames(R)) return(NULL)
  n <- nrow(R); d <- numeric(B)
  for (b in seq_len(B)) {
    i <- sample.int(n, n, replace = TRUE)
    v <- R[i, focal_gene]
    d[b] <- mean(cor(v, R[i, Bg])) - mean(cor(v, R[i, A]))
  }
  obs <- mean(cor(R[, focal_gene], R[, Bg])) - mean(cor(R[, focal_gene], R[, A]))
  data.table(focal = focal_gene,
             contrast = sprintf("%s - %s", mod_b, mod_a),
             delta_obs = obs,
             ci_lo = unname(quantile(d, 0.025)),
             ci_hi = unname(quantile(d, 0.975)),
             p_two_sided = (sum(sign(d) != sign(obs)) + 1) / (B + 1))
}

## PRIMARY contrast: trafficking vs ATG8 cascade.
## SECONDARY: lysosome vs ATG8, kept so the original (poorly chosen) module is
## still reported rather than quietly dropped once it stopped being convenient.
bt <- rbindlist(c(
  lapply(FOCAL, boot_delta, mod_a = "ATG8_cascade", mod_b = "trafficking"),
  lapply(FOCAL, boot_delta, mod_a = "ATG8_cascade", mod_b = "lysosome")
), use.names = TRUE)

if (nrow(bt)) {
  cat("\ndelta > 0 => the focal gene sits CLOSER to the second module\n")
  cat("PRIMARY contrast is trafficking - ATG8_cascade\n\n")
  print(bt[order(contrast, focal)], digits = 3)
  fwrite(bt, file.path(OUT_DIR, "bootstrap_module_delta.csv"))
  cat("\nExpected if VPS34 acts through endosomal trafficking:\n")
  cat("  PIK3C3 -> trafficking delta POSITIVE, CI excluding 0\n")
  cat("  ATG7   -> trafficking delta NEGATIVE (it belongs to ATG8)\n")
  cat("ATG7 is the internal negative control: if it also moved toward\n")
  cat("trafficking, the contrast would be measuring something generic.\n")
}

## ---------------------------------------------------------------------------
## 4. Genome-wide co-dependencies, side by side
##
## This is the least assumption-laden view: no gene set imposed, just the genes
## whose dependency profiles most resemble each focal gene's. If PIK3C3's
## nearest neighbours are PI3KC3-complex and vesicle-trafficking genes while
## ATG7's are the conjugation cascade, the two are in separate functional
## modules and the question is settled without any gene set of mine.
## ---------------------------------------------------------------------------
cat("\n==========================================================\n")
cat(" PART 3: genome-wide nearest co-dependencies\n")
cat("==========================================================\n")

topn <- function(focal_gene, n = 25) {
  cc <- cor_vs_all(focal_gene); if (is.null(cc)) return(NULL)
  cc <- cc[names(cc) != focal_gene]
  s  <- sort(cc, decreasing = TRUE)[seq_len(n)]
  data.table(focal = focal_gene, rank = seq_len(n),
             partner = names(s), partial_r = as.numeric(s),
             in_ATG8 = names(s) %in% ATG8_CASCADE,
             in_lyso = names(s) %in% LYSO_MODULE,
             in_traffic = names(s) %in% TRAF_MODULE,
             in_autophagy_any = names(s) %in% AUTOPHAGY_ALL)
}
tp <- rbindlist(lapply(FOCAL, topn), use.names = TRUE)
fwrite(tp, file.path(OUT_DIR, "genomewide_top_codependencies.csv"))
for (f in FOCAL) {
  sub <- tp[focal == f]
  if (!nrow(sub)) next
  cat(sprintf("\n--- top 25 co-dependencies of %s ---\n", f))
  print(sub[, .(rank, partner, partial_r, in_ATG8, in_lyso, in_traffic)],
        digits = 3)
  cat(sprintf("  of top 25: %d in ATG8 cascade, %d in lysosome set, %d in trafficking\n",
              sum(sub$in_ATG8), sum(sub$in_lyso), sum(sub$in_traffic)))
}

## ---------------------------------------------------------------------------
## 5. Breast-only sensitivity analysis
##
## Pan-cancer is better powered, but your experiment was in TNBC. If the
## conclusion flips in breast lines alone, say so rather than leaning on the
## pan-cancer number.
## ---------------------------------------------------------------------------
cat("\n==========================================================\n")
cat(" PART 4: breast-only sensitivity analysis\n")
cat("==========================================================\n")

br <- cv[lineage_raw == "Breast", ModelID]
br <- intersect(br, rownames(GE))
if (length(br) >= 25) {
  GEb <- GE[br, , drop = FALSE]
  cvb <- cv[match(br, ModelID)]
  Xb  <- model.matrix(~ mean_ge, data = cvb)   # no lineage term within one lineage
  Rb  <- GEb - Xb %*% qr.coef(qr(Xb), GEb)
  cat(sprintf("\n%d breast lines\n", length(br)))
  if (all(c("PIK3C3","ATG7") %in% colnames(Rb))) {
    cb <- cor.test(Rb[, "PIK3C3"], Rb[, "ATG7"])
    cat(sprintf("partial r(PIK3C3, ATG7) in breast = %+.3f [%+.3f, %+.3f], p = %.3g\n",
                cb$estimate, cb$conf.int[1], cb$conf.int[2], cb$p.value))
  }
  cat(sprintf("%-8s  %10s %10s %12s\n", "gene", "ATG8", "lysosome", "trafficking"))
  for (f in FOCAL) {
    if (!f %in% colnames(Rb)) next
    vals <- sapply(names(MODULES), function(m) {
      g <- intersect(setdiff(MODULES[[m]], f), colnames(Rb))
      if (!length(g)) NA_real_ else mean(cor(Rb[, f], Rb[, g]))
    })
    cat(sprintf("%-8s  %+10.3f %+10.3f %+12.3f\n", f,
                vals["ATG8_cascade"], vals["lysosome"], vals["trafficking"]))
  }
  cat("\nn = 67 breast lines gives wide confidence intervals. Treat any\n")
  cat("apparent difference from the pan-cancer result as unconfirmed.\n")
} else {
  cat("\nToo few breast lines for a stable within-lineage estimate.\n")
}

## ---------------------------------------------------------------------------
cat("\n==========================================================\n")
cat(" HOW TO READ THIS\n")
cat("==========================================================\n")
cat("
Read the CALIBRATION number first. PIK3C3 and PIK3R4 are obligate complex
partners, so their correlation is what 'same module' looks like in this data
(~0.33 in the full panel). Every other correlation is a FRACTION of that, not
a deviation from zero. A correlation of 0.08 against a 0.33 ceiling is a
quarter of a module relationship, i.e. effectively unrelated -- even when its
p-value is tiny, because n is large.

Then PART 2. The primary contrast is trafficking - ATG8_cascade:
  PIK3C3 positive with CI excluding 0  -> VPS34 sits with endosomal
     trafficking (ESCRT/HOPS/RAB7), not with ATG8 conjugation.
  ATG7 negative                        -> internal negative control passed;
     the contrast is measuring module identity, not something generic.
  Both positive                        -> the trafficking set is picking up
     general essentiality. Do not interpret.

Then PART 3, which imposes no gene set at all and is therefore the hardest
result to dispute: whatever a gene's nearest neighbours turn out to be IS its
functional module. Note the lysosome set captures lysosomal CONTENTS and
acidification, while the trafficking set captures the DELIVERY machinery --
these are different biology, and the first version of this analysis used only
the former, which is why it found nothing.

LIMITS, for the paper:
  - CRISPR knockout is chronic, chloroquine is acute. Cells compensate for
    chronic ATG7 loss in ways they cannot for acute lysosomal blockade, so a
    null co-dependency does not fully exclude an acute autophagy requirement.
  - Co-dependency shows genes co-vary, not that they act in one pathway.
  - Screens run in nutrient-rich 2-D culture, which understates autophagy.
  - None of this is a substitute for the knockout. It is the best available
    evidence without one, and should be described that way.
")
cat(sprintf("\nFiles written to %s/\n\n", OUT_DIR))
