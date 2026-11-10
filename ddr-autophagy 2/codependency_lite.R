#!/usr/bin/env Rscript
## codependency_lite.R -------------------------------------------------------
## Memory-light version of codependency_module.R. Same headline result,
## ~15 MB instead of ~1 GB.
##
##   Rscript codependency_lite.R                 # uses ~/depmap
##   Rscript codependency_lite.R /path/to/depmap
##
## HOW IT SAVES MEMORY
## The full script loads all ~18,000 genes. This one reads the header, works out
## which columns it actually needs, and asks fread for only those -- about 60
## genes of interest plus 1,500 random genes used solely to estimate each line's
## screen quality. rowMeans over 1,500 random genes correlates >0.99 with the
## mean over all 18,000, so the covariate is essentially unchanged.
##
## The 1,500 random genes do double duty: they also give an empirical null
## distribution of partial correlations, so you can see whether
## r(PIK3C3, ATG7) is unusual or just typical of any two genes.
##
## WHAT IT DROPS
## The genome-wide nearest-co-dependency scan (Part 3 of the full script). Run
## that later on a machine with more RAM, or ask for a chunked version.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))
setDTthreads(1)   # keep memory predictable on a small machine

args     <- commandArgs(trailingOnly = TRUE)
DATA_DIR <- path.expand(if (length(args)) args[1] else "~/depmap")
OUT_DIR  <- "results_codependency"
N_RANDOM <- 1500

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

GE_FILE <- file.path(DATA_DIR, "CRISPRGeneEffect.csv")
MD_FILE <- file.path(DATA_DIR, "Model.csv")
for (f in c(GE_FILE, MD_FILE)) if (!file.exists(f)) stop("Missing: ", f)

ATG8_CASCADE <- AUTOPHAGY_PRIMARY
LYSO_MODULE  <- LYSOSOME_SET
FOCAL        <- c("PIK3C3", "RB1CC1", "BECN1", "ATG7")
WANTED <- unique(c(FOCAL, ATG8_CASCADE, LYSO_MODULE,
                   "ULK1", "ATG13", "ATG101", "ATG14", "PIK3R4",
                   "ATG2A", "ATG9A", "WIPI2", "MAP1LC3B", "GABARAPL1"))

## ---- 1. read the header only, to map gene -> column index ----------------
message("[lite] reading header ...")
hdr <- names(fread(GE_FILE, nrows = 0))
sym <- sub("\\s*\\(\\d+\\)$", "", hdr); sym <- trimws(sym)
message(sprintf("[lite] file has %d columns", length(hdr)))

idx_wanted <- which(sym %in% WANTED)
found <- sym[idx_wanted]
missing <- setdiff(WANTED, found)
if (length(missing))
  message("[lite] not in file: ", paste(missing, collapse = ", "))

set.seed(42)
pool <- setdiff(seq_along(hdr)[-1], idx_wanted)      # exclude col 1 (ModelID)
idx_random <- sort(sample(pool, min(N_RANDOM, length(pool))))

idx <- sort(unique(c(1L, idx_wanted, idx_random)))
message(sprintf("[lite] reading %d of %d columns (%.1f%%) ...",
                length(idx), length(hdr), 100 * length(idx) / length(hdr)))

GE <- fread(GE_FILE, select = idx, showProgress = FALSE)
ids <- as.character(GE[[1]]); GE[, (1) := NULL]
g   <- sub("\\s*\\(\\d+\\)$", "", names(GE)); g <- trimws(g)
M   <- as.matrix(GE); rm(GE); invisible(gc())
storage.mode(M) <- "double"
dimnames(M) <- list(ids, make.unique(g))
message(sprintf("[lite] loaded %d lines x %d genes (%.1f MB)",
                nrow(M), ncol(M), as.numeric(object.size(M)) / 1e6))

## ---- 2. covariates --------------------------------------------------------
rand_cols <- which(dimnames(M)[[2]] %in% sym[idx_random])
if (length(rand_cols) < 200) rand_cols <- seq_len(ncol(M))   # fallback
qc <- rowMeans(M[, rand_cols, drop = FALSE], na.rm = TRUE)

md <- fread(MD_FILE, select = c("ModelID", "OncotreeLineage"), showProgress = FALSE)
cv <- data.table(ModelID = ids, mean_ge = qc)
cv <- merge(cv, md, by = "ModelID", sort = FALSE)
cv[, lineage := factor(ifelse(
  OncotreeLineage %in% names(which(table(OncotreeLineage) >= 10)),
  OncotreeLineage, "Other"))]
cv <- cv[!is.na(mean_ge)]
M  <- M[cv$ModelID, , drop = FALSE]
message(sprintf("[lite] %d lines, %d lineages after filtering",
                nrow(cv), nlevels(cv$lineage)))

## ---- 3. residualize on lineage + screen quality --------------------------
nas <- which(is.na(M), arr.ind = TRUE)
if (nrow(nas)) { cm <- colMeans(M, na.rm = TRUE); M[nas] <- cm[nas[, "col"]] }
X <- model.matrix(~ lineage + mean_ge, data = cv)
R <- M - X %*% qr.coef(qr(X), M)
rm(M); invisible(gc())
message("[lite] residualization done\n")

## ---- 4. HEADLINE ----------------------------------------------------------
cat("==========================================================\n")
cat(" HEADLINE: is VPS34 dependency autophagy-mediated?\n")
cat("==========================================================\n\n")
if (all(c("PIK3C3", "ATG7") %in% colnames(R))) {
  ct <- cor.test(R[, "PIK3C3"], R[, "ATG7"])
  cat(sprintf("partial r(PIK3C3 dep, ATG7 dep) = %+.3f  [95%% CI %+.3f, %+.3f]\n",
              ct$estimate, ct$conf.int[1], ct$conf.int[2]))
  cat(sprintf("p = %.3g,  n = %d cell lines\n\n", ct$p.value, nrow(R)))
  cat(if (abs(ct$estimate) < 0.15)
        ">>> |r| < 0.15: VPS34 and ATG8 conjugation are SEPARATE modules.\n>>> Your combination is endolysosomal, not autophagic.\n"
      else if (ct$estimate > 0.30)
        ">>> r > 0.30: they DO co-vary. An autophagy interpretation stays open.\n"
      else ">>> Intermediate. Read the module comparison below before concluding.\n")
} else cat("PIK3C3 or ATG7 absent from the file.\n")

## ---- 5. module comparison -------------------------------------------------
cat("\n==========================================================\n")
cat(" MODULE COMPARISON (partial r, lineage + QC adjusted)\n")
cat("==========================================================\n")
res <- rbindlist(lapply(FOCAL, function(f) {
  if (!f %in% colnames(R)) return(NULL)
  a <- intersect(setdiff(ATG8_CASCADE, f), colnames(R))
  l <- intersect(setdiff(LYSO_MODULE,  f), colnames(R))
  data.table(focal = f,
             mean_r_ATG8 = mean(cor(R[, f], R[, a])),
             mean_r_lyso = mean(cor(R[, f], R[, l])),
             delta = mean(cor(R[, f], R[, l])) - mean(cor(R[, f], R[, a])))
}))
cat("\ndelta > 0 means the gene behaves like the LYSOSOME module\n\n")
print(res, digits = 3)
fwrite(res, file.path(OUT_DIR, "lite_module_comparison.csv"))

## per-gene detail
det <- rbindlist(lapply(FOCAL, function(f) {
  if (!f %in% colnames(R)) return(NULL)
  mk <- function(gs, lab) {
    gs <- intersect(setdiff(gs, f), colnames(R)); if (!length(gs)) return(NULL)
    data.table(focal = f, module = lab, gene = gs,
               partial_r = as.numeric(cor(R[, f], R[, gs])))
  }
  rbindlist(list(mk(ATG8_CASCADE, "ATG8_cascade"), mk(LYSO_MODULE, "lysosome")))
}))
fwrite(det, file.path(OUT_DIR, "lite_module_detail.csv"))
cat("\n--- per-gene detail for PIK3C3 ---\n")
print(det[focal == "PIK3C3"][order(module, -abs(partial_r)),
                             .(module, gene, partial_r)], digits = 3)

## ---- 6. empirical null from the random genes ------------------------------
cat("\n==========================================================\n")
cat(" EMPIRICAL NULL (random gene pairs)\n")
cat("==========================================================\n")
rc <- intersect(dimnames(R)[[2]][dimnames(R)[[2]] %in% sym[idx_random]], colnames(R))
if (length(rc) > 300 && "PIK3C3" %in% colnames(R)) {
  set.seed(7)
  s  <- sample(rc, 300)
  nullr <- as.numeric(cor(R[, "PIK3C3"], R[, s]))
  obs <- if ("ATG7" %in% colnames(R)) cor(R[, "PIK3C3"], R[, "ATG7"]) else NA
  cat(sprintf("\nPIK3C3 vs 300 random genes: mean r = %+.3f, SD = %.3f\n",
              mean(nullr), sd(nullr)))
  cat(sprintf("PIK3C3 vs ATG7            : r = %+.3f  (%.2f SD from the null mean)\n",
              obs, (obs - mean(nullr)) / sd(nullr)))
  cat(sprintf("PIK3C3 vs lysosome module : mean r = %+.3f  (%.2f SD from null)\n",
              res[focal=="PIK3C3", mean_r_lyso],
              (res[focal=="PIK3C3", mean_r_lyso] - mean(nullr)) / sd(nullr)))
  cat("\nIf ATG7 sits inside the random-gene null but the lysosome module sits\n")
  cat("well outside it, that is the result -- VPS34 travels with the lysosome,\n")
  cat("not with ATG8 conjugation.\n")
}

cat(sprintf("\nFiles written to %s/\n", OUT_DIR))
cat("\nLimits: knockout is chronic, chloroquine acute; co-dependency shows\n")
cat("co-variation, not shared pathway; screens run nutrient-rich, which\n")
cat("understates autophagy. This is the best evidence available without the\n")
cat("knockout -- describe it that way.\n\n")
