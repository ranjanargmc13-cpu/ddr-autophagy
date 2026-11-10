#!/usr/bin/env Rscript
## prism_drugs.R -------------------------------------------------------------
## DOES VPS34 DEPENDENCY PREDICT CHLOROQUINE SENSITIVITY?
##
##   Rscript prism_drugs.R [/path/to/depmap]
##
## This is the direct test of your bench result. Your three drugs are all in
## the PRISM Repurposing screen, so you can ask across ~500-900 cell lines:
##
##   Do the lines that DEPEND on VPS34 happen to be the lines chloroquine KILLS?
##
##   yes -> your drug effect and your genetic finding are the same phenomenon
##   no  -> chloroquine is doing something else, and the VPS34 story is
##          a separate observation that does not explain your synergy
##
## The ATG7 contrast is the point of the whole thing: if chloroquine
## sensitivity tracks PIK3C3 dependency but NOT ATG7 dependency, that is
## evidence chloroquine acts through the VPS34/endolysosomal route rather than
## canonical ATG8-dependent autophagy -- the claim you cannot make at the bench.
##
## DOWNLOAD (https://depmap.org/portal/data_page/?tab=allData, "Drug screens")
## Either of these works; the script detects which you have:
##   A) Repurposing_Public_*_Extended_Primary_Data_Matrix.csv
##      + Repurposing_Public_*_Extended_Primary_Compound_List.csv
##   B) secondary-screen-dose-response-curve-parameters.csv   (older, simplest)
## Put it in the same folder as the other DepMap files.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))
setDTthreads(1)

args     <- commandArgs(trailingOnly = TRUE)
DATA_DIR <- path.expand(if (length(args)) args[1] else "~/depmap")
OUT_DIR  <- "results_prism"
dir.create(OUT_DIR, showWarnings = FALSE)

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
clean <- function(x) trimws(sub("\\s*\\(\\d+\\)$", "", x))

DRUGS <- c("chloroquine", "hydroxychloroquine", "olaparib", "metformin",
           "bafilomycin", "niclosamide", "cisplatin", "talazoparib")
GENES <- c("PIK3C3", "ATG7", "ATG5", "ATG3", "RB1CC1", "BECN1", "PIK3R4",
           "RAB7A", "VPS33A", "HGS", "PARP1", "ATR")

## ---------------------------------------------------------------------------
## dependencies
## ---------------------------------------------------------------------------
GE_FILE <- file.path(DATA_DIR, "CRISPRGeneEffect.csv")
if (!file.exists(GE_FILE)) stop("Missing CRISPRGeneEffect.csv in ", DATA_DIR)
message("[1/3] reading dependencies ...")
hdr <- names(fread(GE_FILE, nrows = 0)); sym <- clean(hdr)
iw  <- which(sym %in% GENES)
set.seed(42)
ir  <- sort(sample(setdiff(seq_along(hdr)[-1], iw), min(1500, length(hdr) - 2)))
GE  <- fread(GE_FILE, select = sort(unique(c(1L, iw, ir))), showProgress = FALSE)
ids <- as.character(GE[[1]]); GE[, (1) := NULL]
M <- as.matrix(GE); rm(GE); invisible(gc())
colnames(M) <- make.unique(clean(colnames(M))); rownames(M) <- ids
qc <- rowMeans(M, na.rm = TRUE)

md <- fread(file.path(DATA_DIR, "Model.csv"), showProgress = FALSE)
md <- md[, intersect(c("ModelID","StrippedCellLineName","OncotreeLineage"),
                     names(md)), with = FALSE]

## ---------------------------------------------------------------------------
## drug response: detect which PRISM format is present
## ---------------------------------------------------------------------------
message("[2/3] locating PRISM data ...")
fs <- list.files(DATA_DIR, full.names = TRUE)
f_sec <- grep("secondary-screen-dose-response-curve-parameters", fs,
              ignore.case = TRUE, value = TRUE)
f_mat <- grep("Repurposing.*Data_Matrix", fs, ignore.case = TRUE, value = TRUE)
f_cmp <- grep("Repurposing.*Compound_List",  fs, ignore.case = TRUE, value = TRUE)

drug <- NULL   # long: ModelID, drug, sensitivity (LOWER = more sensitive)

if (length(f_sec)) {
  message("       format B: ", basename(f_sec[1]))
  sc <- fread(f_sec[1], showProgress = FALSE)
  nm <- names(sc)
  idc <- intersect(c("depmap_id","DepMap_ID","ModelID"), nm)[1]
  dnc <- intersect(c("name","Drug.Name","drug_name"), nm)[1]
  vlc <- intersect(c("auc","AUC","ic50","IC50"), nm)[1]
  if (any(is.na(c(idc, dnc, vlc)))) stop("Unexpected columns in secondary-screen file.")
  sc <- sc[, .(ModelID = get(idc), drug = tolower(get(dnc)),
               sensitivity = as.numeric(get(vlc)))]
  drug <- sc[grepl(paste(DRUGS, collapse = "|"), drug)]
  message(sprintf("       using '%s' as sensitivity (lower = more sensitive)", vlc))

} else if (length(f_mat) && length(f_cmp)) {
  message("       format A: ", basename(f_mat[1]))
  cmp <- fread(f_cmp[1], showProgress = FALSE)
  nmc <- intersect(c("Drug.Name","name","drug_name"), names(cmp))[1]
  idc <- intersect(c("IDs","broad_id","BroadID","column_name"), names(cmp))[1]
  if (any(is.na(c(nmc, idc)))) stop("Unexpected columns in compound list.")
  hit <- cmp[grepl(paste(DRUGS, collapse = "|"), tolower(get(nmc)))]
  if (!nrow(hit)) stop("None of the drugs found in the compound list.")
  message(sprintf("       %d matching compound entries", nrow(hit)))
  mh  <- names(fread(f_mat[1], nrows = 0))
  sel <- which(sapply(mh, function(h)
            any(sapply(hit[[idc]], function(b) grepl(b, h, fixed = TRUE)))))
  if (!length(sel)) stop("Compound IDs not found among matrix columns.")
  DM  <- fread(f_mat[1], select = sort(unique(c(1L, sel))), showProgress = FALSE)
  dmi <- as.character(DM[[1]]); DM[, (1) := NULL]
  lf  <- melt(cbind(data.table(ModelID = dmi), DM), id.vars = "ModelID",
              variable.name = "col", value.name = "sensitivity", na.rm = TRUE)
  lf[, drug := tolower(sapply(as.character(col), function(cc) {
        m <- hit[sapply(hit[[idc]], function(b) grepl(b, cc, fixed = TRUE))]
        if (nrow(m)) as.character(m[[nmc]][1]) else NA_character_ }))]
  drug <- lf[!is.na(drug), .(ModelID, drug, sensitivity)]
  message("       NOTE: matrix values are log fold-change (lower = more killing)")

} else {
  stop("No PRISM file found in ", DATA_DIR, "\n",
       "Download from https://depmap.org/portal/data_page/?tab=allData ",
       "(Drug screens):\n",
       "  secondary-screen-dose-response-curve-parameters.csv   (simplest)\n",
       "or Repurposing_Public_*_Extended_Primary_Data_Matrix.csv\n",
       "   + Repurposing_Public_*_Extended_Primary_Compound_List.csv")
}

drug <- drug[is.finite(sensitivity)]
drug <- drug[, .(sensitivity = median(sensitivity)), by = .(ModelID, drug)]
cat("\ndrugs found:\n"); print(drug[, .(n_lines = .N), by = drug][order(-n_lines)])

## ---------------------------------------------------------------------------
## correlate dependency with drug sensitivity
## ---------------------------------------------------------------------------
message("\n[3/3] testing dependency vs sensitivity ...")
cat("\n==========================================================\n")
cat(" DEPENDENCY vs DRUG SENSITIVITY\n")
cat("==========================================================\n")
cat("\nBoth scales: LOWER = more dependent / more sensitive.\n")
cat("So a POSITIVE r means dependent lines are also sensitive lines,\n")
cat("i.e. the drug and the gene are hitting the same vulnerability.\n")

res <- rbindlist(lapply(unique(drug$drug), function(dd) {
  dsub <- drug[drug == dd]
  rbindlist(lapply(intersect(GENES, colnames(M)), function(g) {
    d <- data.table(ModelID = rownames(M), dep = M[, g], mean_ge = qc)
    d <- merge(d, dsub[, .(ModelID, sensitivity)], by = "ModelID")
    d <- merge(d, md, by = "ModelID")
    d <- d[is.finite(dep) & is.finite(sensitivity)]
    if (nrow(d) < 40) return(NULL)
    d[, lin := factor(ifelse(OncotreeLineage %in%
         names(which(table(OncotreeLineage) >= 8)), OncotreeLineage, "Other"))]
    X <- model.matrix(~ lin + mean_ge, data = d)
    a <- residuals(lm(d$dep ~ X - 1)); b <- residuals(lm(d$sensitivity ~ X - 1))
    ct <- suppressWarnings(cor.test(a, b, method = "spearman"))
    data.table(drug = dd, gene = g, n = nrow(d),
               rho = unname(ct$estimate), p = ct$p.value)
  }))
}))
if (!nrow(res)) stop("No drug/gene pair had enough overlapping lines.")
res[, fdr := p.adjust(p, "BH")]
setorder(res, drug, -rho)
fwrite(res, file.path(OUT_DIR, "dependency_vs_drug.csv"))

for (dd in unique(res$drug)) {
  cat(sprintf("\n--- %s ---\n", dd))
  print(res[drug == dd, .(gene, n, rho, p, fdr)], digits = 3)
}

## ---------------------------------------------------------------------------
## THE KEY CONTRAST
## ---------------------------------------------------------------------------
cat("\n==========================================================\n")
cat(" KEY CONTRAST: chloroquine vs PIK3C3 and ATG7\n")
cat("==========================================================\n\n")
cq <- res[grepl("chloroquine", drug) & gene %in% c("PIK3C3", "ATG7")]
if (nrow(cq)) {
  print(cq[, .(drug, gene, n, rho, p, fdr)], digits = 3)
  for (dd in unique(cq$drug)) {
    a <- cq[drug == dd & gene == "PIK3C3", rho]
    b <- cq[drug == dd & gene == "ATG7",   rho]
    if (length(a) && length(b)) {
      cat(sprintf("\n%s:  PIK3C3 rho = %+.3f   ATG7 rho = %+.3f\n", dd, a, b))
      cat(if (a > 0.1 && abs(b) < 0.1)
        "  >>> Sensitivity tracks VPS34 dependency but NOT ATG7 dependency.\n  >>> Supports an endolysosomal rather than ATG8-autophagy mechanism.\n"
       else if (a > 0.1 && b > 0.1)
        "  >>> Tracks both. Cannot separate the two routes with this test.\n"
       else if (abs(a) < 0.1 && abs(b) < 0.1)
        "  >>> Tracks neither. Chloroquine sensitivity across lines is driven\n  >>> by something other than either dependency -- an important negative.\n"
       else "  >>> Mixed. Read the full table before concluding.\n")
    }
  }
} else cat("Chloroquine not found in this PRISM release.\n")

cat("\nLIMITS: PRISM is a single-agent screen, so it cannot show synergy -- it\n")
cat("tests whether the same cells are vulnerable to both. Sensitivity is also\n")
cat("measured at fixed doses that may not match yours. And chloroquine has\n")
cat("off-target lysosomal effects independent of any of these genes.\n")
cat(sprintf("\nFiles written to %s/\n\n", OUT_DIR))
