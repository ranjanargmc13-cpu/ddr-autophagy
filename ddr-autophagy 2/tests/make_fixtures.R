## make_fixtures.R -----------------------------------------------------------
## Synthetic DepMap-shaped data with a KNOWN planted effect.
##
## Why this exists: if you run the pipeline on real data first and it returns a
## significant result, you cannot tell whether the biology is real or the code
## is wrong. Running it on synthetic data where you already know the answer
## separates those two possibilities. Run --selftest before --data_dir, always.
##
## The fixture reproduces DepMap's schema exactly: ACH-###### model IDs,
## "SYMBOL (ENTREZ)" column naming, Chronos-like value distributions, lineage
## imbalance, and per-line screen-quality variation.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

make_synthetic_depmap <- function(n_lines = 900,
                                  n_bg_genes = 1200,
                                  planted_rho = -0.35,
                                  seed = 42) {
  set.seed(seed)

  ids <- sprintf("ACH-%06d", seq_len(n_lines))

  ## ---- lineages, deliberately imbalanced like the real thing -------------
  lineages <- c("Lung", "Breast", "Colorectal", "Ovary/Fallopian Tube",
                "Pancreas", "Skin", "CNS/Brain", "Lymphoid", "Myeloid",
                "Esophagus/Stomach", "Kidney", "Bladder/Urinary Tract",
                "Prostate", "Soft Tissue", "Thyroid", "Eye")
  probs <- c(.14,.10,.09,.07,.06,.07,.08,.08,.06,.06,.05,.04,.03,.04,.02,.01)
  lin <- sample(lineages, n_lines, replace = TRUE, prob = probs)

  ## ---- per-line screen quality (the key confounder) ----------------------
  ## Real screens vary in efficiency. This shifts ALL genes in a line together
  ## and will manufacture co-dependency if the analysis fails to control for it.
  ## We plant it deliberately so the self-test proves the control works.
  qc <- rnorm(n_lines, 0, 0.12)

  ## ---- latent biology ----------------------------------------------------
  ## A true latent "autophagy reliance" factor, plus a lineage offset so that
  ## some tumour types are genuinely more autophagy-dependent.
  lin_offset <- setNames(rnorm(length(lineages), 0, 0.25), lineages)
  latent_auto <- rnorm(n_lines, 0, 1) + lin_offset[lin]
  latent_auto <- as.numeric(scale(latent_auto))

  ## ---- gene sets we must emit -------------------------------------------
  auto_genes <- unique(c(
    c("ULK1","ATG13","RB1CC1","ATG101"),
    c("PIK3C3","PIK3R4","BECN1","ATG14","NRBF2"),
    c("ATG9A","ATG2A","ATG2B","WIPI2"),
    c("ATG7","ATG3","ATG10","ATG12","ATG5","ATG16L1"),
    c("ATG4B","MAP1LC3B","GABARAPL1")))
  lyso_genes <- c("ATP6V1A","ATP6V1B2","ATP6V0C","ATP6AP1","ATP6AP2","LAMP1",
                  "LAMP2","CTSB","CTSD","CTSL","TFEB","MCOLN1","GBA","NPC1","PSAP")
  ddr_genes <- unique(c(
    "PARP1","PARP2","PARG","XRCC1","POLB","FEN1","LIG3",
    "ATR","CHEK1","WEE1","PKMYT1","TOPBP1","CLSPN","RPA1","RPA2",
    "BRCA1","BRCA2","PALB2","RAD51","RAD51C","BARD1","RAD54L","BRIP1",
    "PRKDC","LIG4","XRCC4","XRCC5","XRCC6","NHEJ1","POLL","POLQ",
    "ATM","CHEK2","MDC1","NBN","MRE11","RAD50","FANCD2","FANCI",
    "USP1","RAD18","UBE2T","SLFN11","TP53","CDKN1A","TP53BP1","SHLD2"))
  neg_genes <- c("RPL13A","ACTB","GAPDH","TUBB","PPIA","B2M")
  rep_genes <- c("MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","PCNA","CDC6","CDT1",
                 "CDC45","GINS1","GINS2","GINS3","GINS4","RRM1","RRM2","TYMS",
                 "TK1","DHFR","E2F1","CCNE1","CCNE2","CDK2","CLSPN","POLA1","PRIM1")
  hrr_genes <- c("BRCA1","BRCA2","PALB2","BARD1","BRIP1","RAD51B","RAD51C",
                 "RAD51D","RAD50","RAD54L","XRCC2","XRCC3","NBN","MRE11","ATM",
                 "ATR","CHEK1","CHEK2","FANCA","FANCB","FANCC","FANCD2","FANCE",
                 "FANCF","FANCG","FANCI","FANCL","FANCM","PPP2R2A","ABRAXAS1","GEN1")
  extra <- c("ARID1A","KRAS","STK11","KEAP1","SMARCA4","RB1","MYC","MYCN",
             "ERBB2","CDKN2A","PTEN")
  bg <- sprintf("BG%04d", seq_len(n_bg_genes))

  all_genes <- unique(c(auto_genes, lyso_genes, ddr_genes, neg_genes,
                        rep_genes, hrr_genes, extra, bg))
  ng <- length(all_genes)

  ## ---- build the gene-effect matrix -------------------------------------
  GE <- matrix(rnorm(n_lines * ng, 0, 0.18), n_lines, ng,
               dimnames = list(ids, all_genes))
  GE <- GE + qc                      # screen-quality shift, all genes together

  ## autophagy genes: loaded on the latent factor.
  ## Chronos is negative for dependency, so a HIGH latent value must push the
  ## score DOWN. Loadings are heterogeneous, as in reality.
  for (g in auto_genes) {
    load_g <- runif(1, 0.10, 0.30)
    base   <- runif(1, -0.18, -0.02)
    GE[, g] <- GE[, g] + base - load_g * latent_auto
  }

  ## lysosome genes: partly shared latent factor, partly their own.
  ## This mirrors the real partial overlap between autophagy and lysosomal
  ## dependency, which is what makes CQ/HCQ interpretation hard.
  latent_lyso <- as.numeric(scale(0.45 * latent_auto + rnorm(n_lines, 0, 1)))
  for (g in lyso_genes) {
    GE[, g] <- GE[, g] + runif(1, -0.35, -0.10) - runif(1, 0.10, 0.25) * latent_lyso
  }

  ## PLANTED EFFECT: a subset of DDR genes is coupled to the latent autophagy
  ## factor with the requested correlation. These are the true positives the
  ## pipeline must recover.
  planted_targets <- c("ATR", "CHEK1", "WEE1", "POLQ", "PARP1")
  for (g in planted_targets) {
    ## planted_rho is the target correlation between autophagy VULNERABILITY
    ## (= -score, so ~ +latent) and this gene's Chronos effect.
    ## Chronos_g = rho * latent + sqrt(1-rho^2) * noise, scaled to Chronos units.
    r <- planted_rho
    z <- r * latent_auto + sqrt(max(0, 1 - r^2)) * rnorm(n_lines)
    GE[, g] <- runif(1, -0.40, -0.15) + 0.30 * as.numeric(scale(z)) + qc
  }

  ## common essentials, so the matrix has a realistic dependency distribution
  for (g in c(neg_genes, "RPA1", "PCNA")) {
    GE[, g] <- GE[, g] + runif(1, -1.4, -0.9)
  }

  ## ---- mutation matrix (0/1/2 damaging counts) --------------------------
  mut_genes <- unique(c(hrr_genes, "TP53","ATM","ATR","BRCA1","BRCA2","ARID1A",
                        "KRAS","STK11","KEAP1","SMARCA4","RB1"))
  MUT <- matrix(0L, n_lines, length(mut_genes),
                dimnames = list(ids, mut_genes))
  for (g in mut_genes) {
    ## TP53 is mutated in roughly half of cell lines; others are rarer.
    pr <- if (g == "TP53") 0.55 else if (g %in% c("ATM","BRCA2","ARID1A")) 0.08 else 0.03
    MUT[, g] <- rbinom(n_lines, 1, pr)
  }
  ## Make HRR mutation mildly predictive of autophagy reliance, so TEST B has
  ## something recoverable too (weaker than the planted co-dependency).
  hrr_present <- intersect(hrr_genes, mut_genes)
  boost <- which(latent_auto > quantile(latent_auto, 0.70))
  for (g in sample(hrr_present, 6)) {
    MUT[boost, g] <- pmax(MUT[boost, g],
                          rbinom(length(boost), 1, 0.22))
  }

  ## ---- copy number (log2 relative, neutral ~ 1.0) ----------------------
  cn_genes <- unique(c("CCNE1","MYC","MYCN","ERBB2","CDKN2A","PTEN","ATM",
                       "RB1","SMARCA4", bg[1:200]))
  CN <- matrix(rnorm(n_lines * length(cn_genes), 1.0, 0.25),
               n_lines, length(cn_genes), dimnames = list(ids, cn_genes))
  amp_lines <- sample(n_lines, round(0.12 * n_lines))
  CN[amp_lines, "CCNE1"] <- CN[amp_lines, "CCNE1"] + runif(length(amp_lines), 0.8, 2.0)

  ## ---- expression (log2 TPM+1) -----------------------------------------
  expr_genes <- unique(c(rep_genes, auto_genes, lyso_genes, bg[1:400]))
  EX <- matrix(rnorm(n_lines * length(expr_genes), 5, 1.5),
               n_lines, length(expr_genes), dimnames = list(ids, expr_genes))
  ## replication-stress genes track the latent factor a little
  for (g in rep_genes) EX[, g] <- EX[, g] + 0.45 * latent_auto
  EX[EX < 0] <- 0

  ## ---- assemble in load_depmap()'s output shape -------------------------
  model <- data.table(
    ModelID = ids,
    CellLineName = paste0("SYN-", seq_len(n_lines)),
    StrippedCellLineName = paste0("SYN", seq_len(n_lines)),
    OncotreeLineage = lin,
    OncotreePrimaryDisease = paste(lin, "carcinoma"),
    OncotreeSubtype = paste(lin, "subtype")
  )

  out <- list(
    gene_effect = GE,
    model = model,
    mut = MUT,
    cn = CN,
    expr = EX,
    line_qc = data.table(ModelID = ids,
                         mean_gene_effect = rowMeans(GE, na.rm = TRUE),
                         n_measured = ncol(GE))
  )
  attr(out, "planted_rho") <- planted_rho
  attr(out, "planted_targets") <- planted_targets
  attr(out, "latent_auto") <- latent_auto

  message(sprintf(paste0(
    "[fixture] synthetic DepMap: %d lines x %d genes\n",
    "[fixture]   planted coupling rho = %+.2f on: %s\n",
    "[fixture]   screen-quality confounder planted (SD 0.12)\n",
    "[fixture]   TEST A must recover the planted targets near the top,\n",
    "[fixture]   and the negative-control set must recover nothing."),
    n_lines, ng, planted_rho, paste(planted_targets, collapse = ", ")))

  out
}

#' Write the fixture out as real CSV files, so you can also test the LOADERS
#' (schema handling, "SYMBOL (ENTREZ)" parsing, ModelID detection) and not just
#' the statistics.
write_synthetic_csvs <- function(dir = "synthetic_depmap", ...) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  dm <- make_synthetic_depmap(...)

  ## re-attach fake entrez IDs to mimic DepMap's column naming
  entrezify <- function(m) {
    cn <- sprintf("%s (%d)", colnames(m), seq_len(ncol(m)) + 1000L)
    dt <- data.table(ModelID = rownames(m))
    body <- as.data.table(m); setnames(body, cn)
    cbind(dt, body)
  }
  fwrite(entrezify(dm$gene_effect), file.path(dir, "CRISPRGeneEffect.csv"))
  fwrite(entrezify(dm$mut),  file.path(dir, "OmicsSomaticMutationsMatrixDamaging.csv"))
  fwrite(entrezify(dm$cn),   file.path(dir, "OmicsCNGene.csv"))
  fwrite(entrezify(dm$expr), file.path(dir, "OmicsExpressionProteinCodingGenesTPMLogp1.csv"))
  fwrite(dm$model, file.path(dir, "Model.csv"))
  message(sprintf("[fixture] CSVs written to %s/ -- test the loaders with:", dir))
  message(sprintf("          Rscript run_all.R --data_dir=%s --null_B=50", dir))
  invisible(dir)
}
