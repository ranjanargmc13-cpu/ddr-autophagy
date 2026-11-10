## 01_genesets.R -------------------------------------------------------------
## Curated gene sets for the DDR x autophagy analysis.
## Every set is annotated with WHY it is included, because reviewers will ask.
## ---------------------------------------------------------------------------

## --- Autophagy machinery ---------------------------------------------------
## Split by functional module so we can test module-level consistency.
## A real synthetic-lethal signal should appear in SEVERAL modules, not one.
AUTOPHAGY_MODULES <- list(
  initiation  = c("ULK1", "ATG13", "RB1CC1", "ATG101"),
  pi3kc3      = c("PIK3C3", "PIK3R4", "BECN1", "ATG14", "NRBF2"),
  atg9_traffic= c("ATG9A", "ATG2A", "ATG2B", "WIPI2"),
  conjugation = c("ATG7", "ATG3", "ATG10", "ATG12", "ATG5", "ATG16L1"),
  lc3_process = c("ATG4B", "MAP1LC3B", "GABARAPL1")
)

## PRIMARY score uses the conjugation cascade only.
## Rationale: ATG7/3/10/12/5/16L1 form the E1-E2-E3-like cascade that is
## essentially dedicated to ATG8 lipidation. They are the most
## autophagy-SPECIFIC nodes available, which matters because the alternatives
## are pleiotropic (see AUTOPHAGY_PLEIOTROPIC below). Using a specific set
## costs statistical power but buys interpretability.
AUTOPHAGY_PRIMARY <- AUTOPHAGY_MODULES$conjugation

## All modules pooled -- used for the sensitivity analysis.
AUTOPHAGY_ALL <- unique(unlist(AUTOPHAGY_MODULES, use.names = FALSE))

## Genes with substantial NON-autophagy functions. If your signal is driven
## by these, it is probably not autophagy. Reported separately, never silently
## pooled into the primary score.
##   RB1CC1  - focal adhesion kinase family interacting protein, TSC/mTOR
##   PIK3C3  - endosomal PI3P, LAP, general lysosome biogenesis
##   ATG9A   - broad membrane trafficking
##   BECN1   - apoptosis crosstalk via BCL2
##   ATG5/7  - LC3-associated phagocytosis, secretory autophagy
AUTOPHAGY_PLEIOTROPIC <- c("RB1CC1", "PIK3C3", "ATG9A", "BECN1")

## Lysosome / autophagosome-lysosome fusion. NOT autophagy-specific, but this
## is what chloroquine and hydroxychloroquine actually hit. Included so you can
## test whether "autophagy" dependency and "lysosome" dependency diverge --
## they often do, and that divergence is the honest explanation for why HCQ
## trials disappoint.
LYSOSOME_SET <- c("ATP6V1A", "ATP6V1B2", "ATP6V0C", "ATP6AP1", "ATP6AP2",
                  "LAMP1", "LAMP2", "CTSB", "CTSD", "CTSL",
                  "TFEB", "MCOLN1", "GBA", "NPC1", "PSAP")

## Endosomal maturation and membrane trafficking.
##
## Added AFTER the genome-wide co-dependency scan, which showed that PIK3C3's
## nearest dependency neighbours are ESCRT, HOPS and RAB7 genes rather than
## either the ATG8 cascade or LYSOSOME_SET. The original LYSOSOME_SET captured
## lysosomal CONTENTS and acidification (cathepsins, v-ATPase, LAMPs); this set
## captures the TRAFFICKING MACHINERY that delivers cargo to lysosomes. They
## are different biology and the distinction turned out to matter.
##
## PIK3R4 and BECN1 are deliberately EXCLUDED: they are obligate PIK3C3 complex
## partners, so including them would guarantee a high correlation for trivial
## reasons and make the module comparison circular.
TRAFFICKING_MODULE <- c(
  ## ESCRT-0/I/II/III and VPS4
  "HGS", "STAM", "TSG101", "VPS28", "VPS37A", "SNF8", "VPS25", "VPS36",
  "CHMP2A", "CHMP4B", "VPS4A",
  ## HOPS / CORVET tethering
  "VPS11", "VPS16", "VPS18", "VPS33A", "VPS39", "VPS41",
  ## endosomal Rabs, tethers and maturation factors
  "RAB5A", "RAB7A", "EEA1", "UVRAG", "EPG5",
  ## retromer / retriever / CCC, endosomal recycling
  "VPS35", "VPS29", "VPS26A", "COMMD3", "COMMD9"
)

## --- DDR dependency targets ------------------------------------------------
## The "partner" side of the co-dependency test: genes whose CRISPR dependency
## reports a cell's reliance on specific DDR routes. These are the druggable or
## near-druggable nodes a combination trial would actually use.
DDR_TARGETS <- c(
  ## PARP / base excision
  "PARP1", "PARP2", "PARG", "XRCC1", "POLB", "FEN1", "LIG3",
  ## replication stress checkpoint (the ATR-CHK1-WEE1 axis)
  "ATR", "CHEK1", "WEE1", "PKMYT1", "TOPBP1", "CLSPN", "RPA1", "RPA2",
  ## homologous recombination
  "BRCA1", "BRCA2", "PALB2", "RAD51", "RAD51C", "BARD1", "RAD54L", "BRIP1",
  ## non-homologous end joining  (the PNAS synthetic-lethal axis)
  "PRKDC", "LIG4", "XRCC4", "XRCC5", "XRCC6", "NHEJ1", "POLL",
  ## microhomology-mediated / backup end joining
  "POLQ", "PARP1",
  ## ATM axis and Fanconi
  "ATM", "CHEK2", "MDC1", "NBN", "MRE11", "RAD50", "FANCD2", "FANCI",
  "USP1", "RAD18", "UBE2T",
  ## damage tolerance / outcome modifiers
  "SLFN11", "TP53", "CDKN1A", "TP53BP1", "SHLD2"
)
DDR_TARGETS <- unique(DDR_TARGETS)

## --- HRR genes scored for damaging mutations -------------------------------
## Used to build the HRD-proxy genotype feature. Note this is a PROXY: DepMap
## does not ship a mutational-signature HRD score, so we are approximating a
## phenotype with a genotype. State this as a limitation in the paper.
HRR_GENES <- c(
  "BRCA1", "BRCA2", "PALB2", "BARD1", "BRIP1",
  "RAD51B", "RAD51C", "RAD51D", "RAD50", "RAD54L",
  "XRCC2", "XRCC3", "NBN", "MRE11",
  "ATM", "ATR", "CHEK1", "CHEK2",
  "FANCA", "FANCB", "FANCC", "FANCD2", "FANCE", "FANCF", "FANCG",
  "FANCI", "FANCL", "FANCM",
  "PPP2R2A", "ABRAXAS1", "GEN1"
)

## --- Replication-stress transcriptional signature --------------------------
## E2F-driven replication/licensing and nucleotide-supply genes. High
## expression marks cells running replication hot -- the context where
## autophagy-supported nucleotide pools should matter most (per the
## replication-stress literature).
REPSTRESS_SIG <- c(
  "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7",
  "PCNA", "CDC6", "CDT1", "CDC45", "GINS1", "GINS2", "GINS3", "GINS4",
  "RRM1", "RRM2", "TYMS", "TK1", "DHFR",
  "E2F1", "CCNE1", "CCNE2", "CDK2", "CLSPN", "POLA1", "PRIM1"
)

## --- Single-gene genotype features ----------------------------------------
GENOTYPE_SNVS <- c("TP53", "ATM", "ATR", "BRCA1", "BRCA2", "ARID1A",
                   "KRAS", "STK11", "KEAP1", "SMARCA4", "RB1")
GENOTYPE_AMPS <- c("CCNE1", "MYC", "MYCN", "ERBB2")
GENOTYPE_DELS <- c("CDKN2A", "PTEN", "ATM", "RB1", "SMARCA4")

## --- Negative-control sets -------------------------------------------------
## Housekeeping-ish, biologically unrelated to DDR or autophagy. If these score
## as strongly as the autophagy set, the pipeline is reporting an artefact.
NEGATIVE_CONTROL_SET <- c("RPL13A", "ACTB", "GAPDH", "TUBB", "PPIA", "B2M")

cat("[genesets] primary autophagy set:", length(AUTOPHAGY_PRIMARY), "genes\n")
cat("[genesets] all autophagy modules:", length(AUTOPHAGY_ALL), "genes\n")
cat("[genesets] DDR dependency targets:", length(DDR_TARGETS), "genes\n")
