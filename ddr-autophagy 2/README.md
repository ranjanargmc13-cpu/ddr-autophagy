# DDR × autophagy dependency analysis

Analysis code for a study of a metformin–chloroquine–olaparib combination in
triple-negative breast cancer, and of whether the autophagy-based rationale for
such combinations is supported by cancer cell line dependency data.

**Headline findings.** Across 1,361 DepMap cell lines: the ATG8 conjugation
cascade is broadly dispensable while PIK3C3/VPS34 is near-common-essential
(77.9% of lines below −0.5); VPS34 is not a member of the ATG8 dependency
module; and autophagy dependency shows no co-dependency with any of 47 DNA
repair genes.

---

## Requirements

```r
install.packages("data.table")
```

R ≥ 4.1. `data.table` is the only hard dependency. Tested on R 4.6.1 (macOS).

## Data

Download from <https://depmap.org/portal/data_page/?tab=allData> into a single
folder (referred to below as `~/depmap`):

| File | Needed by |
|---|---|
| `CRISPRGeneEffect.csv` | everything (required) |
| `Model.csv` | everything (required) |
| `OmicsSomaticMutationsMatrixDamaging.csv` | `run_all.R`, `biomarker_scan.R` |
| `OmicsCNGene.csv` *(or `OmicsCNGeneWGS.csv`, symlinked to that name)* | `run_all.R` |
| `OmicsExpressionProteinCodingGenesTPMLogp1.csv` *or* `OmicsExpressionRawReadCountHumanProteinCodingGenesStranded.csv` | `biomarker_scan.R` |

File names are declared once, in `DEPMAP_FILES` at the top of
`R/02_load_depmap.R`. Change them there if your release differs.

The DepMap release used for the reported results was **[FILL IN RELEASE]**.
State it in any publication — gene effect values shift between quarters.

## Reproducing the analysis

**Always run the self-test first.** It uses synthetic data with a known planted
effect and a planted confounder, and verifies the pipeline recovers one and
controls for the other:

```bash
Rscript run_all.R --selftest --null_B=1000
```

Three criteria must pass: the five planted targets (ATR, CHEK1, WEE1, POLQ,
PARP1) appear in the top five with negative rho; the negative-control gene set
returns ~0 hits; and the planted targets survive the empirical null. Do not run
on real data until they do.

Then:

```bash
# DDR × autophagy co-dependency, lineage ranking, empirical null
Rscript run_all.R --data_dir=~/depmap --no_expression --null_B=1000

# dependency profile of specific cell lines
Rscript check_lines.R ~/depmap

# module membership of VPS34 vs the ATG8 cascade
Rscript codependency_module.R ~/depmap

# memory-light variant of the above (~15 MB instead of ~1 GB)
Rscript codependency_lite.R ~/depmap

# what context creates VPS34 dependency
Rscript biomarker_scan.R ~/depmap
```

### `--null_B` matters

The smallest achievable empirical p-value is `1/(B+1)`. With `B = 50` and 47
targets, no result can reach FDR < 0.05 however strong the effect — the first
run of this analysis produced exactly that artefact. Use `B ≥ 1000`. The code
now warns when `B` is too small for the test to succeed.

---

## Scripts

| Script | Purpose | Status |
|---|---|---|
| `R/01_genesets.R` | Gene set definitions with provenance | ✅ used |
| `R/02_load_depmap.R` | Schema-validating loaders | ✅ used |
| `R/03_scores.R` | Dependency scores + coherence diagnostics | ✅ used |
| `R/04_stats.R` | Co-dependency and genotype association tests | ✅ used |
| `R/05_lineage.R` | Per-lineage aggregation | ✅ used |
| `run_all.R` | Main pipeline | ✅ used |
| `tests/make_fixtures.R` | Synthetic data with planted effect | ✅ used |
| `check_lines.R` | Per-cell-line dependency lookup | ✅ used |
| `codependency_module.R` | Module membership, bootstrap, genome-wide scan | ✅ used |
| `codependency_lite.R` | Low-memory variant | ✅ used |
| `biomarker_scan.R` | Lineage, mutation and expression correlates | ✅ used |
| `selftest_interactive.R` | Self-test from an interactive R session | ✅ used |
| `prism_drugs.R` | Drug sensitivity vs dependency | ⚠️ **never executed** |

`prism_drugs.R` was written but never run, because the PRISM drug-response
matrix was not available. It is included for completeness and should be treated
as untested code.

---

## Methodological notes

These are design decisions that materially affect the results. They are
documented in the source as well.

**Sign conventions.** Chronos gene effect: 0 = no effect, negative = dependent.
`build_dependency_score(orient="vulnerability")` returns a sign-flipped score
where higher = more dependent. In `test_codependency` the autophagy score is
vulnerability-oriented while the DDR gene is raw Chronos, so **negative rho is
hypothesis-consistent**. In `test_ddr_features`, **positive beta** is. In
`rank_lineages` both sides are flipped, so **positive `coupling_rho`** is.
Getting these backwards inverts the conclusion.

**Covariates are not optional.** All correlations are partial, adjusted for
lineage and for each line's mean gene effect across all genes. Tissue of origin
dominates cancer cell line biology, and poor-quality screens make every gene
look essential at once — either will manufacture correlation between arbitrary
gene pairs. The self-test plants a screen-quality confounder specifically to
verify the adjustment works.

**Effect sizes are calibrated, not compared to zero.** PIK3C3 and PIK3R4 are
obligate complex partners and correlate at r = 0.331; that is what "same
module" looks like in this data. A correlation of 0.08 against that ceiling is
a quarter of a module relationship, regardless of its p-value. With n > 1,300,
trivial correlations are highly significant.

**Parametric p-values from composite scores are anti-conservative.** The genes
making up a composite are correlated, so effective degrees of freedom are below
n. Empirical p-values from variance-matched random gene sets are reported
alongside and should be preferred. In the real-data run, a housekeeping
negative-control set produced 3/47 nominal hits while the autophagy set
produced 0/47 — direct evidence the parametric test was not conservative.

**Gene set choice drove an early false negative.** The original `LYSOSOME_SET`
captured lysosomal contents and acidification (cathepsins, v-ATPase, LAMPs) and
found nothing. The genome-wide scan, which imposes no gene set, showed VPS34's
nearest neighbours are ESCRT, HOPS and RAB7 — trafficking machinery, absent
from that set. `TRAFFICKING_MODULE` was added afterwards. The unsupervised
analysis is therefore the primary evidence and the module comparison secondary;
refining a gene set until it produces a result is not a safe procedure.

**Common essentials compress correlations.** A gene essential in most lines
varies little, which mechanically pulls its correlations toward zero. This
applies to VPS34 itself and limits the sensitivity of the module and expression
analyses.

---

## Known issues

- `test_ddr_features` (TEST B in `run_all.R`) returned `n_pos = 0` for all HRR
  mutation features on the real data — the mutation matrix did not align with
  the dependency table. **Not fixed.** The lineage ranking in `run_all.R`
  depends on those features and is therefore unreliable on real data; it was
  excluded from the reported results.
- Two earlier versions of `codependency_module.R` printed a fixed verdict line
  regardless of the computed value. Fixed — the interpretation is now
  conditional and calibrated.
- `p_two_sided` in `empirical_null_codependency` counts opposite-sign bootstrap
  replicates without doubling and is therefore one-sided. Prefer the confidence
  intervals, which are computed correctly.

---

## Limitations of this approach

CRISPR knockout is chronic and profiles untreated cells in nutrient-replete 2-D
culture. The therapeutic hypothesis concerns autophagy *induced* by drug
exposure, which baseline dependency data cannot address. The claim supported
here is the absence of *constitutive* coupling. Co-dependency also shows
co-variation, not shared pathway membership.

## License

MIT — see `LICENSE`.

## Citation

[Add once published.]
