# Study 3 — Phase 1 missing-data strategy (decision record)

_Last updated 2026-06-04._

## The problem

Study 3 fits, per brain region, `log(region) ~ log(rest-of-brain)` by PGLS, predicts
the human value, and compares regions' human prediction errors against rCMRGlc.
Each region is measured in a **different subset of species**, and those subsets
differ in **primate-grade composition**. Phylogenetic GLS accounts for shared
ancestry but does **not** fully remove a sampling confound: if the insula slope
is estimated from an ape-heavy/prosimian-light sample and the cerebellum slope
from the opposite, the two slopes (and their human residuals) are not strictly
comparable.

Grade composition per region (from `checks/s3/phase1/s3_region_missingness_clade.csv`):

| Region | Strepsirrhine | Tarsier | NW monkey | OW monkey | Ape+human | N |
|---|---|---|---|---|---|---|
| Cerebellum / LGN / Striatum / Hippocampus | 17 | 1 | 12 | 8–10 | 5–7 | 45–49 |
| Insula (grey) | 6 | 0 | 8 | 7 | **9** | 30 |
| Amygdala | 8 | 1 | 5 | 2 | 2 | 18 |
| Neocortex grey (Frahm) | 9 | 1 | 7 | 3 | 4 | 24 |

The concern is real and is exactly the insula-vs-cerebellum contrast.

## Why strict listwise deletion is not viable

Restricting to species complete across **all 13** modelled regions leaves **4
species** (Cebus albifrons, Galagoides demidoff, Homo sapiens, Piliocolobus
badius). The "common-N cliff" (`s3_common_N_cliff.csv`) shows where the sample
collapses as regions are added in order of decreasing coverage:

```
LGN (49) → +Striatum/Mesenc/Cerebellum/Hippocampus → 45
        → +Area striata (V1)                       → 42
        → +Insula                                  → 18   <- first big cliff
        → +Pallidum/SubthN/Capsula                 → 14
        → +Neocortex white/grey                    → 7
        → +Amygdala                                → 4
```

There is, however, a **natural 6-region core (LGN, striatum, mesencephalon,
cerebellum, hippocampus, V1) with 42 species in common** — a genuinely usable
balanced set if a fully comparable analysis is wanted.

## Decision

**Primary analysis: max-data per-region fitting (status quo), with the confound
made transparent.** This keeps statistical power and is what the main script does
by default (`CONFIG$restrict_species = NULL`).

**Supported alternatives, all wired in:**

1. **Diagnostic (always report).** `s3_0_missingness_clade_diagnostic_04062026.R`
   produces the grade-composition table/figure and the common-N cliff. Include in
   the supplement so reviewers see the sampling structure.

2. **Balanced common-set option.** Pass a species vector to
   `CONFIG$restrict_species` (e.g. the 42-species 6-region core) to refit all
   regions on identical taxa. Fully comparable; lower power.

3. **Phylogenetic multiple imputation (sensitivity).**
   `s3_1_phylo_multiple_imputation_04062026.R` puts every region on the same
   species set by imputing missing values under a multivariate Brownian model on
   the tree (`Rphylopars`), generates _m_ imputations, refits the BM PGLS on each,
   and pools by **Rubin's rules** — preserving the between-imputation variance.
   Use this to show the cross-region pattern is not an artefact of differential
   sampling.

### On imputation (the points that matter)

- **Use phylogenetic imputation, not plain `mice`.** Standard MI assumes
  independent rows; these species are not independent. `Rphylopars` borrows
  strength from phylogeny + cross-trait correlations. (`mice` + phylogenetic
  eigenvectors is a turnkey alternative if Rubin pooling out-of-the-box is
  preferred.)
- **Keeping all imputations and pooling is correct** — that is what preserves the
  variance you (rightly) did not want to throw away.
- **Homo sapiens is observed for all 13 regions, so it is never imputed.** We only
  impute non-human species to balance the comparison sample; the human value being
  tested is untouched. This makes the imputation far safer than typical use.
- **Imputation cannot create information.** Amygdala (~70% missing) and the Frahm
  neocortex pair (~60%) are heavily imputed; their imputed slopes are
  model-driven. `frac_imputed` is reported per region so heavily-imputed cells can
  be flagged or discounted. For these regions, transparent reporting may be more
  defensible than imputation.

## Phase 1.2 — neocortex grey de-overlap

`Neocortex grey (remaining) = NeoG_Frahm − V1(ASG_Sousa) − Insula(grey)` (grey
matter only; white matter handled separately).

Key data fact: **all 24 neocortex-grey species have V1, but only 9 also have
insula**, and insula is only **~0.5–1.1%** of neocortex grey.

- Subtracting **V1** is "free" — keeps n = 24, and V1 is the material term (3.7%
  for humans, up to ~17% in small species).
- Subtracting **insula** is what drops neocortex to n = 9, for a ~1% correction.

`CONFIG$neocortex_subtract_insula` controls this:
- `TRUE` (default) — full definition, n = 9.
- `FALSE` — V1-only remainder, n = 24 (insula ~1% left in).

**Recommended:** keep the full definition but **impute the ~1% insula term** to
recover n = 24 — a clean, low-stakes imputation (contrast with imputing amygdala).
This is the natural bridge between 1.1 and 1.2.

## Files

| File | Purpose |
|---|---|
| `scripts/s3_predicValuesPGLS_04062026.R` | Tidied main analysis + Phase 1 toggles |
| `scripts/s3_0_missingness_clade_diagnostic_04062026.R` | Phase 1.1 diagnostic (tables + figures) |
| `scripts/s3_1_phylo_multiple_imputation_04062026.R` | Phase 1.1 MI sensitivity (Rphylopars + Rubin) |
| `checks/s3/phase1/s3_region_missingness_clade.csv` | Per-region N + grade composition |
| `checks/s3/phase1/s3_common_N_cliff.csv` | Common-N as regions are added |

The previous main script (`s3_predicValuesPGLS_01062026.R`) is left untouched for
diffing; set the de-overlap/restrict toggles off to reproduce its numbers.
