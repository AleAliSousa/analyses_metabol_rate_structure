# Borghammer 2012 Table 2: convergence feasibility

Checked 12 November 2026 against the local article, rendered Methods/Fig. 2
(p. 305) and Table 2 (p. 307), Table 1, and the existing Table S1 extraction.
Primary article: https://doi.org/10.1111/j.1600-0404.2011.01556.x
Heiss is a coauthor. This is a useful *separate relative-metabolism check*, most
promising for fine-grained s1b composition. It does not supply whole-thalamus or
cerebellar-cortex rates for the existing broad-region comparisons.

## Measurement and population

The controls comprise 11 adults (8 Aarhus, 3 Cologne; age 60.1 ± 7.5 in Table 1).
Use only Controls for normative comparisons; the other 21 participants have PD.
The HRRT scans have approximately 2-mm resolution. Histology-derived basal-ganglia
and thalamic atlas VOIs were warped to individual MRIs, visually checked and
occasionally manually corrected; intensities were extracted in native PET space.
The 10-mm smoothing described later is for the voxel-based analysis, not these
VOI means. Some PD scans were mirrored to align symptoms; do not infer normal
hemispheric biology from their sides.

Table 2 contains ratios under three normalizations: predefined white-matter
reference VOIs, global mean (GM means **global mean**, not grey matter), and a
reference cluster selected using the PD-versus-control comparison. Treat them as
three sensitivity versions of the same subjects, never three replications.
The bilateral control means of all eight subcortical nuclei/structures have
identical rank order across the three normalizations (Spearman rho = 1); this
is descriptive normalization stability, not three independent confirmations.
Prefer WM as the initial predefined-reference check; report GM and reference
cluster results as sensitivities. This is an analysis-design recommendation,
not a claim that WM ratios are absolute metabolism or statistically independent.

Absolute measurements differed significantly between centers, so the authors
abandoned pooled absolute analyses. Table S1 gives global/reference measurements
by center, not center-specific regional absolute rates. It cannot recover them by
multiplying a group-mean ratio by a group-mean denominator: covariance and
individual values are missing. Do not splice Table 2 ratios into Heiss units or
use them for absolute energy budgets.

## Anatomical boundaries

- Putamen and caudate: named structures; retain caudate-body histology explicitly
  as a sample, without claiming complete-mask identity.
- GPe and GPi: external and internal pallidum separately. Do not average them to
  create a whole-pallidum rate without the appropriate tissue weights.
- VA, VL and CM: ventral anterior, ventral lateral and centromedian thalamic
  nuclei. None is medial/mediodorsal thalamus, and their mean is not whole thalamus.
- STN: subthalamic nucleus, outside the thalamus aggregate.
- Pons: selected normalization VOI described among WM structures; not a measured
  complete-pons compartment. Cerebellum: **central cerebellar white matter**, as
  explicitly shown in Fig. 2D and Methods. Neither matches cortical tissue or
  a gross whole-structure volume merely because the table uses a short label.

## Actual overlap with the current data

| Comparison | Defensible candidate regions | Regional n |
|---|---|---:|
| s1b rebuilt at native ROI level | Putamen, caudate body, GPe, GPi, VA, pure CM | 6 |
| s1a stereology | Caudate, putamen | 2 |
| Johansen SV2A | Caudate, putamen | 2 |
| Holmes SV2A | Caudate, putamen, STN | 3 |
| Sjölin proteins | Caudate | 1 |
| s2 ERA / BEIP broad volumes | Caudate, putamen | 2 |
| Current s5a strict volumes | STN | 1 |

The s1b candidates were checked against `roi`, `dissection` and `donor_id` in
both extracted neuronal and nonneuronal metadata. Putamen, caudate, GPe and GPi
have three donors. Pure CM and VA occur in only H18.30.002 (one donor each).
Thus a six-region analysis has uneven donor support; a complete six-region
within-donor check is possible in only that one donor. The four-region,
three-donor core should be reported separately. These n values are anatomical
coverage, not guaranteed complete-case n for every cell measure.

Exclude `Human CM-Pf`: its dissection combines centromedian and parafascicular
nuclei. Do not automatically equate `Human VLN` to Borghammer VL: its source
annotation says ventral group of lateral nucleus, with uncertain extent. Treat
`Human STH` as a separate sensitivity at most: its source explicitly says
subthalamic nucleus **and nearby** (two donors). This newly inspected dissection
wording also warrants revisiting the existing Heiss STH assignment; the present
feasibility check does not silently change previously fitted primary analyses.
Do not confuse thalamic `Human CM` with amygdalar `Human CMN`.

The existing s1b parent summaries pool GPe/GPi and omit VA/CM because Heiss lacks
those nuclear rates. A Borghammer check must therefore rebuild composition from
the native ROIs rather than joining rates to those parent summaries. Most other
studies have too few matches for meaningful regional regression under the
project's four-region minimum; retain these as descriptive comparisons only.

## Recommended separate analysis

1. Use Controls only, with an equal-weight mean of published L/R means to obtain
   one point per bilateral structure. Label this estimator; do not call it a
   volume-weighted mean. Check left/right sensitivity, particularly STN and VA.
2. Rebuild donor-level composition for the six candidate native ROIs. Keep the
   four-region/three-donor core, six-region set, and one complete-donor analysis
   separate. Do not use cells as independent donor replicates.
3. Use WM normalization initially, then repeat with global mean and reference
   cluster on exactly the same regions. Report slopes, Spearman ranks and
   leave-one-region-out results; avoid highly parameterized models at n=6.
4. Test convergence of the composition hypothesis. Anatomical coverage differs
   substantially from the Heiss analysis; numerical slopes have different units.
   Agreement is supporting evidence, not an interchangeable absolute-rate dataset.
5. Check subject reuse before claiming fully independent cohorts. The paper
   explicitly reuses some subjects from Eggers et al. 2009; overlap with Heiss
   2004 is not established here.

No primary fits, input manifests or sibling exports were changed. The standalone
`checks/assess_borghammer_table2.R` reads the existing source and writes
`tables/anatomy/borghammer_controls_summary.csv` and
`tables/anatomy/borghammer_siletti_feasibility.csv`. No bilateral dispersion is
manufactured: Table 2 does not specify its printed dispersion type, and L/R
covariance is unavailable. Existing source extraction still has pending human
visual sign-off; this audit is an AI check of the available primary material.
