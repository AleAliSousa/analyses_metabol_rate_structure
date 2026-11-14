# Primary anatomical matching — 12 November 2026

This policy supersedes the earlier thalamus-only exclusion. It applies before
examining regression results. Heiss 2004 remains the primary absolute-rate source.

## Inclusion rule

Use the same named anatomical structure and compatible tissue compartment.
Histological samples contained in that structure can be compared at the parent
structure level, but they are sampled-structure evidence, not identical whole-ROI
measurements. Atlas versions, sampled location and hemispheric coverage may differ.
Do not claim voxel-identical anatomy across PET, MRI, stereology and tissue assays.
A larger structure cannot inherit a rate from one constituent, or from an arbitrary
mean of selected constituents. Do not treat a parent and its children as independent
regional points. Missing anatomy stays missing; inclusion does not depend on p values.

Every study crosswalk now has `primary_included` and `primary_match_note`.
`derive_study_heiss_rates()` retains rejected source rows but returns NA for their
rate and SD. Explicit `policy = "legacy_proxy"` is available for labelled historical
audits; it is never a fallback. s1b drops rejected ROI assignments before pooling
cells. s1a descriptive source data retain the cortex aggregate; regressions omit it.

## Decisions

- s1a: four cortical lobes, three striatal structures and amygdala remain (8 region
  means for neuron/glia density; only 5 have detailed glial subtype/count data).
  KP neocortex excludes limbic territories. Hanson samples dorsal/medial caudate,
  associative putamen and accumbens; the caudate density is an unweighted mean of
  two samples, not a whole-nucleus estimate. Morgan delineates the amygdaloid
  complex and five subdivisions. These sampling limitations remain explicit.
- s1b: 23 parent groups become 19. Exclude SEP/SI -> unspecified basal forebrain,
  CBL -> unverified cortex-only compartment, CbDN -> dentate, and MD/MD-Re ->
  unspecified medial-thalamic VOI. Exclude BNST from amygdala pooling; other
  amygdala ROIs remain. LG and MG are retained separately. Cortical samples are
  grouped by anatomical lobe; no inference of whole-lobe cell totals is made.
- Johansen: 12 primary regions; exclude aggregate white matter and thalamus.
  Cerebellum is retained because Methods explicitly specifies cerebellum cortex.
- Holmes: 11 primary regions; omit thalamus and conservatively omit cerebellum
  until its aggregate mask is demonstrated to match the Heiss cortical compartment.
  Mixed regional partial-volume correction still requires a separate sensitivity.
- Sjölin: 7 primary regions (four cortical lobes, caudate, hippocampus, internal
  capsule). Exclude pooled lobar white matter, unspecified thalamic tissue,
  unspecified midbrain tissue and cerebellum. Do not call its small samples
  whole-structure homogenates. Five-predictor models at n=7 are exploratory and
  nearly saturated; the regional univariate/LOO results are more informative.
- s2 ERA: 7 primary regions; whole-thalamus volume stays in descriptive tables
  but lacks a compatible Heiss-2004 rate. BEIP: 7 childhood, 6 adolescent regions
  per contrast after also excluding whole cerebellum/whole white matter. Ages
  share one cohort; CAUG/FCG share controls. Do not force s2 to the s1 intersection.
- s4a: zero primary compartment-matched metabolic estimates from endocasts.
  Grey+white parcels and the 50:50 sensorimotor split remain explicitly modelled
  proxies under `s4a_proxy_sensitivity` output paths and script names. They are
  not evidence from a precise grey-matter metabolic match.
- s5a: 8 primary regional volume/rate pairs; omit parent splits, remainder
  cortex, V1 -> occipital lobe, nuclear means -> thalamus/midbrain, and mixed
  cerebellar/white-matter proxies. Report coverage as a partial regional sum (currently only 3.23% of brain volume).
- s5b: use a published whole-brain glucose rate (Clarke & Sokoloff via Karbowski)
  for whole-brain metabolism. Never scale the sparse s5a regional mean to the
  entire brain. s5c and exported s2 regressions inherit the strict s2 map.
- s3 is archived; its dependent network-residual script is excluded from the active runner pending its new anatomy crosswalk. The primary s4b arterial and s5d species models use whole-brain
  measurements. s4b also carries an s4a cross-check; its derived columns and
  overlap table now explicitly say `proxy_sensitivity`, with corrected input paths.

Shared-region comparisons remain within questions. A common parent label permits
an exploratory comparison but does not prove identical samples or atlas boundaries.
Source values are preserved; the audit is in `tables/anatomy/`.

## Metabolic-source convergence

1. **Heiss 1984, not 1991.** Table 1 of Heiss 2004 cites reference 21:
   Heiss et al., JCBFM 4:212–223 (1984), DOI 10.1038/jcbfm.1984.30.
   The 2004 discussion says the older rates account for different lumped constants.
   Our `heiss_1984_mean_umol_100g_min` is this *2004 re-tabulation*, not a fresh
   extraction of native 1984 Table 3. `s1_heiss1984_sensitivity.R` keeps the same
   complete-case regions for both rate sources. Broad old thalamus and mixed
   midbrain do not become medial thalamus, substantia nigra or red nucleus simply
   because the 2004 table places values beside those rows; those are excluded.
2. **A fuller 1984 version is a useful next extraction.** The primary article
   reports 43 side-specific sites, summarized as 22 bilateral/midline regions,
   including cingulate, pons, mixed midbrain, lobar white matter and broad thalamus.
   Seven awake male volunteers; approximately 8-mm PET resolution and 11-mm
   sections. Native kinetic and autoradiographic estimates must be kept distinct.
   Broad sampled thalamus is closer to MRI thalamus than three selected nuclei,
   but complete-mask equivalence is not established. Check sampled volume and
   boundaries before including it in a strict s2 analysis. Lentiform nucleus is
   not putamen alone. No 1984-only regions are spliced into the 2004 regression.
3. **Andersen 2022 is now a separate paired-tracer replication.** Table S1 supplies
   13 regions in 11 young and 15 older healthy adults, with FDG and UCB-J measured
   in the same atlas ROIs and subjects. This allows thalamus, striatum and cingulate
   without borrowing a Heiss rate. Fit age groups separately; report OLS, ranks,
   LOO and omission of influential substantia nigra. GM-normalized uptake is not
   absolute metabolic rate or Bmax. Shared normalization, rounding and regional
   dependence limit inference from group means.
4. **Heiss 1985 is a state-sensitivity candidate, lower priority.** Four evaluated
   volunteers (three nondreamers plus one dreamer), after about 28 h on duty with
   little/no sleep; waking scans are not an ordinary rested normative cohort.
   Keep waking/sleep states separate and do not count repeated people or assume
   independence from the 1984 laboratory cohort without participant records.

Source references: [Heiss 1984 primary article](https://journals.sagepub.com/doi/pdf/10.1038/jcbfm.1984.30),
[Heiss 2004](https://pubmed.ncbi.nlm.nih.gov/15534048/),
[Heiss 1985](https://pubmed.ncbi.nlm.nih.gov/3872693/).
Local PDFs and supplements in Evo-M1 supplied the remaining methods checks.
