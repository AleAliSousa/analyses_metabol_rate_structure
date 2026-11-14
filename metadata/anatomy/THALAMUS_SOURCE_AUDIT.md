> Updated policy: [PRIMARY_MATCH_POLICY.md](PRIMARY_MATCH_POLICY.md). The initial thalamus-only counts below are retained as history, not current primary results.

# Thalamic anatomy: source audit, 12 November 2026

The questions differ: s1 compares tissue composition with metabolic rate; s2
compares regional volume differences with metabolic rate. Do not force s2 onto
the intersection of the s1 region sets. Match each outcome to an appropriate
metabolic measurement, and use shared-region checks within a research question.

Primary articles in Evo-M1-Trait-Data were checked, including Methods, tables,
selected rendered pages, and Sjölin's supplementary Methods. Paper-specific
`THALAMIC_ANATOMY.md` notes were added to these twelve upstream folders:

| Folder | What the source actually measures | Remaining boundary uncertainty |
|---|---|---|
| Heiss_etal_2004 | Medial-thalamic, lateral-geniculate and medial-geniculate VOIs, separately | No complete masks or explicit equivalence of the medial label to a histological MDT delineation |
| Karlsen_etal_2014 | Histologically delineated mediodorsal nucleus, one hemisphere | Its exact overlap with Heiss's spherical medial VOI is unverified |
| Wegiel_etal_2014 | Structure-level thalamic stereology, one hemisphere | Detailed boundaries referred to an earlier report; constituent nuclei not fully established here |
| Johansen_etal_2024 | FreeSurfer 7.2 aggregate thalamic SV2A ROI | Exact nucleus inclusion not enumerated |
| Holmes_etal_2024 | AAL aggregate thalamic BPND ROI; not in the PVC-corrected list | Exact equivalence to FreeSurfer/Hammers not established |
| Varnäs_etal_2020 | Thalamic autoradiography ROI on sampled sections | Table 1's nuclear weighting/boundaries unspecified; figures show nuclear heterogeneity |
| Andersen_etal_2022 | Bilateral Hammers-atlas thalamic ROI for both FDG and UCB-J | FDG is relative SUVR, not absolute rCMRGlc |
| Sjölin_etal_2022 | Small tissue sample labelled thalamus | Nucleus, boundaries, coordinates and hemisphere unspecified in article and supplement |
| Mackes_etal_2020 | Aggregate thalamus volume, FreeSurfer 6.0.0 | No nuclear breakdown; exact boundaries may differ across software versions |
| Walhovd_etal_2011 | Aggregate thalamus volume, FreeSurfer 4.0.1 | Same limitation; do not assume identical segmentation to Mackes |
| McLaughlin_etal_2014 | Aggregate thalamus volume, FreeSurfer 5.0 | Units and side basis are not printed in Table 4 |
| Sheridan_etal_2022 | Left/right aggregate thalamus volumes, FreeSurfer 5.3 follow-up | Bilateral sums are derived; L/R covariance unavailable |

The upstream root used was
`~/Library/CloudStorage/OneDrive-AllenInstitute/Species/Evo-M1-Trait-Data`.
Each note gives exact local source filenames and page/section references.
Only notes were written upstream; source values, registry and snapshots were not changed.
No human visual sign-off or author contact is implied.

## Historical first-pass consequences (superseded)

- Heiss's three rates are 36.6, 23.1 and 20.2 µmol/100 g/min. Their unweighted
  mean, 26.6333 (rounded to 26.6 by the reader), is an analysis-derived proxy.
  It is not a published whole-thalamus rate. It gives geniculate territories
  two-thirds of the weight and does not sample all nuclei.
- Sjölin must be labelled **thalamic tissue, subregion unspecified**, not
  whole-thalamus tissue. The supplement describes cutting roughly 450–550 mg
  from stored samples but does not locate the thalamic dissection. Fig. 1 is
  a concentration schematic, not evidence of a complete-thalamus homogenate.
- Current Heiss-2004 metabolic fits conservatively omit the unmatched broad
  thalamus point, with the three-territory proxy retained in named sensitivity
  fits. This is a provisional anatomical matching choice, not a requirement
  that s1 and s2 have identical regions. The source thalamus volumes remain in
  the data and s2 descriptive volume plots.
- Primary region counts under that choice: Johansen 13, proteins 10, ERA 7,
  Holmes 12, BEIP age 8–10 8, BEIP age 16 7 per contrast. Including the proxy
  adds one region to each. Johansen already excluded thalamus before this session.
- Pairwise composition comparisons currently include Johansen–proteins (8
  shared regions) and Johansen–Holmes (9); the separate s2 comparisons use
  ERA–BEIP shared regions. Cross-question s1–s2 intersection fits were removed
  following the user's clarification. These are descriptive association checks,
  not formal tests that slopes in different units are equal.
- **Additional source correction:** Johansen's Methods explicitly says
  cerebellum cortex. Its Cerebellum row is a cerebellar-cortex match, not the
  whole-cerebellum proxy described in the old map. The map has been corrected.
- s1b already keeps lateral-geniculate, medial-geniculate and medial-thalamic
  groups separately. Its current map combines `Human MD` and `Human MD-Re`:
  MD+Re contains reuniens as well as mediodorsal tissue. That aggregation is
  not a proven exact match to Heiss's medial VOI. No Siletti paper folder was
  found in Evo-M1; this remains a project crosswalk review item. Dataset origin:
  [Human Brain Cell Atlas](https://data.humancellatlas.org/hca-bio-networks/nervous-system/atlases/brain-v1-0).

## Next anatomical work

1. For s2, investigate a whole-thalamus metabolic source with compatible
   boundaries; keep all its regions on the same measurement scale. Do not
   splice Andersen SUVR into Heiss absolute rates.
2. Andersen's paired FDG/SV2A data offer a separate matched-ROI replication
   within the tissue-composition question, including thalamus.
3. Seek Sjölin dissection records before a nucleus-specific match. Determine
   Heiss's medial VOI boundaries before treating Karlsen MDT or Siletti MD/MD+Re
   as exact equivalents. Varnäs figures show relevant nuclei but do not provide
   a ready numeric human MDT row.

Online primary-source corroboration: [Sjölin](https://pmc.ncbi.nlm.nih.gov/articles/PMC9241296/),
[Johansen](https://pmc.ncbi.nlm.nih.gov/articles/PMC11326867/),
[Sheridan](https://pmc.ncbi.nlm.nih.gov/articles/PMC9544316/).
