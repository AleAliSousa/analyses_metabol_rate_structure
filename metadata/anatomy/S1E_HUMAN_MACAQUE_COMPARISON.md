# S1E — Human <-> macaque region-group comparison: receptor classes vs regional CMRglc (PROTOTYPE)

Created 2026-11-17. Script `scripts/s1e_3_human_macaque_groups.R` (tables) + `scripts/s1e_3_fig_human_macaque.py` (figures).
**Descriptive only.** Nothing here enters `run_all.R` / `output_manifest.csv`. Companion to
`human_macaque_region_group_homology.csv` / `HUMAN_MACAQUE_REGION_GROUPS.md` (name crosswalk; read that first) and to
`S1E_RECEPTOR_PROTOTYPE.md` (human-only tiers s1e_0-s1e_2).

## Outputs
- `tables/s1e/human_macaque_group_panel.csv` — 126 rows = 63 region groups x 2 species. Per row: receptor-class means
  (raw, in the species' own units) and within-species z-scores, CMRglc (raw, z, z over cortical groups only, source,
  n datasets), the extended-CMRglc sensitivity fill, dataset-specific CMRglc columns (`value_*`, `structures_*`),
  member areas / structures, and coverage flags (`has_receptors`, `has_cmrglc`, `shared_group_primary`,
  `shared_group_extended`, `shared_receptors_both_species`).
- `tables/s1e/human_macaque_group_coverage.csv` — one row per group: why it is / is not in the shared set.
- `tables/s1e/human_macaque_spearman.csv` — within-species Spearman (class vs CMRglc) per basis and cross-species
  concordance (human z vs macaque z) per variable; `groups` column lists the groups each row uses.
- `figs/s1e/human_macaque_panels.png/.pdf`, `figs/s1e/human_macaque_concordance.png/.pdf`.

## Methods
**Unit of analysis.** The region groups of `human_macaque_region_group_homology.csv`. Every dataset's region names were
taken verbatim from that table; no name was reassigned here. Each Zilles 2017 area (44) and each Froudist-Walsh 2023 area
(109) enters exactly one group (asserted in the script), as does every CMRglc structure of each metabolic table.

**Receptors.** Human: Zilles & Palomero-Gallagher 2017 Supplementary Table 2, `stratum == all_layers`, fmol/mg protein,
15 receptors. Macaque: Froudist-Walsh et al. 2023 Supplementary Table 1, fmol **per neuron**, 14 receptors, *Macaca
fascicularis*. Group value = unweighted mean over member areas of the per-area class sum. Classes as in s1e_1/s1e_2:
iGlu = AMPA+NMDA+kainate; GABA = GABAA+GABAA/BZ+GABAB; monoamine GPCR = alpha1+alpha2+5-HT1A+5-HT2+D1; E/I = iGlu/GABA.
**ACh is harmonised to M1+M2+M3 in both species** because FW2023 has no alpha4beta2; the human 4-receptor version is kept
as `class_ACh_incl_a4b2` (within-human sensitivity only). `total_shared14` sums the 14 receptors both tables share;
`total_15` (human only) adds alpha4beta2.

**CMRglc, human (umol/100 g/min, FDG PET).** Cortical groups: Heiss et al. 1985 Table I, awake non-dreamer arm (n = 3;
the same column choice as s1e_2, `cmrglc_heiss1985_awake_nd3`), mean over member structures where a group has two
(e.g. `inferior_middle_temporal`). Non-cortical groups and lobe-level cortical composites: Heiss et al. 2004 Table 1 via
`data_intermediate/heiss_2004_regions.csv`. `cmrglc_source` records which. Chugani et al. 1987 adult values (n = 7) are
carried as `value_chugani1987_adult` and used only in the extended sensitivity fill (below).

**CMRglc, macaque (umol/100 g/min, 2-DG autoradiography).** Primary = simple mean of the Kennedy et al. 1978 group value
(conscious *M. mulatta*, n = 7; 3 of 39 structures n = 6) and the Palombo et al. 1990 group value (CONTROL animals only,
`Controls (n = 4)`, right and left hemisphere arms averaged per structure; MPTP arms never used; Tables 1-5 treated as one
dataset), each being the mean over member structures; `cmrglc_n_datasets` = 1 or 2. Noda et al. 2002 (FDG PET, conscious
young rhesus, `young_average` arm) and Meguro et al. 1999 (FDG PET, **anaesthetised baboon**, pre-surgery baseline of both
surgery groups averaged) are carried as dataset-specific columns and used only in the extended fill.

**Extended CMRglc (sensitivity only).** Where the primary source lacks a group, the group is filled from a secondary
dataset after a multiplicative scale match k = mean(primary)/mean(secondary) over the groups both have (any tissue
class): k_chugani = 1.336 over 12 groups (human); k_noda = 1.364 over 3 groups, k_meguro = 1.094 over 13 groups
(macaque). A raw fill would mis-rank: Chugani adult values are ~25 % below Heiss for the same structures. The fill adds
`anterior_cingulate` (human side, Chugani) and `inferior_parietal` (macaque side, Meguro baboon). Noda fills no
receptor-bearing group (its ROIs are lobe-level).

**z-scores.** Within species, over all groups of that species that have the variable (receptor classes: the 19 human /
14 macaque groups with receptor areas; CMRglc: all groups with a primary value, cortex and subcortex together).
`cmrglc_z_cortex` repeats the CMRglc z over cortical groups only. Spearman rho is invariant to this choice; only the
scatter positions in the figures depend on it. **Raw receptor values are never compared across species** (fmol/mg protein
vs fmol/neuron); the cross-species rows use z-scores (equivalently ranks) only.

**Statistics.** Spearman rho with naive asymptotic p-values; NA when n < 5. No spatial-autocorrelation null, no model,
no correction for the shared phylogenetic/structural dependence of the groups. At n = 6-13 the p-values are not
informative and are printed only for completeness.

## Group coverage — why the shared set has 6 groups (not 8-9)
`HUMAN_MACAQUE_REGION_GROUPS.md` flagged 9 groups as `fully_covered_both_species`, counting *any* human metabolic dataset
(incl. Andersen 2022 SUVR, Chugani 1987) and *any* macaque metabolic dataset (incl. Meguro baboon). With the inputs this
task names for CMRglc (human: Heiss 1985/2004; macaque: Kennedy 1978 + Palombo 1990 controls) the primary shared set is
**6 groups**, all cortical: dorsolateral_prefrontal, orbitofrontal_ventral_prefrontal, premotor, primary_somatosensory_motor, primary_visual, visual_association.
The extended fill gives **8** (anterior_cingulate, inferior_parietal added); 7 when the baboon dataset is excluded.
`posterior_cingulate_retrosplenial` (the ninth group of the draft count) has no human CMRglc in any absolute-units source —
its only human metabolic entry is Andersen 2022 (CSO-normalised SUVR, not used here). All 13 groups with receptors in
both species, and the reason the 7 non-shared ones drop out:

| group (receptors in both species) | why not in the primary shared set |
|---|---|
| frontal_eye_field_area8 | no human Heiss CMRglc group |
| frontal_pole | no CMRglc in either species (Heiss / Kennedy+Palombo) |
| inferior_frontal_broca_homolog | no macaque Kennedy/Palombo CMRglc structure |
| anterior_cingulate | no human Heiss CMRglc group |
| posterior_cingulate_retrosplenial | no human Heiss CMRglc group |
| inferior_parietal | no macaque Kennedy/Palombo CMRglc structure |
| superior_parietal | no macaque Kennedy/Palombo CMRglc structure |

### Shared groups (extended set) — composition and CMRglc provenance
| group | tissue | homology | human: n Zilles areas | human CMRglc (umol/100 g/min), source | macaque: n FW areas | macaque CMRglc (umol/100 g/min), source | n macaque 2-DG datasets |
|---|---|---|---|---|---|---|---|
| anterior_cingulate | cortex | high | 2 | 39.8, Chugani1987_adult_n7 x 1.336 (scale-matched over 12 groups) | 9 | 40.0, Palombo1990_controls_n4 | 1 |
| dorsolateral_prefrontal | cortex | moderate | 2 | 41.4, Heiss1985_nondreamers_awake_n3 | 10 | 54.5, Palombo1990_controls_n4 | 1 |
| inferior_parietal | cortex | moderate | 4 | 38.2, Heiss1985_nondreamers_awake_n3 | 8 | 48.0, Meguro1999_BABOON_anaesthetised_PET x 1.094 (scale-matched over 13 groups) | 0 |
| orbitofrontal_ventral_prefrontal | cortex | moderate | 2 | 42.1, Heiss1985_nondreamers_awake_n3 | 10 | 52.0, Palombo1990_controls_n4 | 1 |
| premotor | cortex | high | 1 | 42.6, Heiss1985_nondreamers_awake_n3 | 13 | 45.0, Palombo1990_controls_n4 | 1 |
| primary_somatosensory_motor | cortex | high | 5 | 44.2, Heiss1985_nondreamers_awake_n3 | 9 | 45.2, Kennedy1978_conscious_n7+Palombo1990_controls_n4 | 2 |
| primary_visual | cortex | high | 1 | 40.1, Heiss1985_nondreamers_awake_n3 | 1 | 66.0, Kennedy1978_conscious_n7+Palombo1990_controls_n4 | 2 |
| visual_association | cortex | moderate | 6 | 36.1, Heiss1985_nondreamers_awake_n3 | 13 | 58.0, Palombo1990_controls_n4 | 1 |

### Within-species z-scores of the 6 primary shared groups
| group | human iGlu z | macaque iGlu z | human GABA z | macaque GABA z | human ACh (M1+M2+M3) z | macaque ACh (M1+M2+M3) z | human monoamine GPCR z | macaque monoamine GPCR z | human E/I index z | macaque E/I index z | human CMRglc z | macaque CMRglc z |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| dorsolateral_prefrontal | 0.54 | 0.49 | 0.32 | 0.63 | 0.23 | 0.73 | -0.34 | 0.65 | 0.17 | 0.06 | 1.07 | 0.56 |
| orbitofrontal_ventral_prefrontal | 0.79 | 1.27 | 0.73 | 1.09 | 0.49 | 1.14 | -0.04 | 0.82 | 0.01 | 1.00 | 1.16 | 0.42 |
| premotor | -2.84 | -0.26 | -2.28 | -0.13 | -2.50 | -0.31 | 3.13 | 0.46 | -1.02 | -0.45 | 1.23 | 0.03 |
| primary_somatosensory_motor | -0.63 | -1.29 | -0.47 | -1.35 | -1.09 | -1.36 | -0.59 | -1.22 | -0.28 | -0.51 | 1.43 | 0.04 |
| primary_visual | 0.82 | -2.16 | 2.17 | -2.20 | 0.48 | -2.10 | 0.10 | -2.12 | -1.04 | -2.08 | 0.90 | 1.20 |
| visual_association | -0.12 | -1.26 | 0.33 | -1.33 | -1.33 | -1.30 | -0.35 | -1.52 | -0.45 | -0.47 | 0.38 | 0.75 |

## Results (descriptive)
### Cross-species concordance: Spearman rho of human z vs macaque z over the shared groups
| variable | primary (n=6) | +Chugani ACC (n=7) | +Chugani ACC +Meguro IPL (n=8) | all groups with receptors in both species (n=13; CMRglc n=6) |
|---|---|---|---|---|
| iGlu | -0.03 (n=6) | 0.29 (n=7) | 0.29 (n=8) | 0.26 (n=13) |
| GABA | -0.14 (n=6) | -0.04 (n=7) | 0.14 (n=8) | 0.08 (n=13) |
| ACh (M1+M2+M3) | 0.20 (n=6) | 0.32 (n=7) | 0.52 (n=8) | 0.54 (n=13) |
| monoamine GPCR | 0.09 (n=6) | 0.32 (n=7) | 0.38 (n=8) | 0.19 (n=13) |
| E/I index | 0.71 (n=6) | 0.82 (n=7) | 0.67 (n=8) | 0.24 (n=13) |
| total (14 shared receptors) | -0.09 (n=6) | 0.11 (n=7) | 0.24 (n=8) | 0.27 (n=13) |
| CMRglc | -0.89 (n=6) | -0.36 (n=7) | -0.33 (n=8) | -0.89 (n=6) |

The `excluding homology_confidence == low` sensitivity is identical to the primary row: none of the 6 shared groups is
rated low (3 high, 3 moderate). Noda 2002 contributes to no shared group, so "excluding Meguro" is the only macaque
sensitivity that changes the set.

### Within-species Spearman: receptor class vs CMRglc, all groups of the species with both variables (human n = 12, macaque n = 9)
| receptor class | human n | human rho | human p_naive | macaque n | macaque rho | macaque p_naive |
|---|---|---|---|---|---|---|
| iGlu | 12 | -0.47 | 0.124 | 9 | -0.49 | 0.185 |
| GABA | 12 | -0.36 | 0.245 | 9 | -0.33 | 0.391 |
| ACh (M1+M2+M3) | 12 | -0.35 | 0.265 | 9 | -0.33 | 0.391 |
| monoamine GPCR | 12 | -0.34 | 0.286 | 9 | -0.62 | 0.075 |
| E/I index | 12 | -0.03 | 0.931 | 9 | -0.76 | 0.017 |
| total (14 shared receptors) | 12 | -0.42 | 0.175 | 9 | -0.48 | 0.194 |
| ACh incl. alpha4beta2 (human only) | 12 | -0.34 | 0.276 | — | — | — |
| total 15 (human only) | 12 | -0.42 | 0.175 | — | — | — |

### Within-species Spearman over the 6 primary shared groups only
| receptor class | human n | human rho | human p_naive | macaque n | macaque rho | macaque p_naive |
|---|---|---|---|---|---|---|
| iGlu | 6 | -0.54 | 0.266 | 6 | -0.37 | 0.468 |
| GABA | 6 | -0.66 | 0.156 | 6 | -0.37 | 0.468 |
| ACh (M1+M2+M3) | 6 | -0.14 | 0.787 | 6 | -0.37 | 0.468 |
| monoamine GPCR | 6 | -0.03 | 0.957 | 6 | -0.54 | 0.266 |
| E/I index | 6 | 0.20 | 0.704 | 6 | -0.37 | 0.468 |
| total (14 shared receptors) | 6 | -0.54 | 0.266 | 6 | -0.54 | 0.266 |
| ACh incl. alpha4beta2 (human only) | 6 | -0.14 | 0.787 | — | — | — |
| total 15 (human only) | 6 | -0.54 | 0.266 | — | — | — |

### Reading
- Within each species the four receptor classes correlate **negatively** with CMRglc at the group level (human rho -0.34
  to -0.47 over 12 groups; macaque -0.33 to -0.62 over 9). In macaque this is driven by V1 (highest 2-DG rate, lowest
  per-neuron density of every class) and in human mainly by premotor (third-highest Heiss 1985 rate of the 12 groups; lowest Zilles iGlu, GABA and ACh density, but the highest monoamine-GPCR density, hence the weaker monoamine rho).
- Across species the group **ranking of CMRglc is anti-concordant over the 6 shared groups (rho = -0.89)**: human Heiss
  1985 ranks S1/M1 > PM > OFC > DLPFC > V1 > VisAssoc, macaque 2-DG ranks V1 > VisAssoc > DLPFC > OFC > S1/M1 ~ PM. With the
  two filled groups the value is -0.33 to -0.36. The direction of the receptor-class concordances is weak and unstable
  (|rho| <= 0.3 for iGlu, GABA, monoamine GPCR in the primary set; ACh 0.20 -> 0.54 as groups are added). The E/I index is
  the only receptor variable with a consistently positive concordance (0.67-0.82), i.e. the ordering of iGlu/GABA across
  the shared groups is similar in both species while the absolute class levels are not.
- None of these values should be read as evidence for or against a relation; the group set is small, the human CMRglc rests
  on 3 subjects, and the two species' receptor variables are not the same quantity (next section).

## Caveats (read before citing any number)
1. **Receptor units differ in kind, not just scale.** Human densities are per mg protein (tissue basis); macaque densities
   are per neuron (FW2023 divided autoradiographic densities by local neuron density). Areas with high neuron density (V1,
   S1) therefore rank low in the macaque per-neuron measure and need not in the human per-protein measure. Rank
   concordance of a receptor class across species thus mixes receptor distribution with neuron-density distribution.
   FW2023 Supplementary Table 2 is the SD of the same per-neuron values, not a per-tissue table; a per-tissue
   (fmol/mg protein) macaque receptor source — the FW2023 GitHub release or the EBRAINS macaque fingerprints — would be
   the proper counterpart to Zilles 2017 and is the first thing to try.
2. **Physiological state.** Human CMRglc is resting-awake FDG PET (Heiss 1985: sleep-deprived wakefulness, n = 3; Heiss
   2004: resting, eyes closed). Macaque 2-DG values (Kennedy, Palombo) come from conscious, chair-restrained animals a few
   hours after catheter anaesthesia; Noda is conscious PET; Meguro is anaesthetised PET in baboon. Anaesthesia lowers and
   flattens regional rates, which is why Meguro is sensitivity-only.
3. **Species.** Receptors: *Macaca fascicularis* (FW2023). Metabolism: *M. mulatta* (Kennedy, Palombo, Noda); *Papio
   anubis* (Meguro). The homology table counts Meguro among the "macaque" sources; it is a different genus.
4. **Method / resolution.** 2-DG autoradiography samples tissue at sub-mm resolution with the 2-DG lumped constant; FDG PET
   at 8-11 mm (1985) with the FDG lumped constant and partial-volume mixing. Absolute values are not comparable across
   species or methods; only within-species, within-dataset relative profiles were compared, after the simple Kennedy +
   Palombo mean (same unit and method) for the macaque primary.
5. **Group composition is asymmetric.** `premotor` is 1 human area (BA6) vs 13 macaque F-areas; `visual_association` is 6
   human vs 13 macaque areas including MT, DP, LOP; `primary_visual` is 1 vs 1. Unweighted means over member areas weight
   each parcellation's areas equally, not cortical surface.
6. **Homology confidence** of the shared groups: high for primary_visual, primary_somatosensory_motor, premotor; moderate
   for visual_association, dorsolateral_prefrontal, orbitofrontal_ventral_prefrontal (see the homology table's reasons).
7. **Extended fills** rescale a secondary dataset by a single factor over overlapping groups (12 groups for Chugani,
   3 for Noda, 13 for Meguro); regional profile differences between datasets are not corrected.
8. Human CMRglc z-scores mix Heiss 1985 (cortical gyral groups) and Heiss 2004 (lobes, subcortex) sources in one
   within-species z; `cmrglc_z_cortex` and `cmrglc_source` allow restricting to one source.

## Unresolved / next steps
- Owner review of the homology assignments flagged in `HUMAN_MACAQUE_REGION_GROUPS.md` (none of the 6 shared groups is
  among the low-confidence ones, but `visual_association` includes FW DP/LOP and `dorsolateral_prefrontal` maps Zilles 9+46
  to Heiss 'Sup. front. cortex').
- A per-tissue macaque receptor table (FW2023 GitHub data or the EBRAINS macaque fingerprints) to remove the per-neuron
  confound (caveat 1).
- Hansen et al. 2022 PET receptor maps (DK parcels) would add human receptor coverage for cingulate groups, giving
  `anterior_cingulate` / `posterior_cingulate_retrosplenial` a non-filled human receptor value — but the human CMRglc gap for
  those groups would remain unless a cingulate-resolving absolute-units human PET table is added.

## Human subcortical extension (added 2026-11-17; human only)

The EBRAINS receptor series (`Evo-M1-Trait-Data/PalomeroGallagher_etal_2018/`, 41 human
items) adds whole-structure receptor fingerprints for putamen, globus pallidus, the
mediodorsal and anterior thalamic nuclei (2 brains each, v1.0) and hippocampal CA1, CA2,
CA3 and DG (3 brains, v2). No macaque subcortical receptor data exist in either series
(Froudist-Walsh 2023 and the EBRAINS macaque deposits are cortical), so this extension is
human-only and the macaque panel is unchanged.

Tables: `tables/s1e/human_subcortical_receptors_cmrglc.csv` (receptor-class sums, fmol/mg
protein, joined to Heiss 2004 rCMRglc) and
`tables/s1e/human_subcortical_receptor_ratio_to_cortex.csv` (each nucleus / mean of the 44
Zilles 2017 areas, per receptor).

| structure | iGlu | GABA | ACh (M1-M3) | monoamine GPCR | total (14) | Heiss 2004 rCMRglc |
|---|---|---|---|---|---|---|
| Globus pallidus | 652 | 765 | 324 | 538 | 2296 | 21.9 |
| Putamen | 2833 | 3039 | 2249 | 1666 | 9963 | 40.2 |
| Thalamus, mediodorsal | 1564 | 3394 | 1059 | 1316 | 7477 | 36.6 |
| Thalamus, anterior | 1626 | 3634 | 1113 | 1588 | 8153 | — (no Heiss 2004 row) |
| Hippocampus (mean CA1-3, DG) | 1999 | n/a* | 1123 | 2051 | 7945 | 25.7 |
| Cerebral cortex (mean of 44 areas) | 2285 | 6493 | 1486 | 1418 | 11733 | 33.5 |

\* hippocampal datasets lack the GABAA/BZ site, so the GABA class sum is not comparable;
GABAA + GABAB only.

Descriptive Spearman over the five whole structures with both variables (GP, putamen, MD
thalamus, hippocampus, cortex mean): iGlu +0.70, ACh +0.70, total +0.50, monoamine GPCR
+0.30 (n = 5; GABA +0.20 over n = 4). The pattern is carried by the pallidum (lowest
receptor densities of any structure — 6–36 % of the cortical mean for most receptors —
and the lowest rCMRglc) against the putamen (highest rCMRglc, D1 3.2×, M2 3.5×,
alpha4beta2 3.5×, kainate 2.5× cortex). The cortex has the highest total receptor density
but only intermediate rCMRglc, and the thalamic nuclei combine low ionotropic-glutamate
(0.5–0.9× cortex) with high alpha1/alpha4beta2/M2 (2.6–3.8× cortex) and high rCMRglc.
With five points these are observations for the whole-brain design, not results; the
subcortical receptor values rest on two brains per nucleus.

## Rat extension (added 2026-11-17)

**What was added.** A third species, rat, enters the s1e_3 pipeline through two Evo-M1 inputs synced
into `data_raw/` (manifest rows `sokoloff_etal_1977_table3/4` and `palomerogallagher_etal_2018_rat*`):

- **CMRglc**: `Sokoloff_etal_1977_Table3.csv` — the founding [14C]deoxyglucose paper (Sokoloff et al.
  1977, J Neurochem 28: 897-916), Table 3: local cerebral glucose utilization in normal conscious
  male Sprague-Dawley rats, 33 structures (29 gray, 4 white), 30-min arm (n = 10) used; the 45-min
  arm (n = 10) is a replication of the same values (only the hypothalamus differs significantly) and
  Table 4 (conscious control n = 6 vs thiopental n = 8) is synced but not used. Units µmol/100 g/min,
  2-DG autoradiography, conscious animals lightly restrained in loose plaster casts — the same lab and
  method as Kennedy 1978 (macaque), so the macaque-rat pair is the one method-matched comparison in
  the whole design (the human side is FDG PET). The Sokoloff structures were assigned to the region
  groups in `human_macaque_region_group_homology.csv` (column `Sokoloff_etal_1977_Table3`; every
  structure exactly once; 26 groups). Three rat-only groups had to be added — `olfactory_cortex`,
  `thalamus_lateral_nuclei`, `thalamus_ventral_nuclear_group` — all low-confidence and absent from the
  primate datasets; see `HUMAN_MACAQUE_REGION_GROUPS.md`.
- **Receptors**: the 10 rat items of the EBRAINS receptor-fingerprint series (Palomero-Gallagher &
  Zilles; Wistar, n = 5, fmol/mg protein): Oc1B + Oc1M → `primary_visual`, Oc2L + Oc2M →
  `visual_association`, CA1 + CA2 + CA3 + DG → `hippocampus_dentate`; the two laminar items (CA stratum
  cellulare / moleculare) are excluded from the group means because they re-measure the CA fields
  layer by layer. **Only 3 region groups have rat receptor data**, and the hippocampal items carry 10
  of the 15 receptors (no D1, GABAA/BZ, GABAB, M3, alpha4beta2), so every class but iGlu is NA for the
  hippocampus. Receptor classes with CMRglc therefore exist in 2 rat groups (iGlu) or 1 (all other
  classes); the rat within-species rows in `species_group_spearman.csv` report n = 2 / 1 and no rho
  (rho is reported from n ≥ 3, p only from n ≥ 5). This is a data-availability limit, not a result.
  The missing input is the rat-cortex receptor chapter (Zilles & Palomero-Gallagher, "Receptor
  architecture of the rat cerebral cortex" / cortical receptor fingerprints in Paxinos, *The Rat
  Nervous System*), which tabulates the 15 receptors for the full set of rat cortical areas and
  would give the rat the 15-20 cortical groups the human and macaque receptor sets cover.

**Outputs.** `tables/s1e/species_group_panel.csv` (198 rows = 66 groups × 3 species; rat rows carry
`species = rat`, `value_sokoloff1977_30min`, `receptors_missing`, and four CMRglc sharing flags
`shared_cmrglc_human_macaque / human_rat / macaque_rat / all_three`), `species_group_spearman.csv`,
`species_group_coverage.csv`; the old `human_macaque_*` filenames are kept as byte-identical copies so
existing links work. Pre-existing human and macaque rows and all previously reported rho / p values are
unchanged (checked row by row after the run). Figures: `figs/s1e/human_macaque_panels.png` now has a
third (rat) row with an explicit n annotation; new `figs/s1e/cmrglc_concordance_three_species.png`.

**Three-way CMRglc concordance** (Spearman rho of within-species z-scores; z computed over all groups
of the species that have CMRglc — human 41 groups, macaque 45, rat 26):

| pair | shared groups | rho | excl. low-confidence groups | over the 17 groups all three share |
|---|---|---|---|---|
| human vs macaque | 27 | 0.47 | 0.47 (n = 25) | 0.56 |
| human vs rat | 17 | 0.60 | 0.55 (n = 15) | 0.60 |
| macaque vs rat | 22 | 0.91 | 0.93 (n = 19) | 0.91 |

The 17 groups in all three species: primary_visual, primary_somatosensory_motor, frontal_lobe,
parietal_lobe, hippocampus_dentate, amygdala, pallidum, nucleus_accumbens, basal_forebrain_septum,
lateral_geniculate, medial_geniculate, superior_colliculus, inferior_colliculus, substantia_nigra,
cerebellar_cortex, cerebellar_nuclei, white_matter_cerebral. Macaque-rat adds primary_auditory,
hypothalamus, mammillary_body, brainstem_nuclei and white_matter_brainstem_cerebellar_spinal (no human
PET ROI); the rat-only groups (olfactory_cortex, striatum_undivided, the two thalamic groups) are shared
with neither.

**Reading.** The two 2-DG datasets agree closely (rho 0.91 over 22 groups): the rank order of regional
glucose use in conscious rat and conscious macaque is nearly the same, with the inferior colliculus the
most active structure in both (rat z = 2.8, macaque z = 4.0), white matter lowest, and the medial
geniculate (z 0.9 / 1.2) the next-highest shared group in both. Each 2-DG species agrees with the
human FDG-PET profile only moderately (0.47, 0.60), and the disagreement is concentrated in the same
groups: the auditory/visual midbrain-thalamic relays (inferior colliculus, medial geniculate, superior
colliculus, lateral geniculate) and the cerebellar nuclei are high in both rodent and monkey
autoradiography but low in the human PET ROIs (human z −0.2 to −1.7; differences of 1.3-4.3 z). Those
are small nuclei whose PET values are partial-volume mixtures with surrounding white matter and CSF,
so the pattern is more plausibly a resolution/method effect than a human-specific reorganisation; the
method-matched macaque-rat agreement supports that reading. The cortical groups (S1/M1, V1, frontal,
parietal) and the large subcortical nuclei (pallidum, amygdala, hippocampus, septum/basal forebrain,
nucleus accumbens, cerebral white matter) sit near the identity line in all three pairs.

**Caveats.** (i) Rat receptor coverage is 3 groups, 2 with CMRglc; no rat receptor-vs-CMRglc
statement can be made until the rat-cortex receptor table is ingested. (ii) The Sokoloff lobe-level
labels (`Frontal cortex`, `Parietal cortex`) and `Caudate-putamen` are composites; the three new
thalamic/olfactory groups are rat-only. (iii) Sokoloff's `Hippocampus` is split into Ammon's horn
and dentate gyrus and averaged; Kennedy's `Hippocampal gyrus` may be parahippocampal (flagged
earlier). (iv) z-scores depend on which groups a dataset covers: the rat z is over 26 groups
including 4 that no other species has, the macaque over 45, the human over 41, so identical ranks do
not give identical z. (v) Sokoloff's rats were immobilised in plaster casts, Kennedy's monkeys sat
in restraining chairs; both are "conscious, restrained", but stress-related differences cannot be
excluded. (vi) Descriptive only — no spatial-autocorrelation null, no phylogenetic model, naive p
values; nothing here enters `run_all.R`.
