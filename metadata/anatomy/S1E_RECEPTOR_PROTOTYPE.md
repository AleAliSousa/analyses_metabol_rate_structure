# S1E prototype — transmitter receptors, GLI/cell-body fraction and cortical metabolism

**Status:** prototype, 2026-11-15. Nothing here is in `run_all.R` or
`metadata/output_manifest.csv`. All metabolic values are RELATIVE (see §3).
Companion note for the parcellation itself: `S1E_PARCELLATION_PROTOTYPE.md`.

## 1. What was added

### Evo-M1-Trait-Data items (4-file convention; registry rows emitted, owner pastes)

| Item | Source | Rows | Status |
|---|---|---|---|
| `Zilles_PalomeroGallagher_2017_SupplementaryTable2` | Front Neuroanat 11:78 Suppl. Table 2, 24-page PDF (text layer), `_snapshot.py` (pdfplumber; recovers the yellow/grey row highlights = primary sensory / M1) | 2640 = 44 areas × 15 receptors × 4 strata; fmol/mg protein; mean, SD, CV | registered by owner beforehand → public TSV `10.3389%2Ffnana.2017.00078_SupplementaryTable2.tsv` written; 7/7 invariants pass |
| `Zilles_PalomeroGallagher_2017_EBRAINSFingerprints` | EBRAINS receptor fingerprints via siibra 1.0.1a27 (Julich-Brain 3.0.3); 34 datasets, each with its own DOI (10.25493/…) | 531 = 42 fingerprints × 14–16 receptors (incl. mGluR2/3); 21 cortical areas + hippocampal CA1/2/3/DG strata | build stops at `Item encoded` until registered; **owner decision**: keep in the 2017 paper folder or move to a dedicated folder (datasets are separate publications) |
| `Hansen_etal_2022_PET-parcellated-scale033` | authors' repo `netneurolab/hansen_receptors` @ f8b41da, `data/PET_parcellated/scale033` (digital-native tabulation; the paper prints no region table) | 3403 = 83 regions (68 Desikan-Killiany cortical + 15 subcortical, per hemisphere) × 41 tracer maps | build stops at `Item encoded` until registered; units per map in `reference_tables/…_tracer_units.csv` (BPND 14, SUVR 5, Bmax 4, VT 4, **see_source 14**); 82/83 regions join Johansen Table 2 DK labels |

Item names deviate from the plan because the registry `Item name` formula strips
underscores (hyphens survive; precedent Capellini `sleep-dataFemale`).
`_keys/General/species_key.csv` gained one row (Homo sapiens / Zilles_PalomeroGallagher_2017).

QA worth knowing: the Table S2 extraction matches Hansen's independent transcription
(`data/autoradiography/ReceptData*.npy`) in 660/660 all-layer, 660/660 supragranular,
660/660 infragranular and 645/645 printed granular cells (the 15 remaining cells are
area 24 granular, printed `--`, stored as 0 by Hansen; kept blank here). The printed CV
column is **not** SD/mean for 1434 of 2625 cells (median |diff| 0.22); carried as
printed with role `info`, not entered as an erratum (basis unknown, not a transcription error).

### Analysis-repo scripts and tables

- `scripts/00_sync_evo_m1_inputs.R` + `metadata/evo_m1_input_manifest.csv`: three
  `s1e` rows added; the three CSVs now sit in `data_raw/`.
- `helpers/parcellate_metabolic_maps_julich.py` + `scripts/s1e_0_extract_metabolic_maps_julich.R`
  → `data_intermediate/s1e_julich_area_metabolism.csv` (438 rows: 146 Julich labels × L/R/bilateral),
  `tables/s1e/s1e_zilles_area_metabolism.csv`, `s1e_lobe_check_vs_heiss*.csv`,
  `s1e_crosswalk_coverage.csv`; inputs and SHA-256 in `s1e_julich_area_metabolism_manifest.csv`.
- `metadata/anatomy/s1e_zilles2017_julich_area_map.csv`: 44 Zilles areas → Julich-Brain v3.1:
  18 exact, 2 part_of, 6 union, 8 approximate, 10 none (20, 21, 36, 38, 37L, 37M, 9, 46, 23, 31).
- `scripts/s1e_1_join_receptors_metabolism.R` → `tables/s1e/s1e_receptor_metabolism_join.csv`
  (44 rows), `s1e_receptor_metabolism_spearman.csv` (252 rows), `s1e_ebrains_metabolism_join.csv`.
- Atlases: `data_raw/anatomical_atlases/JulichBrain_v31/` (MNI152 1 mm label volumes + LUT;
  fsaverage-164k label GIfTIs with embedded label table); `data_raw/neuromaps_metabolic/`
  (Raichle cmrglc/cmr02/cbf/cbv fsLR-164k, Castrillon 2023 CMRglc MNI 3 mm; manifest).

## 2. Hansen 2022 — decision

The paper's data are voxelwise group-average PET maps (19 receptors/transporters), so
there was nothing to snapshot in the paper itself. The authors' own parcellated CSVs
are the digital-native table and were taken as the frozen source; hemisphere-specific
values are kept, and the n-weighted averaging of same-receptor maps that Hansen applies
(`code/make_receptor_matrix.py`) is deliberately **not** done in the item (aggregation
deferred to analysis time, DeCasien-style). 14 of 41 maps have no stated unit in the
paper (Supplementary Table 3 would resolve most; not in the frozen download). Hansen's
parcellation is Desikan-Killiany, not Julich-Brain, so it is not joined in s1e_1; the
natural pairing is with Johansen Table 2 (DK, 82/83 regions) and with a DK parcellation
of the same neuromaps metabolic maps (not yet run). Hansen's maps cover subcortex
(thalamus, striatum, pallidum, hippocampus, amygdala, accumbens, brainstem), which is the
only whole-brain receptor route for the Heiss 26-region design.

## 3. Join coverage and first descriptive look

- Zilles areas with parcellated metabolism: **26/44 strict** (exact + part_of + union),
  **34/44 any** (incl. approximate). No matched area is flagged for too few vertices.
- EBRAINS fingerprints: 24 Julich anchors, 21 with cortical metabolism (hippocampal
  strata have no surface counterpart).
- Lobe sanity check of the parcellation against Heiss 2004 rCMRglc (from
  `s1e_lobe_check_vs_heiss.csv`): Raichle cmrglc aggregated over all Julich areas
  r = 0.97 across five lobes (rescaled frontal 34.60 vs 35.3, parietal 33.51 vs 35.8,
  temporal 30.27 vs 30.5, occipital 34.47 vs 35.8, insular 29.82 vs 30.3); over the
  Zilles matches r = 0.85 (any) / 0.57 (strict; temporal lobe rests on FG1, FG2, 41).
  The Castrillon 3 mm volume parcellation does **not** reproduce the Heiss ordering
  (negative r) and should not be used at area level until projected to the surface or
  partial-volume corrected.
- Raichle maps as distributed by neuromaps are not in µmol/100 g/min (cortical mean
  6463); treat as relative. Castrillon cortical mean 25.87 is of the absolute order.

Spearman correlations (all-layers receptor density vs. relative metabolism; descriptive,
**no spatial-autocorrelation null, no multiple-comparison control**; strict basis n = 26,
any-match n = 34; full table in `s1e_receptor_metabolism_spearman.csv`):

| outcome | strongest predictors (strict n=26) | rho | any-match n=34 rho |
|---|---|---|---|
| cmrglc | M2 | +0.41 | +0.49 |
| cmrglc | 5-HT1A | −0.37 | −0.50 |
| cmrglc | alpha2 | +0.37 | +0.36 |
| cmrglc | GABAB | −0.36 | −0.37 |
| cmr02 | 5-HT1A | −0.60 | −0.63 |
| cmr02 | M2 | +0.57 | +0.54 |
| cmr02 | GABAA | +0.51 | +0.48 |
| cmr02 | GABAB / M3 / kainate / alpha1 | −0.46 / −0.46 / −0.45 / −0.43 | −0.56 / −0.55 / −0.53 / −0.49 |
| gi_proxy | NMDA | −0.63 | −0.59 |
| gi_proxy | GABAA | −0.56 | −0.47 |
| gi_proxy | D1 | −0.55 | −0.36 |

Class aggregates are weak on cmrglc (iGlu +0.05, GABA −0.11, ACh −0.22, monoamine GPCR
+0.02, receptor E/I index +0.13; strict basis). Reading: the sensory→association
receptor axis (high M2/GABAA in primary sensory areas, high 5-HT1A/GABAB/M3 in
association cortex) is what these correlations track; the energy-budget prediction
(ionotropic glutamate ↔ CMRglc) is not visible at this resolution, and the Castrillon
2023 GPCR/neuromodulator result is only partly echoed (5-HT1A negative, M2 positive).
None of this is a finding until the hierarchy axis is partialled out and spin/BrainSMASH
nulls are applied (the repo already has `01_network_residual_autocorrelation_analysis.R`
for the regional design).

## 4. Limitations carried forward

1. fsLR→fsaverage resampling used an in-house barycentric resampler (Connectome
   Workbench not installable on macOS arm64 sandbox); the helper switches to
   `wb_command` when present. Not yet compared numerically.
2. Template correspondence (Julich fsaverage MPM vs neuromaps fsLR projection; ICBM 2009c
   labels vs 3 mm grid) unvalidated beyond the medial-wall and lobe checks.
3. Zilles 2017 sampling positions inside large unions (areas 6, 8, 24, 32) are not
   encoded; the union pools the whole Julich territory. Approximate matches (V3A, 37B,
   42, 22, 10L, 10M, 11, 47) belong in sensitivity analyses only.
4. Receptor densities are per mg protein; conversion to per-volume / per-neuron needs
   the Sjölin total-protein values (s1d) and regional neuron densities (s1a) — not done.
5. Three brains (72–77 y) for the autoradiography vs. young-adult PET cohorts.

## 5. Follow-up (new analyses folder, layer-resolved)

- EBRAINS/siibra exposes for Julich-Brain v3.1: 7 receptor-density **profiles**
  (laminar, areas incl. 44, 45, 3b, hOc1/hOc2, PFt, PGa), 12 `LayerwiseCellDensity`
  features and 1 `CellDensityProfile` (BigBrain-based cell-body densities by layer —
  the GLI-type quantity that answers the soma-fraction question), plus the strata
  columns already in Table S2 (supragranular / granular / infragranular). A laminar
  study would pair these with laminar metabolic proxies (no laminar PET exists; BigBrain
  intensity profiles, Wagstyl et al. 2020, are the cell-body-fraction side).
- DK-parcellate the same metabolic maps (neuromaps `Parcellater` with the Cammoun
  scale033 atlas Hansen used) so Hansen 2022 + Johansen 2024 + metabolism share one
  83-region frame including subcortex.
- Register the two pending items (owner), then re-run their build scripts to write the
  public TSVs; decide the folder for the EBRAINS item.

## 6. Middle tier (Heiss 1985 groups)

A second, coarser join of the Zilles & Palomero-Gallagher 2017 receptor areas to
**absolute** rCMRglc (µmol/100 g/min) at the gyral-group resolution of
Heiss et al. 1985 (Brain Res 327:362–366, FDG-PET, 4 healthy male volunteers: 3
non-dreamers + 1 dreamer, sleep-deprived waking condition; 7-slice scanner, ~8 mm in-plane
FWHM, 11 mm slices), with the Chugani et al. 1987 adult reference group (19–30 y, n = 7,
Ann Neurol 22:487–497, Table 2) as a second outcome where a Heiss group has a Chugani
counterpart. It sits between the lobe check (s1e_0, 4 Heiss 2004 lobes) and the
vertex-parcellated relative maps (s1e_1, 25–33 areas).

### Files

- `metadata/evo_m1_input_manifest.csv`: rows `heiss_etal_1985_tablei` and
  `chugani_etal_1987_table2` (study s1e); synced to `data_raw/Heiss_etal_1985_TABLEI.csv`
  and `data_raw/Chugani_etal_1987_Table2.csv`.
- `metadata/anatomy/s1e_zilles2017_heiss1985_group_map.csv` — 44 rows (one per Zilles
  area): `heiss1985_structure_as_published`, `heiss1985_group_id`,
  `chugani1987_structure_as_printed` (NA where none), `heiss2004_anatomy_id` (lobe
  vocabulary copied from `s1e_zilles2017_julich_area_map.csv`), `assignment_relationship`
  ∈ {contained, partial, approximate}, `note` (every row).
- `scripts/s1e_2_middle_tier_heiss1985.R` → `tables/s1e/s1e_middle_tier_group_receptors.csv`
  (14 groups × 4 strata = 56 rows; 15 receptors, class sums, EI index, total, outcomes),
  `tables/s1e/s1e_middle_tier_spearman.csv` (252 rows), `tables/s1e/s1e_middle_tier_group_coverage.csv`
  (all 24 Heiss 1985 structures). Not in `run_all.R` / `output_manifest.csv`.

### Heiss 1985 column used

TABLE I prints four arms: non-dreamers awake (n = 3, mean ± SD), non-dreamers % change
asleep, dreamer awake (n = 1), dreamer % change asleep. Only the **awake** arms are used.
Primary outcome `cmrglc_heiss1985_awake_nd3` = non-dreamers awake mean with its printed SD
(the only awake arm with a dispersion). Secondary `cmrglc_heiss1985_awake_pooled4` =
n-weighted mean of all 4 awake subjects, (3 × non-dreamer mean + dreamer value)/4 — the pool
Karbowski 2007 used; no SD can be derived for it. Sleep values exist only as % change and
are not used. Chugani outcome `lcmrglc_chugani1987_adult` = adult-group mean (± SEM).

### Method

Unweighted mean of the member Zilles areas per group × stratum × receptor (`na.rm`; area 24 is
agranular and has no granular-stratum values, so the cingulate granular mean is over 23/31/32 —
recorded as `n_areas_nonmissing_min`). Class aggregates as in s1e_1: iGlu = AMPA+NMDA+kainate,
GABA = GABAA+GABAA/BZ+GABAB, ACh = M1+M2+M3+α4β2, monoamine GPCR = α1+α2+5-HT1A+5-HT2+D1,
EI = iGlu/GABA. Spearman ρ per predictor × outcome × stratum; `p_naive` is reported only when
n ≥ 8 groups (Chugani has 5 groups → ρ only). Descriptive; no spatial null.

### Coverage (44/44 Zilles areas assigned; 14 of 24 Heiss structures receive receptor data)

| group | Heiss 1985 structure | Zilles areas (n, n contained) | Chugani 1987 counterpart | awake nd3 (SD) | pooled4 |
|---|---|---|---|---|---|
| sensorimotor | Sensorimot. cortex | 1, 2, 3a, 3b, 4 (5, 5) | Sensorimotor cortex (match) | 44.2 (3.07) | 42.4 |
| premotor | Premotor cortex | 6 (1, 1) | — | 42.6 (1.83) | 41.55 |
| sup_frontal | Sup. front. cortex | 8, 9, 10L, 10M (4, 0) | — | 41.4 (1.99) | 38.55 |
| inf_frontal | Inf. front. cortex | 44, 45, 46, 47 (4, 2) | — | 44.6 (2.35) | 43.3 |
| orbital | Orbital cortex | 11 (1, 1) | — | 42.1 (1.85) | 40.65 |
| cingulate | Cingulate gyrus | 23, 24, 31, 32 (4, 3) | Anterior cingulate cortex (partial) | 44.6 (2.86) | 43.7 |
| inf_parietal | Inf. par. cortex | PFt, PFm, PGa, PGp (4, 4) | — | 38.2 (2.68) | 36.4 |
| sup_parietal | Sup. par. cortex | 5L, 5M (2, 1) | — | 38.8 (1.75) | 37.525 |
| sup_temporal | Sup. temp. cortex | 22, 41, 42 (3, 3) | Transverse temporal cortex (partial) | 41.3 (2.11) | 40.175 |
| middle_temporal | Middle temp. cortex | 21 (1, 1) | — | 38.8 (3.00) | 37.6 |
| inf_temporal | Inf. temp. cortex | 20, 37L, 37B, 37M (4, 1) | — | 35.4 (2.92) | 34.525 |
| inf_limbic | Inf. limbic struct. | 36, 38 (2, 0) | — | 32.3 (2.50) | 31.2 |
| prim_visual | Prim. vis. cortex | V1 (1, 1) | Calcarine cortex (match) | 40.1 (4.09) | 38.95 |
| vis_assoc | Vis. assoc. cortex | V2d, V2v, V3d, V3A, V3v, V4v, FG1, FG2 (8, 6) | Occipital cortex (approximate) | 36.1 (3.07) | 34.35 |

Relationship counts over the 44 areas: 29 contained, 6 partial, 9 approximate. Heiss
structures without a Zilles area: Insular cortex, Hippocampus, both white-matter regions,
Caudate, Lentiform, Thalamus, Brainstem, Cerebellar gray, Whole brain (10). Decisions
recorded in the crosswalk notes: frontal pole 10L/10M → Sup. front. (not Orbital); FG1/FG2 →
Vis. assoc. (posterior fusiform, "higher visual areas" in Zilles 2017) while 37B/37M → Inf.
temp. (mid-fusiform, temporal-lobe area 37); temporal pole 38 and perirhinal 36 → Inf. limbic
struct. (approximate; composition of that 1985 ROI is not described in the paper); 46 → Inf.
front. (no middle-frontal region in Heiss 1985); 8, 9 → Sup. front. (partial, MFG portion
uncovered).

### Top |ρ| vs `cmrglc_heiss1985_awake_nd3`, all_layers (n = 14 groups, descriptive)

| predictor | ρ | p_naive |
|---|---|---|
| α4β2 (nicotinic) | +0.583 | 0.029 |
| NMDA | −0.564 | 0.036 |
| GABAA | −0.370 | 0.19 |
| α2 | +0.370 | 0.19 |
| 5-HT2 | −0.368 | 0.20 |
| class_ACh | −0.368 | 0.20 |
| D1 | −0.346 | 0.23 |
| M1 | −0.334 | 0.24 |

The same two receptors lead at every stratum (granular: α4β2 ρ = +0.795, p = 0.0007; GABAA
−0.593; NMDA −0.581; supragranular: α4β2 +0.690, NMDA −0.619; infragranular: α4β2 +0.410,
NMDA −0.385). With the pooled 4-subject outcome the ordering is unchanged but weaker
(α4β2 +0.484, NMDA −0.481, both p_naive ≈ 0.08). The α4β2 result is driven by the three
frontal groups (sup_frontal 78.3, cingulate 75.3, inf_frontal 73.0 fmol/mg protein, all
with rCMRglc ≥ 41.4) against temporal/parietal groups near 34–40; NMDA is highest in V1
(1735) and the inf_limbic/vis_assoc groups, which have mid-to-low 1985 rCMRglc. Chugani
(5 groups) yields ρ values of ±0.9–1.0 for several predictors that are not interpretable
at n = 5 and are left without p-values.

### Caveats

1. n = 4 subjects (3 in the primary arm), sleep-deprived waking state, 1985 PET (≈8 mm
   in-plane, 11 mm slices): a 14-point outcome with SDs of 2–4 µmol/100 g/min on
   between-group differences of similar size.
2. Unweighted area means: a group's value is the mean of whatever Zilles areas happen to
   fall in it, not an area- or volume-weighted estimate of the gyral region; groups range
   from 1 area (orbital, premotor, middle temporal, primary visual) to 8 (visual
   association), so single-area groups carry the same weight as 8-area means.
3. Group-size imbalance and 9 approximate assignments (10L, 10M, 46, FG1, FG2, 37B, 37M,
   36, 38) — the inf_limbic group rests entirely on two approximate assignments; drop it in a
   sensitivity run.
4. 14 groups, 21 predictors × 3 outcomes × 4 strata: p_naive values are nominal, no
   multiple-comparison or spatial-autocorrelation correction; neighbouring groups share
   subjects, scanner partial-volume effects and cytoarchitectonic gradients.
5. Receptor densities are per mg protein from 3 brains aged 72–77 y; PET cohorts are young
   adults (Heiss 28–38 y, Chugani 19–30 y).
