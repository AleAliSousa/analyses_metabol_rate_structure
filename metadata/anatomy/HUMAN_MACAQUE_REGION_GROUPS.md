# Human <-> macaque region-group homology table (DRAFT for owner review)

Created 2026-11-15. Companion to `human_macaque_region_group_homology.csv` in this folder.
**Status: draft.** This table maps region *names* to shared region groups. It joins **no numbers**. Nothing here
enters `run_all.R` or `output_manifest.csv`.

## Purpose
Provide one crosswalk from the region vocabularies of the human and macaque metabolic (rCMRglc), receptor and
neuron-density datasets held in Evo-M1-Trait-Data to a common set of region groups, so that later scripts can
decide, group by group, which datasets can be compared within species and which comparisons across species are
anatomically defensible.

## Design
- One row per region group (66 groups: 28 cortex, 33 subcortex, 2 cerebellum,
  3 white matter; 63 original + 3 rat-only groups added 2026-11-17, see "Rat extension" below). Groups were adapted to the data: lobe-level composites (`frontal_lobe`, ...), undivided
  composites (`thalamus_undivided`, `lentiform_nucleus`, `striatum_undivided`, `cingulate_undivided`, `brainstem_undivided`) and
  nucleus-level groups are kept **separate** from the area-level groups, following `DATASET_REGION_COVERAGE.md`
  (no filling from parents, children or neighbours).
- Columns `group_id`, `group_label`, `tissue_class` (cortex / subcortex / white_matter / cerebellum; hippocampus and amygdala
  are filed under subcortex), then **one column per dataset** holding the dataset's region names **verbatim**
  (';'-separated, blank if absent), then coverage counts and flags, `homology_confidence` (high / moderate / low with a
  one-clause reason) and `note`.
- Every name of every vocabulary below is assigned to exactly one group (checked programmatically); no name was edited,
  including the printed misspelling `Dentate nculeus` (Palombo Table 5).

### Dataset columns and their vocabularies
- `Heiss_etal_2004` (26 names)
- `Heiss_etal_1985_TABLEI` (24 names)
- `Chugani_etal_1987_Table2` (14 names)
- `Borghammer_etal_2012_Table2` (10 names)
- `Andersen_etal_2022_TableS2` (14 names)
- `Zilles_PalomeroGallagher_2017_areas` (44 names)
- `Hansen_etal_2022_DK_scale033` (42 names)
- `Kennedy_etal_1978_Table` (39 names)
- `Palombo_etal_1990_Table1` (18 names)
- `Palombo_etal_1990_Table2` (15 names)
- `Palombo_etal_1990_Table3` (14 names)
- `Palombo_etal_1990_Table4` (15 names)
- `Palombo_etal_1990_Table5` (20 names)
- `Noda_etal_2002_Table3` (7 names)
- `Meguro_etal_1999_Table2_baboon` (19 names)
- `Sokoloff_etal_1977_Table3` (33 names; rat, added 2026-11-17)
- `FroudistWalsh_etal_2023_SuppTable1_areas` (109 names)
- `EBRAINS_rat_receptors` (8 names; rat, added 2026-11-17)
- `Collins_etal_2010_DatasetS1_area_pieces_galago_owlmonkey` (12 names)
- `Collins_etal_2016_Table1_chimpanzee` (7 names)

Notes on specific columns:
- `Heiss_etal_2004`: `canonical_label` of `data_intermediate/heiss_2004_regions.csv` (26 rows, lobe-level cortex + subcortex + white matter).
- `Zilles_PalomeroGallagher_2017_areas`: the 44 `area` values of `Zilles_PalomeroGallagher_2017_SupplementaryTable2.csv`;
  `s1e_zilles2017_julich_area_map.csv` gives the corresponding Julich-Brain names (not repeated here).
- `Hansen_etal_2022_DK_scale033`: the 42 distinct `region_label` values (83 rows = left/right x 41 + single `brainstem`); Desikan-Killiany
  cortex + aseg subcortex; no cerebellum.
- `Palombo_etal_1990_Table1..5`: one column per printed table (Table 1 basal ganglia/nigra, Table 2 thalamus, Table 3 cortex +
  hippocampus + two WM tracts, Table 4 basal forebrain/hypothalamus, Table 5 brainstem/cerebellum). The five tables count as **one**
  macaque metabolic dataset in `n_macaque_metabolic_datasets`.
- `Meguro_etal_1999_Table2_baboon`: **baboon (*Papio*), not macaque, and anaesthetised**; counted in `n_macaque_metabolic_datasets`
  because the task lists it with the macaque metabolic sources, but it should be treated as a separate species in any analysis.
- `FroudistWalsh_etal_2023_SuppTable1_areas`: 109 area names extracted with pdfplumber from the first column of Supplementary
  Table 1 (41593_2023_1351_MOESM1_ESM.pdf, pp. 13-16; Julich Brain Macaque Maps). The parcellation covers frontal, parietal,
  occipital and cingulate cortex plus MT; **no temporal (other than MT) and no insular areas**, so the temporal and insular groups
  have `receptor_macaque = n`.
- `Collins_etal_2010_DatasetS1_area_pieces_galago_owlmonkey`: the task lists Collins 2010 as macaque neuron density, but in Dataset S1
  the **macaque hemisphere (08-59) is an unlabelled numbered grid** (pieces 1-41 plus sub-piece `29a`); only the galago #2 and owl-monkey
  hemispheres carry area labels (V1, V2, dV3, vV3, S1, M1, A1, MT, premotor, DL, Remain Ctx, Surr Ctx). Those 12 labels are what the
  column holds; `29a` (a macaque grid sub-piece typed `area` in the CSV) was excluded. The macaque grid cannot be assigned to groups
  without the dissection map.
- `Collins_etal_2016_Table1_chimpanzee`: this item is **one chimpanzee hemisphere**, not macaque (README). Its 7 `structure` values are
  included as the only area-level neuron-density vocabulary in the Collins items; `somatosensory block`, `premotor block` and
  `Prefrontal Cortex` are dissection blocks.

- `Sokoloff_etal_1977_Table3`: the 33 `structure_as_published` values of `Sokoloff_etal_1977_Table3.csv` (Evo-M1; 2-DG, conscious
  Sprague-Dawley rats, 30- and 45-min arms; the analysis uses the 30-min arm). Rat-specific labels and their placement: `Caudate-putamen` ->
  `striatum_undivided` (the rat has no internal capsule separating caudate and putamen; exact fit to the composite group); `Olfactory cortex`
  -> new group `olfactory_cortex`; `Thalamus: lateral nucleus` -> new group `thalamus_lateral_nuclei`; `Thalamus: ventral nucleus` -> new group
  `thalamus_ventral_nuclear_group` (the rat label spans VA/VL and VPL/VPM, which the macaque table keeps as `thalamus_ventral_motor` and
  `thalamus_ventral_posterolateral`; no guess was made); `Hippocampus: ammon's horn` and `Hippocampus: dentate gyrus` -> `hippocampus_dentate`;
  `Septal nucleus` -> `basal_forebrain_septum`; `Lateral lemniscus` and `Pontine gray matter` -> `brainstem_nuclei` with the vestibular, cochlear
  and superior olivary nuclei (Kennedy's pontine reticular nucleus sits there too); `Genu of corpus callosum` -> `white_matter_cerebral` with
  `Corpus callosum` and `Internal capsule`; `Cerebellar white matter` -> `white_matter_brainstem_cerebellar_spinal`; `Parietal cortex` /
  `Frontal cortex` -> the lobe-level composites, as for Kennedy. 26 groups receive at least one Sokoloff structure.
- `EBRAINS_rat_receptors`: the `structure_as_published` values of the 10 rat items of `PalomeroGallagher_etal_2018` (EBRAINS receptor
  fingerprints, Wistar, n = 5, fmol/mg protein): `Area Oc1B;Area Oc1M` -> `primary_visual` (binocular and monocular parts of Oc1 = V1),
  `Area Oc2L;Area Oc2M` -> `visual_association` (lateral and medial secondary visual cortex), `CA1 (Hippocampus);CA2 (Hippocampus);CA3
  (Hippocampus);DG (Hippocampus)` -> `hippocampus_dentate`. The two laminar items `CA, stratum cellulare (Hippocampus)` and `CA, stratum
  moleculare (Hippocampus)` are NOT listed: they are layer-wise re-measurements of the same CA fields and would double-count CA1-3 at a
  different grain in a whole-structure mean. The four hippocampal items carry 10 of the 15 receptors (no D1, GABAA/BZ, GABAB, M3,
  alpha4beta2), so only the iGlu class (AMPA + NMDA + kainate) is complete for `hippocampus_dentate`.

### Coverage flags
- `n_rat_metabolic_datasets`: 1 if Sokoloff 1977 Table 3 has a structure in the group (26 groups), else 0.
- `receptor_rat`: y if an EBRAINS rat item is listed in the group (3 groups).
- `n_human_metabolic_datasets`: count among Heiss 2004, Heiss 1985, Chugani 1987, Borghammer 2012, Andersen 2022 (max 5).
- `n_macaque_metabolic_datasets`: count among Kennedy 1978, Palombo 1990 (all tables = 1), Noda 2002, Meguro 1999 (baboon) (max 4).
- `receptor_human` = y if Zilles 2017 or Hansen 2022 has a name in the group (`receptor_human_source` says which);
  `receptor_macaque` = y if FW2023 has a name in the group.
- `fully_covered_both_species` = y when all four of the above are non-zero.

## Physiological state of each metabolic dataset
| dataset | species | measurement | physiological state |
|---|---|---|---|
| Heiss_etal_2004 (TABLE1 via heiss_2004_regions.csv) | human | FDG PET, 26 ROIs | resting awake; eyes closed, ears unplugged, low light/noise; subjects stayed awake (paper Methods) |
| Heiss_etal_1985_TABLEI | human | FDG PET, 24 ROIs; 4 subjects (3 non-dreamers + 1 dreamer) | awake arm = sleep-deprived wakefulness (~28 h on duty), supine, eyes closed, dimmed room; separate sleep arm (stages I-IV, one REM/dreaming subject) printed as % change (README, ANATOMICAL_MATCHING.md) |
| Chugani_etal_1987_Table2 | human | FDG PET, 14 ROIs; infants-adults | awake, kept alert by EEG-monitored tapping; paediatric subjects with transient neurological events, 18/29 on anticonvulsants (paper Methods) |
| Borghammer_etal_2012_Table2 | human | FDG PET, 10 subcortical ROIs; normalised (WM / GM / reference-cluster) values only | resting awake, supine, eyes closed, quiet room (paper Methods) |
| Andersen_etal_2022_TableS2 | human | FDG PET SUVR (CSO-normalised) + [11C]UCB-J; 14 composite ROIs | resting awake in a quiet, dimly lit room for >=20 min after injection; eyes state not extracted (paper Methods) |
| Kennedy_etal_1978_Table | macaque (Macaca mulatta) | 2-DG autoradiography, 39 structures | conscious; restraining chair, 2-5 h after recovery from halothane/N2O catheter anaesthesia (README title; paper Methods) |
| Palombo_etal_1990 Tables 1-5 | macaque | 2-DG autoradiography; control hemispheres of an MPTP (unilateral intracarotid) study | conscious, chaired, resting; >=3 h recovery from halothane/N2O catheter anaesthesia (paper Methods) |
| Noda_etal_2002_Table3 | macaque (rhesus) | FDG PET, 7 ROIs; young vs aged | conscious; chair-restrained with head fixation after >1 month habituation; catheter placement under transient anaesthesia, discontinued before scans (paper title/Methods) |
| Meguro_etal_1999_Table2 | baboon (Papio), NOT macaque | FDG PET, 19 ROIs; baseline and post-rhinal-lesion/sham time points | anaesthetised (phencyclidine infusion, muscle relaxation, light anaesthesia/analgesia) (paper Methods) |
| Sokoloff_etal_1977_Table3 | rat (Rattus norvegicus; male Sprague-Dawley, 325-450 g) | 2-DG autoradiography, 33 structures; 30-min (n = 10) and 45-min (n = 10) arms | conscious, lightly restrained: femoral catheters under light halothane, then immobilised in loose-fitting, partially bisected plaster casts, >= 2 h recovery before the experiment, kept alert by tactile stimulation as required (paper Methods p. 899). Sokoloff 1977 Table 4 (same lab; conscious control n = 6 vs thiopental n = 8) is synced but not used. |

Human values are resting-awake PET (Heiss 1985 also has a sleep arm, printed as % change); macaque 2-DG values are from conscious,
chair-restrained animals a few hours after catheter anaesthesia; Noda is conscious PET; Meguro (baboon) is anaesthetised PET.
Anaesthesia lowers and flattens rCMRglc, so Meguro values are not state-comparable with the others.

## Comparability caveat (read before joining any values)
Macaque 2-DG autoradiographic rates (Kennedy, Palombo; umol/100 g/min from tissue-level autoradiograms with the 2-DG lumped constant)
and human FDG-PET rates (Heiss, Chugani; umol/100 g/min at ~8-11 mm resolution with the FDG lumped constant and partial-volume
mixing), as well as the normalised-only human tables (Borghammer: WM/GM/reference-cluster ratios; Andersen: CSO-normalised SUVR)
and the Noda/Meguro PET tables, are **not unit-comparable across species or methods**. Any cross-species use of this table must first
normalise within species and within dataset (e.g. ratio to a within-dataset reference group or z-score across groups present in that
dataset) and only then compare the resulting relative profiles across the groups both species share. Receptor densities
(Zilles 2017 fmol/mg protein; Hansen 2022 PET tracer binding in heterogeneous units; FW2023 fmol per neuron) are likewise
dataset-specific scales and need the same within-dataset normalisation.

## Groups flagged for owner review (assignments I was unsure about)
- `temporal_pole_perirhinal_entorhinal`: Heiss 1985 `Inf. limbic struct.` placed here provisionally; the item does not document its
  composition (could include amygdala / parahippocampal gyrus). Meguro `Anterior temporal` placed here because the study's ROIs
  surround the rhinal cortex; could equally be `inferior_middle_temporal`.
- `hippocampus_dentate`: Kennedy `Hippocampal gyrus` may mean the parahippocampal gyrus rather than hippocampus proper.
- `inferior_parietal`: FW2023 `TSA` and `PPt` placed here only because FW2023 Fig. 5 lists them under 'Parietal'; their position in the
  Julich Brain Macaque Maps was not verified. Its macaque metabolic coverage rests on Meguro (baboon) alone.
- `visual_association`: FW2023 `DP` and `LOP` (dorsal prelunate / lateral occipito-parietal) placed with extrastriate visual cortex;
  `V6l`/`V6m` with visual, `V6A*` with `superior_parietal`. Kennedy `Visual cortex` taken as V1-dominated.
- `primary_somatosensory_motor`: Hansen `paracentral` (contains SMA) placed here; `premotor`: Hansen `caudalmiddlefrontal` (BA6/8) placed here
  while Zilles `8` and FW `8A*/8B*` form `frontal_eye_field_area8`.
- `orbitofrontal_ventral_prefrontal`: Hansen `parsorbitalis` filed with Zilles `47` (which maps to Julich Fo4-Fo7) rather than with the
  inferior frontal (Broca) group.
- `anterior_cingulate`: FW midcingulate areas (`a24'ab`, `a24'c`, `p24'ab`) kept with ACC because no human ROI separates MCC.
- `intraparietal_sulcus`: macaque receptor-only; human IPS areas are absent from the Zilles 44-area set and absorbed into DK gyri.
- `cerebellar_cortex`: Borghammer `Cerebellum`, Andersen `CBL`, Noda `Cerebellum`, Meguro `Cerebellar lobe` are undivided cerebellum
  (gray + white matter), filed with cerebellar cortex rather than as a separate `cerebellum_undivided` group.
- `nucleus_accumbens`: Palombo `Olfactory tubercle` included as ventral striatum.

## Rat extension (added 2026-11-17)
Three rat-only groups were added because a Sokoloff 1977 structure fitted none of the 63 groups: `olfactory_cortex` (cortex),
`thalamus_lateral_nuclei` and `thalamus_ventral_nuclear_group` (subcortex); all three are `homology_confidence = low` and have no human or
macaque entry. No existing group definition was changed. Rat coverage: CMRglc in 26 groups (17 shared with the human CMRglc set, 22 with
the macaque primary set, 17 with both); receptors in 3 groups (`primary_visual`, `visual_association`, `hippocampus_dentate`), of which
`visual_association` has no rat CMRglc structure (Sokoloff prints no secondary visual cortex) and `hippocampus_dentate` has a complete
receptor set for the iGlu class only. The rat physiological state is conscious, lightly restrained (plaster cast), 2-DG, i.e. the same
method and a comparable state to Kennedy 1978 / Palombo 1990 (conscious, chair-restrained). Rat receptor densities are fmol/mg protein, the
same unit as Zilles 2017 (human) but not the macaque fmol/neuron; the within-dataset normalisation rule above applies unchanged.

## Summary
- Groups with metabolism AND receptors in BOTH species (`fully_covered_both_species = y`): **9** of 63:
  primary_visual, visual_association, primary_somatosensory_motor, premotor, dorsolateral_prefrontal, orbitofrontal_ventral_prefrontal, anterior_cingulate, posterior_cingulate_retrosplenial, inferior_parietal. All are cortical (FW2023 is cortex-only). One of them (`inferior_parietal`) owes its macaque metabolic
  coverage to the baboon dataset alone; with Meguro excluded the count is **8**.
- Groups with metabolic data in both species (receptors not required): 40 of 63.
- Cortical groups: 27; with human receptors 19; with macaque receptors 14; macaque-receptor-only cortical groups
  lacking any macaque metabolic row: 4.

Draft prepared for owner review; no values were joined or compared.
