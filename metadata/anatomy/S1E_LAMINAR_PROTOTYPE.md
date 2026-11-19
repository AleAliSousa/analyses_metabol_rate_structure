# S1E laminar prototype — receptor strata × BigBrain cell-body density × GLI (human; cross-species GLI)

**Status:** prototype, 2026-11-18. Descriptive only. Nothing here is in `run_all.R` or
`metadata/output_manifest.csv`. Companion notes: `S1E_RECEPTOR_PROTOTYPE.md` (areal receptor join),
`S1E_HUMAN_MACAQUE_COMPARISON.md` (region groups). All numbers below are read from the tables named.

## 1. Inputs and files

Evo-M1 items synced through `metadata/evo_m1_input_manifest.csv` (21 `s1e` rows added; `Rscript scripts/00_sync_evo_m1_inputs.R`
refreshed them, provenance in `data_intermediate/evo_m1_input_provenance.csv`):

| Input (data_raw/) | Content | Unit of observation |
|---|---|---|
| `Zilles_PalomeroGallagher_2017_SupplementaryTable2.csv` | 15 receptors × 44 areas × 4 strata (supragranular I–III, granular IV, infragranular V–VI, all layers), fmol/mg protein | mean over 3 brains |
| `Dickscheid_etal_2021_Area{44,4a,7P,FG1,FG3,FG4,Fp1,hOc1,hOc2,hOc3d,hOc4d,hOc5}.csv` | BigBrain layer-wise cell-body density (# detected cells / 0.1 mm³, all cell types) and layer thickness (µm, shrunken tissue), layers I–VI + I–VI | one donor (BigBrain, 65-y male), 9–11 patches on 7–11 sections |
| `PalomeroGallagher_Zilles_2018_TableS3.csv` | GLI % by stratum, areas 44/45, 7 primate species (human n=4 brains, macaque n=2 (44) / 3 (45), others n=1) | brain |
| `Zilles_etal_1986_Table1.csv`, `_Table2.csv` | posterior cingulate areas 29, 30, 23, 31: laminar % thickness (T1) and GLI (T2) by sublayer (outer main = II–III, granular = IV for 23/31 only, inner main = V–VI) | one row per species, **17 species** (the task brief said 12; the item holds 17) |
| `Preuss_etal_1999_Fig2.csv`, `_MainText.csv`, `Preuss_Coleman_2002_Fig5.csv`, `_Fig7.csv`, `WongRiley_etal_1993_Text.csv`, `Hevner_WongRiley_1990_Table1.csv` | ordinal state rows (Preuss schema) / text rows for the layer-4A alignment | taxon statements |

Laminar receptor **profiles** (not items): 128 TSVs from the EBRAINS fingerprint datasets, downloaded untouched to
`<Evo-M1>/PalomeroGallagher_etal_2018/ebrains_downloads/profiles/` with `profiles_manifest.csv` (sha256). 8 human areas
(3b, 4p, 7A, 44d, 45, 47, PFm, hOc1) × 16 receptors (the 15 of Table S2 with `BZ` = GABAA/BZ, plus mGlu2/3); 101 samples at
0–100 % cortical depth (0 = pial surface, verified from the laminar shape of GABAA/AMPA in hOc1), fmol/mg. The EBRAINS folder is
`*_pr_examples/<receptor>/`: these are single **example** profiles (one brain/section), not dataset means. Anomaly: the 3b `M3/`
file is byte-identical to the 3b `M2/` file (same sha256) → the 3b M3 profile is treated as missing (127 usable profiles).
Profile items are pending (series README updated). `scripts/s1e_4_laminar_prototype.R` reads the profiles in place from the
Evo-M1 root (`EVO_M1_TRAIT_DATA` or the standard macOS path) and stops on any sha256 mismatch.

Scripts and outputs:

- `scripts/s1e_4_laminar_prototype.R` → `tables/s1e/s1e_laminar_strata_panel.csv`, `s1e_laminar_stratum_contrasts.csv`,
  `s1e_laminar_spearman.csv`, `s1e_laminar_profile_rebinned.csv`, `s1e_laminar_profile_vs_tableS2.csv` (+ `_summary.csv`),
  `s1e_laminar_gli_species.csv`, `s1e_laminar_gli_rank_concordance.csv`, `s1e_laminar_gli_stratum_order_by_species.csv`,
  `s1e_laminar_preuss_alignment.csv`.
- `scripts/s1e_4_fig_laminar_prototype.py` → `figs/s1e/laminar_prototype.png` / `.pdf`.
- `metadata/anatomy/s1e_laminar_area_overlap.csv` — the crosswalk (20 rows).

## 2. Area overlap (`s1e_laminar_area_overlap.csv`)

Zilles-2017 areas with a BigBrain Dickscheid partner: **10** (exact/part_of: 44↔Area44, V1↔hOc1, V2d & V2v↔hOc2 (pooled),
V3d↔hOc3d, FG1↔FG1, 4↔4a (agranular); approximate: V3A↔hOc4d, 37B↔FG3+FG4, 10L↔Fp1). With a laminar profile: **7** (3b, 4↔4p,
44↔44d, 45, 47, PFm, V1↔hOc1; 7A has no Zilles partner). With a GLI item: **4** (44, 45 from Table S3; 23, 31 from Zilles
1986). All three laminar sources meet only in **area 44** (Table S2 strata + BigBrain layers + 44d profile + 7-species GLI).
Unmatched: BigBrain 7P and hOc5; profile 7A; GLI areas 29/30.

## 3. Methods

**Strata aggregation.** BigBrain layer densities were aggregated to the three Zilles strata as thickness-weighted means
(Σ density·thickness / Σ thickness). Check: the weighted mean over the six layers reproduces each area's published I–VI density
to within 2.3 % (max |rel. diff| 0.0226, from the same patches). Cumulative layer thicknesses also give the depth fractions used
to re-bin the profiles. Zilles 2017 delimited strata cytoarchitectonically on neighbouring cell-body sections and published
no depth fractions; for agranular areas 4 and 6 they defined layer IV as a 3 % depth stripe below layer III. The BigBrain 4a
layer IV is 45 µm = 1.2 % of depth, so the area-4 granular cell rests on a 45-µm stripe and is flagged (`agranular_flag`).

**Receptor classes** (pharmacology, as in s1e_1/s1e_3): iGlu = AMPA+NMDA+kainate; GABA = GABAA+GABAA/BZ+GABAB;
ACh3 = M1+M2+M3 (α4β2 kept separately); monoamineGPCR = α1+α2+5-HT1A+5-HT2+D1; total_15; EI = iGlu/GABA.
`*_per_celldensity_ratio` = (fmol receptor / mg protein) ÷ (cells / 0.1 mm³): a **ratio of two differently normalised
quantities** (protein mass in 3 unfixed frozen brains vs tissue volume in one fixed, shrunken brain) — a relative laminar
index only, not receptors per cell.

**Panel sets.** `primary`: 6 area keys × 3 strata = 18 cells (V2d+V2v averaged, unweighted, onto hOc2). `unpooled_V2`: 21 cells
(V2d and V2v share the hOc2 cell density; not independent). `sensitivity_incl_approximate`: + V3A, 37B (FG3+FG4
thickness-weighted), 10L = 27 cells. Spearman ρ are descriptive, uncorrected, on 17–27 non-independent cells.

## 4. Results (descriptive)

**Within-area stratum contrasts (primary, `s1e_laminar_stratum_contrasts.csv`).** Every receptor class in every one of the 6
areas is highest in the supragranular stratum and lower in the infragranular stratum (24/24 class×area cells have
granular/supra < 1 and infra/supra < 1; median granular/supra: iGlu 0.86, GABA 0.80, ACh3 0.86, monoamine 0.64, total 0.80;
median infra/supra: iGlu 0.66, GABA 0.49, ACh3 0.60, monoamine 0.52, total 0.53). EI rises with depth in all 6 areas
(infra > gran > supra; median infra/supra 1.36). BigBrain cell-body density peaks in the granular stratum in 5/6 areas
(gran > supra > infra; median granular/supra 1.36, infra/supra 0.94); the exception is agranular area 4 (supra > gran > infra).
Example values (primary panel): V1 cell density supra 78.0 / gran 98.2 / infra 73.7 cells per 0.1 mm³, total receptors
16 466 / 14 263 / 6 687 fmol/mg; area 44: 51.5 / 64.7 / 49.7 and 10 882 / 9 281 / 8 163; human GLI (Table S3, n=4)
14.4 / 15.6 / 13.7 %.

**Receptor class vs cell density across cells (`s1e_laminar_spearman.csv`).** Primary (n = 18): iGlu ρ = 0.55, GABA 0.50,
ACh3 0.44, monoamine 0.43, total 0.51, EI −0.24. Excluding the area-4 granular cell (n = 17): 0.48, 0.46, 0.36, 0.33, 0.47,
−0.30. Sensitivity set (n = 27): 0.40, 0.44, 0.35, 0.25, 0.42, −0.24. The positive ρ is driven by the shared supra+granular
> infragranular pattern (the between-area component — e.g. V1 having both the densest cells and the most receptors — and the
within-area component point the same way); the within-area 3-stratum orders disagree (receptors: supra > gran > infra; cells:
gran > supra > infra), giving a mean per-area ρ of 0.58 for the classes and −0.58 for EI (n = 6 areas). The laminar receptor
pattern therefore does not follow cell-body density within areas; it tracks the neuropil-rich supragranular layers.

**Laminar profiles vs Table S2 (`s1e_laminar_profile_vs_tableS2_summary.csv`).** Expected agreement within rounding was
**not** found. All-layer profile means vs Table S2 all-layer means: median |% diff| 33–55 % per area (3b 38 %, 44d 33 %,
45 34 %, 47 33 %, 4p 34 %, PFm 56 %, hOc1 52 %); 1–7 of 15 receptors within 25 %; receptor rank orders are preserved
(Spearman across receptors 0.74–0.96). Strata re-binned with BigBrain depth fractions: 44d median |% diff| 30–36 %
(ρ 0.85–0.89), hOc1 47–51 % (ρ 0.71–0.79), 4p 20–25 % (ρ 0.88–0.97; the Zilles 3 % rule and the BigBrain 4a layer-IV
boundary give near-identical bins). Some receptors match exactly (44d AMPA 492.8 vs 492, NMDA 967.8 vs 964, α1 125.6 vs 126,
GABAA +4 %), others differ by >100 % (α2 is 2–4.7× higher in the profiles for 5 of 7 areas; kainate +165 % and D1 +118 % in
hOc1). Interpretation: the profiles are single-section "examples", the EBRAINS fingerprints (= Table S2 means; checked for 44d)
are 3-brain means; the α2 discrepancy is systematic and unexplained (ligand/protocol or brain mismatch) — see §7. The profile
re-binning is therefore not a validation of Table S2 strata, and the profile strata should not be substituted for Table S2.

**Cross-species GLI (`s1e_laminar_gli_species.csv`, `s1e_laminar_gli_rank_concordance.csv`).** Areas 44/45, 7 species: the
granular stratum is densest in 4/7 species in both areas; the modal order is tied (44: gran > supra > infra 3, supra > gran >
infra 3, gran > infra > supra 1). Humans have the highest GLI in every stratum of 45 and in supra/granular of 44 (2nd in
infragranular 44, after *Pan troglodytes* 13.68 vs 13.66). Species rank orders between strata: area 44 ρ = 0.71 (supra–gran),
0.75 (supra–infra), 0.93 (gran–infra); area 45 only 0.25, 0.29, 0.68 (n = 7; uncorrected p 0.05–0.59). Posterior cingulate,
17 species: granular stratum densest in 9/17 (23) and 13/17 (31) species; rank concordance 0.92–0.98 (23), 0.93–0.96 (31),
0.81 (29), 0.93 (30) (n = 17). Human cingulate GLI is low in absolute terms (area 31: 10.1 / 12.4 / 10.1 %) against a
17-species median of 14.9 / 24.0 / 18.8 %; the laminar *shape* (granular peak) is conserved, the level is not.

## 5. What layers-and-metabolism can and cannot be said with the data in hand

- There is **no numeric laminar metabolic measure** in the repository. Laminar metabolism is represented only by ordinal
  cytochrome-oxidase (CO) states (Preuss 1999/2002 schema: present/absent/reduced; Wong-Riley 1993 descriptive reactivity by
  layer, with one numeric contrast: IVcα 5–14 % lower optical density than IVcβ in adult human V1). Numeric CO exists only
  regionally and for rat (Hevner et al. 1995 Fig. 31: 11 regions, mU/mg tissue and mU/mg protein) and for macaque V1
  compartments (Hevner & Wong-Riley 1990 Table 1: puff/interpuff in layers 2–3 = 1.32/1.00 relative OD; 4C dark/light columns
  under TTX = 1.00/0.80). None of these is a layer-by-layer profile that could be joined to the strata panel.
- What can be said: the laminar **receptor** gradient (supragranular-heavy, E/I rising with depth) and the laminar **cell-body**
  gradient (granular peak) are both quantified in human for 6 areas; they are not parallel within areas. The CO literature puts
  the most CO-reactive human V1 stratum in layer IVc (Wong-Riley 1993), i.e. with the cell-body peak rather than the receptor
  peak — a qualitative statement only, not testable here beyond ordinal consistency.
- What cannot be said: any ratio of receptors or cells to laminar metabolic rate, any cross-species laminar metabolic
  comparison, and anything about layer 4A from BigBrain (one layer IV).

## 6. Preuss alignment — human V1 layer 4A (`s1e_laminar_preuss_alignment.csv`, 41 rows)

| Source | Observation | Consistency |
|---|---|---|
| Preuss et al. 1999 Fig. 2 / main text | CO-dense band in layer 4A **present** in *Saimiri*, *Ateles*, *Macaca*; **absent** in *Pongo*, *Pan*, *Homo* (Fig. 6: hominoid loss; 4B CO band present in hominoids, absent in monkeys) | — |
| Wong-Riley et al. 1993 text rows 40, 58 | "lamina IVa was again not discernable in CO-reacted sections" (4th postnatal month); "lack of a distinct lamina IVa in CO preparations" obscures the upper boundary of IVb (adult, 0.11–0.16 mm) | consistent with Preuss: no CO-dense 4A in human |
| Dickscheid et al. 2021 hOc1 | layer IV is a single segmented layer: 98.2 cells / 0.1 mm³, 852 µm = 33.1 % of depth, 36 229 cells; IVa/IVb/IVc are **not** subdivided | BigBrain cannot show or test a 4A feature |
| Preuss & Coleman 2002 Fig. 5 / Fig. 7 | Cat-301 mesh in 4A **present only in humans** (absent in *Pan*, *Pongo*, *Macaca*, *Ateles*, *Saimiri*); calbindin: dense irregular territories in humans, uniform moderate density in *Pan*/*Pongo*, label-sparse band in monkeys | the human-specific 4A compartmentation |
| Hevner & Wong-Riley 1990 Table 1 | macaque CO quantified for puff/interpuff (layers 2–3) and 4C columns only; **no 4A row** in the item | macaque CO-dense 4A rests on Preuss 1999 (ordinal) in this repository |

## 7. Unresolved

1. Profile vs Table S2 disagreement (median |% diff| 30–55 %; α2 2–4.7× higher in profiles for 5/7 areas): check the EBRAINS
   dataset descriptions for the profile ligand/brain and whether the "pr_examples" are from one of the three Table S2 brains.
2. The 3b M3 profile is an upload duplicate of M2 on EBRAINS; report to the dataset custodian.
3. Area 4 granular stratum: 45-µm BigBrain layer IV vs Zilles' 3 % stripe — the area-4 granular cell is not comparable to the
   granular cells of granular areas (flagged, and the Spearman is reported with and without it).
4. V2d/V2v → hOc2 pooling is unweighted (no dorsal/ventral split in BigBrain hOc2).
5. Zilles 1986 species values for n = 2–3 brains are not stated to be means (item README).

## 8. Prioritised next data

1. **Macaque laminar receptor profiles** (EBRAINS, Rapan et al. 2021/2023 macaque receptor atlas; profiles per area) — the only
   route to a cross-species laminar receptor comparison matching Table S2 strata.
2. **Wagstyl et al. 2020 BigBrain layer segmentation** — layer thicknesses and staining intensities for all cortex, extractable
   for all 44 Zilles areas (replaces the 12-area Dickscheid coverage; intensity ≠ cell density, so pair with Dickscheid for calibration).
3. **A numeric laminar CO source** (e.g. optical-density profiles of CO histochemistry in human/macaque V1 and association
   cortex; Wong-Riley lab quantitative studies) — without it "layers and metabolism" stays ordinal.
4. **Neuron vs glia split of BigBrain cell counts** (cell-type classification of the Merker-stained detections) — needed before
   the cell-density-normalised receptor ratio can be read as a per-neuron quantity.
5. Register the laminar profile items in Evo-M1 (owner decision on granularity), then replace the in-place read in s1e_4 with
   manifest-synced inputs.
