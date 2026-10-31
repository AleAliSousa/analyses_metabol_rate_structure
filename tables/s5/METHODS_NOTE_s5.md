# Study s5 — Coto-style analyses on the Heiss et al. 2004 human dataset

_Methods note, 2026-10-31. Scripts `scripts/s5a_*.R` … `s5e_*.R`; helper `helpers/s5_helpers.R`;
anatomy map `metadata/anatomy/s5_heiss_region_volume_map.csv`. All numbers below are read from the
tables in this folder._

## Purpose
Reproduce, on humans, the three analyses of Coto et al. 2026 (PNAS; *Tetramorium immigrans*):
(i) scaling of metabolic rate on mass, (ii) relative brain mass vs relative brain metabolic rate,
(iii) preservation of brain energy under deprivation.

## s5a — within-brain scaling of regional glucose use on regional volume
* **Unit of analysis:** Heiss 2004 region (rCMRGlc, µmol/100 g/min, both hemispheres) paired with a
  human volume (mm³, both hemispheres). Pairings follow `tables/Table 1 … .xlsx` and the s4a
  Kochiyama crosswalk; volumes are the `Homo sapiens` row of `data_raw/volumes_wide.csv`
  (Stephan/Zilles-collection lineage) and `data_intermediate/s4a_endocranial_region_budgets.csv`.
* **Non-overlap rules (primary set, 18 regions, 96.5 % of total brain volume 1 257 139 mm³):**
  neocortex grey minus area striata grey minus insula carries the cortex global average (33.5);
  V1 grey stands for the occipital lobe (Table 1 simplification); striatum split 0.45/0.49/0.07
  (de Jong 2012 proportions from the Table 1 note, normalised); mesencephalon carried as one point
  with the Table 1 composite rate 25.3 (weights 0.08/0.11/0.53/0.16); thalamus minus LGN with the
  unweighted mean of the two remaining Heiss thalamic terms (LGN is its own point); cerebellar
  hemispheres (grey + white) carry the cerebellar-cortex rate (same device as s4a); neocortex white
  carries the centrum-semiovale rate.
* **Left out (listed in `s5a_excluded_regions.csv`):** frontal, parietal, temporal lobes (no
  non-overlapping Stephan-basis volume; represented inside the remaining neocortex), basal
  forebrain (no defensible volume), medial geniculate and medial thalamic nuclei and the four
  midbrain nuclei individually (carried only inside composites).
* **Sensitivity set `kochiyama_lobes` (19 regions):** the four Kochiyama grey + white lobes and
  Kochiyama cerebellar cortex / vermis replace the neocortex-grey, V1, cerebellar-hemisphere and vermis
  rows; the Stephan neocortical white-matter row is dropped to avoid double counting.
* **Conversion:** mass = volume × 1 mg/mm³ (task convention; s4a used 1.036 g/cc, +3.6 %);
  glucose (µmol/min) = rCMRGlc × mass(g)/100; energy = × 2.87 J/µmol / 60 → W.
* **Fits:** `smatr::sma` (slope test vs 1) and OLS of log10(glucose use) on log10(volume); also
  log10(rCMRGlc) on log10(volume) (Karbowski-style mass-specific exponent; OLS is the meaningful
  estimator there because r² ≈ 0 makes the SMA slope a pure SD ratio).
* **Results:** all-region SMA slope **1.002 [0.935, 1.074]**, OLS 0.993 [0.924, 1.063], r² 0.98,
  n = 18 — isometric. rCMRGlc on volume: OLS slope −0.007 [−0.076, 0.063], p = 0.85. Excluding white
  matter: SMA 1.055 [1.009, 1.103] (n = 16); grey-only volumes: 1.062 [1.006, 1.121] (n = 13);
  Kochiyama set: 1.061 [1.027, 1.096] (n = 19). Coto's cross-species slope 0.87 and Karbowski's
  cross-species regional exponents (median −0.15, range −0.24 to 0.05) are *between-species* exponents
  on brain mass; within one human brain there is no systematic size dependence of rCMRGlc, so the
  budget scales with volume with slope ≈ 1 (slightly > 1 once the low-rate white matter is removed,
  because the large grey structures — cortex, putamen — have above-median rates).

## s5b — relative brain mass and relative brain metabolic rate
* Human whole-brain rCMRGlc = volume-weighted mean over the s5a primary set = **25.4 µmol/100 g/min**
  (non-white-matter regions only 32.5; cortex global average 33.5; Clarke & Sokoloff 1994 whole brain
  31.0; Evo-M1 merge 28.9). The low whole-brain value is driven by white matter (35 % of volume at
  12.3). Scaled to brain mass 1300.4 g (Evo-M1 brain_mass merge): 330.6 µmol/min = **15.8 W**
  (Clarke & Sokoloff route 403 µmol/min = 19.3 W; Karbowski S1 O₂ route 48.6 ml/min = **16.3 W**).
* Body MR = BMR 12 925 ml O₂/h × 20.1 J/ml = 72.2 W; body mass 61 535 g (Evo-M1 body_ecology).
* **Human: relative brain mass 0.0211; relative brain MR 0.219 (Heiss glucose), 0.226 (O₂),
  0.267 (Clarke & Sokoloff glucose).**
* Mouse (Evo-M1 whole-brain CMRgl 74.8, brain 0.534 g, BMR 50.3 ml O₂/h): relative brain mass
  0.0184, relative brain MR 0.068 (the shared species table gives 0.055 because it uses the
  Karbowski S2 total 0.32 µmol/min instead of 74.8 × mass = 0.40 µmol/min).
* Ants (Coto Dataset S1 treatment means, 31 °C): relative brain mass 0.042 / 0.054 / 0.052
  (control / food-deprived / queenright); relative brain MR as *ratio of treatment means*
  0.211 / 0.367 / 0.276, vs the per-worker means carried in the species table (0.250 / 0.450 /
  0.309) which match the published 24 % / 45 % / 30 %. Both forms are in the table.
* Brain mass-specific O₂: human 0.0374 ml/min/g (Coto quotes ~0.035); ant queenright 0.0593 at
  31 °C → **0.0899 at 37 °C with Q10 = 2** (Coto: ~0.102) → **2.4× human** (Coto: ~3×). Direction
  reproduced, magnitude lower; not adjusted.

## s5c — budget preservation under deprivation (Mackes 2020 adoptees)
* Inputs are the s2 tables: per-hemisphere volumes (cm³) of 8 regions with a Heiss rate (cerebral
  cortex, thalamus [s2 three-term composite], caudate, putamen, pallidum, hippocampus, amygdala,
  accumbens); normative = Walhovd 2011 (n = 262), English adoptees n = 21, Romanian adoptees n = 67.
* Budget = rCMRGlc × volume (g at 1 g/cm³)/100, summed over the 8 regions. % change =
  (adoptee − normative)/normative (sign flipped relative to s2, which reports normative − adoptee).
  Parametric bootstrap of each group mean from mean, SD/√n (2000 draws, seed 20261031).
* **English adoptees:** volume −2.84 % [−4.33, −1.32]; budget −2.94 % [−4.44, −1.40];
  budget − volume −0.10 pp [−0.14, −0.06]. **Romanian adoptees:** volume −11.52 % [−12.81, −10.23];
  budget −11.59 % [−12.89, −10.30]; budget − volume −0.07 pp [−0.10, −0.05]. Volume-weighted mean
  rCMRGlc falls by 0.1 %. There is **no sparing of expensive regions**: the summed budget falls
  marginally *more* than volume, and the s2 regression (standardised form) has a negative slope of
  % volume change on rCMRGlc: English −1.02 [−1.88, −0.16] %/unit, R² 0.59, p = 0.027;
  Romanian −0.88 [−1.56, −0.20], R² 0.62, p = 0.020; all adoptees −0.95, R² 0.62, p = 0.020
  (n = 8 regions; Spearman ρ −0.64 to −0.69, p 0.06–0.09).
* The totals are dominated by cerebral cortex (259 of 284 cm³), so the summed budget essentially
  tracks the cortical volume change.
* **Caveat (stressor):** early institutional deprivation measured as adult volume ≠ acute food
  deprivation with a metabolic-rate outcome (ants) ≠ perinatal protein restriction (mice). Only the
  direction of budget preservation is comparable.

## s5d — cross-species scaling (Coto Fig. 7)
* (a) Coto's own table (`coto_fig7_points.csv`, 7 mammals + ant Q10-corrected): SMA slope mammals
  only **0.875 [0.807, 0.948]**, r² 0.995; with ant **0.864 [0.843, 0.886]**, r² 0.999 — matches the
  published 0.87 ± 0.07, r² 0.99.
* (b) Evo-M1 species table: O₂ species only (Karbowski S1, n = 7) 0.830 [0.670, 1.028]; all mammals
  with any route (n = 10; O₂ where available, else glucose µW) **0.819 [0.690, 0.974]**, r² 0.955;
  adding the Evo-M1 queenright ant (31 °C, not Q10-corrected) 0.967 [0.877, 1.067].
* Prediction bands (mammal-only OLS): human within (observed/predicted 1.09 Coto table; 1.17 Evo-M1;
  1.29 leave-one-out); mouse within (0.56; 0.39 leave-one-out); Coto's Q10-corrected ant within
  (1.15), but the Evo-M1 ant row (2.69 × 10⁻⁶ ml O₂/min) is **outside/below** (0.06×). Coto's ant
  point (1.53 × 10⁻⁵ ml/min) is 5.7× the Evo-M1 Dataset-S1 queenright mean at 31 °C (3.8× after Q10 = 2 correction); the
  discrepancy is inherited from the inputs and is not resolved here.

## Deviations from the source paper / task
* Volumes are literature means (not the Heiss subjects); no subject-level uncertainty for s5a.
* Tissue density 1 mg/mm³ (task) rather than s4a's 1.036 g/cc.
* The thalamus-minus-LGN rate rule is new (documented in the map); all other pairings reuse Table 1 /
  s2 / s4a.
* `human_comparison_rows.csv` carries no "all adoptees" bootstrap (no pooled SD in the s2 inputs);
  the regression row for all adoptees uses s2's `norm_all` column.
