# s1e_0 prototype: Julich-Brain v3.1 parcellation of cortical metabolic maps — 15 November 2026

Status: **prototype**. `scripts/s1e_0_extract_metabolic_maps_julich.R` is not in
`scripts/run_all.R` and its outputs are not in `metadata/output_manifest.csv`.
The purpose is to supply per-area metabolic values for the 44 receptor areas of
Zilles & Palomero-Gallagher (2017, Front. Neuroanat., receptor fingerprints) by
way of Julich-Brain v3.1, and to check that the parcellated group maps reproduce
the Heiss et al. (2004) lobe rCMRglc ordering before any area-level analysis is run.

## Reproduce

```sh
# Python with numpy, scipy, pandas, nibabel, neuromaps (templates are fetched to
# data_raw/neuromaps_metabolic/templates on first run; set NEUROMAPS_DATA to move them)
S1E_PYTHON=/path/to/python Rscript scripts/s1e_0_extract_metabolic_maps_julich.R
# S1E_SKIP_PYTHON=1 reuses data_intermediate/s1e_julich_area_metabolism.csv
```

Inputs (checksums: `data_raw/neuromaps_metabolic/neuromaps_download_manifest.csv`
and `data_intermediate/s1e_julich_area_metabolism_manifest.csv`):

| File | Content |
|---|---|
| `data_raw/neuromaps_metabolic/source-raichle_desc-{cmrglc,cmr02,cbf,cbv}_space-fsLR_den-164k_hemi-{L,R}_feature.func.gii` | Vaishnavi et al. 2010 (PNAS) group maps, as distributed by neuromaps 0.0.7 |
| `data_raw/neuromaps_metabolic/source-castrillon2023_desc-cmrglc_space-MNI152_res-3mm_feature.nii.gz` | Castrillon et al. 2023 (Sci. Adv.) group CMRglc, MNI152 3 mm |
| `data_raw/anatomical_atlases/JulichBrain_v31/JulichBrain_v31_fsaverage_labelled_{left,right}_hemisphere.label.gii` | Julich-Brain v3.1 maximum-probability map on fsaverage-164k; labels 1–146 cortical areas and GapMaps, 147 CorpusCallosumMask, 148 BrainStemMask; every vertex carries a label |
| `data_raw/anatomical_atlases/JulichBrain_v31/JulichBrain_v31_MNI152_labelled_{merged,left_hemisphere,right_hemisphere}.nii.gz` + `_lut.csv` | Julich-Brain v3.1 MPM in ICBM 2009c asymmetric space, 1 mm; labels 1–207 shared across hemispheres |
| `data_intermediate/heiss_2004_regions.csv` | Heiss et al. 2004 rCMRglc (µmol/100 g/min); lobe rows used here |

The GIfTI label table is treated as authoritative. The right-hemisphere table
writes label 145 as `Area MFG5 (MFG,IFG)` (no space); the helper normalises
comma spacing so both hemispheres use the exported CSV spelling `Area MFG5 (MFG, IFG)`.

## Methods

### Crosswalk (`metadata/anatomy/s1e_zilles2017_julich_area_map.csv`)

One row per Zilles-2017 area: `zilles_area`, `lobe` (frontal / parietal /
temporal / occipital / cingulate; used only for the Heiss lobe check),
`heiss_lobe_anatomy_id`, `julich_regions` (`;`-separated GIfTI names, all
validated against the label table), `n_julich_regions`, `match_relationship`,
`note`. Relationship codes:

* `exact` — same cytoarchitectonic entity (Jülich-lab nomenclature, e.g. V3v = hOc3v).
* `part_of` — the Zilles area is a part of one Julich area (V2d, V2v ⊂ hOc2; no vertex-level split attempted, so both receive the hOc2 value).
* `union` — the Zilles area is represented by the union of Julich subdivisions (4 = 4a+4p; 41 = Te1.0+1.1+1.2; 6, 8, 24, 32 — see notes for incomplete coverage).
* `approximate` — a Julich area occupies the same territory but is not defined as the Brodmann area (V3A≈hOc4d, 42≈Te2, 22≈Te3, 10L/10M≈Fp1/Fp2, 11≈Fo1–3, 47≈Fo4–7, 37B≈FG3+FG4).
* `none` — no Julich v3.1 counterpart; the territory is a GapMap or absent (20, 21, 36, 38, 37L, 37M, 9, 46, 23, 31).

### Coverage of the 44 areas

| match_relationship | n areas | areas |
|---|---|---|
| exact | 18 | 1, 2, 3a, 3b, V1, V3d, V3v, V4v, FG1, FG2, 5L, 5M, PGa, PGp, PFt, PFm, 44, 45 |
| part_of | 2 | V2d, V2v |
| union | 6 | 41, 24, 32, 4, 6, 8 |
| approximate | 8 | V3A, 37B, 42, 22, 10L, 10M, 11, 47 |
| none | 10 | 37L, 37M, 20, 21, 36, 38, 23, 31, 9, 46 |

26 of 44 areas have a strict (exact/part_of/union) match; 34 of 44 have any
match. Of the 10 unmatched, 8 are temporal (20, 21, 36, 38, 37L, 37M) or
posterior cingulate (23, 31), so the temporal lobe is thinly represented in
the Zilles-based lobe aggregates (3 strict areas: FG1, FG2, 41).

### fsLR-164k → fsaverage-164k

neuromaps' `fslr_to_fsaverage` requires Connectome Workbench (`wb_command`),
which is not installable on this macOS arm64 sandbox (no conda-forge build).
`helpers/parcellate_metabolic_maps_julich.py` therefore uses `wb_command` when
it is on PATH and otherwise a pure-Python barycentric resampler on the same
neuromaps template spheres (`tpl-fsLR_space-fsaverage_den-164k_hemi-*_sphere`
→ `tpl-fsaverage_den-164k_hemi-*_sphere`): each fsaverage vertex is located in a
triangle of the fsLR sphere deformed to fsaverage (containing triangle found
for 163,842/163,842 vertices in both hemispheres; no nearest-neighbour
fallback was needed), values are interpolated with barycentric weights,
source vertices outside the fsLR no-medial-wall ROI get zero weight, and the
fsaverage no-medial-wall mask is applied. For equal densities this is the
barycentric step of Workbench's `ADAP_BARY_AREA` without the vertex-area
adaptation. The manifest column `resample_method` records which path ran
(`python_barycentric` in the shipped outputs). Alignment check: after
resampling, 100 % of vertices labelled CorpusCallosumMask (147) and
BrainStemMask (148) carry zero data; 5.2 % of vertices in cortical labels are
zero, concentrated in hippocampal subfields, area 25, area 33 and
GapMap-Frontal-to-Temporal-II (medial-wall borders).

### Parcellation

Surface: per hemisphere, mean (and SD) over vertices of each Julich label with
finite, non-zero data; `n_vertices`, `n_vertices_valid` and
`frac_vertices_excluded` are reported; labels 0, 147, 148 excluded. The
bilateral row (`hemisphere = B`) is the valid-vertex-weighted mean of L and R.
Volume: the three Julich MNI label volumes are resampled to the Castrillon
3 mm grid by nearest neighbour through the two affines (target voxel → world
mm → source voxel; no re-registration); mean over voxels with finite non-zero
data; `n_voxels`, `n_voxels_valid` reported. Relative columns
(`*_rel_cortex`) divide by the vertex- or voxel-weighted mean over labels 1–146
of the same hemisphere set.

### Glycolytic-index proxy

`gi_proxy` is the vertex-wise residual of cmrglc regressed on cmr02 (OLS over
296,832 valid vertices, both hemispheres pooled; fitted
`cmrglc = −77.07 + 1.2582·cmr02`, r = 0.847). Positive values mean more
glucose use than expected from oxygen use. Vaishnavi et al. 2010 defined GI on
quantile-normalised volumetric images of individual subjects; this proxy uses
the distributed group surface maps in neuromaps' units and is **not** the
published GI. It is parcellated like the other measures.

### Units

Cortical means (Julich labels 1–146, bilateral, vertex-weighted) are 6463
(cmrglc), 5198 (cmr02), 5995 (cbf), 4770 (cbv) in the units neuromaps
distributes; these are not µmol/100 g/min (Heiss 2004 cortex: 33.5) and no
documented scale factor accompanies the files, so the Raichle maps are treated
as normalised/arbitrary and used only as relative values. The Castrillon
cortical mean is 25.87, of the order of absolute µmol/100 g/min (lower than
Heiss's 33.5, as expected for a 3 mm group map with partial-volume dilution and
a different cohort/quantification).

## Heiss 2004 lobe check (`tables/s1e/s1e_lobe_check_vs_heiss*.csv`)

Three aggregation bases, all vertex-weighted (voxel-weighted for Castrillon):
`zilles44_strict_matches` (exact/part_of/union areas), `zilles44_incl_approximate`,
and `all_julich_by_gyrus` (every Julich cortical area assigned a lobe from the
gyrus token in its name; allocortex, cingulate, POS and cross-lobe GapMaps
excluded — `tables/s1e/s1e_julich_region_lobe_assignment.csv`). Raichle values
are rescaled to the Heiss cortex mean (lobe/cortex × 33.5) only for display.

Raichle cmrglc, rescaled, versus Heiss (µmol/100 g/min):

| Lobe | Heiss 2004 | all_julich_by_gyrus (n areas) | zilles44_incl_approximate | zilles44_strict |
|---|---|---|---|---|
| Frontal | 35.3 | 34.60 (50) | 33.90 (9) | 33.94 (5) |
| Parietal | 35.8 | 33.51 (30) | 32.77 (10) | 32.45 (10) |
| Temporal | 30.5 | 30.27 (21) | 31.23 (6) | 32.37 (3) |
| Occipital | 35.8 | 34.47 (9) | 34.58 (7) | 34.29 (6) |
| Insular | 30.3 | 29.82 (16) | — | — |

Castrillon CMRglc (µmol/100 g/min, voxel-weighted) for the same bases: frontal
26.74 / 26.03 / 26.87, parietal 25.00 / 24.31 / 24.31, temporal 26.40 / 26.74 /
28.86, occipital 25.12 / 25.74 / 25.92, insular 27.38 / — / —.

Pearson r with Heiss across lobes (descriptive, n = 4–5 lobes):

| Measure | all_julich_by_gyrus (5 lobes) | zilles44_incl_approximate (4) | zilles44_strict (4) |
|---|---|---|---|
| cmrglc (Raichle) | 0.97 | 0.85 | 0.57 |
| cmr02 | 0.80 | 0.32 | 0.00 |
| cbf | 0.62 | 0.70 | 0.37 |
| cbv | 0.10 | 0.09 | −0.10 |
| gi_proxy | 0.57 | 0.46 | 0.39 |
| castrillon_cmrglc | −0.72 | −0.72 | −0.87 |

Reading: the surface-based Raichle CMRglc reproduces the Heiss lobe pattern
(temporal and insular ≈ 30, the other three ≈ 34–35) when all Julich areas are
used; with only the Zilles-matched areas the parietal lobe (dominated by
somatosensory and IPL areas) comes out low and the temporal lobe rests on
FG1, FG2 and 41 (primary auditory cortex pulls the strict temporal value up to
32.4, above Heiss's 30.5), so r drops. The Castrillon volume parcellation does **not**
reproduce the Heiss ordering (insula and temporal highest, parietal lowest);
with 3 mm voxels and 1 mm labels nearest-neighbour-resampled across two MNI
template variants, thick insular/temporal cortex is less partial-volume
diluted than thin occipital/parietal cortex, and the 44-area subset gives
only 418 valid voxels for the temporal lobe. The Castrillon column should be
regarded as a secondary, volume-space estimate until a surface-projected or
partial-volume-corrected version is used.

## Flags

* Bilateral regions with < 200 vertices (of 146): CA3 (154), TrS (178),
  Ig3 (184), IFS2 (188). None of the 44-area matches falls below 200 vertices
  bilaterally; per hemisphere, Te 1.0 left has only 31 vertices and Te 1.2
  151/126 (41 is therefore carried mainly by Te 1.1 and the right hemisphere),
  Te 2.1 right 190, 7M left 97, Fo4 left 171, 6mp left 198.
* Vertex exclusion > 10 % among crosswalk regions: p24ab 24 %, s24 15 %
  (medial-wall border); the area-24 union has 10.5 % excluded. Hippocampal
  subfields, area 25 and area 33 are mostly or wholly masked and are not used.
* Castrillon voxel counts < 20 per hemisphere occur in 42 region×hemisphere
  rows (small insular, IFS, Heschl, hippocampal labels); `few_voxels_flag` in
  `tables/s1e/s1e_zilles_area_metabolism.csv` marks the 44-area table (none
  bilaterally).
* Left–right difference in cmrglc across the 34 matched areas: median −0.3 %
  of the bilateral mean, range −7.9 % to +4.3 %.

## Limitations

1. Template correspondence is unvalidated: Julich-Brain fsaverage labels are
   the maximum-probability projection of a volumetric MPM; the Raichle maps
   were projected to fsLR by the neuromaps authors from volumetric PET; and
   the MNI step mixes an ICBM 2009c-asym label volume with a 3 mm grid in the
   FSL MNI152 bounding box. Only medial-wall agreement was checked.
2. Group maps from different cohorts (Vaishnavi et al. 2010 young adults;
   Castrillon et al. 2023) with unknown smoothing; SD columns are spatial SDs
   across vertices/voxels, not standard errors.
3. Raichle units are arbitrary/normalised; only ratios are interpretable.
   The `gi_proxy` is a residual of group maps, not the published GI.
4. Hemisphere handling: bilateral rows are vertex-weighted L/R means; the
   crosswalk does not separate hemispheres; Zilles 2017 fingerprints are
   mostly from one hemisphere per brain.
5. The crosswalk is a nomenclature mapping, not a boundary comparison; `union`
   rows for 6, 8, 24, 32 do not cover the full Brodmann area, and
   `approximate` rows should be used only in sensitivity analyses.
6. Lobe assignment of FG1/FG2 (ventral occipitotemporal) to temporal and of
   cingulate areas to no Heiss lobe are choices recorded in the crosswalk.
7. Pure-Python resampling replaces Workbench `ADAP_BARY_AREA`; differences are
   expected to be far below the regional averaging scale but were not
   quantified against a Workbench run.

## Outputs

* `metadata/anatomy/s1e_zilles2017_julich_area_map.csv` — crosswalk (44 rows).
* `data_intermediate/s1e_julich_area_metabolism.csv` — 146 Julich regions × {L, R, B}
  = 438 rows: cmrglc, cmr02, cbf, cbv, gi_proxy (+ SD, relative-to-cortex),
  castrillon_cmrglc, n_vertices, n_vertices_valid, frac_vertices_excluded,
  n_voxels, n_voxels_valid, frac_voxels_excluded, resample_method.
* `data_intermediate/s1e_julich_area_metabolism_manifest.csv` — input/output paths,
  bytes, SHA-256, software versions, method and unit notes.
* `data_intermediate/s1e_fsaverage164k_metabolic_maps.npz` — resampled vertex data
  and labels (regenerable; not a primary output).
* `tables/s1e/s1e_zilles_area_metabolism.csv` (34 matched areas × {L, R, B}),
  `s1e_zilles_area_unmatched.csv` (10), `s1e_crosswalk_coverage.csv`,
  `s1e_julich_region_lobe_assignment.csv`, `s1e_lobe_check_vs_heiss.csv`,
  `s1e_lobe_check_vs_heiss_summary.csv`.
