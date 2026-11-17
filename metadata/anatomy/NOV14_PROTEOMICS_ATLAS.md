# Tüshaus proteomics, Johansen atlas sampling, and BEIP — 14 November 2026

The three requested arms are implemented and executed. Tüshaus is a separate
protein-abundance analysis. Johansen sampling is an exploratory extraction with
explicit unresolved registration and comparability flags. BEIP remains a separate
s2 replication. Existing strict primary analyses retain their anatomical scope.

## Reproduce

From the project root, with the existing R dependencies (`readxl`, `dplyr`,
`tidyr`, `ggplot2`, `ggrepel`) installed:

```sh
# Optional public-source acquisition; verifies SHA-256 of existing files too.
Rscript checks/fetch_nov14_inputs.R
python3 -m venv .venv-atlas
.venv-atlas/bin/python -m pip install -r helpers/requirements-atlas.txt
Rscript scripts/s1d_2a_extract_tushaus_proteomics.R
Rscript scripts/s1d_2b_tushaus_protein_metabolic_analysis.R
Rscript scripts/s1c_sample_johansen_atlas.R
Rscript scripts/s2_beip_replication.R
Rscript checks/test_nov14_proteomics_atlas.R
Rscript checks/test_strict_anatomical_matching.R
Rscript checks/test_nov12_analysis_expansion.R
Rscript scripts/verify_outputs.R
```

These commands assume the existing prepared Heiss rates, primary Johansen table,
and synced BEIP inputs. `scripts/run_all.R` orders those upstream dependencies
and includes all three arms. It uses cached inputs, not network downloads.
`TUSHAUS_EV1` can override the workbook path. `JOHANSEN_PYTHON` can select an
alternative Python interpreter with the listed libraries. The local virtual
environment is ignored by Git. `metadata/nov14_source_downloads.csv` records
download URLs, archive members, retrieval date, and immutable SHA-256 checksums.
Do not regenerate checksums merely to accept a changed upstream download.

## Tüshaus 2023: extraction and anatomy

Source: [Tüshaus et al., EMBO Journal, DOI 10.15252/embj.2023114665](https://link.springer.com/article/10.15252/embj.2023114665),
Dataset EV1, locally `data_raw/Tushaus_2023/Dataset_EV1.xlsx`.
Sheets `SampleID`, `LFQ`, and `iBAQ` supply 38 tissue samples from 13 territories
in **one donor**, a 56-year-old man. There are three samples per territory except
olives (two). SN has published replicate labels 2, 3, 3; the extractor preserves
these labels, flags the duplication, and uses the distinct sample IDs.

LFQ contains 303,170 protein/sample measurements; iBAQ contains 305,816. The
extractor preserves their published log-intensity units. The EV1 headers do not
unambiguously specify the logarithm base for both measures, so it neither takes
another logarithm nor exponentiates values. It performs no new imputation.
The regional enrichment `class` annotation is retained as published; it is
**not a cell class**. These data cannot directly yield Siletti-style cell fractions
or protein mass fractions.

The workflow parallels s1b's extraction → explicit anatomical mapping → regional
summary → metabolic association, with a different biological outcome. The shared
anatomy map is `metadata/anatomy/s1d_tushaus_region_map.csv`; it is checked by the
same strict-policy reader and audit as the other studies.

| Source code | Territory | Metabolic analysis |
|---|---|---|
| CBM | Cerebellar gray matter | Cerebellar cortex; separated from white matter |
| CN | Caudate nucleus | Included |
| HPC | Hippocampus | Included; subfield mixture unspecified |
| NAC | Nucleus accumbens | Included; core/shell mixture unspecified |
| NR | Red nucleus | Included; small-nucleus sampling limitation |
| PUT | Putamen | Included |
| SN | Substantia nigra | Included; pars mixture unspecified |
| CBw | Cerebellar white matter | No corresponding Heiss rate |
| CC | Cerebral cortical gray matter | Unspecified cortical cube is not global cortex |
| CCw | Cerebral white matter | Not established to be centrum semiovale |
| CL | Claustrum | No corresponding Heiss rate |
| OB | Olives | No corresponding Heiss rate |
| THL | Thalamus | Cannot assign unspecified tissue to a Heiss nucleus |

Included matches establish named territory/tissue agreement, not identical
sampling boundaries across donors and studies. All 13 territories remain in the
source summaries and marker-profile figure. Excluded rows receive no primary
metabolic rate.

The screen requires at least two tissue samples per protein/region and complete
coverage of the **same seven regions**: 6,810 LFQ and 6,797 iBAQ protein IDs.
OLS and Spearman results use BH correction separately within each measurement
family. LFQ and iBAQ measure the same material and are not independent replications.
The GFAP, MBP, NEFL, MAPT, UCHL1, SV2A panel also has leave-one-region-out
checks and model diagnostics. Its six tests are corrected separately within
each measure. A centered log-abundance heatmap describes regional profiles,
not composition percentages.

No panel association survives BH correction. The large screen has zero OLS
BH discoveries in either measure and one Spearman BH discovery in LFQ, none
in iBAQ. These are exploratory spatial associations with n = 7 regional points
and one donor; tissue cubes do not increase donor sample size. Region dependence
and uncertain sampling boundaries limit inferential interpretation.

Outputs: two RDS files under `data_intermediate/tushaus_2023_*`, eight CSVs under
`tables/s1d/tushaus_*`, two PNG/PDF figures under `figs/s1d/tushaus_*`, and marker
diagnostics under `figs/s1d/diagnostics/Tushaus_*`.

## Johansen atlas: coverage and limits

Source: [Johansen atlas release](https://xtra.nru.dk/SV2A-atlas/), associated with
[Johansen et al. 2024](https://doi.org/10.1523/JNEUROSCI.1750-23.2024).
The distributed MNI Bmax images are 5-mm-smoothed maps; published regional
Table 2 estimates use a different, nonsmoothed regional procedure. Cortical
volumetric data are explicitly provided for visualization only. The analysis
does not derive new cortical estimates from these volumes.

`scripts/s1c_sample_johansen_atlas.R` invokes
`helpers/sample_johansen_atlas.py`. Affine-aware resampling handles voxel sizes,
origins, and axis flips. Label images use nearest-neighbor resampling; the
CIT168 probability image uses linear interpolation followed by a 0.5 threshold.
Mean Bmax is calculated over the bilateral mask union, weighting voxels equally.
Outputs retain mask labels, volumes, map coverage, left/right means and
centroids, spatial SD, and the mean voxelwise between-subject SD. Neither SD
is the between-subject SD or standard error of a regional mean.

| New territory | Main mask | Small-nucleus flag |
|---|---|---|
| Lateral geniculate nucleus | MASSP2 | Yes |
| Medial geniculate nucleus | MASSP2 | Yes |
| Substantia nigra | MASSP2 | Yes |
| Red nucleus | MASSP2 | Yes |
| Subthalamic nucleus | MASSP2 | Yes |
| Inferior colliculus | MASSP2 | Yes |
| Superior colliculus | MASSP2 | Yes |
| Dentate nucleus | Diedrichsen | Yes |
| Internal capsule | MASSP2 | No |
| Vermis | AAL3 lobular union | No; descriptive only, tissue match unresolved |

Mask sources:

- [MASSP2 release](https://doi.org/10.21942/uva.27291579.v1), MNI2009b labels.
- [AAL3 release](https://www.gin.cnrs.fr/aal-3/), nominal MNI/Colin space.
- [CIT168 subcortical atlas](https://evendim.caltech.edu/subcortical-atlas),
  MNI2009c STN probability mask as an alternative.
- [Diedrichsen cerebellar atlas](https://github.com/DiedrichsenLab/cerebellar_atlases),
  MNI152NLin6AsymC dentate labels, not dentate gyrus.

All ten new masks have complete nonzero map coverage, but **grid resampling
does not establish anatomical registration between MNI template variants**.
Every extraction has `template_registration_validated = FALSE` and
`primary_included = FALSE`. Mask overlays provide QC, not registration proof.
Small-nucleus estimates are especially vulnerable to smoothing, partial-volume
effects, and boundary differences. Six ±1-mm translations quantify location
sensitivity. Alternative masks quantify definition sensitivity for five nuclei:
relative ranges are about 41% for LGN, 68% for MGN, 8% for red nucleus,
6% for SN, and 1% for STN.

**The comparability check fails.** In five overlapping structures, map-derived
means are approximately 30–65% below published Table 2: caudate −60%, putamen
−30%, hippocampus −49%, amygdala −65%, accumbens −38%. This does not by itself
identify the cause; template alignment, sampling/weighting, and smoothing all
require review. Do not rescale these differences away or assume smoothing alone
explains them. The script fits the 12 published regions and nine eligible new
atlas regions separately; it does **not** fit a pooled 21-region model. The
21-row inventory explicitly says `combine_for_inference = FALSE` and the figure
separates the sources. The atlas-only Spearman correlation is 0.333 (p ≈ 0.381),
which remains exploratory and unvalidated.

Global cerebral cortex, medial thalamic nuclei, basal forebrain, and centrum
semiovale remain unsampled because a justified source-matched mask is unavailable.
No whole-thalamus or three-nucleus proxy is substituted. The original 12-region
Johansen primary analysis is unchanged. Promotion of any new estimate requires
verified template correspondence and an explanation of the overlap discrepancy.

The seven atlas CSV/JSON outputs under `tables/s1c/johansen_atlas_*` document
samples, coverage, geometry/checksums, source-separated fits, overlap checks,
and mask sensitivity. Figures include mask overlays and the separate-source
comparison.

Atlas attribution/licensing: NRU, Copenhagen University Hospital, Rigshospitalet;
Johansen et al. 2024. The downloaded atlas and derived Johansen tables/figures
are distributed under [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/),
as specified by the release. Modifications here comprise mask resampling,
regional sampling, sensitivity calculations, and new figures; the original
NIfTI images are unchanged. This notice concerns these derived data, not the
entire repository. Retain each anatomical atlas's own attribution and terms.

## BEIP replication arm for s2

The existing Nov 12 arm `scripts/s2_beip_replication.R` was rerun and validated,
not duplicated. It reads McLaughlin 2014 Tables 1/4 for the childhood EIG–NIG
contrast and Sheridan 2022 Table S3 for adolescent CAUG–NIG and FCG–NIG
contrasts. These are separate from ERA. The strict map retains seven childhood
regions and six adolescent regions per contrast, excluding unmatched whole
thalamus from Heiss-based fits.

The three OLS p-values are approximately 0.580, 0.241, and 0.251; Spearman
p-values are 1.000, 0.156, and 0.257. The matched regional analyses do not provide
clear convergent support. They use published raw means, not covariate-adjusted
effects. Childhood and adolescence belong to the same longitudinal cohort, and
adolescent contrasts share controls; these are not three independent studies.

Outputs remain `tables/s2/beip_replication_*` and `figs/s2/beip_replication.*`.
See `metadata/anatomy/PRIMARY_MATCH_POLICY.md` and the source maps for
the prior policy decisions. No sibling-project exports are changed by this work.
