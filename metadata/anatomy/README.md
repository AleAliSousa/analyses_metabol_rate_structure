# Anatomical mapping architecture

The project uses `anatomy_id` as its stable anatomical key. Source labels and
plot labels are retained for provenance and display. Persistent mapping
decisions use IDs; `join_heiss_rates()` is the one compatibility boundary for
analysis tables that have already aggregated rows under a display label.

## Authoritative layers

1. `data_raw/Heiss_etal_2004_TABLE1.csv` is the authoritative local Heiss
   measurement snapshot. The pipeline refreshes it from the matching source in
   Evo-M1-Trait-Data; it is never edited or supplemented with calculated rows.
2. `anatomy_registry.csv` defines stable project anatomy IDs, hierarchy,
   preferred labels, plotting colours, and display order.
3. `heiss_2004_region_map.csv` is the single primary map for all 26 native
   Heiss rows. It maps each published source label to one registry ID.
4. A comparison paper gets one additional map only when its source terms do
   not map directly through the primary Heiss map.

`scripts/00_prepare_heiss_rates.R` validates these layers and writes
`data_intermediate/heiss_2004_regions.csv`. The same preparation step validates
and stages the already separated s1a and s2 comparison inputs as:

- `data_intermediate/s1a_stereology_comparison.csv`
- `data_intermediate/s2_stress_volume_comparison.csv`

The one-time split is documented by the archived
`scripts/archive/prepare_legacy_stereology_stress_tables.R`, which produced
`data_raw/stereology_comparison.csv` and
`data_raw/stress_volume_comparison.csv`. The `rCMRGlc` columns still present in
the old mixed files are legacy audit values only. The archived migration checks
them against the central Heiss source (including the documented s2 composite)
before discarding them. No active preparation or analysis reads those embedded
rates. The clean comparison inputs retain source labels but deliberately do not
store `anatomy_id`; the active preparation attaches IDs from the authoritative
central or study map.

## Active comparison maps

| Study | Authoritative map | Reason for a study layer |
|---|---|---|
| s1a stereology | Primary Heiss map only | All retained comparison rows use native Heiss labels. |
| s1b transcriptomics | `s1b_linnarsson_roi_map.csv` | Fine Linnarsson ROIs are explicitly grouped into Heiss regions. Source-specific choices such as FI/Idg/Ig -> insular lobe and TH-TL -> temporal lobe are recorded here. The same rows carry the two explicit cortical-scope flags described below. |
| s2 stress volume | `s2_stress_volume_region_map.csv` | English synonyms are resolved and the whole-thalamus rate is explicitly defined as the unweighted mean of the three Heiss thalamic measurements, rounded to the source precision (26.6). |
| s1c Johansen synaptic density | `s1c_johansen_region_map.csv` | Johansen labels are resolved; the whole-cerebellum -> cerebellar-cortex and aggregate-white-matter -> centrum-semiovale comparisons are marked as proxies. |
| s4a Kochiyama | `s4a_kochiyama_region_map.csv` | Kochiyama parcels are allocated to native Heiss lobes/tissues, including the 50/50 sensorimotor split. |

s1c deliberately excludes Johansen's whole-thalamus measurement because Heiss
reports only selected thalamic nuclei; it also excludes source regions without
a defensible Heiss counterpart. Those exclusions are analysis scope, not
alternative mappings.

The local s1c input, `data_raw/Johansen_etal_2024_Table2.csv`, is an exact copy
of `Evo-M1-Trait-Data/Johansen_etal_2024/Johansen_etal_2024_Table2.csv` and
retains the paper's uppercase `SV2A_*` source fields. The s1c script translates
the required fields to `synaptic_density_mean` and `synaptic_density_sd` at the
import boundary; lowercase `sv2a` is not used as an analysis variable or output
name.

The s1b mapping-validation script reads
`s1b_linnarsson_roi_map.csv` directly. The former generated
`data_intermediate/rcmr_roi_relationship.csv` is no longer part of the mapping
system, so there is no second editable s1b map.

### s1b anatomical scopes

The s1b map distinguishes three nested analysis scopes instead of inferring
an ambiguous bare "cortex" category from free-text dissection descriptions.
The documentation and active output fields always say `neocortex`, `cerebral
cortex`, or `cerebellar cortex` explicitly.

- `is_telencephalon` is the broadest scope. It includes the cerebral-cortex
  scope plus basal forebrain, caudate, nucleus accumbens, pallidum, and
  putamen. The default s1b_3 scripts analyze all mapped brain regions; their
  telencephalon counterparts use this flag for the planned sensitivity split.
  The telencephalon analyses in s1b_4 and s1b_6 consume the same flag.
- `is_cerebral_cortex` is the project's broad cerebral-cortex scope. It
  includes neocortical ROIs plus the ROIs mapped to `Corpus amygdaloideum` and
  `Hippocampus`.
- `is_neocortex` is the narrower scope used for the Jorstad-like E:I
  comparison. It includes ROIs whose source `ROIGroupCoarse` is `Cerebral
  cortex` and excludes the amygdaloid complex and hippocampus.

Human A25 and A32 are explicitly retained in the neocortex scope. Their
agranular cytoarchitecture does not make them hippocampal or amygdaloid cortex,
and the s1b validation script protects this decision from accidental removal.
The Jorstad-like analysis is a cell-subclass crosswalk applied across the
available s1b neocortical ROIs grouped to Heiss regions; it is not a claim that
the dataset exactly replicates [Jorstad et al.'s eight-area neocortical
sampling design](https://doi.org/10.1126/science.adf6812). Retaining anterior
cingulate cortex follows the interpretation of ACC as a neocortical
specialization described by [Allman et al.
2001](https://doi.org/10.1111/j.1749-6632.2001.tb03476.x).

The E:I scope comparison keeps anatomical and cell-class provenance separate:

Anatomical membership is defined here by `s1b_linnarsson_roi_map.csv`. The
independent Jorstad-to-Siletti cell taxonomy is defined in
`metadata/cell_taxonomy/s1b_jorstad_to_siletti_cell_class_map.csv`; anatomical
scope and cell-class membership must not be inferred from one another.

- **Jorstad-aligned neocortex:** the neocortical scope and neocortical-class
  crosswalk are used together.
- **Siletti-derived cerebral cortex:** the expanded cerebral-cortex scope is
  paired with the broad Siletti neuronal classes, retaining amygdaloid and
  hippocampal excitatory classes.
- **Cross-study extrapolation:** a separately labelled diagnostic applies the
  neocortical-class crosswalk to the expanded cerebral-cortex scope. This
  reproduces the earlier calculation, but it is attributed to neither
  Jorstad's anatomical design nor the native Siletti broad-class definition.

The three directly comparable analyses are kept in separate, explicitly named
scripts:

- `scripts/s1b_5a_neocortical_cell_classes_x_neocortical_regions.R` applies
  the neocortical cell-class crosswalk only to neocortical regions.
- `scripts/s1b_5b_neocortical_cell_classes_x_cerebral_cortex_regions.R`
  applies that same crosswalk to the broader cerebral-cortex scope, including
  the amygdaloid complex and hippocampus. Its output is therefore labelled as
  an anatomical-scope extrapolation rather than a Jorstad replication.
- `scripts/s1b_5c_cerebral_cortex_cell_classes_x_cerebral_cortex_regions.R`
  matches the broader cerebral-cortex scope to the native broad Siletti cell
  classes, including amygdalar and hippocampal excitatory classes. MSN-excluded
  E:I is primary; the same script produces a clearly labelled MSN-included
  sensitivity figure.

All three scripts use `helpers/s1b_EI.R`, which is the single source
for the cell-class definitions, donor-balanced calculations, statistics, and
plot construction. The generated
`data_analysis/s1b_5_EI_cell_class_x_anatomical_scope_definition_table.csv`
places cell-class scope, anatomical scope, MSN treatment, linear-model results,
multiple-testing results, and Spearman results in separate columns. It is
generated by
`scripts/s1b_5d_EI_cell_class_x_anatomical_scope_definition_table.R`.

## Study-map schema

Every study map uses these columns:

- `source_region`: the exact term in the comparison dataset.
- `analysis_anatomy_id`: the registry ID used to group and join the analysis.
- `heiss_anatomy_id`: the native Heiss measurement supplying the rate.
- `component_weight`: a positive weight for a direct or composite Heiss rate.
- `relationship`: a short controlled description such as `exact`, `proxy`, or
  `composite_unweighted_mean`.
- `mapping_note`: the anatomical rationale.

Maps may add study-specific columns. s1b uses `is_telencephalon`,
`is_cerebral_cortex`, and `is_neocortex` for its nested analysis scopes. s4a
uses `allocation_weight` to divide a Kochiyama source parcel between analysis
regions. Multiple rows are permitted only when the file explicitly defines a
composite or split.

## Future s3

s3 is not active and has no current study map. Reintroduction should add one
file using the same schema. If the Heiss occipital-lobe rate is paired with V1,
that row must use distinct registry IDs for `analysis_anatomy_id` and
`heiss_anatomy_id`, set `relationship` to `proxy`, and explain the approximation
in `mapping_note`. Plotting helpers intentionally contain no hidden
occipital-lobe/V1 alias.

## Editing workflow

1. To change a preferred term or colour, edit its registry row without changing
   `anatomy_id`.
2. To add a genuinely new anatomical concept, add one unique `anatomy_id` to
   the registry, then reference that ID from a source map.
3. To add or rename a Heiss source row, update `heiss_2004_region_map.csv` and
   run `Rscript scripts/00_prepare_heiss_rates.R`.
4. Put approximations, parent-child substitutions, composite definitions, and
   weights in the one map belonging to the relevant comparison paper. They do
   not belong in the registry, plotting helpers, or the pure Heiss map.
5. Run the preparation step and the relevant study script. Active mapping
   validation fails on unknown IDs, duplicate assignments, missing source rows,
   or invalid weights. The archived migration additionally checks the old
   embedded rates before reproducing the clean comparison inputs.

The mixed-table extraction is already archived and is not part of `run_all.R`.
Publication analyses use the clean comparison inputs, the central Heiss table,
and the version-controlled study maps.
