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
| s1b transcriptomics | `s1b_linnarsson_roi_map.csv` | Fine Linnarsson ROIs are explicitly grouped into Heiss regions. Source-specific choices such as FI/Idg/Ig -> insular lobe and TH-TL -> temporal lobe are recorded here. |
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

## Study-map schema

Every study map uses these columns:

- `source_region`: the exact term in the comparison dataset.
- `analysis_anatomy_id`: the registry ID used to group and join the analysis.
- `heiss_anatomy_id`: the native Heiss measurement supplying the rate.
- `component_weight`: a positive weight for a direct or composite Heiss rate.
- `relationship`: a short controlled description such as `exact`, `proxy`, or
  `composite_unweighted_mean`.
- `mapping_note`: the anatomical rationale.

Maps may add study-specific columns. For example, s4a uses
`allocation_weight` to divide a Kochiyama source parcel between analysis
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
