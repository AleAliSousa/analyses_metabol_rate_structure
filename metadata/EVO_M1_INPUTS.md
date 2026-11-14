# Evo-M1 input synchronization

The publication project keeps version-controlled snapshots of the Evo-M1 data
it consumes under `data_raw/`. At the beginning of `scripts/run_all.R`,
`scripts/00_sync_evo_m1_inputs.R` refreshes those snapshots from the external
Evo-M1-Trait-Data repository.

Set its root once before running the pipeline:

```sh
export EVO_M1_TRAIT_DATA=/path/to/Evo-M1-Trait-Data
Rscript scripts/run_all.R
```

On macOS, the synchronizer also recognizes the existing standard OneDrive
location when the environment variable is unset. On a machine without the
external repository, it uses the committed `data_raw/` snapshots, keeping these
dependencies available offline. It stops if neither an upstream file nor its
local snapshot is available.

All source-to-local relationships live in
`metadata/evo_m1_input_manifest.csv`. Add a row there when an analysis begins
using another Evo-M1 file. Do not add paper-specific copy commands to analysis
scripts.

The synchronizer validates that every declared upstream file exists before it
changes any local file. It copies only files whose checksums differ and writes
the resulting checksums to `data_intermediate/evo_m1_input_provenance.csv`.
That provenance file contains no machine-specific absolute paths or transient
sync status, so its contents remain stable across online and offline runs.

The manifest deliberately points the Boyer and Harrington inputs to
`__merging_fossil_brain_glucose/inputs/`. Those are the normalized schemas
consumed by s4b. The Kochiyama, Johansen, and Heiss inputs come directly from
their named source directories, while `volumes_wide.csv` comes from the current
comparative-volume merge.

`data_raw/species.nwk` is not in the refresh manifest. It is the fixed,
project-specific 59-species primate tree used to delimit the taxonomic scope of
the traits plot, rather than a direct copy of a current Evo-M1 file.


## Study 1a and Study 2 assembly

Study 1a/2 now use the same public synchronization path as every other analysis
input. The manifest pulls the canonical per-paper analysis-ready tables (standard
`<Paper>_<Item>.csv` names; no `_derived` layer) into `data_raw/`.

`scripts/00_assemble_comparison_inputs.R` then rebuilds:

- `data_raw/stereology_comparison.csv` from Morgan, Hanson, and Karlsen & Pakkenberg;
- `data_raw/stress_volume_comparison.csv` from Walhovd and Mackes/ERA.

`scripts/run_all.R` runs the sequence:

```text
00_sync_evo_m1_inputs.R
  -> 00_assemble_comparison_inputs.R
  -> 00_prepare_heiss_rates.R
  -> Study 1a / Study 2 analyses
```

The independent Heiss Table 1 source remains separate and is joined only during
`00_prepare_heiss_rates.R`. Author-provided source snapshots are built and documented
in `Evo-M1-Trait-Data`; this analysis repository does not duplicate those raw sources
or maintain a second restricted manifest.


## 12 November 2026 expansion: sync status versus analytical use

All 44 manifest sources were checked byte-for-byte against their local snapshots
on 12 November. The 25 recently added rows were synced before their analytical
integration; synchronization alone does not mean an analysis consumes a table.

Four of those 25 inputs are now consumed:

- Holmes Supplementary Table 2 by `scripts/s1c_holmes_replication.R` (healthy
  controls, BPND on its native scale; separate from Johansen).
- McLaughlin Tables 1 and 4 and Sheridan Table S3 by
  `scripts/s2_beip_replication.R` (BEIP separate from ERA, with time points and
  contrasts kept separate). The reports share a cohort, and the adolescent
  contrasts share controls. They are not independent replications of each other.

The other 21 new inputs remain synced but are not yet integrated into these
analyses. This includes the remaining stereology expansion and s5 species
additions. The source-snapshot human review remains pending; the current AI
check does not fill the human-signoff registry column.

The existing Karlsen 2011 author data now supplies **both** bilateral cell
counts and uncorrected volume (twice the author table's hemisphere volume) for
s1a. Density is the ratio of group-mean count to group-mean volume, not the mean
of individual densities. `data_intermediate/s1a_volume_basis_audit.csv` retains
the previous corrected denominator and change factor. This corrects a tissue-
basis mismatch; it does not harmonize all inter-study fixation or sampling effects.

Replication outputs include OLS, Spearman and leave-one-region-out checks.
These are exploratory, unadjusted region-level associations, not subject-level
inference or the papers' covariate-adjusted group effects. The age-16 BEIP source
has coarse 0.1-cm3 rounding and unspecified dispersion; no bilateral uncertainty
was manufactured from missing L/R covariance.

Thalamus anatomy was checked in the primary sources and documented in twelve
Evo-M1 paper folders. See [the source audit](anatomy/THALAMUS_SOURCE_AUDIT.md).
The metabolic fits omit the anatomically unmatched broad thalamus point as a
conservative current choice and retain a named proxy sensitivity. Source volumes
remain intact. Shared-region comparisons are kept within s1 composition questions
or within s2 volume-response questions; s2 is not forced onto the s1 intersection.

`run_all.R` explicitly orders the two new replication scripts, the regional
sensitivity script, and the Heiss-1984 sensitivity script (attribution corrected against 2004 reference 21). New outputs are in
`metadata/output_manifest.csv`. The Heiss sensitivity compares both rate sources
on identical complete-case regions and reports a separate white-matter exclusion.

## Strict matching continuation, 12 November

See `metadata/anatomy/PRIMARY_MATCH_POLICY.md`. Andersen Table S1 is now consumed by a separate native-scale paired-tracer replication (5 of the 25 new inputs consumed; 20 not yet integrated). No proxy rate is used to fill a missing primary anatomical match. Earlier counts and source-year labels in the handoff are superseded by the new policy and generated audit.
