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
