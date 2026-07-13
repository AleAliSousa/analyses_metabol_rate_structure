# Stephan_primates — per-datapoint reference resolution

Built by `scripts/build_stephan_primates_reference_sheet.R` on 2026-07-13.

## What this does

The old metadata (`Stephan_primates_metadata.xlsx`) records references **per
variable** — a column can list several candidate sources (e.g. `Body_weight` →
*Stephan et al 1981; Zilles and Rehkamper 1988; Bauernfeind et al 2013*). This
pipeline resolves, for **each data point** (species × column), the *specific*
source table the value came from, by value-matching every cell of
`data_raw/Stephan_primates.csv` against the per-source provenance exported from
the organised **Evo-M1-Trait-Data** project.

## Inputs

- `data_raw/Stephan_primates.csv` — the compiled data (60 species × 77 columns).
- `metadata/Stephan_primates_metadata.xlsx` — per-variable reference groupings.
- `data_intermediate/volumes_unfiltered.csv` — Evo-M1 per-source values
  (mirror of `Evo-M1-Trait-Data/__merging_volumes/volumes_unfiltered.csv`).
- `data_intermediate/Stephan_etal_1981_TablesI-VI.csv` — Evo-M1 mirror of the
  Stephan 1981 combined table, used as a fallback (see finding 1).

> The two `data_intermediate/` files are mirrors from Evo-M1; refresh them when
> the Evo-M1 merge is re-run.

## Outputs

- **`Stephan_primates_references_long.csv`** — one row per non-empty data cell,
  with the resolved specific reference(s). Key columns: `preferred_reference`,
  `preferred_citation`, `all_matching_sources`, `status`, `metadata_reference`,
  `unit_rescale_flag`, `merge_deviation_flag`, `closest_nonmatching`.
- **`Stephan_primates_reference_mismatches.csv`** — cells to review:
  `value_differs`, `matched_outside_metadata`, unit-rescale, and merge-deviation
  cells.
- **`Stephan_primates_references_by_column.csv`** — per-column coverage.
- **`Stephan_primates_source_hunt.csv`** — for every unresolved cell
  (`value_differs`, `matched_outside_metadata`, `provenance_gap`), the file and
  row in the raw cut-and-paste directory (`datasets brain and regions species/
  Stephan and Frahm`) where the value actually appears — i.e. where the compiler
  copied it from.

## Coverage (2,268 data cells)

| status | n | meaning |
|---|---|---|
| resolved_unique | 1,389 | one source carries the value |
| resolved_multiple | 429 | several sources share it (preferred picked by the rule below) |
| metadata_only | 240 | no Evo-M1 volumetric equivalent (see below) |
| provenance_gap | 105 | Evo-M1 has the variable but not this species |
| matched_outside_metadata | 68 | value matches an Evo-M1 source that is *not* the metadata source — **review** |
| value_differs | 37 | value matches no Evo-M1 source — **review** |

`metadata_only` columns (reference taken straight from the xlsx): the seven
Smaers 2017 cortical grey/white columns, `Meninges_hypophysis_nerves_etc`,
`Brainvol` (de Sousa 2010), and `Brain_volume` (Bauernfeind 2013).

## How the preferred reference is chosen

This follows the universal rule of the Evo-M1 volumes merge (`volumes_compiled.R`
§5–6), **not** a hand-picked priority:

1. **Scope to the metadata source.** Candidate sources are first restricted to
   the paper(s) the metadata lists for that column, so a coincidental identical
   number in another team or structure is not mistaken for the source. (This is
   also why the older Stephan 1970 monograph is not chosen where the metadata
   lists Stephan 1981 — 1970 is simply not in scope, not a special case.)
2. **Within scope, apply the merge rule:** per team, take the **most recent**
   source; the Stephan_collection team is used first; body/brain **mass →
   Stephan 1981** (as the merge does).
3. **Deviations** flagged by the merge (newest-vs-next within a team, i.e. the
   typo/anomaly cases in `volumes_flags.csv`) are carried through in
   `merge_deviation_flag`.

`matched_outside_metadata` means the Stephan_primates value matched an Evo-M1
source that is not the metadata's stated source (e.g. `ASG_Sousa`/`LGN_Sousa`,
which are de Sousa's columns but whose numbers coincide with the Stephan/Frahm
data de Sousa built on; or a body/brain-mass value that coincidentally equals a
Stephan volume). These are attributed to the metadata source and flagged so you
can confirm; the `source_hunt` file shows where the value really came from.

A value that matches **no** source (`value_differs`) is the real surprise to
watch for — that is where unpublished/raw data would show up.

## Matching rules

A Stephan value matches a source value if they agree within 1% (relative), or
agree after rounding the source to the Stephan value's decimal places (handles
2-sig-fig rounding of small volumes), or both are zero. Body/brain-mass columns
are additionally tried at ×1000 / ×0.001 and the rescale is flagged. Vestibular
columns are matched against Evo-M1's `*_unilateral_*` variants because the
Stephan file reports one-side volumes (Stephan 1981 Table XIII).

---

## Findings that need your decision (raw vs published)

**1. Gorilla uses Stephan 1981 *revised* values — and Evo-M1's merge lost them.**
The Gorilla structural columns (Total_brain_net_volume = 470359, Telencephalon =
369878, Body_weight = 105000, Brain_weight = 500000, and all fundamental parts,
~7% larger than the 1970 values) match `Stephan_etal_1981_TablesI-VI.csv`
**exactly**. But Evo-M1's `volumes_unfiltered.csv` carries only the older
Stephan **1970** Gorilla row. So these are *not* raw/unpublished numbers — they
are the published 1981 revision, and the discrepancy is an **Evo-M1 merge gap**:
the 1981 Gorilla (and Homo) rows were dropped from the compilation.
*Action:* reference = Stephan et al 1981; and the Evo-M1 merge should ingest the
1981 Gorilla/Homo rows. (The pipeline now resolves these via a 1981 fallback.)

**2. Brain-mass unit inconsistency in Evo-M1.** 38 `Brain_weight`/`Brain_mass`
cells only match Evo-M1 at ×1000: Stephan_primates stores brain mass in **mg**
(internally consistent), whereas Evo-M1's `Brain_Mass.mg` holds **grams** for
these primates (e.g. *Alouatta caraya* 55.8 in Evo-M1 vs 55800 in Stephan).
*Action:* fix the Evo-M1 unit (or its label); Stephan_primates is correct.

**3. Orangutan split (14 cells).** Stephan_primates separates *Pongo abelii* and
*Pongo pygmaeus*, but Evo-M1 pools both as *Pongo sp.* Body/brain mass, de Sousa
brain size, and Bauernfeind insula subfields differ between the two files.
*Action:* reconcile the taxon (see Evo-M1 `specimen-taxon-tracking`); decide
which orangutan values belong to which species in both projects.

**4. de Sousa's own area-striata / LGN (6 cells).** `ASG_Sousa` and `LGN_Sousa`
for Homo, Pan and Hylobates differ (3–22%) from Evo-M1's Frahm / Bush-Allman /
deSousa-2013 values for the same structures. These columns are de Sousa's own
measurements (metadata already attributes them to de Sousa 2010/2013), so the
difference is expected — different measurement, not an error. No change needed
unless you want to align on one measurement.

**5. Minor precision (4 cells).** *Cebus albifrons* subcommissurale (1.2%),
*Lepilemur* lateral vestibular nucleus (3.4%), *Nycticebus* subfornicale (2.8%),
*Callithrix pygmaea* brain mass (1.2%) — effectively the same value at coarser
rounding; the source is still identifiable in `closest_nonmatching`.

**6. Genuinely different, worth a look (small set).** *Tarsius syrichta*
neocortex grey/white vs Bush-Allman (~9%); *Avahi laniger* / *Callithrix
pygmaea* body mass vs Stephan/Frahm (8–17%, Bauernfeind column).

## Metadata note

The xlsx lists `NeoG_Frahm`/`NeoW_Frahm` as *"Frahm et al 1992"* — this is a
typo for **Frahm et al 1982** (Frahm, Stephan & Stephan, *J Hirnforsch* 23:375),
which is the source that actually carries the neocortex grey/white values. The
script treats "Frahm 199x" as Frahm 1982 when choosing the preferred reference.

## Important

The script only **reads** the data; it never edits `Stephan_primates.csv`.
After deciding each case above, change the value in both projects by hand.
