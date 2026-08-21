# Stephan_primates — per-datapoint reference resolution

Built by `scripts/qc_stephan/build_stephan_primates_reference_sheet.R`; last
verified 2026-08-20.

## Purpose

`metadata/Stephan_primates_metadata.xlsx` records references per variable, but
one variable can list several papers. The script resolves a source for each
numeric data point in `data_raw/Stephan_primates.csv` and now broadens the
search when the metadata-listed value is not a near-perfect match.

The search is intentionally ordered:

1. Match the species, standardized structure, and value in the selected
   Evo-M1 provenance ledger.
2. For a weak cell, search every table belonging to the metadata-listed
   paper(s), including other printed tables, while retaining the documented
   structure key.
3. Search the canonical all-table Evo-M1 volume ledger.
4. If mounted, search dated volume inputs in
   `Evo-M1-Trait-Data-restricted`.

Only sources dated 2019 or earlier are eligible. The cutoff uses publication or
manually curated source dates, not filesystem modification times, because the
files have been moved and re-synced since the compilation was assembled.

## Inputs

- `data_raw/Stephan_primates.csv`
- `metadata/Stephan_primates_metadata.xlsx`
- `metadata/qc_stephan/stephan_wide_crosswalk.csv` — authoritative anatomy and
  unit key shared with the comparison scripts
- `Evo-M1-Trait-Data/__merging_volumes/volumes_unfiltered_select.csv` —
  selected direct-match ledger
- `Evo-M1-Trait-Data/__merging_volumes/volumes_flags_select.csv`
- `Evo-M1-Trait-Data/__merging_volumes/volumes_unfiltered.csv` — canonical
  inventory of all final brain-volume tables
- `Evo-M1-Trait-Data/Stephan_etal_1981/Stephan_etal_1981_TablesI-VI.csv`
- finalized de Sousa 2010 Supplementary Table 2 and Smaers 2017 Table S1 public
  CSVs. They are valid provenance tables even though the modern merge excludes
  some secondary/derived rows from its pooled estimate.
- `Evo-M1-Trait-Data-restricted/MIGRATED_INDEX.csv` and the explicitly listed
  unpublished volume folders, when mounted

Set `EVOM1_RESTRICTED=/path/to/Evo-M1-Trait-Data-restricted` to override the
default restricted-repo location. If the private repo is absent, the script
warns and completes the public search.

## Outputs

- `Stephan_primates_references_long.csv` — one row per numeric cell. New audit
  fields include `structure_key_rule`, `resolution_method`, `fallback_candidate_sources`,
  `fallback_match_detail`, `preferred_source_prior_contributions`, and
  `source_cutoff_year`.
- `Stephan_primates_reference_mismatches.csv` — unresolved, likely, fallback,
  rescaled, and merge-deviation cases for review.
- `Stephan_primates_references_by_column.csv` — coverage by source column.
- `Stephan_primates_candidate_source_inventory.csv` — the 58 public source
  tables and 50 restricted documents eligible under the cutoff in this run.
- `Stephan_primates_fallback_candidates.csv` — all fallback candidates with
  tier, difference, contribution count, and a source locator.
- `Stephan_primates_source_hunt.csv` — remaining unresolved cells traced into
  the older local cut-and-paste directory.

Restricted rows are inspected in memory but are not copied into these outputs;
the audit records only the private file/sheet/cell locator and match flags.

## Structure keys

`Stephan_primates_structure_keys.md` documents how anatomy, composition,
laterality, and units establish every key. Source definitions come before
numeric similarity: an exact value elsewhere in a metadata-listed paper is not
eligible unless the table heading/methods support the same structure. The
source-specific `*_definitions.csv` files in Evo-M1 retain the detailed printed
definitions.

## Resolution and ranking rules

A near-perfect match agrees within 1% relative difference, agrees at the
precision printed in `Stephan_primates.csv`, or is zero in both files. Mass
columns also try ×1000 and ×0.001. Restricted workbooks try those scale factors
for all structures because their raw units vary.

For cells with no near-perfect direct match, a value up to 5% away can be
reported as `likely_fallback_*`, but only when the same structure and a
metadata-listed/eligible source support it. It is not treated as resolved.

Within the anatomically eligible candidates in a search tier, a table already preferred for other clean cells ranks
ahead of a numerical one-off. This is recorded in
`preferred_source_prior_contributions`. Metadata-listed sources rank ahead of
other public volume sources, and restricted matches require species and
structure evidence from the same row/column context.

`Pons` is source-conditioned. Forty rows are Matano et al. 1985b ventral pons
and map to `Ventral_pons_Vol.mm3`; the orangutan row is Zilles & Rehkamper's
whole pons and maps to `Pons_Vol.mm3`. The two terms are not global aliases.

Seven Smaers cortical columns use the explicit Table S1 definitions and a
cm3-to-mm3 scale. Prefrontal and frontal-motor volumes are the anterior and
posterior five-section proxies, not a two-part decomposition of frontal lobe;
association cortex is the source's derived residual.

## Coverage from the 2026-08-20 run

| status | n | meaning |
|---|---:|---|
| `resolved_unique` | 1,601 | one direct selected-ledger/source-key source |
| `resolved_multiple` | 299 | multiple direct sources; one preferred |
| `resolved_fallback_metadata` | 168 | exact match in another metadata-listed public table |
| `resolved_fallback_restricted` | 53 | exact/near-exact private-document evidence with structure context |
| `resolved_fallback_evom1` | 2 | exact match in another final Evo-M1 table |
| `likely_fallback_metadata` | 14 | same-source/structure candidate 1–5% away |
| `likely_fallback_restricted` | 1 | private metadata-source candidate 1–5% away |
| `metadata_only` | 118 | one of three columns with no defensible volume key |
| `value_differs` | 11 | no supported candidate |
| `provenance_gap` | 1 | no supported candidate row |

The 12 unresolved cells are now tightly concentrated: 11 orangutan mass/insula
cells caused by the *Pongo abelii* / *Pongo pygmaeus* split versus pooled
`Pongo sp.` sources, plus `Macaca fascicularis` `Body_weight = 2900`.

## Main findings

- The finalized public de Sousa Supplementary Table 2 now supplies V1/LGN
  provenance directly rather than through a restricted working copy.
- All available values in the seven Smaers cortical columns resolve to the
  finalized 2017 Table S1 using their documented bilateral cm3 definitions.
- The restricted fallback retains the Avahi/Callithrix Bauernfeind rows, the
  missing Avahi Stephan structures, and other source-level locators, but no
  longer promotes exact numbers that lack structure context.
- Fourteen small discrepancies (1–5%) retain a most-likely metadata source but
  remain explicitly `likely`, including several printed-rounding cases and
  pooled orangutan means.
- Gorilla's revised fundamental-part values continue to resolve to Stephan et
  al. 1981 via the combined-table fallback.
- The existing brain-mass unit inconsistency remains visible through
  `unit_rescale_flag`; `Stephan_primates` uses mg while several Evo-M1 source
  rows carry gram-scale numbers under the mg term.
- The metadata label “Frahm et al 1992” for `NeoG_Frahm`/`NeoW_Frahm` is treated
  as the known typo for Frahm et al. 1982.
- The only metadata-only fields are `Brain_volume`, `Brainvol`, and
  `Meninges_hypophysis_nerves_etc`; generic whole brain is not silently equated
  with Stephan's explicitly net-brain construct.

The script only reads `Stephan_primates.csv`. Any decision to change a data
value still needs to be made explicitly in the source datasets.
