# Primary-data compilation for the four gap-filling sources

_Drafted 2026-06-05. Goal: extend the Stephan composite with PRIMARY (histology / histology+MRI)
volumes from Smaers, MacLeod, Barger, Semendeferi, recording provenance so same-specimen vs
separate-specimen merges are explicit. Schema below is a proposal that maps 1:1 onto the
Evo-M1-Trait-Data columns once you share a sample row — it is long-format, so remapping is a rename._

## Principle: primary data, with provenance

Prioritise measurements made directly on specimens, and record **which collection** each came from.
Per your note, the **Zilles and Stephan collections are now housed together and combined across many
papers**; Smaers, MacLeod (Düsseldorf portion) and Semendeferi measured *different slides / scans of
the same brains*. So those three can be merged **specimen-to-specimen** with Stephan (ideal), while
Barger's apes are **separate brains** and merge only at the **species-mean** level.

## Proposed schema (long format — one row per species × structure × source)

| column | meaning |
|---|---|
| `species` | binomial matching `species.nwk` / `Stephan_primates.csv` |
| `heiss_region` | the Heiss 2004 Table-1 region this value serves |
| `structure_measured` | exact structure named in the source (e.g. "cerebellar vermis") |
| `hemisphere` | both / left / right / unilateral×2 |
| `value_mm3` | volume in mm³ (convert at entry; note original units) |
| `value_sd`, `n_individuals` | dispersion + specimen count for that species |
| `measurement_type` | histology / MRI / histology+MRI / stereology |
| `data_type` | **primary** / secondary-derived |
| `source_short`, `source_doi_pmid` | citation + ID |
| `specimen_collection` | e.g. "Düsseldorf C&O Vogt (Zilles/Stephan)", "Yerkes", "ape-aging network" |
| `same_specimens_as_stephan` | yes / partly / no  → controls specimen vs species-mean merge |
| `table_ref` | table/figure/page the value came from |
| `notes` | caveats, conversions |

## Sources × coverage (pre-filled; values pending the source tables)

| Source | Heiss region(s) filled | Structures measured | Method | Species / N | Collection | Same brains as Stephan? | Already in your data? |
|---|---|---|---|---|---|---|---|
| **Smaers et al. 2017** (Curr Biol) + 2011/2018 | Frontal lobe (+ assoc. → parietal/temporal lump) | prefrontal & motor frontal grey+white; non-frontal association grey+white | histology | ~18–19 anthropoids incl. 6 apes + Homo | Düsseldorf C&O Vogt (Zilles) | **yes** (likely same brains) | **partly** — `Prefrontal.Gray`, `Frontal.motor.Gray`, `Other.cortical.association.areas.Gray` (+ white) already in `Stephan_primates.csv` |
| **MacLeod et al. 2003** (J Hum Evol) | Cerebellar cortex, **Vermis** | whole cerebellum, vermis, lateral hemispheres | MRI + histology | 97 (42 apes, 14 humans, 41 monkeys) | Yerkes **+** Düsseldorf "Hirnforschung" | **partly** (Düsseldorf cases yes; Yerkes no) | no per-species columns — only used as rCMRGlc weights |
| **Barger et al. 2007, 2012** (J Comp Neurol) | Corpus amygdaloideum | whole amygdala + lateral/basal/accessory-basal/central etc. nuclei | stereology | apes + humans (great-ape sample) | ape-aging / comparative collections | **no** (separate brains → species-mean merge) | no — Stephan amygdala only |
| **Semendeferi et al. 1997, 2000, 2002** | **Parietal lobe**, **Temporal lobe**, Frontal lobe (+ limbic) | frontal/parietal/temporal/occipital lobes; limbic | MRI (on Zilles-collection brains) | apes + humans | Zilles collection (per your note: same brains) | **yes** | no |

## What each source closes (vs the Heiss inventory gaps)

- **Parietal lobe & Temporal lobe** (the Tier-4 cortical-lobe gaps) → **Semendeferi** lobe volumes.
- **Vermis / cerebellar cortex** as real per-species volumes (not just proportions) → **MacLeod**.
- **Ape amygdala** (Stephan amygdala is ape-poor) → **Barger** (species-mean merge).
- **Frontal lobe** is already largely covered by your **Smaers 2017** columns; use Smaers white/grey
  directly rather than re-deriving.

## Merge rules

1. **Specimen-level (Smaers, MacLeod-Düsseldorf, Semendeferi):** where a source reports the same
   individual brains as Stephan, join on specimen/`Code_number` if available, else species mean, and
   set `same_specimens_as_stephan = yes/partly`. These are internally consistent with Stephan volumes.
2. **Species-mean (Barger, MacLeod-Yerkes):** merge as a species-level augmentation; flag `no`. Use
   only to fill species Stephan lacks, and check scaling consistency before pooling.
3. **Primary beats derived:** if a structure is available both as a primary measurement and as a
   Stephan-derived proportion (e.g. vermis), keep the primary value and retain the derived one only as
   a fallback (`data_type` distinguishes them).

## What I still need from you

1. **The Evo-M1-Trait-Data template** — a sample row + column list (or point Cowork at the folder).
   I'll rename my columns to match; nothing else changes.
2. **The source tables/PDFs** for the actual numbers. MacLeod 2003, Barger 2012, Semendeferi are
   paywalled and not web-fetchable here; if you drop the PDFs/supp spreadsheets in the workspace (or
   mount the OneDrive Stephan-team folder), I'll extract the per-species volumes into the template.
   Smaers values you largely already have.

## Files
- `data_raw/primary_volume_compilation_long_TEMPLATE.csv` — the empty scaffold to fill.
- This plan.
