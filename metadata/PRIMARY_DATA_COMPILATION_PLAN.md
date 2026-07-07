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

## Data-quality caveats on candidate sources (per DeCasien & Higham 2019)

DeCasien & Higham (2019, *Nat. Ecol. Evol.*), in their comparative brain-region
compilation, flagged specific problems with two sources adjacent to the ones planned
above. Recording these here so we apply the same caution before pooling:

> "Specifically, although Semendeferi and Damasio (ref. 66) presents whole brain
> volumes, this measure excludes the medulla, the pons and most of the midbrain.
> Navarrete et al. (ref. 67) explicitly compare their data to others (refs 24,51,52,59)
> and note several regions for which there are marked, statistically significant
> differences in average brain region volumes (for example, the hippocampus is up to
> 60% smaller in their dataset). On investigation of these data, we found further
> inconsistencies with earlier data, which may be based on differences in the regional
> boundaries used for measurement. Since regional boundary information is not available,
> even after contacting the authors, we are not using these data at this time — please
> see the authors' published erratum (ref. 68)."

**Implications for this compilation and for Study 3:**

1. **Semendeferi whole-brain ≠ our `Total_brain_net_volume`.** Their whole-brain
   measure omits medulla, pons and most of the midbrain, so it is *not* interchangeable
   with the Stephan `Total_brain_net_volume` used as the Study 3 predictor / rest-of-brain
   base. If Semendeferi parietal/temporal **lobe** volumes are added, use only the lobe
   values and keep the Stephan total as the predictor; do **not** substitute or pool their
   whole-brain figure, and check lobe-to-total scaling against Stephan before merging
   specimen-to-specimen.

2. **Navarrete et al. data: exclude, following DeCasien & Higham.** The marked,
   unexplained regional discrepancies (e.g. hippocampus up to ~60% smaller) and the
   absence of regional-boundary definitions make these values unsafe to pool. We adopt
   the same position and do not ingest Navarrete volumes; see the authors' published
   erratum (ref. 68).

3. Both points reinforce the merge rule above ("Primary beats derived; check scaling
   consistency before pooling") — differences in **regional boundary definitions** are
   the likely driver of cross-dataset disagreement, so any new lobe/region source must be
   boundary-checked against Stephan before it enters `volumes_wide` / `volumes_long`.

*Refs as numbered in DeCasien & Higham 2019:*
- *66. Semendeferi, K. & Damasio, H. The brain and its main anatomical subdivisions in living
  hominoids using magnetic resonance imaging. J. Hum. Evol. 38, 317–332 (2000).*
- *67. Navarrete, A. F. et al. Primate brain anatomy: new volumetric MRI measurements for
  neuroanatomical studies. Brain Behav. Evol. 91, 1–9 (2018).*
- *68. Erratum. Brain Behav. Evol. 92, 182–184 (2019).*

## Files
- `data_raw/primary_volume_compilation_long_TEMPLATE.csv` — the empty scaffold to fill.
- This plan.
