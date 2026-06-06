# Comparative primate volume data for the Heiss et al. (2004) rCMRGlc regions

_Inventory + equivalencies, drafted 2026-06-05. "Higher primate" = monkeys, apes, humans; the
central question is which structures have **great-ape + human** volume data, not just strepsirrhine-
heavy Stephan coverage._

## The big-picture constraint

There is still **no second multi-region comparative volumetric dataset** to rival Stephan, Frahm &
Baron (1981) — DeCasien & Higham (2019) analysed 33 regions drawn essentially *from* the Stephan/
Frahm data, and Barton & Harvey (2000) did the same. So the realistic strategy is the one you have
already been following and the one Smaers pioneered: **take Stephan as the backbone and graft on
targeted, region-specific datasets that add the apes (Pan, Gorilla, Pongo, Hylobates) and refine
humans.** The honest weak point is *subcortical* great-ape coverage: for nuclei like pallidum, STN,
thalamus, mesencephalon, capsula interna, the only apes in the classic tables are *Pongo* (via Zilles
& Rehkämper 1988) plus *Homo*; Pan/Gorilla volumes for those nuclei are largely absent. Cortex, insula,
V1, LGN and cerebellum are where the ape-inclusive add-ons exist.

## What you already have (multi-source composite)

From `metadata/Stephan_primates_metadata.xlsx`, your `Stephan_primates.csv` is already assembled from:

| Source | Adds | Ape/human reach |
|---|---|---|
| Stephan, Frahm & Baron 1981 (+ Stephan 1984) | backbone: all fundamental parts + telencephalic/diencephalic subdivisions, striatum, hippocampus, amygdala(+subnuclei), pallidum, STN, capsula interna, thalamus, etc. | Homo; apes essentially absent |
| Zilles & Rehkämper 1988 | *Pongo* values across the Stephan structures | the main great-ape graft for subcortex |
| Frahm & Stephan 1982 (+ Frahm 1984) | neocortex grey vs white; area striata (V1) | as Stephan |
| de Sousa et al. 2010, 2013 | V1 (area striata grey), LGN | adds apes incl. Pan, Gorilla, Hylobates, Pongo + Homo |
| Bauernfeind et al. 2013 | insula (total, FI, granular/dysgranular/agranular) | anthropoid-wide incl. all great apes + Homo |
| Smaers et al. 2017 | frontal-lobe grey+white along the prefrontal↔motor axis, association cortex | 19 anthropoids incl. 6 apes (with Homo) |
| Matano et al. 1985 (pons), Matano 2001 (dentate) | pons; dentate nucleus | apes + human |

So you are not starting from scratch — you are at ~16 of the 26 Heiss regions with usable, mostly
ape-aware equivalencies. The exercise below is really "where can the remaining ~10 be filled, and
which of my existing equivalencies are approximations to flag."

## Master equivalency table (Heiss 2004 → primate volume data)

Legend for **Status**: ✅ have it, ape-aware · 🟡 have a proxy/approximation · 🔶 combined only (no
clean subdivision) · ⛔ genuine gap for apes.

### Cerebral cortex
| Heiss region | rCMRGlc | Volume equivalency | Build / source | Status |
|---|---|---|---|---|
| Cerebral cortex (global average) | 33.5 | Neocortex grey (total) | Frahm 1982 `NeoG_Frahm`; = volume-weighted mean of the lobes, so it *contains* occipital(V1) + insula | ✅ |
| Frontal lobe | 35.3 | Prefrontal + frontal-motor grey | `Prefrontal.Gray + Frontal.motor.Gray` (Smaers 2017); cross-check Semendeferi 1997/2002, Bush & Allman 2004, Rilling & Insel 1999 | ✅ (apes) |
| Parietal lobe | 35.8 | — | no clean parietal-only volume in Stephan; partition Smaers "Other cortical association" or use Rilling & Insel 1999 MRI lobes | ⛔ |
| Temporal lobe | 30.5 | (neocortical temporal + hippocampus + entorhinal) | Stephan has hippocampus + schizocortex(entorhinal); neocortical temporal only via Rilling & Insel 1999 / Semendeferi MRI | 🔶 |
| Occipital lobe | 35.8 | V1 / area striata grey | `ASG_Sousa` (de Sousa 2010) — **V1 only; occipital lobe also has extrastriate**, so this underestimates the lobe | 🟡 |
| Insular lobe | 30.3 | Insula total | `Total_insula_volume_L` (Bauernfeind 2013) | ✅ (apes) |
| Hippocampus | 25.7 | Hippocampus | `Hippocampus` (Stephan; Zilles & Rehkämper); cf. DeCasien & Higham 2019 (hippocampus), Barger et al. | ✅ |

### Forebrain nuclei
| Heiss region | rCMRGlc | Volume equivalency | Build / source | Status |
|---|---|---|---|---|
| Caudatum | 32.6 | part of Striatum | Stephan `Striatum` is combined; split via de Jong 2012 proportions, or individual MRI (sparse for apes) | 🔶 |
| Putamen | 40.2 | part of Striatum | as above | 🔶 |
| Nucleus accumbens | 32.3 | part of Striatum | as above | 🔶 |
| Pallidum | 21.9 | Pallidum | `Pallidum` (Stephan; Pongo via Zilles) | ✅ (Pongo+Homo only) |
| Nucleus medial thalami | 36.6 | part of Thalamus | Stephan `Thalamus` (whole); MD-nucleus-specific only via Armstrong 1980s | 🔶 |
| Corpus geniculatum laterale | 23.1 | LGN | `LGN_Sousa` (de Sousa 2010/2013; Stephan 1984; Barton 1998) | ✅ (apes) |
| Corpus geniculatum mediale | 20.2 | — | MGN volume rarely measured comparatively | ⛔ |
| Nucleus subthalamicus | 25.2 | STN | `Nucleus_subthalamicus` (Stephan) | ✅ (Pongo+Homo) |
| Basal forebrain | 21.4 | (Septum as partial proxy) | Stephan `Septum`; nucleus basalis-specific comparative data sparse (Raghanti et al.) | 🔶 |
| Corpus amygdaloideum | 22.2 | Amygdala | `Amygdala` (Stephan + subnuclei); Barger et al. 2007/2012 add apes+human subnuclei | ✅ |

### Brain stem nuclei  (Heiss measured 4 nuclei; Stephan has only the whole mesencephalon)
| Heiss region | rCMRGlc | Volume equivalency | Build / source | Status |
|---|---|---|---|---|
| Colliculus inferior | 31.4 | part of Mesencephalon | Stephan `Mesencephalon` whole; individual ape volumes sparse (you already weight rCMRGlc via Garcia-Gomar 2019) | 🔶 |
| Colliculus superior | 24.4 | part of Mesencephalon | as above | 🔶 |
| Substantia nigra | 22.8 | part of Mesencephalon | as above; Eapen 2011 (human) used for rCMRGlc weighting | 🔶 |
| Nucleus ruber | 31.0 | part of Mesencephalon | as above | 🔶 |

### Cerebellum
| Heiss region | rCMRGlc | Volume equivalency | Build / source | Status |
|---|---|---|---|---|
| Cerebellar cortex | 29.8 | part of Cerebellum | Stephan `Cerebellum`; cortex/vermis split via MacLeod 2003 (apes+humans); Rilling & Insel 1998 | ✅ (apes) |
| Vermis | 30.1 | part of Cerebellum | MacLeod et al. 2003 (vermis vs hemispheres, apes+humans) | ✅ (apes) |
| Nucleus dentatus cerebelli | 24.2 | part of Cerebellum | Matano 2001 | 🟡 |

### White matter
| Heiss region | rCMRGlc | Volume equivalency | Build / source | Status |
|---|---|---|---|---|
| Centrum semiovale | 12.3 | Neocortex white | `NeoW_Frahm` (Frahm 1982); cf. Smaers white, Rilling & Insel 1999 | ✅ |
| Capsula interna | 25.5 | Capsula interna | `Capsula_interna` (Stephan) | ✅ (Pongo+Homo) |

## "From scratch": how to assemble each region, by tractability

**Tier 1 — clean, ape-aware equivalency already available (8):** Cerebral cortex global (NeoG),
Insular lobe (Bauernfeind), Hippocampus, LGN (de Sousa), Amygdala (Stephan/Barger), Cerebellar
cortex + Vermis (MacLeod), Neocortex white, plus Frontal lobe (Smaers prefrontal+motor).

**Tier 2 — present but as an approximation to flag (4):** Occipital lobe ≈ V1 only (add extrastriate);
Pallidum / STN / Capsula interna are fine *measures* but ape coverage = Pongo + Homo only;
Dentate nucleus from a single proportion source.

**Tier 3 — only as a combined parent, needs a split rule (8):** Caudatum/Putamen/Nucleus accumbens
(→ Striatum, split via de Jong 2012); the four mesencephalic nuclei (→ Mesencephalon); Nucleus medial
thalami (→ Thalamus); Cerebellar dentate. These can be carried at the *parent* level with a
literature split — exactly the weighted-average device you already use for the rCMRGlc side.

**Tier 4 — genuine gaps for apes (≈3):** Parietal lobe, Temporal lobe (as clean neocortical lobes),
Corpus geniculatum mediale (MGN). For the lobes, the only realistic ape+human source is MRI
parcellation (Rilling & Insel 1999; Heuer et al. 2019 open MRI; Semendeferi/Sherwood).

## How the field built these equivalencies

- **Barton & Harvey (2000)** and **Barton (1998)** used Stephan 1981 (+ Stephan 1984 for the visual
  system: LGN, V1, optic tract) more or less directly, with whole-brain / rest-of-brain as the size
  covariate — the same part-whole logic as your Study 3 — to argue mosaic vs. concerted evolution.
- **DeCasien & Higham (2019)** assembled **33 regions** from the Stephan/Frahm corpus with an updated
  phylogeny and PGLS; their dataset is the closest thing to a ready-made, curated version of yours and
  is worth pulling as a cross-check.
- **Smaers et al. (2011, 2017, 2018; Smaers & Soligo 2013)** is the template for *adding apes*: rather
  than rely on Stephan's ape-poor sample, they generated histological frontal-lobe and cerebellar
  volumes for ~19 anthropoids including the great apes + humans. Your `Smaers 2017` columns already
  embed this.

## Recommended additions (and which gaps they close)

1. **Rilling & Insel 1999** (MRI, 11 species incl. apes+humans) + **Rilling & Insel 1998** (cerebellum):
   gives lobe-level neocortex (frontal, **temporal**, **parietal**-ish), grey/white, and cerebellum for
   apes+humans → closes Tier-4 lobes and strengthens Temporal.
2. **Semendeferi et al. (1997, 2002), Sherwood et al.**: ape/human frontal & temporal lobe and limbic
   volumes → corroborate Frontal/Temporal.
3. **Barger et al. (2007, 2012)**: amygdala + subnuclei in apes+humans → strengthens Amygdala beyond
   Pongo+Homo.
4. **Heuer et al. (2019)**: open MRI for **34 primate species** (Paris collection, via BrainBox) —
   cerebral volume + folding; broad species net incl. apes, useful for neocortex/cerebral size and as
   an independent check, though it is not finely parcellated subcortically.
5. **DeCasien & Higham (2019) supplementary data**: a curated 33-region Stephan-based table to diff
   against your composite.
6. Human anchors you already hold — **Karlsen & Pakkenberg 2011**, **Morgan et al. 2014** — pin the
   *Homo* values for the subcortical/striatal subdivisions where comparative apes are missing.

## Caveats that matter for the equivalencies

- **"Cerebral cortex (global average)" is the volume-weighted mean of the cortical lobes**, so it
  already includes occipital(V1) and insula. This is the formal justification for the Study-3 v2
  de-overlap (remaining neocortex = neocortex − V1 − insula) and for reweighting its rCMRGlc.
- **Occipital lobe ≠ V1**: V1/area striata is the most metabolically distinctive part but the lobe also
  contains extrastriate cortex; matching `ASG_Sousa` to "Occipital lobe" is a deliberate
  (defensible) simplification.
- **Frahm neocortex** explicitly *excludes* internal capsule and *includes* the corpus callosum and
  some transitional cortices (insular, subgenual, cingulate, retrosplenial), and measures entorhinal/
  presubicular separately (per your metadata Notes). That boundary definition should drive the
  V2 correspondence, not Heiss's vaguer lobe labels.
- **Subcortical ape coverage is thin** (mostly Pongo + Homo). This is the same clade-imbalance issue
  flagged in `PHASE1_missing_data_strategy.md`, and it is *worse* for nuclei than for cortex — which is
  precisely where phylogenetic imputation on the raw data (before any correspondence step) would help.

---

### Sources
- DeCasien & Higham 2019, *Nat Ecol Evol* — primate mosaic brain evolution (33 regions): https://www.nature.com/articles/s41559-019-0969-0
- DeCasien, Williams & Higham 2017, *Nat Ecol Evol* — brain size, diet vs sociality: https://www.nature.com/articles/s41559-017-0112
- Barton & Harvey 2000, *Nature* — mosaic evolution of brain structure: https://www.nature.com/articles/35016580
- Rilling & Insel 1999, *J Hum Evol* — primate neocortex via MRI (11 species): https://pubmed.ncbi.nlm.nih.gov/10444351/
- Smaers et al. 2017, *Current Biology* — prefrontal expansion in great apes & humans: https://pubmed.ncbi.nlm.nih.gov/28162899/
- Smaers et al. 2018, *PNAS* — quantitative prefrontal assessment: https://www.pnas.org/doi/10.1073/pnas.1721653115
- Heuer et al. 2019, *Cortex* — neocortical folding, MRI of 34 primate species: https://pubmed.ncbi.nlm.nih.gov/31235272/
- Local: `metadata/Stephan_primates_metadata.xlsx`, `metadata/matching brain structures to metabolic.xlsx`, `data_raw/Heiss_etal_2004_TABLE1.csv`
