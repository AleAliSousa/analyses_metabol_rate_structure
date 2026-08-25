# Citations & data provenance — s4a endocranial energy-budget analysis

This file records where each number in the `s4a_endocranial*` analyses comes from,
and — importantly — keeps the **metabolic-rate inputs** (Heiss) separate from the
**whole-brain benchmark values** used only as an external sanity check (which are
NOT from Heiss).

## 1. Metabolic rates (rCMRGlc) — the model input

**Heiss W-D, Habedank B, Klein JC, Herholz K, Wienhard K, Lenox M, Nutt R.
"Metabolic Rates in Small Brain Nuclei Determined by High-Resolution PET."
J Nucl Med. 2004;45(11):1811–1815.** (PMID 15534048)

- PDF in the repo: `Evo-M1-Trait-Data/Heiss_etal_2004/Heiss_etal_2004.pdf`
- Local table: `data_raw/Heiss_etal_2004_TABLE1.csv`
- Method: ¹⁸F-FDG PET on a high-resolution research tomograph (HRRT, 2.2 mm
  FWHM) in 9 healthy volunteers (8 M, 1 F; mean age 35 y, range 25–65).
- **Units: µmol glucose / 100 g / min.** Regional means (both hemispheres),
  e.g. cerebral cortex 33.5, frontal 35.3, parietal 35.8, temporal 30.5,
  occipital 35.8, cerebellar cortex 29.8, vermis 30.1, white matter 12.3.
- **What this paper does NOT contain:** any whole-brain total glucose use, any
  grams-of-glucose-per-day figure, or any "% of body metabolism" statement. It
  is a table of regional *rates* only. Every absolute budget in this project is
  those rates × Kochiyama regional volumes — see §2.
- The Heiss Table 1 also carries a historical comparison column ("Data of Heiss
  et al., ref 21") from an **earlier Heiss study**; that earlier work is also
  regional-rate data (used only for cross-reference), not a whole-brain total.
  Its full citation is not pinned here — cite from the 2004 paper's reference
  list (ref 21) if needed.

## 2. Regional brain volumes — the other model input

**Kochiyama T, Ogihara N, Tanabe HC, et al. "Reconstructing the Neanderthal
brain using computational anatomy." Sci Rep. 2018;8:6296.**
DOI 10.1038/s41598-018-24331-0.

- Fig 3A: ICV-size-adjusted *relative* parcel volumes (shape), NT/EH/MH.
- Fig 3 legend: modern-human absolute parcel volumes (cc), n = 1185.
- Cerebrum/cerebellum group means from the fossil-specimen supplementary text.

## 3. Whole-brain benchmark values — external sanity check ONLY (not Heiss)

These were used in conversation to check that our measured-region total
(≈ 329 µmol/min for the 6 modern-human regions, ≈ 85 g glucose/day for ~75 %
of the brain) lands in a physiologically sensible place. They are **not** part
of the calculation and are **not** from Heiss.

- **Clarke DD, Sokoloff L. "Circulation and Energy Metabolism of the Brain."
  In: Siegel GJ et al. (eds), Basic Neurochemistry. Raven Press, New York,
  1994, pp. 645–680.** — human whole-brain glucose: CMRglc ≈ 0.31 µmol/g/min,
  total ≈ 428.6 µmol/min (≈ 111 g glucose/day for a ~1382 g brain), CMRO2
  ≈ 0.035 ml/g/min. This is the whole-brain benchmark used in the summing
  check; the value is recorded in the Karbowski (2007) energetics compilation
  in this repo (see below), which is where we read it.
  NOTE: an earlier in-conversation figure of "~120 g/day" was an over-estimate
  from extrapolating the high *cortical* rate to the whole brain (ignoring
  low-rate white matter); the grounded whole-brain value is ~111 g/day.
- **Karbowski J. "Global and regional brain metabolic scaling and its
  functional consequences." BMC Biol. 2007;5:18.** — a multi-species compilation
  of regional/whole-brain CMRgl, CMRO2 and CBF. Used here only as the source
  that tabulates the Clarke & Sokoloff 1994 human whole-brain values.
  Local: `Evo-M1-Trait-Data/Karbowski__2007/`.
- Foundational methods behind whole-brain glucose measurement: **Kety & Schmidt
  (1948)** (whole-brain blood flow / metabolism) and **Sokoloff et al. (1977)**
  (2-deoxyglucose method). Not verified against a retrieved source in this
  session — confirm the full citations before using them in a paper.

The whole-brain summing check that uses these values lives in
`Evo-M1-Trait-Data/__energetics_comparison/heiss_wholebrain_check.R`.

## 3b. Phylogeny — the PGLS covariance structure

- **Upham NS, Esselstyn JA, Jetz W. "Inferring the mammal tree: species-level
  sets of phylogenies for questions in ecology, evolution, and conservation."
  PLoS Biol. 2019;17(12):e3000494. doi:10.1371/journal.pbio.3000494** — source
  of every tree used here. Two items from it are in play:
  - *Single tree*: the DNA-only MCC (4,098 sp, `topoFree NDexp`), which is what
    the legacy `data_raw/species.nwk` descends from.
  - *100-tree sample*: `MamPhy_BDvr_Completed_5911sp_topoCons_NDexp`, 100 trees
    drawn from the published 10k pseudoposterior via vertlife.org
    (download job recorded in the source folder's `config.yaml`). Used for
    PGLS-across-trees sensitivity analysis, so topological and branch-length
    uncertainty enters the result rather than being fixed by one point estimate.
    In the *Completed* set, species without DNA are placed by the **authors'**
    taxonomic imputation, differently in each posterior tree.
  Local: `Evo-M1-Trait-Data/Upham_etal_2019/`.
- Tip-name reconciliation (e.g. `Callithrix pygmaea` → `Cebuella pygmaea`,
  `Plecturocebus moloch` → `Callicebus moloch`) uses only synonym tables already
  recorded in that repo — `__merging_trees/tree_tip_crosswalk.csv`,
  `_keys/species_display_aliases.csv`, `_keys/species_reference.csv`. Nothing is
  grafted or imputed locally: a species no published tip carries is reported
  absent in the coverage report and left off the trees.

`scripts/s3_prepare_tree_sample.R` builds the sample and writes
`data_intermediate/s3_tree_coverage_report.csv`, which states how every input
species was handled.

## 4. Unit conversion used throughout

- Volume → mass: cc × 1.036 g/cc (brain tissue density).
- Rate is per 100 g, so region budget = rate × (mass_g / 100) × shape × size.
- Glucose molar mass 180.16 g/mol for any µmol/min → g/day conversion.

_Written 2026-07-10. §1–2 verified against the repo PDFs/CSVs; §3 whole-brain
values are external literature benchmarks, cited for the sanity check only._
