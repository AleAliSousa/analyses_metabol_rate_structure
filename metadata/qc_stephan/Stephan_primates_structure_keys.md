# Stephan primates structure keys

This note documents how columns in `data_raw/Stephan_primates.csv` are keyed to
Evo-M1 standardized variables. The executable one-row-per-column map is
`stephan_wide_crosswalk.csv`; `build_stephan_primates_reference_sheet.R` and
`compare_stephan.R` both read that file.

## Keying rule

A key is accepted from anatomical description, laterality, composition, and
units in the source table caption, footnotes, or methods. Similar values are a
check on an established key, not evidence that two structures are equivalent.
Consequently:

- a match in another table must still use the same documented structure;
- left, right, unilateral, author-estimated bilateral, and measured bilateral
  terms remain distinct;
- components and totals are not aliases, even when one can be derived from the
  other;
- a source-conditioned exception is written explicitly rather than made a
  global term alias; and
- `unit_scale` converts the Stephan column into the standardized Evo unit
  before comparison (`1000` for cm3 to mm3).

The crosswalk contains 65 defensible keys. Three numeric columns remain
`metadata_only`: `Brain_volume`, `Brainvol`, and
`Meninges_hypophysis_nerves_etc`. The first two are generic whole-brain fields
whose relationship to Stephan's explicitly net brain volume is not sufficiently
defined; the third is not a brain-structure volume.

## Definition-sensitive keys

| Stephan column(s) | Evo-M1 key | Anatomical basis and constraint | Source evidence |
|---|---|---|---|
| `Pons` | `Ventral_pons_Vol.mm3` for Matano; `Pons_Vol.mm3` for Zilles | This column is heterogeneous by source. Matano defines and measures only ventral/basilar pons (VPo). The orangutan row is the broader pons printed by Zilles. These terms are never global aliases. | Matano et al. 1985b title, Introduction, and Table I; Zilles & Rehkamper 1988 Table 12-2. |
| `ASG_Sousa` | `Area_striata_grey_matter_Vol.mm3` | Grey matter only. de Sousa Supplementary Table 2 is an author-side bilateral estimate, obtained by doubling measured left V1. It is not Stephan 1981 `Area_striata`, which includes underlying white matter. | de Sousa et al. 2010 Table 1 and Methods; Stephan et al. 1981 Table IX footnote 27. |
| `LGN_Sousa` | `Corpus_geniculatum_laterale_Vol.mm3` | Entire LGN; Supplementary Table 2 publishes the author-estimated bilateral value (2 x left). The measured left Table 1 term remains separate. | de Sousa et al. 2010 Table 1 and Methods, p. 4. |
| `NeoWG` | `Neocortex_Vol.mm3` | Stephan code 18 includes neocortical grey, underlying white matter, corpus callosum, and the insular, subgenual, cingular, and retrosplenial transitional cortices. | Stephan et al. 1981 Tables IV-VI, footnote 18. |
| `NeoG_Frahm`, `NeoW_Frahm` | grey/white neocortex terms | Frahm's neocortex includes grey laminae 1-6 and white matter with corpus callosum, excludes internal capsule, and includes some transitional cortex but excludes presubicular and entorhinal cortex. | Frahm et al. 1982 Methods, p. 376, and Table 2. |
| `Total_brain_net_volume` | `Total_brain_net_volume_Vol.mm3` | Stephan code 5 is explicitly the sum of codes 6-10. A generic whole-brain volume is not accepted as an anatomical alias. | Stephan et al. 1981 Tables I-III and footnote 5. |
| `Cerebellum` | `Cerebellum_Vol.mm3` | Stephan code 7 includes brachium and nuclei pontis. This does not make a separately published pons/VPo value part of the `Pons` key. | Stephan et al. 1981 footnote 7. |
| `Mesencephalon`, `Medulla_oblongata` | corresponding terms | Mesencephalon excludes substantia reticularis; the inseparable reticular portion is assigned to medulla oblongata. | Stephan et al. 1981 footnotes 6 and 8. |
| `Diencephalon` | `Diencephalon_Vol.mm3` | Includes pallidum; excludes hypophysis, which is assigned to code 4. | Stephan et al. 1981 footnote 9. |
| `Lobus_piriformis` | `Lobus_piriformis_Vol.mm3` | Includes palaeocortex and amygdala; excludes schizocortex. | Stephan et al. 1981 footnote 13. |
| `Schizocortex` | `Schizo_cortex_Vol.mm3` | Entorhinal, perirhinal, presubicular, and parasubicular cortices plus underlying white matter. | Stephan et al. 1981 footnote 16. |
| insula subdivision columns ending `_L` | corresponding `*_left_Vol.mm3` terms | Bauernfeind Table 1 is measured left hemisphere and undoubled. `FI` is the VEN-containing frontoinsular part; total left insula is the sum of the printed subdivisions. | Bauernfeind et al. 2013 Table 1, Methods pp. 265-266, and Table 3. |
| vestibular complex and nuclei | corresponding `*_unilateral_Vol.mm3` terms | Stephan 1981 values were measured on one side; the side is not identified. They must not be pooled with Baron 1988 measured-bilateral values. | Stephan 1981 Tables XII-XIII, clarified by Baron et al. 1988. |
| cerebellar nuclei columns | corresponding nucleus terms | Matano 1985a cerebellar-nucleus measures; the historical `lateral` nucleus is the nucleus dentatus cerebelli. | Matano et al. 1985a table caption and definitions. |
| seven Smaers cortical columns | `Area_striata_white_matter`, `PrefrontalCortex_*`, `FrontalMotor_*`, and `AssociationCortex_*` | Published values are both-hemisphere cm3. Prefrontal and frontal-motor measures are opposite-end, five-of-20-section proxies and do not partition the whole frontal lobe. Other association cortex is the published residual `neocortex - whole frontal lobe - primary visual`. | Smaers et al. 2017 Table S1 and Supplemental Experimental Procedures; Smaers et al. 2011 Supplementary Tables 1-2. |

## Straight canonical keys

The remaining mappings use the same source-defined structure on both sides,
with spelling normalized by the Evo standard term. These include the Stephan
1981 fundamental parts; olfactory, septal, striatal, hippocampal, diencephalic,
amygdaloid, and periventricular structures; body/brain mass; and Matano 1985a
cerebellar nuclei. Their detailed inclusions remain in the source-specific
`*_definitions.csv` files. A matching label does not erase those definitions:
for example, Stephan `Striatum` includes the internal-capsule portions running
through the striatum, and `Palaeocortex` includes the structures enumerated in
the 1981 footnote.

## How the reference build uses the keys

1. Read the primary standardized-variable key and `unit_scale` from
   `stephan_wide_crosswalk.csv`.
2. Apply only explicit source-conditioned alternatives (`Pons` is currently the
   sole case).
3. Search the selected Evo ledger, then other finalized public volume tables,
   then the restricted document inventory.
4. Require the candidate to satisfy the structure key. Restricted candidates
   additionally require species and structure words in the same row/column
   context. An exact number elsewhere in a listed paper is not enough.
5. Use numerical agreement and a source's prior contributions to rank the
   anatomically eligible candidates.

When a source definition or key changes, update together:

- the source folder's `*_definitions.csv` and table README;
- `stephan_wide_crosswalk.csv`;
- this note; and
- any source-aware rule in `build_stephan_primates_reference_sheet.R`.
