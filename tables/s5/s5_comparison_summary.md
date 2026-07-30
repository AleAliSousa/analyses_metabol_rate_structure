# Fossil-hominin brain glucose metabolism: three-method comparison

Comparison of three independent estimates of fossil-hominin brain glucose
metabolism — a regional volumetric budget (s4), a carotid-blood-flow estimate
(Seymour), and an arterial-canal regression (Boyer) — produced by
`scripts/s5_arterial_canal.R`.

Sources: `data_raw/Boyer_Harrington_2018_Table2.csv`,
`data_raw/Seymour_etal_2017_TableS1.csv`,
`data_intermediate/s4_specimen_budgets.csv`.

---

## Methods

Brain glucose metabolism of fossil hominins was estimated by three independent
methods and, for cross-method comparison, expressed as a ratio to the
modern-human value to remove differences in anatomical scope and volume
definition. First, a regional volumetric budget partitioned each endocast into
six regions matched to the resting cerebral metabolic rates of glucose
(rCMRGlc) measured by ^18F-FDG PET in Heiss et al. (2004); the grey-plus-white-
matter regional volumes reconstructed by Kochiyama et al. (2018) were converted
to mass (tissue density 1.036 g cm^-3), multiplied by the corresponding
rCMRGlc, and summed. Because Heiss et al. report only the cortical lobes and
cerebellum, this budget represents ~77% of whole-brain glucose use (328.5
versus 428.6 umol min^-1 in the modern human). Second, a cerebral-blood-flow
estimate followed Seymour et al. (2016, 2017): total internal-carotid flow rate
(Q), derived from the carotid foramen radius, was treated as proportional to
brain metabolic rate, and each specimen's whole-brain glucose use was scaled
from the modern-human reference (specimen AS8078; Q = 7.34 cm^3 s^-1; 428.6
umol min^-1) as BGU = 428.6 x (Q / 7.34). Third, a regression estimate followed
Boyer & Harrington (2018): their calibration of whole-brain glucose utilization
on total arterial canal area (ACA) and endocranial volume (ECV) was refit across
the seven extant taxa with directly measured metabolism (Karbowski 2007;
ln BGU = -0.139 + 0.440*ln ACA + 0.541*ln ECV; R^2 = 0.9997; with a Duan (1983)
smearing correction for log-retransformation bias) and applied to each fossil.
Because fossil crania preserve only the carotid canal, total ACA was
reconstructed by scaling the modern-human ACA (159.8 mm^2) by each fossil's
carotid cross-sectional area relative to the modern human, both taken from
Seymour's foramen radii; an ACA predicted from ECV through the catarrhine
ACA-ECV allometry was retained as an upper bound. Note that the endocranial
volumes used in methods two and three are cranial capacities (the whole
braincase cavity), distinct from the grey-plus-white-matter brain-tissue volumes
used in method one.

## Results

On the three specimens common to the flow and volumetric datasets
(Gibraltar/Forbes' Quarry, La Chapelle-aux-Saints, and Skhul 5), the
arterial-canal regression and the regional volumetric budget agreed closely,
differing by a mean of only 0.043 in modern-human-ratio units, whereas the
carotid-flow estimate was more variable (mean absolute difference 0.142 against
the volumetric budget), largely because Skhul 5's unusually small carotid
foramen produced an anomalously low flow-based value (0.65 versus 0.91-0.96
modern-human units from the other methods). Across all 30 fossil specimens the
flow and regression estimates were strongly correlated (Pearson r = 0.977,
p ~ 2x10^-20; Spearman rho = 0.967), indicating concordant rank ordering,
although the regression estimate was systematically higher in absolute terms
(by a mean factor of 1.47x). Group means relative to the modern human
(flow / regression) declined from recent and early *Homo sapiens* (0.93 / 0.98)
and Neanderthals (0.89 / 0.96, n = 2) through early *Homo* (0.41 / 0.56) to
*Australopithecus* (0.18 / 0.30). The methods diverged, however, on the central
question of scaling: regressing estimated glucose use on endocranial volume
across the hominin series returned a hyperallometric exponent of 1.43 for the
carotid-flow method — consistent with Seymour et al.'s inference that brain
metabolism increased faster than brain size during human evolution — but a
near-isometric 1.01 for the arterial-canal regression, so the evidence for
metabolic acceleration is contingent on whether metabolism is inferred from
carotid flow alone or from total arterial canal area together with brain volume.

---

## Supporting statistics

### Overlap specimens (ratio to modern human)

| Specimen | Seymour flow | Boyer ACA | s4 volume |
|---|---|---|---|
| Forbes'/Gibraltar | 0.83 | 0.87 | 0.82 |
| La Chapelle-aux-Saints | 0.94 | 1.06 | 1.04 |
| Skhul 5 | 0.65 | 0.91 | 0.96 |

Mean absolute difference (ratio units): Boyer vs s4 = 0.043; Seymour vs s4 =
0.142; Seymour vs Boyer = 0.139.

### Group means (ratio to modern human; whole-brain scope for arterial methods)

| Group | n | Endocranial capacity (cc) | Seymour flow | Boyer ACA |
|---|---|---|---|---|
| Recent + early *H. sapiens* | 6 | 1479 | 0.93 | 0.98 |
| Neanderthal | 2 | 1413 | 0.89 | 0.96 |
| Early *Homo* | 12 | 839 | 0.41 | 0.56 |
| *Australopithecus* | 10 | 458 | 0.18 | 0.30 |

### Early vs modern *Homo sapiens* (`sapiens_grade`)

The six arterial *H. sapiens* above pool early fossils with recent and living
humans. The `sapiens_grade` field (in `s5_fossil_estimates.csv`; means in
`s5_sapiens_grade_means.csv`) separates them:

| grade | n | specimens | Endocranial (cc) | Seymour flow | Boyer ACA |
|---|---|---|---|---|---|
| early | 3 | BC1, LH 18, Skhul 5 | 1466 | 0.84 | 0.95 |
| recent | 2 | M3-A343, M4-A344 | 1493 | 1.03 | 1.01 |
| modern reference | 1 | AS8078 | 1493 | 1.00 | 1.00 |

By the arterial methods the early group averages slightly below the recent and
modern humans, but this is carried by two individuals rather than a grade-level
effect: Skhul 5's anomalously small carotid foramen depresses its flow estimate
(0.65, the flagged outlier), and LH 18 has the smallest endocranial volume in the
*H. sapiens* set (1367 cc). Excluding the Skhul 5 flow anomaly, the early-grade
Seymour mean rises to 0.93, and the Boyer estimate (0.95) already overlaps the
modern range. On the volumetric (s4) budget — which carries all four Kochiyama
early *H. sapiens* fossils (Qafzeh 9, Skhul 5, Mladeč 1, Cro-Magnon 1) — the early
grade is indistinguishable from living humans: EH mean ratio 1.04 (± 0.07) versus
the modern (MH) baseline 1.00 and Neanderthals 1.05
(`tables/s4/species_absolute_budgets.csv`). The early-versus-modern *H. sapiens*
difference in brain glucose use is therefore small and method-dependent —
negligible on the tissue-volumetric budget, marginally lower on the arterial
estimate, and there attributable to specimen-specific vascular anatomy (Skhul 5)
and brain size (LH 18) rather than to the fossil-versus-living distinction itself.

Only Skhul 5 of Kochiyama's early-*H. sapiens* set enters the arterial (Seymour)
sample; Qafzeh 9, Mladeč 1 and Cro-Magnon 1 are volumetric-only. Fossils are never
pooled with the living-human reference in either method — the grade tag keeps them
separable (cf. the fossil-vs-extant note in the Evo-M1 specimen crosswalk).

### Cross-method statistics (all 30 Seymour specimens)

- Seymour flow vs Boyer ACA: Pearson r = 0.977 (p ~ 2x10^-20); Spearman rho = 0.967
- Boyer / Seymour absolute whole-brain BGU: mean ratio 1.47x
- Scaling of estimated glucose use on endocranial volume: exponent 1.43
  (Seymour flow, hyperallometric) vs 1.01 (Boyer ACA, near-isometric)

### Anchors and calibration

- Modern-human whole-brain glucose: 428.6 umol min^-1 (Clarke & Sokoloff 1994 /
  Boyer Table 2)
- s4 modern-human budget: 328.5 umol min^-1 (6-region cortical+cerebellar =
  76.7% of whole brain)
- Boyer calibration (n = 7): ln BGU = -0.139 + 0.440*ln ACA + 0.541*ln ECV;
  R^2 = 0.9997; smearing 1.001

## Caveats

- The Neanderthal arterial sample is only n = 2 (Gibraltar, La Chapelle).
- Seymour and s4 use different specimen sets and volume conventions (cranial
  capacity vs grey+white-matter brain tissue); only the three overlap
  individuals are strict like-for-like.
- The ECV-predicted-ACA variant of the Boyer branch is an upper bound: hominins
  fall below the general euarchontan ACA-ECV allometry.
- Absolute umol min^-1 values are comparable only within a single metabolic
  scope; all cross-method comparisons above use the ratio to modern human.

## References

- Boyer DM, Harrington AR (2018). *J. Hum. Evol.* 114, 85-101.
- Clarke DD, Sokoloff L (1994). Circulation and energy metabolism of the brain.
- Duan N (1983). Smearing estimate: a nonparametric retransformation method.
  *J. Am. Stat. Assoc.* 78, 605-610.
- Heiss W-D, et al. (2004). regional cerebral metabolic rates of glucose.
- Karbowski J (2007). Global and regional brain metabolic scaling. *BMC Biol.* 5, 18.
- Kochiyama T, et al. (2018). Reconstructing the Neanderthal brain. *Sci. Rep.* 8, 6296.
- Seymour RS, Bosiocic V, Snelling EP (2016; Correction 2017). *R. Soc. Open Sci.*
  3:160305 / 4:170846.
