# Kochiyama_etal_2018_Figure3

## Source

PDF: `Kochiyama-2018-Reconstructing the Neanderthal.pdf`

Paper: Kochiyama, T., Tanabe, H. C., Sawada, R., Kondo, O., Amano, H., Suzuki, H., Homma, S., Ogihara, N., et al. (2018). *Reconstructing the Neanderthal brain using computational anatomy.* Scientific Reports, 8, 6296. https://doi.org/10.1038/s41598-018-24331-0

Figure: **Figure 3. Comparisons of the relative volumes of the parcellated brain regions among NT, EH and MH.**

> **Citation / naming note.** The task brief cited the 2018 book chapter (Kochiyama et al., in *Digital Endocasts*, Springer) and the token `Kochiyama_etal_2019`. However, the repo registry (`__ReadMe.xlsx`), the existing folder, and the PDF all correspond to the **Scientific Reports 2018** paper (DOI 10.1038/s41598-018-24331-0), which is the source that contains "Figure 3 ... relative volumes of the parcellated brain regions among NT, EH and MH" with the data in the legend. This dataset was therefore built as **`Kochiyama_etal_2018_Figure3`** to match the registry (so the DOI-encoded TSV lookup resolves). NT = Neanderthal, EH = early Homo sapiens, MH = modern Homo sapiens.

## What the figure reports

The numeric data are printed in the **figure legend**: the mean (± s.d.) MH (modern human) volumes of the **13 parcellated brain regions** in cc, plus the ANOVA F/p for the four regions with a significant relative-volume difference among the three groups (Pa SI, Oc SM, Ce V, Ce P; df = 2,1190, Bonferroni p < 0.003). Panel (a) shows the NT/EH/MH relative-volume bars, but the NT and EH values are **not** given numerically in the text and are not transcribed (they would require figure digitization). Panel (b) shows absolute Ce A / Ce P volumes by hemisphere and group as bars — also not numerically in the text.

## Files

| Path | Role |
|---|---|
| `Kochiyama_etal_2018_Figure3_snapshot.csv` | Faithful capture of the legend data (caption + 13 region rows: MH mean/sd cc, ANOVA F/p). Source of truth. |
| `Kochiyama_etal_2018_Figure3.R` | Preparation: snapshot -> CSV (+ DOI-named TSV). Maps region codes to canonical structure/subregion (from the legend key) and converts cc -> mm3. |
| `Kochiyama_etal_2018_Figure3.csv` | Analysis-ready data, 13 rows (one per region). "use this". |
| `reference_tables/Kochiyama_etal_2018_Figure3_definitions.csv` | Data dictionary. |

No `comparison/` folder: there is no independent curated copy of this figure's data to audit against (allowed by the pipeline — comparison is only required when such a source exists).

## Units & conversions

Region volumes are printed in **cc (cm³)**. The analysis CSV keeps the printed cc values (`MH_mean_Vol.cc`, `MH_sd_Vol.cc`) and adds project-unit columns in **mm³** (`× 1000`).

## Sample sizes

Reconstructions used 1,185 living-human MR brains and eight fossil endocasts (four NT: Amud 1, La Chapelle-aux-Saints 1, La Ferrassie 1, Forbes' Quarry 1; four early H. sapiens: Qafzeh 9, Skhul 5, Mladeč 1, Cro-Magnon 1). The MH means are over n = 1185 (consistent with the ANOVA within-group df of 1190 for three groups of 1185 / 4 / 4).

## Data role

`primary` for the reconstructed MH region volumes (Kochiyama et al.'s own measurements), but note they are **computationally reconstructed** modern-human region volumes for a single species; add to a merge only deliberately. ANOVA columns are `secondary` (test statistics).
