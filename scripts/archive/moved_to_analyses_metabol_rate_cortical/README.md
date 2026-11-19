# Moved 2026-11-19 → ../../../analyses_metabol_rate_cortical

The cortical and laminar s1e scripts now live (renamed) in the sibling repo `analyses_metabol_rate_cortical`:

| here (frozen copy) | there (maintained) |
|---|---|
| s1e_0_extract_metabolic_maps_julich.R | scripts/c0_extract_metabolic_maps_julich.R |
| s1e_1_join_receptors_metabolism.R | scripts/c1_receptors_metabolism_julich_areas.R |
| s1e_2_middle_tier_heiss1985.R | scripts/c2_receptors_metabolism_heiss1985_groups.R |
| s1e_4_laminar_prototype.R / s1e_4_fig_laminar_prototype.py | scripts/c3_laminar_prototype.R / scripts/c3_fig_laminar_prototype.py |

Output tables/figures here (`tables/s1e/s1e_*`, `figs/s1e/laminar_prototype.*`, `metadata/anatomy/s1e_*`, `S1E_*.md`)
are left in place as the frozen record; the maintained versions are `tables/cortical/cortical_*` etc. in the new repo
(regenerated 2026-11-19 and verified identical cell-for-cell). `s1e_3_human_macaque_groups.R` (three-species
regional comparison, not cortex-only) remains active in this repo; its Evo-M1 manifest rows are unchanged.
