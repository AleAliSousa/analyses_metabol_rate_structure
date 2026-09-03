# Project anatomy registry

`anatomy_registry.csv` is the project-wide source of truth for stable anatomy
IDs, preferred display labels, hierarchy, plotting colours, and display order.
Analysis code should join on `anatomy_id`; labels must not be used as database
keys.

`heiss_2004_region_map.csv` maps the 26 native region labels in Heiss et al.
(2004) Table 1 to registry IDs. It contains only exact regions, source-language
synonyms, and the explicit global-average measurement. It must not contain
comparative-volume proxies or calculated composite rates.

## Editing workflow

1. To change a preferred term or colour, edit its registry row without changing
   `anatomy_id`.
2. To add a genuinely new anatomical concept, add a new unique `anatomy_id` and
   then add source-specific mappings separately.
3. To add or rename a Heiss source row, update `heiss_2004_region_map.csv` and
   run `Rscript scripts/0_prepare_heiss_rates.R`.
4. Record approximations, parent-child substitutions, composite definitions,
   and weights in study-specific crosswalks. They do not belong in the registry
   or in the pure Heiss preparation.

Related concepts may currently share a colour to preserve the established
project palette, but they retain different `anatomy_id` values. For example,
`occipital_lobe` and `primary_visual_cortex_grey_matter` are not synonyms even
though both currently use the established occipital/visual blue.
