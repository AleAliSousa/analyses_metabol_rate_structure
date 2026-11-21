############################################################
# mouse_s1b_3_nn (+ s1b_6/s1b_7 analogues) -- Non-neuronal cell composition
# per brain structure vs mouse regional glucose utilisation. Proportions are
# fractions of non-neuronal cells, so the (sorting-dependent) neuron:non-neuron
# ratio does not enter them. Astrocytes are split into telencephalic (TE) and
# non-telencephalic (NT) types, the mouse counterpart of the human type 1/2
# astrocyte split; OPC and oligodendrocytes are separate predictors.
# Outputs: tables/mouse_s1b/nonneuronal_*.csv, figs/mouse_s1b/nonneuronal_vs_rcmr.*
############################################################
.a <- commandArgs(FALSE); .f <- sub("^--file=", "", .a[grepl("^--file=", .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else "analyses_metabol_rate_structure_nonhuman/scripts", "common.R"))
source(nh("scripts/mouse_s1b_helpers.R"))

counts <- read_mouse_counts(); roi_map <- read_roi_map(); rates <- read_mouse_rates()
sub <- read_csv(nh("metadata/mouse_subclass_to_nonneuronal_category_map.csv"), show_col_types = FALSE)
nn <- counts %>% filter(grepl("^3[0-4] ", class)) %>%
  left_join(select(sub, subclass = allen_subclass, nn_category = nonneuronal_category), by = "subclass")
stopifnot(!anyNA(nn$nn_category))

EXCLUDE <- c("Ependymal_tanycyte_other", "Other_nonneuronal", "Immune_other", "Astrocyte_other")

dp <- roi_donor_props(nn, "nn_category")
rp <- region_props(dp, roi_map, primary_only = TRUE)
wide <- rp %>% select(cmrgl_region, category, prop) %>%
  pivot_wider(names_from = category, values_from = prop, values_fill = 0, names_prefix = "p_") %>%
  inner_join(rates, by = "cmrgl_region")
write_csv(left_join(wide, distinct(rp, cmrgl_region, rois, n_donors, n_cells), by = "cmrgl_region"),
          nh("tables/mouse_s1b/nonneuronal_category_p_by_region_primary.csv"))
preds <- setdiff(grep("^p_", names(wide), value = TRUE), paste0("p_", EXCLUDE))
st <- spearman_table(wide, preds)
write_csv(st, nh("tables/mouse_s1b/nonneuronal_spearman_vs_rcmr_primary.csv"))
print(st)

# Whole-population non-neuronal fraction per region (flag: depends on Allen's neuron enrichment)
tot <- counts %>% mutate(nonneuronal = grepl("^3[0-4] ", class)) %>%
  group_by(roi, donor) %>% summarise(n = sum(n_cells), f = sum(n_cells[nonneuronal]) / n, .groups = "drop") %>%
  filter(n >= MIN_CELLS_DONOR_ROI) %>% group_by(roi) %>% summarise(nonneuronal_frac = mean(f), .groups = "drop")
write_csv(tot, nh("tables/mouse_s1b/nonneuronal_fraction_of_all_cells_by_roi.csv"))

plot_vs_rate(wide, preds, st, nh("figs/mouse_s1b/nonneuronal_vs_rcmr"),
             "Fraction of non-neuronal cells (mean across donors)", ncol = 3)
message("Wrote non-neuronal tables and figure")
