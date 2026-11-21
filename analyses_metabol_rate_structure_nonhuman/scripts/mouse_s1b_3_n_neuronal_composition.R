############################################################
# mouse_s1b_3_n -- Neuronal cell-class composition per brain structure vs
# mouse regional glucose utilisation. Mouse analogue of
# s1b_3_n_transcriptomic_neuronal (human). Proportions are mean donor x ROI
# fractions of neurons, averaged over the ROIs of each CMRgl region.
# Outputs: tables/mouse_s1b/neuronal_*.csv, figs/mouse_s1b/neuronal_vs_rcmr.*
############################################################
.a <- commandArgs(FALSE); .f <- sub("^--file=", "", .a[grepl("^--file=", .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else "analyses_metabol_rate_structure_nonhuman/scripts", "common.R"))
source(nh("scripts/mouse_s1b_helpers.R"))

counts <- read_mouse_counts(); roi_map <- read_roi_map(); rates <- read_mouse_rates()
cls <- read_csv(nh("metadata/mouse_class_to_s1b_category_map.csv"), show_col_types = FALSE)

neurons <- counts %>% filter(!grepl("^3[0-4] ", class)) %>%
  left_join(select(cls, class = allen_class, s1b_category), by = "class")
stopifnot(!anyNA(neurons$s1b_category))

# Exclusions parallel the human script: MSNs are basal-ganglia specific, and
# neuromodulatory/other subcortical GABA classes are regional rather than graded.
EXCLUDE <- c("Inhibitory_projection_MSN", "Neuromodulatory")

run <- function(primary_only, tag) {
  dp <- roi_donor_props(neurons, "s1b_category")
  rp <- region_props(dp, roi_map, primary_only)
  wide <- rp %>% select(cmrgl_region, category, prop) %>%
    pivot_wider(names_from = category, values_from = prop, values_fill = 0, names_prefix = "p_") %>%
    inner_join(rates, by = "cmrgl_region")
  info <- rp %>% distinct(cmrgl_region, rois, n_donors, n_cells)
  write_csv(left_join(wide, info, by = "cmrgl_region"), nh(sprintf("tables/mouse_s1b/neuronal_category_p_by_region_%s.csv", tag)))
  preds <- setdiff(grep("^p_", names(wide), value = TRUE), paste0("p_", EXCLUDE))
  st <- spearman_table(wide, preds)
  write_csv(st, nh(sprintf("tables/mouse_s1b/neuronal_spearman_vs_rcmr_%s.csv", tag)))
  list(wide = wide, preds = preds, stats = st)
}
prim <- run(TRUE, "primary"); sens <- run(FALSE, "with_partial_regions")
print(prim$stats); cat("n regions (primary):", nrow(prim$wide), "\n")

# Finer-grained view: cortical subclasses (IT/ET/NP/CT/L6b and interneuron groups) in S1 vs rest
plot_vs_rate(prim$wide, prim$preds, prim$stats, nh("figs/mouse_s1b/neuronal_vs_rcmr"),
             "Neuronal class proportion (mean across donors)")
message("Wrote neuronal tables and figure")
