############################################################
# mouse_s1b_2 -- Validate the Allen ROI -> Evo-M1 CMRgl region map against
# the extracted cell counts and the mouse rate table (mouse analogue of
# s1b_2_mapping_*). Fails on any inconsistency; never rewrites the map.
# Output: data_intermediate/mouse_s1b_roi_mapping_coverage.csv
############################################################
.a <- commandArgs(FALSE); .f <- sub("^--file=", "", .a[grepl("^--file=", .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else "analyses_metabol_rate_structure_nonhuman/scripts", "common.R"))
source(nh("scripts/mouse_s1b_helpers.R"))

counts <- read_mouse_counts(); map <- read_roi_map(); rates <- read_mouse_rates()

if (length(miss <- setdiff(unique(counts$roi), map$allen_roi)))
  stop("Allen ROI(s) not in the map: ", paste(miss, collapse = ", "), call. = FALSE)
if (length(orph <- setdiff(map$allen_roi, unique(counts$roi))))
  stop("Mapped ROI(s) absent from the data: ", paste(orph, collapse = ", "), call. = FALSE)
if (length(bad <- setdiff(na.omit(map$cmrgl_region), rates$cmrgl_region)))
  stop("Mapped region(s) without a mouse rate: ", paste(bad, collapse = ", "), call. = FALSE)
if (any(map$primary_included & is.na(map$cmrgl_region)))
  stop("primary_included ROIs must map to a rate region.", call. = FALSE)

cls <- read_csv(nh("metadata/mouse_class_to_s1b_category_map.csv"), show_col_types = FALSE)
sub <- read_csv(nh("metadata/mouse_subclass_to_nonneuronal_category_map.csv"), show_col_types = FALSE)
nn_class <- grepl("^3[0-4] ", unique(counts$class))
if (length(u <- setdiff(unique(counts$class[!grepl("^3[0-4] ", counts$class)]), cls$allen_class)))
  stop("Neuronal class(es) missing from class map: ", paste(u, collapse = ", "), call. = FALSE)
if (length(u <- setdiff(unique(counts$subclass[grepl("^3[0-4] ", counts$class)]), sub$allen_subclass)))
  stop("Non-neuronal subclass(es) missing from map: ", paste(u, collapse = ", "), call. = FALSE)

cov <- counts %>% group_by(roi) %>%
  summarise(n_cells = sum(n_cells), n_donors = n_distinct(donor),
            nonneuronal_frac = sum(n_cells[grepl("^3[0-4] ", class)]) / sum(n_cells), .groups = "drop") %>%
  left_join(map, by = c("roi" = "allen_roi")) %>%
  left_join(rates, by = "cmrgl_region") %>%
  arrange(desc(primary_included), cmrgl_region, roi)
write_csv(cov, nh("data_intermediate/mouse_s1b_roi_mapping_coverage.csv"))
print(as.data.frame(select(cov, roi, n_cells, n_donors, nonneuronal_frac, cmrgl_region, primary_included, rcmr_value)), digits = 3)
message("Validated ", nrow(map), " ROI mappings; ", sum(map$primary_included), " primary")
