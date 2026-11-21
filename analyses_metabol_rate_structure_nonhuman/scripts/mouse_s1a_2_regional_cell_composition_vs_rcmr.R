############################################################
# mouse_s1a_2 -- Regional cell-type density and ratios vs mouse regional
# glucose utilisation (CMRgl). Mouse analogue of s1a_1_stereology_cell_types
# (neurons, glia, astrocytes, oligodendrocytes, microglia, glia:neuron).
#
# Cell counts are in-situ MERFISH counts per CCF structure (mouse_s1a_1),
# pooled over the structures of each CMRgl region
# (metadata/mouse_ccf_structure_to_cmrgl_region_map.csv).
# "Density" = cells per mm^3 of CCF volume. This is a relative measure (one
# donor, sections sampled every 200 um, segmentation losses); use ranks only.
#
# Outputs: tables/mouse_s1a/*.csv, figs/mouse_s1a/*
############################################################
.a <- commandArgs(FALSE); .f <- sub("^--file=", "", .a[grepl("^--file=", .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else "analyses_metabol_rate_structure_nonhuman/scripts", "common.R"))
source(nh("scripts/mouse_s1b_helpers.R"))
for (d in c("tables/mouse_s1a", "figs/mouse_s1a")) dir.create(nh(d), recursive = TRUE, showWarnings = FALSE)

counts <- read_csv(nh("data_intermediate/allen_merfish_counts_by_ccf_structure.csv"), show_col_types = FALSE, na = character())
vol    <- read_csv(nh("data_intermediate/allen_ccf_structure_volumes.csv"), show_col_types = FALSE)
map    <- read_csv(nh("metadata/mouse_ccf_structure_to_cmrgl_region_map.csv"), show_col_types = FALSE)
cls    <- read_csv(nh("metadata/mouse_class_to_s1b_category_map.csv"), show_col_types = FALSE)
nnmap  <- read_csv(nh("metadata/mouse_subclass_to_nonneuronal_category_map.csv"), show_col_types = FALSE)
rates  <- read_mouse_rates()

# ---- checks: every mapped structure has data and a volume; rate regions exist
if (length(m <- setdiff(map$ccf_structure, counts$structure))) stop("Mapped structure(s) without cells: ", paste(m, collapse = ", "))
if (length(m <- setdiff(map$ccf_structure, vol$acronym))) stop("Mapped structure(s) without CCF volume: ", paste(m, collapse = ", "))
if (length(m <- setdiff(map$cmrgl_region, rates$cmrgl_region))) stop("Mapped region(s) without a mouse rate: ", paste(m, collapse = ", "))
stopifnot(!anyDuplicated(vol$acronym))

# ---- cell groups
is_nn <- grepl("^3[0-4] ", counts$class)
counts <- counts %>%
  left_join(select(cls, class = allen_class, ei_class), by = "class") %>%
  left_join(select(nnmap, subclass = allen_subclass, nn_category = nonneuronal_category), by = "subclass") %>%
  mutate(group = case_when(
    !is_nn & ei_class == "E" ~ "neuron_E",
    !is_nn & ei_class == "I" ~ "neuron_I",
    !is_nn ~ "neuron_other",
    nn_category %in% c("Astrocyte_TE", "Astrocyte_NT", "Astrocyte_other") ~ "astrocyte",
    nn_category == "Oligodendrocyte" ~ "oligodendrocyte",
    nn_category == "OPC" ~ "opc",
    nn_category == "Microglia" ~ "microglia",
    nn_category == "Vascular" ~ "vascular",
    TRUE ~ "other_nonneuronal"))
stopifnot(!anyNA(counts$group))

by_struct <- counts %>% group_by(structure, group) %>% summarise(n = sum(n_cells), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = n, values_fill = 0)
need <- c("neuron_E", "neuron_I", "neuron_other", "astrocyte", "oligodendrocyte", "opc", "microglia", "vascular", "other_nonneuronal")
for (g in need) if (!g %in% names(by_struct)) by_struct[[g]] <- 0

region_table <- function(map_rows, with_rates = TRUE) {
  out <- map_rows %>% inner_join(by_struct, by = c("ccf_structure" = "structure")) %>%
    inner_join(select(vol, acronym, volume_mm3), by = c("ccf_structure" = "acronym")) %>%
    group_by(cmrgl_region) %>%
    summarise(across(all_of(need), sum), volume_mm3 = sum(volume_mm3),
              n_structures = n(), structures = paste(ccf_structure, collapse = "+"), .groups = "drop") %>%
    mutate(neurons = neuron_E + neuron_I + neuron_other,
           glia = astrocyte + oligodendrocyte + opc + microglia,
           total_cells = neurons + glia + vascular + other_nonneuronal,
           neuron_density = neurons / volume_mm3, glia_density = glia / volume_mm3,
           astrocyte_density = astrocyte / volume_mm3, oligodendrocyte_density = oligodendrocyte / volume_mm3,
           microglia_density = microglia / volume_mm3,
           neuron_fraction = neurons / total_cells,
           glia_to_neuron_ratio = glia / neurons,
           astrocyte_to_neuron_ratio = astrocyte / neurons,
           oligodendrocyte_to_neuron_ratio = oligodendrocyte / neurons,
           microglia_to_neuron_ratio = microglia / neurons,
           E_fraction_of_neurons = neuron_E / neurons)
  if (with_rates) inner_join(out, rates, by = "cmrgl_region") else out
}

PRED <- c("neuron_density", "glia_density", "astrocyte_density", "oligodendrocyte_density", "microglia_density",
          "glia_to_neuron_ratio", "astrocyte_to_neuron_ratio", "oligodendrocyte_to_neuron_ratio",
          "microglia_to_neuron_ratio", "neuron_fraction")

sets <- list(
  primary_all            = filter(map, primary_included, !nested_in_other_region),
  primary_grey_matter    = filter(map, primary_included, !nested_in_other_region, cmrgl_region != "Corpus_callosum"),
  with_partial_regions   = filter(map, (primary_included | cmrgl_region %in% c("Parietal_cortex", "Basal_ganglia", "White_matter")),
                                  cmrgl_region != "Neocortex"))
res <- lapply(names(sets), function(nm) {
  tb <- region_table(sets[[nm]]); write_csv(tb, nh(sprintf("tables/mouse_s1a/regional_cell_composition_%s.csv", nm)))
  st <- spearman_table(tb, PRED); write_csv(st, nh(sprintf("tables/mouse_s1a/spearman_vs_rcmr_%s.csv", nm)))
  cat("\n==", nm, "(n regions =", nrow(tb), ")\n"); print(as.data.frame(st), digits = 3)
  list(tb = tb, st = st)
}); names(res) <- names(sets)
plot_vs_rate(res$primary_grey_matter$tb, PRED, res$primary_grey_matter$st,
             nh("figs/mouse_s1a/cell_composition_vs_rcmr_grey_matter"), "Regional cell measure (MERFISH, CCF-resolved)", ncol = 4)

# ---- S1 (SSp-*) profile vs the rest of isocortex (nested pool shown for reference)
iso <- map %>% filter(cmrgl_region == "Neocortex") %>% pull(ccf_structure)
ssp <- grep("^SSp-", iso, value = TRUE)
prof <- bind_rows(
  mutate(region_table(tibble(ccf_structure = ssp, cmrgl_region = "S1_SSp"), with_rates = FALSE), group = "S1 (SSp-*)"),
  mutate(region_table(tibble(ccf_structure = setdiff(iso, ssp), cmrgl_region = "Neocortex"), with_rates = FALSE), group = "other isocortex")) %>%
  select(group, volume_mm3, neurons, glia, neuron_density, glia_density, glia_to_neuron_ratio,
         astrocyte_to_neuron_ratio, oligodendrocyte_to_neuron_ratio, microglia_to_neuron_ratio, E_fraction_of_neurons)
write_csv(prof, nh("tables/mouse_s1a/S1_vs_other_isocortex_cell_composition.csv")); print(as.data.frame(prof), digits = 3)

# ---- validity check against isotropic fractionator (Herculano-Houzel 2011)
hh_path <- nh("data_raw/Evo_M1_mouse_isotropic_fractionator.csv")
if (file.exists(hh_path)) {
  hh <- read_csv(hh_path, show_col_types = FALSE); g <- function(v) hh$Value[hh$Variable == v]
  hh_tab <- tibble(region = c("Neocortex", "Cerebellum"),
    iso_frac_neuron_fraction = c(g("CerebralCortex_N.n") / g("CerebralCortex_C.n"), g("Cerebellum_N.n") / g("Cerebellum_C.n")),
    iso_frac_other_to_neuron = c(g("CerebralCortex_O.n") / g("CerebralCortex_N.n"), g("Cerebellum_O.n") / g("Cerebellum_N.n")))
  cc <- region_table(bind_rows(filter(map, cmrgl_region %in% c("Neocortex", "Cerebellum")))) %>%
    transmute(region = cmrgl_region, merfish_neuron_fraction = neuron_fraction,
              merfish_nonneuron_to_neuron = (total_cells - neurons) / neurons)
  write_csv(left_join(hh_tab, cc, by = "region"), nh("tables/mouse_s1a/validation_vs_isotropic_fractionator.csv"))
  print(as.data.frame(left_join(hh_tab, cc, by = "region")), digits = 3)
}
message("Wrote mouse s1a tables and figure")
