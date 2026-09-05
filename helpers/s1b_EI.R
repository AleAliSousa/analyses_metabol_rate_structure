# Shared calculation and plotting functions for the s1b E:I analyses.
# Anatomical scope and cell-class definition remain separate and explicit.

s1b_jorstad_cell_class_map <- function(
  path = "metadata/cell_taxonomy/s1b_jorstad_to_siletti_cell_class_map.csv"
) {
  if (!file.exists(path)) {
    stop("Missing s1b Jorstad-to-Siletti cell-class map: ", path, call. = FALSE)
  }
  map <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA")
  )
  required <- c(
    "jorstad_class", "jorstad_neighborhood", "jorstad_subclass",
    "jorstad_description", "EI_class", "siletti_supercluster_term",
    "siletti_cluster_ids", "mapping_level", "mapping_status", "mapping_note"
  )
  missing <- setdiff(required, names(map))
  if (length(missing)) {
    stop("Cell-class map is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(map$jorstad_subclass)) {
    stop("Jorstad subclasses must be unique in the cell-class map.", call. = FALSE)
  }
  expected_subclasses <- c(
    "L2/3 IT", "L4 IT", "L5 IT", "L6 IT", "L6 IT Car3",
    "L5 ET", "L5/6 NP", "L6 CT", "L6b",
    "Lamp5", "Lamp5 Lhx6", "Pax6", "Sncg", "Vip",
    "Chandelier", "Pvalb", "Sst", "Sst Chodl"
  )
  if (!setequal(map$jorstad_subclass, expected_subclasses)) {
    stop(
      "Cell-class map must cover every neuronal subclass in Jorstad Table S1.",
      call. = FALSE
    )
  }
  if (!all(map$EI_class %in% c("E", "I"))) {
    stop("EI_class must contain only E or I.", call. = FALSE)
  }
  if (!all(map$mapping_level %in% c("supercluster", "cluster"))) {
    stop("mapping_level must contain only supercluster or cluster.", call. = FALSE)
  }
  if (any(map$mapping_level == "cluster" & is.na(map$siletti_cluster_ids))) {
    stop("Cluster-level mappings require siletti_cluster_ids.", call. = FALSE)
  }
  l5et <- map[map$jorstad_subclass == "L5 ET", , drop = FALSE]
  if (
    nrow(l5et) != 1L ||
    l5et$siletti_supercluster_term != "Miscellaneous" ||
    l5et$mapping_level != "cluster"
  ) {
    stop("L5 ET must use a cluster-specific selector within Miscellaneous.", call. = FALSE)
  }
  map
}

s1b_neocortical_ei_crosswalk <- function() {
  s1b_jorstad_cell_class_map() %>%
    dplyr::group_by(
      siletti_supercluster_term, siletti_cluster_ids,
      mapping_level, EI_class
    ) %>%
    dplyr::summarise(
      class_equivalent = paste(jorstad_subclass, collapse = " + "),
      .groups = "drop"
    ) %>%
    dplyr::rename(supercluster_term = siletti_supercluster_term)
}

s1b_siletti_l5et_cluster_ids <- function() {
  ids <- s1b_jorstad_cell_class_map() %>%
    dplyr::filter(jorstad_subclass == "L5 ET") %>%
    dplyr::pull(siletti_cluster_ids)
  unique(trimws(unlist(strsplit(ids, ";", fixed = TRUE))))
}

attach_s1b_jorstad_ei_class <- function(obs) {
  required <- c("supercluster_term", "cluster_id")
  missing <- setdiff(required, names(obs))
  if (length(missing)) {
    stop("Siletti observations are missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  crosswalk <- s1b_neocortical_ei_crosswalk()
  indexed <- obs %>%
    dplyr::mutate(
      .s1b_ei_row_id = dplyr::row_number(),
      .s1b_cluster_id = as.character(cluster_id)
    )

  supercluster_map <- crosswalk %>%
    dplyr::filter(mapping_level == "supercluster") %>%
    dplyr::select(supercluster_term, EI_class, class_equivalent, mapping_level)
  cluster_map <- crosswalk %>%
    dplyr::filter(mapping_level == "cluster") %>%
    tidyr::separate_rows(siletti_cluster_ids, sep = ";") %>%
    dplyr::transmute(
      supercluster_term,
      .s1b_cluster_id = trimws(siletti_cluster_ids),
      EI_class,
      class_equivalent,
      mapping_level
    )

  hits <- dplyr::bind_rows(
    indexed %>% dplyr::inner_join(supercluster_map, by = "supercluster_term"),
    indexed %>%
      dplyr::inner_join(
        cluster_map,
        by = c("supercluster_term", ".s1b_cluster_id")
      )
  )
  if (anyDuplicated(hits$.s1b_ei_row_id)) {
    stop("A Siletti observation maps to more than one Jorstad E:I class.", call. = FALSE)
  }

  hits %>%
    dplyr::arrange(.s1b_ei_row_id) %>%
    dplyr::select(-.s1b_ei_row_id, -.s1b_cluster_id)
}

s1b_siletti_broad_e_terms <- function() {
  c(
    "Amygdala excitatory",
    "Deep-layer corticothalamic and 6b",
    "Deep-layer intratelencephalic",
    "Deep-layer near-projecting",
    "Upper-layer intratelencephalic",
    "Thalamic excitatory",
    "Hippocampal CA1-3",
    "Hippocampal CA4",
    "Hippocampal dentate gyrus",
    "Mammillary body",
    "Lower rhombic lip",
    "Upper rhombic lip"
  )
}

s1b_siletti_broad_i_interneuron_terms <- function() {
  c(
    "CGE interneuron",
    "MGE interneuron",
    "LAMP5-LHX6 and Chandelier",
    "Cerebellar inhibitory",
    "Midbrain-derived inhibitory"
  )
}

s1b_siletti_broad_i_msn_terms <- function() {
  c("Medium spiny neuron", "Eccentric medium spiny neuron")
}

# Cell classes appropriate to the project's cerebral-cortex anatomical scope.
# These are intentionally narrower than the whole-brain broad E/I lists above:
# hippocampal and amygdaloid excitatory classes are retained, but superclusters
# associated with thalamus, hypothalamus, midbrain, hindbrain, or cerebellum are
# excluded even when a small number of their cells occur in a cortical ROI.
s1b_siletti_cerebral_cortex_e_terms <- function() {
  c(
    "Amygdala excitatory",
    "Deep-layer corticothalamic and 6b",
    "Deep-layer intratelencephalic",
    "Deep-layer near-projecting",
    "Upper-layer intratelencephalic",
    "Hippocampal CA1-3",
    "Hippocampal CA4",
    "Hippocampal dentate gyrus"
  )
}

s1b_siletti_cerebral_cortex_i_terms <- function() {
  c(
    "CGE interneuron",
    "MGE interneuron",
    "LAMP5-LHX6 and Chandelier"
  )
}

calculate_s1b_neocortical_class_ei <- function(
  scope_flag,
  scope_label,
  interpretation,
  obs_path = "data_intermediate/linnarsson_adult_human_brain_obs_metadata_neuronal.rds"
) {
  allowed_scopes <- c("is_neocortex", "is_cerebral_cortex")
  if (!scope_flag %in% allowed_scopes) {
    stop(
      "scope_flag must be one of: ", paste(allowed_scopes, collapse = ", "),
      call. = FALSE
    )
  }
  if (!exists("attach_s1b_anatomy", mode = "function")) {
    stop("Source helpers/read_heiss_rates.R before this helper.", call. = FALSE)
  }

  obs <- readRDS(obs_path) %>%
    dplyr::filter(cell_type == "neuron") %>%
    dplyr::mutate(
      roi = stringr::str_squish(as.character(roi)),
      supercluster_term = stringr::str_squish(as.character(supercluster_term)),
      donor_id = stringr::str_squish(as.character(donor_id))
    )
  obs <- attach_s1b_anatomy(obs)

  obs_use <- obs %>%
    dplyr::filter(
      .data[[scope_flag]] %in% TRUE,
      !is.na(anatomy_id),
      anatomy_group != "",
      anatomy_group != "Unmapped"
    ) %>%
    attach_s1b_jorstad_ei_class()

  class_levels <- tibble::tibble(EI_class = c("E", "I"))
  donor_region <- obs_use %>%
    dplyr::distinct(anatomy_group, donor_id)
  counts <- obs_use %>%
    dplyr::count(anatomy_group, donor_id, EI_class, name = "n_cells")

  long <- donor_region %>%
    tidyr::crossing(class_levels) %>%
    dplyr::left_join(
      counts,
      by = c("anatomy_group", "donor_id", "EI_class")
    ) %>%
    dplyr::mutate(n_cells = tidyr::replace_na(n_cells, 0)) %>%
    dplyr::group_by(anatomy_group, donor_id) %>%
    dplyr::mutate(
      donor_total_cells = sum(n_cells),
      donor_prop = dplyr::if_else(
        donor_total_cells > 0,
        n_cells / donor_total_cells,
        NA_real_
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(anatomy_group, EI_class) %>%
    dplyr::summarise(
      proportion = mean(donor_prop, na.rm = TRUE),
      n_donors = dplyr::n_distinct(donor_id[donor_total_cells > 0]),
      n_cells = sum(n_cells, na.rm = TRUE),
      .groups = "drop"
    )

  scope_lookup <- obs_use %>%
    dplyr::distinct(
      anatomy_id, anatomy_group,
      is_neocortex, is_cerebral_cortex, is_telencephalon
    )
  if (anyDuplicated(scope_lookup$anatomy_group)) {
    stop("An anatomy group has inconsistent s1b scope flags.", call. = FALSE)
  }

  out <- long %>%
    dplyr::select(anatomy_group, EI_class, proportion) %>%
    tidyr::pivot_wider(
      names_from = EI_class,
      values_from = proportion,
      values_fill = 0,
      names_prefix = "p_"
    ) %>%
    dplyr::left_join(
      long %>%
        dplyr::group_by(anatomy_group) %>%
        dplyr::summarise(
          n_donors = max(n_donors),
          n_cells_EI = sum(n_cells),
          .groups = "drop"
        ),
      by = "anatomy_group"
    ) %>%
    join_heiss_rates() %>%
    dplyr::left_join(
      scope_lookup,
      by = c("anatomy_id", "anatomy_group")
    ) %>%
    dplyr::mutate(
      EI_ratio = dplyr::if_else(p_I > 0, p_E / p_I, NA_real_),
      cell_class_definition = "Neocortical E:I cell-class crosswalk",
      anatomical_analysis_scope = scope_label,
      interpretation = interpretation
    ) %>%
    dplyr::filter(is.finite(EI_ratio), is.finite(rcmr_value)) %>%
    dplyr::arrange(anatomy_group)

  out
}

calculate_s1b_cerebral_cortex_class_ei <- function(
  obs_path = "data_intermediate/linnarsson_adult_human_brain_obs_metadata_neuronal.rds"
) {
  if (!exists("attach_s1b_anatomy", mode = "function")) {
    stop("Source helpers/read_heiss_rates.R before this helper.", call. = FALSE)
  }

  obs <- readRDS(obs_path) %>%
    dplyr::filter(cell_type == "neuron") %>%
    dplyr::mutate(
      roi = stringr::str_squish(as.character(roi)),
      supercluster_term = stringr::str_squish(as.character(supercluster_term)),
      donor_id = stringr::str_squish(as.character(donor_id))
    )
  obs <- attach_s1b_anatomy(obs)

  obs_use <- obs %>%
    dplyr::filter(
      is_cerebral_cortex %in% TRUE,
      !is.na(anatomy_id),
      anatomy_group != "",
      anatomy_group != "Unmapped"
    ) %>%
    dplyr::mutate(
      EI_class = dplyr::case_when(
        supercluster_term == "Miscellaneous" &
          as.character(cluster_id) %in% s1b_siletti_l5et_cluster_ids() ~ "E",
        supercluster_term %in% s1b_siletti_cerebral_cortex_e_terms() ~ "E",
        supercluster_term %in% s1b_siletti_cerebral_cortex_i_terms() ~
          "I_interneuron",
        supercluster_term %in% s1b_siletti_broad_i_msn_terms() ~ "I_MSN",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(EI_class))

  # Calculate all three classes together. This preserves the donor-balancing
  # used by the original s1b_5 analysis before the two I denominators diverge.
  class_levels <- tibble::tibble(
    EI_class = c("E", "I_interneuron", "I_MSN")
  )
  donor_region <- obs_use %>%
    dplyr::distinct(anatomy_group, donor_id)
  counts <- obs_use %>%
    dplyr::count(anatomy_group, donor_id, EI_class, name = "n_cells")

  long <- donor_region %>%
    tidyr::crossing(class_levels) %>%
    dplyr::left_join(
      counts,
      by = c("anatomy_group", "donor_id", "EI_class")
    ) %>%
    dplyr::mutate(n_cells = tidyr::replace_na(n_cells, 0)) %>%
    dplyr::group_by(anatomy_group, donor_id) %>%
    dplyr::mutate(
      donor_total_cells = sum(n_cells),
      donor_prop = dplyr::if_else(
        donor_total_cells > 0,
        n_cells / donor_total_cells,
        NA_real_
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(anatomy_group, EI_class) %>%
    dplyr::summarise(
      proportion = mean(donor_prop, na.rm = TRUE),
      n_donors = dplyr::n_distinct(donor_id[donor_total_cells > 0]),
      n_cells = sum(n_cells, na.rm = TRUE),
      .groups = "drop"
    )

  scope_lookup <- obs_use %>%
    dplyr::distinct(
      anatomy_id, anatomy_group,
      is_neocortex, is_cerebral_cortex, is_telencephalon
    )
  if (anyDuplicated(scope_lookup$anatomy_group)) {
    stop("An anatomy group has inconsistent s1b scope flags.", call. = FALSE)
  }

  wide <- long %>%
    dplyr::select(anatomy_group, EI_class, proportion) %>%
    tidyr::pivot_wider(
      names_from = EI_class,
      values_from = proportion,
      values_fill = 0,
      names_prefix = "p_"
    ) %>%
    dplyr::left_join(
      long %>%
        dplyr::group_by(anatomy_group) %>%
        dplyr::summarise(
          n_donors = max(n_donors),
          n_cells_EI = sum(n_cells),
          .groups = "drop"
        ),
      by = "anatomy_group"
    ) %>%
    join_heiss_rates() %>%
    dplyr::left_join(
      scope_lookup,
      by = c("anatomy_id", "anatomy_group")
    ) %>%
    dplyr::mutate(
      p_I_with_MSN = p_I_interneuron + p_I_MSN,
      EI_ratio_no_MSN = dplyr::if_else(
        p_I_interneuron > 0,
        p_E / p_I_interneuron,
        NA_real_
      ),
      EI_ratio_with_MSN = dplyr::if_else(
        p_I_with_MSN > 0,
        p_E / p_I_with_MSN,
        NA_real_
      )
    )

  no_msn <- wide %>%
    dplyr::mutate(
      p_I = p_I_interneuron,
      EI_ratio = EI_ratio_no_MSN,
      MSN_included_in_denominator = FALSE,
      cell_class_definition =
        "Siletti cerebral-cortex E:I classes (off-scope classes and MSN excluded)",
      interpretation =
        "Cerebral-cortex Siletti classes within cerebral-cortex scope"
    )
  with_msn <- wide %>%
    dplyr::mutate(
      p_I = p_I_with_MSN,
      EI_ratio = EI_ratio_with_MSN,
      MSN_included_in_denominator = TRUE,
      cell_class_definition =
        "Siletti cerebral-cortex E:I classes (off-scope classes excluded; MSN included)",
      interpretation =
        "Cerebral-cortex Siletti classes within cerebral-cortex scope; MSN sensitivity"
    )

  dplyr::bind_rows(no_msn, with_msn) %>%
    dplyr::mutate(
      anatomical_analysis_scope =
        "Cerebral cortex (including amygdaloid complex and hippocampus)"
    ) %>%
    dplyr::filter(is.finite(EI_ratio), is.finite(rcmr_value)) %>%
    dplyr::arrange(MSN_included_in_denominator, anatomy_group)
}

summarize_s1b_ei <- function(data) {
  fit <- stats::lm(rcmr_value ~ EI_ratio, data = data)
  fit_summary <- summary(fit)
  spearman <- suppressWarnings(stats::cor.test(
    data$EI_ratio,
    data$rcmr_value,
    method = "spearman",
    exact = FALSE
  ))

  tibble::tibble(
    cell_class_definition = data$cell_class_definition[[1]],
    anatomical_analysis_scope = data$anatomical_analysis_scope[[1]],
    interpretation = data$interpretation[[1]],
    n_regions = nrow(data),
    lm_beta = unname(stats::coef(fit)[["EI_ratio"]]),
    lm_p_value = stats::coef(fit_summary)["EI_ratio", "Pr(>|t|)"],
    lm_r_squared = fit_summary$r.squared,
    spearman_rho = unname(spearman$estimate),
    spearman_p_value = spearman$p.value
  )
}

summarize_s1b_neocortical_class_ei <- function(data) {
  summarize_s1b_ei(data)
}

plot_s1b_ei <- function(
  data,
  stats,
  title,
  subtitle,
  x_label = "E:I ratio"
) {
  check_region_palette(data, region_col = "anatomy_group")
  plot_data <- set_region_order(data, region_col = "anatomy_group")
  stats_label <- sprintf(
    "n = %d regions; linear model R2 = %.2f, p = %.3g; Spearman rho = %.2f, p = %.3g",
    stats$n_regions,
    stats$lm_r_squared,
    stats$lm_p_value,
    stats$spearman_rho,
    stats$spearman_p_value
  )

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = EI_ratio, y = rcmr_value, color = anatomy_group)
  ) +
    ggplot2::geom_smooth(
      ggplot2::aes(group = 1),
      method = "lm",
      se = TRUE,
      color = "grey35",
      fill = "grey85",
      linewidth = 0.8
    ) +
    ggplot2::geom_point(size = 3.2) +
    scale_color_regions(guide = "none") +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      caption = stats_label,
      x = x_label,
      y = "rCMRGlc (umol/100 g/min.)"
    ) +
    theme_project(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.caption = ggplot2::element_text(hjust = 0)
    )

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      ggplot2::aes(label = anatomy_group),
      show.legend = FALSE,
      size = 3.2,
      max.overlaps = Inf,
      min.segment.length = 0,
      box.padding = 0.35
    )
  }
  p
}

plot_s1b_neocortical_class_ei <- function(data, stats, title, subtitle) {
  plot_s1b_ei(
    data = data,
    stats = stats,
    title = title,
    subtitle = subtitle,
    x_label = "E:I ratio from neocortical cell-class crosswalk"
  )
}

run_s1b_neocortical_class_ei_scope <- function(
  scope_flag,
  scope_label,
  interpretation,
  title,
  subtitle,
  output_stem
) {
  dir.create("data_analysis", showWarnings = FALSE, recursive = TRUE)
  dir.create("figs/s1b", showWarnings = FALSE, recursive = TRUE)

  data <- calculate_s1b_neocortical_class_ei(
    scope_flag = scope_flag,
    scope_label = scope_label,
    interpretation = interpretation
  )
  stats <- summarize_s1b_neocortical_class_ei(data)
  plot <- plot_s1b_neocortical_class_ei(data, stats, title, subtitle)

  utils::write.csv(data, paste0("data_analysis/", output_stem, ".csv"), row.names = FALSE)
  utils::write.csv(stats, paste0("data_analysis/", output_stem, "_statistics.csv"), row.names = FALSE)
  ggplot2::ggsave(
    paste0("figs/s1b/", output_stem, ".pdf"),
    plot,
    width = 8.5,
    height = 6.8,
    units = "in"
  )
  ggplot2::ggsave(
    paste0("figs/s1b/", output_stem, ".jpg"),
    plot,
    width = 8.5,
    height = 6.8,
    units = "in",
    dpi = 300
  )

  print(stats)
  invisible(list(data = data, stats = stats, plot = plot))
}
