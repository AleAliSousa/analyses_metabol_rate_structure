# Shared functions for the mouse s1b scripts (loaded after common.R).
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr); library(ggplot2)})

MIN_CELLS_DONOR_ROI <- 100   # donor x ROI units with fewer cells of the relevant universe are dropped

read_mouse_counts <- function() {
  x <- read_csv(nh("data_intermediate/allen_wmb10x_counts_roi_donor_celltype.csv"),
                show_col_types = FALSE, na = character())
  rename(x, roi = region_of_interest_acronym, donor = donor_label)
}

read_roi_map <- function() {
  m <- read_csv(nh("metadata/mouse_roi_to_cmrgl_region_map.csv"), show_col_types = FALSE, na = c("", "NA"))
  stopifnot(!anyDuplicated(m$allen_roi))
  m
}

read_mouse_rates <- function() {
  r <- read_csv(nh("data_raw/Evo_M1_mouse_CMRgl_by_region.csv"), show_col_types = FALSE)
  stopifnot(!anyDuplicated(r$Region))
  select(r, cmrgl_region = Region, rcmr_value = Value, n_rate_studies = n_studies,
         rate_references = Primary_reference_ids)
}

# Mean donor-level proportions per ROI (zeros filled), universe = rows of `x`
# (already filtered), categories in column `cat_col`.
roi_donor_props <- function(x, cat_col, min_cells = MIN_CELLS_DONOR_ROI) {
  x <- rename(x, category = all_of(cat_col))
  cats <- distinct(x, category)
  units <- x %>% group_by(roi, donor) %>% summarise(total = sum(n_cells), .groups = "drop") %>%
    filter(total >= min_cells)
  x %>% group_by(roi, donor, category) %>% summarise(n = sum(n_cells), .groups = "drop") %>%
    right_join(crossing(units, cats), by = c("roi", "donor", "category")) %>%
    mutate(n = replace_na(n, 0), prop = n / total)
}

# Average to Allen ROI, then to CMRgl region (equal weight per ROI), keep chosen ROIs
region_props <- function(donor_props, roi_map, primary_only = TRUE) {
  m <- filter(roi_map, !is.na(cmrgl_region), if (primary_only) primary_included else TRUE)
  donor_props %>%
    group_by(roi, category) %>%
    summarise(prop = mean(prop), n_donors = n_distinct(donor), n_cells = sum(n), .groups = "drop") %>%
    inner_join(select(m, roi = allen_roi, cmrgl_region), by = "roi") %>%
    group_by(cmrgl_region, category) %>%
    summarise(prop = mean(prop), n_donors = sum(n_donors), n_cells = sum(n_cells),
              rois = paste(sort(unique(roi)), collapse = "+"), .groups = "drop")
}

spearman_table <- function(df, predictors, response = "rcmr_value") {
  bind_rows(lapply(predictors, function(v) {
    ok <- complete.cases(df[[response]], df[[v]])
    if (sum(ok) < 4 || sd(df[[v]][ok]) == 0)
      return(data.frame(predictor = v, rho = NA_real_, p_value = NA_real_, n = sum(ok)))
    ct <- suppressWarnings(cor.test(df[[response]][ok], df[[v]][ok], method = "spearman", exact = FALSE))
    data.frame(predictor = v, rho = unname(ct$estimate), p_value = ct$p.value, n = sum(ok))
  })) %>% mutate(p_adj_BH = p.adjust(p_value, "BH")) %>% arrange(p_value)
}

# Facetted scatter: proportion (x) vs rate (y), one panel per predictor
plot_vs_rate <- function(df, predictors, stats, file, xlab, ncol = 3) {
  pd <- df %>% select(cmrgl_region, rcmr_value, all_of(predictors)) %>%
    pivot_longer(all_of(predictors), names_to = "predictor", values_to = "prop") %>%
    left_join(stats, by = "predictor") %>%
    mutate(panel = sprintf("%s\nrho=%.2f, p=%.3f, n=%d", predictor, rho, p_value, n))
  p <- ggplot(pd, aes(prop, rcmr_value)) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "steelblue", linewidth = 0.6) +
    geom_point(aes(colour = cmrgl_region), size = 2.4, alpha = 0.9) +
    facet_wrap(~panel, scales = "free_x", ncol = ncol) +
    labs(x = xlab, y = "Mouse rCMRgl (µmol/100 g/min)", colour = "Region") +
    theme_bw(base_size = 11) + theme(strip.text = element_text(size = 8))
  nrow_f <- ceiling(length(predictors) / ncol)
  ggsave(paste0(file, ".pdf"), p, width = 3.2 * ncol + 1.5, height = 2.8 * nrow_f)
  ggsave(paste0(file, ".jpg"), p, width = 3.2 * ncol + 1.5, height = 2.8 * nrow_f, dpi = 200)
  invisible(p)
}

