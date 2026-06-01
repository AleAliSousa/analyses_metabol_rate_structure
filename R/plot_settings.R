############################################################
# Shared plotting settings
# Source this file at the top of plotting scripts:
#   source("R/01_plot_settings.R")
############################################################

# -----------------------------
# Region order
# -----------------------------
# Keep this as the single source of truth for region display order.
region_order <- c(
  "Cerebral cortex (global average)",
  "Frontal lobe",
  "Parietal lobe",
  "Temporal lobe",
  "Occipital lobe",
  "Insular lobe",
  "Hippocampus",
  "Corpus amygdaloideum",
  "Basal forebrain",
  "Caudatum",
  "Putamen",
  "Nucleus accumbens",
  "Pallidum",
  "Nucleus subthalamicus",
  "Nucleus medial thalami",
  "Corpus geniculatum laterale",
  "Corpus geniculatum mediale",
  "Substantia nigra",
  "Nucleus ruber",
  "Colliculus superior",
  "Colliculus inferior",
  "Cerebellar cortex",
  "Vermis",
  "Nucleus dentatus cerebelli",
  "Centrum semiovale",
  "Capsula interna"
)

# -----------------------------
# Region color palette
# -----------------------------
# Named vector: colors remain stable even when a plot contains only a subset
# of regions. Add new regions here rather than defining palettes inside plots.
region_palette <- c(
  "Cerebral cortex (global average)" = "#000000",
  "Frontal lobe"                    = "#56B4E9",
  "Parietal lobe"                   = "#D55E00",
  "Temporal lobe"                   = "#CC79A7",
  "Occipital lobe"                  = "#0072B2",
  "Insular lobe"                    = "#009E73",
  "Hippocampus"                     = "#882255",
  "Corpus amygdaloideum"            = "#E69F00",
  "Basal forebrain"                 = "#AA4499",
  "Caudatum"                        = "#44AA99",
  "Putamen"                         = "#117733",
  "Nucleus accumbens"               = "#88CCEE",
  "Pallidum"                        = "#999933",
  "Nucleus subthalamicus"           = "#332288",
  "Nucleus medial thalami"          = "#DDCC77",
  "Corpus geniculatum laterale"     = "#6699CC",
  "Corpus geniculatum mediale"      = "#661100",
  "Substantia nigra"                = "#A6761D",
  "Nucleus ruber"                   = "#E7298A",
  "Colliculus superior"             = "#66A61E",
  "Colliculus inferior"             = "#7570B3",
  "Cerebellar cortex"               = "#1B9E77",
  "Vermis"                          = "#B2DF8A",
  "Nucleus dentatus cerebelli"      = "#33A02C",
  "Centrum semiovale"               = "#BDBDBD",
  "Capsula interna"                 = "#636363"
)

stopifnot(all(region_order %in% names(region_palette)))
stopifnot(!anyDuplicated(names(region_palette)))

# -----------------------------
# Helper functions
# -----------------------------

# Check that every region in a data frame has a defined color.
check_region_palette <- function(data, region_col = "anatomy_group", palette = region_palette) {
  regions <- sort(unique(as.character(data[[region_col]])))
  regions <- regions[!is.na(regions)]
  missing_regions <- setdiff(regions, names(palette))

  if (length(missing_regions) > 0) {
    stop(
      "Missing colors for regions: ",
      paste(missing_regions, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

# Apply the shared region order to a plotting data frame.
set_region_order <- function(data, region_col = "anatomy_group") {
  data[[region_col]] <- factor(
    data[[region_col]],
    levels = region_order
  )
  data
}

# ggplot2 scales using the shared region palette.
scale_color_regions <- function(...) {
  ggplot2::scale_color_manual(values = region_palette, drop = FALSE, ...)
}

scale_fill_regions <- function(...) {
  ggplot2::scale_fill_manual(values = region_palette, drop = FALSE, ...)
}

# Shared base theme for figures.
theme_project <- function(base_size = 14) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = NA, color = "black"),
      strip.text = ggplot2::element_text(size = base_size * 0.85),
      legend.title = ggplot2::element_text(size = base_size * 0.9),
      legend.text = ggplot2::element_text(size = base_size * 0.8)
    )
}
