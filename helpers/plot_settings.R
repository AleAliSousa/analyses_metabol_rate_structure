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
  "Capsula interna"                 = "#636363",

  # Coarser regions used in Study 2 / Study 3 that have no single fine-grained
  # equivalent above. Given their own stable colors so they are consistent
  # across those studies' figures. (Regions that DO have a fine equivalent are
  # handled by region_aliases below rather than duplicated here.)
  "Thalamus"                        = "#117777",
  "Striatum"                        = "#7B3294",
  "Neocortex white"                 = "#878787",
  "Mesencephalon"                   = "#B15928",
  "Cerebellum"                      = "#1F78B4"
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

# -----------------------------
# Faceted multi-panel figure helpers
# -----------------------------
# These keep multi-panel (facet_wrap / facet_grid) figures from being squashed
# into a short, wide canvas. Use facet_ncol() to set the panel grid, the
# compact theme to shrink text, and facet_dims() to size the saved file so each
# panel stays tall enough to read.

# Sensible number of columns for n facet panels.
facet_ncol <- function(n) {
  if (n <= 3) n else if (n <= 8) 3L else 4L
}

# Compact theme for faceted scatter grids: smaller strip/axis text and more
# space between panels so equation labels and axis numbers do not overlap.
theme_facet_compact <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      panel.spacing    = ggplot2::unit(1.2, "lines"),
      strip.background = ggplot2::element_rect(fill = NA, color = "black"),
      strip.text       = ggplot2::element_text(size = base_size * 0.80),
      axis.text        = ggplot2::element_text(size = base_size * 0.70),
      plot.title       = ggplot2::element_text(size = base_size * 1.05, face = "bold"),
      plot.subtitle    = ggplot2::element_text(size = base_size * 0.85),
      plot.caption     = ggplot2::element_text(size = base_size * 0.70, hjust = 0)
    )
}

# Recommended saved size (inches) for an n-panel facet grid. Width grows with
# the number of columns, height with the number of rows, so panels keep a
# readable aspect ratio instead of being compressed. Returns ncol/width/height.
facet_dims <- function(n, ncol = facet_ncol(n),
                       panel_w = 3.3, panel_h = 3.0,
                       legend_w = 2.6, title_h = 1.4) {
  nrow <- ceiling(n / ncol)
  list(
    ncol   = ncol,
    nrow   = nrow,
    width  = ncol * panel_w + legend_w,
    height = nrow * panel_h + title_h
  )
}

# -----------------------------
# Cross-study region label aliases
# -----------------------------
# Different studies label the same brain region differently (Study 1 uses the
# Heiss/Stephan terms that key region_palette; Studies 2 and 3 use coarser or
# English names). Map those alternates onto the canonical palette key so a given
# region gets the SAME colour in every figure. Each entry asserts that the two
# labels denote the same region -- review before trusting cross-study colours.
region_aliases <- c(
  # Study 2 (stress_volume) coarse English labels
  "Cerebral cortex"             = "Cerebral cortex (global average)",
  "Caudate"                     = "Caudatum",
  "Amygdala"                    = "Corpus amygdaloideum",
  "Accumbens"                   = "Nucleus accumbens",
  # Study 3 (PGLS) labels
  "Neocortex grey"              = "Cerebral cortex (global average)",
  "Insular cortex (grey)"       = "Insular lobe",
  "Nucleus subthalamicus Luysi" = "Nucleus subthalamicus",
  "Area striata grey"           = "Occipital lobe"   # area striata = V1 (occipital)
)

# Normalise region labels: trim whitespace and apply region_aliases so they line
# up with region_palette / region_order. Base R only (safe before tidyverse is
# loaded). Vectorised.
canonical_region <- function(x) {
  x <- trimws(as.character(x))
  hit <- !is.na(x) & x %in% names(region_aliases)
  x[hit] <- region_aliases[x[hit]]
  x
}

# Colour vector for a set of region labels given in display order, for colouring
# discrete axis tick labels (e.g. Study 3 plots where region is the x/y axis
# rather than a colour aesthetic). Unmatched labels fall back to `default`.
region_axis_colors <- function(levels_in_order, default = "grey30") {
  canon <- canonical_region(levels_in_order)
  cols <- unname(region_palette[canon])
  cols[is.na(cols)] <- default
  cols
}
