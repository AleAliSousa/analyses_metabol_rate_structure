# =====================================================================
# 0_prepare_heiss_rates.R
#
# Prepare the 26 published Heiss et al. (2004) Table 1 region rows.
# This script deliberately contains no Stephan-derived volume terminology,
# comparative-volume crosswalks, or calculated composite rates.
# =====================================================================

# Locate the repository from either Rscript or an interactive/source call.
.script_args <- commandArgs(trailingOnly = FALSE)
.file_arg <- .script_args[grepl("^--file=", .script_args)]
.start_dir <- if (length(.file_arg)) {
  # Rscript can encode spaces as "~+~" in commandArgs() on macOS.
  .script_path <- sub("^--file=", "", .file_arg[1])
  .script_path <- gsub("~+~", " ", .script_path, fixed = TRUE)
  dirname(normalizePath(.script_path, mustWork = FALSE))
} else {
  getwd()
}

.repo_root <- normalizePath(.start_dir, mustWork = FALSE)
while (!file.exists(file.path(.repo_root, ".git")) &&
       dirname(.repo_root) != .repo_root) {
  .repo_root <- dirname(.repo_root)
}
if (!file.exists(file.path(.repo_root, ".git"))) {
  stop("Could not locate the project root from: ", .start_dir, call. = FALSE)
}
setwd(.repo_root)

raw_path <- "data_raw/Heiss_etal_2004_TABLE1.csv"
registry_path <- "metadata/anatomy/anatomy_registry.csv"
map_path <- "metadata/anatomy/heiss_2004_region_map.csv"
output_path <- "data_intermediate/heiss_2004_regions.csv"

required_files <- c(raw_path, registry_path, map_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop(
    "Missing required input file(s): ",
    paste(missing_files, collapse = ", "),
    call. = FALSE
  )
}

read_project_csv <- function(path, check.names = TRUE) {
  read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = check.names,
    na.strings = c("", "NA")
  )
}

require_columns <- function(data, columns, object_name) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) {
    stop(
      object_name, " is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

assert_unique_nonmissing <- function(x, label) {
  if (anyNA(x) || any(trimws(x) == "")) {
    stop(label, " contains missing or blank values.", call. = FALSE)
  }
  duplicates <- unique(x[duplicated(x)])
  if (length(duplicates)) {
    stop(
      label, " contains duplicate value(s): ",
      paste(duplicates, collapse = ", "),
      call. = FALSE
    )
  }
}

registry <- read_project_csv(registry_path)
require_columns(
  registry,
  c(
    "anatomy_id", "canonical_label", "parent_id", "anatomy_level",
    "tissue_class", "color_hex", "display_order", "definition_note"
  ),
  registry_path
)

assert_unique_nonmissing(registry$anatomy_id, "anatomy_registry$anatomy_id")
assert_unique_nonmissing(
  registry$canonical_label,
  "anatomy_registry$canonical_label"
)

bad_ids <- registry$anatomy_id[
  !grepl("^[a-z][a-z0-9_]*$", registry$anatomy_id)
]
if (length(bad_ids)) {
  stop(
    "Registry anatomy_id values must be lower snake case: ",
    paste(bad_ids, collapse = ", "),
    call. = FALSE
  )
}

parent_ids <- registry$parent_id[!is.na(registry$parent_id)]
missing_parents <- setdiff(parent_ids, registry$anatomy_id)
if (length(missing_parents)) {
  stop(
    "Registry parent_id value(s) are not defined: ",
    paste(missing_parents, collapse = ", "),
    call. = FALSE
  )
}

bad_colors <- registry$color_hex[
  is.na(registry$color_hex) |
    !grepl("^#[0-9A-Fa-f]{6}$", registry$color_hex)
]
if (length(bad_colors)) {
  stop(
    "Registry color_hex values must use #RRGGBB format: ",
    paste(unique(bad_colors), collapse = ", "),
    call. = FALSE
  )
}

if (anyNA(registry$display_order) ||
    anyDuplicated(registry$display_order)) {
  stop("Registry display_order values must be present and unique.", call. = FALSE)
}

heiss_map <- read_project_csv(map_path)
require_columns(
  heiss_map,
  c("heiss_region", "anatomy_id", "relationship", "mapping_note"),
  map_path
)
assert_unique_nonmissing(heiss_map$heiss_region, "heiss_map$heiss_region")

unknown_ids <- setdiff(heiss_map$anatomy_id, registry$anatomy_id)
if (length(unknown_ids)) {
  stop(
    "Heiss mapping uses anatomy_id value(s) absent from the registry: ",
    paste(unknown_ids, collapse = ", "),
    call. = FALSE
  )
}

allowed_relationships <- c("exact", "source_synonym", "source_measurement")
bad_relationships <- setdiff(unique(heiss_map$relationship), allowed_relationships)
if (length(bad_relationships)) {
  stop(
    "Unsupported Heiss mapping relationship(s): ",
    paste(bad_relationships, collapse = ", "),
    ". Allowed values: ", paste(allowed_relationships, collapse = ", "),
    call. = FALSE
  )
}

heiss_raw <- read_project_csv(raw_path, check.names = FALSE)
raw_columns <- c(
  "category",
  "Region",
  "Both hemispheres Mean",
  "Both hemispheres SD",
  "Left minus right hemisphere Difference",
  "Left minus right hemisphere SD",
  "P",
  "Data of Heiss et al. (21)"
)
require_columns(heiss_raw, raw_columns, raw_path)
assert_unique_nonmissing(heiss_raw$Region, "Heiss raw Region")

unmapped_regions <- setdiff(heiss_raw$Region, heiss_map$heiss_region)
orphan_map_rows <- setdiff(heiss_map$heiss_region, heiss_raw$Region)
if (length(unmapped_regions) || length(orphan_map_rows)) {
  details <- c(
    if (length(unmapped_regions)) {
      paste0("Unmapped Heiss region(s): ", paste(unmapped_regions, collapse = ", "))
    },
    if (length(orphan_map_rows)) {
      paste0("Mapping row(s) absent from raw Heiss data: ",
             paste(orphan_map_rows, collapse = ", "))
    }
  )
  stop(paste(details, collapse = "\n"), call. = FALSE)
}

# Preserve the published Heiss row order. No aggregation or new rows occur here.
map_index <- match(heiss_raw$Region, heiss_map$heiss_region)
registry_index <- match(heiss_map$anatomy_id[map_index], registry$anatomy_id)

heiss_prepared <- data.frame(
  source_row = seq_len(nrow(heiss_raw)),
  source_dataset = "Heiss et al. 2004 Table 1",
  heiss_category = heiss_raw$category,
  heiss_region = heiss_raw$Region,
  anatomy_id = heiss_map$anatomy_id[map_index],
  canonical_label = registry$canonical_label[registry_index],
  parent_id = registry$parent_id[registry_index],
  mapping_relationship = heiss_map$relationship[map_index],
  rcmrglc_mean_umol_100g_min = as.numeric(heiss_raw[["Both hemispheres Mean"]]),
  rcmrglc_sd_umol_100g_min = as.numeric(heiss_raw[["Both hemispheres SD"]]),
  left_minus_right_mean = as.numeric(
    heiss_raw[["Left minus right hemisphere Difference"]]
  ),
  left_minus_right_sd = as.numeric(
    heiss_raw[["Left minus right hemisphere SD"]]
  ),
  left_minus_right_p_value = as.numeric(heiss_raw[["P"]]),
  heiss_1991_mean_umol_100g_min = as.numeric(
    heiss_raw[["Data of Heiss et al. (21)"]]
  ),
  color_hex = registry$color_hex[registry_index],
  display_order = registry$display_order[registry_index],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (nrow(heiss_prepared) != nrow(heiss_raw)) {
  stop("Preparation changed the number of Heiss rows.", call. = FALSE)
}
if (anyNA(heiss_prepared$anatomy_id) || anyNA(heiss_prepared$canonical_label)) {
  stop("Preparation produced an unmapped anatomy row.", call. = FALSE)
}
if (anyNA(heiss_prepared$rcmrglc_mean_umol_100g_min) ||
    any(heiss_prepared$rcmrglc_mean_umol_100g_min <= 0)) {
  stop("Every Heiss row must have a positive mean rCMRGlc value.", call. = FALSE)
}
if (anyNA(heiss_prepared$rcmrglc_sd_umol_100g_min) ||
    any(heiss_prepared$rcmrglc_sd_umol_100g_min < 0)) {
  stop("Every Heiss row must have a non-negative rCMRGlc SD.", call. = FALSE)
}

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write.csv(heiss_prepared, output_path, row.names = FALSE, na = "")

message(
  "Prepared ", nrow(heiss_prepared),
  " published Heiss region rows: ", output_path
)
