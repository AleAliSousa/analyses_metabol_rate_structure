# =====================================================================
# ARCHIVED ONE-TIME MIGRATION
#
# Split comparison measurements out of the legacy mixed tables:
#   data_raw/stereology.csv
#   data_raw/stress_volume.csv
#
# Outputs (now version-controlled publication inputs):
#   data_raw/stereology_comparison.csv
#   data_raw/stress_volume_comparison.csv
#
# Active analyses do not run this script. It is retained only to document and
# reproduce the one-time separation. The active preparation workflow validates
# and copies the clean comparison inputs to data_intermediate/.
# =====================================================================

setwd(local({
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d)
  d
}))

read_project_csv <- function(path) {
  read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA")
  )
}

heiss_raw <- read_project_csv("data_raw/Heiss_etal_2004_TABLE1.csv")
heiss_map <- read_project_csv("metadata/anatomy/heiss_2004_region_map.csv")
heiss_id <- heiss_map$anatomy_id[match(heiss_raw$Region, heiss_map$heiss_region)]
rate_by_id <- setNames(heiss_raw$`Both hemispheres Mean`, heiss_id)

# ---- s1a stereology -------------------------------------------------
stereology_raw <- read_project_csv("data_raw/stereology.csv")
measurement_columns <- setdiff(names(stereology_raw), c("Region", "rCMRGlc"))
keep <- rowSums(!is.na(stereology_raw[measurement_columns])) > 0
stereology <- stereology_raw[keep, , drop = FALSE]
stereology$Region <- trimws(stereology$Region)
stereology_id <- heiss_map$anatomy_id[
  match(stereology$Region, heiss_map$heiss_region)
]

if (anyNA(stereology_id) ||
    any(abs(stereology$rCMRGlc - rate_by_id[stereology_id]) > 1e-8)) {
  stop("Legacy stereology anatomy/rate values do not match Heiss Table 1.")
}

stereology_out <- data.frame(
  source_region = stereology$Region,
  stereology[measurement_columns],
  check.names = FALSE
)
write.csv(
  stereology_out,
  "data_raw/stereology_comparison.csv",
  row.names = FALSE,
  na = ""
)

# ---- s2 stress volume -----------------------------------------------
stress_raw <- read_project_csv("data_raw/stress_volume.csv")
stress_raw$anatomy_group <- trimws(stress_raw$anatomy_group)
stress_map <- read_project_csv(
  "metadata/anatomy/s2_stress_volume_region_map.csv"
)

source_terms <- unique(stress_map$source_region)
resolved <- lapply(source_terms, function(term) {
  rows <- stress_map[stress_map$source_region == term, , drop = FALSE]
  components <- rows[!is.na(rows$heiss_anatomy_id), , drop = FALSE]
  rate <- if (!nrow(components)) {
    NA_real_
  } else {
    weights <- components$component_weight / sum(components$component_weight)
    sum(rate_by_id[components$heiss_anatomy_id] * weights)
  }
  data.frame(
    source_region = term,
    anatomy_id = unique(rows$analysis_anatomy_id),
    rcmr_value = rate,
    stringsAsFactors = FALSE
  )
})
resolved <- do.call(rbind, resolved)
idx <- match(stress_raw$anatomy_group, resolved$source_region)

legacy_rate <- stress_raw$rcmr_value
central_rate <- resolved$rcmr_value[idx]
mismatch <- xor(is.na(legacy_rate), is.na(central_rate)) |
  (!is.na(legacy_rate) & abs(legacy_rate - central_rate) > 0.05)
if (anyNA(idx) || any(mismatch)) {
  stop("Legacy stress-volume anatomy/rate values do not match the s2 map.")
}

stress_columns <- setdiff(
  names(stress_raw), c("anatomy_group", "rcmr_value")
)
stress_out <- data.frame(
  source_region = stress_raw$anatomy_group,
  stress_raw[stress_columns],
  check.names = FALSE
)
write.csv(
  stress_out,
  "data_raw/stress_volume_comparison.csv",
  row.names = FALSE,
  na = ""
)

message(
  "Archived migration reproduced ", nrow(stereology_out),
  " stereology and ", nrow(stress_out), " stress-volume comparison rows."
)
