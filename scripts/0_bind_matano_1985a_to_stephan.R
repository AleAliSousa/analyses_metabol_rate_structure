# Bind Matano 1985a cerebellar nuclei data onto Stephan_primates.csv.
# The join is by species.  Name changes/synonyms are handled with species_key.csv
# plus a few manual mappings for names used in the current Stephan_primates.csv.

setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))  # repo root (portable; replaces hardcoded path -- see R/project_root.R)

# ---------- Helpers ----------
first_existing <- function(paths, label) {
  checked_paths <- paths
  paths <- checked_paths[file.exists(checked_paths)]
  if (length(paths) == 0L) {
    stop("Could not find ", label, ". Checked: ", paste(checked_paths, collapse = ", "))
  }
  paths[[1L]]
}

trim_character_columns <- function(dat) {
  dat[] <- lapply(dat, function(x) {
    if (is.character(x)) {
      x <- trimws(x)
      x[x == ""] <- NA_character_
    }
    x
  })
  dat
}

read_csv_plain <- function(path, encoding = "UTF-8") {
  dat <- read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA", "NaN"),
    fileEncoding = encoding
  )
  trim_character_columns(dat)
}

parse_number_plain <- function(x) {
  # Handles values such as "6,400" and "1,166.90".
  suppressWarnings(as.numeric(gsub(",", "", as.character(x), fixed = TRUE)))
}

normalise_for_matching <- function(x) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- NA_character_
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", " ", x)
  x <- trimws(gsub("\\s+", " ", x))
  x[x == ""] <- NA_character_
  x
}

species_to_file_name <- function(x) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- NA_character_
  gsub("\\s+", "_", trimws(x))
}

coalesce_character <- function(x, y) {
  out <- x
  take_y <- is.na(out) | out == ""
  out[take_y] <- y[take_y]
  out
}

# ---------- Input paths ----------
stephan_path <- "data_raw/Stephan_primates.csv"
matano_path <- first_existing(
  c("data_raw/to_add/Matano_1985_a.csv", "data_raw/to_add/Matano_1985a.csv"),
  "Matano 1985a file"
)
species_key_path <- first_existing(
  c("data_raw/species_key.csv", "data_raw/to_add/species_key.csv", "species_key.csv"),
  "species_key.csv"
)

# ---------- Read data ----------
stephan <- read_csv_plain(stephan_path)
# Matano file contains a Latin-1 character in at least one non-primate species name.
matano_raw <- read_csv_plain(matano_path, encoding = "latin1")
species_key <- read_csv_plain(species_key_path)

stopifnot("Species" %in% names(stephan))
stopifnot("Species" %in% names(matano_raw))
stopifnot(all(c("accepted_name", "source_publication", "variant_name") %in% names(species_key)))

# ---------- Remove Matano metadata rows and parse numeric columns ----------
matano_source_cols <- c(
  "Species_Matano1985a",
  "code_Matano1985a",
  "Body_weight_1985a",
  "Cerebellar_nuclei_total",
  "Interpositus_cerebellar_nuclei",
  "Lateral_cerebellar_nuclei",
  "Medial_cerebellar_nuclei",
  "Number_cerebellar_nuclei"
)
matano_numeric_cols <- c(
  "Body_weight_1985a",
  "Cerebellar_nuclei_total",
  "Interpositus_cerebellar_nuclei",
  "Lateral_cerebellar_nuclei",
  "Medial_cerebellar_nuclei",
  "Number_cerebellar_nuclei"
)

missing_cols <- setdiff(matano_source_cols, names(matano_raw))
if (length(missing_cols) > 0L) {
  stop("Matano file is missing expected columns: ", paste(missing_cols, collapse = ", "))
}

matano <- matano_raw[!grepl("^AAAA_", matano_raw$Species), , drop = FALSE]

has_source_data <- Reduce(
  `|`,
  lapply(matano_source_cols, function(col) !is.na(matano[[col]]) & matano[[col]] != "")
)
matano <- matano[has_source_data, , drop = FALSE]

for (col in matano_numeric_cols) {
  matano[[col]] <- parse_number_plain(matano[[col]])
}
if ("Body_weight_current" %in% names(matano)) {
  matano[["Body_weight_current"]] <- parse_number_plain(matano[["Body_weight_current"]])
}

# ---------- Species-name matching ----------
# Use only Matano 1985a entries from the key.  The source labels sometimes appear
# with or without underscores, so strip punctuation before filtering.
species_key$source_publication_clean <- gsub(
  "[^a-z0-9]+", "", tolower(species_key$source_publication)
)
key_matano <- species_key[species_key$source_publication_clean == "matano1985a", , drop = FALSE]
key_matano$variant_clean <- normalise_for_matching(key_matano$variant_name)

# Check that duplicated Matano variant names do not map to conflicting accepted names.
key_conflicts <- aggregate(
  accepted_name ~ variant_clean,
  data = key_matano[!is.na(key_matano$variant_clean), c("variant_clean", "accepted_name")],
  FUN = function(z) paste(sort(unique(z)), collapse = " | ")
)
key_conflicts$n_accepted <- vapply(strsplit(key_conflicts$accepted_name, " \\| "), length, integer(1))
key_conflicts <- key_conflicts[key_conflicts$n_accepted > 1L, , drop = FALSE]
if (nrow(key_conflicts) > 0L) {
  stop("Conflicting species_key mappings for Matano1985a. Inspect species_key.csv before joining.")
}

lookup_names <- key_matano$variant_clean
lookup_values <- key_matano$accepted_name
keep_lookup <- !is.na(lookup_names) & !duplicated(lookup_names)
lookup_names <- lookup_names[keep_lookup]
lookup_values <- lookup_values[keep_lookup]

accepted_from_species <- lookup_values[match(normalise_for_matching(matano$Species), lookup_names)]
accepted_from_matano_name <- lookup_values[match(normalise_for_matching(matano$Species_Matano1985a), lookup_names)]

matano$matano_accepted_name <- coalesce_character(accepted_from_species, accepted_from_matano_name)
matano$matano_accepted_name <- coalesce_character(matano$matano_accepted_name, matano$Species)
matano$target_species_raw <- species_to_file_name(matano$matano_accepted_name)
matano$target_species <- matano$target_species_raw

# Manual mappings where the provided key uses older/lumped names, while the
# current Stephan file uses newer names or a slightly different spelling.
manual_name_map <- c(
  "Callicebus_moloch" = "Plecturocebus_moloch",
  "Lagothrix_lagothricha" = "Lagothrix_lagotricha",
  "Gorilla_sp." = "Gorilla_gorilla",
  "Gorilla_sp" = "Gorilla_gorilla",
  "Tarsius_sp." = "Tarsius_syrichta",
  "Tarsius_sp" = "Tarsius_syrichta"
)

not_direct <- !(matano$target_species %in% stephan$Species)
can_manual <- not_direct & matano$target_species %in% names(manual_name_map)
manual_targets <- unname(manual_name_map[matano$target_species[can_manual]])
# Only apply manual maps when the mapped name is actually present in Stephan_primates.csv.
manual_targets[!(manual_targets %in% stephan$Species)] <- NA_character_
idx_manual <- which(can_manual)
idx_manual <- idx_manual[!is.na(manual_targets)]
matano$target_species[idx_manual] <- manual_targets[!is.na(manual_targets)]

matano$matched_in_Stephan <- matano$target_species %in% stephan$Species

# ---------- Keep source columns to add ----------
matano_keep_cols <- c("target_species", matano_source_cols)
matano_for_join <- matano[matano$matched_in_Stephan, matano_keep_cols, drop = FALSE]

# Avoid accidental many-to-one joins.
duplicated_targets <- unique(matano_for_join$target_species[duplicated(matano_for_join$target_species)])
if (length(duplicated_targets) > 0L) {
  stop(
    "More than one Matano row maps to the same Stephan species: ",
    paste(duplicated_targets, collapse = ", "),
    ". Resolve duplicates before joining."
  )
}

# ---------- Unit checks ----------
# Body weights: Matano 1985a body weight is in the same gram scale as Stephan's
# Body_weight in the supplied files.  This check catches kg/g mismatches.
stephan_qc <- data.frame(
  Species = stephan$Species,
  Body_weight = if ("Body_weight" %in% names(stephan)) parse_number_plain(stephan$Body_weight) else NA_real_,
  Cerebellum = if ("Cerebellum" %in% names(stephan)) parse_number_plain(stephan$Cerebellum) else NA_real_,
  stringsAsFactors = FALSE
)
qc <- merge(
  matano_for_join,
  stephan_qc,
  by.x = "target_species",
  by.y = "Species",
  all.x = TRUE,
  sort = FALSE
)

qc$Body_weight_1985a_to_Stephan <- qc$Body_weight_1985a / qc$Body_weight
bw_ratio_median <- median(qc$Body_weight_1985a_to_Stephan, na.rm = TRUE)

if (is.finite(bw_ratio_median)) {
  if (bw_ratio_median > 500 && bw_ratio_median < 2000) {
    matano_for_join$Body_weight_1985a <- matano_for_join$Body_weight_1985a / 1000
    message("Converted Body_weight_1985a by /1000 to match Stephan Body_weight scale.")
  } else if (bw_ratio_median > 0.0005 && bw_ratio_median < 0.002) {
    matano_for_join$Body_weight_1985a <- matano_for_join$Body_weight_1985a * 1000
    message("Converted Body_weight_1985a by *1000 to match Stephan Body_weight scale.")
  } else if (bw_ratio_median < 0.1 || bw_ratio_median > 10) {
    warning(
      "Body_weight_1985a does not look like the same unit as Stephan Body_weight. ",
      "Median ratio = ", signif(bw_ratio_median, 4), ". Inspect before using."
    )
  } else {
    message("Body_weight_1985a is on the same scale as Stephan Body_weight; no conversion applied.")
  }
}

# Volume check: Matano cerebellar nuclei volumes should be much smaller than
# Stephan's total Cerebellum volume if both are in mm^3.  In the supplied files,
# median Cerebellar_nuclei_total/Cerebellum is about 0.022, so no cm/mm
# conversion is applied.
qc <- merge(
  matano_for_join,
  stephan_qc,
  by.x = "target_species",
  by.y = "Species",
  all.x = TRUE,
  sort = FALSE
)
qc$Cerebellar_nuclei_total_to_Cerebellum <- qc$Cerebellar_nuclei_total / qc$Cerebellum
vol_ratio_median <- median(qc$Cerebellar_nuclei_total_to_Cerebellum, na.rm = TRUE)
vol_ratio_max <- max(qc$Cerebellar_nuclei_total_to_Cerebellum, na.rm = TRUE)

if (is.finite(vol_ratio_max) && vol_ratio_max > 0.5) {
  stop(
    "Matano cerebellar nucleus volumes look too large relative to Stephan Cerebellum. ",
    "This suggests a bad species join or a volume-unit mismatch. Inspect match report before writing."
  )
}
if (is.finite(vol_ratio_median)) {
  message(
    "Cerebellar nuclei volumes look compatible with Stephan volume scale; no cm/mm conversion applied. ",
    "Median TCN/Cerebellum = ", signif(vol_ratio_median, 4), "."
  )
}

# ---------- Join and write outputs ----------
# Remove any existing Matano columns so the script can be rerun safely.
columns_to_add <- setdiff(names(matano_for_join), "target_species")
stephan_base <- stephan[, setdiff(names(stephan), columns_to_add), drop = FALSE]

stephan_base$.row_id_before_join <- seq_len(nrow(stephan_base))
stephan_updated <- merge(
  stephan_base,
  matano_for_join,
  by.x = "Species",
  by.y = "target_species",
  all.x = TRUE,
  sort = FALSE
)
stephan_updated <- stephan_updated[order(stephan_updated$.row_id_before_join), , drop = FALSE]
stephan_updated$.row_id_before_join <- NULL
stephan_updated <- stephan_updated[, c(names(stephan_base)[names(stephan_base) != ".row_id_before_join"], columns_to_add), drop = FALSE]

backup_path <- sub("\\.csv$", "_before_Matano1985a.csv", stephan_path)
if (!file.exists(backup_path)) {
  file.copy(stephan_path, backup_path)
}

write.csv(stephan_updated, stephan_path, row.names = FALSE, na = "")

report_dir <- dirname(matano_path)
match_report_path <- file.path(report_dir, "Matano_1985a_match_report.csv")
unmatched_path <- file.path(report_dir, "Matano_1985a_unmatched_rows.csv")
qc_path <- file.path(report_dir, "Matano_1985a_unit_check.csv")

report_cols <- intersect(
  c(
    "Species", "Species_Matano1985a", "matano_accepted_name", "target_species_raw",
    "target_species", "matched_in_Stephan", "Body_weight_current",
    "Body_weight_current_source", matano_source_cols
  ),
  names(matano)
)
write.csv(matano[, report_cols, drop = FALSE], match_report_path, row.names = FALSE, na = "")
write.csv(matano[!matano$matched_in_Stephan, report_cols, drop = FALSE], unmatched_path, row.names = FALSE, na = "")
write.csv(qc, qc_path, row.names = FALSE, na = "")

message("Matched and joined ", nrow(matano_for_join), " Matano 1985a rows into ", stephan_path, ".")
message("Backup written to: ", backup_path)
message("Match report written to: ", match_report_path)
message("Unmatched rows written to: ", unmatched_path)
message("Unit-check report written to: ", qc_path)
