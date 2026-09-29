# =====================================================================
# 00_sync_evo_m1_inputs-restricted.R
#
# OPT-IN provenance layer. Pulls the RESTRICTED raw source files (author /
# ERA-team workbooks) that the public group-mean inputs were derived from,
# into a GITIGNORED data_restricted/ area, so that
# 00_derive_from_restricted.R can regenerate and verify the public derived
# files for peer review. This is NOT part of the default public pipeline and
# only runs when Evo-M1-Trait-Data-restricted is available.
#
# EVO_M1_TRAIT_DATA_RESTRICTED is the single configurable upstream root.
# Every source/destination pair is declared in:
#   metadata/evo_m1_input_manifest-restricted.csv
# =====================================================================

# Locate the repository from either Rscript or an interactive/source call.
.script_args <- commandArgs(trailingOnly = FALSE)
.file_arg <- .script_args[grepl("^--file=", .script_args)]
.start_dir <- if (length(.file_arg)) {
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

manifest_path <- "metadata/evo_m1_input_manifest-restricted.csv"
provenance_path <- "data_intermediate/evo_m1_input_provenance-restricted.csv"

if (!file.exists(manifest_path)) {
  stop("Restricted input manifest not found: ", manifest_path, call. = FALSE)
}

manifest <- read.csv(
  manifest_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
required_columns <- c(
  "dependency_id", "study", "source_relative_path", "local_path"
)
missing_columns <- setdiff(required_columns, names(manifest))
if (length(missing_columns)) {
  stop(
    manifest_path, " is missing required column(s): ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

for (column in required_columns) {
  manifest[[column]] <- trimws(as.character(manifest[[column]]))
  if (anyNA(manifest[[column]]) || any(manifest[[column]] == "")) {
    stop(manifest_path, "$", column, " contains missing values.",
         call. = FALSE)
  }
}
if (anyDuplicated(manifest$dependency_id)) {
  stop(manifest_path, " contains duplicate dependency_id values.",
       call. = FALSE)
}
if (anyDuplicated(manifest$local_path)) {
  stop(manifest_path, " maps more than one source to the same local path.",
       call. = FALSE)
}

# Restricted raw files must land in the gitignored data_restricted/ area so
# they are never committed to this public analysis repository.
unsafe_source <- grepl("^(/|[A-Za-z]:)", manifest$source_relative_path) |
  grepl("(^|/)\\.\\.(/|$)", manifest$source_relative_path)
unsafe_local <- !startsWith(manifest$local_path, "data_restricted/") |
  grepl("(^|/)\\.\\.(/|$)", manifest$local_path)
if (any(unsafe_source) || any(unsafe_local)) {
  stop(
    "Manifest paths must be safe relative source paths and data_restricted/ targets.",
    call. = FALSE
  )
}

configured_root <- trimws(Sys.getenv("EVO_M1_TRAIT_DATA_RESTRICTED", unset = ""))
standard_root <- path.expand(paste0(
  "~/Library/CloudStorage/OneDrive-AllenInstitute/",
  "Evo-M1-Trait-Data-restricted"
))

if (nzchar(configured_root)) {
  if (!dir.exists(configured_root)) {
    stop(
      "EVO_M1_TRAIT_DATA_RESTRICTED does not identify an existing directory: ",
      configured_root,
      call. = FALSE
    )
  }
  evo_m1_root <- normalizePath(configured_root)
  root_source <- "EVO_M1_TRAIT_DATA_RESTRICTED"
} else if (dir.exists(standard_root)) {
  evo_m1_root <- normalizePath(standard_root)
  root_source <- "standard macOS location"
} else {
  evo_m1_root <- NA_character_
  root_source <- NA_character_
}

md5_file <- function(path) {
  unname(as.character(tools::md5sum(path)))
}

atomic_copy <- function(source, destination) {
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(
    pattern = paste0(".", basename(destination), ".sync-"),
    tmpdir = dirname(destination)
  )
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)

  if (!file.copy(
    source, temporary,
    overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE
  )) {
    stop("Could not stage restricted input: ", source, call. = FALSE)
  }
  if (!identical(md5_file(source), md5_file(temporary))) {
    stop("Checksum mismatch while staging: ", source, call. = FALSE)
  }

  backup <- NA_character_
  if (file.exists(destination)) {
    backup <- tempfile(
      pattern = paste0(".", basename(destination), ".backup-"),
      tmpdir = dirname(destination)
    )
    if (!file.rename(destination, backup)) {
      stop("Could not stage existing local input: ", destination,
           call. = FALSE)
    }
  }

  if (!file.rename(temporary, destination)) {
    if (!is.na(backup) && file.exists(backup)) {
      file.rename(backup, destination)
    }
    stop("Could not install refreshed restricted input: ", destination,
         call. = FALSE)
  }
  if (!is.na(backup) && file.exists(backup)) unlink(backup)
  invisible(TRUE)
}

if (is.na(evo_m1_root)) {
  # This is expected on any machine without the private repo. The restricted
  # provenance layer is optional; the default public pipeline does not need it.
  stop(
    "Evo-M1-Trait-Data-restricted is not available, so restricted raw sources ",
    "cannot be synced. This is expected outside peer-review reproduction. ",
    "Set EVO_M1_TRAIT_DATA_RESTRICTED to the private repository root to run ",
    "the provenance layer.",
    call. = FALSE
  )
}

source_paths <- file.path(evo_m1_root, manifest$source_relative_path)
missing_upstream <- manifest$source_relative_path[!file.exists(source_paths)]
if (length(missing_upstream)) {
  stop(
    "No files were changed. Restricted source file(s) are missing: ",
    paste(missing_upstream, collapse = ", "),
    call. = FALSE
  )
}

source_hashes <- vapply(source_paths, md5_file, character(1))
local_hashes <- vapply(manifest$local_path, function(path) {
  if (file.exists(path)) md5_file(path) else NA_character_
}, character(1))
changed <- is.na(local_hashes) | source_hashes != local_hashes

for (i in which(changed)) {
  atomic_copy(source_paths[i], manifest$local_path[i])
}

hashes <- vapply(manifest$local_path, md5_file, character(1))
if (!identical(unname(source_hashes), unname(hashes))) {
  stop("One or more refreshed restricted inputs failed checksum validation.",
       call. = FALSE)
}
message(
  "Restricted root (", root_source, "): ", evo_m1_root,
  "\nRefreshed ", sum(changed), " restricted raw file(s); ", sum(!changed),
  " already current (in gitignored data_restricted/)."
)

provenance <- data.frame(
  dependency_id = manifest$dependency_id,
  study = manifest$study,
  source_relative_path = manifest$source_relative_path,
  local_path = manifest$local_path,
  md5 = hashes,
  stringsAsFactors = FALSE
)
dir.create(dirname(provenance_path), recursive = TRUE, showWarnings = FALSE)
write.csv(provenance, provenance_path, row.names = FALSE, na = "")
message("Wrote restricted dependency provenance: ", provenance_path)
