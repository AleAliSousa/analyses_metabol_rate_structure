# ============================================================
# Compare the two saved "checks" datasets side by side  -- DIAGNOSTIC VERSION
# ============================================================
# The two s3 scripts read different source files:
#   s3_..._UPDATED_MU_PLOTS_1_PATCHED.R  <- data_raw/Stephan_primates.csv
#   s3_..._VOLUMES_WIDE_SELECT.R         <- data_raw/volumes_wide_select.csv
# and each writes the exact matrix it modelled to checks/.
#
# The two matrices LOOK similar cell-by-cell, but the fits differ because of
# three things a plain side-by-side view hides:
#   (1) Preferred_brain_volume is DEFINED differently in the two scripts
#       (coalesce(Brain_volume, Brainvol, Total_brain_net_volume) vs
#        Total_brain_net_volume_Vol.mm3 alone).  It is the x-axis of every
#        PGLS, so small offsets propagate into every slope.
#   (2) Species labels in volumes_wide_select.csv do not all match the tree
#       tips, so those rows are silently pruned before fitting.
#   (3) A handful of region cells genuinely disagree (some by >25%), and they
#       are concentrated in the high-leverage apes.
#
# This script surfaces all three.
#
# NOTE ON OUT-OF-SCOPE ROWS.  Stephan_primates.csv is primates-only (59 spp,
# exactly the tree tips).  volumes_wide_select.csv draws on the wider
# Stephan-Frahm-Baron mammal series, so it carries ~51 extra rows that are not
# tips: some are non-primates (Chrysochloris, Sorex, Tenrec, Talpa, Tupaia...),
# some are primates whose LABEL does not match a tip (abbreviations like
# "C._medius", "sp."/"spp." placeholders, genus synonymy, subspecies trinomials).
#
# This repo is DOWNSTREAM and does not own taxonomy, so nothing here asserts
# which of those rows is a primate.  The analysis phylogeny is the scope
# definition: a row is in scope iff its label is a tip.  Off-tree rows are
# reported below as a single list to be resolved UPSTREAM, in the compilation
# project's species key.
#
# Inputs   : checks/Stephan_data_used.csv
#            checks/volumes_wide_select_data_used.csv
#            data_raw/species.nwk                (optional but recommended)
# Outputs  : checks/Stephan_vs_volumes_wide_side_by_side.csv   (wide view + diffs)
#            checks/Stephan_vs_volumes_wide_cell_discrepancies.csv (long, only problems)
#            checks/Stephan_vs_volumes_wide_variable_summary.csv   (per-variable N)
#            checks/Stephan_vs_volumes_wide_species_coverage.csv   (who drops out)
# ============================================================

stephan_file <- "checks/Stephan_data_used.csv"
wide_file    <- "checks/volumes_wide_select_data_used.csv"
tree_file    <- "data_raw/species.nwk"

out_wide     <- "checks/Stephan_vs_volumes_wide_side_by_side.csv"
out_cells    <- "checks/Stephan_vs_volumes_wide_cell_discrepancies.csv"
out_varsum   <- "checks/Stephan_vs_volumes_wide_variable_summary.csv"
out_species  <- "checks/Stephan_vs_volumes_wide_species_coverage.csv"

# Relative tolerance below which two values count as "same" (rounding only).
REL_TOL <- 1e-6
# Relative difference above which a cell is called a MATERIAL disagreement.
REL_FLAG <- 0.01
# Restrict the main side-by-side view to species that can actually reach a model
# (i.e. tree tips).  Set FALSE to dump every row from both files.
RESTRICT_TO_TREE <- TRUE

# ---------- Read saved check files ----------
stephan <- read.csv(stephan_file, check.names = FALSE, stringsAsFactors = FALSE)
wide    <- read.csv(wide_file,    check.names = FALSE, stringsAsFactors = FALSE)

# ---------- Tree tips: the species that actually survive into the models ----------
tree_tips <- NULL
if (file.exists(tree_file)) {
  if (requireNamespace("ape", quietly = TRUE)) {
    tree_tips <- ape::read.tree(tree_file)$tip.label
  } else {
    warning("Package 'ape' not available; species-coverage section will be skipped.")
  }
} else {
  warning("Tree file not found at ", tree_file,
          "; species-coverage section will be skipped.")
}

# ---------- Exact variable equivalences used by the two analysis scripts ----------
# Rows are in the same order as region_labels in the scripts.
variable_map <- data.frame(
  label = c(
    "Preferred brain volume",
    "Corpus geniculatum laterale",
    "Amygdala",
    "Pallidum",
    "Neocortex white",
    "Insular cortex (grey)",
    "Nucleus subthalamicus Luysi",
    "Capsula interna",
    "Striatum",
    "Area striata grey",
    "Neocortex grey",
    "Mesencephalon",
    "Cerebellum",
    "Hippocampus",
    "Nucleus dentatus cerebelli"
  ),
  stephan = c(
    "Preferred_brain_volume",
    "LGN_Sousa",
    "Amygdala",
    "Pallidum",
    "NeoW_Frahm",
    "Total_insula_volume_L",
    "Nucleus_subthalamicus",
    "Capsula_interna",
    "Striatum",
    "ASG_Sousa",
    "NeoG_Frahm",
    "Mesencephalon",
    "Cerebellum",
    "Hippocampus",
    "Lateral_cerebellar_nuclei"
  ),
  volumes_wide = c(
    "Preferred_brain_volume",
    "Corpus_geniculatum_laterale_Vol.mm3",
    "Amygdala_Vol.mm3",
    "Pallidum_Vol.mm3",
    "Neocortex_white_matter_Vol.mm3",
    "Insula_left_Vol.mm3",
    "Nucleus_subthalamicus_Vol.mm3",
    "Capsula_interna_Vol.mm3",
    "Striatum_Vol.mm3",
    "Area_striata_grey_matter_Vol.mm3",
    "Neocortex_grey_matter_Vol.mm3",
    "Mesencephalon_Vol.mm3",
    "Cerebellum_Vol.mm3",
    "Hippocampus_Vol.mm3",
    "Lateral_cerebellar_nuclei_Vol.mm3"
  ),
  stringsAsFactors = FALSE
)

# ---------- Structural checks ----------
required_stephan <- c("Species", variable_map$stephan)
required_wide    <- c("Species", variable_map$volumes_wide)

missing_stephan <- setdiff(required_stephan, names(stephan))
missing_wide    <- setdiff(required_wide, names(wide))

if (length(missing_stephan) > 0) {
  stop("Missing column(s) from Stephan check file: ",
       paste(missing_stephan, collapse = ", "))
}
if (length(missing_wide) > 0) {
  stop("Missing column(s) from volumes_wide check file: ",
       paste(missing_wide, collapse = ", "))
}

if (anyDuplicated(stephan$Species)) {
  warning("Duplicate Species rows found in Stephan check file: ",
          paste(unique(stephan$Species[duplicated(stephan$Species)]), collapse = ", "))
}
if (anyDuplicated(wide$Species)) {
  warning("Duplicate Species rows found in volumes_wide check file: ",
          paste(unique(wide$Species[duplicated(wide$Species)]), collapse = ", "))
}

# ---------- Match species, retaining species from either dataset ----------
all_species_raw <- unique(c(stephan$Species, wide$Species))

# Optionally drop everything that cannot reach a model. This is what removes the
# 31 non-primate rows (Chrysochloris, Sorex, Tenrec, ...) from the main view.
if (RESTRICT_TO_TREE && !is.null(tree_tips)) {
  all_species <- all_species_raw[all_species_raw %in% tree_tips]
} else {
  all_species <- all_species_raw
}

comparison  <- data.frame(Species = all_species, stringsAsFactors = FALSE)

idx_s <- match(all_species, stephan$Species)
idx_w <- match(all_species, wide$Species)

# Does this species reach the models at all?
if (!is.null(tree_tips)) {
  comparison$on_tree <- ifelse(all_species %in% tree_tips, "yes", "PRUNED")
}
comparison$in_Stephan      <- ifelse(all_species %in% stephan$Species, "yes", "-")
comparison$in_volumes_wide <- ifelse(all_species %in% wide$Species,    "yes", "-")

# ---------- Cell-level classifier ----------
classify_cell <- function(a, b) {
  both_na <- is.na(a) & is.na(b)
  only_a  <- !is.na(a) &  is.na(b)
  only_b  <-  is.na(a) & !is.na(b)
  both    <- !is.na(a) & !is.na(b)

  rel <- rep(NA_real_, length(a))
  denom <- pmax(abs(a), abs(b))
  rel[both] <- ifelse(denom[both] == 0, 0,
                      abs(a[both] - b[both]) / denom[both])

  status <- rep(NA_character_, length(a))
  status[both_na] <- "both_missing"
  status[only_a]  <- "Stephan_only"
  status[only_b]  <- "volumes_wide_only"
  status[both & rel <= REL_TOL]                    <- "same"
  status[both & rel >  REL_TOL & rel <= REL_FLAG]  <- "rounding"
  status[both & rel >  REL_FLAG]                   <- "DIFFERENT"

  list(status = status, rel = rel)
}

# ---------- Add each equivalent pair, plus diff + % + status ----------
cell_rows <- list()

for (i in seq_len(nrow(variable_map))) {
  base_name <- gsub("[^A-Za-z0-9]+", "_", variable_map$label[i])
  base_name <- gsub("^_|_$", "", base_name)

  a <- suppressWarnings(as.numeric(stephan[[variable_map$stephan[i]]][idx_s]))
  b <- suppressWarnings(as.numeric(wide[[variable_map$volumes_wide[i]]][idx_w]))

  cc <- classify_cell(a, b)

  comparison[[paste0(base_name, "_Stephan")]]      <- a
  comparison[[paste0(base_name, "_volumes_wide")]] <- b
  comparison[[paste0(base_name, "_pct_diff")]]     <-
    round(100 * (b - a) / a, 2)
  comparison[[paste0(base_name, "_status")]]       <- cc$status

  keep <- cc$status %in% c("Stephan_only", "volumes_wide_only", "DIFFERENT")
  if (any(keep)) {
    cell_rows[[length(cell_rows) + 1]] <- data.frame(
      Species      = all_species[keep],
      variable     = variable_map$label[i],
      col_Stephan  = variable_map$stephan[i],
      col_wide     = variable_map$volumes_wide[i],
      Stephan      = a[keep],
      volumes_wide = b[keep],
      pct_diff     = round(100 * (b[keep] - a[keep]) / a[keep], 2),
      status       = cc$status[keep],
      on_tree      = if (is.null(tree_tips)) NA_character_
                     else ifelse(all_species[keep] %in% tree_tips, "yes", "PRUNED"),
      stringsAsFactors = FALSE
    )
  }
}

cell_discrepancies <- if (length(cell_rows)) do.call(rbind, cell_rows) else
  data.frame(Species = character(0))

# Problems on pruned species cannot affect any fit; sort so the ones that
# actually matter are at the top.
if (nrow(cell_discrepancies) > 0 && !is.null(tree_tips)) {
  cell_discrepancies <- cell_discrepancies[
    order(cell_discrepancies$on_tree != "yes",
          cell_discrepancies$status  != "DIFFERENT",
          -abs(ifelse(is.na(cell_discrepancies$pct_diff), 0,
                      cell_discrepancies$pct_diff))), ]
}

# ---------- Per-variable summary: effective N in each script ----------
# N is counted ON TREE TIPS ONLY, because that is the set each PGLS sees.
in_model <- if (is.null(tree_tips)) rep(TRUE, length(all_species)) else
  all_species %in% tree_tips

var_summary <- do.call(rbind, lapply(seq_len(nrow(variable_map)), function(i) {
  a <- suppressWarnings(as.numeric(stephan[[variable_map$stephan[i]]][idx_s]))
  b <- suppressWarnings(as.numeric(wide[[variable_map$volumes_wide[i]]][idx_w]))
  cc <- classify_cell(a, b)
  data.frame(
    variable        = variable_map$label[i],
    N_Stephan       = sum(!is.na(a) & in_model),
    N_volumes_wide  = sum(!is.na(b) & in_model),
    N_shared        = sum(!is.na(a) & !is.na(b) & in_model),
    n_Stephan_only  = sum(cc$status == "Stephan_only"      & in_model, na.rm = TRUE),
    n_wide_only     = sum(cc$status == "volumes_wide_only" & in_model, na.rm = TRUE),
    n_DIFFERENT     = sum(cc$status == "DIFFERENT"         & in_model, na.rm = TRUE),
    max_abs_pct_diff = {
      d <- abs(100 * (b - a) / a)[in_model]
      if (all(is.na(d))) NA_real_ else round(max(d, na.rm = TRUE), 2)
    },
    stringsAsFactors = FALSE
  )
}))

# ---------- Species coverage: who is lost, and from where ----------
# Built on the UNRESTRICTED species list, so the pruned rows are still audited.
species_coverage <- data.frame(
  Species         = all_species_raw,
  in_Stephan      = all_species_raw %in% stephan$Species,
  in_volumes_wide = all_species_raw %in% wide$Species,
  on_tree         = if (is.null(tree_tips)) NA else all_species_raw %in% tree_tips,
  stringsAsFactors = FALSE
)

# Verdicts, in descending order of how much they should worry you.
# No clade is asserted: off-tree rows are handed upstream for name resolution.
species_coverage$verdict <- with(species_coverage, ifelse(
  !is.na(on_tree) & !on_tree,
    "off-tree - out of scope, resolve upstream",
  ifelse(in_Stephan & !in_volumes_wide,
    "on tree but absent from volumes_wide - LOST",
  ifelse(!in_Stephan & in_volumes_wide,
    "extra in volumes_wide",
    "in both"))))

sev <- match(species_coverage$verdict,
             c("on tree but absent from volumes_wide - LOST",
               "extra in volumes_wide",
               "in both",
               "off-tree - out of scope, resolve upstream"))
species_coverage <- species_coverage[order(sev, species_coverage$Species), ]

# ---------- Save ----------
dir.create("checks", showWarnings = FALSE, recursive = TRUE)
write.csv(comparison,         out_wide,    row.names = FALSE, na = "")
write.csv(cell_discrepancies, out_cells,   row.names = FALSE, na = "")
write.csv(var_summary,        out_varsum,  row.names = FALSE, na = "")
write.csv(species_coverage,   out_species, row.names = FALSE, na = "")

# ---------- Console report ----------
cat("\n============================================================\n")
cat("Stephan  vs  volumes_wide_select : why the two s3 scripts differ\n")
cat("============================================================\n\n")

cat("Rows read      : Stephan", nrow(stephan), "| volumes_wide", nrow(wide), "\n")
if (!is.null(tree_tips)) {
  cat("Tree tips      :", length(tree_tips), "(primates only)\n")
  cat("Modelled       : Stephan", sum(stephan$Species %in% tree_tips),
      "| volumes_wide", sum(wide$Species %in% tree_tips), "\n")
  if (RESTRICT_TO_TREE)
    cat("Side-by-side view restricted to tree tips (RESTRICT_TO_TREE = TRUE):",
        nrow(comparison), "rows;",
        length(all_species_raw) - length(all_species), "off-tree rows moved to",
        basename(out_species), "\n")
}
cat("\n")

lost <- species_coverage$Species[
  species_coverage$verdict == "on tree but absent from volumes_wide - LOST"]
if (length(lost)) {
  cat("--- Tree species present in Stephan but ABSENT from volumes_wide ---\n")
  cat("    (dropped from every volumes_wide fit)\n")
  for (sp in lost) cat("      ", sp, "\n")
  cat("\n")
}

offtree <- species_coverage$Species[
  species_coverage$verdict == "off-tree - out of scope, resolve upstream"]
if (length(offtree)) {
  cat("--- Off-tree rows in volumes_wide:", length(offtree), "---\n")
  cat("    Not tips of the analysis phylogeny, so never modelled. A mix of\n")
  cat("    genuinely out-of-clade taxa and labels that simply do not match a\n")
  cat("    tip (abbreviations, sp./spp. placeholders, synonymy, trinomials).\n")
  cat("    Deciding which is which is an UPSTREAM job for the compilation\n")
  cat("    project's species key -- not something this repo should encode.\n")
  for (sp in offtree) cat("      ", sp, "\n")
  cat("\n")
}

cat("--- Per-variable effective N (tree tips only) ---\n")
print(var_summary, row.names = FALSE)
cat("\n")

if (nrow(cell_discrepancies) > 0) {
  onmodel <- cell_discrepancies[
    is.na(cell_discrepancies$on_tree) | cell_discrepancies$on_tree == "yes", ]
  cat("--- Cell-level discrepancies affecting the models:",
      nrow(onmodel), "---\n")
  cat("    (full list in ", out_cells, ")\n", sep = "")
  print(utils::head(onmodel[onmodel$status == "DIFFERENT", ], 20), row.names = FALSE)
  cat("\n")
}

cat("Wrote:\n  ", out_wide, "\n  ", out_cells, "\n  ",
    out_varsum, "\n  ", out_species, "\n", sep = "")

# Optional: view in RStudio
# View(comparison); View(cell_discrepancies); View(var_summary)
