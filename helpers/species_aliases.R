# =============================================================================
# species_aliases.R -- the ONE species-label alias table for this repo.
#
# WHY THIS FILE EXISTS
#   The same handful of synonyms were historically redeclared in several
#   comparative-volume scripts, in opposite directions, and they disagreed:
#     compare_stephan_vs_volumes_wide_merged_v3.R  taxon_synonyms  (merge -> modern)
#     build_stephan_primates_reference_sheet.R     sp_xwalk        (Stephan -> merge)
#     0_bind_matano_1985a_to_stephan.R             manual_name_map + an inline
#                                                  rename that WROTE BACK to the
#                                                  source csv (now archived)
#   Declaring them once, in one direction, is the whole point of this file.
#
# THE CANONICAL NAMESPACE IS data_raw/species.nwk TIP LABELS.
#   Not "the modern name" -- the tree tip. This repo is downstream: it does not
#   decide taxonomy, it consumes a fixed tip set. Anything that does not reduce
#   to a tip cannot be modelled, so the tip labels ARE the scope. As of this
#   writing, the project tree contains 59 tips. Comparative datasets such as
#   Evo-M1 volumes_wide.csv are standardized to this fixed namespace.
#
#   Consequence worth stating plainly: some tip labels are NOT the currently
#   accepted binomial. Callithrix_pygmaea is a tip; the accepted name for the
#   pygmy marmoset is Cebuella pygmaea. (It is NOT Callicebus -- that genus is
#   the titi monkeys, and the confusion is easy to make because the
#   Callicebus/Plecturocebus moloch case sits two rows below in this same
#   table.) Renaming the tip is an upstream decision about species.nwk, not
#   something to fix by editing a data file. The alias below absorbs the modern
#   name so that a file arriving with Cebuella_pygmaea still lands on the tip.
#
# TWO NAMESPACES, ONE TABLE
#   tree_label   underscore form, the tip label            e.g. Plecturocebus_moloch
#   merge_label  space form, what Evo-M1-Trait-Data uses   e.g. "Callicebus moloch"
#   Both directions are derived from the single table, so they cannot drift.
#
# Sourced by scripts/02_traits_neocortex_grey_white.R and available to other
# comparative-volume analyses.
# =============================================================================

## --- the table --------------------------------------------------------------
## kind:
##   spelling   same taxon, different spelling in the merge
##   genus      merge keeps a superseded genus
##   pooled     merge pools to "<Genus> sp."; unambiguous, one tip only
##   ambiguous  merge pools to "<Genus> sp." but SEVERAL tips fall under it --
##              never auto-resolved; see pongo_pairs below
species_alias_table <- data.frame(
  tree_label  = c("Callithrix_pygmaea", "Lagothrix_lagotricha", "Plecturocebus_moloch",
                  "Gorilla_gorilla", "Tarsius_syrichta", "Pongo_abelii", "Pongo_pygmaeus"),
  merge_label = c("Callithrix pygmaea", "Lagothrix lagothricha", "Callicebus moloch",
                  "Gorilla sp.", "Tarsius sp.", "Pongo sp.", "Pongo sp."),
  kind        = c("spelling", "spelling", "genus",
                  "pooled", "pooled", "ambiguous", "ambiguous"),
  note        = c("tip and merge both keep Callithrix; accepted name is Cebuella pygmaea",
                  "merge spells the epithet lagothricha",
                  "merge keeps the old genus Callicebus",
                  "merge pools gorillas as 'Gorilla sp.'",
                  "merge pools tarsiers as 'Tarsius sp.'",
                  "merge 'Pongo sp.' is a sensu lato mean (Bornean + Sumatran + hybrid)",
                  "merge 'Pongo sp.' is a sensu lato mean (Bornean + Sumatran + hybrid)"),
  stringsAsFactors = FALSE
)

## --- extra inbound labels ---------------------------------------------------
## Labels that are neither a tip nor a merge label but have turned up in a file
## and must still land on a tip. Left side is matched in SPACE form.
## Cebuella pygmaea is here because an archived data-bind script once rewrote
## that name in a comparative source file, silently knocking the row off the
## project tree and out of its own join. Cheap insurance against a repeat.
species_extra_aliases <- c(
  "Cebuella pygmaea"     = "Callithrix_pygmaea",
  "Gorilla sp"           = "Gorilla_gorilla",
  "Tarsius sp"           = "Tarsius_syrichta",
  "Lagothrix lagothricha" = "Lagothrix_lagotricha",
  "Callicebus moloch"    = "Plecturocebus_moloch"
)

## --- the explicit one-to-many pairing ---------------------------------------
## "Pongo sp." cannot be resolved to a single tip, so it is never guessed.
## Compare it against BOTH tips, flagged, and let the reader decide.
pongo_pairs <- data.frame(
  merge_label = c("Pongo sp.", "Pongo sp."),
  tree_label  = c("Pongo_abelii", "Pongo_pygmaeus"),
  note        = c("merge row is a sensu lato mean; paired tree tip is Sumatran",
                  "merge row is a sensu lato mean; paired tree tip is Bornean"),
  stringsAsFactors = FALSE
)

## --- helpers ----------------------------------------------------------------
.sp_space <- function(x) gsub("\\s+", " ", trimws(gsub("_", " ", as.character(x))))
.sp_under <- function(x) gsub(" ", "_", .sp_space(x))

#' Any species label -> tree-tip label (underscore form).
#'
#' Unknown labels pass through unchanged rather than being guessed at, so an
#' unmapped name surfaces as an unpaired species in the QC report instead of
#' quietly attaching itself to the wrong row. Ambiguous merge labels
#' ("Pongo sp.") are deliberately left alone -- use pongo_pairs.
canon_species <- function(x) {
  sp  <- .sp_space(x)
  out <- .sp_under(sp)
  unamb <- species_alias_table[species_alias_table$kind != "ambiguous", ]
  hit <- match(sp, unamb$merge_label)
  out[!is.na(hit)] <- unamb$tree_label[hit[!is.na(hit)]]
  hit2 <- match(sp, names(species_extra_aliases))
  out[!is.na(hit2)] <- unname(species_extra_aliases[hit2[!is.na(hit2)]])
  out
}

#' Tree-tip label -> the label Evo-M1-Trait-Data uses (space form).
#' Used to join against volumes_unfiltered_select.csv, which is keyed on the
#' merge's own names.
to_merge_species <- function(x) {
  tl  <- .sp_under(x)
  out <- .sp_space(tl)
  hit <- match(tl, species_alias_table$tree_label)
  out[!is.na(hit)] <- species_alias_table$merge_label[hit[!is.na(hit)]]
  out
}
