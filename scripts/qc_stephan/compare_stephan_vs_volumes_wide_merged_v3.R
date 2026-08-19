# =============================================================================
# SUPERSEDED 2026-08-18 -- do not run. Use scripts/qc_stephan/compare_stephan.R
#
# The full 565-line original is preserved at
#   scripts/archive/compare_stephan_vs_volumes_wide_merged_v3.R
# (and in git history). It is kept because its mode "model" is still the only
# code that compares the two s3 FITTED matrices -- tree restriction, effective N
# per variable, which species drop out of which analysis. The replacement is
# raw-mode only, by choice: raw source files are the only level at which a 2x
# laterality or 1000x unit slip is still visible.
#
# What moved where:
#   the 60-row variable crosswalk -> metadata/qc_stephan/stephan_wide_crosswalk.csv
#   the species synonyms          -> R/species_aliases.R
#   everything else               -> scripts/qc_stephan/compare_stephan.R
# =============================================================================

stop("compare_stephan_vs_volumes_wide_merged_v3.R is superseded.\n",
     "  Run instead : scripts/qc_stephan/compare_stephan.R\n",
     "  Original at : scripts/archive/compare_stephan_vs_volumes_wide_merged_v3.R",
     call. = FALSE)
