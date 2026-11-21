############################################################
# mouse_s1a_0 -- Snapshot mouse isotropic-fractionator cell counts
# (Herculano-Houzel et al. 2011; cerebral cortex, cerebellum, olfactory bulb,
# rest of brain, whole brain) from Evo-M1-Trait-Data. Used only as an
# independent validity check of the MERFISH neuron fractions (n = 2 rate-matched
# regions, too few for a regression).
# Upstream: __merging_cell_counts_isotropic_fractionator/cell_counts_isotropic_fractionator_long.csv
# Output:   data_raw/Evo_M1_mouse_isotropic_fractionator.csv
############################################################
.a <- commandArgs(FALSE); .f <- sub("^--file=", "", .a[grepl("^--file=", .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else "analyses_metabol_rate_structure_nonhuman/scripts", "common.R"))
SRC_REL <- "__merging_cell_counts_isotropic_fractionator/cell_counts_isotropic_fractionator_long.csv"
OUT <- nh("data_raw/Evo_M1_mouse_isotropic_fractionator.csv")
cands <- c(Sys.getenv("EVO_M1_TRAIT_DATA", unset = NA), "/home/user/alealisousa/evo-m1-trait-data",
           "../Evo-M1-Trait-Data", "../../Evo-M1-Trait-Data")
root <- Filter(function(p) !is.na(p) && file.exists(file.path(p, SRC_REL)), cands)[1]
if (is.na(root)) {
  if (!file.exists(OUT)) stop("Evo-M1-Trait-Data not found and no local snapshot: ", OUT, call. = FALSE)
  message("Upstream unavailable; using existing snapshot")
} else {
  x <- read.csv(file.path(root, SRC_REL), stringsAsFactors = FALSE)
  x <- x[x$Species == "Mus musculus", c("Species", "Variable", "Value", "Source")]
  stopifnot(nrow(x) > 0)
  write.csv(x[order(x$Variable), ], OUT, row.names = FALSE)
  message("Wrote ", nrow(x), " mouse isotropic-fractionator values")
}
