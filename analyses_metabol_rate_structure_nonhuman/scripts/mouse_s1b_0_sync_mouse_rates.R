############################################################
# mouse_s1b_0 -- Snapshot the mouse regional glucose-utilisation rates
# (CMRgl, umol/100 g/min) from Evo-M1-Trait-Data into this project.
#
# Upstream: Evo-M1-Trait-Data/__merging_cerebral_metabolic_rate/
#           cerebral_metabolic_rate_long.csv  (Species == "Mus musculus")
# Set EVO_M1_TRAIT_DATA to the upstream clone; if it is unavailable the
# existing local snapshot is used. Output:
#   data_raw/Evo_M1_mouse_CMRgl_by_region.csv
#   data_intermediate/evo_m1_mouse_rates_provenance.csv
############################################################
.a <- commandArgs(FALSE); .f <- sub("^--file=", "", .a[grepl("^--file=", .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else "analyses_metabol_rate_structure_nonhuman/scripts", "common.R"))

SRC_REL <- "__merging_cerebral_metabolic_rate/cerebral_metabolic_rate_long.csv"
OUT  <- nh("data_raw/Evo_M1_mouse_CMRgl_by_region.csv")
PROV <- nh("data_intermediate/evo_m1_mouse_rates_provenance.csv")

candidates <- c(Sys.getenv("EVO_M1_TRAIT_DATA", unset = NA),
                "/home/user/alealisousa/evo-m1-trait-data",
                "../Evo-M1-Trait-Data", "../../Evo-M1-Trait-Data")
root <- Filter(function(p) !is.na(p) && file.exists(file.path(p, SRC_REL)), candidates)[1]

if (is.na(root)) {
  if (!file.exists(OUT)) stop("Evo-M1-Trait-Data not found and no local snapshot at ", OUT, call. = FALSE)
  message("Upstream unavailable; using existing snapshot ", OUT)
} else {
  long <- read.csv(file.path(root, SRC_REL), stringsAsFactors = FALSE, check.names = FALSE)
  mouse <- subset(long, Species == "Mus musculus" & Measure == "CMRgl")
  stopifnot(nrow(mouse) > 0, all(mouse$Units == "umol/100g/min"), !anyDuplicated(mouse$Region))
  keep <- c("Species", "Region", "Measure", "Units", "Value", "n_studies",
            "Primary_reference_ids", "Compilations")
  mouse <- mouse[order(mouse$Region), keep]
  write.csv(mouse, OUT, row.names = FALSE)
  head_sha <- tryCatch(system2("git", c("-C", root, "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
                       error = function(e) NA_character_)
  write.csv(data.frame(source = file.path(root, SRC_REL), upstream_commit = head_sha[1],
                       n_regions = nrow(mouse), synced = format(Sys.time(), "%Y-%m-%d")),
            PROV, row.names = FALSE)
  message("Wrote ", nrow(mouse), " mouse CMRgl regions to ", OUT)
}
