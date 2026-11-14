# Feasibility audit only: no change to primary Heiss fits or source measurements.
source("helpers/regional_robustness.R")
root <- Sys.getenv("EVO_M1_ROOT", "/Users/crossmodal/Library/CloudStorage/OneDrive-AllenInstitute/Species/Evo-M1-Trait-Data")
path <- file.path(root, "Borghammer_etal_2012", "Borghammer_etal_2012_Table2.csv")
d <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(nrow(d) == 108L)
c <- d[d$group_as_published == "Controls", ]
stopifnot(nrow(c) == 54L, !anyDuplicated(c[c("structure_as_published", "side_as_published", "normalization_as_published")]))
# Equal-weight hemisphere summaries, not volume-weighted means. No bilateral
# dispersion: the printed type and left/right covariance are unavailable.
groups <- split(c, interaction(c$structure_as_published, c$normalization_as_published, drop = TRUE))
b <- do.call(rbind, lapply(groups, function(z) data.frame(
 structure = z$structure_as_published[1], normalization = z$normalization_as_published[1],
 normalized_rate = mean(z$normalized_CMRglc_mean), n_printed_sides = nrow(z),
 hemisphere_summary = if(nrow(z)==2) "equal mean of published L/R means" else "single published row",
 n_controls = 11L)))
stopifnot(nrow(b)==30L)
write_analysis_csv(b, "tables/anatomy/borghammer_controls_summary.csv")
# Native ROI labels must be used before the existing Heiss parent aggregation.
map <- data.frame(
 borghammer_structure = c("Putamen","Caudate","GPe","GPi","VA","CM","CM","VL","STN"),
 siletti_roi = c("Human Pu","Human CaB","Human GPe","Human GPi","Human VA","Human CM","Human CM-Pf","Human VLN","Human STH"),
 proposed_status = c(rep("candidate_named_structure",6),"exclude_mixed_nuclei","boundary_unresolved","sensitivity_nearby_tissue"),
 reason = c("Putamen samples", "Caudate body is a sampled subregion; no identical-mask claim",
            "External pallidum kept separate", "Internal pallidum kept separate",
            "Ventral anterior nucleus", "Pure centromedian sample only",
            "Combined centromedian and parafascicular is broader than CM",
            "Source says ventral group of lateral nucleus; equivalence to Borghammer VL unproven",
            "Source says subthalamic nucleus AND nearby"))
coverage <- list()
for (kind in c("neuronal","nonneuronal")) {
 obs <- readRDS(paste0("data_intermediate/linnarsson_adult_human_brain_obs_metadata_",kind,".rds"))
 a <- unique(obs[obs$roi %in% map$siletti_roi, c("roi","donor_id","dissection")])
 coverage[[kind]] <- do.call(rbind,lapply(seq_len(nrow(map)),function(i) {
  z <- a[a$roi==map$siletti_roi[i], ]
  cbind(map[i,], dataset=kind, n_donors=length(unique(z$donor_id)),
        donors=paste(sort(unique(z$donor_id)),collapse=" | "),
        source_dissection=paste(unique(z$dissection),collapse=" | "))
 }))
}
write_analysis_csv(do.call(rbind,coverage), "tables/anatomy/borghammer_siletti_feasibility.csv")
# Rank stability is descriptive: three normalizations are the same people.
structs <- c("Putamen","Caudate","GPe","GPi","VA","VL","CM","STN")
w <- reshape(b[b$structure %in% structs,c("structure","normalization","normalized_rate")],
             idvar="structure",timevar="normalization",direction="wide")
print(w)
print(cor(w[-1], method="spearman"))
stopifnot(all(coverage$neuronal$n_donors[coverage$neuronal$borghammer_structure %in% c("VA","CM") &
          coverage$neuronal$proposed_status == "candidate_named_structure"] == 1L))
cat("Borghammer feasibility: six named-structure candidates; VA and pure CM each one donor. No primary fits changed.\n")
