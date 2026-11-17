# Scientific checks on actual extracted outputs, including the published sample
# metadata inconsistency and known anatomy exclusions.
library(dplyr)
rd <- function(p) read.csv(p,check.names=FALSE,stringsAsFactors=FALSE)
s <- rd("tables/s1d/tushaus_sample_inventory.csv")
stopifnot(nrow(s)==38L,!anyDuplicated(s$sample_id),length(unique(s$donor_id))==1L,
          sum(s$duplicate_replicate_label)==2L,
          sum(s$source_region=="SN")==3L)
d <- readRDS("data_intermediate/tushaus_2023_protein_sample_long.rds")
r <- readRDS("data_intermediate/tushaus_2023_protein_regional_summary.rds")
stopifnot(nrow(d)==608986L,!anyDuplicated(d[c("measure","protein_id","sample_id")]),
          length(unique(r$anatomy_id[r$primary_included]))==7L,
          all(is.na(r$rcmr_value[!r$primary_included])))
# Recompute one marker from individual tissue samples, rather than comparing
# two exports of the same summary.
z <- d[d$measure=="LFQ" & d$protein_id=="P14136" & d$source_region=="CN",]
q <- r[r$measure=="LFQ" & r$protein_id=="P14136" & r$source_region=="CN",]
stopifnot(nrow(z)==3L,nrow(q)==1L,
          abs(mean(z$published_log_intensity)-q$mean_published_log_intensity)<1e-10,
          abs(sd(z$published_log_intensity)-q$sd_between_tissue_samples)<1e-10)
screen <- rd("tables/s1d/tushaus_proteome_screen.csv")
stopifnot(all(screen$n==7L),all(screen$n_donors==1L))
for (measure in unique(screen$measure)) {
  a<-screen[screen$measure==measure,]
  stopifnot(max(abs(a$q_spearman_BH-p.adjust(a$p_spearman,"BH")),na.rm=TRUE)<1e-10)
}
a <- rd("tables/s1c/johansen_atlas_samples.csv")
ext <- a[a$role=="extension",]
stopifnot(nrow(ext)==10L,!anyDuplicated(ext$anatomy_id),
          all(ext$map_coverage_fraction>=0.95),all(is.finite(ext$bmax_pmol_ml)),
          !any(as.logical(a$primary_included)),
          !any(as.logical(a$template_registration_validated)),
          all(ext$left_centroid_x[ext$anatomy_id!="vermis"]<0),
          all(ext$right_centroid_x[ext$anatomy_id!="vermis"]>0))
audit <- rd("tables/s1c/johansen_atlas_coverage_audit.csv")
stopifnot(nrow(audit)==26L,sum(audit$status=="not_sampled_anatomy_unresolved")==4L)
expanded <- rd("tables/s1c/johansen_atlas_exploratory_expansion.csv")
stopifnot(nrow(expanded)==21L,!anyDuplicated(expanded$anatomy_id),
          !"vermis" %in% expanded$anatomy_id,
          !"medial_thalamic_nuclei" %in% expanded$anatomy_id,
          !any(expanded$combine_for_inference))
af <- rd("tables/s1c/johansen_atlas_sensitivity_fits.csv")
stopifnot(setequal(af$analysis,c("published_table_primary","atlas_only_exploratory")),
          identical(as.integer(af$n),c(12L,9L)))
b <- rd("tables/s2/beip_replication_fits.csv")
stopifnot(nrow(b)==3L,identical(as.integer(b$n),c(7L,6L,6L)))
cat("Nov 14 protein extraction, atlas coverage and BEIP checks passed.\n")
