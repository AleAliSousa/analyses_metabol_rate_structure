# R entry point: sample images with affine-aware Python NIfTI libraries, then
# produce regional sensitivity analyses using the common R statistics helper.
library(dplyr)
library(ggplot2)
source("helpers/regional_robustness.R")
local_python <- file.path(getwd(), ".venv-atlas", "bin", "python")
python <- Sys.getenv("JOHANSEN_PYTHON", if (file.exists(local_python)) local_python else Sys.which("python3"))
if (!nzchar(python)) stop("Python required; install helpers/requirements-atlas.txt")
status <- system2(python, shQuote("helpers/sample_johansen_atlas.py"))
if (status != 0L) stop("Atlas sampler failed. Install helpers/requirements-atlas.txt in JOHANSEN_PYTHON environment.")
a <- read.csv("tables/s1c/johansen_atlas_samples.csv",stringsAsFactors=FALSE)
a <- a %>% mutate(across(c(small_nucleus,named_anatomy_match,primary_included,
                            template_registration_validated),as.logical))
stopifnot(!anyNA(a$small_nucleus),!anyNA(a$named_anatomy_match))
base <- read.csv("tables/s1c/synaptic_density_rcmrglc_matched_regions.csv") %>%
  transmute(anatomy_id,rcmr_value=rcmrglc_mean,bmax_pmol_ml=synaptic_density_mean,
            source="Published Table 2",small_nucleus=anatomy_id=="nucleus_accumbens")
new <- a %>% filter(role=="extension",named_anatomy_match,map_coverage_fraction>=0.95) %>%
  transmute(anatomy_id,rcmr_value=rcMRGlc,bmax_pmol_ml,
            source="Atlas extraction (template review pending)",small_nucleus)
stopifnot(!anyDuplicated(c(base$anatomy_id,new$anatomy_id)))
expanded <- bind_rows(base,new) %>%
  mutate(combine_for_inference=FALSE,
    comparability_note="Do not pool: map estimates are 30-65% below Table 2 in five overlap regions; calibration and template review pending")
write_analysis_csv(expanded,"tables/s1c/johansen_atlas_exploratory_expansion.csv")
sets <- list(published_table_primary=base,
             atlas_only_exploratory=new)
# Do not fit the pooled table/map data: the overlap check fails comparability.
# Excluding all new small nuclei leaves only internal capsule, too few to fit.
fits <- bind_rows(lapply(names(sets),function(nm)cbind(analysis=nm,
  regional_fit(sets[[nm]],"rcmr_value","bmax_pmol_ml"))))
fits$interpretation <- ifelse(fits$analysis=="published_table_primary",
  "Existing table-based reference", "Exploratory atlas-only fit: smoothed map and template correspondence unvalidated; overlap values differ 30-65% from Table 2")
write_analysis_csv(fits,"tables/s1c/johansen_atlas_sensitivity_fits.csv")
overlap <- a %>% filter(role=="table_overlap_check") %>%
  left_join(select(base,anatomy_id,table_bmax=bmax_pmol_ml),by="anatomy_id") %>%
  mutate(difference_pct=100*(bmax_pmol_ml-table_bmax)/table_bmax)
write_analysis_csv(overlap,"tables/s1c/johansen_atlas_table_overlap_check.csv")
mask_comparison <- a %>% filter(role %in% c("extension","alternative_mask")) %>%
  group_by(anatomy_id) %>% summarise(n_mask_definitions=n(),
    min_bmax=min(bmax_pmol_ml),max_bmax=max(bmax_pmol_ml),
    relative_range_pct=100*(max_bmax-min_bmax)/mean(bmax_pmol_ml),
    .groups="drop")
write_analysis_csv(mask_comparison,"tables/s1c/johansen_atlas_mask_sensitivity.csv")
p <- ggplot(expanded,aes(rcmr_value,bmax_pmol_ml,color=source,shape=small_nucleus)) +
  geom_point(size=2.6) + facet_wrap(~source,nrow=1) + ggrepel::geom_text_repel(aes(label=anatomy_id),size=2.6,seed=14,max.overlaps=Inf) +
  theme_classic(base_size=11) + labs(title="Johansen atlas extension: exploratory regional estimates",
    subtitle="Separate sources: overlap checks do not support pooling table and map estimates",
    x="Heiss rCMRGlc (µmol / 100 g / min)",y="SV2A Bmax (pmol/mL)",color=NULL,shape="Small nucleus",
    caption="Map estimates are 30–65% below Table 2 in five overlap regions. MNI template variants require registration review. Spatial SD and mean voxelwise subject SD are not regional standard errors.\nVermis is sampled descriptively but excluded here because its tissue match remains uncertain.") +
  theme(legend.position="bottom")
ggsave("figs/s1c/johansen_atlas_extension.png",p,width=12,height=7,dpi=180)
ggsave("figs/s1c/johansen_atlas_extension.pdf",p,width=12,height=7)
print(fits)
