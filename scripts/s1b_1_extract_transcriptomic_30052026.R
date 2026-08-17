setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))  # repo root (portable; replaces hardcoded path -- see R/project_root.R)

############################################################
# s1b_1 — One-time extraction of cell metadata from the
# Linnarsson adult human brain h5ad files (neuronal +
# nonneuronal). Merges the former s1b_1_n / s1b_1_nn pair,
# which differed only in dataset URL and output filename.
#
# Inputs : the two h5ad files below (local OneDrive copies of
#          https://github.com/linnarsson-lab/adult-human-brain/)
# Outputs: data_intermediate/linnarsson_adult_human_brain_obs_metadata_neuronal.rds
#          data_intermediate/linnarsson_adult_human_brain_obs_metadata_nonneuronal.rds
# Consumers: s1b_2_mapping_*, s1b_2_check_dissection_roi.R
############################################################

## Install and Load up packages
library(anndataR)

DATASETS <- list(
  neuronal = list(
    # https://datasets.cellxgene.cziscience.com/a71efd3c-765c-466b-8eca-0b29024094d4.h5ad
    h5ad = "~/Library/CloudStorage/OneDrive-AllenInstitute/Analysis_region_celltype_human/a71efd3c-765c-466b-8eca-0b29024094d4.h5ad",
    out  = "data_intermediate/linnarsson_adult_human_brain_obs_metadata_neuronal.rds"
  ),
  nonneuronal = list(
    # https://datasets.cellxgene.cziscience.com/c1d05de1-d442-48b1-a32c-86f4f0dc5f82.h5ad
    h5ad = "~/Library/CloudStorage/OneDrive-AllenInstitute/Analysis_region_celltype_human/c1d05de1-d442-48b1-a32c-86f4f0dc5f82.h5ad",
    out  = "data_intermediate/linnarsson_adult_human_brain_obs_metadata_nonneuronal.rds"
  )
)

for (nm in names(DATASETS)) {
  d <- DATASETS[[nm]]
  message("Extracting obs metadata: ", nm)
  adata <- read_h5ad(d$h5ad, as = "HDF5AnnData")
  obs <- as.data.frame(adata$obs)
  saveRDS(obs, d$out, compress = "xz")
  message("  wrote ", d$out, " (", nrow(obs), " cells)")
}
