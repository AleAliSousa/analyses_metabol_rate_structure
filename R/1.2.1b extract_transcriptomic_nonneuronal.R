setwd("~/Library/CloudStorage/Dropbox/COLLABORATIVE/Do expensive brain regions increase less in humans/analyses_metabol_rate_structure")

## Install and Load up packages
library(anndataR)

# https://github.com/linnarsson-lab/adult-human-brain/
# https://datasets.cellxgene.cziscience.com/c1d05de1-d442-48b1-a32c-86f4f0dc5f82.h5ad
# read anndata object from h5ad file
adata <- read_h5ad(
  "~/Library/CloudStorage/OneDrive-AllenInstitute/Analysis_region_celltype_human/c1d05de1-d442-48b1-a32c-86f4f0dc5f82.h5ad",
  as = "HDF5AnnData"
)

######################################
## One-time extraction from original h5ad
######################################
obs <- as.data.frame(adata$obs)

saveRDS(
  obs,
  "data/linnarsson_adult_human_brain_obs_metadata_nonneuronal.rds",
  compress = "xz"
)
