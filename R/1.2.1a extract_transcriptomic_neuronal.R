setwd("~/Library/CloudStorage/Dropbox/COLLABORATIVE/Do expensive brain regions increase less in humans/analyses_metabol_rate_structure")

## Install and Load up packages
library(anndataR)

# https://github.com/linnarsson-lab/adult-human-brain/
# https://datasets.cellxgene.cziscience.com/a71efd3c-765c-466b-8eca-0b29024094d4.h5ad
# read anndata object from h5ad file
adata <- read_h5ad(
  "~/Library/CloudStorage/OneDrive-AllenInstitute/Analysis_region_celltype_human/a71efd3c-765c-466b-8eca-0b29024094d4.h5ad",
  as = "HDF5AnnData"
)

######################################
## One-time extraction from original h5ad
######################################
obs <- as.data.frame(adata$obs)

saveRDS(
  obs,
  "data/linnarsson_adult_human_brain_obs_metadata_neuronal.rds",
  compress = "xz"
)
