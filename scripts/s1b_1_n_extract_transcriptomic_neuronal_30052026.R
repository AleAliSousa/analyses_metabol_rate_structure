setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))  # repo root (portable; replaces hardcoded path -- see R/project_root.R)

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
  "data_intermediate/linnarsson_adult_human_brain_obs_metadata_neuronal.rds",
  compress = "xz"
)
