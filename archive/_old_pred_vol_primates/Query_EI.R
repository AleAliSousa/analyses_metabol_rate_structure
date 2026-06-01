library(cellxgene.census)
library(dplyr)
library(SingleCellExperiment)

census <- open_soma(census_version = "2025-11-08")  # optional but reproducible

# 1) Read dataset metadata (this is where collection_id lives)
datasets <- as.data.frame(
  census$get("census_info")$get("datasets")$read()$concat()
)

hbca_collection_id <- "283d65eb-dd53-496d-adb7-7570c7caa443"

# 2) Find the dataset_id(s) for HBCA v1.0 "All neurons"
hbca_neuron_ds <- datasets %>%
  filter(collection_id == hbca_collection_id,
         dataset_title == "All neurons")

hbca_neuron_ids <- unique(hbca_neuron_ds$dataset_id)

# 3) Query expression for exactly those dataset(s)
# (dataset_id is an obs column, so this filter is valid)
cell_filter <- paste0(
  "dataset_id %in% c(",
  paste0("'", hbca_neuron_ids, "'", collapse = ","),
  ")"
)

sce_neurons <- get_single_cell_experiment(
  census = census,
  organism = "Homo sapiens",
  measurement_name = "RNA",
  obs_value_filter = cell_filter,
  # request only useful obs columns; all must exist in obs
  obs_column_names = c("dataset_id", "cell_type", "assay", "tissue", "donor_id", "is_primary_data")
)

# optional: drop duplicates if desired (Census supports is_primary_data)
sce_neurons <- sce_neurons[, colData(sce_neurons)$is_primary_data]