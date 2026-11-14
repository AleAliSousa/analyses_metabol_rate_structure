# Siletti et al. 2023 (CELLxGENE revision of 2026-08-10) cell_type handling.
#
# In the original release every cell in the "All neurons" dataset carried
# cell_type == "neuron" (2,477,203 cells) or "leukocyte" (3,753 cells, all in
# supercluster "Miscellaneous"). The revised release replaces these with
# specific Cell Ontology terms, so cell_type == "neuron" now matches almost
# nothing. Neurons are therefore defined as every cell that is not one of the
# four immune terms that replaced "leukocyte".
SILETTI_NEURONAL_FILE_IMMUNE_TYPES <- c(
  "T cell", "natural killer cell", "monocyte", "B cell"
)

is_siletti_neuron <- function(cell_type) {
  !(as.character(cell_type) %in% SILETTI_NEURONAL_FILE_IMMUNE_TYPES)
}
