############################################################
# mouse_s1a_1 -- In-situ cell counts per CCF brain structure from the Allen
# whole-mouse-brain MERFISH atlas (Yao et al. 2023; C57BL6J-638850,
# release 20231215; 4.0 M cells in 59 coronal sections, one donor).
# Mouse stand-in for the human s1a stereology (counts of neurons, astrocytes,
# oligodendrocytes, microglia per structure). MERFISH is in situ, so unlike the
# 10x dissections (s1b) the neuron:glia balance is not altered by cell sorting.
#
# Inputs (downloaded if absent from the Allen ABC Atlas S3 bucket):
#   metadata/MERFISH-C57BL6J-638850-CCF/20231215/views/cell_metadata_with_parcellation_annotation.csv (1.6 GB)
#   metadata/Allen-CCF-2020/20230630/views/parcellation_term_with_counts.csv   (structure volumes)
# Outputs (small, committed):
#   data_intermediate/allen_merfish_counts_by_ccf_structure.csv
#   data_intermediate/allen_ccf_structure_volumes.csv
############################################################
.a <- commandArgs(FALSE); .f <- sub("^--file=", "", .a[grepl("^--file=", .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else "analyses_metabol_rate_structure_nonhuman/scripts", "common.R"))
suppressPackageStartupMessages({library(dplyr); library(readr)})

BASE <- "https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com"
RAW <- nh("data_raw/allen_abc_mouse"); dir.create(RAW, showWarnings = FALSE, recursive = TRUE)
files <- c(
  "MERFISH-CCF_cell_metadata_with_parcellation_annotation.csv" =
    "metadata/MERFISH-C57BL6J-638850-CCF/20231215/views/cell_metadata_with_parcellation_annotation.csv",
  "CCF_parcellation_term_with_counts.csv" =
    "metadata/Allen-CCF-2020/20230630/views/parcellation_term_with_counts.csv")
for (f in names(files)) if (!file.exists(file.path(RAW, f))) {
  message("Downloading ", files[[f]])
  download.file(file.path(BASE, files[[f]]), file.path(RAW, f), mode = "wb", quiet = TRUE)
}

cells <- read_csv(file.path(RAW, "MERFISH-CCF_cell_metadata_with_parcellation_annotation.csv"),
  col_select = c(parcellation_category, parcellation_division, parcellation_structure,
                 neurotransmitter, class, subclass), show_col_types = FALSE, na = character())
counts <- cells %>%
  count(category = parcellation_category, division = parcellation_division,
        structure = parcellation_structure, neurotransmitter, class, subclass, name = "n_cells")
stopifnot(sum(counts$n_cells) == nrow(cells))
write_csv(counts, nh("data_intermediate/allen_merfish_counts_by_ccf_structure.csv"))

vol <- read_csv(file.path(RAW, "CCF_parcellation_term_with_counts.csv"), show_col_types = FALSE) %>%
  select(acronym, name, identifier, voxel_count, volume_mm3) %>%
  distinct()   # some terms are listed under two parents with identical volumes
stopifnot(!anyDuplicated(vol$acronym))
write_csv(vol, nh("data_intermediate/allen_ccf_structure_volumes.csv"))
message("Wrote ", nrow(counts), " count rows (", nrow(cells), " cells); ", nrow(vol), " CCF terms")
