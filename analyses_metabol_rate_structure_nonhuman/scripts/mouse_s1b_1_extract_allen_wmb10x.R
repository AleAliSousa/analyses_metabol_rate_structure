############################################################
# mouse_s1b_1 -- Extract cell-type counts per brain structure from the
# Allen Brain Cell Atlas whole-mouse-brain 10x dataset (Yao et al. 2023,
# Nature 624:317; WMB-10X, release 20231215). Mouse analogue of
# s1b_1_extract_transcriptomic (human, Siletti/Linnarsson).
#
# Inputs (downloaded if absent from
#   https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/):
#   metadata/WMB-10X/20231215/cell_metadata.csv          (~860 MB)
#   metadata/WMB-10X/20231215/region_of_interest_metadata.csv
#   metadata/WMB-taxonomy/20231215/cluster_to_cluster_annotation_membership.csv
# Output (small, committed):
#   data_intermediate/allen_wmb10x_counts_roi_donor_celltype.csv
#     one row per ROI x donor x cluster-level class/subclass/neurotransmitter
############################################################
.a <- commandArgs(FALSE); .f <- sub("^--file=", "", .a[grepl("^--file=", .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else "analyses_metabol_rate_structure_nonhuman/scripts", "common.R"))
suppressPackageStartupMessages({library(dplyr); library(readr)})

BASE <- "https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com"
RAW  <- nh("data_raw/allen_abc_mouse")
dir.create(RAW, showWarnings = FALSE, recursive = TRUE)
files <- c(
  "WMB-10X_cell_metadata.csv" = "metadata/WMB-10X/20231215/cell_metadata.csv",
  "WMB-10X_region_of_interest_metadata.csv" = "metadata/WMB-10X/20231215/region_of_interest_metadata.csv",
  "WMB-taxonomy_cluster_to_cluster_annotation_membership.csv" =
    "metadata/WMB-taxonomy/20231215/cluster_to_cluster_annotation_membership.csv"
)
for (f in names(files)) {
  dest <- file.path(RAW, f)
  if (!file.exists(dest)) {
    message("Downloading ", files[[f]])
    download.file(file.path(BASE, files[[f]]), dest, mode = "wb", quiet = TRUE)
  }
}

cm <- read_csv(file.path(RAW, "WMB-10X_cell_metadata.csv"),
  col_select = c(region_of_interest_acronym, donor_label, library_method, cluster_alias),
  show_col_types = FALSE)
mem <- read_csv(file.path(RAW, "WMB-taxonomy_cluster_to_cluster_annotation_membership.csv"),
  show_col_types = FALSE)

# Wide cluster -> class / subclass / neurotransmitter lookup
lvl <- c("neurotransmitter", "class", "subclass")
tax <- mem %>%
  filter(cluster_annotation_term_set_name %in% lvl) %>%
  select(cluster_alias, level = cluster_annotation_term_set_name,
         term = cluster_annotation_term_name) %>%
  tidyr::pivot_wider(names_from = level, values_from = term)
stopifnot(!anyNA(tax$class), !anyNA(tax$subclass), !anyDuplicated(tax$cluster_alias))

counts <- cm %>%
  inner_join(tax, by = "cluster_alias") %>%
  count(region_of_interest_acronym, donor_label, library_method,
        neurotransmitter, class, subclass, name = "n_cells")
stopifnot(sum(counts$n_cells) == nrow(cm))   # every cell has a taxonomy assignment

write_csv(counts, nh("data_intermediate/allen_wmb10x_counts_roi_donor_celltype.csv"))
message("Wrote ", nrow(counts), " rows covering ", sum(counts$n_cells), " cells")
