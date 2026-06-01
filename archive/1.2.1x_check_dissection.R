library(dplyr)
library(stringr)
library(tidyr)
library(purrr)

# ------------------------------------------------------------
# Investigate rows where dissection is blank
# Focus on ROIs with any blank-dissection rows
# ------------------------------------------------------------

# 1) define blank dissection rows
obs2 <- obs %>%
  mutate(
    dissection_chr = as.character(dissection),
    tissue_chr = as.character(tissue),
    roi_chr = as.character(roi),
    dissection_blank = is.na(dissection_chr) | str_trim(dissection_chr) == "",
    tissue_blank = is.na(tissue_chr) | str_trim(tissue_chr) == ""
  )

cat("Total rows:", nrow(obs2), "\n")
cat("Blank dissection rows:", sum(obs2$dissection_blank), "\n")
cat("Nonblank dissection rows:", sum(!obs2$dissection_blank), "\n\n")

# 2) basic completeness checks for key columns
check_cols <- c(
  "dissection", "tissue", "roi", "ROIGroup", "ROIGroupCoarse", "ROIGroupFine",
  "sample_id", "supercluster_term", "cluster_id", "subcluster_id",
  "cell_type", "obs_dataset"
)
check_cols <- check_cols[check_cols %in% colnames(obs2)]

cat("=== Completeness checks ===\n")
for (col in check_cols) {
  x <- as.character(obs2[[col]])
  cat("\n---", col, "---\n")
  cat("all rows: missing =", sum(is.na(x)),
      " blank =", sum(!is.na(x) & str_trim(x) == ""), "\n")
  cat("blank dissection rows only: missing =", sum(is.na(x[obs2$dissection_blank])),
      " blank =", sum(!is.na(x[obs2$dissection_blank]) & str_trim(x[obs2$dissection_blank]) == ""), "\n")
}

# 3) sample blank-dissection rows
cat("\n=== Sample blank-dissection rows ===\n")
print(
  obs2 %>%
    filter(dissection_blank) %>%
    select(any_of(check_cols)) %>%
    head(30)
)

# 4) ROIs with blank dissection
cat("\n=== ROI counts among blank-dissection rows ===\n")
roi_blank_counts <- obs2 %>%
  filter(dissection_blank) %>%
  count(roi, sort = TRUE)

print(roi_blank_counts, n = 100)

cat("\n=== ROI summary: fraction blank dissection ===\n")
roi_blank_summary <- obs2 %>%
  group_by(roi) %>%
  summarise(
    n = n(),
    n_blank = sum(dissection_blank),
    frac_blank = n_blank / n,
    tissue_values = n_distinct(tissue, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_blank > 0) %>%
  arrange(desc(frac_blank), desc(n_blank))

print(roi_blank_summary, n = 100)

# 5) define problematic ROIs automatically as any ROI with blank dissection rows
problem_rois <- roi_blank_summary$roi

cat("\n=== Problematic ROIs ===\n")
print(problem_rois)

# 6) tissue terms among blank-dissection rows
if ("tissue" %in% colnames(obs2)) {
  cat("\n=== Tissue values among blank-dissection rows ===\n")
  tissue_blank_counts <- obs2 %>%
    filter(dissection_blank) %>%
    count(tissue, sort = TRUE)
  print(tissue_blank_counts, n = 100)
  
  cat("\n=== Cross-tab: blank dissection vs blank tissue ===\n")
  print(table(dissection_blank = obs2$dissection_blank,
              tissue_blank = obs2$tissue_blank,
              useNA = "ifany"))
}

# 7) ROIGroupCoarse among blank-dissection rows
if ("ROIGroupCoarse" %in% colnames(obs2)) {
  cat("\n=== ROIGroupCoarse among blank-dissection rows ===\n")
  print(
    obs2 %>%
      filter(dissection_blank) %>%
      count(ROIGroupCoarse, sort = TRUE),
    n = 100
  )
}

# 8) sample_id / dataset / supercluster / cell_type distributions among blank-dissection rows
for (col in c("sample_id", "obs_dataset", "supercluster_term", "cell_type")) {
  if (col %in% colnames(obs2)) {
    cat("\n===", col, "among blank-dissection rows ===\n")
    print(
      obs2 %>%
        filter(dissection_blank) %>%
        count(.data[[col]], sort = TRUE),
      n = 100
    )
  }
}

# 9) columns nearly constant among blank-dissection rows
near_constant_blank <- map_dfr(names(obs2), function(col) {
  x <- as.character(obs2[[col]][obs2$dissection_blank])
  x[is.na(x)] <- "<NA>"
  tab <- sort(table(x), decreasing = TRUE)
  tibble(
    column = col,
    n_unique = length(tab),
    top_value = names(tab)[1],
    top_prop = as.numeric(tab[1]) / length(x)
  )
})

cat("\n=== Columns nearly constant among blank-dissection rows ===\n")
print(
  near_constant_blank %>%
    arrange(desc(top_prop), n_unique),
  n = 100
)

# 10) numeric summaries by blank/nonblank dissection
num_cols <- intersect(
  c("fraction_mitochondrial", "fraction_unspliced", "cell_cycle_score", "total_genes", "total_UMIs"),
  colnames(obs2)
)

if (length(num_cols) > 0) {
  cat("\n=== Numeric summaries by blank vs nonblank dissection ===\n")
  print(
    obs2 %>%
      group_by(dissection_blank) %>%
      summarise(
        across(
          all_of(num_cols),
          list(
            mean = ~mean(.x, na.rm = TRUE),
            median = ~median(.x, na.rm = TRUE)
          ),
          .names = "{.col}_{.fn}"
        )
      )
  )
}

# ------------------------------------------------------------
# Detailed investigation of ROIs with blank dissection
# ------------------------------------------------------------

# 11) compare blank vs nonblank within each problematic ROI
cat("\n=== ROI-level summary: blank vs nonblank dissection ===\n")
roi_summary <- obs2 %>%
  filter(roi %in% problem_rois) %>%
  group_by(roi, dissection_blank) %>%
  summarise(
    n = n(),
    n_tissue = n_distinct(tissue),
    n_sample = n_distinct(sample_id),
    n_supercluster = n_distinct(supercluster_term),
    n_cell_type = n_distinct(cell_type),
    .groups = "drop"
  )

print(roi_summary, n = 100)

# 12) sample-level concentration of blanks within each ROI
cat("\n=== Sample-level blank fraction within each problematic ROI ===\n")
roi_sample_blank <- obs2 %>%
  filter(roi %in% problem_rois) %>%
  count(roi, sample_id, dissection_blank, sort = TRUE) %>%
  pivot_wider(
    names_from = dissection_blank,
    values_from = n,
    values_fill = 0,
    names_prefix = "blank_"
  ) %>%
  mutate(frac_blank = blank_TRUE / (blank_TRUE + blank_FALSE)) %>%
  arrange(roi, desc(frac_blank), desc(blank_TRUE))

print(roi_sample_blank, n = 200)

# 13) supercluster-level concentration of blanks within each ROI
if ("supercluster_term" %in% colnames(obs2)) {
  cat("\n=== Supercluster-level blank fraction within each problematic ROI ===\n")
  roi_supercluster_blank <- obs2 %>%
    filter(roi %in% problem_rois) %>%
    count(roi, supercluster_term, dissection_blank, sort = TRUE) %>%
    pivot_wider(
      names_from = dissection_blank,
      values_from = n,
      values_fill = 0,
      names_prefix = "blank_"
    ) %>%
    mutate(frac_blank = blank_TRUE / (blank_TRUE + blank_FALSE)) %>%
    arrange(roi, desc(frac_blank), desc(blank_TRUE))
  
  print(roi_supercluster_blank, n = 200)
}

# 14) cell_type-level concentration of blanks within each ROI
if ("cell_type" %in% colnames(obs2)) {
  cat("\n=== Cell-type-level blank fraction within each problematic ROI ===\n")
  roi_celltype_blank <- obs2 %>%
    filter(roi %in% problem_rois) %>%
    count(roi, cell_type, dissection_blank, sort = TRUE) %>%
    pivot_wider(
      names_from = dissection_blank,
      values_from = n,
      values_fill = 0,
      names_prefix = "blank_"
    ) %>%
    mutate(frac_blank = blank_TRUE / (blank_TRUE + blank_FALSE)) %>%
    arrange(roi, desc(frac_blank), desc(blank_TRUE))
  
  print(roi_celltype_blank, n = 200)
}

# 15) what nonblank dissection values exist for the same problematic ROIs?
cat("\n=== Nonblank dissection values within problematic ROIs ===\n")
roi_nonblank_dissection <- obs2 %>%
  filter(roi %in% problem_rois, !dissection_blank) %>%
  count(roi, dissection, sort = TRUE)

print(roi_nonblank_dissection, n = 200)

# 16) compare tissue vs dissection within problematic ROIs
cat("\n=== ROI x tissue x dissection combinations ===\n")
roi_tissue_dissection <- obs2 %>%
  filter(roi %in% problem_rois) %>%
  count(roi, tissue, dissection, sort = TRUE)

print(roi_tissue_dissection, n = 300)

# 17) compact per-sample summary
cat("\n=== Per-sample summary within problematic ROIs ===\n")
roi_sample_summary <- obs2 %>%
  filter(roi %in% problem_rois) %>%
  group_by(roi, sample_id) %>%
  summarise(
    n = n(),
    n_blank = sum(dissection_blank),
    frac_blank = n_blank / n,
    tissue = paste(unique(as.character(tissue)), collapse = "; "),
    dissection_values = paste(unique(as.character(dissection)), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(roi, desc(frac_blank), desc(n_blank))

print(roi_sample_summary, n = 200)

# 18) suspicious variant-name inspection
variant_rois <- c("Human Gpe", "Human GPe", "Human A35-36", "Human A35-A36")

cat("\n=== Variant ROI inspection ===\n")
variant_summary <- obs2 %>%
  filter(roi %in% variant_rois) %>%
  group_by(roi) %>%
  summarise(
    n = n(),
    n_blank = sum(dissection_blank),
    frac_blank = n_blank / n,
    tissues = paste(unique(as.character(tissue)), collapse = "; "),
    dissection_values = paste(unique(as.character(dissection)), collapse = "; "),
    .groups = "drop"
  )

print(variant_summary, n = Inf)

# 19) export useful tables
write.csv(
  obs2 %>% filter(dissection_blank),
  "blank_dissection_rows.csv",
  row.names = FALSE
)

write.csv(
  roi_blank_summary,
  "roi_blank_summary.csv",
  row.names = FALSE
)

write.csv(
  roi_sample_blank,
  "roi_sample_blank_summary.csv",
  row.names = FALSE
)

write.csv(
  roi_nonblank_dissection,
  "roi_nonblank_dissection_values.csv",
  row.names = FALSE
)

write.csv(
  roi_tissue_dissection,
  "roi_tissue_dissection_combinations.csv",
  row.names = FALSE
)

write.csv(
  roi_sample_summary,
  "roi_sample_summary.csv",
  row.names = FALSE
)

write.csv(
  variant_summary,
  "roi_variant_summary.csv",
  row.names = FALSE
)

cat("\nWrote:\n")
cat(" - blank_dissection_rows.csv\n")
cat(" - roi_blank_summary.csv\n")
cat(" - roi_sample_blank_summary.csv\n")
cat(" - roi_nonblank_dissection_values.csv\n")
cat(" - roi_tissue_dissection_combinations.csv\n")
cat(" - roi_sample_summary.csv\n")
cat(" - roi_variant_summary.csv\n")