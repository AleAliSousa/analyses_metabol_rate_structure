# Tüshaus 2023 Dataset EV1: extraction analogous to s1b obs preparation.
# One donor. Tissue samples are not independent donors; retain unique sample IDs.
library(readxl)
library(dplyr)
source("helpers/read_heiss_rates.R")
source("helpers/regional_robustness.R")
f <- Sys.getenv("TUSHAUS_EV1", "data_raw/Tushaus_2023/Dataset_EV1.xlsx")
if (!file.exists(f)) stop("Missing Tüshaus Dataset EV1: ", f,
  ". See metadata/anatomy/NOV14_PROTEOMICS_ATLAS.md for download and provenance.")
map <- read_study_anatomy_map("metadata/anatomy/s1d_tushaus_region_map.csv")
rates <- derive_study_heiss_rates(map) %>% select(source_region, anatomy_id, rcmr_value)
stopifnot(!anyDuplicated(map$source_region))
# Read as text: mixed annotation columns otherwise generate coercion warnings.
samples <- read_excel(f, sheet = "SampleID", col_types = "text") %>%
  transmute(sample_id = ID, source_region = `Brain region`,
            replicate_as_published = `Biological replicate`, donor_id = "Tushaus_single_donor")
stopifnot(nrow(samples) == 38L, !anyDuplicated(samples$sample_id),
          setequal(samples$source_region, map$source_region))
samples <- samples %>% group_by(source_region, replicate_as_published) %>%
  mutate(duplicate_replicate_label = n() > 1L) %>% ungroup() %>%
  left_join(map, by = "source_region")
write_analysis_csv(samples, "tables/s1d/tushaus_sample_inventory.csv")

extract_sheet <- function(sheet) {
  raw <- suppressMessages(read_excel(f, sheet = sheet, col_types = "text"))
  required <- c("ids", "ID", "Brain region", "Gene.names", "intensity", "Biological replicate", "class")
  stopifnot(all(required %in% names(raw)))
  d <- raw %>% transmute(protein_id = ids, gene_names = Gene.names, sample_id = ID,
    sheet_region = `Brain region`, sheet_replicate = `Biological replicate`,
    published_log_intensity = as.numeric(intensity),
    regional_enrichment_class_as_published = class, measure = sheet)
  stopifnot(!anyDuplicated(d[c("protein_id", "sample_id")]),
            all(is.finite(d$published_log_intensity)),
            all(d$sample_id %in% samples$sample_id))
  d <- left_join(d, samples, by = "sample_id")
  stopifnot(all(d$sheet_region == d$source_region),
            all(d$sheet_replicate == d$replicate_as_published))
  # No new imputation, gene reassignment, normalization, or exponentiation.
  # Preserve the published log scale: EV1 does not explicitly label its base.
  d
}
long <- bind_rows(extract_sheet("LFQ"), extract_sheet("iBAQ"))
stopifnot(sum(long$measure == "LFQ") == 303170L,
          sum(long$measure == "iBAQ") == 305816L)
dir.create("data_intermediate", showWarnings = FALSE)
saveRDS(long, "data_intermediate/tushaus_2023_protein_sample_long.rds", compress = "gzip")
region <- long %>% group_by(measure, protein_id, source_region) %>%
  summarise(gene_names = paste(sort(unique(na.omit(gene_names))), collapse = ";"),
    n_tissue_samples = n_distinct(sample_id),
    mean_published_log_intensity = mean(published_log_intensity),
    sd_between_tissue_samples = if (n() > 1L) sd(published_log_intensity) else NA_real_,
    .groups = "drop") %>% left_join(map, by = "source_region") %>%
  left_join(rates, by = c("source_region", "anatomy_id")) %>%
  mutate(rcmr_value = ifelse(primary_included, rcmr_value, NA_real_),
         n_donors = 1L, eligible_region_mean = n_tissue_samples >= 2L)
stopifnot(all(is.finite(region$rcmr_value[region$primary_included])),
          length(unique(region$anatomy_id[region$primary_included])) == 7L)
saveRDS(region, "data_intermediate/tushaus_2023_protein_regional_summary.rds", compress = "gzip")
coverage <- region %>% group_by(source_region, region_label, measure, primary_included, match_note) %>%
  summarise(n_proteins = n(), n_proteins_at_least_two_samples = sum(eligible_region_mean),
            .groups = "drop")
write_analysis_csv(coverage, "tables/s1d/tushaus_region_coverage.csv")
write_analysis_csv(data.frame(source_file = f, md5 = unname(tools::md5sum(f)),
  doi = "10.15252/embj.2023114665", sheets = "SampleID; LFQ; iBAQ",
  scale = "Published log-intensity units; no second logarithm or back-transformation",
  donor_n = 1L, tissue_sample_n = nrow(samples)), "tables/s1d/tushaus_input_provenance.csv")
message("Tüshaus extraction: ", nrow(long), " values; 38 tissue samples; 13 regions; 7 metabolic matches; one donor.")
