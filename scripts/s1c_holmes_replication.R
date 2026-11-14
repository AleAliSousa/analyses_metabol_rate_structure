# Independent healthy-control SV2A BPND replication; never calibrated onto Bmax.
setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))
library(ggplot2)
source("helpers/read_heiss_rates.R")
source("helpers/regional_robustness.R")
source("helpers/model_diagnostics.R"); diag_stage("s1c")

raw <- read.csv("data_raw/Holmes_etal_2024_SupplementaryTable2.csv",
                check.names = FALSE, stringsAsFactors = FALSE)
d <- raw[raw$is_healthy_control == 1 & raw$measure == "BPND", ]
stopifnot(nrow(d) == 26L, all(d$group == "HC"), all(d$n == 30))
map <- read_study_anatomy_map("metadata/anatomy/s1c_holmes_region_map.csv")
rates <- derive_study_heiss_rates(map)
d <- merge(d, rates, by.x = "region_as_published", by.y = "source_region")
metadata <- unique(map[c("source_region", "relationship", "mg_pvc_applied")])
d <- merge(d, metadata, by.x = "region_as_published", by.y = "source_region")
stopifnot(!anyDuplicated(d$anatomy_id), nrow(d) == 13L)
stopifnot(all(!d$primary_included | is.finite(d$rcmr_value)))
d$analysis_note <- "HC only; BPND is not Bmax; mixed regional partial-volume correction"
write_analysis_csv(d, "tables/s1c/holmes_rcmrglc_matched_regions.csv")
unmapped <- raw[raw$is_healthy_control == 1 & raw$measure == "BPND" &
                 !raw$region_as_published %in% d$region_as_published, ]
unmapped$exclusion_reason <- "No equivalent Heiss ROI or nested cortical parcel; ventral striatum is broader than accumbens"
write_analysis_csv(unmapped, "tables/s1c/holmes_unmatched_regions.csv")

sets <- list(primary = primary_region_rows(d),
             pvc_corrected_only = subset(primary_region_rows(d), mg_pvc_applied == "yes"))
fits <- do.call(rbind, lapply(names(sets), function(nm) {
  z <- sets[[nm]]
  diag_lm(lm(mean ~ rcmr_value, data = z), paste("Holmes BPND", nm),
          labels = z$region_as_published)
  cbind(analysis = nm, regional_fit(z, "rcmr_value", "mean"))
}))
write_analysis_csv(fits, "tables/s1c/holmes_replication_fits.csv")
p <- ggplot(sets$primary, aes(rcmr_value, mean)) +
  geom_smooth(method = "lm", formula = y ~ x, color = "steelblue", se = TRUE) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0, alpha = .4) +
  geom_point(aes(shape = mg_pvc_applied), size = 2.8) +
  ggrepel::geom_text_repel(aes(label = region_as_published), seed = 12,
                          max.overlaps = Inf, size = 3) +
  labs(title = "Independent SV2A replication: Holmes healthy controls",
       subtitle = "30 participants; one point per region; strict anatomical matches",
       x = "Heiss 2004 rCMRGlc (µmol / 100 g / min)", y = "SV2A BPND (unitless)",
       shape = "Partial-volume\ncorrection",
       caption = "BPND remains on its native scale. Thalamus and cerebellum excluded for uncertain matching.\nError bars: between-participant SD; fitted band: region-level OLS 95% CI.") +
  theme_classic(base_size = 12)
dir.create("figs/s1c", recursive = TRUE, showWarnings = FALSE)
ggsave("figs/s1c/holmes_replication.png", p, width = 9, height = 6.5, dpi = 180)
ggsave("figs/s1c/holmes_replication.pdf", p, width = 9, height = 6.5, device = cairo_pdf)
print(fits)
