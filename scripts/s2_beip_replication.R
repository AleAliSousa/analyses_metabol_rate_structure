# BEIP is independent of ERA, but its 2014 and 2022 reports share a cohort.
# Within-report percent differences cancel unspecified volume units/side basis.
# These are descriptive raw-mean contrasts, not the papers' adjusted effects.
setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))
library(dplyr)
library(ggplot2)
source("helpers/read_heiss_rates.R")
source("helpers/regional_robustness.R")
source("helpers/model_diagnostics.R"); diag_stage("s2")
rd <- function(name) read.csv(paste0("data_raw/", name, ".csv"), check.names = FALSE,
                              stringsAsFactors = FALSE, na.strings = c("", "NA"))
t1 <- rd("McLaughlin_etal_2014_Table1")
t4 <- rd("McLaughlin_etal_2014_Table4")
s3 <- rd("Sheridan_etal_2022_TableS3")
early <- bind_rows(
  transmute(t4, source_region = structure_as_published, group_code, n_group,
            volume = volume_mean, source = "McLaughlin_etal_2014_Table4"),
  transmute(filter(t1, variable_as_published == "Cortical Gray Matter"),
            source_region = variable_as_published, group_code, n_group,
            volume = mean, source = "McLaughlin_etal_2014_Table1")
) %>% mutate(timepoint = "BEIP age 8-10", volume_unit = "not printed; native units")
# One bilateral value per region; do not count L, R and derived L+R together.
late <- s3 %>%
  filter(hemisphere_as_published == "Left+Right (computed)" |
           is.na(hemisphere_as_published)) %>%
  transmute(source_region = structure_as_published, group_code, n_group,
            volume = mean.cm3, source = "Sheridan_etal_2022_TableS3",
            timepoint = "BEIP age 16", volume_unit = "cm3")
stopifnot(nrow(early) == 18L, nrow(late) == 24L)
means <- bind_rows(early, late)
stopifnot(all(is.finite(means$volume)), all(means$volume > 0),
          !anyDuplicated(means[c("timepoint", "source_region", "group_code")]))
control <- means %>% filter(group_code == "NIG") %>%
  select(timepoint, source_region, control_volume = volume, n_control = n_group)
d <- means %>% filter(group_code != "NIG") %>%
  left_join(control, by = c("timepoint", "source_region")) %>%
  mutate(volume_deficit_pct = 100 * (control_volume - volume) / control_volume,
         contrast = paste(timepoint, group_code, "vs NIG"),
         cohort = "BEIP; repeated cohort across time points",
         outcome_note = "Positive = smaller exposed-group volume; raw means, not covariate-adjusted effects")
rates <- derive_study_heiss_rates("metadata/anatomy/s2_beip_region_map.csv")
relationships <- read_study_anatomy_map("metadata/anatomy/s2_beip_region_map.csv") %>%
  distinct(source_region, relationship)
d <- d %>% left_join(rates, by = "source_region") %>%
  left_join(relationships, by = "source_region")
stopifnot(all(is.finite(d$rcmr_value[d$primary_included])), all(is.finite(d$volume_deficit_pct)),
          !anyDuplicated(d[c("contrast", "anatomy_id")]))
write_analysis_csv(d, "tables/s2/beip_volume_contrasts.csv")
fits <- list()
for (contrast_name in unique(d$contrast)) {
  z <- d[d$contrast == contrast_name, ]
  sets <- list(primary = primary_region_rows(z))
  for (nm in names(sets)) {
    a <- sets[[nm]]
    fits[[paste(contrast_name, nm)]] <- cbind(contrast = contrast_name, analysis = nm,
      regional_fit(a, "rcmr_value", "volume_deficit_pct"))
    diag_lm(lm(volume_deficit_pct ~ rcmr_value, data = a),
            paste("BEIP", contrast_name, nm), labels = a$source_region)
  }
}
fits <- bind_rows(fits)
write_analysis_csv(fits, "tables/s2/beip_replication_fits.csv")
plot_data <- primary_region_rows(d) %>% mutate(
  contrast = factor(contrast, levels = unique(d$contrast)),
  region_label = case_when(
    source_region %in% c("Cortical Gray Matter", "Total grey matter volume (cubic centimeters)") ~ "Cortical grey matter",
    source_region == "Total white matter volume" ~ "White matter",
    TRUE ~ source_region
  )
)
p <- ggplot(plot_data, aes(rcmr_value, volume_deficit_pct)) +
  geom_hline(yintercept = 0, color = "grey70", linetype = 2) +
  geom_smooth(method = "lm", formula = y ~ x, color = "steelblue") +
  geom_point(size = 2.5) +
  ggrepel::geom_text_repel(aes(label = region_label), seed = 12, size = 2.7,
                          max.overlaps = Inf) +
  facet_wrap(~ contrast, nrow = 1) +
  labs(title = "BEIP replication, separate from ERA",
       subtitle = "Time points share a cohort; CAUG and FCG share NIG controls",
       x = "Heiss 2004 rCMRGlc (µmol / 100 g / min)",
       y = "Volume deficit relative to NIG (%)",
       caption = "Positive = smaller volume. Raw group means; strict anatomical matches only.\nAge-16 means rounded to 0.1 cm3; no uncertainty inferred from unspecified dispersion or unknown L/R covariance.") +
  theme_classic(base_size = 11)
dir.create("figs/s2", recursive = TRUE, showWarnings = FALSE)
ggsave("figs/s2/beip_replication.png", p, width = 13, height = 6, dpi = 180)
ggsave("figs/s2/beip_replication.pdf", p, width = 13, height = 6, device = cairo_pdf)
print(fits)
