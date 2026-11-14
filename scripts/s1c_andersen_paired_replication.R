# Convergence with both tracers measured in the same subjects and atlas ROIs.
# Native GM-normalized uptake only; these values are not absolute rCMRGlc/Bmax.
library(ggplot2)
source("helpers/regional_robustness.R")
d <- read.csv("data_raw/Andersen_etal_2022_TableS1.csv", check.names = FALSE)
stopifnot(all(d$is_healthy_control == 1), nrow(d) == 26,
          !anyDuplicated(d[c("group_as_published", "region_as_published")]))
d$anatomy_id <- d$region_as_published
labels <- c(FCX = "Frontal cortex", Amg = "Amygdala", MTCx = "Medial temporal cortex",
            LTCx = "Lateral temporal cortex", PCX = "Parietal cortex", OCX = "Occipital cortex",
            STR = "Striatum", SN = "Substantia nigra", ICX = "Insular cortex",
            ACC = "Anterior cingulate", PCCx = "Posterior cingulate", CBL = "Cerebellum", Thal = "Thalamus")
d$region_label <- unname(labels[d$region_as_published])
stopifnot(!anyNA(d$region_label))
d$analysis_note <- "Same-subject paired tracer ROIs; GM-normalized regional means; ages fitted separately; exploratory across-region association"
write_analysis_csv(d, "tables/s1c/andersen_paired_regions.csv")
fits <- list()
for (g in unique(d$group_as_published)) {
 z <- d[d$group_as_published == g, ]
 for (scope in c("all_paired_regions", "without_substantia_nigra")) {
  a <- if (scope == "all_paired_regions") z else z[z$anatomy_id != "SN", ]
  fits[[paste(g, scope)]] <- cbind(group = g, scope = scope,
    regional_fit(a, "FDG_GMnorm", "UCBJ_GMnorm"))
 }
}
f <- do.call(rbind, fits)
write_analysis_csv(f, "tables/s1c/andersen_paired_fits.csv")
p <- ggplot(d, aes(FDG_GMnorm, UCBJ_GMnorm)) +
 geom_smooth(method = "lm", formula = y ~ x, colour = "steelblue") +
 geom_point(size = 2.5) +
 ggrepel::geom_text_repel(aes(label = region_label), seed = 12, size = 2.8, max.overlaps = Inf) +
 facet_wrap(~group_as_published) +
 labs(title = "Paired FDG–SV2A convergence: Andersen 2022",
      subtitle = "Same atlas regions for both tracers, including thalamus; 11 young and 15 older adults",
      x = "FDG uptake (GM normalized)", y = "UCB-J uptake (GM normalized)",
      caption = "Table S1 rounded group means; 13 regions per age group. Both scans partial-volume corrected.\nRegional associations, not subject-level tests; shared GM normalization and spatial dependence limit inference.") +
 theme_classic(base_size = 12)
ggsave("figs/s1c/andersen_paired_replication.png", p, width = 12, height = 6, dpi = 180)
ggsave("figs/s1c/andersen_paired_replication.pdf", p, width = 12, height = 6, device = cairo_pdf)
print(f)
