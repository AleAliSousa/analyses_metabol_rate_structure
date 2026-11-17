# Across-region protein associations, on native published scales.
# Functional enrichment classes in EV1 are NOT cell classes or cell fractions.
library(dplyr)
library(tidyr)
library(ggplot2)
source("helpers/regional_robustness.R")
source("helpers/model_diagnostics.R"); diag_stage("s1d")
f <- "data_intermediate/tushaus_2023_protein_regional_summary.rds"
if (!file.exists(f)) stop("Run s1d_2a_extract_tushaus_proteomics.R first")
all_regions <- readRDS(f)
d <- all_regions %>% filter(primary_included, eligible_region_mean,
                            is.finite(mean_published_log_intensity))
# Every primary protein uses the SAME seven regions. This avoids changing
# anatomy across proteins. Partial coverage is retained in the coverage audit.
coverage <- d %>% group_by(measure, protein_id, gene_names) %>%
  summarise(n_matched_regions = n_distinct(anatomy_id), .groups = "drop") %>%
  mutate(primary_complete = n_matched_regions == 7L)
write_analysis_csv(coverage, "tables/s1d/tushaus_protein_coverage.csv")
complete <- d %>% semi_join(filter(coverage, primary_complete), by = c("measure", "protein_id"))
screen_one <- function(z) {
  x <- z$rcmr_value; y <- z$mean_published_log_intensity
  if (length(unique(y)) < 2L) return(data.frame(n = length(y), slope = NA_real_,
    p_ols = NA_real_, rho = NA_real_, p_spearman = NA_real_))
  fit <- summary(lm(y ~ x))$coefficients
  r <- suppressWarnings(cor.test(x,y,method = "spearman",exact = FALSE))
  data.frame(n = length(y), slope = fit["x","Estimate"], p_ols = fit["x","Pr(>|t|)"],
    rho = unname(r$estimate), p_spearman = r$p.value)
}
screen <- complete %>% group_by(measure, protein_id, gene_names) %>%
  group_modify(~screen_one(.x)) %>% ungroup() %>% group_by(measure) %>%
  mutate(q_ols_BH = p.adjust(p_ols,"BH"), q_spearman_BH = p.adjust(p_spearman,"BH"),
         n_tests_in_family = n(), n_donors = 1L) %>% ungroup()
write_analysis_csv(screen, "tables/s1d/tushaus_proteome_screen.csv")
# Panel selected from the existing Sjölin comparison plus the SV2A bridge.
markers <- c("GFAP", "MBP", "NEFL", "MAPT", "UCHL1", "SV2A")
panel <- d %>% filter(gene_names %in% markers)
stopifnot(!anyDuplicated(panel[c("measure","gene_names","anatomy_id")]));
write_analysis_csv(panel, "tables/s1d/tushaus_marker_regional_summary.csv")
fits <- panel %>% group_by(measure,gene_names,protein_id) %>% group_modify(~
  regional_fit(.x, "rcmr_value", "mean_published_log_intensity")) %>% ungroup() %>%
  group_by(measure) %>% mutate(q_ols_BH = p.adjust(p_ols,"BH"),
    q_spearman_BH = p.adjust(p_spearman,"BH")) %>% ungroup()
write_analysis_csv(fits,"tables/s1d/tushaus_marker_fits.csv")
for (key in unique(paste(panel$measure,panel$gene_names))) {
  z <- panel[paste(panel$measure,panel$gene_names) == key,]
  diag_lm(lm(mean_published_log_intensity ~ rcmr_value,data=z),
          paste("Tushaus", key),labels=z$region_label)
}
# Center each protein WITHIN each measure for a relative regional profile.
# This is a log-abundance deviation, not a cell fraction or protein mass share.
profile <- all_regions %>% filter(gene_names %in% markers,eligible_region_mean) %>%
  group_by(measure,protein_id) %>% mutate(centered_log_abundance =
    mean_published_log_intensity - mean(mean_published_log_intensity)) %>% ungroup()
write_analysis_csv(profile,"tables/s1d/tushaus_marker_regional_profiles.csv")
dir.create("figs/s1d",showWarnings=FALSE,recursive=TRUE)
p <- ggplot(panel,aes(rcmr_value,mean_published_log_intensity)) +
  geom_smooth(method="lm",formula=y~x,se=FALSE,color="#43677D") +
  geom_point(aes(color=region_label),size=2) +
  facet_grid(measure ~ gene_names,scales="free_y") + theme_classic(base_size=10) +
  labs(title="Tüshaus protein abundance and regional metabolism",
    subtitle="Seven anatomically matched regions; tissue-sample means from one donor",
    x="Heiss rCMRGlc (µmol / 100 g / min)",y="Published log-intensity units",color="Region",
    caption="LFQ and iBAQ are separate measurements of the same material. Regional inference is exploratory; no donor-population inference.")
q <- ggplot(profile,aes(region_label,gene_names,fill=centered_log_abundance)) +
  geom_tile() + facet_wrap(~measure,ncol=1) +
  scale_fill_gradient2(low="#3969AC",mid="white",high="#B84B3E") +
  theme_minimal(base_size=10) + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  labs(title="Protein profiles across all 13 sampled regions",x=NULL,y=NULL,
    fill="Centered\nlog abundance",caption="Each protein is centered across regions within its measurement scale. These are not cell proportions.")
for (ext in c("png","pdf")) {
  ggsave(paste0("figs/s1d/tushaus_marker_metabolism.",ext),p,width=14,height=7,dpi=180)
  ggsave(paste0("figs/s1d/tushaus_marker_profiles.",ext),q,width=12,height=7,dpi=180)
}
print(fits %>% select(measure,gene_names,n,rho,p_spearman,q_spearman_BH,loo_max_p))
message("Screened ",nrow(screen)," protein-by-measure combinations on identical seven-region coverage.")
