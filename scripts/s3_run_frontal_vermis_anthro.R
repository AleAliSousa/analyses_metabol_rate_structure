# ============================================================
# Study 3 — extended run: add Frontal lobe (grey) + Cerebellar vermis as
# regions 15 & 16, ANTHROPOIDS ONLY.
#
# Why anthropoids only: frontal-lobe grey (Smaers) and vermis are measured
# almost exclusively in anthropoids, so an all-primate frame would force a
# prosimian-empty sample. See metadata/PRIMARY_DATA_COMPILATION_PLAN.md for the
# source caveats (esp. Semendeferi whole-brain exclusion; Navarrete exclusion).
#
# rCMRGlc values (Heiss inventory, both hemispheres):
#   Frontal lobe = 35.3   |   Vermis = 30.1
# These are injected into the rCMRGlc join via CONFIG$extra_regions.
#
# Outputs go to run_tag = "merged_all_anthro_frontalvermis"
#   figs/s3/<tag>/  checks/s3/<tag>/  tables/s3/<tag>/
# so nothing overwrites the core 14-region runs.
#
# Run from the project root.
# ============================================================

setwd("~/Library/CloudStorage/Dropbox/COLLABORATIVE/Do expensive brain regions increase less in humans/analyses_metabol_rate_structure")

library(tidyverse)

VARIANT_SCRIPT <- "scripts/s3_predicValuesPGLS_MERGED_variant.R"

extra <- list(
  list(raw = "FrontalCortex_grey", merged_col = "FrontalCortex_grey_matter_Vol.mm3",
       label = "Frontal lobe grey",  rCMRGlc = 35.3),
  list(raw = "CerebellarVermis",   merged_col = "CerebellarVermis_Vol.mm3",
       label = "Cerebellar vermis",  rCMRGlc = 30.1)
)

# ---- Primary run: anthropoids only, 16 regions ----
.s3_config_override <<- list(
  data_source   = "merged",
  source_filter = "all",
  clade_restrict = "anthropoids",
  extra_regions = extra,
  extra_tag     = "frontalvermis"
)
# The sourced variant crops primate -> tree -> anthropoid at its first step and
# hard-stops on any non-primate. assert_primates_only() re-confirms it at this entry
# point (is_primate/data live in the global env after source()).
assert_primates_only <- function(tag) {
  if (!exists("data") || !exists("is_primate"))
    stop("[driver] variant did not populate data/is_primate for ", tag)
  bad <- unique(data$Species[!is_primate(data$Species)])
  if (length(bad) > 0)
    stop("[driver] NON-PRIMATE leaked in ", tag, ": ", paste(bad, collapse = ", "))
  message(sprintf("[driver] %s: %d primate species confirmed (0 non-primates).", tag, nrow(data)))
}

source(VARIANT_SCRIPT, local = FALSE)
anthro_tag <- CONFIG$run_tag                       # "merged_all_anthro_frontalvermis"
assert_primates_only(anthro_tag)                   # consistency guard

# ---- Optional comparator: all primates, same 16 regions ----
# (Frontal/vermis are anthropoid-dominated so their rows change little; the OTHER
#  14 regions shift. Comment out if you only want the anthropoid frame.)
.s3_config_override <<- list(
  data_source   = "merged",
  source_filter = "all",
  clade_restrict = "none",
  extra_regions = extra,
  extra_tag     = "frontalvermis"
)
source(VARIANT_SCRIPT, local = FALSE)
allprim_tag <- CONFIG$run_tag                      # "merged_all_frontalvermis"
assert_primates_only(allprim_tag)                  # consistency guard

if (exists(".s3_config_override")) rm(.s3_config_override)

# ---- Focused summary for the two NEW regions across the two frames ----
new_labels <- c("Frontal lobe grey", "Cerebellar vermis")
read_core <- function(tag) {
  f <- file.path("checks", "s3", tag, "core_with_rCMRGlc_predicted_volumes.csv")
  if (!file.exists(f)) { warning("missing: ", f); return(NULL) }
  readr::read_csv(f, show_col_types = FALSE) %>% mutate(run_tag = tag)
}
num <- function(x) readr::parse_number(as.character(x))
focus <- bind_rows(read_core(anthro_tag), read_core(allprim_tag)) %>%
  filter(Structure %in% new_labels) %>%
  mutate(across(c(Predicted, Observed, `95% CI min`, `95% CI max`,
                  Diff.pre, lambda), num),
         frame = ifelse(grepl("anthro", run_tag), "anthropoids", "all primates")) %>%
  select(Structure, frame, Model, N, lambda,
         `95% CI min`, Predicted, `95% CI max`, Observed, Diff.pre, rCMRGlc)

dir.create("tables/s3/_frontal_vermis", showWarnings = FALSE, recursive = TRUE)
write_csv(focus, "tables/s3/_frontal_vermis/frontal_vermis_predicted_vs_observed.csv")

message("\n==== Frontal lobe & Vermis: predicted vs observed (human) ====")
print(as.data.frame(focus), row.names = FALSE)
message("\nFigures/tables per frame in:")
message("  figs/s3/",   anthro_tag,  " , figs/s3/",   allprim_tag)
message("Focused table: tables/s3/_frontal_vermis/frontal_vermis_predicted_vs_observed.csv")
