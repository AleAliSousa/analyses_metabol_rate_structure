# ============================================================
# Study 3 — run the merged-data variant across configurations and
# compare human predicted-vs-observed volumes + per-region N.
#
# Runs s3_predicValuesPGLS_MERGED_variant.R once per configuration by setting
# .s3_config_override before each source(). Each run writes its own tagged
# outputs under figs/checks/tables/s3/<run_tag>/. This driver then reads back
# each run's predicted-volume table and builds a side-by-side comparison.
#
# Configurations:
#   stephan            data_source=stephan                     (baseline / prior results)
#   merged_all         data_source=merged, all sources         (max coverage)
#   merged_histol      data_source=merged, histological only   (Stephan-comparable; human kept)
#   merged_all_anthro  data_source=merged, all sources, anthropoids only
#
# INCLUDE_FRONTAL_VERMIS (below):
#   TRUE  -> adds Frontal lobe grey (rCMRGlc 35.3) + Cerebellar vermis (30.1) as
#            regions 15-16 to EVERY config, so they appear alongside the core 14 in
#            the comparison table. In the stephan config those columns don't exist,
#            so they drop out (shown as NA). NOTE: this also adds the two structures
#            to each config's internal rCMRGlc 12-way regression (Table_S_linear_*),
#            i.e. that regression is now fit on 16 (minus always-excluded) structures.
#   FALSE -> canonical core-14 comparison only.
# Per-region predicted/observed/N rows are independent of region count either way;
# only the pooled rCMRGlc regression is affected.
#
# Run from the project root.
# ============================================================

setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))  # repo root (portable; replaces hardcoded path -- see R/project_root.R)

library(tidyverse)
library(writexl)

VARIANT_SCRIPT <- "scripts/s3_predicValuesPGLS_MERGED_variant.R"

INCLUDE_FRONTAL_VERMIS <- TRUE

extra_fv <- list(
  list(raw = "FrontalCortex_grey", merged_col = "FrontalCortex_grey_matter_Vol.mm3",
       label = "Frontal lobe grey",  rCMRGlc = 35.3),
  list(raw = "CerebellarVermis",   merged_col = "CerebellarVermis_Vol.mm3",
       label = "Cerebellar vermis",  rCMRGlc = 30.1)
)

configs <- list(
  stephan           = list(data_source = "stephan"),
  merged_all        = list(data_source = "merged", source_filter = "all"),
  merged_histol     = list(data_source = "merged", source_filter = "histological"),
  merged_all_anthro = list(data_source = "merged", source_filter = "all",
                           clade_restrict = "anthropoids")
)

# Attach frontal/vermis to every config (own run_tag suffix keeps outputs separate)
if (isTRUE(INCLUDE_FRONTAL_VERMIS)) {
  configs <- lapply(configs, function(cfg) {
    cfg$extra_regions <- extra_fv
    cfg$extra_tag     <- "frontalvermis"
    cfg
  })
}

# ---- Run each configuration ----
# The sourced variant performs the primate -> tree -> anthropoid crop at its first
# step and hard-stops on any non-primate. We re-assert here so this entry point
# independently confirms primates-only (is_primate/grade_of/data live in the global
# env after each source()).
assert_primates_only <- function(tag) {
  if (!exists("data") || !exists("is_primate"))
    stop("[driver] variant did not populate data/is_primate for ", tag)
  bad <- unique(data$Species[!is_primate(data$Species)])
  if (length(bad) > 0)
    stop("[driver] NON-PRIMATE leaked in ", tag, ": ", paste(bad, collapse = ", "))
  message(sprintf("[driver] %s: %d primate species confirmed (0 non-primates).", tag, nrow(data)))
}

run_tags <- character(0)
for (nm in names(configs)) {
  message("\n=================  RUNNING CONFIG: ", nm, "  =================")
  .s3_config_override <<- configs[[nm]]
  source(VARIANT_SCRIPT, local = FALSE)         # populates CONFIG + crops data, writes tagged outputs
  assert_primates_only(nm)                      # consistency guard
  run_tags[nm] <- CONFIG$run_tag
}
if (exists(".s3_config_override")) rm(.s3_config_override)

# ---- Collect each run's predicted-volume table ----
read_run <- function(nm, tag) {
  f <- file.path("checks", "s3", tag, "core_with_rCMRGlc_predicted_volumes.csv")
  if (!file.exists(f)) { warning("Missing output for ", nm, ": ", f); return(NULL) }
  readr::read_csv(f, show_col_types = FALSE) %>%
    mutate(config = nm, run_tag = tag)
}
all_runs <- imap_dfr(run_tags, ~ read_run(.y, .x))

# Normalise numeric columns that may have been written as text
num <- function(x) readr::parse_number(as.character(x))
all_runs <- all_runs %>%
  mutate(across(c(Predicted, Observed, `95% CI min`, `95% CI max`,
                  Diff.pre, Diff.min, Diff.max, lambda), num),
         N = as.integer(N))

dir.create("tables/s3/_comparison", showWarnings = FALSE, recursive = TRUE)
write_csv(all_runs, "tables/s3/_comparison/all_configs_long.csv")

# ---- Side-by-side comparison (Brownian model) ----
make_wide <- function(df, model_key, value_col) {
  df %>%
    filter(grepl(model_key, Model, ignore.case = TRUE)) %>%
    select(Structure, config, !!value_col) %>%
    pivot_wider(names_from = config, values_from = !!value_col)
}

bm_predicted <- make_wide(all_runs, "brownian|\\bbm\\b", "Predicted")
bm_N         <- make_wide(all_runs, "brownian|\\bbm\\b", "N")
bm_diffpre   <- make_wide(all_runs, "brownian|\\bbm\\b", "Diff.pre")

ml_predicted <- make_wide(all_runs, "pagel|lambda|\\bml\\b", "Predicted")
ml_diffpre   <- make_wide(all_runs, "pagel|lambda|\\bml\\b", "Diff.pre")

observed_tbl <- all_runs %>%
  filter(grepl("brownian|\\bbm\\b", Model, ignore.case = TRUE)) %>%
  distinct(Structure, config, Observed) %>%
  pivot_wider(names_from = config, values_from = Observed)

write_xlsx(
  list(
    "BM_Predicted"   = bm_predicted,
    "BM_N"           = bm_N,
    "BM_DiffPre"     = bm_diffpre,
    "PagelML_Pred"   = ml_predicted,
    "PagelML_DiffPre"= ml_diffpre,
    "Observed"       = observed_tbl,
    "all_long"       = all_runs
  ),
  "tables/s3/_comparison/Study3_config_comparison.xlsx"
)

# ---- Comparison figure: prediction error by structure, per config ----
plot_df <- all_runs %>%
  filter(grepl("brownian|\\bbm\\b", Model, ignore.case = TRUE)) %>%
  mutate(config = factor(config, levels = names(configs)))

p_cmp <- ggplot(plot_df, aes(x = reorder(Structure, Diff.pre), y = Diff.pre,
                             color = config, group = config)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 2.6, position = position_dodge(width = 0.5)) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_bw() +
  labs(title = "Study 3 — BM prediction error by data configuration",
       subtitle = "Diff.pre = (Observed - Predicted)/Observed for the human. Positive = human larger than predicted.",
       x = "Brain structure", y = "Prediction error (BM)", color = "Configuration") +
  theme(legend.position = "bottom")

ggsave("tables/s3/_comparison/Study3_config_comparison_BM.png", p_cmp, width = 9, height = 6, dpi = 300)
ggsave("tables/s3/_comparison/Study3_config_comparison_BM.pdf", p_cmp, width = 9, height = 6)

# ---- Per-region N summary across configs ----
N_summary <- all_runs %>%
  filter(grepl("brownian|\\bbm\\b", Model, ignore.case = TRUE)) %>%
  select(Structure, config, N) %>%
  pivot_wider(names_from = config, values_from = N)
write_csv(N_summary, "tables/s3/_comparison/per_region_N_by_config.csv")

message("\nComparison written to tables/s3/_comparison/ (xlsx, csv, png/pdf).")
print(N_summary)
