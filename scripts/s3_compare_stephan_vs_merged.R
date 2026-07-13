# ============================================================
# QC: compare the SAME structure values between Stephan_primates.csv and the
# merged volumes_wide.csv, cell by cell (species x structure).
#
# For each mapped structure it reports Stephan value, merged value, % difference,
# and the merged team(s) that supplied the value (from volumes_long.csv), so you
# can see where the merge reproduces Stephan exactly vs where another source /
# a different regional boundary / laterality changed the number.
#
# Outputs:
#   checks/s3/stephan_vs_merged/value_comparison_long.csv     (every cell)
#   checks/s3/stephan_vs_merged/value_comparison_summary.csv  (per-structure)
#
# Known definitional note: Stephan Total_insula_volume_L is LEFT insula; merged
# Insula_Vol.mm3 is BILATERAL (~2x). Insula_left_Vol.mm3 matches 1:1 and is what
# the Study 3 merged variant uses. This QC compares against the LEFT column.
#
# Run from the project root.
# ============================================================

setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))  # repo root (portable; replaces hardcoded path -- see R/project_root.R)
library(tidyverse)

out_dir <- "checks/s3/stephan_vs_merged"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

xwalk_species <- c("Lagothrix_lagothricha" = "Lagothrix_lagotricha",
                   "Pongo_sp."             = "Pongo_pygmaeus")
norm_sp <- function(x) {
  x <- gsub(" ", "_", x)
  hit <- x %in% names(xwalk_species); x[hit] <- unname(xwalk_species[x[hit]]); x
}

S <- read.csv("data_raw/Stephan_primates.csv", stringsAsFactors = FALSE)
S <- S[!is.na(S$Species) & trimws(S$Species) != "", ]
S$Species <- norm_sp(S$Species)
W <- read.csv("data_intermediate/volumes_wide.csv", check.names = FALSE, stringsAsFactors = FALSE)
W$Species <- norm_sp(W$Species)
L <- read.csv("data_intermediate/volumes_long.csv", check.names = FALSE, stringsAsFactors = FALSE)
L$Species <- norm_sp(L$Species)

# Stephan column -> merged column crosswalk (structures present in both)
cross <- c(
  Total_brain_net_volume = "Total_brain_net_volume_Vol.mm3",
  Medulla_oblongata = "Medulla_oblongata_Vol.mm3", Cerebellum = "Cerebellum_Vol.mm3",
  Mesencephalon = "Mesencephalon_Vol.mm3", Diencephalon = "Diencephalon_Vol.mm3",
  Telencephalon = "Telencephalon_Vol.mm3", Bulbus_olfactorius = "Bulbus_olfactorius_Vol.mm3",
  Bulbus_olfactorius_accessorius = "Bulbus_olfactorius_accessorius_Vol.mm3",
  Lobus_piriformis = "Lobus_piriformis_Vol.mm3", Septum = "Septum_Vol.mm3",
  Striatum = "Striatum_Vol.mm3", Schizocortex = "Schizo_cortex_Vol.mm3",
  Hippocampus = "Hippocampus_Vol.mm3", Epithalamus = "Epithalamus_Vol.mm3",
  Thalamus = "Thalamus_Vol.mm3", Hypothalamus = "Hypothalamus_Vol.mm3",
  Subthalamus = "Subthalamus_Vol.mm3", Pallidum = "Pallidum_Vol.mm3",
  Nucleus_subthalamicus = "Nucleus_subthalamicus_Vol.mm3", Capsula_interna = "Capsula_interna_Vol.mm3",
  Tractus_opticus = "Tractus_opticus_Vol.mm3", Palaeocortex = "Palaeocortex_Vol.mm3",
  Amygdala = "Amygdala_Vol.mm3", NeoG_Frahm = "Neocortex_grey_matter_Vol.mm3",
  NeoW_Frahm = "Neocortex_white_matter_Vol.mm3", ASG_Sousa = "Area_striata_grey_matter_Vol.mm3",
  LGN_Sousa = "Corpus_geniculatum_laterale_Vol.mm3",
  Lateral_cerebellar_nuclei = "Lateral_cerebellar_nuclei_Vol.mm3",
  Interpositus_cerebellar_nuclei = "Interpositus_cerebellar_nuclei_Vol.mm3",
  Medial_cerebellar_nuclei = "Medial_cerebellar_nuclei_Vol.mm3",
  Total_insula_volume_L = "Insula_left_Vol.mm3",           # LEFT insula (1:1 with Stephan)
  Complexus_centromedialis = "Complexus_centromedialis_Vol.mm3",
  Nucleus_tractus_olfactorius = "Nucleus_tractus_olfactorius_Vol.mm3"
)

team_lookup <- L %>% distinct(Species, Variable, Teams)

long <- imap_dfr(cross, function(wcol, scol) {
  if (!scol %in% names(S) || !wcol %in% names(W)) return(NULL)
  sdat <- tibble(Species = S$Species, Stephan_value = suppressWarnings(as.numeric(S[[scol]])))
  wdat <- tibble(Species = W$Species, Merged_value  = suppressWarnings(as.numeric(W[[wcol]])))
  inner_join(sdat, wdat, by = "Species") %>%
    filter(!is.na(Stephan_value), !is.na(Merged_value)) %>%
    mutate(Stephan_col = scol, Merged_col = wcol,
           abs_diff = Merged_value - Stephan_value,
           pct_diff = 100 * (Merged_value - Stephan_value) / Stephan_value) %>%
    left_join(team_lookup %>% filter(Variable == wcol) %>% select(Species, merged_teams = Teams),
              by = "Species")
})

long <- long %>%
  mutate(match_cat = cut(abs(pct_diff), c(-Inf, 0.5, 2, 5, 100, Inf),
                         labels = c("exact(<0.5%)","close(0.5-2%)","minor(2-5%)",
                                    "notable(5-100%)","gross(>100%)"))) %>%
  select(Species, Stephan_col, Merged_col, Stephan_value, Merged_value,
         abs_diff, pct_diff, match_cat, merged_teams)

write_csv(long, file.path(out_dir, "value_comparison_long.csv"))

summary_tbl <- long %>%
  group_by(Stephan_col, Merged_col) %>%
  summarise(n = n(),
            n_exact_lt0.5pct = sum(abs(pct_diff) < 0.5),
            n_gt5pct = sum(abs(pct_diff) > 5),
            median_abs_pct = median(abs(pct_diff)),
            max_abs_pct = max(abs(pct_diff)),
            pearson_r = cor(Stephan_value, Merged_value),
            .groups = "drop") %>%
  arrange(desc(median_abs_pct))

write_csv(summary_tbl, file.path(out_dir, "value_comparison_summary.csv"))

message("cells compared: ", nrow(long),
        " | exact(<0.5%): ", sum(abs(long$pct_diff) < 0.5),
        " | >5%: ", sum(abs(long$pct_diff) > 5))
print(summary_tbl, n = 40)
