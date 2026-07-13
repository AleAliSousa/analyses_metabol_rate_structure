# ============================================================
# Study 3 — Phase 1.1 diagnostic: missingness + clade composition
#
# Purpose: quantify the slope-comparability concern. Different regions are
# measured in different species subsets, and those subsets differ in primate
# grade composition (e.g. insula is ape-heavy / prosimian-light, cerebellum the
# opposite). PGLS does not fully remove this sampling confound, so we (a) report
# per-region N + grade composition and (b) show how a common (listwise) species
# set shrinks as regions are added — the "common-N cliff".
#
# Outputs:
#   checks/s3/phase1/s3_region_missingness_clade.csv
#   checks/s3/phase1/s3_common_N_cliff.csv
#   figs/s3/phase1/plot_clade_composition_by_region.png/.pdf
#   figs/s3/phase1/plot_common_N_cliff.png/.pdf
# ============================================================

setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))  # repo root (portable; replaces hardcoded path -- see R/project_root.R)

library(tidyverse)

dir.create("checks/s3/phase1", showWarnings = FALSE, recursive = TRUE)
dir.create("figs/s3/phase1",   showWarnings = FALSE, recursive = TRUE)

# --- Data ---
Stephan <- read.csv("data_raw/Stephan_primates.csv") %>%
  filter(!is.na(Species), trimws(Species) != "") %>%   # drop trailing blank row (59 real species)
  mutate(Preferred_brain_volume = coalesce(Brain_volume, Brainvol, Total_brain_net_volume))

# Modelled regions: raw column -> display label (same order as the main script)
region_labels <- c(
  LGN_Sousa = "Corpus geniculatum laterale", Amygdala = "Amygdala", Pallidum = "Pallidum",
  NeoW_Frahm = "Neocortex white", Total_insula_volume_L = "Insular cortex (grey)",
  Nucleus_subthalamicus = "Nucleus subthalamicus Luysi", Capsula_interna = "Capsula interna",
  Striatum = "Striatum", ASG_Sousa = "Area striata grey", NeoG_Frahm = "Neocortex grey",
  Mesencephalon = "Mesencephalon", Cerebellum = "Cerebellum", Hippocampus = "Hippocampus"
)
target_cols <- names(region_labels)

# --- Primate grade lookup (genus -> grade) ---
grade_of <- function(species) {
  genus <- sub("_.*$", "", species)
  strepsirrhini <- c("Avahi","Cheirogaleus","Daubentonia","Eulemur","Galago","Galagoides",
                     "Indri","Lepilemur","Loris","Microcebus","Nycticebus","Otolemur",
                     "Perodicticus","Propithecus","Varecia")
  platyrrhini   <- c("Alouatta","Aotus","Ateles","Callithrix","Cebus","Lagothrix",
                     "Leontopithecus","Pithecia","Plecturocebus","Saguinus","Saimiri")
  cercopith     <- c("Cercocebus","Cercopithecus","Colobus","Erythrocebus","Lophocebus",
                     "Macaca","Miopithecus","Nasalis","Papio","Piliocolobus","Pygathrix",
                     "Trachypithecus")
  hominoidea    <- c("Gorilla","Homo","Hylobates","Pan","Pongo","Symphalangus")
  dplyr::case_when(
    genus == "Tarsius"        ~ "Tarsiiformes",
    genus %in% strepsirrhini  ~ "Strepsirrhini",
    genus %in% platyrrhini    ~ "Platyrrhini",
    genus %in% cercopith      ~ "Cercopithecoidea",
    genus %in% hominoidea     ~ "Hominoidea",
    TRUE                      ~ "Unassigned"
  )
}
clade_levels <- c("Strepsirrhini","Tarsiiformes","Platyrrhini","Cercopithecoidea","Hominoidea")
Stephan$grade <- factor(grade_of(Stephan$Species), levels = clade_levels)
N_total <- nrow(Stephan)

# ============================================================
# Table 1 — per-region missingness + clade composition
# ============================================================
region_clade <- map_dfr(target_cols, function(col) {
  present <- Stephan[!is.na(Stephan[[col]]), ]
  counts  <- as.list(table(factor(present$grade, levels = clade_levels)))
  tibble(raw_column = col, Structure = region_labels[[col]],
         N_present = nrow(present), N_missing = N_total - nrow(present),
         pct_missing = round(100 * (N_total - nrow(present)) / N_total, 1)) %>%
    bind_cols(as_tibble(counts))
})
write_csv(region_clade, "checks/s3/phase1/s3_region_missingness_clade.csv")

# Stacked-bar figure (clade composition per region, ordered by N)
clade_long <- region_clade %>%
  pivot_longer(all_of(clade_levels), names_to = "grade", values_to = "n") %>%
  mutate(grade = factor(grade, levels = clade_levels),
         Structure = fct_reorder(Structure, N_present))

p_clade <- ggplot(clade_long, aes(x = Structure, y = n, fill = grade)) +
  geom_col() +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  theme_bw(base_size = 13) +
  labs(title = "Species sampling differs by region (Phase 1.1 concern)",
       subtitle = "Counts of each primate grade contributing to each region's slope",
       x = NULL, y = "Number of species", fill = "Primate grade")
ggsave("figs/s3/phase1/plot_clade_composition_by_region.png", p_clade, width = 9, height = 5.5, dpi = 300)
ggsave("figs/s3/phase1/plot_clade_composition_by_region.pdf", p_clade, width = 9, height = 5.5)

# ============================================================
# Table 2 — common-N cliff (regions added in order of decreasing N)
# ============================================================
ordered_cols <- region_clade %>% arrange(desc(N_present)) %>% pull(raw_column)

common <- Stephan$Species
cliff <- map_dfr(seq_along(ordered_cols), function(i) {
  col <- ordered_cols[i]
  present <- Stephan$Species[!is.na(Stephan[[col]])]
  common <<- intersect(common, present)
  comp <- table(factor(grade_of(common), levels = clade_levels))
  tibble(step = i, region_added = region_labels[[col]],
         region_N = length(present), cumulative_common_N = length(common),
         common_species_clades = paste(sprintf("%s:%d", names(comp), as.integer(comp))[comp > 0],
                                       collapse = "; "))
})
write_csv(cliff, "checks/s3/phase1/s3_common_N_cliff.csv")

p_cliff <- ggplot(cliff, aes(reorder(region_added, step), cumulative_common_N, group = 1)) +
  geom_line(color = "steelblue4", linewidth = 1) +
  geom_point(size = 3, color = "steelblue4") +
  geom_text(aes(label = cumulative_common_N), vjust = -0.8, size = 3.3) +
  coord_cartesian(ylim = c(0, max(cliff$cumulative_common_N) + 5)) +
  theme_bw(base_size = 13) +
  labs(title = "Common-N cliff: listwise-complete species as regions are added",
       subtitle = "Regions added in order of decreasing coverage. Big drops mark the binding regions.",
       x = "Region added (cumulative)", y = "Species complete across all added regions") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave("figs/s3/phase1/plot_common_N_cliff.png", p_cliff, width = 9, height = 5.5, dpi = 300)
ggsave("figs/s3/phase1/plot_common_N_cliff.pdf", p_cliff, width = 9, height = 5.5)

message("Phase 1.1 diagnostic written to checks/s3/phase1/ and figs/s3/phase1/")
print(region_clade)
print(cliff)
