## Install and Load up packages
library(here)
library(anndataR)
library(ggplot2)
library(ggpmisc)
library(readxl)
library(tidyverse)

# Paths relative to project root:
r_dir      <- here("R")
data_dir   <- here("data")
output_dir <- here("output")

# https://github.com/linnarsson-lab/adult-human-brain/
# https://datasets.cellxgene.cziscience.com/a71efd3c-765c-466b-8eca-0b29024094d4.h5ad
# read anndata object from h5ad file
adata <- read_h5ad(
  "~/Library/CloudStorage/OneDrive-AllenInstitute/Analysis_region_celltype_human/a71efd3c-765c-466b-8eca-0b29024094d4.h5ad",
  as = "HDF5AnnData"
)

obs <- as.data.frame(adata$obs)

# inspect the anndata object
colnames(obs)

# filter to neurons only
obs_neuron <- obs %>%
  filter(cell_type == "neuron")

#####################################
# Anatomical grouping of dissections
#####################################

# --- NEW: derive cortical lobe for cerebral cortex dissections ---
anatomy_rules <- tribble(
  ~pattern, ~anatomy_group,
  
  # ===================== CORTEX =====================
  
  "^Cerebral cortex \\(Cx\\).*cuneus|lingual|occipital|peristriate|prostriata|visual|V1|V2|A19|MT",
  "Occipital lobe",
  
  "^Cerebral cortex \\(Cx\\).*temporal|parahippocampal|entorhinal|perirhinal|fusiform|occipitotemporal|hippocamp",
  "Temporal lobe",
  
  "^Cerebral cortex \\(Cx\\).*frontal|precentral|orbit|gyrus rectus|subcallosal|A44|A45|A46|A13|A14|A24|A25|A32",
  "Frontal lobe",
  
  "^Cerebral cortex \\(Cx\\).*parietal|postcentral|supramarginal|supraparietal|precuneus|operculum|A40|A43|A5|A7|S1",
  "Parietal lobe",
  
  "^Cerebral cortex \\(Cx\\).*insula|insular",
  "Insular lobe",
  
  "^Cerebral cortex \\(Cx\\).*cingulate|retrosplenial|A23|A29|A30",
  "Cingulate / limbic",
  
  # ===================== LIMBIC =====================
  
  "Amygdaloid|Amygdala|Basolateral nuclear group|Central nuclear group|Corticomedial nuclear group",
  "Corpus amygdaloideum",
  
  "Hippocampal CA|Dentate gyrus",
  "Hippocampus",
  
  # ===================== BASAL GANGLIA =====================
  
  "Caudate",
  "Caudatum",
  
  "Putamen",
  "Putamen",
  
  "Accumbens",
  "Nucleus accumbens",
  
  "Pallidum|Globus pallidus",
  "Pallidum",
  
  # ===================== THALAMUS =====================
  
  "Medial thalam",
  "Nucleus medial thalami",
  
  "Lateral geniculate",
  "Corpus geniculatum laterale",
  
  "Medial geniculate",
  "Corpus geniculatum mediale",
  
  "Subthalamic",
  "Nucleus subthalamicus",
  
  # ===================== MIDBRAIN =====================
  
  "Inferior colliculus",
  "Colliculus inferior",
  
  "Superior colliculus",
  "Colliculus superior",
  
  "Substantia nigra",
  "Substantia nigra",
  
  "Red nucleus",
  "Nucleus ruber",
  
  # ===================== CEREBELLUM =====================
  
  "Cerebellar cortex",
  "Cerebellar cortex",
  
  "Vermis",
  "Vermis",
  
  "Dentate nucleus",
  "Nucleus dentatus cerebelli",
  
  # ===================== WHITE MATTER =====================
  
  "Centrum semiovale",
  "Centrum semiovale",
  
  "Internal capsule",
  "Capsula interna"
)

# ---- Apply rules (first match wins) ----
obs_neuron <- obs_neuron %>%
  mutate(anatomy_group = NA_character_)

for (i in seq_len(nrow(anatomy_rules))) {
  obs_neuron$anatomy_group[
    is.na(obs_neuron$anatomy_group) &
      str_detect(obs_neuron$dissection,
                 regex(anatomy_rules$pattern[i], ignore_case = TRUE))
  ] <- anatomy_rules$anatomy_group[i]
}

obs_neuron$anatomy_group[is.na(obs_neuron$anatomy_group)] <- "Unmapped"

# inspect the distribution of anatomy groups
unmapped_df <- obs_neuron %>%
  filter(anatomy_group == "Unmapped") %>%
  select(
    dissection,
    ROIGroupFine,
    ROIGroupCoarse,
    roi
  ) %>%
  distinct()

unmapped_df

###################################################################
# Cell type proportion calculations (ALL major neuronal categories)
###################################################################
# ---- Define major neuronal categories ----
# These are mutually exclusive and exhaustive for neurons
obs_neuron <- obs_neuron %>%
  mutate(
    cell_category = case_when(
      # ---- Excitatory projection neurons ----
      supercluster_term %in% c(
        "Amygdala excitatory",
        "Deep-layer corticothalamic and 6b",
        "Deep-layer intratelencephalic",
        "Deep-layer near-projecting",
        "Upper-layer intratelencephalic",
        "Thalamic excitatory",
        "Hippocampal CA1-3",
        "Hippocampal CA4",
        "Hippocampal dentate gyrus",
        "Mammillary body",
        "Lower rhombic lip",
        "Upper rhombic lip"
      ) ~ "Excitatory_projection",
      # ---- Inhibitory interneurons ----
      supercluster_term %in% c(
        "CGE interneuron",
        "MGE interneuron",
        "LAMP5-LHX6 and Chandelier",
        "Cerebellar inhibitory",
        "Midbrain-derived inhibitory"
      ) ~ "Inhibitory_interneuron",
      # ---- Inhibitory projection neurons (basal ganglia principal cells) ----
      supercluster_term %in% c(
        "Medium spiny neuron",
        "Eccentric medium spiny neuron"
      ) ~ "Inhibitory_projection_MSN",
      # ---- Other projection neurons (rare / region-specific) ----
      supercluster_term %in% c(
        "Deep-layer corticothalamic and 6b"
      ) ~ "Other_projection",
      # ---- QC / unresolved neuronal clusters ----
      supercluster_term %in% c(
        "Splatter",
        "Miscellaneous"
      ) ~ "Unresolved_or_QC",
      # ---- Fallback (should be rare) ----
      TRUE ~ "Other_neuron"
    )
  )
# ---- Sanity check: every neuron assigned ----
table(obs_neuron$cell_category)
# ---- Region × cell-category counts and proportions ----

celltype_table_long <- obs_neuron %>%
  filter(anatomy_group != "", anatomy_group != "Unmapped") %>%
  count(anatomy_group, cell_category, name = "n_cells") %>%
  group_by(anatomy_group) %>%
  mutate(
    n_all_neurons = sum(n_cells),
    p_cells = n_cells / n_all_neurons
  ) %>%
  ungroup()

celltype_proportion_table <- celltype_table_long %>%
  # Pivot counts
  select(anatomy_group, cell_category, n_cells) %>%
  tidyr::pivot_wider(
    names_from  = cell_category,
    values_from = n_cells,
    values_fill = 0,
    names_prefix = "n_"
  ) %>%
  # Join proportions
  left_join(
    celltype_table_long %>%
      select(anatomy_group, cell_category, p_cells) %>%
      tidyr::pivot_wider(
        names_from  = cell_category,
        values_from = p_cells,
        values_fill = 0,
        names_prefix = "p_"
      ),
    by = "anatomy_group"
  ) %>%
  # Add total neuron count only (no forced proportions)
  left_join(
    celltype_table_long %>%
      distinct(anatomy_group, n_all_neurons),
    by = "anatomy_group"
  )

# ---- QC: do neuron proportions sum to ~1 by region? ----

qc_prop_sum <- celltype_proportion_table %>%
  mutate(
    p_sum = rowSums(
      dplyr::select(., dplyr::starts_with("p_")),
      na.rm = TRUE
    )
  ) %>%
  select(anatomy_group, p_sum) %>%
  arrange(desc(abs(p_sum - 1)))

qc_prop_sum

######################################################
# Read and Bind rCMRGlc values from Heiss et al. 2004
######################################################

# Read in rCMRGlc values from Excel file (from Sup Table 1 of Study 1)
xlsx_path <- "/Users/crossmodal/Library/CloudStorage/Dropbox/COLLABORATIVE/Do expensive brain regions increase less in humans/MS Is human brain organization economical/Supplementary Figs and Tables/Sup Table 1 Regional cerebral metabolic rates and data sources using in Study 1.xlsx"
sheet_name <- 1

# Skip the title row(s) so the next row becomes the header
raw <- read_excel(xlsx_path, sheet = sheet_name, skip = 1, col_names = TRUE)

# Select and clean the relevant columns (3rd = anatomy term, 4th = rCMRGlc value)
rcmr <- raw %>%
  select(3, 4) %>%                                 # 3rd = anatomy term, 4th = rCMRGlc value
  setNames(c("rcmr_term", "rcmr_value")) %>%        # rename by position
  mutate(
    rcmr_term  = str_trim(as.character(rcmr_term)),
    rcmr_value = suppressWarnings(as.numeric(rcmr_value))
  ) %>%
  filter(!is.na(rcmr_term), rcmr_term != "", !is.na(rcmr_value))

head(rcmr)

# join to rCMRGlc table
analysis_df <- celltype_proportion_table %>%
  inner_join(rcmr, by = c("anatomy_group" = "rcmr_term"))

###########
# Analyses
###########

#############################
## AUTOMATIC PREDICTOR SETUP
#############################

exclude <- c("p_all_neurons", "p_all_nonneurons", "p_sum")

# Detect all p_ columns and optionally exclude some
all_p <- names(analysis_df)[grepl("^p_", names(analysis_df))]
predictors <- setdiff(all_p, exclude)

cat("Using predictors:\n")
print(predictors)

# Check that predictors match all p_ columns except those explicitly excluded

cat("\nCheck: predictors match all p_ columns except excluded:\n")

if (setequal(predictors, setdiff(all_p, exclude))) {
  cat(" ✓ OK\n")
} else {
  cat(" ⚠ Mismatch detected\n")
}

###########################
## MULTIVARIATE REGRESSION
###########################

fit <- lm(
  reformulate(predictors, response = "rcmr_value"),
  data = analysis_df
)
summary(fit)

###########################
## CORRELATIONS (LOOP)
###########################

cor_results <- do.call(
  rbind,
  lapply(predictors, function(v) {
    ct <- cor.test(
      analysis_df$rcmr_value,
      analysis_df[[v]],
      method = "spearman",
      exact = FALSE
    )
    data.frame(
      variable = v,
      r        = unname(ct$estimate),
      p_value  = ct$p.value,
      n        = sum(complete.cases(analysis_df$rcmr_value, analysis_df[[v]]))
    )
  })
) %>%
  dplyr::arrange(p_value) %>%
  dplyr::mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

cor_results

###########################
## FACETED SCATTERPLOT
###########################

plot_df <- analysis_df %>%
  dplyr::select(anatomy_group, rcmr_value, all_of(predictors)) %>%
  tidyr::pivot_longer(
    cols = all_of(predictors),
    names_to = "predictor",
    values_to = "prop"
  )

# Build high-contrast palette
regions <- sort(unique(plot_df$anatomy_group))
base_cols <- c(
  "#000000", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)
pal <- grDevices::colorRampPalette(base_cols)(length(regions))
names(pal) <- regions

ggplot(plot_df, aes(x = prop, y = rcmr_value, color = anatomy_group)) +
  geom_point(size = 2.8, alpha = 0.85) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "steelblue") +
  stat_poly_eq(
    aes(label = paste(after_stat(rr.label), after_stat(p.value.label), sep = "*\", \"*")),
    formula = y ~ x,
    parse = TRUE,
    label.x = "right",
    label.y = "top",
    size = 4,
    color = "black"
  ) +
  facet_wrap(~ predictor, scales = "free_x") +
  scale_color_manual(values = pal) +
  labs(
    x = "Cell-type proportion",
    y = "rCMRGlc (µmol/100 g/min.)",
    color = "Region"
  ) +
  theme_classic(base_size = 14)