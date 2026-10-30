options(device = "png")
library(tidyverse)
source("helpers/plot_settings.R")

## ============================================================
## Study s1d_2 — Regional composition of the five neuroglial proteins
##
## Replicates the design of s1a_2 (per-region pie charts of cell-type
## proportions) with the five Sjölin et al. 2022 proteins (GFAP, MBP,
## NFL, tau, UCHL1) in place of cell types. Shares are of the summed
## mass of the five proteins per g wet tissue, from the regional
## geometric means written by s1d_1.
##
## Two variants per region, as s1a_2 has Counts and Density pies:
##   AllFive  - all five proteins (MBP dominates by mass everywhere)
##   ExclMBP  - the four non-myelin proteins, so their relative
##              composition is visible
## All mapped analysis regions are drawn, including the three without a
## Heiss rate (pons, medulla, cervical spinal cord); the combined facet
## figure orders regions by rCMRGlc with the unrated regions last.
## ============================================================

# -----------------------------
# Paths
# -----------------------------
input_file <- "data_intermediate/s1d_neuroglial_protein_regional_summary.csv"
outdir <- "figs/s1d/protein_pies"

if (!file.exists(input_file)) {
  stop("Run s1d_1_neuroglial_proteins_metabolic_rate.R first; missing ",
       input_file, call. = FALSE)
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Load data
# -----------------------------
d <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE,
              na.strings = c("", "NA"))

protein_columns <- c(
  GFAP  = "gfap_ug_per_g",
  MBP   = "mbp_ug_per_g",
  NFL   = "nfl_ug_per_g",
  Tau   = "tau_ug_per_g",
  UCHL1 = "uchl1_ug_per_g"
)
missing_columns <- setdiff(c("anatomy_id", "anatomy_group", "rcmr_value",
                             protein_columns), names(d))
if (length(missing_columns)) {
  stop(input_file, " is missing column(s): ",
       paste(missing_columns, collapse = ", "), call. = FALSE)
}

# Stable protein colours (Okabe-Ito), shared by every pie.
protein_palette <- c(
  GFAP  = "#E69F00",
  MBP   = "#56B4E9",
  NFL   = "#009E73",
  Tau   = "#D55E00",
  UCHL1 = "#CC79A7"
)

# Region display order: by rCMRGlc, unrated regions last by registry order.
registry <- read.csv("metadata/anatomy/anatomy_registry.csv",
                     stringsAsFactors = FALSE, check.names = FALSE)
d <- d %>%
  left_join(registry %>% select(anatomy_id, display_order), by = "anatomy_id") %>%
  arrange(is.na(rcmr_value), rcmr_value, display_order) %>%
  mutate(
    region_label = ifelse(
      is.na(rcmr_value),
      paste0(anatomy_group, "\n(no rCMRGlc)"),
      paste0(anatomy_group, "\n", formatC(rcmr_value, format = "f", digits = 1),
             " µmol/100 g/min")
    ),
    region_label = factor(region_label, levels = region_label)
  )

if (anyDuplicated(d$anatomy_id)) {
  stop("Regional summary has more than one row per anatomy_id.", call. = FALSE)
}

# -----------------------------
# Long format: one row per region x protein, with shares
# -----------------------------
long <- d %>%
  select(anatomy_id, anatomy_group, region_label, rcmr_value,
         all_of(unname(protein_columns))) %>%
  pivot_longer(cols = all_of(unname(protein_columns)),
               names_to = "column", values_to = "ug_per_g") %>%
  mutate(protein = names(protein_columns)[match(column, protein_columns)],
         protein = factor(protein, levels = names(protein_palette))) %>%
  select(-column)

shares <- bind_rows(
  long %>%
    group_by(anatomy_id) %>%
    mutate(variant = "AllFive", share = ug_per_g / sum(ug_per_g)) %>%
    ungroup(),
  long %>%
    filter(protein != "MBP") %>%
    group_by(anatomy_id) %>%
    mutate(variant = "ExclMBP", share = ug_per_g / sum(ug_per_g)) %>%
    ungroup()
)

shares %>%
  select(variant, anatomy_id, anatomy_group, rcmr_value, protein,
         ug_per_g, share) %>%
  arrange(variant, rcmr_value, anatomy_group, protein) %>%
  write.csv("tables/s1d/neuroglial_protein_composition_shares.csv",
            row.names = FALSE, na = "")

# -----------------------------
# Helper to save pie charts (as in s1a_2)
# -----------------------------
save_pie <- function(values, labels, title, filename, outdir) {

  filename <- basename(filename)

  df <- data.frame(
    label = labels,
    value = as.numeric(values),
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(value), value > 0) %>%
    mutate(label = factor(label, levels = names(protein_palette)))

  if (nrow(df) == 0) {
    warning("Skipping empty plot: ", title)
    return(invisible(NULL))
  }

  p <- ggplot(df, aes(x = "", y = value, fill = label)) +
    geom_col(color = "white", width = 1) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = protein_palette, drop = TRUE) +
    theme_void() +
    labs(
      title = title,
      fill = "Protein"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5)
    )

  outfile <- file.path(outdir, filename)

  ggsave(
    filename = outfile,
    plot = p,
    width = 5,
    height = 5,
    dpi = 300,
    create.dir = TRUE
  )

  message("Saved: ", outfile)

  invisible(outfile)
}

# -----------------------------
# Generate pies (one per region and variant)
# -----------------------------
for (i in seq_len(nrow(d))) {

  row <- d[i, ]
  r <- row$anatomy_group
  file_stem <- gsub("[^A-Za-z0-9]+", "_", r)

  save_pie(
    values = unlist(row[protein_columns]),
    labels = names(protein_columns),
    title = paste(r, "(five proteins, µg/g share)"),
    filename = paste0("Pie_AllFive_", file_stem, ".png"),
    outdir = outdir
  )

  non_mbp <- protein_columns[names(protein_columns) != "MBP"]
  save_pie(
    values = unlist(row[non_mbp]),
    labels = names(non_mbp),
    title = paste(r, "(excluding MBP, µg/g share)"),
    filename = paste0("Pie_ExclMBP_", file_stem, ".png"),
    outdir = outdir
  )
}

# -----------------------------
# Combined facet figures: all regions on one page, ordered by rCMRGlc
# -----------------------------
make_combined <- function(variant_id, title) {
  plot_df <- shares %>%
    filter(variant == variant_id) %>%
    group_by(region_label) %>%
    arrange(desc(protein), .by_group = TRUE) %>%
    mutate(label_pos = cumsum(share) - share / 2,
           pct_label = ifelse(share >= 0.06,
                              sprintf("%.0f%%", 100 * share), "")) %>%
    ungroup()

  ggplot(plot_df, aes(x = "", y = share, fill = protein)) +
    geom_col(color = "white", width = 1) +
    geom_text(aes(y = label_pos, label = pct_label), size = 2.6) +
    coord_polar(theta = "y") +
    facet_wrap(~ region_label, ncol = 5) +
    scale_fill_manual(values = protein_palette, drop = TRUE) +
    theme_void(base_size = 11) +
    labs(
      title = title,
      subtitle = "Share of summed protein mass per g wet tissue (regional geometric means, Sjölin et al. 2022); regions ordered by rCMRGlc (Heiss et al. 2004)",
      fill = "Protein"
    ) +
    theme(
      strip.text = element_text(size = 8),
      plot.subtitle = element_text(size = 8, color = "grey40"),
      legend.position = "right"
    )
}

p_all <- make_combined("AllFive",
                       "Regional composition of the five neuroglial proteins")
p_excl <- make_combined("ExclMBP",
                        "Regional composition of GFAP, NFL, tau and UCHL1 (MBP excluded)")

for (nm in c("png", "pdf")) {
  ggsave(paste0("figs/s1d/protein_composition_pies_all_five.", nm),
         plot = p_all, width = 11, height = 8, units = "in", dpi = 300,
         bg = "white", device = if (nm == "pdf") cairo_pdf else nm)
  ggsave(paste0("figs/s1d/protein_composition_pies_excl_mbp.", nm),
         plot = p_excl, width = 11, height = 8, units = "in", dpi = 300,
         bg = "white", device = if (nm == "pdf") cairo_pdf else nm)
}

message("Study s1d_2 done. Pies in ", outdir, "/, combined figures in figs/s1d/.")
