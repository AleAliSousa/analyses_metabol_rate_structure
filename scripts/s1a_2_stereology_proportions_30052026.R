options(device = "png")
library(tidyverse)
source("helpers/read_heiss_rates.R")

# -----------------------------
# Paths
# -----------------------------
input_file <- "data_intermediate/s1a_stereology_comparison.csv"
outdir <- "figs/s1a/stereology_pies"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Load data
# -----------------------------
d <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE) %>%
  inner_join(
    as_tibble(read_heiss_rcmr(include_global_average = TRUE)) %>%
      select(anatomy_id, rcmr_term),
    by = "anatomy_id"
  )

# -----------------------------
# Standardize region names
# -----------------------------
d <- d %>%
  mutate(
    Region = case_when(
      rcmr_term == "Cerebral cortex (global average)" ~ "Cerebral Cortex",
      rcmr_term == "Frontal lobe" ~ "Frontal",
      rcmr_term == "Parietal lobe" ~ "Parietal",
      rcmr_term == "Temporal lobe" ~ "Temporal",
      rcmr_term == "Occipital lobe" ~ "Occipital",
      rcmr_term == "Corpus amygdaloideum" ~ "Amygdala",
      TRUE ~ rcmr_term
    )
  )

regions <- c(
  "Cerebral Cortex",
  "Frontal",
  "Parietal",
  "Temporal",
  "Occipital",
  "Amygdala"
)

d <- d %>%
  filter(Region %in% regions)

# -----------------------------
# Helper to save pie charts
# -----------------------------
save_pie <- function(values, labels, title, filename, outdir) {
  
  # Prevent duplicated directory paths if filename accidentally contains a path
  filename <- basename(filename)
  
  df <- data.frame(
    label = labels,
    value = as.numeric(values),
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(value), value > 0)
  
  if (nrow(df) == 0) {
    warning("Skipping empty plot: ", title)
    return(invisible(NULL))
  }
  
  p <- ggplot(df, aes(x = "", y = value, fill = label)) +
    geom_col(color = "white", width = 1) +
    coord_polar(theta = "y") +
    theme_void() +
    labs(
      title = title,
      fill = "Cell type"
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
# Generate pies
# -----------------------------
for (r in regions) {
  
  row <- d %>%
    filter(Region == r)
  
  if (nrow(row) == 0) {
    warning("No data found for region: ", r)
    next
  }
  
  if (nrow(row) > 1) {
    warning("Multiple rows found for region: ", r, ". Using the first row.")
    row <- row[1, ]
  }
  
  save_pie(
    values = c(
      row$Neuron_N,
      row$Astro_N,
      row$Oligo_N,
      row$Microglia_N
    ),
    labels = c(
      "Neuron",
      "Astro",
      "Oligo",
      "Microglia"
    ),
    title = paste(r, "(Cell Number)"),
    filename = paste0("Pie_Counts_", gsub(" ", "_", r), ".png"),
    outdir = outdir
  )
  
  save_pie(
    values = c(
      row$NeurDensity,
      row$AstroDensity,
      row$OligoDensity,
      row$MicroDensity
    ),
    labels = c(
      "Neuron",
      "Astro",
      "Oligo",
      "Microglia"
    ),
    title = paste(r, "(Cell Density)"),
    filename = paste0("Pie_Density_", gsub(" ", "_", r), ".png"),
    outdir = outdir
  )
}
