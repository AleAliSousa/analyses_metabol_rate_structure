# ============================================================
# Heiss <-> Stephan correspondence + Table 1  —  PIPELINE v2
# Copy of 0_Heiss_Stephan_and_table1_30052026.R, feeding the v2 ("replace
# neocortex grey") analysis s3v2_predicValuesPGLS_05062026.R.
#   * writes data_intermediate/Heiss_Stephan_data_v2.csv (NOT the v1 table)
#   * contains a PHASE 3.2 STUB (see section 5b) that will reweight the
#     Neocortex grey rCMRGlc to the de-overlapped remainder; pass-through for now.
# ============================================================

setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))  # repo root (portable; replaces hardcoded path -- see R/project_root.R)

library(tidyverse)
library(writexl)
library(openxlsx)

heiss <- read.csv("data_raw/Heiss_etal_2004_TABLE1.csv", stringsAsFactors = FALSE, check.names = FALSE)

# -----------------------------
# 1. Helper functions
# -----------------------------

get_mean <- function(region_name) {
  heiss$`Both hemispheres Mean`[heiss$Region == region_name]
}

weighted_mean_available <- function(values, weights) {
  ok <- !is.na(values) & !is.na(weights)
  sum(values[ok] * weights[ok]) / sum(weights[ok])
}

# -----------------------------
# 2. Add equivalence columns
# -----------------------------

heiss$term_1 <- NA
heiss$term_2 <- NA
heiss$term_3 <- heiss$Region

heiss$volume_term <- NA
heiss$volume_source <- NA
heiss$calculation <- NA
heiss$is_calculated_row <- FALSE

# Direct equivalences to volume terms
heiss$volume_term[heiss$Region == "Cerebral cortex (global average)"] <- "Neocortex grey"
heiss$volume_source[heiss$Region == "Cerebral cortex (global average)"] <- 
  "Frahm et al. 1982; Zilles and Rehkamper 1988"

heiss$volume_term[heiss$Region == "Occipital lobe"] <- "Area striata grey"
heiss$volume_source[heiss$Region == "Occipital lobe"] <- 
  "de Sousa et al. 2010; Frahm et al. 1984"

heiss$volume_term[heiss$Region == "Insular lobe"] <- "Insular cortex (grey)"
heiss$volume_source[heiss$Region == "Insular lobe"] <- 
  "Bauernfeind et al. 2013"

heiss$volume_term[heiss$Region == "Hippocampus"] <- "Hippocampus"
heiss$volume_source[heiss$Region == "Hippocampus"] <- 
  "Stephan et al. 1981; Zilles and Rehkamper 1988"

heiss$volume_term[heiss$Region == "Pallidum"] <- "Pallidum"
heiss$volume_source[heiss$Region == "Pallidum"] <- 
  "Stephan et al. 1981; Zilles and Rehkamper 1988"

heiss$volume_term[heiss$Region == "Corpus geniculatum laterale"] <- 
  "Corpus geniculatum laterale"
heiss$volume_source[heiss$Region == "Corpus geniculatum laterale"] <- 
  "de Sousa et al. 2010; Stephan et al. 1984"

heiss$volume_term[heiss$Region == "Nucleus subthalamicus"] <- 
  "Nucleus subthalamicus Luysi"
heiss$volume_source[heiss$Region == "Nucleus subthalamicus"] <- 
  "Stephan et al. 1981"

heiss$volume_term[heiss$Region == "Corpus amygdaloideum"] <- "Amygdala"
heiss$volume_source[heiss$Region == "Corpus amygdaloideum"] <- 
  "Stephan et al. 1981; Zilles and Rehkamper 1988"

heiss$volume_term[heiss$Region == "Centrum semiovale"] <- "Neocortex white"
heiss$volume_source[heiss$Region == "Centrum semiovale"] <- 
  "Frahm et al. 1982; Zilles and Rehkamper 1988"

heiss$volume_term[heiss$Region == "Capsula interna"] <- "Capsula interna"
heiss$volume_source[heiss$Region == "Capsula interna"] <- 
  "Stephan et al. 1981"

# -----------------------------
# 3. Calculate composite rows
# -----------------------------
# These are the weighted averages described in the Supplementary Table 1 note.
# For Cerebellum and Brain stem nuclei, "Other" refers to structures that are unknown and/or not present
# in Heiss. Therefore this script calculates the weighted
# average using the available Heiss regions and renormalizes the available weights.

get_mean <- function(region) {
  x <- heiss$`Both hemispheres Mean`[heiss$Region == region]
  if (length(x) == 0) return(NA_real_)
  x
}

weighted_mean_available <- function(values, weights) {
  ok <- !is.na(values)
  sum(values[ok] * weights[ok]) / sum(weights[ok])
}

calc_mean <- function(weights) {
  values <- sapply(names(weights), get_mean)
  weighted_mean_available(values, weights)
}

calc_text <- function(weights, source) {
  paste0(
    paste0(weights, "*", names(weights), collapse = " + "),
    " (", source, ")"
  )
}

striatum_weights <- c(
  "Caudatum" = 0.45,
  "Putamen" = 0.49,
  "Nucleus accumbens" = 0.07
)

cerebellum_weights <- c(
  "Cerebellar cortex" = 0.56, # Sullivan et al 2000 (as proportion cerebellum)
  "Vermis" = 0.08, # MacLeod et al. 2003
  "Nucleus dentatus cerebelli" = 0.01, # Matano 2001
  "Other" = 0.45
)

brainstem_weights <- c(
  "Colliculus inferior" = 0.08, # Garcia-Gomar et al. 2019
  "Colliculus superior" = 0.11, # Garcia-Gomar et al. 2019
  "Substantia nigra" = 0.53, # Eapen et al. 2011
  "Nucleus ruber" = 0.16, # Eapen et al. 2011
  "Other" = 0.12
)

striatum_mean <- calc_mean(striatum_weights)
cerebellum_mean <- calc_mean(cerebellum_weights)
brainstem_mean <- calc_mean(brainstem_weights)

striatum_calculation <- calc_text(
  striatum_weights,
  "de Jong et al. 2012"
)

cerebellum_calculation <- calc_text(
  cerebellum_weights,
  "Sullivan et al. 2000, MacLeod et al. 2003, Matano 2001"
)

brainstem_calculation <- calc_text(
  brainstem_weights,
  "Eapen et al. 2011, Garcia-Gomar et al. 2019"
)

# -----------------------------
# 4. Add calculated composite rows
# -----------------------------

striatum_row <- data.frame(
  category = "Forebrain nuclei",
  Region = "Striatum",
  `Both hemispheres Mean` = striatum_mean,
  `Both hemispheres SD` = NA,
  `Left minus right hemisphere Difference` = NA,
  `Left minus right hemisphere SD` = NA,
  P = NA,
  `Data of Heiss et al. (21)` = NA,
  term_1 = NA,
  term_2 = "Striatum1",
  term_3 = NA,
  volume_term = "Striatum",
  volume_source = "Stephan et al. 1981; Zilles and Rehkamper 1988",
  calculation = striatum_calculation,
  is_calculated_row = TRUE,
  check.names = FALSE
)

cerebellum_row <- data.frame(
  category = "Cerebellum",
  Region = "Cerebellum",
  `Both hemispheres Mean` = cerebellum_mean,
  `Both hemispheres SD` = NA,
  `Left minus right hemisphere Difference` = NA,
  `Left minus right hemisphere SD` = NA,
  P = NA,
  `Data of Heiss et al. (21)` = NA,
  term_1 = "Cerebellum1",
  term_2 = NA,
  term_3 = NA,
  volume_term = "Cerebellum",
  volume_source = "Stephan et al. 1981; Zilles and Rehkamper 1988",
  calculation = cerebellum_calculation,
  is_calculated_row = TRUE,
  check.names = FALSE
)

brainstem_row <- data.frame(
  category = "Brain stem nuclei",
  Region = "Brain stem nuclei",
  `Both hemispheres Mean` = brainstem_mean,
  `Both hemispheres SD` = NA,
  `Left minus right hemisphere Difference` = NA,
  `Left minus right hemisphere SD` = NA,
  P = NA,
  `Data of Heiss et al. (21)` = NA,
  term_1 = "Brain stem nuclei1",
  term_2 = NA,
  term_3 = NA,
  volume_term = "Mesencephalon",
  volume_source = "Stephan et al. 1981; Zilles and Rehkamper 1988",
  calculation = brainstem_calculation,
  is_calculated_row = TRUE,
  check.names = FALSE
)

data_table <- rbind(
  heiss,
  striatum_row,
  cerebellum_row,
  brainstem_row
)

# -----------------------------
# 5. Order table
# -----------------------------

# Determine order from data_table itself
data_table$original_row_order <- seq_len(nrow(data_table))

# Use the original non-calculated rows to define the group order
category_order <- unique(data_table$category[!data_table$is_calculated_row])

data_table$category_order <- match(data_table$category, category_order)

# Within each category:
#   calculated rows first
#   then original rows in their original order
order_table <- data_table[order(
  data_table$category_order,
  !data_table$is_calculated_row,
  data_table$original_row_order
), ]

# Rename columns and omit helper columns
cols <- c(
  category = "category",
  term_1 = "term_1",
  term_2 = "term_2",
  term_3 = "term_3",
  Heiss_region = "Region",
  rCMRGlc_mean_both_hemispheres = "Both hemispheres Mean",
  volume_term = "volume_term",
  volume_source = "volume_source",
  calculation = "calculation"
)

order_table <- setNames(order_table[unname(cols)], names(cols))

# ============================================================
# PHASE 3.2 STUB — reweight Neocortex grey rCMRGlc for the de-overlapped remainder
# ------------------------------------------------------------
# In pipeline v2, the modelled "Neocortex grey" is the REMAINDER after removing
# V1/occipital and insula. Its rCMRGlc should therefore be the cortex-wide value
# with the V1 and insula glucose contributions removed by volume weighting —
# exactly the approach already used above for Cerebellum (cortex + vermis) and
# Brain stem nuclei.
#
#   rCMR(neocortex global) is a volume-weighted mean that INCLUDES occipital + insula:
#     rCMR_global * V_neo = rCMR_rem * V_rem + rCMR_occ * V_V1 + rCMR_ins * V_ins
#   => rCMR_rem = (rCMR_global * V_neo - rCMR_occ * V_V1 - rCMR_ins * V_ins) / V_rem
#            = (rCMR_global - rCMR_occ * f_V1 - rCMR_ins * f_ins) / (1 - f_V1 - f_ins)
#     where f_V1 = V_V1 / V_neo, f_ins = V_ins / V_neo.
#
#   Inputs (TODO Phase 3.2):
#     rCMR_global = Heiss "Cerebral cortex (global average)" mean   (in this table)
#     rCMR_occ    = Heiss "Occipital lobe" mean                     (in this table)
#     rCMR_ins    = Heiss "Insular lobe" mean                       (in this table)
#     f_V1, f_ins = volume fractions of V1 (ASG_Sousa) and insula (Total_insula_volume_L)
#                   within neocortex grey (NeoG_Frahm), from data_raw/Stephan_primates.csv
#                   or literature (mirror the cerebellum literature-fraction approach).
#
# STUB: pass-through for now, so the v2 table is structurally identical to v1 and
# the v2 pipeline runs end-to-end. Phase 3.2 fills in the function body.
# ============================================================
reweight_neocortex_grey <- function(rCMR_global,
                                    rCMR_occ = NA, rCMR_ins = NA,
                                    f_V1 = NA, f_ins = NA) {
  # TODO Phase 3.2: return (rCMR_global - rCMR_occ*f_V1 - rCMR_ins*f_ins) / (1 - f_V1 - f_ins)
  rCMR_global   # placeholder: no reweighting yet
}

i_neo <- which(order_table$volume_term == "Neocortex grey")
if (length(i_neo) == 1) {
  order_table$rCMRGlc_mean_both_hemispheres[i_neo] <-
    reweight_neocortex_grey(order_table$rCMRGlc_mean_both_hemispheres[i_neo])
}

# -----------------------------
# 6. Save intermediate final table  (v2 — see Phase 3.2 stub above)
# -----------------------------

write.csv(
  order_table,
  "data_intermediate/Heiss_Stephan_data_v2.csv",
  row.names = FALSE,
  na = ""
)

# -----------------------------
# 7. Publication formatting
# -----------------------------

pub <- order_table[, c(
  "term_1",
  "term_2",
  "term_3",
  "rCMRGlc_mean_both_hemispheres",
  "volume_term",
  "volume_source",
  "Heiss_region"
)]

i <- pub$Heiss_region == "Cerebral cortex (global average)"
pub$term_1[i] <- pub$term_3[i]
pub$term_3[i] <- NA

section <- pub[NA_integer_, ]

out <- pub[0, ]

for (i in seq_len(nrow(pub))) {
  if (pub$Heiss_region[i] == "Striatum") {
    section[1, ] <- NA
    section$term_1 <- "Forebrain nuclei"
    out <- rbind(out, section)
  }
  
  if (pub$Heiss_region[i] == "Centrum semiovale") {
    section[1, ] <- NA
    section$term_1 <- "White matter"
    out <- rbind(out, section)
  }
  
  out <- rbind(out, pub[i, ])
}

# Drop helper column
out <- out[, c("term_1", "term_2", "term_3", "rCMRGlc_mean_both_hemispheres", "volume_term", "volume_source")]

# Make workbook
wb <- createWorkbook()
addWorksheet(wb, "Sheet1")

# Title
writeData(
  wb, "Sheet1",
  "Table 1. Regional cerebral metabolic rates, corresponding volume terms and data sources.",
  startRow = 1,
  startCol = 1,
  colNames = FALSE
)

# Headers
writeData(
  wb, "Sheet1",
  matrix(c(
    "rCMRGlc term (Heiss et al. 2004)", NA, NA,
    "rCMRGlc (µmol/100 g/min.)",
    "Volume term",
    "Region volume source2"
  ), nrow = 1),
  startRow = 3,
  startCol = 1,
  colNames = FALSE
)

writeData(
  wb, "Sheet1",
  matrix(c(NA, NA, NA, "mean (both hemispheres)", NA, NA), nrow = 1),
  startRow = 4,
  startCol = 1,
  colNames = FALSE
)

# Body
writeData(
  wb, "Sheet1",
  out,
  startRow = 5,
  startCol = 1,
  colNames = FALSE,
  na.string = ""
)

# Note
note_row <- 5 + nrow(out) + 1

note_text <- paste0(
  "Note. 1. rCMRGlc calculated in the current study as weighted averages based on the relative size of the subregion (sources in parentheses). ",
  "Striatum = ", striatum_calculation, "; ",
  "Cerebellum = ", cerebellum_calculation, "; ",
  "Brain stem nuclei = ", brainstem_calculation, ".\n",
  "2. All region volume source studies included phylogenetically diverse primate brain specimens from the Duesseldorf Stephan and Zilles collections, ",
  "and obtained by overlapping research teams."
)

writeData(wb, "Sheet1", note_text, startRow = note_row, startCol = 1, colNames = FALSE)
mergeCells(wb, "Sheet1", cols = 1:6, rows = note_row)

# Formatting
header_style <- createStyle(
  border = c("top", "bottom"),
  borderStyle = "thin",
  halign = "center",
  valign = "center"
)

number_style <- createStyle(numFmt = "0.0")
note_style <- createStyle(wrapText = TRUE, valign = "top")

addStyle(wb, "Sheet1", header_style, rows = 3:4, cols = 1:6, gridExpand = TRUE)
addStyle(wb, "Sheet1", number_style, rows = 5:(4 + nrow(out)), cols = 4)
addStyle(wb, "Sheet1", note_style, rows = note_row, cols = 1)
 
setColWidths(wb, "Sheet1", cols = 1, widths = 22)
setColWidths(wb, "Sheet1", cols = 2, widths = 14)
setColWidths(wb, "Sheet1", cols = 3, widths = 28)
setColWidths(wb, "Sheet1", cols = 4, widths = 22)
setColWidths(wb, "Sheet1", cols = 5, widths = 26)
setColWidths(wb, "Sheet1", cols = 6, widths = 45)

setRowHeights(wb, "Sheet1", rows = note_row, heights = 90)

saveWorkbook(wb, "tables/Table 1 v2 Regional cerebral metabolic rates corresponding volume terms and data sources.xlsx", overwrite = TRUE)