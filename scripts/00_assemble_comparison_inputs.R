# =====================================================================
# 00_assemble_comparison_inputs.R
#
# Build the Study 1a and Study 2 comparison inputs directly from the
# per-paper PUBLIC source files that metadata/evo_m1_input_manifest.csv
# syncs into data_raw/. This REPLACES the previously hand-assembled
# data_raw/stereology_comparison.csv and data_raw/stress_volume_comparison.csv
# with a reproducible build, and runs BEFORE 00_prepare_heiss_rates.R.
#
# CORRECTED BUILD (see RECONSTRUCTION_FINDINGS.txt). Differences vs the old
# hand-made files, all deliberate:
#   [1] Thalamus norm_english now computed by formula (old file had a 10x
#       decimal error: -0.273 -> -2.727).
#   [2] Total Brain SD kept on the whole-brain basis, consistent with its
#       whole-brain mean (old file halved only the SD: 62.589 -> 125.178).
#   [3] Cerebral-cortex adoptee SD uses independent-variance propagation
#       sqrt(sum(parcel_sd^2)) instead of sum(parcel_sd) (old file inflated
#       it ~7x by summing SDs: 34.77 -> ~5.09 UK, 32.35 -> ~4.69 Romanian).
#   [4] Stereology cell-count SD columns are set NA: they are unused by any
#       Study 1a script and have no source in the public per-paper files.
#   [minor] Nucleus accumbens glia:neuron ratio uses Hanson's published value
#       (2.03) rather than the recomputed 2.04.
#
# All numeric outputs otherwise reproduce the previous files within rounding.
# =====================================================================

# ---- locate repo root (portable: Rscript or interactive) ------------
.script_args <- commandArgs(trailingOnly = FALSE)
.file_arg <- .script_args[grepl("^--file=", .script_args)]
.start_dir <- if (length(.file_arg)) {
  .p <- sub("^--file=", "", .file_arg[1])
  .p <- gsub("~+~", " ", .p, fixed = TRUE)
  dirname(normalizePath(.p, mustWork = FALSE))
} else getwd()
.repo_root <- normalizePath(.start_dir, mustWork = FALSE)
while (!file.exists(file.path(.repo_root, ".git")) &&
       dirname(.repo_root) != .repo_root) {
  .repo_root <- dirname(.repo_root)
}
if (!file.exists(file.path(.repo_root, ".git"))) {
  stop("Could not locate the project root from: ", .start_dir, call. = FALSE)
}
setwd(.repo_root)

rd <- function(p) {
  if (!file.exists(p)) {
    stop("Missing per-paper input (run 00_sync_evo_m1_inputs.R first): ", p,
         call. = FALSE)
  }
  read.csv(p, stringsAsFactors = FALSE, check.names = FALSE,
           na.strings = c("", "NA"))
}

# 95% CI half-width using a t-distribution (matches the original build).
ci_half <- function(sd, n) stats::qt(0.975, n - 1) * sd / sqrt(n)

# =====================================================================
# STUDY 1a: stereology_comparison.csv
# =====================================================================
au <- rd("data_raw/Karlsen_Pakkenberg_2011_authordata_groupmeans.csv")
au <- au[au$Group == "Control", ]
t2 <- rd("data_raw/Karlsen_Pakkenberg_2011_Table2_derived.csv")
t2 <- t2[t2$Group == "Control", ]
mo <- rd("data_raw/Morgan_etal_2014_DataS1_derived.csv")
moC <- mo[mo$Diagnosis == "Control", ]
ha <- rd("data_raw/Hanson_etal_2018_Table2_derived.csv")

stereo_cols <- c(
  "source_region", "Volume", "VolumeSD",
  "Neuron_N", "NeuronSD", "NeurDensity",
  "Glia_N", "GliaSD", "GliaDensity",
  "Astro_N", "AstroSD", "AstroDensity",
  "Oligo_N", "OligoSD", "OligoDensity",
  "Microglia_N", "MicroSD", "MicroDensity",
  "GliaNeurDensityRatio"
)
new_stereo_row <- function() {
  r <- as.list(rep(NA_real_, length(stereo_cols)))
  names(r) <- stereo_cols
  r
}

# KP Control cortex + lobes: cell numbers (authordata) / volume (Table2).
kp_rows <- list(
  "Cerebral cortex (global average)" = c(au = "Neocortex total", t2 = "Total neocortex"),
  "Frontal lobe"  = c(au = "Frontal",   t2 = "Frontal"),
  "Parietal lobe" = c(au = "Parietal",  t2 = "Parietal"),
  "Temporal lobe" = c(au = "Temporal",  t2 = "Temporal"),
  "Occipital lobe"= c(au = "Occipital", t2 = "Occipital")
)
stereo_list <- list()
for (lbl in names(kp_rows)) {
  a <- au[au$Region == kp_rows[[lbl]][["au"]], ]
  v <- t2[t2$Region == kp_rows[[lbl]][["t2"]], ]
  vol <- v$Volume.cm3
  nN  <- a$`Neuron_N.x10.9`        * 1e9
  gN  <- a$`Glia_N.x10.9`          * 1e9
  asN <- a$`Astrocyte_N.x10.9`     * 1e9
  olN <- a$`Oligodendrocyte_N.x10.9` * 1e9
  miN <- a$`Microglia_N.x10.9`     * 1e9
  r <- new_stereo_row()
  r$source_region <- lbl
  r$Volume <- vol
  r$Neuron_N <- nN;     r$NeurDensity  <- nN / vol
  r$Glia_N   <- gN;     r$GliaDensity  <- gN / vol
  r$Astro_N  <- asN;    r$AstroDensity <- asN / vol
  r$Oligo_N  <- olN;    r$OligoDensity <- olN / vol
  r$Microglia_N <- miN; r$MicroDensity <- miN / vol
  r$GliaNeurDensityRatio <- gN / nN
  stereo_list[[lbl]] <- as.data.frame(r, stringsAsFactors = FALSE,
                                      check.names = FALSE)
}

# Nucleus accumbens: Hanson (typically-developing group).
hn <- ha[ha$Region == "Nucleus Accumbens" & ha$Group == "TD", ]
r <- new_stereo_row()
r$source_region <- "Nucleus accumbens"
r$NeurDensity <- hn$`Neuron_density.per.mm3`
r$GliaDensity <- hn$`Glia_density.per.mm3`
r$GliaNeurDensityRatio <- hn$Glia_Neuron_ratio   # Hanson published (2.03)
stereo_list[["Nucleus accumbens"]] <- as.data.frame(r, stringsAsFactors = FALSE,
                                                    check.names = FALSE)

# Corpus amygdaloideum: Morgan Control-group mean x2 (bilateral).
volA <- mean(moC$`Amyg Volume`, na.rm = TRUE) * 2 / 1000  # mm3 -> cm3
nN  <- mean(moC$Whole_Num_Neuron, na.rm = TRUE) * 2
gN  <- mean(moC$Whole_Num_Glia,   na.rm = TRUE) * 2
asN <- mean(moC$Whole_Num_Astro,  na.rm = TRUE) * 2
olN <- mean(moC$Whole_Num_Olig,   na.rm = TRUE) * 2
miN <- mean(moC$Whole_Num_Micro,  na.rm = TRUE) * 2
r <- new_stereo_row()
r$source_region <- "Corpus amygdaloideum"
r$Volume <- volA
r$Neuron_N <- nN;     r$NeurDensity  <- nN / volA
r$Glia_N   <- gN;     r$GliaDensity  <- gN / volA
r$Astro_N  <- asN;    r$AstroDensity <- asN / volA
r$Oligo_N  <- olN;    r$OligoDensity <- olN / volA
r$Microglia_N <- miN; r$MicroDensity <- miN / volA
r$GliaNeurDensityRatio <- gN / nN
stereo_list[["Corpus amygdaloideum"]] <- as.data.frame(r, stringsAsFactors = FALSE,
                                                       check.names = FALSE)

stereology <- do.call(rbind, stereo_list)[, stereo_cols]
dir.create("data_raw", showWarnings = FALSE)
write.csv(stereology, "data_raw/stereology_comparison.csv",
          row.names = FALSE, na = "")
message("Assembled data_raw/stereology_comparison.csv (",
        nrow(stereology), " rows) from KP/Hanson/Morgan public sources.")

# =====================================================================
# STUDY 2: stress_volume_comparison.csv
# =====================================================================
wa <- rd("data_raw/Walhovd_etal_2011_derived.csv")
sc <- rd("data_raw/Mackes_etal_2020_subcortical_derived.csv")
er <- rd("data_raw/Mackes_etal_2020_ERABIS_derived.csv")
er_parcels <- er[!is.na(er$Hemisphere) & er$Hemisphere %in% c("L", "R"), ]

N_NORM <- 262; N_UK <- 21; N_ROM <- 67

# region label -> (Walhovd structure, Mackes subcortical region, whole_brain?)
stress_defs <- list(
  "Cerebral cortex" = list(wa = "Cerebral Cortex (GM)", sc = NA,               kind = "cortex"),
  "Thalamus"        = list(wa = "Thalamus",             sc = "Thalamus",       kind = "sub"),
  "Caudate"         = list(wa = "Caudate",              sc = "Caudate",        kind = "sub"),
  "Putamen"         = list(wa = "Putamen",              sc = "Putamen",        kind = "sub"),
  "Pallidum"        = list(wa = "Pallidum",             sc = "Pallidum",       kind = "sub"),
  "Hippocampus"     = list(wa = "Hippocampus",          sc = "Hippocampus",    kind = "sub"),
  "Amygdala"        = list(wa = "Amygdala",             sc = "Amygdala",       kind = "sub"),
  "Accumbens"       = list(wa = "Accumbens",            sc = "Nucleus Accumbens", kind = "sub"),
  "Total Brain"     = list(wa = "Total volume",         sc = NA,               kind = "brain")
)

stress_rows <- list()
for (lbl in names(stress_defs)) {
  d <- stress_defs[[lbl]]
  w <- wa[wa$Structure == d$wa, ]

  # --- normative (Walhovd) ---
  if (d$kind == "brain") {
    norm_mean <- w$Overall_N262_mean.mm3 / 1000        # whole brain, cm3
    norm_sd   <- w$Overall_N262_SD.mm3   / 1000        # [2] consistent whole-brain SD
  } else {
    norm_mean <- w$Overall_N262_mean.mm3 / 2000        # per hemisphere, cm3
    norm_sd   <- w$Overall_N262_SD.mm3   / 2000
  }

  # --- adoptees (Mackes) ---
  if (d$kind == "sub") {
    su <- sc[sc$Region == d$sc & sc$Group == "UK", ]
    sr <- sc[sc$Region == d$sc & sc$Group == "Romanian", ]
    eng_mean <- su$Volume_per_hemisphere.cm3; eng_sd <- su$SD.cm3
    rom_mean <- sr$Volume_per_hemisphere.cm3; rom_sd <- sr$SD.cm3
  } else if (d$kind == "cortex") {
    puk <- er_parcels[er_parcels$Group == "UK", ]
    pro <- er_parcels[er_parcels$Group == "Romanian", ]
    eng_mean <- sum(puk$mean.mm3) / 2000               # per hemisphere, cm3
    rom_mean <- sum(pro$mean.mm3) / 2000
    eng_sd <- sqrt(sum(puk$sd.mm3^2)) / 2000           # [3] independent propagation
    rom_sd <- sqrt(sum(pro$sd.mm3^2)) / 2000
  } else { # brain: Mackes BrainSegVolNotVent (whole)
    bu <- er[er$Group == "UK" & er$Structure == "BrainSegVolNotVent.x", ]
    br <- er[er$Group == "Romanian" & er$Structure == "BrainSegVolNotVent.x", ]
    eng_mean <- bu$mean.mm3 / 1000; eng_sd <- bu$sd.mm3 / 1000
    rom_mean <- br$mean.mm3 / 1000; rom_sd <- br$sd.mm3 / 1000
  }

  # --- % volume change vs normative (not defined for whole brain) ---
  if (d$kind == "brain") {
    norm_all <- NA_real_; norm_english <- NA_real_; norm_romanian <- NA_real_
  } else {
    norm_all      <- (norm_mean - (eng_mean + rom_mean) / 2) / norm_mean * 100
    norm_english  <- (norm_mean - eng_mean) / norm_mean * 100   # [1] correct formula
    norm_romanian <- (norm_mean - rom_mean) / norm_mean * 100
  }

  stress_rows[[lbl]] <- data.frame(
    source_region = lbl,
    norm_all = norm_all, norm_english = norm_english, norm_romanian = norm_romanian,
    english_mean = eng_mean, english_sd = eng_sd,
    english_lci = eng_mean - ci_half(eng_sd, N_UK),
    english_uci = eng_mean + ci_half(eng_sd, N_UK),
    romanian_mean = rom_mean, romanian_sd = rom_sd,
    romanian_lci = rom_mean - ci_half(rom_sd, N_ROM),
    romanian_uci = rom_mean + ci_half(rom_sd, N_ROM),
    norm_mean = norm_mean, norm_sd = norm_sd,
    norm_lci = norm_mean - ci_half(norm_sd, N_NORM),
    norm_uci = norm_mean + ci_half(norm_sd, N_NORM),
    stringsAsFactors = FALSE, check.names = FALSE
  )
}
stress <- do.call(rbind, stress_rows)
write.csv(stress, "data_raw/stress_volume_comparison.csv",
          row.names = FALSE, na = "")
message("Assembled data_raw/stress_volume_comparison.csv (",
        nrow(stress), " rows) from Walhovd/Mackes public sources.")
