# ============================================================
# s4_endocranial.R
# Endocast-scaled regional brain energy budgets: Neanderthal (NT),
# early modern human (EH) and modern human (MH).
#
# DATA PROVENANCE & CITATIONS: see ../CITATIONS.md
#   Rates  = Heiss et al. 2004 (J Nucl Med 45:1811), regional rCMRGlc only --
#            umol glucose/100 g/min; NO whole-brain total in that paper.
#   Volumes= Kochiyama et al. 2018 (Sci Rep 8:6296).
#   Whole-brain benchmark (~429 umol/min; ~111 g glucose/day) is Clarke &
#   Sokoloff 1994 (via the Karbowski 2007 compilation) -- external sanity
#   check, NOT Heiss. See ../CITATIONS.md and the summing check in
#   Evo-M1-Trait-Data/__energetics_comparison/heiss_wholebrain_check.R.
#
# QUESTION
#   Do fossil hominins have similar regional brain energy budgets to
#   modern humans? Kochiyama et al. (2018) give, for each parcellated
#   brain region, the volume of NT and EH RELATIVE to modern human
#   (MH = 1). Heiss et al. (2004) give the modern-human regional glucose
#   metabolic rate (rCMRGlc). Combining them:
#
#       budget_region(species) = rCMRGlc_region            # per-mass rate, MH
#                              x  MH_mass_region            # MH regional mass
#                              x  rel_volume_region(species)# NT/EH/MH size ratio
#
#   i.e. we hold the per-gram metabolic intensity fixed at the modern-human
#   value (the only species measured) and let the endocast size ratios scale
#   the regional mass. The species contrast is therefore driven entirely by
#   the Kochiyama relative volumes; rCMRGlc and any density/unit constant
#   cancel in the NT/MH and EH/MH budget RATIOS.
#
# METHOD (base R)
#   1. Build a region crosswalk between Heiss regions and Kochiyama
#      parcels. Kochiyama parcels are cortical/cerebellar sub-regions;
#      several aggregate up to one Heiss lobe. Aggregation of relative
#      volumes to the Heiss level is VOLUME-WEIGHTED using the MH absolute
#      sub-region volumes from the Kochiyama Figure-3 legend, so the lobe
#      relative volume is the true volume ratio of the summed sub-regions.
#      The sensorimotor parcel ("Sm") straddles the frontal/parietal border
#      and is split 50/50 between the frontal and parietal lobes.
#   2. Matched regions only: the 6 Heiss regions with Kochiyama coverage
#      (frontal, parietal, temporal, occipital lobes; cerebellar cortex;
#      vermis). Heiss subcortical / brain-stem / white-matter regions have
#      no endocast size data and are excluded (no size assumption imposed).
#   3. For each matched region and each species (NT, EH, MH), multiply the
#      MH glucose budget by the species relative volume, then total across
#      regions and take the NT/MH and EH/MH ratios.
#
# INPUTS
#   data_raw/Heiss_etal_2004_TABLE1.csv            (modern-human rCMRGlc)
#   data_raw/Kochiyama_etal_2018_Figure3.csv       (NT/EH/MH relative volumes)
#   data_raw/Kochiyama_etal_2018_Figure3legend.csv (MH absolute sub-region vols)
#   data_raw/Kochiyama_etal_2018_Figure3B.csv       (DIRECT absolute cerebellar
#     cortex volume, cc, L/R x group; Ce A + Ce P). Fig 3B ratio NT 0.925 /
#     EH 1.004. *** USED AS AN INDEPENDENT CHECK ONLY *** -- it does NOT enter
#     any species or specimen budget that feeds the figures/tables (see below).
#   data_raw/Kochiyama_etal_2018_Figure3A.csv       (ICV-size-adjusted relative
#     volumes; same 13-region content as Figure3.csv). *** THE SINGLE SHAPE
#     SOURCE for every budget: species averages AND specimens use Fig3A shape,
#     scaled to absolute size via the Fig3-legend MH volumes. ***
#     -- Kochiyama files copied from
#        Evo-M1-Trait-Data/Kochiyama_etal_2018/ on 2026-07-07.
#
# VOLUME BASIS (CSF)
#   The Fig3-legend MH volumes are PARCELLATED BRAIN-TISSUE region volumes
#   (grey+white; frontal, parietal, ... cerebellar cortex, vermis), NOT
#   intracranial volume -- cerebrospinal fluid is already EXCLUDED. The density
#   constant below (1.036 g/cc) therefore correctly treats them as tissue. The
#   group cerebral/cerebellar volumes used for absolute-size scaling are the
#   Kochiyama reconstructed brain (tissue) volumes on the same CSF-excluded basis.
#
# OUTPUTS
#   data_intermediate/s4_region_crosswalk.csv        region correspondence + weights
#   data_intermediate/s4_endocranial_region_budgets.csv  per-region budgets by species
#   tables/s4/endocranial_budget_summary.csv         totals and NT/MH, EH/MH ratios
# ============================================================

setwd("~/Library/CloudStorage/Dropbox/COLLABORATIVE/Do expensive brain regions increase less in humans/analyses_metabol_rate_structure")

# Assumed constants (only affect ABSOLUTE budget units, NOT the species ratios)
BRAIN_DENSITY_G_PER_CC <- 1.036   # assumed fresh brain-tissue density (~1.03-1.04 g/cc); cancels in the species ratios
# Heiss rCMRGlc is per 100 g/min, so mass enters in units of (100 g):
#   regional glucose (umol/min) = rCMRGlc * (mass_g / 100)

# -----------------------------
# 1. Read inputs
# -----------------------------

heiss <- read.csv("data_raw/Heiss_etal_2004_TABLE1.csv",
                  stringsAsFactors = FALSE, check.names = FALSE)

koch <- read.csv("data_raw/Kochiyama_etal_2018_Figure3A.csv",
                 stringsAsFactors = FALSE, check.names = FALSE)

leg <- read.csv("data_raw/Kochiyama_etal_2018_Figure3legend.csv",
                stringsAsFactors = FALSE, check.names = FALSE)

# Convenience lookups keyed by Kochiyama region code
koch_rel <- koch[, c("Region_code", "NT_rel", "EH_rel", "MH_rel")]
mh_vol_cc <- setNames(leg$MH_mean_Vol.cc, leg$Region_code)   # MH sub-region vol (cc)

# rCMRGlc lookup keyed by Heiss region name
heiss_rcmr <- setNames(heiss$`Both hemispheres Mean`, heiss$Region)

# -----------------------------
# 2. Region crosswalk: Heiss region <- Kochiyama parcels (+ split weight)
# -----------------------------
# split_weight = fraction of a Kochiyama parcel assigned to the Heiss region.
# All are 1 except the sensorimotor parcel "Sm", split 0.5 frontal / 0.5 parietal.

crosswalk <- data.frame(
  Heiss_region  = c("Frontal lobe","Frontal lobe","Frontal lobe","Frontal lobe",
                    "Parietal lobe","Parietal lobe","Parietal lobe",
                    "Temporal lobe","Temporal lobe",
                    "Occipital lobe","Occipital lobe",
                    "Cerebellar cortex","Cerebellar cortex",
                    "Vermis"),
  Koch_code     = c("Fr SM","Fr I","Fr O","Sm",
                    "Pa SI","Pa TP","Sm",
                    "Te SM","Te I",
                    "Oc SM","Oc I",
                    "Ce A","Ce P",
                    "Ce V"),
  split_weight  = c(1, 1, 1, 0.5,
                    1, 1, 0.5,
                    1, 1,
                    1, 1,
                    1, 1,
                    1),
  stringsAsFactors = FALSE
)

# Attach MH sub-region volume and relative volumes to each crosswalk row
crosswalk$MH_vol_cc <- mh_vol_cc[crosswalk$Koch_code]
crosswalk <- merge(crosswalk, koch_rel,
                   by.x = "Koch_code", by.y = "Region_code",
                   all.x = TRUE, sort = FALSE)

# Effective (split-adjusted) MH volume contributed by each parcel to its lobe
crosswalk$eff_vol_cc <- crosswalk$MH_vol_cc * crosswalk$split_weight

stopifnot(!any(is.na(crosswalk$MH_vol_cc)),   # every parcel found in legend
          !any(is.na(crosswalk$NT_rel)))      # every parcel found in Figure 3

write.csv(crosswalk[order(crosswalk$Heiss_region), ],
          "data_intermediate/s4_region_crosswalk.csv", row.names = FALSE)

# -----------------------------
# 3. Aggregate to Heiss-region level
# -----------------------------
# For each Heiss region:
#   MH_vol_cc(region)     = sum of effective sub-region MH volumes
#   rel_volume(species)   = volume-weighted mean of parcel relative volumes,
#                           weighting by the parcel's effective MH volume.
#     (This equals the true summed-volume ratio: because a parcel's species
#      volume = MH_vol * rel, the lobe ratio = sum(MH_vol*rel)/sum(MH_vol).)

regions <- unique(crosswalk$Heiss_region)

agg <- data.frame(Heiss_region = regions, stringsAsFactors = FALSE)

wmean <- function(x, w) sum(x * w) / sum(w)

agg$MH_vol_cc <- sapply(regions, function(r) {
  s <- crosswalk[crosswalk$Heiss_region == r, ]
  sum(s$eff_vol_cc)
})
for (sp in c("NT_rel","EH_rel","MH_rel")) {
  agg[[sp]] <- sapply(regions, function(r) {
    s <- crosswalk[crosswalk$Heiss_region == r, ]
    wmean(s[[sp]], s$eff_vol_cc)
  })
}

# Attach modern-human rCMRGlc for each Heiss region
agg$rCMRGlc <- heiss_rcmr[agg$Heiss_region]
stopifnot(!any(is.na(agg$rCMRGlc)))   # every matched Heiss region has a rate

# -----------------------------
# 4. Regional glucose budgets by species
# -----------------------------
#   MH regional mass (g)      = MH_vol_cc * density
#   MH regional glucose       = rCMRGlc * mass_g / 100      (umol glucose / min)
#   species regional glucose  = MH regional glucose * rel_volume(species)

agg$MH_mass_g       <- agg$MH_vol_cc * BRAIN_DENSITY_G_PER_CC
agg$budget_MH_base  <- agg$rCMRGlc * agg$MH_mass_g / 100   # umol/min, MH size

agg$budget_NT <- agg$budget_MH_base * agg$NT_rel
agg$budget_EH <- agg$budget_MH_base * agg$EH_rel
agg$budget_MH <- agg$budget_MH_base * agg$MH_rel           # == budget_MH_base

region_budgets <- agg[, c("Heiss_region","rCMRGlc","MH_vol_cc",
                          "NT_rel","EH_rel","MH_rel",
                          "budget_NT","budget_EH","budget_MH")]
region_budgets[ , -1] <- lapply(region_budgets[ , -1], round, 4)
write.csv(region_budgets,
          "data_intermediate/s4_endocranial_region_budgets.csv", row.names = FALSE)

# -----------------------------
# 5. Totals and species ratios
# -----------------------------

tot_NT <- sum(agg$budget_NT)
tot_EH <- sum(agg$budget_EH)
tot_MH <- sum(agg$budget_MH)

summary_tab <- data.frame(
  species              = c("Neanderthal (NT)", "Early modern human (EH)", "Modern human (MH)"),
  code                 = c("NT", "EH", "MH"),
  total_budget_umol_min = round(c(tot_NT, tot_EH, tot_MH), 2),
  ratio_to_MH          = round(c(tot_NT, tot_EH, tot_MH) / tot_MH, 4),
  pct_diff_from_MH     = round((c(tot_NT, tot_EH, tot_MH) / tot_MH - 1) * 100, 2),
  stringsAsFactors = FALSE
)

write.csv(summary_tab, "tables/s4/endocranial_budget_summary.csv", row.names = FALSE)

# -----------------------------
# 5b. Per-specimen and species-average ABSOLUTE-SIZE budgets
# -----------------------------
# Sections 1-5 use the Kochiyama relative volumes, which are ICV-SIZE-ADJUSTED:
# they describe regional SHAPE (proportional allocation at matched overall
# brain size), so their budget ratios are ~1 for every group. Here we restore
# absolute SIZE using each specimen's own cerebral & cerebellar volume.
#
#   budget(specimen) = sum_region [ budget_MH_base(region)          # MH per-region glucose
#                                   x rel_group(region)             # group SHAPE (Kochiyama)
#                                   x size_factor(region, specimen) ]
#   size_factor = cerebral_cc  / CEREBRAL_MH   for cortical regions
#               = cerebellar_cc / CEREBELLAR_MH for cerebellar regions
#
# Because rel_group is size-free, multiplying by the specimen's measured size
# does not double-count. A specimen is assigned its GROUP's shape (Skhul 5 uses
# the EH pattern) but its OWN size. Budgets are expressed as a ratio to the MH
# species reference (= sum of budget_MH_base, i.e. shape=MH, size=1).

CEREBRAL_MH   <- 1097   # Kochiyama MH mean cerebral volume (cc)
CEREBELLAR_MH <-  149   # Kochiyama MH mean cerebellar volume (cc)

# Classify each matched Heiss region as cortical or cerebellar tissue
agg$tissue <- ifelse(agg$Heiss_region %in% c("Cerebellar cortex", "Vermis"),
                     "cerebellar", "cortical")

rel_by_group <- list(NT = agg$NT_rel, EH = agg$EH_rel, MH = agg$MH_rel)
budget_MH_ref <- sum(agg$budget_MH_base)     # MH species reference (umol/min)

specimen_budget <- function(cerebral_cc, cerebellar_cc, group) {
  size_factor <- ifelse(agg$tissue == "cortical",
                        cerebral_cc  / CEREBRAL_MH,
                        cerebellar_cc / CEREBELLAR_MH)
  sum(agg$budget_MH_base * rel_by_group[[group]] * size_factor)
}

# Read the authoritative fossil specimen table (compiled from Kochiyama et al.
# 2018 main-text; per-specimen dates in years BP + reconstructed cerebral and
# cerebellar volumes + cerebellum:cerebrum ratio). Normalise the column names
# the rest of this script uses (ages in kya, *_cc volumes, group code).
spec_raw <- read.csv("data_raw/Kochiyama_etal_2018_FossilSpecimensText.csv",
                     stringsAsFactors = FALSE, check.names = FALSE,
                     encoding = "UTF-8")
spec <- data.frame(
  specimen      = spec_raw$Specimen,
  group         = spec_raw$Taxon_code,                 # NT / EH
  age_min_kya   = spec_raw$date_min_yBP / 1000,
  age_max_kya   = spec_raw$date_max_yBP / 1000,
  age_mid       = spec_raw$date_mean_yBP / 1000,       # paper's own midpoint
  cerebral_cc   = spec_raw$Cerebrum_Vol.cc,
  cerebellar_cc = spec_raw$Cerebellum_Vol.cc,
  cc_ratio      = spec_raw$Cerebellum_Cerebrum_ratio,  # cerebellum : cerebrum
  stringsAsFactors = FALSE)
spec$total_cc <- spec$cerebral_cc + spec$cerebellar_cc

# --- specimen sex estimates -------------------------------------------------
# Sex estimates compiled from the dedicated demographics lookup
# (material_other_sources/Kochiyama_fossil_demographics_journal_sources.xlsx),
# one journal citation per specimen (see that file's "Preferred journal
# citation" column). Codes: M / F when a sex is reported; "M?" / "F?" when the
# assignment is contested; "?" when unknown. Qafzeh 9 is coded "M?" because the
# 2021 reassessment reports male but older assessments were debated/female.
spec_sex <- c(
  "Amud 1"                   = "M",
  "La Chapelle-aux-Saints 1" = "M",
  "La Ferrassie 1"           = "M",
  "Forbes' Quarry 1"         = "F",
  "Qafzeh 9"                 = "M?",
  "Skhul 5"                  = "M",
  "Mladeč 1"                 = "F",
  "Cro-Magnon 1"             = "M")
spec$sex <- unname(spec_sex[spec$specimen])
spec$sex[is.na(spec$sex)] <- "?"

# --- per-specimen budgets ---
spec$budget_umol_min  <- mapply(specimen_budget,
                                spec$cerebral_cc, spec$cerebellar_cc, spec$group)
spec$budget_ratio_MH  <- spec$budget_umol_min / budget_MH_ref

specimen_budgets <- spec[, c("specimen","group","age_min_kya","age_max_kya",
                             "cerebral_cc","cerebellar_cc","total_cc",
                             "budget_umol_min","budget_ratio_MH")]
specimen_budgets$budget_umol_min <- round(specimen_budgets$budget_umol_min, 2)
specimen_budgets$budget_ratio_MH <- round(specimen_budgets$budget_ratio_MH, 4)
# write.csv escapes UTF-8 strings to <U+010D> when the R session locale is not
# UTF-8, mangling the c-caron in "Mladec 1". Build the CSV lines by hand and
# write them as raw bytes (useBytes = TRUE) so the UTF-8 stays intact.
write_csv_utf8 <- function(df, path) {
  q <- function(x) ifelse(grepl("[,\"]", x), paste0("\"", gsub("\"","\"\"",x), "\""), x)
  body <- apply(df, 1, function(r) paste(q(as.character(r)), collapse = ","))
  writeLines(c(paste(names(df), collapse = ","), body), path, useBytes = TRUE)
}
write_csv_utf8(specimen_budgets, "data_intermediate/s4_specimen_budgets.csv")

# --- species-average budgets from GROUP-MEAN volumes (+ SD propagation) ---
# IDENTICAL MODEL to the per-specimen function above: Fig3A shape everywhere,
# scaled to the group's absolute cerebral/cerebellar size. The budget is linear
# in the group volumes, so it splits into three additive tissue terms:
#   budget_g = A_g  *(cerebral_g  /CEREBRAL_MH)   cortical lobes    (Fig3A shape x cerebral size)
#            + CeC_g*(cerebellar_g/CEREBELLAR_MH) cerebellar cortex (Fig3A shape x cerebellar size)
#            + V_g  *(cerebellar_g/CEREBELLAR_MH) vermis            (Fig3A shape x cerebellar size)
# Figure 3B is NOT used here -- it is computed later purely as an independent
# check (see "Fig3B independent check" block below).
gm <- read.csv("data_raw/Kochiyama_etal_2018_FossilSpecimensText_groupmeans.csv",
               stringsAsFactors = FALSE, check.names = FALSE)
gm <- gm[match(c("NT","EH","MH"), gm$Group), ]
grp_vol <- data.frame(
  group      = gm$Group, n = gm$n,
  cerebral   = gm$Cerebrum_mean_cc,   cerebral_sd   = gm$Cerebrum_sd_cc,
  cerebellar = gm$Cerebellum_mean_cc, cerebellar_sd = gm$Cerebellum_sd_cc,
  stringsAsFactors = FALSE)

# Per-tissue MH glucose bases and Fig3A shape sums (per group):
A_g <- sapply(grp_vol$group, function(g)                      # cortical lobes
  sum(agg$budget_MH_base[agg$tissue == "cortical"] * rel_by_group[[g]][agg$tissue == "cortical"]))
is_verm <- agg$Heiss_region == "Vermis"
V_g <- sapply(grp_vol$group, function(g)                      # vermis
  sum(agg$budget_MH_base[is_verm] * rel_by_group[[g]][is_verm]))
is_cec_reg <- agg$Heiss_region == "Cerebellar cortex"
CeC_g <- sapply(grp_vol$group, function(g)                    # cerebellar cortex (Fig3A shape)
  sum(agg$budget_MH_base[is_cec_reg] * rel_by_group[[g]][is_cec_reg]))

# Species budget = cortical (Fig3A shape x cerebral size)
#               + cerebellar cortex (Fig3A shape x cerebellar size)
#               + vermis           (Fig3A shape x cerebellar size)
grp_vol$budget_umol_min <-
  A_g   * (grp_vol$cerebral   / CEREBRAL_MH) +
  CeC_g * (grp_vol$cerebellar / CEREBELLAR_MH) +
  V_g   * (grp_vol$cerebellar / CEREBELLAR_MH)

# Population SD (individual spread): the cortical term varies with cerebral
# volume; both cerebellar terms (cortex + vermis) vary with cerebellar volume.
sd_cort  <- (A_g / CEREBRAL_MH)               * grp_vol$cerebral_sd
sd_cec   <- (CeC_g / CEREBELLAR_MH)           * grp_vol$cerebellar_sd
sd_verm  <- (V_g / CEREBELLAR_MH)             * grp_vol$cerebellar_sd
grp_vol$budget_sd <- sqrt(sd_cort^2 + sd_cec^2 + sd_verm^2)
# 95% CI on the MEAN budget: SEM = SD/sqrt(n), half-width = t(n-1,0.975)*SEM.
# CI answers "how well is the group mean pinned down"; fossil n = 4 -> wide,
# MH n = 1185 -> tight.
sem_cort <- sd_cort / sqrt(grp_vol$n)
sem_cec  <- sd_cec  / sqrt(grp_vol$n)
sem_verm <- sd_verm / sqrt(grp_vol$n)
grp_vol$budget_sem <- sqrt(sem_cort^2 + sem_cec^2 + sem_verm^2)
tcrit <- qt(0.975, df = pmax(grp_vol$n - 1, 1))
grp_vol$budget_ci95 <- tcrit * grp_vol$budget_sem      # half-width, umol/min

grp_vol$budget_ratio_MH     <- grp_vol$budget_umol_min / budget_MH_ref
grp_vol$budget_ratio_MH_sd  <- grp_vol$budget_sd       / budget_MH_ref
grp_vol$budget_ratio_MH_ci95 <- grp_vol$budget_ci95    / budget_MH_ref

species_budgets <- data.frame(
  group                  = grp_vol$group,
  n                      = grp_vol$n,
  cerebral_cc            = grp_vol$cerebral,
  cerebellar_cc          = grp_vol$cerebellar,
  budget_umol_min        = round(grp_vol$budget_umol_min, 2),
  budget_sd              = round(grp_vol$budget_sd, 2),
  budget_ci95_halfwidth  = round(grp_vol$budget_ci95, 2),
  budget_ratio_MH        = round(grp_vol$budget_ratio_MH, 4),
  budget_ratio_MH_sd     = round(grp_vol$budget_ratio_MH_sd, 4),
  budget_ratio_MH_ci95   = round(grp_vol$budget_ratio_MH_ci95, 4),
  stringsAsFactors = FALSE)
write.csv(species_budgets,
          "tables/s4/species_absolute_budgets.csv", row.names = FALSE)

# -----------------------------
# 5c. Fig3B INDEPENDENT CHECK (does NOT feed any budget/figure)
# -----------------------------
# Figure 3B gives the DIRECT digitized absolute cerebellar-cortex volume
# (Ce A + Ce P, L+R summed) per GROUP. We use it ONLY to sanity-check the
# Fig3A-derived cerebellar-cortex term against an independent measurement.
# It is deliberately kept out of grp_vol$budget_umol_min above.
f3b <- read.csv("data_raw/Kochiyama_etal_2018_Figure3B.csv", stringsAsFactors = FALSE)
is_cec <- f3b$Region_code %in% c("Ce A", "Ce P")
f3b_cortex <- tapply(f3b$Volume_cc[is_cec], f3b$Group_code[is_cec], sum)[c("NT","EH","MH")]
f3b_ratio  <- f3b_cortex / f3b_cortex["MH"]      # NT 0.925, EH 1.004, MH 1.000

# Fig3A-implied cerebellar-cortex ratio to MH (shape x cerebellar size), for
# the same three groups -- this is what actually enters the budget.
cebvol <- setNames(grp_vol$cerebellar, grp_vol$group)
f3a_cec_ratio <- setNames(
  (CeC_g * (cebvol[grp_vol$group] / CEREBELLAR_MH)) /
    (CeC_g["MH"] * (cebvol["MH"] / CEREBELLAR_MH)),
  grp_vol$group)[c("NT","EH","MH")]

fig3b_check <- data.frame(
  group                 = c("NT","EH","MH"),
  fig3A_cerebellarcortex_ratio_MH = round(as.numeric(f3a_cec_ratio), 4),  # used in budget
  fig3B_cerebellarcortex_ratio_MH = round(as.numeric(f3b_ratio), 4),      # independent check
  abs_diff              = round(abs(as.numeric(f3a_cec_ratio) - as.numeric(f3b_ratio)), 4),
  stringsAsFactors = FALSE)
write.csv(fig3b_check, "tables/s4/fig3B_cerebellar_cortex_check.csv", row.names = FALSE)

# -----------------------------
# 6. Console report
# -----------------------------

cat("\n=== s4_endocranial: matched regions (Heiss <-> Kochiyama) ===\n")
print(region_budgets, row.names = FALSE)

cat("\n=== Total brain energy budget across matched regions ===\n")
print(summary_tab, row.names = FALSE)

cat("\nNT is", round((tot_NT/tot_MH - 1)*100, 2),
    "% relative to modern human across the matched regions.\n")
cat("EH is", round((tot_EH/tot_MH - 1)*100, 2),
    "% relative to modern human across the matched regions.\n")
cat("\n(Budget ratios are invariant to brain density and rCMRGlc units;\n",
    "the species contrast is driven by the Kochiyama endocast size ratios.)\n")

cat("\n=== Per-specimen ABSOLUTE-SIZE budgets (group shape x specimen size) ===\n")
print(specimen_budgets, row.names = FALSE)

cat("\n=== Species-average ABSOLUTE-SIZE budgets (from group-mean volumes) ===\n")
print(species_budgets, row.names = FALSE)
cat("\nUnlike the size-adjusted ratios above (all ~1.00), restoring absolute\n",
    "brain size gives NT ~", round((species_budgets$budget_ratio_MH[1]-1)*100,1),
    "% and EH ~+", round((species_budgets$budget_ratio_MH[2]-1)*100,1),
    "% vs MH, with wide overlapping s.d. (larger fossil brains, same shape).\n")

# -----------------------------
# 7. Figure: per-region budget by species (base R barplot)
# -----------------------------

mat <- t(as.matrix(agg[, c("budget_NT","budget_EH","budget_MH")]))
colnames(mat) <- agg$Heiss_region
rownames(mat) <- c("NT","EH","MH")

# ---------------------------------------------------------------------------
# Figure export helper: render one base-graphics figure to BOTH a raster PNG
# (for slides / PowerPoint) and a vector PDF (for print). `draw` is a function
# holding the plotting calls; it is replayed once per device so the two files
# are identical. cairo_pdf sizes are in inches = pixels / res.
# ---------------------------------------------------------------------------
save_png_pdf <- function(stem, draw, width, height, res = 220) {
  png(paste0(stem, ".png"), width = width, height = height, res = res)
  draw(); dev.off()
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width / res, height = height / res)
  draw(); dev.off()
}

.draw_region_budgets <- function() {
op <- par(mar = c(7.5, 4.6, 3.2, 1.2), mgp = c(3, 0.7, 0))
cols <- c(NT = "#4C72B0", EH = "#DD8452", MH = "#55A868")
bp <- barplot(mat, beside = TRUE, col = cols, border = NA,
              las = 2, ylab = "Regional glucose budget (umol glucose / min)",
              cex.names = 0.9, ylim = c(0, max(mat) * 1.12))
title(main = "Endocast-scaled regional brain energy budgets by hominin group",
      cex.main = 1.05, font.main = 1)
legend("topright",
       legend = c("Neanderthal (NT)","Early modern human (EH)","Modern human (MH)"),
       fill = cols, border = NA, bty = "n", cex = 0.9)
mtext(paste0("Total budget ratio to MH:  NT = ", summary_tab$ratio_to_MH[1],
             "   EH = ", summary_tab$ratio_to_MH[2],
             "   (MH = 1)"),
      side = 1, line = 6.2, cex = 0.8, adj = 0)
par(op)
}
save_png_pdf("figs/s4/endocranial_region_budgets",
             .draw_region_budgets, width = 2100, height = 1350, res = 220)

# -----------------------------
# 7b. Figure: regions as parts of total brain VOLUME (donut)
# -----------------------------
# Shows each region sized by its modern-human volume as a fraction of the whole
# brain. The 6 regions that carry BOTH a Kochiyama volume and a Heiss metabolic
# rate (i.e. that enter the energy budget) are coloured; regions for which no
# metabolic rate is used -- subcortical nuclei, white matter, insula,
# hippocampus, brain stem, cerebellar white matter/dentate -- are shown GREY:
# they contribute to brain volume but not to the metabolic-rate budget.
CEREBRUM_MH_TOTAL   <- 1097   # MH cerebrum group mean (cc)
CEREBELLUM_MH_TOTAL <- 149    # MH cerebellum group mean (cc)
BRAIN_MH_TOTAL      <- CEREBRUM_MH_TOTAL + CEREBELLUM_MH_TOTAL   # 1246 cc

# matched (budgeted) regions, split into cerebral vs cerebellar
is_cbll  <- grepl("[Cc]erebell|[Vv]ermis", agg$Heiss_region)
cort_matched_vol <- sum(agg$MH_vol_cc[!is_cbll])
cbll_matched_vol <- sum(agg$MH_vol_cc[ is_cbll])
cerebral_grey  <- CEREBRUM_MH_TOTAL   - cort_matched_vol   # subcortical, WM, insula, hippocampus
cerebellar_grey<- CEREBELLUM_MH_TOTAL - cbll_matched_vol   # cerebellar WM + dentate

pie_df <- data.frame(
  label = c(agg$Heiss_region,
            "Other cerebral tissue\n(subcortical, white matter,\ninsula, hippocampus)",
            "Other cerebellar tissue\n(white matter, dentate)"),
  vol   = c(agg$MH_vol_cc, cerebral_grey, cerebellar_grey),
  budgeted = c(rep(TRUE, nrow(agg)), FALSE, FALSE),
  stringsAsFactors = FALSE)
# order: budgeted regions first (largest -> smallest), then grey remainders
ord <- c(order(-pie_df$vol[pie_df$budgeted]),
         nrow(agg) + order(-pie_df$vol[!pie_df$budgeted]))
pie_df <- pie_df[ord, ]
pie_df$frac <- pie_df$vol / sum(pie_df$vol)

reg_pal <- c("Frontal lobe"="#4C72B0","Parietal lobe"="#55A868",
             "Temporal lobe"="#C44E52","Occipital lobe"="#8172B3",
             "Cerebellar cortex"="#DD8452","Vermis"="#937860")
pie_df$col <- ifelse(pie_df$budgeted, reg_pal[pie_df$label], NA)
pie_df$col[!pie_df$budgeted] <- c("#B9B9B9","#DBDBDB")[seq_len(sum(!pie_df$budgeted))]

# annular-wedge donut (manual; no extra package)
ann_wedge <- function(cx, cy, r_in, r_out, a0, a1, col, border="white", lwd=1.5) {
  tt <- seq(a0, a1, length.out = max(2, round((a1-a0)/(pi/180))))
  xx <- c(cx + r_out*cos(tt), cx + r_in*rev(cos(tt)))
  yy <- c(cy + r_out*sin(tt), cy + r_in*rev(sin(tt)))
  polygon(xx, yy, col = col, border = border, lwd = lwd)
}

png("figs/s4/endocranial_region_volume_pie.png", width = 2100, height = 1500, res = 220)
op <- par(mar = c(2.2, 1, 3.4, 1))
plot(NA, xlim = c(-1.72, 1.72), ylim = c(-1.35, 1.35), asp = 1,
     axes = FALSE, xlab = "", ylab = "")
r_in <- 0.55; r_out <- 1.05
ang <- pi/2 - 2*pi*c(0, cumsum(pie_df$frac))     # start at top, go clockwise
for (i in seq_len(nrow(pie_df))) {
  ann_wedge(0, 0, r_in, r_out, ang[i+1], ang[i], col = pie_df$col[i])
  mid <- (ang[i] + ang[i+1]) / 2
  lab_r <- r_out + 0.10
  lx <- lab_r*cos(mid); ly <- lab_r*sin(mid)
  adj_x <- if (cos(mid) < -0.10) 1 else if (cos(mid) > 0.10) 0 else 0.5
  txt <- sprintf("%s\n%.0f cc  (%.0f%%)",
                 gsub("\n", " ", pie_df$label[i]), pie_df$vol[i], 100*pie_df$frac[i])
  # keep the two long grey labels multi-line for readability
  if (!pie_df$budgeted[i]) txt <- sprintf("%s\n%.0f cc  (%.0f%%)",
                                          pie_df$label[i], pie_df$vol[i], 100*pie_df$frac[i])
  text(lx, ly, txt, cex = 0.62, adj = c(adj_x, 0.5),
       col = if (pie_df$budgeted[i]) "grey15" else "grey40", xpd = NA)
  # leader line
  segments(r_out*cos(mid), r_out*sin(mid), (r_out+0.06)*cos(mid), (r_out+0.06)*sin(mid),
           col = "grey70", lwd = 0.8)
}
matched_pct <- 100*sum(pie_df$frac[pie_df$budgeted])
text(0,  0.14, sprintf("%.0f cc", BRAIN_MH_TOTAL), cex = 1.5, font = 2, col = "grey15")
text(0, -0.06, "total MH brain", cex = 0.8, col = "grey40")
text(0, -0.24, sprintf("%.0f%% of volume\nhas a rate", matched_pct), cex = 0.66, col = "grey30")
title(main = "Modern-human brain by VOLUME: which regions enter the energy budget",
      cex.main = 1.02, font.main = 1)
legend("bottomleft",
       legend = c("In energy budget (Heiss rate available)",
                  "Volume only -- no metabolic rate used"),
       fill = c("#4C72B0", "#B9B9B9"), border = "white", bty = "n", cex = 0.72)
mtext("Wedge size and all percentages are shares of brain VOLUME (Kochiyama Fig 3 legend + cerebrum/cerebellum group means), NOT of energy budget. Grey regions have no metabolic rate, so a whole-brain budget denominator is undefined -- % of total budget cannot be computed from these data.",
      side = 1, line = 0.9, cex = 0.53, col = "grey40")
par(op)
dev.off()

# -----------------------------
# 7c. Figure: three-species comparison of brain volume composition (donuts)
# -----------------------------
# Same donut as 7b, drawn for NT, EH and MH side by side. Each species' regional
# volume = MH parcel volume x Fig 3A relative (shape) x that species' size factor
# (cerebrum ratio for cortical regions, cerebellum ratio for cerebellar regions);
# the grey "other tissue" wedge is the per-species remainder of cerebrum/cerebellum.
sp_tot <- data.frame(sp = c("NT","EH","MH"),
                     cerebrum = c(1161,1135,1097), cerebellum = c(149,153,149),
                     stringsAsFactors = FALSE)
sp_size <- list(NT = c(cort=1161/1097, cbll=149/149),
                EH = c(cort=1135/1097, cbll=153/149),
                MH = c(cort=1097/1097, cbll=149/149))
cw2 <- crosswalk
cw2$is_cbll <- grepl("[Cc]erebell|[Vv]ermis", cw2$Heiss_region)
region_vol_sp <- function(sp) {
  rel <- cw2[[paste0(sp,"_rel")]]
  sf  <- ifelse(cw2$is_cbll, sp_size[[sp]]["cbll"], sp_size[[sp]]["cort"])
  tapply(cw2$MH_vol_cc * cw2$split_weight * rel * sf, cw2$Heiss_region, sum)
}
reg_order <- c("Frontal lobe","Parietal lobe","Temporal lobe","Occipital lobe",
               "Cerebellar cortex","Vermis")
reg_pal <- c("Frontal lobe"="#4C72B0","Parietal lobe"="#55A868",
             "Temporal lobe"="#C44E52","Occipital lobe"="#8172B3",
             "Cerebellar cortex"="#DD8452","Vermis"="#937860")

build_pie <- function(sp) {
  Vr <- region_vol_sp(sp)[reg_order]
  cort_m <- sum(Vr[1:4]); cbll_m <- sum(Vr[5:6])
  cerebral_grey  <- sp_tot$cerebrum[sp_tot$sp==sp]  - cort_m
  cerebellar_grey<- sp_tot$cerebellum[sp_tot$sp==sp]- cbll_m
  d <- data.frame(label = c(reg_order,"Other cerebral","Other cerebellar"),
                  vol = c(Vr, cerebral_grey, cerebellar_grey),
                  budgeted = c(rep(TRUE,6), FALSE, FALSE), stringsAsFactors = FALSE)
  d$col <- c(reg_pal[reg_order], "#B9B9B9", "#DBDBDB")
  d$frac <- d$vol/sum(d$vol)
  d
}
if (!exists("ann_wedge")) {
  ann_wedge <- function(cx, cy, r_in, r_out, a0, a1, col, border="white", lwd=1.2) {
    tt <- seq(a0, a1, length.out = max(2, round((a1-a0)/(pi/180))))
    polygon(c(cx + r_out*cos(tt), cx + r_in*rev(cos(tt))),
            c(cy + r_out*sin(tt), cy + r_in*rev(sin(tt))),
            col = col, border = border, lwd = lwd)
  }
}
draw_one <- function(d, sp, brain_cc) {
  plot(NA, xlim = c(-1.25,1.25), ylim = c(-1.25,1.25), asp = 1, axes = FALSE, xlab="", ylab="")
  r_in <- 0.52; r_out <- 1.02
  ang <- pi/2 - 2*pi*c(0, cumsum(d$frac))
  for (i in seq_len(nrow(d))) {
    ann_wedge(0,0,r_in,r_out, ang[i+1], ang[i], col = d$col[i])
    if (d$budgeted[i] && d$frac[i] > 0.03) {
      mid <- (ang[i]+ang[i+1])/2
      text((r_in+r_out)/2*cos(mid), (r_in+r_out)/2*sin(mid),
           sprintf("%.0f", d$vol[i]), cex = 0.6, col = "white", font = 2)
    }
  }
  matched_pct <- 100*sum(d$frac[d$budgeted])
  text(0, 0.10, sprintf("%d cc", round(brain_cc)), cex = 1.15, font = 2, col = "grey15")
  text(0,-0.12, sprintf("%.0f%% of vol.", matched_pct), cex = 0.64, col = "grey35")
  sp_full <- c(NT="Neanderthal", EH="Early modern human", MH="Modern human")[sp]
  title(main = sprintf("%s (%s)", sp_full, sp), cex.main = 0.95, font.main = 1, line = 0.2)
}

png("figs/s4/endocranial_region_volume_pie_3species.png", width = 2550, height = 1150, res = 220)
op <- par(mfrow = c(1,3), mar = c(1.4,0.6,2.6,0.6), oma = c(4.2,0,3.0,0))
for (sp in c("NT","EH","MH")) {
  d <- build_pie(sp)
  draw_one(d, sp, sp_tot$cerebrum[sp_tot$sp==sp] + sp_tot$cerebellum[sp_tot$sp==sp])
}
mtext("Brain volume composition across hominin groups (wedge = modern-human volume x Fig 3A shape x species size)",
      outer = TRUE, cex = 0.92, font = 2, line = 0.6)
par(op); par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = TRUE)
plot(NA, xlim=0:1, ylim=0:1, axes=FALSE, xlab="", ylab="")
legend("bottom", legend = c(reg_order, "Other cerebral (volume only)", "Other cerebellar (volume only)"),
       fill = c(reg_pal[reg_order], "#B9B9B9", "#DBDBDB"), border = "white",
       bty = "n", cex = 0.68, ncol = 4, xpd = NA, inset = c(0, 0.005))
mtext("Wedge labels = cc; centre % = share of brain VOLUME with a metabolic rate, NOT of energy budget (grey has no rate, so total-budget share is undefined).",
      side = 1, outer = TRUE, cex = 0.55, col = "grey40", line = 2.6)
dev.off()

# -----------------------------
# 7d. Figure: donut sized by VOLUME, labelled with ABSOLUTE metabolic cost
# -----------------------------
# Same volume-proportioned wedges as 7b/7c, but each measured wedge is labelled
# with its ABSOLUTE glucose cost (umol glucose/min) = rate x (volume x 1.036/100),
# instead of a percentage. Grey regions have no Heiss rate, so they are labelled
# "no rate" -- we deliberately do NOT invent a cost for them. The centre shows the
# summed cost over the MEASURED regions (this is NOT a whole-brain total).
rate_by_region <- setNames(agg$rCMRGlc, agg$Heiss_region)
sp_size2 <- list(NT = c(cort=1161/1097, cbll=149/149),
                 EH = c(cort=1135/1097, cbll=153/149),
                 MH = c(cort=1097/1097, cbll=149/149))
cw3 <- crosswalk; cw3$is_cbll <- grepl("[Cc]erebell|[Vv]ermis", cw3$Heiss_region)
region_vol_sp2 <- function(sp) {
  rel <- cw3[[paste0(sp,"_rel")]]
  sf  <- ifelse(cw3$is_cbll, sp_size2[[sp]]["cbll"], sp_size2[[sp]]["cort"])
  tapply(cw3$MH_vol_cc * cw3$split_weight * rel * sf, cw3$Heiss_region, sum)
}
reg_order2 <- c("Frontal lobe","Parietal lobe","Temporal lobe","Occipital lobe",
                "Cerebellar cortex","Vermis")

build_cost_pie <- function(sp) {
  Vr <- region_vol_sp2(sp)[reg_order2]
  cost <- rate_by_region[reg_order2] * (Vr * BRAIN_DENSITY_G_PER_CC / 100)
  cort_m <- sum(Vr[1:4]); cbll_m <- sum(Vr[5:6])
  cerebral_grey  <- sp_tot$cerebrum[sp_tot$sp==sp]  - cort_m
  cerebellar_grey<- sp_tot$cerebellum[sp_tot$sp==sp]- cbll_m
  d <- data.frame(label = c(reg_order2,"Other cerebral","Other cerebellar"),
                  vol = c(Vr, cerebral_grey, cerebellar_grey),
                  cost = c(cost, NA, NA),
                  budgeted = c(rep(TRUE,6), FALSE, FALSE), stringsAsFactors = FALSE)
  d$col <- c(reg_pal[reg_order2], "#B9B9B9", "#DBDBDB")
  d$frac <- d$vol/sum(d$vol)
  d
}
draw_cost <- function(d, sp) {
  plot(NA, xlim = c(-1.55,1.55), ylim = c(-1.3,1.3), asp = 1, axes = FALSE, xlab="", ylab="")
  r_in <- 0.52; r_out <- 1.02
  ang <- pi/2 - 2*pi*c(0, cumsum(d$frac))
  for (i in seq_len(nrow(d))) {
    ann_wedge(0,0,r_in,r_out, ang[i+1], ang[i], col = d$col[i])
    mid <- (ang[i]+ang[i+1])/2
    lab_r <- r_out + 0.08
    adj_x <- if (cos(mid) < -0.10) 1 else if (cos(mid) > 0.10) 0 else 0.5
    if (d$budgeted[i]) {
      txt <- sprintf("%s\n%.0f umol/min", gsub("\n"," ",d$label[i]), d$cost[i])
      colr <- "grey15"
    } else {
      txt <- sprintf("%s\n(no rate; %.0f cc)", d$label[i], d$vol[i]); colr <- "grey45"
    }
    text(lab_r*cos(mid), lab_r*sin(mid), txt, cex = 0.6, adj = c(adj_x,0.5), col = colr, xpd = NA)
  }
  tot_cost <- sum(d$cost, na.rm = TRUE)
  text(0, 0.12, sprintf("%.0f", tot_cost), cex = 1.35, font = 2, col = "grey15")
  text(0,-0.06, "umol glucose/min", cex = 0.62, col = "grey40")
  text(0,-0.22, "(measured regions)", cex = 0.58, col = "grey45")
  sp_full <- c(NT="Neanderthal", EH="Early modern human", MH="Modern human")[sp]
  title(main = sprintf("%s (%s)", sp_full, sp), cex.main = 0.95, font.main = 1, line = 0.2)
}

png("figs/s4/endocranial_region_cost_pie_3species.png", width = 2700, height = 1200, res = 220)
op <- par(mfrow = c(1,3), mar = c(1.2,2.0,2.6,2.0), oma = c(5.2,0,3.2,0))
for (sp in c("NT","EH","MH")) draw_cost(build_cost_pie(sp), sp)
mtext("Absolute regional glucose cost across hominin groups  (wedge SIZE = brain volume; wedge LABEL = rate x mass, umol glucose/min)",
      outer = TRUE, cex = 0.9, font = 2, line = 0.7)
par(op); par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = TRUE)
plot(NA, xlim=0:1, ylim=0:1, axes=FALSE, xlab="", ylab="")
legend("bottom", legend = c(reg_order2, "Other cerebral (no rate)", "Other cerebellar (no rate)"),
       fill = c(reg_pal[reg_order2], "#B9B9B9", "#DBDBDB"), border = "white",
       bty = "n", cex = 0.66, ncol = 4, xpd = NA, inset = c(0, 0.005))
mtext("Wedge area is proportional to VOLUME; the label converts that volume to an absolute glucose cost via the Heiss rate. Centre total sums only the 6 measured regions -- it is NOT a whole-brain metabolic total (grey regions have no rate).",
      side = 1, outer = TRUE, cex = 0.53, col = "grey40", line = 2.7)
mtext("Point estimates from group-mean sizes; between-species differences are not statistically significant (n = 4 fossils/group).",
      side = 1, outer = TRUE, cex = 0.53, col = "grey40", line = 3.4)
dev.off()

# -----------------------------
# 7e. Figure: 7d + volume context annotations
# -----------------------------
# Identical to 7d (wedge size = volume, coloured labels = absolute glucose cost),
# but adds three volume references requested for teaching:
#   (1) under the centre total: the measured regions as a % of TOTAL brain volume;
#   (2) under each species title: that species' estimated TOTAL brain volume (cc);
#   (3) inside each coloured wedge: the % of total brain volume it represents.
draw_cost2 <- function(d, sp, brain_cc) {
  plot(NA, xlim = c(-1.55,1.55), ylim = c(-1.35,1.35), asp = 1, axes = FALSE, xlab="", ylab="")
  r_in <- 0.52; r_out <- 1.02
  ang <- pi/2 - 2*pi*c(0, cumsum(d$frac))
  for (i in seq_len(nrow(d))) {
    ann_wedge(0,0,r_in,r_out, ang[i+1], ang[i], col = d$col[i])
    mid <- (ang[i]+ang[i+1])/2
    lab_r <- r_out + 0.08
    adj_x <- if (cos(mid) < -0.10) 1 else if (cos(mid) > 0.10) 0 else 0.5
    if (d$budgeted[i]) {
      txt <- sprintf("%s\n%.0f umol/min", gsub("\n"," ",d$label[i]), d$cost[i]); colr <- "grey15"
    } else {
      txt <- sprintf("%s\n(no rate; %.0f cc)", d$label[i], d$vol[i]); colr <- "grey45"
    }
    text(lab_r*cos(mid), lab_r*sin(mid), txt, cex = 0.6, adj = c(adj_x,0.5), col = colr, xpd = NA)
    # (3) percent of total brain volume, written INSIDE the coloured wedge
    if (d$budgeted[i] && d$frac[i] > 0.03) {
      text((r_in+r_out)/2*cos(mid), (r_in+r_out)/2*sin(mid),
           sprintf("%.0f%%", 100*d$frac[i]), cex = 0.62, col = "white", font = 2)
    }
  }
  tot_cost   <- sum(d$cost, na.rm = TRUE)
  matched_pct<- 100*sum(d$frac[d$budgeted])
  text(0, 0.20, sprintf("%.0f", tot_cost), cex = 1.3, font = 2, col = "grey15")
  text(0, 0.03, "umol glucose/min", cex = 0.60, col = "grey40")
  text(0,-0.11, "(measured regions)", cex = 0.56, col = "grey45")
  # (1) measured regions as % of total brain volume
  text(0,-0.28, sprintf("= %.0f%% of total\nbrain volume", matched_pct),
       cex = 0.56, col = "grey30")
  sp_full <- c(NT="Neanderthal", EH="Early modern human", MH="Modern human")[sp]
  title(main = sprintf("%s (%s)", sp_full, sp), cex.main = 0.95, font.main = 1, line = 1.1)
  # (2) species total brain volume, under the title
  mtext(sprintf("total brain volume = %d cc", round(brain_cc)),
        side = 3, line = 0.05, cex = 0.66, col = "grey30")
}

.draw_cost_pie_3species <- function() {
op <- par(mfrow = c(1,3), mar = c(1.2,2.0,3.4,2.0), oma = c(5.2,0,3.2,0))
for (sp in c("NT","EH","MH")) {
  brain_cc <- sp_tot$cerebrum[sp_tot$sp==sp] + sp_tot$cerebellum[sp_tot$sp==sp]
  draw_cost2(build_cost_pie(sp), sp, brain_cc)
}
mtext("Absolute regional glucose cost across hominin groups  (wedge SIZE = brain volume; wedge LABEL = rate x mass, umol glucose/min; wedge % = share of total brain volume)",
      outer = TRUE, cex = 0.86, font = 2, line = 0.7)
par(op); par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = TRUE)
plot(NA, xlim=0:1, ylim=0:1, axes=FALSE, xlab="", ylab="")
legend("bottom", legend = c(reg_order2, "Other cerebral (no rate)", "Other cerebellar (no rate)"),
       fill = c(reg_pal[reg_order2], "#B9B9B9", "#DBDBDB"), border = "white",
       bty = "n", cex = 0.66, ncol = 4, xpd = NA, inset = c(0, 0.005))
mtext("Wedge area proportional to VOLUME; label = absolute glucose cost via Heiss rate; in-wedge % = that region's share of total brain volume. Centre total sums only the 6 measured regions -- NOT a whole-brain metabolic total (grey regions have no rate).",
      side = 1, outer = TRUE, cex = 0.53, col = "grey40", line = 2.7)
mtext("Point estimates from group-mean sizes; between-species differences are not statistically significant (n = 4 fossils/group).",
      side = 1, outer = TRUE, cex = 0.53, col = "grey40", line = 3.4)
}
save_png_pdf("figs/s4/endocranial_region_cost_pie_3species_volctx",
             .draw_cost_pie_3species, width = 2700, height = 1250, res = 220)

# =============================================================================
# PLOTTING  (reduced set -- 3 figures)
# =============================================================================
# The energy budget is encoded as MARKER SIZE on the volume-through-time
# scatter (bigger budget = bigger point), not as a separate budget-vs-age plot.
# Figures produced:
#   Fig A  volume_timeline.png         brain volume vs age; point size = budget
#   Fig A2 volume_timeline_logage.png  same, LOG10 age axis (spreads the fossils)
#   Fig C  species_budget_ci.png       species means +/- 95% CI (true scale)
# The earlier linear-"AMPLIFY" deviation variants have been removed: rescaling
# deviations about 1.00 while relabelling ticks with the true numbers adds no
# information. A log AGE axis (Fig A2) is the one genuinely useful re-scaling.
# All figures share ONE specimen marker scheme (colour = group, shape = specimen).

grp_col <- c(NT = "#4C72B0", EH = "#DD8452", MH = "#55A868")

# ---- specimen plotting attributes -----------------------------------------
nt_dated_mid  <- mean(spec$age_mid[spec$group == "NT" & !is.na(spec$age_mid)])
spec$age_plot <- ifelse(is.na(spec$age_mid), nt_dated_mid, spec$age_mid)
spec$undated  <- is.na(spec$age_mid)
spec$col      <- grp_col[spec$group]
fill_pch <- c(21, 22, 23, 24)                     # circle, square, diamond, tri
spec$pch <- ave(seq_len(nrow(spec)), spec$group,
                FUN = function(ix) fill_pch[seq_along(ix)])
mh_total_cc <- CEREBRAL_MH + CEREBELLAR_MH
MH_PCH <- 23

# ---- marker-size encoding: area scaled to the energy-budget ratio ----------
# Points are sized by their energy-budget ratio to the modern-human mean
# (bigger budget = bigger marker). Because the ratios span only a narrow band
# (~0.81-1.19), a strict area-proportional map is nearly invisible; instead
# marker AREA is mapped LINEARLY across the observed budget range onto a legible
# size band. The map is monotonic and shown in a size legend, so ordering and
# relative magnitude read correctly. MH (=1.00) is sized on the same scale.
SIZE_CEX  <- c(1.5, 4.6)                          # min / max marker cex
size_vals <- c(spec$budget_ratio_MH, 1.0)         # include MH reference
size_rng  <- range(size_vals)
cex_from_budget <- function(r)                    # budget ratio -> plotting cex
  sqrt(SIZE_CEX[1]^2 + (SIZE_CEX[2]^2 - SIZE_CEX[1]^2) *
         (r - size_rng[1]) / diff(size_rng))
spec$cex_budget <- cex_from_budget(spec$budget_ratio_MH)
mh_cex_budget   <- cex_from_budget(1.0)

# ---- shared "famous-figure" scatter ---------------------------------------
# Identity in a LEFT legend (shape=specimen, colour=group); age left->right with
# PRESENT on the right; y-axis on the RIGHT; dashed group-mean lines; shaded
# NT-MH divergence window. Point AREA encodes the energy budget (see above);
# a size legend at bottom-left decodes it. xlog = TRUE draws age on a LOG10 axis:
# age 0 is undefined on a log scale, so the modern-human marker and the "present"
# line move to a nominal recent age MH_LOG_AGE (kya); dated fossils then spread
# across the axis instead of bunching in the last ~135 ky.
MH_LOG_AGE <- 5     # kya; nominal "recent" position for MH on a log axis

famous_scatter <- function(yvar, ylab, main, mh_y, y_pad = 0.06,
                           note = "", xlog = FALSE) {
  yv <- spec[[yvar]]; mh <- mh_y
  yr <- range(c(yv, mh)); pad <- diff(yr) * y_pad
  ylim <- c(yr[1] - pad, yr[2] + pad)
  mh_x <- if (xlog) MH_LOG_AGE else 0                 # MH marker x-position
  xlim <- if (xlog) c(1000, 3) else c(850, -55)       # log: high->low kya
  logarg <- if (xlog) "x" else ""

  op <- par(mar = c(5.6, 1.2, 3.2, 6.4), mgp = c(3, 0.7, 0)); on.exit(par(op))
  plot(NA, xlim = xlim, ylim = ylim, xlab = "Age before present (thousand years)",
       ylab = "", xaxt = "n", yaxt = "n", log = logarg)
  if (xlog) {
    xt <- c(1000, 300, 100, 30, MH_LOG_AGE)
    axis(1, at = xt, labels = c("1000","300","100","30", paste0("~",MH_LOG_AGE,"*")))
  } else {
    axis(1, at = seq(800, 0, by = -200))
  }
  axis(4)
  mtext(ylab, side = 4, line = 3.6, cex = 1.0)
  title(main = main, cex.main = 1.1, font.main = 1)

  rect(800, ylim[1], 600, ylim[2], col = "#00000010", border = NA)
  div_x <- if (xlog) 470 else 560       # just right of the band, in open space
  div_y <- ylim[1] + diff(ylim) * 0.24  # low-middle band, clear of the mean lines
  text(div_x, div_y, "NT-MH divergence\n(~600-800 kya)",
       cex = 0.72, col = "grey30", adj = c(0.5, 0.5))
  abline(v = mh_x, col = "grey80", lty = 3)

  lab_x <- if (xlog) 2.4 else -52                      # right-margin mean label x
  gmean <- tapply(yv, spec$group, mean)
  for (g in c("NT","EH")) {
    abline(h = gmean[g], col = grp_col[g], lty = 2, lwd = 1.2)
    text(lab_x, gmean[g], sprintf("%.0f", gmean[g]), col = grp_col[g],
         cex = 0.72, adj = c(1, -0.4), xpd = NA)
  }
  abline(h = mh_y, col = grp_col["MH"], lty = 2, lwd = 1.2)

  for (i in seq_len(nrow(spec))) {
    if (!spec$undated[i] && spec$age_min_kya[i] != spec$age_max_kya[i]) {
      segments(spec$age_min_kya[i], yv[i], spec$age_max_kya[i], yv[i],
               col = spec$col[i], lwd = 1.8)
      segments(spec$age_min_kya[i], yv[i] - diff(ylim)*0.008,
               spec$age_min_kya[i], yv[i] + diff(ylim)*0.008, col = spec$col[i], lwd = 1.8)
      segments(spec$age_max_kya[i], yv[i] - diff(ylim)*0.008,
               spec$age_max_kya[i], yv[i] + diff(ylim)*0.008, col = spec$col[i], lwd = 1.8)
    }
    points(spec$age_plot[i], yv[i], pch = spec$pch[i],
           bg = ifelse(spec$undated[i], NA, paste0(spec$col[i], "CC")),
           col = spec$col[i], cex = spec$cex_budget[i], lwd = 1.8)
  }
  points(mh_x, mh, pch = MH_PCH, bg = paste0(grp_col["MH"], "CC"),
         col = grp_col["MH"], cex = mh_cex_budget, lwd = 1.8)

  # identity legend (shape = specimen, colour = group), with budget ratio
  leg_lab <- sprintf("%s [%s]  (budget %.2f)", spec$specimen, spec$sex, spec$budget_ratio_MH)
  legend("topleft",
         legend = c(leg_lab, "Modern human (present, 1.00)"),
         pch    = c(spec$pch, MH_PCH),
         pt.bg  = c(ifelse(spec$undated, NA, paste0(spec$col, "CC")), paste0(grp_col["MH"],"CC")),
         col    = c(spec$col, grp_col["MH"]),
         pt.cex = 1.3, bty = "n", cex = 0.72, y.intersp = 1.05)
  # size legend: decode marker area -> energy-budget ratio
  s_ref <- pretty(size_rng, n = 3); s_ref <- s_ref[s_ref >= size_rng[1] & s_ref <= size_rng[2]]
  legend("bottomleft", title = "Energy budget (ratio to MH)",
         legend = sprintf("%.2f", s_ref), pch = 21, col = "grey40", pt.bg = "grey80",
         pt.cex = cex_from_budget(s_ref), bty = "n", cex = 0.72,
         y.intersp = 1.5, x.intersp = 1.4)
  if (nzchar(note))
    mtext(note, side = 1, line = 4.1, cex = 0.60, adj = 0, col = "grey40")
  invisible(gmean)
}

# ---- Fig A: brain volume through time (true scale), point size = budget ----
png("figs/s4/volume_timeline.png", width = 2100, height = 1450, res = 220)
famous_scatter("total_cc",
               ylab = "Total brain volume (cerebral + cerebellar, cc)",
               main = "Fossil brain volume through time (point size = energy budget)",
               mh_y = mh_total_cc,
               note = "Colour = group; shape = specimen; point AREA = energy-budget ratio to MH (legend). Whiskers = dating range; hollow = undated (Forbes' Quarry 1). Dashed = group means.")
dev.off()

# ---- Fig A2: same, LOG10 age axis (spreads the fossils out) ----------------
png("figs/s4/volume_timeline_logage.png", width = 2100, height = 1450, res = 220)
famous_scatter("total_cc",
               ylab = "Total brain volume (cerebral + cerebellar, cc)",
               main = "Fossil brain volume through time (log age; point size = energy budget)",
               mh_y = mh_total_cc, xlog = TRUE,
               note = "Log10 age axis. *Modern human plotted at a nominal ~5 kya (age 0 undefined on log scale). Point AREA = energy-budget ratio to MH; colour = group; shape = specimen. Markers sit at the reported mean date; on a log axis a linearly-symmetric date range is not symmetric, so points appear slightly toward the older (left) end of their whisker.")
dev.off()

# ---- Fig A2-flip: glucose utilization through time; point size = brain volume
# Second version of Fig A2 with the two data encodings SWAPPED:
#   Y axis     = glucose utilization (ratio to MH)   [was the point size]
#   point AREA = total brain volume (cc)             [was the y axis]
# "Energy budget" is relabelled to the more specific "glucose utilization"
# throughout this figure (glucose utilization is what the budget measures).
# marker AREA mapped LINEARLY across the observed brain-volume range (as before)
VOL_SIZE_CEX <- c(1.5, 4.6)
vol_vals <- c(spec$total_cc, mh_total_cc)
vol_rng  <- range(vol_vals)
cex_from_vol <- function(v)                         # brain volume (cc) -> cex
  sqrt(VOL_SIZE_CEX[1]^2 + (VOL_SIZE_CEX[2]^2 - VOL_SIZE_CEX[1]^2) *
         (v - vol_rng[1]) / diff(vol_rng))
spec$cex_vol <- cex_from_vol(spec$total_cc)
mh_cex_vol   <- cex_from_vol(mh_total_cc)

famous_scatter_flip <- function(main, note = "", xlog = FALSE) {
  yv <- spec$budget_ratio_MH; mh <- 1.0            # glucose utilization; MH = 1.00
  yr <- range(c(yv, mh)); pad <- diff(yr) * 0.06
  ylim <- c(yr[1] - pad, yr[2] + pad)
  mh_x <- if (xlog) MH_LOG_AGE else 0
  xlim <- if (xlog) c(1000, 3) else c(850, -55)
  logarg <- if (xlog) "x" else ""

  op <- par(mar = c(5.6, 1.2, 3.2, 6.4), mgp = c(3, 0.7, 0)); on.exit(par(op))
  plot(NA, xlim = xlim, ylim = ylim, xlab = "Age before present (thousand years)",
       ylab = "", xaxt = "n", yaxt = "n", log = logarg)
  if (xlog) {
    xt <- c(1000, 300, 100, 30, MH_LOG_AGE)
    axis(1, at = xt, labels = c("1000","300","100","30", paste0("~",MH_LOG_AGE,"*")))
  } else {
    axis(1, at = seq(800, 0, by = -200))
  }
  axis(4)
  mtext("Glucose utilization (ratio to MH)", side = 4, line = 3.6, cex = 1.0)
  title(main = main, cex.main = 1.1, font.main = 1)

  rect(800, ylim[1], 600, ylim[2], col = "#00000010", border = NA)
  div_x <- if (xlog) 470 else 560
  div_y <- ylim[1] + diff(ylim) * 0.24
  text(div_x, div_y, "NT-MH divergence\n(~600-800 kya)",
       cex = 0.72, col = "grey30", adj = c(0.5, 0.5))
  abline(v = mh_x, col = "grey80", lty = 3)

  lab_x <- if (xlog) 2.4 else -52
  gmean <- tapply(yv, spec$group, mean)
  for (g in c("NT","EH")) {
    abline(h = gmean[g], col = grp_col[g], lty = 2, lwd = 1.2)
    text(lab_x, gmean[g], sprintf("%.2f", gmean[g]), col = grp_col[g],
         cex = 0.72, adj = c(1, -0.4), xpd = NA)
  }
  abline(h = mh, col = grp_col["MH"], lty = 2, lwd = 1.2)

  for (i in seq_len(nrow(spec))) {
    if (!spec$undated[i] && spec$age_min_kya[i] != spec$age_max_kya[i]) {
      segments(spec$age_min_kya[i], yv[i], spec$age_max_kya[i], yv[i],
               col = spec$col[i], lwd = 1.8)
      segments(spec$age_min_kya[i], yv[i] - diff(ylim)*0.008,
               spec$age_min_kya[i], yv[i] + diff(ylim)*0.008, col = spec$col[i], lwd = 1.8)
      segments(spec$age_max_kya[i], yv[i] - diff(ylim)*0.008,
               spec$age_max_kya[i], yv[i] + diff(ylim)*0.008, col = spec$col[i], lwd = 1.8)
    }
    points(spec$age_plot[i], yv[i], pch = spec$pch[i],
           bg = ifelse(spec$undated[i], NA, paste0(spec$col[i], "CC")),
           col = spec$col[i], cex = spec$cex_vol[i], lwd = 1.8)
  }
  points(mh_x, mh, pch = MH_PCH, bg = paste0(grp_col["MH"], "CC"),
         col = grp_col["MH"], cex = mh_cex_vol, lwd = 1.8)

  # identity legend (shape = specimen, colour = group), with brain volume
  leg_lab <- sprintf("%s [%s]  (%.0f cc)", spec$specimen, spec$sex, spec$total_cc)
  legend("topleft",
         legend = c(leg_lab, sprintf("Modern human (present, %.0f cc)", mh_total_cc)),
         pch    = c(spec$pch, MH_PCH),
         pt.bg  = c(ifelse(spec$undated, NA, paste0(spec$col, "CC")), paste0(grp_col["MH"],"CC")),
         col    = c(spec$col, grp_col["MH"]),
         pt.cex = 1.3, bty = "n", cex = 0.72, y.intersp = 1.05)
  # size legend: decode marker area -> total brain volume
  v_ref <- pretty(vol_rng, n = 3); v_ref <- v_ref[v_ref >= vol_rng[1] & v_ref <= vol_rng[2]]
  legend("bottomleft", title = "Total brain volume (cc)",
         legend = sprintf("%.0f", v_ref), pch = 21, col = "grey40", pt.bg = "grey80",
         pt.cex = cex_from_vol(v_ref), bty = "n", cex = 0.72,
         y.intersp = 1.5, x.intersp = 1.4)
  if (nzchar(note))
    mtext(note, side = 1, line = 4.1, cex = 0.60, adj = 0, col = "grey40")
  invisible(gmean)
}

save_png_pdf("figs/s4/volume_timeline_logage_flipped",
  function() famous_scatter_flip(
    main = "Fossil glucose utilization through time (log age; point size = brain volume)",
    xlog = TRUE,
    note = "Log10 age axis. *Modern human plotted at a nominal ~5 kya (age 0 undefined on log scale). Y = glucose utilization (ratio to MH); point AREA = total brain volume (cc); colour = group; shape = specimen. Markers sit at the reported mean date; on a log axis a linearly-symmetric date range is not symmetric, so points appear slightly toward the older (left) end of their whisker."),
  width = 2100, height = 1450, res = 220)

# ---- species-average helper: means with 95% CI ---------------------------
# CI (not SD) answers "how well is each group mean pinned down". Fossil n = 4
# -> wide intervals; modern-human n = 1185 -> tight.
sp_age <- c(NT = mean(spec$age_plot[spec$group == "NT"]),
            EH = mean(spec$age_plot[spec$group == "EH"]), MH = 0)
sp <- species_budgets
sp$age <- sp_age[sp$group]

species_ci_plot <- function(ytick = NULL,
                            tick_fmt = function(v) sprintf("%.2f", v),
                            main = "Species-average brain energy budget", note = "") {
  lo <- sp$budget_ratio_MH - sp$budget_ratio_MH_ci95
  hi <- sp$budget_ratio_MH + sp$budget_ratio_MH_ci95
  ctr <- sp$budget_ratio_MH
  yr <- range(c(lo, hi, 1)); pad <- diff(yr) * 0.10
  ylim <- c(yr[1] - pad, yr[2] + pad); xlim <- c(135, -15)
  op <- par(mar = c(5.4, 5.2, 3.2, 1.4), mgp = c(3, 0.7, 0)); on.exit(par(op))
  plot(NA, xlim = xlim, ylim = ylim, xlab = "Mean age before present (thousand years)",
       ylab = "Energy budget (ratio to modern-human mean)", yaxt = "n")
  if (is.null(ytick)) axis(2) else axis(2, at = ytick, labels = tick_fmt(ytick))
  title(main = main, cex.main = 1.05, font.main = 1)
  abline(h = 1, col = "grey70", lty = 2); abline(v = 0, col = "grey85", lty = 3)
  arrows(sp$age, lo, sp$age, hi, code = 3, angle = 90, length = 0.06,
         col = grp_col[sp$group], lwd = 2)
  points(sp$age, ctr, pch = c(19,19,18)[match(sp$group, c("NT","EH","MH"))],
         col = grp_col[sp$group], cex = 1.9)
  text(sp$age, hi, sprintf("%s\n%.3f (95%% CI +/- %.3f, n=%d)",
                           sp$group, sp$budget_ratio_MH, sp$budget_ratio_MH_ci95, sp$n),
       pos = 3, offset = 0.5, cex = 0.62, col = "grey20", xpd = NA)
  legend("topright", legend = c("Neanderthal (NT)","Early modern human (EH)","Modern human (MH)"),
         col = grp_col[c("NT","EH","MH")], pch = c(19,19,18), pt.cex = 1.3, bty = "n", cex = 0.8)
  if (nzchar(note)) mtext(note, side = 1, line = 3.6, cex = 0.64, adj = 0, col = "grey40")
}

# ---- Fig C: species means +/- 95% CI (true scale) -------------------------
png("figs/s4/species_budget_ci.png", width = 2000, height = 1300, res = 220)
species_ci_plot(ytick = seq(0.85, 1.20, 0.05),
                main = "Species-average brain energy budget (95% CI, true scale)",
                note = "Error bars = 95% CI on the mean (t-based, from group volume s.d. and n). Wide fossil CIs reflect n=4.")
dev.off()

# =============================================================================
# 8. Does regional metabolic rate predict evolutionary change in region SHAPE?
# =============================================================================
# DIFFERENT QUESTION from the budget sections above. Everything above holds
# rCMRGlc FIXED at the modern-human value and asks how endocast SIZE scales the
# regional energy cost. Here we instead ask whether the modern-human regional
# metabolic INTENSITY itself is associated with how much each region's SHARE of
# the brain has changed over hominin evolution.
#
#   x = rCMRGlc                  modern-human regional glucose rate
#                                (Heiss et al. 2004; umol glucose/100 g/min)
#   y = Fig 3A relative volume   region's ICV-size-adjusted volume in a fossil
#                                group RELATIVE to modern human (Kochiyama et al.
#                                2018, Fig 3A; MH = 1). This is the "deviation
#                                from MH" metric: y>1 = larger share than MH,
#                                y<1 = smaller share than MH. y is size-free, so
#                                it isolates SHAPE (regional reallocation), not
#                                overall brain enlargement.
#
# Two contrasts, both against modern human:  1. Neanderthal (NT/MH),
#                                            2. early modern human (EH/MH).
# A positive correlation would mean metabolically expensive regions expanded
# their share toward the modern-human condition (or contracted in the fossil);
# a null result means metabolic intensity does NOT structure which regions were
# evolutionarily reallocated. n = 6 matched regions (same regions that carry
# both a Heiss rate and a Kochiyama volume), so this is a low-power descriptive
# correlation, reported with its p-value and CI, not a significance claim.
#
# OUTPUTS
#   tables/s4/rcmrglc_vs_fig3a_deviation.csv     per-region + fitted stats
#   figs/s4/rcmrglc_vs_fig3a_deviation.png       scatter + OLS fit + 95% CI band

corr_df <- data.frame(
  Heiss_region = agg$Heiss_region,
  rCMRGlc      = agg$rCMRGlc,
  NT_rel       = agg$NT_rel,
  EH_rel       = agg$EH_rel,
  stringsAsFactors = FALSE)

# --- correlations + OLS fits (Neanderthal and early modern human vs MH) ------
corr_stats <- do.call(rbind, lapply(c("NT","EH"), function(sp) {
  y  <- corr_df[[paste0(sp,"_rel")]]; x <- corr_df$rCMRGlc
  ct <- cor.test(x, y)                          # Pearson (+ 95% CI on r)
  sr <- suppressWarnings(cor.test(x, y, method = "spearman"))
  fit <- lm(y ~ x); sm <- summary(fit)
  data.frame(
    contrast    = paste0(sp, " / MH"),
    n           = length(x),
    pearson_r   = unname(ct$estimate),
    pearson_p   = ct$p.value,
    r_CI_low    = ct$conf.int[1],
    r_CI_high   = ct$conf.int[2],
    spearman_rho= unname(sr$estimate),
    spearman_p  = sr$p.value,
    slope       = unname(coef(fit)[2]),
    slope_se    = sm$coefficients[2,2],
    intercept   = unname(coef(fit)[1]),
    r_squared   = sm$r.squared,
    stringsAsFactors = FALSE)
}))
corr_stats[ , sapply(corr_stats, is.numeric)] <-
  lapply(corr_stats[ , sapply(corr_stats, is.numeric)], round, 4)

write.csv(corr_df,    "tables/s4/rcmrglc_vs_fig3a_deviation.csv", row.names = FALSE)
write.csv(corr_stats, "tables/s4/rcmrglc_vs_fig3a_deviation_stats.csv", row.names = FALSE)

cat("\n=== rCMRGlc vs Fig-3A deviation-from-MH (per matched region) ===\n")
print(corr_df, row.names = FALSE)
cat("\n=== Correlation / regression summary ===\n")
print(corr_stats, row.names = FALSE)

# --- figure OPTION 1: base-R, Pearson + OLS fit + 95% CI ribbon -------------
# Descriptive linear view: shows the OLS fit and its 95% CI band, annotated with
# the Pearson r, slope and n. Written to *_fit.png. (The fit is not significant;
# this panel is kept only as an alternative descriptive view -- see OPTION 2 for
# the significance-neutral rank-correlation version with no line.)
corr_panel <- function(sp, col) {
  y <- corr_df[[paste0(sp,"_rel")]]; x <- corr_df$rCMRGlc
  ct <- cor.test(x, y); fit <- lm(y ~ x)
  xr <- seq(min(x) - 0.5, max(x) + 0.5, length.out = 100)
  pr <- predict(fit, newdata = data.frame(x = xr), interval = "confidence", level = 0.95)
  yr <- range(c(y, pr)); pad <- diff(yr) * 0.12
  plot(NA, xlim = range(x) + c(-0.6, 0.6), ylim = c(yr[1] - pad, yr[2] + pad),
       xlab = "Modern-human rCMRGlc  (umol glucose / 100 g / min)",
       ylab = sprintf("Fig 3A relative volume  (%s / MH)", sp), las = 1)
  abline(h = 1, col = "grey70", lty = 2)                         # 1.0 = no difference vs MH
  polygon(c(xr, rev(xr)), c(pr[,"lwr"], rev(pr[,"upr"])),
          col = paste0(col, "33"), border = NA)                  # 95% CI band
  abline(fit, col = col, lwd = 2)
  points(x, y, pch = 21, bg = paste0(col, "CC"), col = col, cex = 1.7, lwd = 1.6)
  text(x, y, corr_df$Heiss_region, pos = 3, offset = 0.6, cex = 0.62, col = "grey25", xpd = NA)
  title(main = sprintf("%s vs modern human",
                       c(NT = "Neanderthal", EH = "Early modern human")[sp]),
        cex.main = 1.05, font.main = 1)
  legend("topleft", bty = "n", cex = 0.74,
         legend = c(sprintf("Pearson r = %.2f  (p = %.2f)", ct$estimate, ct$p.value),
                    sprintf("slope = %.4f per unit rCMRGlc", coef(fit)[2]),
                    sprintf("n = %d matched regions", nrow(corr_df))))
}

png("figs/s4/rcmrglc_vs_fig3a_deviation_fit.png", width = 2500, height = 1150, res = 220)
op <- par(mfrow = c(1,2), mar = c(4.8, 4.8, 3.0, 1.2), mgp = c(2.9, 0.7, 0),
          oma = c(2.2, 0, 2.4, 0))
corr_panel("NT", grp_col["NT"]); corr_panel("EH", grp_col["EH"])
mtext("Does modern-human regional metabolic rate predict evolutionary change in regional volume?",
      outer = TRUE, cex = 1.0, font = 2, line = 0.5)
mtext(paste0("x = MH regional glucose rate (Heiss et al. 2004);  y = ICV-size-adjusted relative volume vs MH ",
             "(Kochiyama et al. 2018, Fig 3A).  Dashed line = 1.0 (no shape difference).  Shaded band = 95% CI of the OLS fit."),
      side = 1, outer = TRUE, cex = 0.58, col = "grey40", line = 0.6)
par(op); dev.off()

# --- figure OPTION 2: two panels (NT, EH) in the s3 house style -------------
# Matches s3_predicValuesPGLS_..._PATCHED.R (p_scatter): theme_bw with no grid
# and a black panel border, black filled points with structure labels, and the
# Spearman rho / p / n annotation in navy at top-left. NO fitted line and NO
# error bars -- the correlation is not significant, so only the rank statistic
# and the raw scatter are shown (a fitted line would over-state the pattern).
library(ggplot2)
library(patchwork)

x_lab_gg <- "rCMRGlc (\u00b5mol/100 g/min)"      # micro sign, matching s3
theme_paper <- theme_bw(base_size = 16) +
  theme(panel.grid   = element_blank(),
        axis.title   = element_text(size = 18), axis.text = element_text(size = 15),
        plot.title   = element_text(size = 18, face = "bold"),
        plot.subtitle= element_text(size = 15),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        plot.margin  = margin(8, 8, 8, 8))

fmt_p <- function(p) format.pval(p, digits = 2, eps = 1e-4)

s4_corr_scatter <- function(sp, sp_full) {
  d <- data.frame(
    rCMRGlc   = corr_df$rCMRGlc,
    y         = corr_df[[paste0(sp, "_rel")]],
    Structure = corr_df$Heiss_region
  )
  
  d <- d[complete.cases(d), ]
  
  ct <- suppressWarnings(
    cor.test(d$rCMRGlc, d$y, method = "spearman", exact = FALSE)
  )
  
  cor_stats <- paste0(
    "Spearman \u03c1 = ", sprintf("%.3f", unname(ct$estimate)),
    "\np = ", fmt_p(ct$p.value),
    "\nn = ", nrow(d)
  )
  
  xr <- range(d$rCMRGlc)
  yr <- range(d$y)
  
  ggplot(d, aes(rCMRGlc, y)) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      color = "gray60"
    ) +
    geom_point(pch = 16, size = 2.2) +
    ggrepel::geom_text_repel(
      aes(label = Structure),
      size = 2.8,
      direction = "both",
      box.padding = 0.25,
      point.padding = 0.15,
      min.segment.length = 0,
      max.overlaps = Inf,
      seed = 1
    ) +
    annotate(
      "text",
      x = xr[1] - diff(xr) * 0.02,
      y = yr[2] + diff(yr) * 0.30,
      label = cor_stats,
      hjust = 0,
      vjust = 1,
      size = 3.2,
      color = "navy"
    ) +
    coord_cartesian(
      xlim = c(
        xr[1] - diff(xr) * 0.05,
        xr[2] + diff(xr) * 0.55
      ),
      ylim = c(
        yr[1] - diff(yr) * 0.10,
        yr[2] + diff(yr) * 0.38
      ),
      clip = "on"
    ) +
    labs(
      title = sprintf(
        "Scatterplot: rCMRGlc and %s volume difference",
        sp
      ),
      subtitle = sprintf(
        "%s relative volume vs modern human (Fig 3A)",
        sp_full
      ),
      x = x_lab_gg,
      y = sprintf(
        "Fig 3A relative volume (%s / MH)",
        sp
      )
    ) +
    theme_paper
}

p_corr_nt <- s4_corr_scatter("NT", "Neanderthal")
p_corr_eh <- s4_corr_scatter("EH", "Early modern human")
p_corr    <- p_corr_nt + p_corr_eh
# type = "cairo" renders the Greek rho (\u03c1) correctly regardless of the
# session locale (base png() substitutes a dot under a non-UTF-8 locale, e.g.
# when the script is run head-less via Rscript). cairo_pdf() does the same for pdf.
ggsave("figs/s4/rcmrglc_vs_fig3a_deviation.png", p_corr, width = 15, height = 6.2,
       dpi = 220, type = "cairo")
ggsave("figs/s4/rcmrglc_vs_fig3a_deviation.pdf", p_corr, width = 15, height = 6.2,
       device = cairo_pdf)
