# =====================================================================
# s4b_arterial_canal.R
# ---------------------------------------------------------------------
# Brain glucose metabolism of fossil hominins from ARTERIAL CANAL size,
# combining the rationales of Boyer & Harrington (2018) and Seymour et al.
# (2016/2017), for comparison with the VOLUME-based (Kochiyama shape x Heiss
# regional rCMRGlc) budgets produced by s4a_endocranial.R.
#
# ---------------------------------------------------------------------
# WHY TWO METHODS (and what each assumes)
# ---------------------------------------------------------------------
# Fossil crania reliably preserve only the INTERNAL CAROTID canal/foramen.
# The vertebral-artery contribution (Boyer's doubled transverse-foramen area,
# DTFA, measured on cervical vertebrae) is essentially never available. So per
# fossil we have: internal-carotid foramen radius (Seymour) + endocranial
# volume. Two published rationales turn that into brain metabolism:
#
#   [A] SEYMOUR flow-anchored.  Seymour et al. predict total internal-carotid
#       blood-flow rate Q (cm^3 s^-1) from foramen radius (Hagen-Poiseuille,
#       Q ~ r^3) and argue perfusion is proportional to the tissue's metabolic
#       rate. Hence relative flow ~ relative metabolism. We anchor the modern
#       human at the same whole-brain glucose value s4a uses (Clarke & Sokoloff
#       1994 / Boyer Table 2 ~= 429 umol/min) and scale each fossil by its flow
#       ratio to the modern-human reference specimen:
#            BGU_fossil = BGU_modH * (Q_fossil / Q_modH)
#       This uses ONLY the fossil data that actually exists (carotid), needs no
#       total-ACA reconstruction, and lands in the same units as s4a.
#
#   [B] BOYER regression.  Boyer & Harrington calibrated whole-Brain Glucose
#       Utilization (BGU) on TOTAL arterial canal area (ACA) + ECV across 7
#       extant taxa with directly measured BGU (their Table 2), and showed ACA
#       predicts metabolism better than ECV or neuron count. They never applied
#       it to fossils, precisely because fossils lack the vertebral canal. To
#       use their (better-calibrated) equation we must reconstruct a TOTAL ACA
#       for each fossil from its carotid measurement. There is no single agreed
#       way to do this, so we compute the Boyer branch TWO ways and report the
#       spread as a sensitivity band:
#         B1 (primary) carotid-ratio scaling. Assuming the carotid share of
#            total ACA is conserved within hominins, scale the modern-human
#            total ACA (Boyer Table 1 Homo = 159.81 mm^2) by the fossil's
#            carotid CSA relative to the modern-human carotid CSA, both taken
#            from Seymour's OWN foramen radii (measurement-consistent, avoiding
#            the foramen-vs-canal landmark mismatch):
#                 ACA_fossil = ACA_Homo * (r_fossil / r_modH)^2
#         B2 (bound/check) ACA predicted from ECV via Boyer's own ACA~ECV
#            allometry. NOTE: hominins sit BELOW the general euarchontan
#            ACA~ECV line (small carotids for their brain size), so this
#            over-predicts ACA and therefore BGU for hominins -- it is reported
#            as an upper bound, not an independent estimate.
#       Then  BGU = exp(b0 + b_ACA*ln(ACA) + b_ECV*ln(ECV)) * smearing
#       with coefficients refit here (log-log OLS) from Boyer Table 2 and a
#       Duan (1983) smearing factor correcting log-back-transform bias.
#
# ---------------------------------------------------------------------
# VOLUME CONVENTIONS  --  three DIFFERENT quantities; do not conflate
# ---------------------------------------------------------------------
#   * Cranial / endocranial CAPACITY = the whole braincase cavity (brain
#     tissue + CSF + ventricles + meninges + vessels). This is what Seymour's
#     "Brain_volume_cm3" and Boyer's "ECV" actually are (verified: Seymour
#     La Chapelle 1625 ~= Hawks & Wolpoff cranial capacity 1626; Skhul 5
#     1520 = 1520). LARGE.
#   * Brain-TISSUE volume (GM+WM) = metabolically active parenchyma, from
#     Kochiyama's MRI-based endocast reconstruction. This is what s4a uses
#     (cerebrum + cerebellum). SMALLER: e.g. Skhul 5 tissue 1199 cc vs
#     capacity 1520 cc; modern-human tissue ~1246 cc vs capacity ~1450-1493.
#   The Boyer branch here is INTERNALLY consistent -- Boyer calibrated BGU on
#   endocranial capacity (Table 2 Homo ECV 1422) and we feed fossils Seymour's
#   endocranial capacity, so no tissue<->capacity conversion is done or needed.
#   The Seymour branch uses flow only, no volume. But when COMPARING against
#   s4a, the volumes are different conventions: output columns are labelled by
#   provenance (endocranial_cc_arterial vs brain_tissue_cc_s4a) and are NOT the
#   same measurement.
#
# ---------------------------------------------------------------------
# METABOLIC SCOPE  --  whole brain vs partial regional sum
# ---------------------------------------------------------------------
#   Boyer/Seymour BGU is WHOLE-BRAIN glucose (~429 umol/min in modern humans).
#   s4a's budget is the SUM OF 6 Heiss regions only (frontal, parietal, temporal,
#   occipital lobes + cerebellar cortex + vermis) = 328.51 umol/min for the
#   modern human = 76.7% of whole brain (it omits deep grey, brainstem,
#   diencephalon, etc.). So the ABSOLUTE umol/min values are NOT directly
#   comparable across the two families of method.
#
#   -> PRIMARY cross-method comparison is therefore the RATIO to the modern
#      human within each method (dimensionless), which cancels both the
#      capacity-vs-tissue and whole-brain-vs-partial offsets with no extra
#      assumption. We additionally emit a scope-adjusted arterial absolute
#      (BGU x 0.767) on s4a's cortical+cerebellar footing, flagged as assuming
#      the covered-region fraction is constant across taxa.
#
# CROSS-METHOD VALIDATION.  La Chapelle-aux-Saints, Forbes'/Gibraltar and
# Skhul 5 appear in BOTH Seymour and s4a, giving direct convergence checks
# (compare on RATIO, not absolute umol/min).
#
# ---------------------------------------------------------------------
# INPUTS  (all relative to repo root)
#   data_raw/Boyer_Harrington_2018_Table2.csv   7-taxon BGU calibration
#   data_raw/Boyer_Harrington_2018_Table1.csv   extant euarchontans: DPA/ACA/ECV
#   data_raw/Seymour_etal_2017_TableS1.csv       30 fossil hominin specimens
#   data_intermediate/s4a_specimen_budgets.csv    s4a volume-based per-specimen
#   tables/s4a/species_absolute_budgets.csv       s4a volume-based group means
#
# OUTPUTS
#   data_intermediate/s4b_boyer_calibration.csv   fitted coefficients + fit stats
#   data_intermediate/s4b_fossil_estimates.csv    per-specimen, all methods (wide)
#   data_intermediate/s4b_fossil_estimates_long.csv  tidy long form (one row/est.)
#   tables/s4b/s4b_specimen_crosswalk.csv           Seymour<->s4a name matching
#   tables/s4b/s4b_method_comparison_overlap.csv    overlap specimens across methods
#   tables/s4b/s4b_group_means.csv                  NT / EH / MH group summaries
#   tables/s4b/s4b_sapiens_grade_means.csv          H. sapiens early/recent/modern split
#
# All glucose values are umol glucose / min (same as s4a). No figures are
# produced (per request) -- these tables feed downstream plotting.
#
# References
#   Boyer DM, Harrington AR (2018). Scaling of bony canals for encephalic
#     vessels in euarchontans. J. Hum. Evol. 114, 85-101.
#   Seymour RS, Bosiocic V, Snelling EP (2016; Correction 2017). Fossil skulls
#     reveal that blood flow rate to the brain increased faster than brain
#     volume during human evolution. R. Soc. Open Sci. 3:160305 / 4:170846.
#   Karbowski J (2007). Global and regional brain metabolic scaling. BMC Biol.
#   Clarke DD, Sokoloff L (1994). Circulation and energy metabolism of the brain.
#   Duan N (1983). Smearing estimate: a nonparametric retransformation method.
# =====================================================================

suppressPackageStartupMessages({
  library(stats)
})

## ---- paths: self-contained (Rscript or RStudio; full repo or lone folder) ----
setwd(local({ d <- normalizePath(getwd())
              while (!file.exists(file.path(d, ".git")) &&
                     dirname(d) != d) d <- dirname(d); d }))  # repo root

dir.create("tables/s4b",         showWarnings = FALSE, recursive = TRUE)
dir.create("data_intermediate", showWarnings = FALSE, recursive = TRUE)

rd <- function(p) read.csv(p, stringsAsFactors = FALSE, check.names = FALSE,
                           comment.char = "#", encoding = "UTF-8")

# UTF-8-safe writer (specimen names carry diacritics, e.g. "Mladec 1")
write_csv_utf8 <- function(df, path) {
  con <- file(path, open = "wb", encoding = "UTF-8")
  on.exit(close(con))
  writeLines(paste(names(df), collapse = ","), con, useBytes = TRUE)
  apply(df, 1, function(r)
    writeLines(paste(ifelse(is.na(r), "", r), collapse = ","), con, useBytes = TRUE))
  invisible(path)
}

## =====================================================================
## 1. Boyer & Harrington (2018) BGU calibration  (Table 2, 7 extant taxa)
## =====================================================================
cal <- rd("data_raw/Boyer_Harrington_2018_Table2.csv")
cal$lnBGU <- log(cal$BGU_umol_min)
cal$lnACA <- log(cal$ACA_mm2)
cal$lnECV <- log(cal$ECV_cc)

fit_full <- lm(lnBGU ~ lnACA + lnECV, data = cal)   # preferred: ACA + ECV
fit_aca  <- lm(lnBGU ~ lnACA,          data = cal)   # ACA only
fit_ecv  <- lm(lnBGU ~ lnECV,          data = cal)   # ECV only

# Duan (1983) smearing factor: mean(exp(residuals)) corrects the bias from
# predicting on the log scale and exponentiating back to umol/min.
smear <- function(m) mean(exp(residuals(m)))
smear_full <- smear(fit_full)

cat("=== Boyer BGU calibration (log-log OLS, n =", nrow(cal), ") ===\n")
cat(sprintf("  ln(BGU) = %.4f + %.4f*ln(ACA) + %.4f*ln(ECV)   R2=%.4f  smear=%.4f\n",
            coef(fit_full)[1], coef(fit_full)[2], coef(fit_full)[3],
            summary(fit_full)$r.squared, smear_full))
cat(sprintf("  ACA-only: ln(BGU)=%.4f + %.4f*ln(ACA)   R2=%.4f\n",
            coef(fit_aca)[1], coef(fit_aca)[2], summary(fit_aca)$r.squared))
cat(sprintf("  ECV-only: ln(BGU)=%.4f + %.4f*ln(ECV)   R2=%.4f\n",
            coef(fit_ecv)[1], coef(fit_ecv)[2], summary(fit_ecv)$r.squared))
# calibration self-check: predicted Homo should match observed 428.55
homo_pred <- exp(predict(fit_full, data.frame(lnACA = log(159.81),
                                              lnECV = log(1422.33)))) * smear_full
cat(sprintf("  self-check: predicted Homo BGU = %.1f (observed 428.55)\n\n", homo_pred))

# save calibration
calib_out <- data.frame(
  model      = c("ACA+ECV", "ACA_only", "ECV_only"),
  intercept  = c(coef(fit_full)[1], coef(fit_aca)[1], coef(fit_ecv)[1]),
  b_lnACA    = c(coef(fit_full)[2], coef(fit_aca)[2], NA),
  b_lnECV    = c(coef(fit_full)[3], NA,               coef(fit_ecv)[2]),
  r_squared  = c(summary(fit_full)$r.squared, summary(fit_aca)$r.squared,
                 summary(fit_ecv)$r.squared),
  smearing   = c(smear(fit_full), smear(fit_aca), smear(fit_ecv)),
  n          = nrow(cal),
  stringsAsFactors = FALSE)
calib_out[ , sapply(calib_out, is.numeric)] <-
  round(calib_out[ , sapply(calib_out, is.numeric)], 4)
write.csv(calib_out, "data_intermediate/s4b_boyer_calibration.csv", row.names = FALSE)

## =====================================================================
## 2. ACA ~ ECV allometry (for the B2 ECV-predicted-ACA bound)
##    Fit on catarrhines (Hominoidea + Cercopithecoidea) -- the appropriate
##    phylogenetic bracket for hominins -- and also on all euarchontans for
##    context. Report both; use the catarrhine fit for the bound.
## =====================================================================
t1 <- rd("data_raw/Boyer_Harrington_2018_Table1.csv")
t1$ACA_mm2 <- as.numeric(t1$ACA_mm2)
t1$ECV_cc  <- as.numeric(t1$ECV_cc)
t1$DPA_mm2 <- as.numeric(t1$DPA_mm2)
t1 <- t1[is.finite(t1$ACA_mm2) & is.finite(t1$ECV_cc), ]

catarrhine <- t1[t1$Taxonomic_group %in% c("Hominoidea", "Cercopithecoidea"), ]
fit_aca_ecv_cat <- lm(log(ACA_mm2) ~ log(ECV_cc), data = catarrhine)
fit_aca_ecv_all <- lm(log(ACA_mm2) ~ log(ECV_cc), data = t1)
smear_cat <- mean(exp(residuals(fit_aca_ecv_cat)))

cat("=== ACA ~ ECV allometry (for ECV-predicted-ACA bound) ===\n")
cat(sprintf("  catarrhine (n=%d): ln(ACA)=%.4f + %.4f*ln(ECV)  R2=%.3f  smear=%.4f\n",
            nrow(catarrhine), coef(fit_aca_ecv_cat)[1], coef(fit_aca_ecv_cat)[2],
            summary(fit_aca_ecv_cat)$r.squared, smear_cat))
cat(sprintf("  all euarchontan (n=%d): ln(ACA)=%.4f + %.4f*ln(ECV)  R2=%.3f\n\n",
            nrow(t1), coef(fit_aca_ecv_all)[1], coef(fit_aca_ecv_all)[2],
            summary(fit_aca_ecv_all)$r.squared))

## Extant hominoid carotid share of ACA (DPA:ACA) -- provenance / uncertainty of
## the "carotid fraction conserved" assumption behind the B1 scaling.
hom <- t1[t1$Taxonomic_group == "Hominoidea" & is.finite(t1$DPA_mm2), ]
hom$dpa_aca <- hom$DPA_mm2 / hom$ACA_mm2
cat("=== extant hominoid DPA:ACA (carotid share of ACA) ===\n")
print(data.frame(Species = hom$Species, dpa_aca = round(hom$dpa_aca, 3)),
      row.names = FALSE)
cat(sprintf("  mean = %.3f  sd = %.3f  CV = %.1f%%\n\n",
            mean(hom$dpa_aca), sd(hom$dpa_aca),
            100 * sd(hom$dpa_aca) / mean(hom$dpa_aca)))

## =====================================================================
## 3. Constants / modern-human anchor
## =====================================================================
BGU_MODH  <- 428.55   # umol/min, WHOLE-BRAIN glucose, Clarke & Sokoloff 1994 / Boyer Table 2 Homo
ACA_HOMO  <- 159.81   # mm^2, Boyer Table 1 total ACA for Homo sapiens
# Modern-human reference from Seymour: specimen "AS8078" = "Mean modern H. sapiens"
Q_MODH    <- 7.34     # cm^3/s, total Q_ICA of AS8078
R_MODH    <- 0.302    # cm, internal-carotid foramen radius of AS8078
ECV_MODH  <- 1493     # cc, endocranial CAPACITY of AS8078 (NOT brain tissue)

# Metabolic-scope bridge to s4a: s4a's 6-region sum captures this fraction of
# whole-brain glucose in the modern human (328.51 / 428.55). Used only to put
# the whole-brain arterial BGU onto s4a's cortical+cerebellar footing.
S4A_MODH_REGIONAL <- 328.51
S4A_SCOPE_FRAC    <- S4A_MODH_REGIONAL / BGU_MODH   # ~0.767

## =====================================================================
## 4. Fossil specimens (Seymour 2017 Table S1)
## =====================================================================
sey <- rd("data_raw/Seymour_etal_2017_TableS1.csv")
sey$Foramen_radius_cm <- as.numeric(sey$Foramen_radius_cm)
sey$Total_QICA_cm3_s  <- as.numeric(sey$Total_QICA_cm3_s)
sey$Brain_volume_cm3  <- as.numeric(sey$Brain_volume_cm3)

## -- [A] Seymour flow-anchored -------------------------------------------------
sey$BGU_seymour <- BGU_MODH * (sey$Total_QICA_cm3_s / Q_MODH)

## -- [B1] Boyer, carotid-ratio-scaled ACA -------------------------------------
sey$ACA_scaled  <- ACA_HOMO * (sey$Foramen_radius_cm / R_MODH)^2
## -- [B2] Boyer, ECV-predicted ACA (catarrhine allometry) ---------------------
sey$ACA_ecvpred <- exp(predict(fit_aca_ecv_cat,
                        data.frame(ECV_cc = sey$Brain_volume_cm3))) * smear_cat

# Boyer prediction helper: returns point estimate + 95% prediction interval
# (propagates calibration uncertainty via lm predict se; smearing on point est).
boyer_pred <- function(aca, ecv) {
  nd <- data.frame(lnACA = log(aca), lnECV = log(ecv))
  p  <- predict(fit_full, nd, se.fit = TRUE)
  tcrit <- qt(0.975, df = fit_full$df.residual)
  data.frame(fit = exp(p$fit) * smear_full,
             lwr = exp(p$fit - tcrit * p$se.fit) * smear_full,
             upr = exp(p$fit + tcrit * p$se.fit) * smear_full)
}
b1 <- boyer_pred(sey$ACA_scaled,  sey$Brain_volume_cm3)
b2 <- boyer_pred(sey$ACA_ecvpred, sey$Brain_volume_cm3)
sey$BGU_boyer_scaled     <- b1$fit
sey$BGU_boyer_scaled_lwr <- b1$lwr
sey$BGU_boyer_scaled_upr <- b1$upr
sey$BGU_boyer_ecvpred    <- b2$fit   # upper-bound branch

# Boyer sensitivity band = range across the two ACA-reconstruction routes
sey$BGU_boyer_lo <- pmin(sey$BGU_boyer_scaled, sey$BGU_boyer_ecvpred)
sey$BGU_boyer_hi <- pmax(sey$BGU_boyer_scaled, sey$BGU_boyer_ecvpred)

## -- RATIO to modern human (the scope/convention-invariant comparison) --------
# Seymour: flow ratio IS the metabolism ratio.
sey$ratioMH_seymour      <- sey$Total_QICA_cm3_s / Q_MODH
# Boyer: divide by the Boyer prediction for the modern human (same pipeline:
# carotid at r_MODH -> ACA_HOMO, endocranial capacity ECV_MODH).
BGU_boyer_modH <- boyer_pred(ACA_HOMO, ECV_MODH)$fit
sey$ratioMH_boyer_scaled <- sey$BGU_boyer_scaled  / BGU_boyer_modH
sey$ratioMH_boyer_ecvpred<- sey$BGU_boyer_ecvpred / BGU_boyer_modH

## -- scope-adjusted arterial absolutes (put whole-brain BGU on s4a's 6-region
##    cortical+cerebellar footing; assumes covered fraction constant across taxa)
sey$BGU_seymour_s4ascope      <- sey$BGU_seymour      * S4A_SCOPE_FRAC
sey$BGU_boyer_scaled_s4ascope <- sey$BGU_boyer_scaled * S4A_SCOPE_FRAC

## =====================================================================
## 5. Taxon-group tag (for group means / comparison with s4a)
## =====================================================================
grp_of <- function(sp) {
  s <- sp
  ifelse(grepl("neanderthal", s, ignore.case = TRUE), "NT",
  ifelse(grepl("^H\\. sapiens", s),                    "EH/MH",
  ifelse(grepl("erectus|heidelbergensis|rudolfensis|habilis|georgicus|naledi|floresiensis",
               s, ignore.case = TRUE),                 "early Homo",
  ifelse(grepl("^A\\.", s),                            "Australopithecus", "other"))))
}
sey$group_arterial <- grp_of(sey$Species)

# Intra-species temporal grade (only H. sapiens is graded; NA elsewhere). This
# splits the lumped "EH/MH" arterial group into early vs recent vs the modern
# reference. early = the 4 Kochiyama EH fossils PLUS BC1 (Border Cave ~80 ka) and
# LH 18 (Ngaloba ~120 ka), which are early H. sapiens carried by Seymour/Boyer but
# not by Kochiyama; recent = the two Bushman crania; modern_reference = AS8078.
# Note: of the Kochiyama EH set only Skhul 5 is in the arterial (Seymour) sample;
# Qafzeh 9, Mladeč 1, Cro-Magnon 1 are volumetric-only (s4a) — see tables/s4a.
grade_by_specimen <- c("AS8078" = "modern_reference", "M3-A343" = "recent",
  "M4-A344" = "recent", "BC1" = "early", "LH 18" = "early", "Skhul 5" = "early")
sey$sapiens_grade <- ifelse(sey$Specimen %in% names(grade_by_specimen),
                            grade_by_specimen[sey$Specimen], "NA")

## =====================================================================
## 6. Crosswalk Seymour <-> s4a / Kochiyama specimen names, and merge
## =====================================================================
s4aspec <- rd("data_intermediate/s4a_specimen_budgets.csv")
s4aspec$budget_umol_min <- as.numeric(s4aspec$budget_umol_min)
s4aspec$total_cc        <- as.numeric(s4aspec$total_cc)         # Kochiyama GM+WM tissue
s4aspec$budget_ratio_MH <- as.numeric(s4aspec$budget_ratio_MH)  # s4a's own ratio to MH

# Seymour label -> s4a specimen label (only the overlapping individuals)
xwalk <- data.frame(
  seymour_specimen = c("Gibraltar (Forbes Quarry)", "La Chapelle-aux-Saints", "Skhul 5"),
  s4a_specimen      = c("Forbes' Quarry 1", "La Chapelle-aux-Saints 1", "Skhul 5"),
  stringsAsFactors = FALSE)
write_csv_utf8(xwalk, "tables/s4b/s4b_specimen_crosswalk.csv")

sey$s4a_specimen <- xwalk$s4a_specimen[match(sey$Specimen, xwalk$seymour_specimen)]
# s4a quantities (brain-TISSUE based, 6-region cortical+cerebellar scope)
sey$brain_tissue_cc_s4a  <- s4aspec$total_cc[match(sey$s4a_specimen, s4aspec$specimen)]
sey$BGU_s4a_regional     <- s4aspec$budget_umol_min[match(sey$s4a_specimen, s4aspec$specimen)]
sey$ratioMH_s4a          <- s4aspec$budget_ratio_MH[match(sey$s4a_specimen, s4aspec$specimen)]

## =====================================================================
## 7. Assemble outputs
## =====================================================================
rnd <- function(x, d = 2) ifelse(is.na(x), NA, round(x, d))

wide <- data.frame(
  Species            = sey$Species,
  Specimen           = sey$Specimen,
  group_arterial     = sey$group_arterial,
  sapiens_grade      = sey$sapiens_grade,
  s4a_specimen        = sey$s4a_specimen,
  foramen_radius_cm  = sey$Foramen_radius_cm,
  QICA_cm3_s         = sey$Total_QICA_cm3_s,
  # --- volumes, labelled by convention (NOT interchangeable) ---
  endocranial_cc_arterial = sey$Brain_volume_cm3,   # cranial CAPACITY (Seymour/Boyer)
  brain_tissue_cc_s4a      = rnd(sey$brain_tissue_cc_s4a, 0),  # GM+WM tissue (Kochiyama/s4a)
  ACA_scaled_mm2     = rnd(sey$ACA_scaled, 2),
  ACA_ecvpred_mm2    = rnd(sey$ACA_ecvpred, 2),
  # --- whole-brain glucose, umol/min (arterial methods) ---
  BGU_seymour_wholebrain    = rnd(sey$BGU_seymour, 1),
  BGU_boyer_scaled_wholebrain = rnd(sey$BGU_boyer_scaled, 1),
  BGU_boyer_scaled_lwr = rnd(sey$BGU_boyer_scaled_lwr, 1),
  BGU_boyer_scaled_upr = rnd(sey$BGU_boyer_scaled_upr, 1),
  BGU_boyer_ecvpred_wholebrain = rnd(sey$BGU_boyer_ecvpred, 1),  # upper bound
  # --- scope-adjusted to s4a's 6-region cortical+cerebellar footing ---
  BGU_seymour_s4ascope      = rnd(sey$BGU_seymour_s4ascope, 1),
  BGU_boyer_scaled_s4ascope = rnd(sey$BGU_boyer_scaled_s4ascope, 1),
  BGU_s4a_regional          = rnd(sey$BGU_s4a_regional, 1),
  # --- ratio to modern human (PRIMARY comparison; unit/scope invariant) ---
  ratioMH_seymour       = rnd(sey$ratioMH_seymour, 3),
  ratioMH_boyer_scaled  = rnd(sey$ratioMH_boyer_scaled, 3),
  ratioMH_boyer_ecvpred = rnd(sey$ratioMH_boyer_ecvpred, 3),
  ratioMH_s4a            = rnd(sey$ratioMH_s4a, 3),
  stringsAsFactors = FALSE)
write_csv_utf8(wide, "data_intermediate/s4b_fossil_estimates.csv")

# tidy long form (one row per specimen x method); scale = metabolic scope,
# quantity = ratio-to-MH (comparable) plus native absolute (scope-specific).
mk_long <- function(method, scale, val, ratio, lo = NA, hi = NA) data.frame(
  Species = sey$Species, Specimen = sey$Specimen,
  group_arterial = sey$group_arterial, sapiens_grade = sey$sapiens_grade,
  method = method, scope = scale,
  BGU_umol_min = rnd(val, 1), ratio_MH = rnd(ratio, 3),
  lwr = rnd(lo, 1), hi = rnd(hi, 1),
  stringsAsFactors = FALSE)
long <- rbind(
  mk_long("Seymour_flow",      "whole_brain",        sey$BGU_seymour,      sey$ratioMH_seymour),
  mk_long("Boyer_ACA_scaled",  "whole_brain",        sey$BGU_boyer_scaled, sey$ratioMH_boyer_scaled,
          sey$BGU_boyer_scaled_lwr, sey$BGU_boyer_scaled_upr),
  mk_long("Boyer_ACA_ecvpred", "whole_brain",        sey$BGU_boyer_ecvpred, sey$ratioMH_boyer_ecvpred),
  mk_long("s4a_volume",         "cortical+cerebellar", sey$BGU_s4a_regional,  sey$ratioMH_s4a))
long <- long[!is.na(long$BGU_umol_min), ]
write_csv_utf8(long, "data_intermediate/s4b_fossil_estimates_long.csv")

## -- overlap specimens: all methods side by side ------------------------------
## Compare on ratio_MH (invariant). Absolutes shown in each method's native
## scope; volumes labelled by convention.
ov <- wide[!is.na(wide$s4a_specimen),
           c("Specimen", "s4a_specimen", "QICA_cm3_s",
             "endocranial_cc_arterial", "brain_tissue_cc_s4a",
             "ratioMH_seymour", "ratioMH_boyer_scaled", "ratioMH_boyer_ecvpred",
             "ratioMH_s4a",
             "BGU_seymour_s4ascope", "BGU_boyer_scaled_s4ascope", "BGU_s4a_regional")]
write_csv_utf8(ov, "tables/s4b/s4b_method_comparison_overlap.csv")
cat("=== overlap specimens (in both Seymour and s4a); compare on ratio_MH ===\n")
print(ov, row.names = FALSE)
cat("\n")

## -- group means (arterial methods) + s4a group means for comparison -----------
grp_summ <- do.call(rbind, lapply(split(sey, sey$group_arterial), function(d) data.frame(
  group_arterial          = d$group_arterial[1], n = nrow(d),
  endocranial_cc_mean     = round(mean(d$Brain_volume_cm3), 0),
  ratioMH_seymour_mean    = round(mean(d$ratioMH_seymour), 3),
  ratioMH_boyer_scaled_mean = round(mean(d$ratioMH_boyer_scaled), 3),
  BGU_seymour_wholebrain_mean = round(mean(d$BGU_seymour), 1),
  BGU_boyer_scaled_wholebrain_mean = round(mean(d$BGU_boyer_scaled), 1),
  stringsAsFactors = FALSE)))
grp_summ <- grp_summ[order(-grp_summ$endocranial_cc_mean), ]
write_csv_utf8(grp_summ, "tables/s4b/s4b_group_means.csv")
cat("=== arterial-method group means (endocranial capacity; whole-brain BGU) ===\n")
print(grp_summ, row.names = FALSE)

## -- H. sapiens grade means: split the lumped EH/MH group ---------------------
## early vs recent vs modern_reference (arterial methods). See sapiens_grade above.
hs <- sey[sey$group_arterial == "EH/MH", ]
grade_summ <- do.call(rbind, lapply(c("early", "recent", "modern_reference"), function(g) {
  d <- hs[hs$sapiens_grade == g, ]
  if (nrow(d) == 0) return(NULL)
  data.frame(sapiens_grade = g, n = nrow(d),
    specimens = paste(d$Specimen, collapse = "; "),
    endocranial_cc_mean = round(mean(d$Brain_volume_cm3), 0),
    ratioMH_seymour_mean = round(mean(d$ratioMH_seymour), 3),
    ratioMH_boyer_scaled_mean = round(mean(d$ratioMH_boyer_scaled), 3),
    BGU_seymour_wholebrain_mean = round(mean(d$BGU_seymour), 1),
    BGU_boyer_scaled_wholebrain_mean = round(mean(d$BGU_boyer_scaled), 1),
    stringsAsFactors = FALSE)
}))
write_csv_utf8(grade_summ, "tables/s4b/s4b_sapiens_grade_means.csv")
cat("=== H. sapiens grade means (arterial; early vs recent vs modern reference) ===\n")
print(grade_summ, row.names = FALSE)

cat(sprintf("\n[note] Boyer modern-human reference BGU = %.1f umol/min (whole brain);\n",
            BGU_boyer_modH))
cat(sprintf("       s4a modern-human budget = %.1f umol/min (6-region cortical+cerebellar,\n",
            S4A_MODH_REGIONAL))
cat(sprintf("       = %.1f%% of whole brain). Absolute umol/min compare only within a scope;\n",
            100 * S4A_SCOPE_FRAC))
cat("       cross-method comparison uses ratio_MH.\n")

s4agrp_path <- "tables/s4a/species_absolute_budgets.csv"
if (file.exists(s4agrp_path)) {
  s4agrp <- rd(s4agrp_path)
  cat("\n=== s4a volume-based group means (brain TISSUE; cortical+cerebellar scope) ===\n")
  print(s4agrp[, c("group", "n", "cerebral_cc", "budget_umol_min",
                  "budget_ci95_halfwidth", "budget_ratio_MH")], row.names = FALSE)
}

cat("\n[s4b] done. Outputs in data_intermediate/s4b_* and tables/s4b/.\n")
