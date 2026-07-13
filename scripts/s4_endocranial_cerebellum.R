# =====================================================================
# s4_endocranial_cerebellum.R
# ---------------------------------------------------------------------
# CEREBELLUM-ONLY energy budget for fossil hominins.
#
# Purpose: test whether the elevated whole-brain budget of Neanderthals
# (NT) and early modern humans (EH) relative to modern humans (MH), seen
# in s4_endocranial.R, is driven by the CEREBELLUM. This script is a
# self-contained sibling of s4_endocranial.R that uses ONLY cerebellar
# data and reports ONLY the cerebellar budget.
#
# METHOD (identical model to s4_endocranial.R, restricted to cerebellum):
#   budget_part = rCMRGlc_part x (MH_mass_part / 100 g) x shape_part(group)
#                                                       x cerebellar_size(indiv)
#   where
#     rCMRGlc_part  = modern-human glucose rate, per 100 g/min (Heiss 2004)
#     MH_mass_part  = MH sub-region volume (cc) x 1.036 g/cc  (Fig 3 legend)
#     shape_part    = Fig 3A size-adjusted relative volume for the group
#     cerebellar_size = (individual cerebellar volume) / (MH cerebellar volume)
#
# The three cerebellar parcels (Kochiyama Fig 3A / legend):
#     Ce A  anterior cerebellar cortex   rate 29.8
#     Ce P  posterior cerebellar cortex  rate 29.8
#     Ce V  vermis                       rate 30.1
#   Cerebellar cortex = Ce A + Ce P.
#
# SHAPE vs SIZE (as in the main script):
#   Fig 3A ratios are ICV-size-adjusted, so they describe cerebellar
#   SHAPE (how the cerebellum is apportioned) at matched brain size. The
#   individual's own cerebellar VOLUME restores absolute SIZE.
#
# VOLUME BASIS (CSF): the Fig-3-legend volumes are parcellated brain
# tissue (grey + white), CSF already excluded; density 1.036 g/cc treats
# them as tissue.
#
# NOTE on the size reference: CEREBELLAR_MH = 149 cc is the modern-human
# whole-cerebellum group mean (FossilSpecimensText groupmeans, n = 1185),
# the same reference used in s4_endocranial.R. The digitized parcels sum
# to ~140.6 cc; because every quantity is expressed as a RATIO to MH, the
# choice of reference cancels and does not affect the group/specimen ratios.
# =====================================================================

suppressWarnings(suppressMessages({
  # base R only
}))

## ---- paths --------------------------------------------------------
# Resolve project root: prefer the sourced script's location, else cwd.
this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (!is.null(this_file)) {
  root <- normalizePath(file.path(dirname(this_file), ".."), mustWork = FALSE)
} else {
  root <- "."
}
if (!dir.exists(file.path(root, "data_raw"))) root <- "."   # fall back to cwd
raw_dir   <- file.path(root, "data_raw")
tab_dir   <- file.path(root, "tables", "s4")
fig_dir   <- file.path(root, "figs",   "s4")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Figure export helper: render a base-graphics figure to BOTH a raster PNG
# (slides) and a vector PDF (print). `draw` is replayed once per device so the
# files are identical. cairo_pdf sizes are in inches = pixels / res.
save_png_pdf <- function(stem, draw, width, height, res = 220) {
  png(paste0(stem, ".png"), width = width, height = height, res = res)
  draw(); dev.off()
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width / res, height = height / res)
  draw(); dev.off()
}

## ---- constants ----------------------------------------------------
BRAIN_DENSITY_G_PER_CC <- 1.036
CEREBELLAR_MH          <- 149      # MH whole-cerebellum mean (cc), n = 1185

## ---- cerebellar parcels (rate, MH volume, Fig 3A shape) -----------
# rates: Heiss et al. 2004 TABLE1 (per 100 g/min)
# MH volumes: Kochiyama Fig 3 legend (cc); shape: Kochiyama Fig 3A
parcel <- data.frame(
  code   = c("Ce A", "Ce P", "Ce V"),
  label  = c("Anterior cerebellar cortex",
             "Posterior cerebellar cortex",
             "Vermis"),
  rate   = c(29.8, 29.8, 30.1),
  mh_cc  = c(13.86, 114.41, 12.38),
  NT     = c(0.990, 0.940, 0.965),
  EH     = c(1.000, 0.995, 1.010),
  MH     = c(1.000, 1.000, 1.000),
  stringsAsFactors = FALSE
)
# per-part MH budget base (umol/min): rate x mass(g)/100
parcel$mass_g    <- parcel$mh_cc * BRAIN_DENSITY_G_PER_CC
parcel$base_umol <- parcel$rate * parcel$mass_g / 100

## ---- size-adjusted cerebellar SHAPE factor per group --------------
# rate-and-mass-weighted mean of the parcel shape ratios; = cerebellar
# budget ratio at matched (MH) size.
groups <- c("NT", "EH", "MH")
shape_factor <- sapply(groups, function(g)
  sum(parcel$base_umol * parcel[[g]]) / sum(parcel$base_umol))

## ---- MH reference cerebellar budget (umol/min) --------------------
budget_MH_cerebellum <- sum(parcel$base_umol)   # shape = 1, size = 1

## ================================================================
## 1. SPECIES cerebellum budget (group shape x group size)
## ================================================================
gm <- read.csv(file.path(raw_dir,
        "Kochiyama_etal_2018_FossilSpecimensText_groupmeans.csv"),
        check.names = FALSE)
gm_cb <- data.frame(
  group        = gm$Group,
  cerebellar_cc= gm$Cerebellum_mean_cc,
  cerebellar_sd= gm$Cerebellum_sd_cc,
  n            = gm$n
)
gm_cb$size_ratio    <- gm_cb$cerebellar_cc / CEREBELLAR_MH
gm_cb$shape_factor  <- shape_factor[gm_cb$group]
gm_cb$budget_umol   <- budget_MH_cerebellum * gm_cb$shape_factor * gm_cb$size_ratio
gm_cb$budget_ratio  <- gm_cb$budget_umol / budget_MH_cerebellum
# propagate the (dominant) size uncertainty into the budget ratio
gm_cb$budget_ratio_sd <- gm_cb$shape_factor * gm_cb$cerebellar_sd / CEREBELLAR_MH
gm_cb$ci95_halfwidth  <- 1.96 * gm_cb$budget_ratio_sd / sqrt(gm_cb$n)

species_out <- gm_cb[, c("group","n","cerebellar_cc","cerebellar_sd",
                         "size_ratio","shape_factor",
                         "budget_umol","budget_ratio",
                         "budget_ratio_sd","ci95_halfwidth")]
write.csv(species_out,
          file.path(tab_dir, "cerebellum_species_budget.csv"),
          row.names = FALSE)

## ================================================================
## 2. PER-SPECIMEN cerebellum budget (group shape x specimen size)
## ================================================================
sp <- read.csv(file.path(raw_dir,
        "Kochiyama_etal_2018_FossilSpecimensText.csv"),
        check.names = FALSE, encoding = "UTF-8")
grp_map <- c("Neanderthal" = "NT", "early Homo sapiens" = "EH")
spec <- data.frame(
  specimen      = sp$Specimen,
  group         = grp_map[sp$Taxon_group],
  cerebellar_cc = sp$`Cerebellum_Vol.cc`,
  age_min_kya   = sp$date_min_yBP / 1000,
  age_max_kya   = sp$date_max_yBP / 1000,
  age_mid       = sp$date_mean_yBP / 1000,     # paper's own midpoint
  stringsAsFactors = FALSE
)
spec$size_ratio   <- spec$cerebellar_cc / CEREBELLAR_MH
spec$shape_factor <- shape_factor[spec$group]
spec$budget_umol  <- budget_MH_cerebellum * spec$shape_factor * spec$size_ratio
spec$budget_ratio <- spec$budget_umol / budget_MH_cerebellum

# --- specimen sex estimates -------------------------------------------------
# Sex estimates compiled from the dedicated demographics lookup
# (material_other_sources/Kochiyama_fossil_demographics_journal_sources.xlsx),
# one journal citation per specimen (see that file's "Preferred journal
# citation" column). Codes: M / F when reported; "M?" / "F?" when contested;
# "?" when unknown. Qafzeh 9 is "M?" (2021 reassessment male; historically
# debated/female).
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

spec <- spec[order(spec$group, -spec$budget_ratio), ]
write.csv(spec[, c("specimen","group","cerebellar_cc",
                   "size_ratio","shape_factor",
                   "budget_umol","budget_ratio")],
          file.path(tab_dir, "cerebellum_specimen_budget.csv"),
          row.names = FALSE)

## ================================================================
## 3. PARCEL breakdown table (which sub-region carries the effect)
## ================================================================
parcel_out <- parcel[, c("code","label","rate","mh_cc","mass_g",
                         "base_umol","NT","EH","MH")]
names(parcel_out)[7:9] <- paste0("shape_", names(parcel_out)[7:9])
write.csv(parcel_out,
          file.path(tab_dir, "cerebellum_parcel_breakdown.csv"),
          row.names = FALSE)

## ================================================================
## 4. FIGURES
## ================================================================
grp_col <- c(NT = "#4C72B0", EH = "#DD8452", MH = "#55A868")

## --- 4a. species cerebellum budget +/- 95% CI ---------------------
png(file.path(fig_dir, "cerebellum_species_budget.png"),
    width = 1900, height = 1300, res = 200)
op <- par(mar = c(5, 6, 4, 2))
o  <- match(c("NT","EH","MH"), species_out$group)
b  <- species_out$budget_ratio[o]
h  <- species_out$ci95_halfwidth[o]
xpos <- 1:3
plot(xpos, b, pch = 19, cex = 2.4, col = grp_col[species_out$group[o]],
     xlim = c(0.5, 3.5), ylim = range(c(b - h, b + h, 1)) + c(-0.02, 0.02),
     xaxt = "n", xlab = "", ylab = "Cerebellum energy budget (ratio to MH)",
     main = "Cerebellum-only energy budget by species", cex.lab = 1.15)
arrows(xpos, b - h, xpos, b + h, angle = 90, code = 3,
       length = 0.08, lwd = 2, col = grp_col[species_out$group[o]])
axis(1, at = xpos, labels = c("Neanderthal","early H. sapiens","modern human"))
abline(h = 1, lty = 2, col = "grey50")
text(xpos, b, sprintf("%.3f", b), pos = 4, offset = 1.1, cex = 0.95)
par(op); dev.off()

## --- 4b. per-specimen cerebellum budget vs volume -----------------
png(file.path(fig_dir, "cerebellum_specimen_budget.png"),
    width = 2000, height = 1300, res = 200)
op <- par(mar = c(5, 6, 4, 2))
sp2 <- spec[order(spec$budget_ratio), ]
cols <- grp_col[sp2$group]
bp <- barplot(sp2$budget_ratio, horiz = TRUE, col = cols, border = NA,
              xlim = c(0, max(sp2$budget_ratio) * 1.15),
              names.arg = sp2$specimen, las = 1,
              xlab = "Cerebellum energy budget (ratio to MH)",
              main = "Cerebellum-only budget per fossil specimen",
              cex.names = 0.85, cex.lab = 1.1)
abline(v = 1, lty = 2, col = "grey40")
text(sp2$budget_ratio, bp, sprintf("%.2f", sp2$budget_ratio),
     pos = 4, cex = 0.8, xpd = NA)
legend("bottomright", legend = c("Neanderthal","early H. sapiens"),
       fill = grp_col[c("NT","EH")], border = NA, bty = "n", cex = 0.9)
par(op); dev.off()

## --- 4c/4d. cerebellum volume through time (point size = cerebellum budget)
## Mirror of s4_endocranial.R's volume_timeline figures, restricted to the
## cerebellum: y = cerebellar volume (cc); point AREA = cerebellum budget
## ratio to MH; colour = group; shape = specimen. Linear + log10 age axes.
spec$age_plot <- ifelse(is.na(spec$age_mid),
                        mean(spec$age_mid[spec$group == "NT"], na.rm = TRUE),
                        spec$age_mid)
spec$undated  <- is.na(spec$age_mid)
spec$col      <- grp_col[spec$group]
fill_pch      <- c(21, 22, 23, 24)                # circle, square, diamond, tri
spec$pch      <- ave(seq_len(nrow(spec)), spec$group,
                     FUN = function(ix) fill_pch[seq_along(ix)])
MH_PCH        <- 23
MH_LOG_AGE    <- 5                                # kya; nominal recent pos. for MH

# marker AREA mapped linearly across the observed budget range (as in main script)
SIZE_CEX  <- c(1.5, 4.6)
size_vals <- c(spec$budget_ratio, 1.0)
size_rng  <- range(size_vals)
cex_from_budget <- function(r)
  sqrt(SIZE_CEX[1]^2 + (SIZE_CEX[2]^2 - SIZE_CEX[1]^2) *
         (r - size_rng[1]) / diff(size_rng))
spec$cex_budget <- cex_from_budget(spec$budget_ratio)
mh_cex_budget   <- cex_from_budget(1.0)

cb_scatter <- function(main, note = "", xlog = FALSE) {
  yv <- spec$cerebellar_cc; mh <- CEREBELLAR_MH
  yr <- range(c(yv, mh)); pad <- diff(yr) * 0.08
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
  } else axis(1, at = seq(800, 0, by = -200))
  axis(4); mtext("Cerebellar volume (cc)", side = 4, line = 3.6, cex = 1.0)
  title(main = main, cex.main = 1.1, font.main = 1)

  rect(800, ylim[1], 600, ylim[2], col = "#00000010", border = NA)
  div_x <- if (xlog) 470 else 560
  div_y <- ylim[1] + diff(ylim) * 0.24
  text(div_x, div_y, "NT-MH divergence\n(~600-800 kya)",
       cex = 0.72, col = "grey30", adj = c(0.5, 0.5))
  abline(v = mh_x, col = "grey80", lty = 3)

  lab_x <- if (xlog) 2.4 else -52
  # dashed group-mean lines use the SAME published cerebellar group means
  # (FossilSpecimensText groupmeans) that the budget table uses -- not the
  # 4-specimen recomputed mean -- so figure labels match cerebellum_species_budget.csv.
  gmean <- setNames(gm_cb$cerebellar_cc, gm_cb$group)
  means <- c(NT = unname(gmean["NT"]), EH = unname(gmean["EH"]), MH = mh)
  # NT and MH cerebellar means are identical (149 cc), so their lines coincide.
  # Draw MH DOTTED (lty 3) over NT DASHED (lty 2) so both colours remain visible
  # on the shared line; draw EH first (it is distinct).
  abline(h = means["EH"], col = grp_col["EH"], lty = 2, lwd = 1.2)
  abline(h = means["NT"], col = grp_col["NT"], lty = 2, lwd = 1.4)
  abline(h = means["MH"], col = grp_col["MH"], lty = 3, lwd = 1.4)
  # numeric labels: stack vertically when two group means fall within 1.5 cc
  ord_lab <- order(means)
  offs <- rep(0, 3); names(offs) <- names(means)
  for (k in seq_along(ord_lab)[-1]) {
    a <- names(means)[ord_lab[k]]; b <- names(means)[ord_lab[k-1]]
    if (abs(means[a] - means[b]) < 1.5) offs[a] <- offs[b] + 0.9
  }
  for (g in names(means))
    text(lab_x, means[g], sprintf("%.0f", means[g]), col = grp_col[g],
         cex = 0.72, adj = c(1, -0.4 - offs[g]), xpd = NA)

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

  leg_lab <- sprintf("%s [%s]  (budget %.2f)", spec$specimen, spec$sex, spec$budget_ratio)
  legend("topleft", legend = c(leg_lab, "Modern human (present, 1.00)"),
         pch = c(spec$pch, MH_PCH),
         pt.bg = c(ifelse(spec$undated, NA, paste0(spec$col, "CC")), paste0(grp_col["MH"],"CC")),
         col = c(spec$col, grp_col["MH"]),
         pt.cex = 1.3, bty = "n", cex = 0.72, y.intersp = 1.05)
  s_ref <- pretty(size_rng, n = 3); s_ref <- s_ref[s_ref >= size_rng[1] & s_ref <= size_rng[2]]
  legend("bottomleft", title = "Cerebellum budget (ratio to MH)",
         legend = sprintf("%.2f", s_ref), pch = 21, col = "grey40", pt.bg = "grey80",
         pt.cex = cex_from_budget(s_ref), bty = "n", cex = 0.72,
         y.intersp = 1.5, x.intersp = 1.4)
  if (nzchar(note))
    mtext(note, side = 1, line = 4.1, cex = 0.60, adj = 0, col = "grey40")
  invisible(gmean)
}

png(file.path(fig_dir, "cerebellum_volume_timeline.png"),
    width = 2100, height = 1450, res = 220)
cb_scatter(main = "Fossil cerebellum volume through time (point size = energy budget)",
           note = "Colour = group; shape = specimen; point AREA = cerebellum energy-budget ratio to MH (legend). Whiskers = dating range; hollow = undated (Forbes' Quarry 1). Dashed = group means; the NT mean (149 cc) equals the MH mean, so those two lines coincide (MH drawn dotted over NT).")
dev.off()

save_png_pdf(file.path(fig_dir, "cerebellum_volume_timeline_logage"),
  function() cb_scatter(main = "Fossil cerebellum volume through time (log age; point size = energy budget)",
           xlog = TRUE,
           note = "Log10 age axis. *Modern human plotted at a nominal ~5 kya (age 0 undefined on log scale). Point AREA = cerebellum energy-budget ratio to MH; colour = group; shape = specimen. The NT mean (149 cc) equals the MH mean, so those dashed lines coincide (MH drawn dotted over NT). Markers sit at the reported mean date; on a log axis a linearly-symmetric date range is not symmetric, so points appear slightly toward the older (left) end of their whisker."),
  width = 2100, height = 1450, res = 220)

## --- 4d-flip. cerebellum glucose utilization through time; point size = volume
## Second version of the log-age figure with the two encodings SWAPPED:
##   Y axis     = cerebellum glucose utilization (ratio to MH)  [was point size]
##   point AREA = cerebellar volume (cc)                        [was the y axis]
## "Energy budget" is relabelled to the more specific "glucose utilization".
## marker AREA mapped LINEARLY across the observed cerebellar-volume range
CBVOL_SIZE_CEX <- c(1.5, 4.6)
cbvol_vals <- c(spec$cerebellar_cc, CEREBELLAR_MH)
cbvol_rng  <- range(cbvol_vals)
cex_from_cbvol <- function(v)
  sqrt(CBVOL_SIZE_CEX[1]^2 + (CBVOL_SIZE_CEX[2]^2 - CBVOL_SIZE_CEX[1]^2) *
         (v - cbvol_rng[1]) / diff(cbvol_rng))
spec$cex_cbvol <- cex_from_cbvol(spec$cerebellar_cc)
mh_cex_cbvol   <- cex_from_cbvol(CEREBELLAR_MH)

cb_scatter_flip <- function(main, note = "", xlog = FALSE) {
  yv <- spec$budget_ratio; mh <- 1.0               # glucose utilization; MH = 1.00
  yr <- range(c(yv, mh)); pad <- diff(yr) * 0.08
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
  } else axis(1, at = seq(800, 0, by = -200))
  axis(4); mtext("Cerebellum glucose utilization (ratio to MH)", side = 4, line = 3.6, cex = 1.0)
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
           col = spec$col[i], cex = spec$cex_cbvol[i], lwd = 1.8)
  }
  points(mh_x, mh, pch = MH_PCH, bg = paste0(grp_col["MH"], "CC"),
         col = grp_col["MH"], cex = mh_cex_cbvol, lwd = 1.8)

  leg_lab <- sprintf("%s [%s]  (%.0f cc)", spec$specimen, spec$sex, spec$cerebellar_cc)
  legend("topleft", legend = c(leg_lab, sprintf("Modern human (present, %.0f cc)", CEREBELLAR_MH)),
         pch = c(spec$pch, MH_PCH),
         pt.bg = c(ifelse(spec$undated, NA, paste0(spec$col, "CC")), paste0(grp_col["MH"],"CC")),
         col = c(spec$col, grp_col["MH"]),
         pt.cex = 1.3, bty = "n", cex = 0.72, y.intersp = 1.05)
  v_ref <- pretty(cbvol_rng, n = 3); v_ref <- v_ref[v_ref >= cbvol_rng[1] & v_ref <= cbvol_rng[2]]
  legend("bottomleft", title = "Cerebellar volume (cc)",
         legend = sprintf("%.0f", v_ref), pch = 21, col = "grey40", pt.bg = "grey80",
         pt.cex = cex_from_cbvol(v_ref), bty = "n", cex = 0.72,
         y.intersp = 1.5, x.intersp = 1.4)
  if (nzchar(note))
    mtext(note, side = 1, line = 4.1, cex = 0.60, adj = 0, col = "grey40")
  invisible(gmean)
}

png(file.path(fig_dir, "cerebellum_volume_timeline_logage_flipped.png"),
    width = 2100, height = 1450, res = 220)
cb_scatter_flip(
  main = "Fossil cerebellum glucose utilization through time (log age; point size = volume)",
  xlog = TRUE,
  note = "Log10 age axis. *Modern human plotted at a nominal ~5 kya (age 0 undefined on log scale). Y = cerebellum glucose utilization (ratio to MH); point AREA = cerebellar volume (cc); colour = group; shape = specimen. Markers sit at the reported mean date; on a log axis a linearly-symmetric date range is not symmetric, so points appear slightly toward the older (left) end of their whisker.")
dev.off()

## ================================================================
## 5. CONSOLE SUMMARY
## ================================================================
cat("=== CEREBELLUM-ONLY BUDGET ===\n")
cat(sprintf("MH reference cerebellum budget: %.2f umol/min\n\n",
            budget_MH_cerebellum))
cat("Size-adjusted cerebellar SHAPE factor (matched-size budget ratio):\n")
print(round(shape_factor, 4))
cat("\nSpecies cerebellum budget (shape x size):\n")
print(species_out[, c("group","cerebellar_cc","size_ratio",
                      "shape_factor","budget_ratio")], row.names = FALSE)
cat("\nPer-specimen cerebellum budget:\n")
print(spec[, c("specimen","group","cerebellar_cc","budget_ratio")],
      row.names = FALSE)

## ================================================================
## 5. Significance summary figure: whole-brain vs cerebellum budgets
##    (mean +/- 95% CI vs MH = 1.0, with two-sided p-values)
## ----------------------------------------------------------------
## Rebuilds figs/s4/budget_significance.{png,pdf}. Whole-brain group budgets
## are read from the main script's output (species_absolute_budgets.csv);
## cerebellum group budgets are `species_out` computed above. The p-value is a
## two-sided one-sample t-test that the group-mean budget ratio differs from
## MH = 1.0, with SEM = sd / sqrt(n) and df = n - 1 (fossil n = 4). Run
## s4_endocranial.R first so the whole-brain table exists.
wb_path <- file.path(tab_dir, "species_absolute_budgets.csv")
if (file.exists(wb_path)) {
  wb <- read.csv(wb_path, stringsAsFactors = FALSE)
  wb <- wb[match(c("NT","EH"), wb$group), ]
  cb <- species_out[match(c("NT","EH"), species_out$group), ]

  sig <- data.frame(
    group = c("NT","EH","NT","EH"),
    label = c("Whole-brain\nNeanderthal", "Whole-brain\nearly H.sapiens",
              "Cerebellum\nNeanderthal",  "Cerebellum\nearly H.sapiens"),
    ratio = c(wb$budget_ratio_MH,      cb$budget_ratio),
    ci    = c(wb$budget_ratio_MH_ci95, cb$ci95_halfwidth),
    sd    = c(wb$budget_ratio_MH_sd,   cb$budget_ratio_sd),
    n     = c(wb$n,                    cb$n),
    y     = c(4.3, 3.5, 1.8, 1.0),
    stringsAsFactors = FALSE)
  sig$sem <- sig$sd / sqrt(sig$n)
  sig$p   <- 2 * pt(-abs((sig$ratio - 1) / sig$sem), df = pmax(sig$n - 1, 1))
  sig$col <- grp_col[sig$group]

  draw_budget_sig <- function() {
    op <- par(mar = c(5.2, 9.0, 3.4, 6.2), mgp = c(3, 0.7, 0)); on.exit(par(op))
    xlim <- range(c(sig$ratio - sig$ci, sig$ratio + sig$ci, 1))
    xlim <- xlim + c(-0.06, 0.06) * diff(xlim)
    plot(NA, xlim = xlim, ylim = c(0.4, 4.9),
         xlab = "Energy budget (ratio to MH)", ylab = "", yaxt = "n", bty = "n")
    rect(par("usr")[1], 3.15, par("usr")[2], 4.65, col = "#00000008", border = NA)
    rect(par("usr")[1], 0.65, par("usr")[2], 2.15, col = "#00000008", border = NA)
    abline(v = 1, lty = 2, col = "grey40")
    arrows(sig$ratio - sig$ci, sig$y, sig$ratio + sig$ci, sig$y,
           code = 3, angle = 90, length = 0.05, col = sig$col, lwd = 2.4)
    points(sig$ratio, sig$y, pch = 19, col = sig$col, cex = 1.9)
    axis(2, at = sig$y, labels = sig$label, las = 1, tick = FALSE, cex.axis = 0.8)
    mtext("Whole-brain budget", side = 2, line = 7.2, at = 3.9, cex = 0.85)
    mtext("Cerebellum budget",  side = 2, line = 7.2, at = 1.4, cex = 0.85)
    text(par("usr")[2], sig$y, sprintf("p = %.2f", sig$p),
         pos = 4, xpd = NA, cex = 0.9, col = "grey20")
    title(main = "Are the budget differences significant? (mean +/- 95% CI, vs MH = 1.0)",
          cex.main = 1.05, font.main = 1)
    mtext(sprintf("all CIs cross 1.0  ->  none significant (n = %d/group)", sig$n[1]),
          side = 1, line = 3.2, adj = 1, cex = 0.7, col = "grey40")
  }
  save_png_pdf(file.path(fig_dir, "budget_significance"),
               draw_budget_sig, width = 2000, height = 1150, res = 220)
  write.csv(sig[, c("group","label","ratio","ci","sd","n","p")],
            file.path(tab_dir, "budget_significance.csv"), row.names = FALSE)
  cat("\nbudget_significance figure written (whole-brain + cerebellum).\n")
} else {
  message("budget_significance skipped: ", wb_path,
          " not found -- run s4_endocranial.R first.")
}
