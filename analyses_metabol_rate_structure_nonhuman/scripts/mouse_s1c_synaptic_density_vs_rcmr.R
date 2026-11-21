############################################################
# mouse_s1c -- Regional synaptic density vs mouse regional glucose
# utilisation (CMRgl). Mouse analogue of s1c_synaptic_density_metabolic_rate
# (human: Johansen 2024 SV2A Bmax vs Heiss 2004 rCMRGlc).
#
# Input : data_raw/mouse_synaptic_density_by_region.csv  (contract in
#         metadata/mouse_synaptic_density_input_spec.md; env MOUSE_SYNAPSE_CSV overrides)
# Output: tables/mouse_s1c/{matched_regions,fits,leave_one_out}.csv,
#         figs/mouse_s1c/synaptic_density_vs_rcmr.{pdf,jpg}
# Stops with a clear message if the input file is absent.
############################################################
.a <- commandArgs(FALSE); .f <- sub("^--file=", "", .a[grepl("^--file=", .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else "analyses_metabol_rate_structure_nonhuman/scripts", "common.R"))
source(nh("scripts/mouse_s1b_helpers.R"))
for (d in c("tables/mouse_s1c", "figs/mouse_s1c")) dir.create(nh(d), recursive = TRUE, showWarnings = FALSE)

path <- Sys.getenv("MOUSE_SYNAPSE_CSV", unset = nh("data_raw/mouse_synaptic_density_by_region.csv"))
if (!file.exists(path))
  stop("Missing ", path, ". See ", nh("metadata/mouse_synaptic_density_input_spec.md"),
       " for the required columns and candidate sources.", call. = FALSE)

syn <- read_csv(path, show_col_types = FALSE)
req <- c("cmrgl_region", "synaptic_density", "measure", "unit", "source")
if (length(m <- setdiff(req, names(syn)))) stop("Input lacks column(s): ", paste(m, collapse = ", "), call. = FALSE)
if (!"analysis_set" %in% names(syn)) syn$analysis_set <- "all"
if (!"primary_included" %in% names(syn)) syn$primary_included <- TRUE
if (!"nested_in_other_region" %in% names(syn)) syn$nested_in_other_region <- FALSE
if (anyNA(syn$synaptic_density) || !is.numeric(syn$synaptic_density)) stop("synaptic_density must be numeric and complete.", call. = FALSE)
if (anyDuplicated(syn[c("analysis_set", "cmrgl_region")])) stop("One row per analysis_set x cmrgl_region is required.", call. = FALSE)
if (any(tapply(paste(syn$measure, syn$unit), syn$analysis_set, function(v) length(unique(v))) != 1L))
  stop("A single measure and unit per analysis_set is required.", call. = FALSE)

rates <- read_mouse_rates()
if (length(m <- setdiff(syn$cmrgl_region, rates$cmrgl_region)))
  stop("Region(s) without a mouse rate: ", paste(m, collapse = ", "), call. = FALSE)
d <- inner_join(syn, rates, by = "cmrgl_region")
write_csv(d, nh("tables/mouse_s1c/matched_regions.csv"))

# Fit every analysis set on (a) primary regions and (b) all non-nested regions (adds partial matches).
# Nested pools (e.g. Neocortex) are never fitted together with their components.
fit_one <- function(z, label, set) {
  n <- nrow(z)
  if (n < 4) return(NULL)
  ct <- suppressWarnings(cor.test(z$rcmr_value, z$synaptic_density, method = "spearman", exact = FALSE))
  fit <- lm(rcmr_value ~ synaptic_density, data = z); sm <- summary(fit)
  loo <- vapply(seq_len(n), function(i) suppressWarnings(
    unname(cor.test(z$rcmr_value[-i], z$synaptic_density[-i], method = "spearman", exact = FALSE)$estimate)), numeric(1))
  tibble(analysis_set = set, region_subset = label, n = n, measure = unique(z$measure), unit = unique(z$unit),
         spearman_rho = unname(ct$estimate), spearman_p = ct$p.value, lm_slope = coef(fit)[[2]],
         lm_r2 = sm$r.squared, lm_p = sm$coefficients[2, 4], loo_rho_min = min(loo), loo_rho_max = max(loo),
         regions = paste(z$cmrgl_region, collapse = "+"))
}
d_fit <- filter(d, !nested_in_other_region)
fits <- bind_rows(lapply(unique(d_fit$analysis_set), function(s) {
  z <- filter(d_fit, analysis_set == s)
  bind_rows(fit_one(filter(z, primary_included), "primary", s), fit_one(z, "primary_plus_partial", s))
}))
if (!nrow(fits)) stop("Fewer than 4 matched regions in every analysis set; no inference attempted.", call. = FALSE)
write_csv(fits, nh("tables/mouse_s1c/fits.csv"))
print(as.data.frame(select(fits, -regions)), digits = 3)

loo <- bind_rows(lapply(unique(d_fit$analysis_set), function(s) {
  z <- filter(d_fit, analysis_set == s, primary_included)
  if (nrow(z) < 5) return(NULL)
  bind_rows(lapply(seq_len(nrow(z)), function(i) {
    r <- suppressWarnings(cor.test(z$rcmr_value[-i], z$synaptic_density[-i], method = "spearman", exact = FALSE))
    tibble(analysis_set = s, left_out = z$cmrgl_region[i], n = nrow(z) - 1, spearman_rho = unname(r$estimate), spearman_p = r$p.value)
  }))
}))
write_csv(loo, nh("tables/mouse_s1c/leave_one_out.csv"))

# Figure: one panel per analysis set (primary regions, partial matches hollow)
pd <- d_fit %>% left_join(distinct(fits, analysis_set) %>% mutate(ok = TRUE), by = "analysis_set") %>% filter(ok)
pf <- fits %>% filter(region_subset == "primary") %>% transmute(analysis_set, lab = sprintf("%s\nrho=%.2f, p=%.3f, n=%d", analysis_set, spearman_rho, spearman_p, n))
pd <- inner_join(pd, pf, by = "analysis_set")
p <- ggplot(pd, aes(synaptic_density, rcmr_value)) +
  geom_smooth(data = filter(pd, primary_included), method = "lm", formula = y ~ x, se = FALSE, colour = "steelblue", linewidth = 0.6) +
  geom_point(aes(shape = primary_included), size = 2.2) + scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = "Primary match") +
  facet_wrap(~lab, scales = "free_x", ncol = 4) +
  labs(x = "Regional SV2A measure (as provided)", y = "Mouse rCMRgl (\u00b5mol/100 g/min)") +
  theme_bw(base_size = 10) + theme(strip.text = element_text(size = 7))
nf <- ceiling(n_distinct(pd$analysis_set) / 4)
ggsave(nh("figs/mouse_s1c/synaptic_density_vs_rcmr.pdf"), p, width = 13, height = 3 * nf + 0.5)
ggsave(nh("figs/mouse_s1c/synaptic_density_vs_rcmr.jpg"), p, width = 13, height = 3 * nf + 0.5, dpi = 150)
