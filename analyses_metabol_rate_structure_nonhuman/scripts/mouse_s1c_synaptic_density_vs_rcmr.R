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
if (anyNA(syn$synaptic_density) || !is.numeric(syn$synaptic_density)) stop("synaptic_density must be numeric and complete.", call. = FALSE)
if (anyDuplicated(syn$cmrgl_region)) stop("One row per cmrgl_region is required; average sub-regions before saving.", call. = FALSE)
if (length(syn_u <- unique(syn$unit)) != 1L || length(unique(syn$measure)) != 1L)
  stop("A single measure and unit per file is required.", call. = FALSE)

rates <- read_mouse_rates()
if (length(m <- setdiff(syn$cmrgl_region, rates$cmrgl_region)))
  stop("Region(s) without a mouse rate: ", paste(m, collapse = ", "), call. = FALSE)
d <- inner_join(syn, rates, by = "cmrgl_region")
write_csv(d, nh("tables/mouse_s1c/matched_regions.csv"))
n <- nrow(d); message("Matched ", n, " regions")
if (n < 4) stop("Fewer than 4 matched regions; no inference attempted.", call. = FALSE)

ct <- suppressWarnings(cor.test(d$rcmr_value, d$synaptic_density, method = "spearman", exact = FALSE))
fit <- lm(rcmr_value ~ synaptic_density, data = d); sm <- summary(fit)
fits <- tibble(n = n, measure = unique(d$measure), unit = syn_u,
               spearman_rho = unname(ct$estimate), spearman_p = ct$p.value,
               lm_slope = coef(fit)[[2]], lm_r2 = sm$r.squared, lm_p = sm$coefficients[2, 4])
write_csv(fits, nh("tables/mouse_s1c/fits.csv")); print(as.data.frame(fits), digits = 3)

loo <- bind_rows(lapply(seq_len(n), function(i) {
  z <- d[-i, ]; r <- suppressWarnings(cor.test(z$rcmr_value, z$synaptic_density, method = "spearman", exact = FALSE))
  tibble(left_out = d$cmrgl_region[i], n = n - 1, spearman_rho = unname(r$estimate), spearman_p = r$p.value)
}))
write_csv(loo, nh("tables/mouse_s1c/leave_one_out.csv"))
message("Leave-one-out rho range: ", sprintf("%.2f to %.2f", min(loo$spearman_rho), max(loo$spearman_rho)))

p <- ggplot(d, aes(synaptic_density, rcmr_value)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "steelblue", linewidth = 0.6) +
  geom_point(size = 2.6) + geom_text(aes(label = cmrgl_region), vjust = -0.8, size = 3) +
  labs(x = sprintf("%s (%s)", unique(d$measure), syn_u), y = "Mouse rCMRgl (µmol/100 g/min)",
       subtitle = sprintf("Spearman rho = %.2f, p = %.3f, n = %d", fits$spearman_rho, fits$spearman_p, n)) +
  theme_bw(base_size = 11)
ggsave(nh("figs/mouse_s1c/synaptic_density_vs_rcmr.pdf"), p, width = 6.5, height = 4.8)
ggsave(nh("figs/mouse_s1c/synaptic_density_vs_rcmr.jpg"), p, width = 6.5, height = 4.8, dpi = 200)
