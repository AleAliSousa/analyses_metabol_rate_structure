############################################################
# mouse_s1b_5 -- Cortical excitatory:inhibitory ratio vs rate across
# neocortical areas. Mouse analogue of s1b_5 (Jorstad-like neocortical E:I).
#   E = Allen classes 01 IT-ET Glut + 02 NP-CT-L6b Glut
#   I = classes 06 CTX-CGE GABA + 07 CTX-MGE GABA
# Computed per donor x ROI (>= MIN_CELLS cortical neurons), then averaged.
# Neocortical areas with a rate are the primary-included isocortex ROIs.
# Also writes the S1 (SSp) profile against the other isocortical ROIs.
# Outputs: tables/mouse_s1b/EI_*.csv, figs/mouse_s1b/EI_ratio_vs_rcmr.*
############################################################
.a <- commandArgs(FALSE); .f <- sub("^--file=", "", .a[grepl("^--file=", .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else "analyses_metabol_rate_structure_nonhuman/scripts", "common.R"))
source(nh("scripts/mouse_s1b_helpers.R"))

counts <- read_mouse_counts(); roi_map <- read_roi_map(); rates <- read_mouse_rates()
E <- c("01 IT-ET Glut", "02 NP-CT-L6b Glut"); I <- c("06 CTX-CGE GABA", "07 CTX-MGE GABA")
iso <- roi_map %>% filter(is_isocortex) %>% pull(allen_roi)

ei <- counts %>% filter(roi %in% iso, class %in% c(E, I)) %>%
  mutate(EI = if_else(class %in% E, "E", "I")) %>%
  group_by(roi, donor, EI) %>% summarise(n = sum(n_cells), .groups = "drop") %>%
  pivot_wider(names_from = EI, values_from = n, values_fill = 0) %>%
  filter(E + I >= MIN_CELLS_DONOR_ROI, I > 0) %>%
  mutate(EI_ratio = E / I, E_frac = E / (E + I))

by_roi <- ei %>% group_by(roi) %>%
  summarise(n_donors = n(), n_E = sum(E), n_I = sum(I),
            EI_ratio = mean(EI_ratio), E_frac = mean(E_frac),
            EI_ratio_pooled = sum(E) / sum(I), .groups = "drop") %>%
  left_join(select(roi_map, roi = allen_roi, cmrgl_region, primary_included), by = "roi")
write_csv(by_roi, nh("tables/mouse_s1b/EI_ratio_by_isocortical_roi.csv"))

# S1 profile: SSp vs mean of the other isocortical ROIs
s1 <- by_roi %>% mutate(group = if_else(roi == "SSp", "S1 (SSp)", "other isocortex")) %>%
  group_by(group) %>% summarise(n_rois = n(), EI_ratio = mean(EI_ratio), E_frac = mean(E_frac), .groups = "drop")
write_csv(s1, nh("tables/mouse_s1b/EI_S1_vs_other_isocortex.csv")); print(s1)

by_region <- by_roi %>% filter(primary_included, !is.na(cmrgl_region)) %>%
  group_by(cmrgl_region) %>% summarise(EI_ratio = mean(EI_ratio), E_frac = mean(E_frac), n_donors = sum(n_donors),
                                       rois = paste(roi, collapse = "+"), .groups = "drop") %>%
  inner_join(rates, by = "cmrgl_region")
write_csv(by_region, nh("tables/mouse_s1b/EI_ratio_by_region_with_rcmr.csv"))
st <- spearman_table(by_region, c("EI_ratio", "E_frac")); write_csv(st, nh("tables/mouse_s1b/EI_spearman_vs_rcmr.csv"))
print(by_region); print(st)

p <- ggplot(by_region, aes(EI_ratio, rcmr_value)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "steelblue", linewidth = 0.6) +
  geom_point(size = 2.6) + geom_text(aes(label = cmrgl_region), vjust = -0.8, size = 3) +
  labs(x = "Cortical E:I ratio (IT-ET + NP-CT-L6b : CGE + MGE), mean across donors",
       y = "Mouse rCMRgl (µmol/100 g/min)",
       subtitle = sprintf("Spearman rho = %.2f, p = %.3f, n = %d", st$rho[st$predictor == "EI_ratio"],
                          st$p_value[st$predictor == "EI_ratio"], st$n[st$predictor == "EI_ratio"])) +
  theme_bw(base_size = 11)
ggsave(nh("figs/mouse_s1b/EI_ratio_vs_rcmr.pdf"), p, width = 6, height = 4.5)
ggsave(nh("figs/mouse_s1b/EI_ratio_vs_rcmr.jpg"), p, width = 6, height = 4.5, dpi = 200)
