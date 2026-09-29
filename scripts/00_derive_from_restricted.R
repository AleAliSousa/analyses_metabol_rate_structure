# =====================================================================
# 00_derive_from_restricted.R
#
# OPT-IN provenance layer (peer-review reproducibility). Regenerates the
# PUBLIC group-mean input files FROM the restricted raw sources that
# 00_sync_evo_m1_inputs-restricted.R places in the gitignored data_restricted/,
# then verifies the regenerated values match the committed public files in
# data_raw/. Runs only when the restricted raw is present; the default public
# pipeline never needs it.
#
# STATUS:
#   [DONE]  Mackes cortical    -> Mackes_etal_2020_ERABIS_derived.csv
#   [DONE]  Mackes subcortical -> Mackes_etal_2020_subcortical_derived.csv
#   [DONE]  KP control + DS    -> Karlsen_Pakkenberg_2011_authordata_groupmeans.csv
# All four public group-mean files regenerate from the restricted raw and are
# verified against the committed public copies.
# =====================================================================

suppressMessages(library(readxl))

# ---- locate repo root (portable: Rscript or interactive) ------------
.script_args <- commandArgs(trailingOnly = FALSE)
.file_arg <- .script_args[grepl("^--file=", .script_args)]
.start_dir <- if (length(.file_arg)) {
  .p <- sub("^--file=", "", .file_arg[1]); .p <- gsub("~+~", " ", .p, fixed = TRUE)
  dirname(normalizePath(.p, mustWork = FALSE))
} else getwd()
.repo_root <- normalizePath(.start_dir, mustWork = FALSE)
while (!file.exists(file.path(.repo_root, ".git")) &&
       dirname(.repo_root) != .repo_root) .repo_root <- dirname(.repo_root)
if (!file.exists(file.path(.repo_root, ".git")))
  stop("Could not locate the project root from: ", .start_dir, call. = FALSE)
setwd(.repo_root)

# Compare a freshly derived table to the committed public file and report.
verify_against_public <- function(derived, public_path, keys, tol = 0.05) {
  if (!file.exists(public_path)) {
    message("  (no public reference at ", public_path, " to verify against)")
    return(invisible(NULL))
  }
  pub <- read.csv(public_path, check.names = FALSE, stringsAsFactors = FALSE)
  key <- function(df) do.call(paste, c(df[keys], sep = "||"))
  dk <- key(derived); pk <- key(pub)
  num_cols <- intersect(names(derived), names(pub))
  num_cols <- num_cols[vapply(num_cols, function(c)
    is.numeric(derived[[c]]) || suppressWarnings(!all(is.na(as.numeric(pub[[c]])))),
    logical(1))]
  mism <- 0L
  for (i in seq_len(nrow(derived))) {
    j <- match(dk[i], pk); if (is.na(j)) { cat("  [missing in public]", dk[i], "\n"); mism <- mism + 1L; next }
    for (cc in setdiff(num_cols, keys)) {
      a <- suppressWarnings(as.numeric(derived[[cc]][i]))
      b <- suppressWarnings(as.numeric(pub[[cc]][j]))
      if (is.na(a) && is.na(b)) next
      if (xor(is.na(a), is.na(b)) || abs(a - b) > tol) {
        cat(sprintf("  [DIFF] %-28s %-10s derived=%s public=%s\n", dk[i], cc, a, b)); mism <- mism + 1L
      }
    }
  }
  if (mism == 0L) message("  VERIFIED: regenerated values match the public file (tol=", tol, ").")
  else message("  ", mism, " difference(s) vs public file -- review above.")
  invisible(mism)
}

# =====================================================================
# Mackes cortical (ERABIS) : raw psych::describe workbook -> ERABIS_derived
# =====================================================================
derive_mackes_cortical <- function(xlsx_path) {
  sheet_group <- c("UK adoptees" = "UK", "Romanian adoptees" = "Romanian")
  # Cells carry trailing non-breaking spaces (\u00a0) that trimws() leaves;
  # normalise all unicode whitespace before parsing.
  clean <- function(x) {
    x <- as.character(x)
    x <- gsub("[\u00a0\u2007\u202f]", " ", x)
    trimws(x)
  }
  numv <- function(x) as.numeric(clean(x))
  out <- list()
  for (sh in names(sheet_group)) {
    d <- suppressMessages(read_excel(xlsx_path, sheet = sh, skip = 2,
                                     col_names = FALSE))
    region_raw_full <- clean(d[[1]])
    keep <- !is.na(region_raw_full) & region_raw_full != ""
    d <- d[keep, ]; region_raw_full <- region_raw_full[keep]
    # strip trailing "_volume"; keep summary rows (BrainSegVolNotVent.x, eTIV.x)
    region_raw <- sub("_volume$", "", region_raw_full)
    hemi <- ifelse(grepl("^lh_", region_raw), "L",
             ifelse(grepl("^rh_", region_raw), "R", NA_character_))
    structure <- sub("^(lh_|rh_)", "", region_raw)
    out[[sh]] <- data.frame(
      Group      = unname(sheet_group[sh]),
      Region_raw = region_raw,
      Structure  = structure,
      Hemisphere = hemi,
      n          = as.integer(round(numv(d[[3]]))),
      mean.mm3   = round(numv(d[[4]]), 3),
      sd.mm3     = round(numv(d[[5]]), 3),
      se.mm3     = round(numv(d[[14]]), 3),
      stringsAsFactors = FALSE, check.names = FALSE
    )
  }
  do.call(rbind, out)
}

erabis_raw <- "data_restricted/Mackes_etal_2020_ERABIS_volume_stats_snapshot.xlsx"
if (file.exists(erabis_raw)) {
  message("Mackes cortical: deriving ERABIS_derived from ", erabis_raw)
  erabis <- derive_mackes_cortical(erabis_raw)
  target <- "Mackes_etal_2020_ERABIS_derived.csv"          # basename of public file
  # locate the public copy currently in data_raw/ (synced by public manifest)
  public_in_dataraw <- file.path("data_raw", target)
  verify_against_public(erabis, public_in_dataraw, keys = c("Group", "Region_raw"))
  write.csv(erabis, public_in_dataraw, row.names = FALSE, na = "")
  message("  wrote ", public_in_dataraw, " (", nrow(erabis), " rows) regenerated from restricted raw.")
} else {
  message("Mackes cortical: restricted raw not present (", erabis_raw,
          "); skipping. Run 00_sync_evo_m1_inputs-restricted.R with the private repo mounted.")
}

# =====================================================================
# Mackes subcortical : ERA descriptives CSV -> subcortical_derived
# Uses the raw bilateral mean_<structure> rows (not the *_reg TBV/sex-adjusted).
# =====================================================================
derive_mackes_subcortical <- function(csv_path) {
  lines <- readLines(csv_path, warn = FALSE)
  struct_map <- c(thalamus = "Thalamus", caudate = "Caudate",
                  putamen = "Putamen", pallidum = "Pallidum",
                  hippocampus = "Hippocampus", amygdala = "Amygdala",
                  accumbens = "Nucleus Accumbens")
  region_order <- c("Thalamus", "Caudate", "Putamen", "Pallidum",
                    "Hippocampus", "Amygdala", "Nucleus Accumbens")
  group <- NA_character_; rows <- list()
  for (ln in lines) {
    if (grepl("UK adoptees", ln)) group <- "UK"
    else if (grepl("Romanian adoptees", ln)) group <- "Romanian"
    f <- strsplit(ln, ",")[[1]]
    key <- trimws(f[1])
    if (grepl("^mean_[a-z]+$", key)) {           # excludes mean_*_reg
      st <- sub("^mean_", "", key)
      if (st %in% names(struct_map) && !is.na(group)) {
        rows[[paste(group, st)]] <- data.frame(
          Region = struct_map[[st]], Group = group,
          n = as.integer(f[2]), taxon = "Homo sapiens",
          Volume_per_hemisphere.cm3 = round(as.numeric(f[3]), 4),
          SD.cm3 = round(as.numeric(f[4]), 4),
          stringsAsFactors = FALSE, check.names = FALSE)
      }
    }
  }
  df <- do.call(rbind, rows)
  df$Group  <- factor(df$Group, levels = c("UK", "Romanian"))
  df$Region <- factor(df$Region, levels = region_order)
  df <- df[order(df$Group, df$Region), ]
  df$Group <- as.character(df$Group); df$Region <- as.character(df$Region)
  rownames(df) <- NULL
  df
}

subcort_raw <- "data_restricted/Mackes_etal_2020_subcortical_snapshot.csv"
if (file.exists(subcort_raw)) {
  message("Mackes subcortical: deriving subcortical_derived from ", subcort_raw)
  subc <- derive_mackes_subcortical(subcort_raw)
  public_in_dataraw <- file.path("data_raw", "Mackes_etal_2020_subcortical_derived.csv")
  verify_against_public(subc, public_in_dataraw, keys = c("Group", "Region"))
  write.csv(subc, public_in_dataraw, row.names = FALSE, na = "")
  message("  wrote ", public_in_dataraw, " (", nrow(subc), " rows) regenerated from restricted raw.")
} else {
  message("Mackes subcortical: restricted raw not present (", subcort_raw, "); skipping.")
}

# =====================================================================
# Karlsen & Pakkenberg : raw counting workbooks -> authordata_groupmeans
# Parsing mirrors the authors' Karlsen_Pakkenberg_2011_authordata.R (Summ sheet,
# per-brain columns under the 'Code' row; cell NUMBERS x2 = bilateral, densities
# and glia:neuron ratio left as recorded). Per-brain values are then aggregated
# to Group x Region means (Control = 6 brains, DS = 4 brains).
# =====================================================================
kp_regions <- data.frame(
  Region = c("Frontal", "Temporal", "Parietal", "Occipital", "Neocortex total"),
  npref  = c("F", "T", "P", "O", "Ne"),        # neuron / volume prefix
  gpref  = c("Fr", "Te", "Pa", "Oc", "Ne"),    # glia-subtype prefix
  stringsAsFactors = FALSE)

kp_parse_file <- function(path, group) {
  M <- as.matrix(suppressWarnings(
    read_excel(path, sheet = "Summ", col_names = FALSE, .name_repair = "minimal")))
  M[is.na(M)] <- ""
  labs <- trimws(paste(M[, 1], M[, 2], M[, 3]))
  code_row <- which(apply(M, 1, function(r) any(r == "Code")) &
                      seq_len(nrow(M)) > 40)[1]
  code_cols <- which(nzchar(M[code_row, ]) & M[code_row, ] != "0" &
                       M[code_row, ] != "Code")
  code_cols <- code_cols[code_cols >= 4]
  codes <- M[code_row, code_cols]
  getrow <- function(label) {
    i <- which(labs == label)[1]
    if (is.na(i)) return(rep(NA_real_, length(code_cols)))
    suppressWarnings(as.numeric(M[i, code_cols]))
  }
  out <- list()
  for (k in seq_len(nrow(kp_regions))) {
    np <- kp_regions$npref[k]; gp <- kp_regions$gpref[k]; rg <- kp_regions$Region[k]
    out[[k]] <- data.frame(
      Group = group, Brain = codes, Region = rg, taxon = "Homo sapiens",
      Volume_hemisphere.cm3 = getrow(paste0(np, " V cm3")),
      Neuron_N.x10.9 = getrow(paste0(np, " N 10^9")) * 2,
      Neuron_density.x10.6.per.cm3 = getrow(paste0(np, " NV Mcm-3")),
      Astrocyte_N.x10.9 = getrow(paste0(gp, "Ast N")) * 2,
      Astrocyte_density.x10.6.per.cm3 = getrow(paste0(gp, "Ast NV")),
      Oligodendrocyte_N.x10.9 = getrow(paste0(gp, "Oli N")) * 2,
      Oligodendrocyte_density.x10.6.per.cm3 = getrow(paste0(gp, "Oli NV")),
      Microglia_N.x10.9 = getrow(paste0(gp, "Mic N")) * 2,
      Microglia_density.x10.6.per.cm3 = getrow(paste0(gp, "Mic NV")),
      Glia_N.x10.9 = getrow(paste0(gp, "_Gl N")) * 2,
      Glia_density.x10.6.per.cm3 = getrow(paste0(gp, "_Gl NV")),
      Glia_Neuron_ratio = getrow(paste0(gp, "_Gl Gl/Ne")),
      stringsAsFactors = FALSE, check.names = FALSE)
  }
  do.call(rbind, out)
}

kp_ctrl <- "data_restricted/Karlsen_Pakkenberg_2011_authordata_control_snapshot.xls"
kp_ds   <- "data_restricted/Karlsen_Pakkenberg_2011_authordata_DS_snapshot.xls"
if (file.exists(kp_ctrl) && file.exists(kp_ds)) {
  message("Karlsen & Pakkenberg: deriving authordata_groupmeans from control + DS workbooks")
  dat <- rbind(kp_parse_file(kp_ctrl, "Control"), kp_parse_file(kp_ds, "DS"))
  dat <- dat[!is.na(dat$Brain) & nzchar(dat$Brain), ]
  measure_cols <- setdiff(names(dat), c("Group", "Brain", "Region", "taxon"))
  agg <- aggregate(dat[measure_cols],
                   by = list(Group = dat$Group, Region = dat$Region),
                   FUN = function(x) mean(x, na.rm = TRUE))
  nbr <- tapply(dat$Brain, dat$Group, function(x) length(unique(x)))
  agg$n_brains <- as.integer(nbr[agg$Group])
  for (cc in measure_cols) agg[[cc]] <- round(agg[[cc]], 4)
  agg$Region <- factor(agg$Region, levels = kp_regions$Region)
  agg$Group  <- factor(agg$Group, levels = c("Control", "DS"))
  agg <- agg[order(agg$Region, agg$Group), ]
  agg$Region <- as.character(agg$Region); agg$Group <- as.character(agg$Group)
  col_order <- c("Group", "Region", "Volume_hemisphere.cm3",
    "Neuron_N.x10.9", "Neuron_density.x10.6.per.cm3",
    "Astrocyte_N.x10.9", "Astrocyte_density.x10.6.per.cm3",
    "Oligodendrocyte_N.x10.9", "Oligodendrocyte_density.x10.6.per.cm3",
    "Microglia_N.x10.9", "Microglia_density.x10.6.per.cm3",
    "Glia_N.x10.9", "Glia_density.x10.6.per.cm3", "Glia_Neuron_ratio", "n_brains")
  agg <- agg[, col_order]
  public_in_dataraw <- file.path("data_raw",
    "Karlsen_Pakkenberg_2011_authordata_groupmeans.csv")
  verify_against_public(agg, public_in_dataraw, keys = c("Group", "Region"))
  write.csv(agg, public_in_dataraw, row.names = FALSE, na = "")
  message("  wrote ", public_in_dataraw, " (", nrow(agg), " rows) regenerated from restricted raw.")
} else {
  message("KP: control/DS workbooks not present in data_restricted/; skipping.")
}
