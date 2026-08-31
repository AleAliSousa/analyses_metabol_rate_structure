# =============================================================================
# build_stephan_primates_reference_sheet.R
#
# PURPOSE
#   The current metadata (metadata/Stephan_primates_metadata.xlsx) records
#   references grouped PER VARIABLE (a column can list several candidate
#   sources, e.g. "Stephan et al 1981; Zilles and Rehkamper 1988; Bauernfeind
#   et al 2013"). This script resolves, for EACH DATA POINT (species x column),
#   the SPECIFIC source table that the value came from, by value-matching every
#   cell of data_raw/Stephan_primates.csv against the per-source provenance
#   table exported from the organised Evo-M1-Trait-Data project
#   (the per-source ledger: Species x Variable x Value x Source-table x Year).
#   Cells without a near-perfect direct match then go through a broader,
#   ordered source search: metadata-listed papers (all of their tables), every
#   final volume table in Evo-M1, then pre-2020 inputs in the restricted repo.
#
# INPUTS  (all read relative to the project root set below)
#   LOCAL to this repo (curated here; no upstream copy exists)
#     data_raw/Stephan_primates.csv                   - the compiled data
#     metadata/Stephan_primates_metadata.xlsx         - per-variable references
#     metadata/qc_stephan/stephan_wide_crosswalk.csv  - anatomy/unit keys
#
#   READ LIVE from Evo-M1-Trait-Data (no local copy is kept, by design)
#     __merging_volumes/volumes_unfiltered_select.csv - per-source values
#     __merging_volumes/volumes_flags_select.csv      - merge deviation flags
#         (outputs of volumes_compiled_select.R, the ONE merge instance this
#          project uses -- see EVOM1_DIR / MERGE_INSTANCE below.)
#     Stephan_etal_1981/Stephan_etal_1981_TablesI-VI.csv - digitised 1981 tables
#         (output of Stephan_etal_1981_TablesI-VI.R in that folder.)
#     __merging_volumes/volumes_unfiltered.csv        - canonical inventory of
#         all final brain-volume tables, used only for fallback matching.
#     deSousa_etal_2010/deSousa_etal_2010_SupTable2.csv and
#     Smaers_etal_2017/Smaers_etal_2017_TableS1part1.csv - final public tables
#         named by the metadata but intentionally not selected into the ledger.
#   READ, WHEN MOUNTED, from Evo-M1-Trait-Data-restricted
#     MIGRATED_INDEX.csv plus the known unpublished volume workbooks. Restricted
#     files are searched only for unresolved cells and never copied here.
#
# OUTPUTS (written to metadata/)
#   Stephan_primates_references_long.csv        - one row per non-empty data
#                                                 cell, with the resolved
#                                                 specific reference(s).
#   Stephan_primates_reference_mismatches.csv   - cells whose value does NOT
#                                                 match any Evo-M1 source
#                                                 (potential raw / pre- or
#                                                 post-publication data) plus
#                                                 cells matched only after a
#                                                 unit rescale. FOR REVIEW.
#   Stephan_primates_references_by_column.csv   - per-column coverage summary.
#   Stephan_primates_candidate_source_inventory.csv - every public/restricted
#                                                 source searched and its date.
#
# MATCHING RULES
#   * A Stephan value matches an Evo-M1 source value if they agree within 1%
#     (relative) OR agree after rounding the source value to the Stephan
#     value's decimal places (handles 2-sig-fig rounding of small volumes) OR
#     both are zero.
#   * For body/brain-mass columns only, matches are also tried at x1000 and
#     x0.001 (kg/g, g/mg) and the rescale is flagged.
#   * Vestibular-complex columns are matched against the *_unilateral_* Evo-M1
#     variants, because the Stephan file reports one-side volumes (Stephan et
#     al 1981 Table XIII), not the bilateral Baron et al 1988 values.
#   * Structure keys come from stephan_wide_crosswalk.csv. A numeric match to a
#     different structure is never promoted merely because it occurs in a
#     metadata-listed paper. Source-aware alternatives are explicit: Pons is
#     ventral pons for Matano 1985b, but generic pons for the Zilles orangutan.
#   * When several sources carry the same value, the "preferred" reference is
#     the one consistent with the metadata's stated source(s) for that column,
#     then a table already known to contribute other cells. Publication year
#     must be <= 2019 because later documents cannot have contributed to the
#     Stephan_primates compilation.
#
# NOTE: this script only READS the source data; it never edits
#       Stephan_primates.csv. Decide per mismatch whether raw data should be
#       adopted, then change it in both projects by hand.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
})

## --- project root -----------------------------------------------------------
## Portable: walk up from the current directory to the repo root (.git).
## All paths below are relative to this. See helpers/project_root.R.
setwd(local({ d <- normalizePath(getwd()); while (!file.exists(file.path(d, ".git")) && dirname(d) != d) d <- dirname(d); d }))

## --- which Evo-M1 merge instance ---------------------------------------------
## __merging_volumes/ holds several merge instances, each written by its own
## script and distinguished by a filename suffix:
##   volumes_compiled_select.R    -> volumes_*_select.csv     <- THE ONE WE USE
##   volumes_compiled.R           -> volumes_*.csv            (no suffix)
##   volumes_compiled_DeCasien.R  -> volumes_*_DeCasien.csv
## The project has standardised on the _select instance; the others are not
## used anywhere in this repo. Everything downstream must come from a single
## instance, so change this in one place if that ever moves.
MERGE_INSTANCE <- "select"
sfx <- if (nzchar(MERGE_INSTANCE)) paste0("_", MERGE_INSTANCE) else ""

## A source published after the compilation existed cannot be its source.
## The compilation was assembled in 2019; include papers from 2019 but exclude
## 2020 onward. This is a publication/source-date cutoff, not a filesystem mtime
## cutoff (OneDrive migrations make mtimes unsuitable evidence).
SOURCE_CUTOFF_YEAR <- 2019L
SEARCH_RESTRICTED <- TRUE

## Read the merge outputs STRAIGHT FROM SOURCE, not from a copy in
## data_intermediate/. Same convention as compare_stephan_vs_volumes_wide.R,
## compare_volumes_files.R and the s3 select script. Rationale: a local copy has
## no script that maintains it, so it drifts from the merge silently, and
## sitting in data_intermediate/ it looks like a regenerable build product and
## gets swept up by cleanups. Reading source means the run either sees the
## current merge or stops.
EVOM1_ROOT <- file.path("/Users/crossmodal/Library/CloudStorage",
                        "OneDrive-AllenInstitute/Species/Evo-M1-Trait-Data")
EVOM1_DIR  <- file.path(EVOM1_ROOT, "__merging_volumes")   # the merge instances
EVOM1_S81  <- file.path(EVOM1_ROOT, "Stephan_etal_1981")   # 1981 table extraction

## Only these two are genuinely local -- curated here, with no upstream copy to
## read. Everything else comes from Evo-M1-Trait-Data.
f_data  <- "data_raw/Stephan_primates.csv"
f_meta  <- "metadata/Stephan_primates_metadata.xlsx"
f_xwalk <- "metadata/qc_stephan/stephan_wide_crosswalk.csv"

f_prov  <- file.path(EVOM1_DIR, sprintf("volumes_unfiltered%s.csv", sfx))
f_flags <- file.path(EVOM1_DIR, sprintf("volumes_flags%s.csv", sfx))
## Canonical all-table ledger. This is the maintained inventory of final
## Evo-M1 brain-volume tables and is broader than the hand-selected ledger.
f_prov_all <- file.path(EVOM1_DIR, "volumes_unfiltered.csv")
## Written by Stephan_etal_1981/Stephan_etal_1981_TablesI-VI.R, which digitises
## the printed tables. Read from source for the same reason as the merge files.
f_s81   <- file.path(EVOM1_S81, "Stephan_etal_1981_TablesI-VI.csv")
f_desousa_supp2 <- file.path(EVOM1_ROOT, "deSousa_etal_2010",
                             "deSousa_etal_2010_SupTable2.csv")
f_smaers_s1 <- file.path(EVOM1_ROOT, "Smaers_etal_2017",
                         "Smaers_etal_2017_TableS1part1.csv")

## The restricted repo normally sits beside the OneDrive "Species" directory;
## EVOM1_RESTRICTED remains the portable override used by the Evo-M1 projects.
restricted_default <- file.path(dirname(dirname(EVOM1_ROOT)),
                                "Evo-M1-Trait-Data-restricted")
EVOM1_RESTRICTED_ROOT <- Sys.getenv("EVOM1_RESTRICTED", unset = restricted_default)

## Outputs live under the qc_stephan/ layer: scaffolding for checking over
## Stephan_primates.csv, meant to be discarded as a unit once that file is
## replaced with a fresh version. Repo-root-relative, so this works from
## wherever the script is invoked.
out_dir <- "metadata/qc_stephan"

## In the _select ledger every row carries Team = "Selected_collection"; the
## unsuffixed instance used "Stephan_collection". The preference rule below
## prefers rows from this collection, so the name must track the instance.
COLLECTION_TEAM <- if (identical(MERGE_INSTANCE, "select"))
                     "Selected_collection" else "Stephan_collection"

## --- required inputs --------------------------------------------------------
## HARD dependency by design: this script assigns a citation to every published
## data point, so a missing input must stop the run rather than quietly produce
## a thinner reference sheet. In particular the flags file and the 1981 tables
## used to be optional; without them the script still ran but silently lost the
## merge-deviation flags and the Stephan-1981 fallback confirmations, changing
## the resolved references with no warning. They are required now.
##
## Two origins, and they fail and repair differently:
##   * Evo-M1 products, read live from OneDrive. Missing almost always means
##     OneDrive is not synced, or the upstream script has not been run -- it
##     does NOT mean anything was lost, so never "restore" these from git.
##   * Stephan_primates.csv and its metadata sheet -- curated in this repo, no
##     upstream copy exists. If these go missing, restore from git.
required_inputs <- c(f_data, f_meta, f_xwalk, f_prov, f_flags, f_s81, f_prov_all,
                     f_desousa_supp2, f_smaers_s1)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) {
  evo_missing   <- missing_inputs[startsWith(missing_inputs, EVOM1_ROOT)]
  local_missing <- setdiff(missing_inputs, evo_missing)
  stop("build_stephan_primates_reference_sheet.R: missing required input(s):\n",
       paste0("  - ", missing_inputs, collapse = "\n"),
       if (length(evo_missing))
         paste0("\n\nRead live from Evo-M1-Trait-Data (volumes_* come from ",
                "volumes_compiled_select.R;\nthe 1981 tables from ",
                "Stephan_etal_1981_TablesI-VI.R). Is OneDrive mounted /\nsynced? ",
                "Has the upstream script been run? Do NOT restore these from git.")
       else "",
       if (length(local_missing))
         paste0("\n\nCurated in this repo -- no upstream copy exists. Restore from git:\n",
                paste0("  git checkout HEAD -- ", local_missing, collapse = "\n"))
       else "",
       call. = FALSE)
}

## Stamp which merge the run actually read, so a result can be traced back to a
## specific state of volumes_compiled_select.R.
message("Evo-M1 source: ", f_prov, "\n               (modified ",
        format(file.mtime(f_prov), "%Y-%m-%d %H:%M"), ")")
message("Source-date cutoff: publication year <= ", SOURCE_CUTOFF_YEAR)

## --- 1. species crosswalk (Stephan space-form -> Evo-M1 name) ---------------
## Label crosswalk only -- these assert nothing about taxonomy, they just map a
## Stephan_primates.csv label onto the label the merge happens to use. Verified
## against the _select ledger's 65 species; without them the cells look like
## provenance gaps.
sp_xwalk <- c(
  "Lagothrix lagotricha" = "Lagothrix lagothricha",  # spelling
  "Pongo pygmaeus"       = "Pongo sp.",              # merge pools orangutans
  "Pongo abelii"         = "Pongo sp.",              #   as "Pongo sp."
  "Gorilla gorilla"      = "Gorilla sp.",            # merge pools as "sp."
  "Gorilla gorilla gorilla" = "Gorilla sp.",          # Smaers 2017 trinomial
  "Pan troglodytes troglodytes" = "Pan troglodytes",  # Smaers 2017 trinomial
  "Tarsius syrichta"     = "Tarsius sp.",            #   likewise (cf. s81_sp)
  "Plecturocebus moloch" = "Callicebus moloch"       # merge keeps the old genus
)
norm_sp <- function(x) {
  x <- gsub("_", " ", trimws(x))
  hit <- x %in% names(sp_xwalk)
  x[hit] <- unname(sp_xwalk[x[hit]])
  x
}

## --- 2. column -> Evo-M1 variable crosswalk ---------------------------------
## Anatomy and unit mappings are data, not code. Keeping the build and the
## comparison scripts on this one CSV prevents their keys from drifting apart.
## Columns absent from it have no defensible Evo-M1 volume equivalent and are
## carried through as metadata_only.
xwalk <- read_csv(f_xwalk, show_col_types = FALSE) %>%
  filter(!is.na(stephan_var), stephan_var != "", !is.na(wide_var), wide_var != "")
if (anyDuplicated(xwalk$stephan_var))
  stop(f_xwalk, ": stephan_var must be unique", call. = FALSE)
if (any(is.na(suppressWarnings(as.numeric(xwalk$unit_scale))) |
        suppressWarnings(as.numeric(xwalk$unit_scale)) <= 0))
  stop(f_xwalk, ": unit_scale must be a positive number", call. = FALSE)
col2var <- setNames(xwalk$wide_var, xwalk$stephan_var)
col_unit_scale <- setNames(as.numeric(xwalk$unit_scale), xwalk$stephan_var)

## Pons is the one source-conditioned anatomy key. The Matano 1985b rows are
## explicitly ventral pons (VPo); the Pongo value is the broader pons printed
## by Zilles & Rehkamper. These are not global aliases.
structure_key_variants <- function(col) {
  base <- tibble(Variable = unname(col2var[[col]]), source_prefix = NA_character_,
                 key_rule = "crosswalk_exact")
  if (identical(col, "Pons")) bind_rows(
    tibble(Variable = "Ventral_pons_Vol.mm3",
           source_prefix = "Matano_etal_1985_b", key_rule = "Matano ventral pons"),
    tibble(Variable = "Pons_Vol.mm3",
           source_prefix = "Zilles_Rehkämper_1988", key_rule = "Zilles whole pons"),
    tibble(Variable = "Pons_Vol.mm3",
           source_prefix = "Zilles_Rehkämper_1988", key_rule = "Zilles whole pons"),
    tibble(Variable = "Pons_Vol.mm3",
           source_prefix = "Zilles__Rehkamper_1988", key_rule = "Zilles whole pons")
  ) else base
}
source_satisfies_key <- function(source, prefix) {
  is.na(prefix) || !nzchar(prefix) || startsWith(source, prefix)
}
candidate_key_row <- function(col, variable, source) {
  kv <- structure_key_variants(col)
  hit <- which(kv$Variable == variable &
                 vapply(kv$source_prefix, function(p)
                   source_satisfies_key(source, p), logical(1)))
  if (length(hit)) kv[hit[1], , drop = FALSE] else kv[0, , drop = FALSE]
}
match_scales_for_col <- function(col, broad = FALSE) {
  canonical <- 1 / unname(col_unit_scale[[col]])
  if (broad || col %in% mass_cols) unique(c(canonical, 1, 1000, 0.001)) else canonical
}
mass_cols <- c("Body_weight", "Brain_weight", "Body_mass", "Brain_mass", "Body_weight_1985a")

## data columns that carry a value but are not meant to resolve to a source
id_cols <- c("Species", "Species_Stephan_file", "Code_number_Stephan", "order",
             "Species_Matano1985a", "code_Matano1985a",
             "Number.of.individuals.vestibular.complex",
             "Number.of.individuals.periventrucular", "Number_cerebellar_nuclei")

## --- 3. citation library (source-table key prefix -> short citation) --------
ref_lib <- tribble(
  ~key,                 ~cite,
  "Stephan_etal_1970",  "Stephan, Bauchot & Andy 1970",
  "Stephan_etal_1981",  "Stephan, Frahm & Baron 1981 (Folia Primatol 35:1-29)",
  "Stephan_etal_1982",  "Stephan et al. 1982 (J Hirnforsch 23:575-591)",
  "Stephan_etal_1984",  "Stephan et al. 1984 (J Hirnforsch 25:385-403)",
  "Stephan_etal_1987",  "Stephan et al. 1987 (J Hirnforsch 28:571-584)",
  "Frahm_etal_1982",    "Frahm, Stephan & Stephan 1982 (J Hirnforsch 23:375-389)",
  "Frahm_etal_1984",    "Frahm et al. 1984 (J Hirnforsch 25:537-557)",
  "Frahm_etal_1997",    "Frahm et al. 1997 (J Hirnforsch 38:209-222)",
  "Frahm_etal_1998",    "Frahm et al. 1998 (J Hirnforsch 39:45-54)",
  "Frahm_Zilles_1994",  "Frahm & Zilles 1994 (J Hirnforsch 35:343-354)",
  "Baron_etal_1983",    "Baron et al. 1983 (J Hirnforsch 24:551-568)",
  "Baron_etal_1987",    "Baron et al. 1987 (J Hirnforsch 28:463-477)",
  "Baron_etal_1988",    "Baron et al. 1988 (J Hirnforsch 29:509-523)",
  "Baron_etal_1990",    "Baron et al. 1990 (J Hirnforsch 31:193-200)",
  "Matano_etal_1985",   "Matano et al. 1985 (Folia Primatol 44:182-203)",
  "deSousa_etal_2010",  "de Sousa et al. 2010 (J Hum Evol 58:281-292)",
  "deSousa_etal_2013",  "de Sousa et al. 2013 (Brain Behav Evol 81:93-108)",
  "Bauernfeind_etal_2013", "Bauernfeind et al. 2013 (J Hum Evol 64:263-279)",
  "Smaers_etal_2017",   "Smaers et al. 2017 (Curr Biol 27:714-720, Table S1)",
  "Sherwood_etal_2005", "Sherwood et al. 2005 (J Hum Evol 48:45-84)",
  "Zilles_Rehkämper_1988", "Zilles & Rehkamper 1988 (Table 12-2)",
  "Zilles__Rehkamper_1988", "Zilles & Rehkamper 1988 (Table 12-2)"
)
cite_of <- function(src) {
  if (is.na(src) || src == "") return(NA_character_)
  keys <- vapply(strsplit(src, ";\\s*")[[1]], function(s) {
    m <- ref_lib$key[startsWith(s, ref_lib$key)]
    if (length(m)) ref_lib$cite[match(m[which.max(nchar(m))], ref_lib$key)] else s
  }, character(1))
  paste(unique(keys), collapse = "; ")
}

## map a metadata citation phrase -> author_year key prefix used in Source IDs
meta_phrase_to_key <- function(phrase) {
  p <- tolower(phrase)
  case_when(
    grepl("bauernfeind", p)                 ~ "Bauernfeind_etal_2013",
    grepl("de sousa.*2013|sousa et al 2013", p) ~ "deSousa_etal_2013",
    grepl("de sousa.*2010|sousa et al 2010", p) ~ "deSousa_etal_2010",
    grepl("smaers.*2017|smaers et al 2017", p) ~ "Smaers_etal_2017",
    grepl("zilles", p)                       ~ "Zilles_Rehkämper_1988",
    grepl("matano", p)                       ~ "Matano_etal_1985",
    grepl("frahm.*199|frahm et al 1992", p)  ~ "Frahm_etal_1982",  # metadata "1992" is a typo for 1982
    grepl("frahm", p)                        ~ "Frahm_etal_1982",
    grepl("stephan et al 1981|stephan 1981", p) ~ "Stephan_etal_1981",
    grepl("stephan", p)                      ~ "Stephan_etal",
    TRUE                                     ~ NA_character_
  )
}

## --- 4. read inputs ---------------------------------------------------------
S <- read.csv(f_data, stringsAsFactors = FALSE, check.names = FALSE)
S <- S[!is.na(S$Species) & trimws(S$Species) != "", ]
S$Species_evo <- norm_sp(S$Species)

# metadata: keep rows where col-1 is a single-token variable name and col-4 has refs
meta_raw <- suppressMessages(read_excel(f_meta, sheet = "Stephan_NHprimates metadata",
                                         col_names = FALSE))
meta_map <- meta_raw %>%
  transmute(var = as.character(.[[1]]), refs = as.character(.[[4]])) %>%
  filter(!is.na(var), !is.na(refs), !grepl("\\s", var))       # variable names have no spaces
# align metadata variable names to the CSV column names (insula/brain cols drop "_L")
meta_lookup <- function(col) {
  v <- meta_map$refs[match(col, meta_map$var)]
  if (is.na(v)) v <- meta_map$refs[match(sub("_L$", "", col), meta_map$var)]
  if (is.na(v)) NA_character_ else v
}

P <- read.csv(f_prov, stringsAsFactors = FALSE, check.names = FALSE)
## Schema guard: the Evo-M1 merge is re-run independently of this repo, so a
## renamed column would otherwise surface as an opaque "undefined columns
## selected" from the split() below.
prov_required <- c("Species", "Variable", "Value", "Source", "Team", "Year")
if (length(setdiff(prov_required, names(P))))
  stop(f_prov, ": expected column(s) not found: ",
       paste(setdiff(prov_required, names(P)), collapse = ", "),
       "\nFound: ", paste(names(P), collapse = ", "),
       "\nThe Evo-M1 merge schema has changed -- update col2var / the reader.",
       call. = FALSE)
P$Value <- suppressWarnings(as.numeric(P$Value))
P$Year  <- suppressWarnings(as.integer(P$Year))
P <- P[!is.na(P$Value) & !is.na(P$Year) & P$Year <= SOURCE_CUTOFF_YEAR, ]
prov <- split(P[, c("Source", "Value", "Year", "Team")], list(P$Species, P$Variable), drop = TRUE)
prov_key <- function(sp, var) paste(sp, var, sep = ".")
direct_candidates <- function(sp, col) {
  kv <- structure_key_variants(col)
  out <- bind_rows(lapply(seq_len(nrow(kv)), function(j) {
    z <- prov[[prov_key(sp, kv$Variable[j])]]
    if (is.null(z)) return(NULL)
    z <- z[vapply(z$Source, source_satisfies_key, logical(1),
                  prefix = kv$source_prefix[j]), , drop = FALSE]
    if (!nrow(z)) return(NULL)
    mutate(as_tibble(z), Candidate_variable = kv$Variable[j],
           key_rule = kv$key_rule[j])
  }))
  if (nrow(out)) out else NULL
}

## Broader fallback ledger: every final source table wired into the canonical
## volume merge. Keep it separate from P so it cannot change direct matches.
P_all <- read.csv(f_prov_all, stringsAsFactors = FALSE, check.names = FALSE)
if (length(setdiff(prov_required, names(P_all))))
  stop(f_prov_all, ": expected column(s) not found: ",
       paste(setdiff(prov_required, names(P_all)), collapse = ", "), call. = FALSE)
P_all$Value <- suppressWarnings(as.numeric(P_all$Value))
P_all$Year  <- suppressWarnings(as.integer(P_all$Year))
P_all <- P_all[!is.na(P_all$Value) & !is.na(P_all$Year) &
                 P_all$Year <= SOURCE_CUTOFF_YEAR, prov_required]
P_all$Candidate_file <- f_prov_all

## Final public tables that are deliberately absent from the canonical ledger
## because the volume merge flags them as secondary/derived. They still matter
## here: these are the tables named by Stephan_primates metadata, and reference
## resolution asks which table supplied a value, not whether the merge selects
## that table for a modern pooled estimate. Values are normalized to Evo units.
public_final_long <- function(file, source, year, species_col, variable_map,
                              value_scale = 1) {
  d <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  missing <- setdiff(c(species_col, names(variable_map)), names(d))
  if (length(missing))
    stop(file, ": expected column(s) not found: ", paste(missing, collapse = ", "),
         call. = FALSE)
  bind_rows(lapply(names(variable_map), function(v) tibble(
    Species = norm_sp(d[[species_col]]), Variable = unname(variable_map[[v]]),
    Value = suppressWarnings(as.numeric(d[[v]])) * value_scale,
    Source = source, Team = "Public_final_table", Year = as.integer(year),
    Candidate_file = file
  ))) %>% filter(!is.na(Value))
}
P_public_extra <- bind_rows(
  public_final_long(
    f_desousa_supp2, "deSousa_etal_2010_SupTable2", 2010L, "species",
    c(V1_area_striata_volume_cm3 = "Area_striata_grey_matter_Vol.mm3",
      LGN_volume_cm3 = "Corpus_geniculatum_laterale_Vol.mm3",
      neocortex_volume_cm3 = "Neocortex_Vol.mm3"), 1000
  ),
  public_final_long(
    f_smaers_s1, "Smaers_etal_2017_TableS1part1", 2017L, "species",
    c(primary_visual_gray = "Area_striata_grey_matter_Vol.mm3",
      primary_visual_white = "Area_striata_white_matter_Vol.mm3",
      prefrontal_gray = "PrefrontalCortex_grey_matter_Vol.mm3",
      prefrontal_white = "PrefrontalCortex_white_matter_Vol.mm3",
      other_association_gray = "AssociationCortex_grey_matter_Vol.mm3",
      other_association_white = "AssociationCortex_white_matter_Vol.mm3",
      frontal_motor_gray = "FrontalMotor_grey_matter_Vol.mm3",
      frontal_motor_white = "FrontalMotor_white_matter_Vol.mm3"), 1000
  )
)
P_all <- bind_rows(P_all, P_public_extra) %>%
  distinct(Species, Variable, Value, Source, Team, Year, .keep_all = TRUE)

## merge deviation flags (newest-vs-next within Stephan_collection; possible typos)
## Path defined and existence-checked at the top; required, not optional.
FL <- read.csv(f_flags, stringsAsFactors = FALSE, check.names = FALSE)
dev_flag <- function(sp, var) {
  h <- which(FL$Species == sp & FL$Variable == var & FL$flag == "deviation")
  if (length(h)) FL$detail[h[1]] else ""
}
MASS_VARS <- c("Body_Mass.g", "Brain_Mass.mg")

## --- fallback provenance: Stephan et al 1981 combined Tables I-VI -----------
## Consulted only when the primary merge (volumes_unfiltered) yields no match.
## The merge drops some 1981 "revised" rows (e.g. Gorilla, Homo), keeping the
## older 1970 values, so those cells otherwise look like differences. Matching
## is value-gated, so the species synonym map below need not be exhaustive:
## wrong guesses simply fail to match and cost nothing.
## f_s81 defined and existence-checked at the top; required, not optional.
s81_sp <- c(
  Callithrix_pygmaea = "Cebuella pygmaea", Eulemur_fulvus = "Lemur fulvus",
  Varecia_variegata = "Lemur variegatus", Otolemur_crassicaudatus = "Galago crassicaudatus",
  Galagoides_demidoff = "Galago demidovii", Homo_sapiens = "Homo sapiens sapiens",
  Miopithecus_talapoin = "Cercopithecus talapoin", Lophocebus_albigena = "Cercocebus albigena",
  Piliocolobus_badius = "Colobus badius", Pithecia_monachus = "Pithecia monacha",
  Daubentonia_madagascariensis = "Daubentonia madagascar.", Plecturocebus_moloch = "Callicebus moloch",
  Avahi_laniger = "Avahi l. laniger", Avahi_occidentalis = "Avahi laniger occidentalis",
  Tarsius_syrichta = "Tarsius sp.", Cebus_albifrons = "Cebus sp.", Alouatta_seniculus = "Alouatta sp."
)
s81_col <- c(
  Body_weight = "Body_weight", Brain_weight = "Brain_weight", Ventricles = "Ventricles",
  Total_brain_net_volume = "Total_brain_net_volume", Medulla_oblongata = "Medulla_oblongata",
  Cerebellum = "Cerebellum", Mesencephalon = "Mesencephalon", Diencephalon = "Diencephalon",
  Telencephalon = "Telencephalon", Bulbus_olfactorius = "Bulbus_olfactorius",
  Bulbus_olfactorius_accessorius = "Bulbus_olfactorius_accessorius", Lobus_piriformis = "Lobus_piriformis",
  Septum = "Septum", Striatum = "Striatum", Schizocortex = "Schizo_cortex", Hippocampus = "Hippocampus",
  NeoWG = "Neocortex", Epithalamus = "Epithalamus", Thalamus = "Thalamus", Hypothalamus = "Hypothalamus",
  Subthalamus = "Subthalamus", Pallidum = "Pallidum", Nucleus_subthalamicus = "Nucleus_subthalamicus",
  Capsula_interna = "Capsula_interna", Tractus_opticus = "Tractus_opticus", Palaeocortex = "Palaeocortex",
  Amygdala = "Amygdala", Complexus_centromedialis = "Complexus_centromedialis",
  Nucleus_tractus_olfactorius = "Nucleus_tractus_olfactorius",
  Complexus_cortico_basolateralis = "Complexus_corticobasolateralis"
)
S81 <- read.csv(f_s81, stringsAsFactors = FALSE, check.names = FALSE)
s81_value <- function(sp_underscore, col) {
  if (is.null(S81) || !col %in% names(s81_col)) return(NA_real_)
  printed <- if (sp_underscore %in% names(s81_sp)) s81_sp[[sp_underscore]] else gsub("_", " ", sp_underscore)
  r <- which(S81$Species == printed); if (!length(r)) return(NA_real_)
  suppressWarnings(as.numeric(S81[[ s81_col[[col]] ]][r[1]]))
}

## intended source for columns not covered by the metadata sheet (Matano 1985a block)
col_primary <- c(
  Body_weight_1985a = "Matano_etal_1985_a", Cerebellar_nuclei_total = "Matano_etal_1985_a",
  Interpositus_cerebellar_nuclei = "Matano_etal_1985_a", Lateral_cerebellar_nuclei = "Matano_etal_1985_a",
  Medial_cerebellar_nuclei = "Matano_etal_1985_a"
)

## --- 5. matching helpers ----------------------------------------------------
rel_pct <- function(a, b) { d <- pmax(abs(a), abs(b)); ifelse(d == 0, 0, 100 * abs(a - b) / d) }
ndec <- function(s) { s <- trimws(s); if (grepl("\\.", s)) nchar(sub("^[^.]*\\.", "", s)) else 0L }
is_match <- function(sval, sstr, v, sc) {
  t <- v * sc
  if (sval == 0 && t == 0) return(TRUE)
  if (t != 0 && rel_pct(sval, t) <= 1) return(TRUE)
  nd <- ndec(sstr)
  if (sval != 0 && round(t, nd) == round(sval, nd)) return(TRUE)
  FALSE
}

## Source/table helpers used by the fallback pass. `source_paper()` removes the
## printed-table suffix while retaining same-author/year distinctions such as
## Matano_etal_1985_a vs _b.
ascii_key <- function(x) {
  y <- iconv(x %||% "", from = "", to = "ASCII//TRANSLIT", sub = "")
  gsub("[^a-z0-9]", "", tolower(y))
}
source_paper <- function(x) {
  sub("_(Table|TABLE|Tables|TABLES|SupplementaryTable|SupTable|SOMTable|Fig).*$",
      "", x)
}
source_is_allowed <- function(src, allow) {
  if (!length(allow)) return(rep(FALSE, length(src)))
  Reduce(`|`, lapply(allow, function(a) startsWith(src, a)), rep(FALSE, length(src)))
}
paper_is_allowed <- function(paper, allow) {
  if (!length(allow)) return(FALSE)
  pk <- ascii_key(paper)
  any(vapply(allow, function(a) {
    ak <- ascii_key(a)
    startsWith(pk, ak) || startsWith(ak, pk)
  }, logical(1)))
}

## Generic document reader for the deliberately looser restricted/local search.
## All sheets are read as text so old Excel workbooks and mixed-type rows retain
## their printed values. It returns NULL rather than aborting the main audit.
read_document <- function(f) {
  tryCatch({
    if (grepl("\\.csv$", f, ignore.case = TRUE)) {
      list(`_` = suppressWarnings(suppressMessages(
        readr::read_csv(f, col_names = FALSE, col_types = cols(.default = "c"),
                        locale = readr::locale(encoding = "latin1"), lazy = FALSE,
                        show_col_types = FALSE, progress = FALSE))))
    } else if (grepl("\\.tsv$", f, ignore.case = TRUE)) {
      list(`_` = suppressWarnings(suppressMessages(
        readr::read_tsv(f, col_names = FALSE, col_types = cols(.default = "c"),
                        locale = readr::locale(encoding = "latin1"), lazy = FALSE,
                        show_col_types = FALSE, progress = FALSE))))
    } else {
      sh <- readxl::excel_sheets(f)
      setNames(lapply(sh, function(s)
        suppressMessages(readxl::read_excel(f, sheet = s, col_names = FALSE,
                                             col_types = "text"))), sh)
    }
  }, error = function(e) NULL)
}

parse_doc_number <- function(x) {
  z <- trimws(as.character(x))
  z <- gsub(",", "", z, fixed = TRUE)
  ok <- grepl("^[+-]?([0-9]+\\.?[0-9]*|\\.[0-9]+)([eE][+-]?[0-9]+)?$", z)
  out <- rep(NA_real_, length(z))
  out[ok] <- suppressWarnings(as.numeric(z[ok]))
  out
}

is_match_vector <- function(sval, sstr, v, sc = 1) {
  t <- v * sc
  d <- pmax(abs(sval), abs(t))
  close <- ifelse(d == 0, TRUE, 100 * abs(sval - t) / d <= 1)
  rounded <- if (sval == 0) rep(FALSE, length(t)) else
    round(t, ndec(sstr)) == round(sval, ndec(sstr))
  !is.na(t) & (close | rounded)
}

## --- 6. walk every data cell ------------------------------------------------
rows <- list(); i <- 0L
for (r in seq_len(nrow(S))) {
  sp_raw <- S$Species[r]; sp_evo <- S$Species_evo[r]
  for (col in names(col2var)) {
    sstr <- as.character(S[[col]][r]); if (is.na(sstr) || trimws(sstr) == "") next
    sval <- suppressWarnings(as.numeric(sstr)); if (is.na(sval)) next
    var  <- col2var[[col]]
    meta_refs <- meta_lookup(col)
    cand <- direct_candidates(sp_evo, col)
    scales <- match_scales_for_col(col)
    m_src <- character(0); m_yr <- integer(0); m_team <- character(0)
    m_var <- character(0); m_rule <- character(0); facs <- numeric(0)
    if (!is.null(cand)) {
      for (k in seq_len(nrow(cand))) {
        for (sc in scales) {
          if (is_match(sval, sstr, cand$Value[k], sc)) {
            m_src <- c(m_src, cand$Source[k]); m_yr <- c(m_yr, cand$Year[k])
            m_team <- c(m_team, cand$Team[k]); m_var <- c(m_var, cand$Candidate_variable[k])
            m_rule <- c(m_rule, cand$key_rule[k]); facs <- c(facs, sc); break
          }
        }
      }
    }
    # confirm against the Stephan 1981 combined table (core Stephan columns): the
    # merge drops some 1981 revised rows (e.g. Gorilla, Homo). Value-gated.
    fv <- s81_value(sp_raw, col)
    fb81 <- !is.na(fv) && is_match(sval, sstr, fv, 1)
    if (fb81 && !any(startsWith(m_src, "Stephan_etal_1981"))) {
      m_src <- c(m_src, "Stephan_etal_1981_TablesI-VI"); m_yr <- c(m_yr, 1981L)
      m_team <- c(m_team, COLLECTION_TEAM)
      m_var <- c(m_var, var); m_rule <- c(m_rule, "Stephan 1981 printed-table confirmation")
      facs <- c(facs, 1 / unname(col_unit_scale[[col]]))
    }
    matched_src_u <- sort(unique(m_src))

    allow <- na.omit(unique(vapply(strsplit(meta_refs %||% "", ";\\s*")[[1]],
                                   meta_phrase_to_key, character(1))))
    if (col %in% names(col_primary)) allow <- c(allow, col_primary[[col]])

    # preferred reference: (1) scope to the metadata-listed source(s) (which
    # paper the compiler used); (2) within scope apply the volumes-merge rule:
    # COLLECTION_TEAM first, most recent year, mass -> Stephan 1981.
    pref <- NA_character_; pref_var <- var; pref_rule <- "crosswalk_exact"
    scope <- "none"; dev <- ""
    if (length(m_src)) {
      in_allow <- if (length(allow))
        Reduce(`|`, lapply(allow, function(a) startsWith(m_src, a)), rep(FALSE, length(m_src))) else rep(FALSE, length(m_src))
      idx <- if (any(in_allow)) which(in_allow) else seq_along(m_src)
      scope <- if (any(in_allow)) "ok" else "outside_metadata"
      sc_i <- idx[m_team[idx] == COLLECTION_TEAM]; if (length(sc_i)) idx <- sc_i
      if (var %in% MASS_VARS) {
        s81 <- idx[grepl("^Stephan_etal_1981_Table(I{1,3})$", m_src[idx])]
        if (length(s81)) idx <- s81[1]
      }
      pick <- idx[order(-m_yr[idx], m_src[idx])][1]
      pref <- m_src[pick]; pref_var <- m_var[pick]; pref_rule <- m_rule[pick]
    }
    dev <- dev_flag(sp_evo, pref_var)

    canonical_scale <- 1 / unname(col_unit_scale[[col]])
    unexpected_facs <- facs[abs(facs - canonical_scale) > .Machine$double.eps^0.5]
    scale_flag <- if (length(unexpected_facs))
      paste(sort(unique(unexpected_facs)), collapse = ";") else ""

    if (length(m_src) && scope != "outside_metadata") {
      status <- if (length(matched_src_u) == 1) "resolved_unique" else "resolved_multiple"
      closest <- if (identical(pref, "Stephan_etal_1981_TablesI-VI"))
                   "confirmed via Stephan 1981 Tables I-VI (absent from Evo-M1 merge)"
                 else if (nzchar(dev)) paste0("merge deviation flag: ", dev) else NA_character_
    } else if (length(m_src) && scope == "outside_metadata") {
      status <- "matched_outside_metadata"
      mp <- na.omit(vapply(strsplit(meta_refs %||% "", ";\\s*")[[1]], meta_phrase_to_key, character(1)))
      if (length(mp)) pref <- mp[1]
      closest <- sprintf("Evo-M1 value matches %s (not the metadata source); attributed to metadata source; see source_hunt",
                         paste(matched_src_u, collapse = "; "))
    } else if (is.null(cand)) {
      status <- "provenance_gap"; closest <- NA_character_
    } else {
      status <- "value_differs"
      j <- which.min(rel_pct(sval, cand$Value))
      closest <- sprintf("%s = %s (%.1f%%)", cand$Source[j], cand$Value[j], rel_pct(sval, cand$Value[j]))
    }

    i <- i + 1L
    rows[[i]] <- tibble(
      Species = sp_raw, Species_evo = sp_evo, Stephan_column = col,
      EvoM1_variable = pref_var, structure_key_rule = pref_rule,
      Stephan_value = sval, status = status,
      preferred_reference = pref, preferred_citation = cite_of(pref),
      all_matching_sources = paste(matched_src_u, collapse = "; "),
      n_candidate_sources = if (is.null(cand)) 0L else nrow(cand),
      metadata_reference = meta_refs, unit_rescale_flag = scale_flag,
      merge_deviation_flag = dev, closest_nonmatching = closest
    )
  }
}
long <- bind_rows(rows)

## --- 7. metadata-only columns (no Evo-M1 variable) --------------------------
meta_only_cols <- setdiff(
  names(S)[sapply(S, function(x) any(!is.na(x) & trimws(as.character(x)) != ""))],
  c(names(col2var), id_cols, "Species_evo")
)
mo <- list(); j <- 0L
for (r in seq_len(nrow(S))) for (col in meta_only_cols) {
  sstr <- as.character(S[[col]][r]); if (is.na(sstr) || trimws(sstr) == "") next
  sval <- suppressWarnings(as.numeric(sstr)); if (is.na(sval)) next
  j <- j + 1L
  mo[[j]] <- tibble(Species = S$Species[r], Species_evo = S$Species_evo[r],
                    Stephan_column = col, EvoM1_variable = NA_character_,
                    structure_key_rule = "no_defensible_volume_key",
                    Stephan_value = sval, status = "metadata_only",
                    preferred_reference = NA_character_, preferred_citation = NA_character_,
                    all_matching_sources = "", n_candidate_sources = 0L,
                    metadata_reference = meta_lookup(col), unit_rescale_flag = "",
                    merge_deviation_flag = "", closest_nonmatching = NA_character_)
}
if (length(mo)) long <- bind_rows(long, bind_rows(mo))

## --- 8. fallback source search ----------------------------------------------
## Direct near-perfect matches above remain untouched. Only weak cells enter
## this two-pass fallback. A 5% window is permitted only when source and
## structure (or a restricted metadata-listed paper) support the attribution;
## this produces a "likely" rather than a "resolved" status.
FALLBACK_MAX_PCT <- 5
weak_status <- c("value_differs", "matched_outside_metadata", "provenance_gap")

long <- long %>% mutate(
  resolution_method = case_when(
    status == "resolved_unique" ~ "direct_same_structure_unique",
    status == "resolved_multiple" ~ "direct_same_structure_multiple",
    status == "metadata_only" ~ "metadata_only_no_volume_crosswalk",
    status == "matched_outside_metadata" ~ "direct_match_outside_metadata",
    TRUE ~ "unresolved_after_direct_match"
  ),
  fallback_candidate_sources = "",
  fallback_match_detail = "",
  preferred_source_prior_contributions = 0L,
  source_cutoff_year = SOURCE_CUTOFF_YEAR
)

## First pass: learn which exact tables already contribute clean cells. This is
## the empirical prior requested for conflicts: a table known to have supplied
## many other cells outranks an otherwise equal, one-off numerical coincidence.
source_contrib <- long %>%
  filter(status %in% c("resolved_unique", "resolved_multiple"),
         !is.na(preferred_reference), preferred_reference != "") %>%
  count(preferred_reference, name = "n_prior_contributions")
source_contrib_lookup <- setNames(source_contrib$n_prior_contributions,
                                  source_contrib$preferred_reference)
contribution_n <- function(src) {
  z <- unname(source_contrib_lookup[src])
  ifelse(is.na(z), 0L, as.integer(z))
}
paper_contrib <- source_contrib %>%
  mutate(paper_key = ascii_key(source_paper(preferred_reference))) %>%
  group_by(paper_key) %>% summarise(n_prior_contributions = sum(n_prior_contributions),
                                    .groups = "drop")
paper_contrib_lookup <- setNames(paper_contrib$n_prior_contributions,
                                 paper_contrib$paper_key)
paper_contribution_n <- function(paper) {
  z <- unname(paper_contrib_lookup[ascii_key(paper)])
  ifelse(is.na(z), 0L, as.integer(z))
}
append_sources <- function(old, new) {
  old <- if (is.na(old) || !nzchar(old)) character(0) else strsplit(old, ";\\s*")[[1]]
  paste(unique(c(old, new[nzchar(new)])), collapse = "; ")
}

fallback_rows <- list(); fallback_i <- 0L

## 8a. Search all public final volume tables. Metadata-listed papers are tier 1
## even when the value lives in a different printed table/standardized column;
## the remaining Evo-M1 volume inventory is tier 2.
public_weak <- which(long$status %in% weak_status)
for (ii in public_weak) {
  sval <- long$Stephan_value[ii]
  sstr <- format(sval, scientific = FALSE, trim = TRUE)
  var <- long$EvoM1_variable[ii]
  if (is.na(var) || !nzchar(var)) next
  cand <- P_all[P_all$Species == long$Species_evo[ii], , drop = FALSE]
  if (!nrow(cand)) next

  refs <- long$metadata_reference[ii]
  refs <- if (is.na(refs)) "" else refs
  allow <- na.omit(unique(vapply(strsplit(refs, ";\\s*")[[1]],
                                 meta_phrase_to_key, character(1))))
  if (long$Stephan_column[ii] %in% names(col_primary))
    allow <- c(allow, col_primary[[long$Stephan_column[ii]]])
  scales <- match_scales_for_col(long$Stephan_column[ii])

  hits <- list(); nh <- 0L
  for (kk in seq_len(nrow(cand))) {
    pct <- vapply(scales, function(sc) rel_pct(sval, cand$Value[kk] * sc), numeric(1))
    near <- vapply(scales, function(sc) is_match(sval, sstr, cand$Value[kk], sc), logical(1))
    best_sc <- which.min(pct)
    metadata_match <- source_is_allowed(cand$Source[kk], allow)
    key_row <- candidate_key_row(long$Stephan_column[ii], cand$Variable[kk],
                                 cand$Source[kk])
    variable_match <- nrow(key_row) > 0
    ## Numeric similarity never licenses a different anatomical structure.
    eligible <- variable_match && (any(near) || pct[best_sc] <= FALLBACK_MAX_PCT)
    if (!eligible) next
    if (any(near)) best_sc <- which(near)[which.min(pct[near])]
    tier_rank <- if (metadata_match) 1L else 3L
    nh <- nh + 1L
    hits[[nh]] <- tibble(
      Species = long$Species[ii], Stephan_column = long$Stephan_column[ii],
      Stephan_value = sval, candidate_reference = cand$Source[kk],
      candidate_source = cand$Source[kk], candidate_variable = cand$Variable[kk],
      candidate_value = cand$Value[kk], candidate_year = cand$Year[kk],
      search_tier = ifelse(tier_rank <= 2, "metadata_tables", "evom1_all_volume_tables"),
      tier_rank = tier_rank,
      match_quality = ifelse(any(near), "near_perfect", "close_same_structure"),
      pct_difference = pct[best_sc], scale = scales[best_sc],
      prior_contributions = contribution_n(cand$Source[kk]),
      candidate_file = cand$Candidate_file[kk], candidate_sheet = "", candidate_cell = "",
      evidence_context = paste0("Species=", cand$Species[kk], "; Variable=", cand$Variable[kk],
                                "; key_rule=", key_row$key_rule[1]),
      structure_context_match = TRUE
    )
  }
  if (!length(hits)) next
  hit <- bind_rows(hits) %>%
    distinct(candidate_source, candidate_variable, candidate_value, scale, .keep_all = TRUE) %>%
    arrange(tier_rank, desc(prior_contributions),
            desc(match_quality == "near_perfect"), pct_difference,
            desc(candidate_year), candidate_source)
  for (jj in seq_len(nrow(hit))) {
    fallback_i <- fallback_i + 1L
    fallback_rows[[fallback_i]] <- hit[jj, ]
  }
  long$fallback_candidate_sources[ii] <- append_sources(
    long$fallback_candidate_sources[ii], unique(hit$candidate_reference))

  best <- hit[1, ]
  if (best$tier_rank <= 3) {
    long$preferred_reference[ii] <- best$candidate_reference
    long$preferred_citation[ii] <- cite_of(best$candidate_reference)
    long$EvoM1_variable[ii] <- best$candidate_variable
    chosen_key <- candidate_key_row(long$Stephan_column[ii], best$candidate_variable,
                                    best$candidate_source)
    long$structure_key_rule[ii] <- chosen_key$key_rule[1]
    long$all_matching_sources[ii] <- append_sources(long$all_matching_sources[ii],
                                                     hit$candidate_source[hit$tier_rank <= 3])
    is_near <- best$match_quality == "near_perfect"
    long$status[ii] <- if (best$search_tier == "metadata_tables") {
      if (is_near) "resolved_fallback_metadata" else "likely_fallback_metadata"
    } else {
      if (is_near) "resolved_fallback_evom1" else "likely_fallback_evom1"
    }
    long$resolution_method[ii] <- paste0(best$search_tier, "_",
                                          best$match_quality)
    long$fallback_match_detail[ii] <- sprintf(
      "%s = %s (%s; %.2f%%; scale %s; %d prior contributions)",
      best$candidate_source, format(best$candidate_value, trim = TRUE),
      best$candidate_variable, best$pct_difference, best$scale,
      best$prior_contributions)
    canonical_scale <- 1 / unname(col_unit_scale[[long$Stephan_column[ii]]])
    if (all(abs(best$scale - c(1, canonical_scale)) > .Machine$double.eps^0.5))
      long$unit_rescale_flag[ii] <- append_sources(long$unit_rescale_flag[ii],
                                                   as.character(best$scale))
  }
}

## 8b. Build a dated restricted-document inventory. MIGRATED_INDEX is the
## maintained list of private comparison inputs; limiting it to paper families
## already present in the public volume ledger prevents unrelated energetic,
## behavioural, and cell-count sheets from creating numerical coincidences.
restricted_inventory <- tibble(
  paper_folder = character(), source_year = integer(), relative_file = character(),
  absolute_file = character(), cutoff_basis = character()
)
if (SEARCH_RESTRICTED && dir.exists(EVOM1_RESTRICTED_ROOT)) {
  migrated_file <- file.path(EVOM1_RESTRICTED_ROOT, "MIGRATED_INDEX.csv")
  if (file.exists(migrated_file)) {
    migrated <- read.csv(migrated_file, stringsAsFactors = FALSE, check.names = FALSE)
    migrated_year <- suppressWarnings(as.integer(str_extract(migrated$paper_folder,
                                                               "(19|20)[0-9]{2}")))
    volume_paper_keys <- unique(ascii_key(source_paper(P_all$Source)))
    restricted_inventory <- migrated %>%
      mutate(source_year = migrated_year,
             paper_key = ascii_key(paper_folder),
             absolute_file = file.path(EVOM1_RESTRICTED_ROOT, moved_to)) %>%
      filter(startsWith(role, "input"),
             grepl("\\.(csv|tsv|xls|xlsx)$", moved_to, ignore.case = TRUE),
             paper_key %in% volume_paper_keys,
             !is.na(source_year), source_year <= SOURCE_CUTOFF_YEAR,
             file.exists(absolute_file)) %>%
      transmute(paper_folder, source_year, relative_file = moved_to,
                absolute_file, cutoff_basis = "publication year parsed from paper folder")
  }

  ## These volume-bearing private collections are not migration/check inputs,
  ## so list them explicitly with a defensible pre-compilation date basis.
  manual_specs <- tribble(
    ~relative_dir, ~paper_folder, ~source_year, ~cutoff_basis,
    "unpublished_data/____Unpublished__Frahm_Stephan_individuals",
      "Unpublished_Frahm_Stephan_individuals", 2019L,
      "manually classified as an input to the 2019 Stephan compilation",
    "unpublished_data/____Unpublished__MacLeod_neocortex",
      "MacLeod_etal_2003", 2003L, "associated publication year",
    "_deSousaDissertation_staging",
      "deSousa_dissertation", 2010L, "dissertation/publication predates compilation"
  )
  manual_inventory <- pmap_dfr(manual_specs, function(relative_dir, paper_folder,
                                                       source_year, cutoff_basis) {
    d <- file.path(EVOM1_RESTRICTED_ROOT, relative_dir)
    ff <- if (dir.exists(d)) list.files(d, pattern = "\\.(csv|tsv|xls|xlsx)$",
                                       full.names = TRUE, recursive = TRUE,
                                       ignore.case = TRUE) else character(0)
    if (!length(ff)) return(tibble())
    tibble(paper_folder = paper_folder, source_year = source_year,
           relative_file = sub(paste0("^", EVOM1_RESTRICTED_ROOT, "/?"), "", ff),
           absolute_file = ff, cutoff_basis = cutoff_basis)
  })
  restricted_inventory <- bind_rows(restricted_inventory, manual_inventory) %>%
    distinct(absolute_file, .keep_all = TRUE) %>%
    filter(source_year <= SOURCE_CUTOFF_YEAR)
} else if (SEARCH_RESTRICTED) {
  warning("Restricted fallback skipped: repo not mounted at ", EVOM1_RESTRICTED_ROOT,
          ". Set EVOM1_RESTRICTED to override.", call. = FALSE)
}

## Cell-level index preserves both row context and the first five cells above a
## numeric cell. The latter makes transposed sheets searchable without treating
## every species appearing anywhere in a column as evidence for every value.
index_documents <- function(inv) {
  out <- list(); oi <- 0L
  if (!nrow(inv)) return(tibble())
  for (dd in seq_len(nrow(inv))) {
    frames <- read_document(inv$absolute_file[dd])
    if (is.null(frames)) {
      message("  (skipped unreadable restricted file: ", inv$relative_file[dd], ")")
      next
    }
    for (shn in names(frames)) {
      df <- frames[[shn]]
      if (is.null(df) || !nrow(df) || !ncol(df)) next
      cells <- lapply(df, function(z) iconv(as.character(z), to = "UTF-8", sub = ""))
      for (cc in seq_along(cells)) {
        vals <- parse_doc_number(cells[[cc]])
        rr_hit <- which(!is.na(vals))
        if (!length(rr_hit)) next
        header_rows <- seq_len(min(5L, nrow(df)))
        header_text <- paste(na.omit(cells[[cc]][header_rows]), collapse = " | ")
        row_text <- vapply(rr_hit, function(rr) {
          z <- vapply(cells, `[`, character(1), rr)
          z <- z[!is.na(z) & nzchar(trimws(z))]
          substr(paste(z, collapse = " | "), 1, 500)
        }, character(1))
        oi <- oi + 1L
        out[[oi]] <- tibble(
          value = vals[rr_hit], paper_folder = inv$paper_folder[dd],
          source_year = inv$source_year[dd], relative_file = inv$relative_file[dd],
          sheet = ifelse(shn == "_", "", shn), row = rr_hit, column = cc,
          row_text = row_text, column_header = substr(header_text, 1, 250)
        )
      }
    }
  }
  bind_rows(out)
}

restricted_index <- tibble()
if (nrow(restricted_inventory)) {
  message("restricted source fallback: indexing ", nrow(restricted_inventory),
          " pre-2020 volume document(s) ...")
  restricted_index <- index_documents(restricted_inventory)
  message("restricted source fallback: ", nrow(restricted_index), " numeric cells indexed")
}

species_match_strength <- function(context, species_raw, species_evo) {
  context <- tolower(gsub("_", " ", context))
  raw_space <- gsub("_", " ", species_raw)
  names_to_try <- unique(c(raw_space, species_evo,
    if (species_raw %in% names(s81_sp)) s81_sp[[species_raw]] else character(0)))
  names_to_try <- tolower(names_to_try[!is.na(names_to_try) & nzchar(names_to_try)])
  if (any(vapply(names_to_try, grepl, logical(1), x = context, fixed = TRUE))) return(2L)
  genera <- unique(sub(" .*", "", names_to_try))
  if (any(vapply(genera, grepl, logical(1), x = context, fixed = TRUE))) return(1L)
  0L
}

structure_context_match <- function(col, context) {
  z <- tolower(gsub("[_\\.]", " ", iconv(context, to = "ASCII//TRANSLIT", sub = "")))
  pat <- switch(col,
    "Body_weight" = "body *(weight|mass)", "Body_mass" = "body *(weight|mass)",
    "Body_weight_1985a" = "body *(weight|mass)",
    "Brain_weight" = "brain *(weight|mass)", "Brain_mass" = "brain *(weight|mass)",
    "ASG_Sousa" = "area *striata|(^|[^a-z])v1([^a-z]|$)",
    "LGN_Sousa" = "lgn|geniculatum",
    "NeoG_Frahm" = "neocortex.*grey|grey.*neocortex",
    "NeoW_Frahm" = "neocortex.*white|white.*neocortex",
    "Granular_volume_L" = "(^|[^a-z])granular([^a-z]|$)",
    "Dysgranular_volume_L" = "dysgranular",
    "Agranular_volume_L" = "agranular",
    "FI_volume_L" = "fronto *insular|(^|[^a-z])fi([^a-z]|$)",
    "Total_insula_volume_L" = "insula|(^|[^a-z])total([^a-z]|$)",
    "Pons" = "ventral *pons|(^|[^a-z])vpo([^a-z]|$)|(^|[^a-z])pons([^a-z]|$)",
    "Primary.visual.White.Smaers" = "primary *visual.*white|white.*primary *visual|area *striata.*white",
    "Prefrontal.Gray" = "prefrontal.*(gray|grey)|(gray|grey).*prefrontal",
    "Prefrontal.White" = "prefrontal.*white|white.*prefrontal",
    "Frontal.motor.Gray" = "frontal *motor.*(gray|grey)|(gray|grey).*frontal *motor",
    "Frontal.motor" = "frontal *motor.*white|white.*frontal *motor",
    "Other.cortical.association.areas.Gray" = "association.*(gray|grey)|(gray|grey).*association",
    "Other.cortical.association.areas.White" = "association.*white|white.*association",
    NA_character_
  )
  if (!is.na(pat)) return(grepl(pat, z, perl = TRUE))
  tokens <- strsplit(tolower(gsub("[_\\.]", " ", col)), " +")[[1]]
  tokens <- tokens[!tokens %in% c("volume", "vol", "left", "right", "mm3", "sousa", "frahm") &
                     nchar(tokens) >= 4]
  if (!length(tokens)) return(FALSE)
  grepl(paste(tokens, collapse = ".*"), z, perl = TRUE)
}

## Search restricted documents only for cells still weak after the public pass,
## including public fallbacks: an exact hit in a metadata-listed restricted
## document may supersede a public match outside the metadata source.
restricted_weak <- public_weak
if (nrow(restricted_index) && length(restricted_weak)) for (ii in restricted_weak) {
  sval <- long$Stephan_value[ii]
  sstr <- format(sval, scientific = FALSE, trim = TRUE)
  scales <- c(1, 1000, 0.001)
  dh <- list(); dn <- 0L
  for (sc in scales) {
    pct <- rel_pct(sval, restricted_index$value * sc)
    near <- is_match_vector(sval, sstr, restricted_index$value, sc)
    hh <- which(near | pct <= FALLBACK_MAX_PCT)
    if (!length(hh)) next
    for (kk in hh) {
      context <- paste(restricted_index$row_text[kk],
                       restricted_index$column_header[kk], sep = " | ")
      sp_strength <- species_match_strength(context, long$Species[ii], long$Species_evo[ii])
      if (!sp_strength) next
      refs <- long$metadata_reference[ii]
      refs <- if (is.na(refs)) "" else refs
      allow <- na.omit(unique(vapply(strsplit(refs, ";\\s*")[[1]],
                                     meta_phrase_to_key, character(1))))
      if (long$Stephan_column[ii] %in% names(col_primary))
        allow <- c(allow, col_primary[[long$Stephan_column[ii]]])
      meta_match <- paper_is_allowed(restricted_index$paper_folder[kk], allow)
      prior <- paper_contribution_n(restricted_index$paper_folder[kk])
      structure_match <- structure_context_match(
        long$Stephan_column[ii], paste(restricted_index$row_text[kk],
                                       restricted_index$column_header[kk], sep = " | "))
      ## Structure evidence is mandatory even for an exact number in a listed
      ## paper. This blocks same-row coincidences such as an insula value that
      ## happens to equal Septum; values corroborate a key but never create it.
      if (!structure_match) next
      if (!near[kk] && !meta_match) next
      if (!meta_match && (prior == 0 || !structure_match)) next
      ref <- paste0(restricted_index$paper_folder[kk], " [restricted: ",
                    basename(restricted_index$relative_file[kk]), "]")
      dn <- dn + 1L
      dh[[dn]] <- tibble(
        Species = long$Species[ii], Stephan_column = long$Stephan_column[ii],
        Stephan_value = sval, candidate_reference = ref,
        candidate_source = restricted_index$paper_folder[kk],
        candidate_variable = NA_character_, candidate_value = restricted_index$value[kk],
        candidate_year = restricted_index$source_year[kk], search_tier = "restricted_documents",
        tier_rank = ifelse(meta_match, 5L, 6L),
        match_quality = ifelse(near[kk], "near_perfect", "close_metadata_source"),
        pct_difference = pct[kk], scale = sc, prior_contributions = prior,
        candidate_file = restricted_index$relative_file[kk],
        candidate_sheet = restricted_index$sheet[kk],
        candidate_cell = paste0("R", restricted_index$row[kk], "C", restricted_index$column[kk]),
        ## Do not copy a restricted row into this repo. The private file/sheet/
        ## cell locator is sufficient for a curator to inspect it in place.
        evidence_context = paste0("species_match_strength=", sp_strength,
                                  "; structure_context_match=", structure_match),
        species_match_strength = sp_strength,
        metadata_source_match = meta_match, structure_context_match = structure_match
      )
    }
  }
  if (!length(dh)) next
  hit <- bind_rows(dh) %>%
    distinct(candidate_reference, candidate_file, candidate_sheet, candidate_cell,
             scale, .keep_all = TRUE) %>%
    arrange(tier_rank, desc(structure_context_match), desc(prior_contributions),
            desc(species_match_strength),
            desc(match_quality == "near_perfect"), pct_difference,
            candidate_file, candidate_sheet, candidate_cell)
  for (jj in seq_len(nrow(hit))) {
    fallback_i <- fallback_i + 1L
    fallback_rows[[fallback_i]] <- hit[jj, setdiff(names(hit),
                                                    c("species_match_strength",
                                                      "metadata_source_match"))]
  }
  long$fallback_candidate_sources[ii] <- append_sources(
    long$fallback_candidate_sources[ii], unique(hit$candidate_reference))
  best <- hit[1, ]
  currently_likely <- startsWith(long$status[ii], "likely_")
  outside_metadata_public <- startsWith(long$status[ii], "resolved_fallback_evom1")
  can_override <- long$status[ii] %in% weak_status ||
    (currently_likely && best$match_quality == "near_perfect") ||
    (outside_metadata_public && best$metadata_source_match &&
       best$match_quality == "near_perfect")
  if (can_override) {
    long$preferred_reference[ii] <- best$candidate_reference
    long$preferred_citation[ii] <- cite_of(best$candidate_reference)
    long$all_matching_sources[ii] <- append_sources(long$all_matching_sources[ii],
                                                     best$candidate_reference)
    long$status[ii] <- if (best$match_quality == "near_perfect")
      "resolved_fallback_restricted" else "likely_fallback_restricted"
    long$resolution_method[ii] <- paste0("restricted_documents_", best$match_quality)
    long$fallback_match_detail[ii] <- sprintf(
      "%s = %s (%.2f%%; scale %s; %s %s; %d prior contributions)",
      best$candidate_reference, format(best$candidate_value, trim = TRUE),
      best$pct_difference, best$scale, best$candidate_file, best$candidate_cell,
      best$prior_contributions)
    canonical_scale <- 1 / unname(col_unit_scale[[long$Stephan_column[ii]]])
    if (all(abs(best$scale - c(1, canonical_scale)) > .Machine$double.eps^0.5))
      long$unit_rescale_flag[ii] <- append_sources(long$unit_rescale_flag[ii],
                                                   as.character(best$scale))
  }
}

fallback_candidates <- if (length(fallback_rows)) bind_rows(fallback_rows) else tibble()

## Stamp the learned contribution prior for every final preferred source,
## including restricted labels (looked up by their paper-family prefix).
long$preferred_source_prior_contributions <- vapply(long$preferred_reference, function(src) {
  if (is.na(src) || !nzchar(src)) return(0L)
  exact <- contribution_n(src)
  if (exact) exact else paper_contribution_n(source_paper(sub(" \\[restricted:.*$", "", src)))
}, integer(1))

source_inventory <- bind_rows(
  P_all %>% distinct(Source, Team, Year, Candidate_file) %>%
    transmute(search_tier = "evom1_all_volume_tables", source = Source,
              source_year = Year, team = Team,
              file = Candidate_file,
              cutoff_basis = ifelse(Candidate_file == f_prov_all,
                "publication year in canonical volume ledger",
                "publication year of finalized public table named by metadata")),
  restricted_inventory %>%
    transmute(search_tier = "restricted_documents", source = paper_folder,
              source_year, team = NA_character_, file = relative_file, cutoff_basis)
) %>% arrange(search_tier, source_year, source, file)

long <- long %>% arrange(Species, Stephan_column)

## --- 9. write outputs -------------------------------------------------------
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write_csv(long, file.path(out_dir, "Stephan_primates_references_long.csv"))
write_csv(source_inventory,
          file.path(out_dir, "Stephan_primates_candidate_source_inventory.csv"))
if (nrow(fallback_candidates))
  write_csv(fallback_candidates,
            file.path(out_dir, "Stephan_primates_fallback_candidates.csv"))

mism <- long %>%
  filter(status %in% c("value_differs", "matched_outside_metadata") |
           startsWith(status, "likely_") | grepl("^resolved_fallback", status) |
           unit_rescale_flag != "" | merge_deviation_flag != "") %>%
  select(Species, Stephan_column, EvoM1_variable, structure_key_rule,
         Stephan_value, status,
         resolution_method, preferred_reference, preferred_source_prior_contributions,
         unit_rescale_flag, merge_deviation_flag, closest_nonmatching,
         fallback_match_detail, fallback_candidate_sources,
         all_matching_sources, metadata_reference) %>%
  arrange(desc(status == "value_differs"), desc(status == "matched_outside_metadata"),
          desc(merge_deviation_flag != ""), Species, Stephan_column)
write_csv(mism, file.path(out_dir, "Stephan_primates_reference_mismatches.csv"))

by_col <- long %>%
  group_by(Stephan_column, EvoM1_variable) %>%
  summarise(n_cells = n(),
            n_resolved = sum(startsWith(status, "resolved")),
            n_unique   = sum(status == "resolved_unique"),
            n_multiple = sum(status == "resolved_multiple"),
            n_resolved_fallback = sum(startsWith(status, "resolved_fallback")),
            n_likely_fallback = sum(startsWith(status, "likely_fallback")),
            n_value_differs = sum(status == "value_differs"),
            n_gap = sum(status == "provenance_gap"),
            n_metadata_only = sum(status == "metadata_only"),
            metadata_reference = dplyr::first(metadata_reference),
            .groups = "drop") %>%
  arrange(desc(n_value_differs), Stephan_column)
write_csv(by_col, file.path(out_dir, "Stephan_primates_references_by_column.csv"))

## --- 10. source hunt: trace unresolved cells to the raw cut-and-paste dir ----
## For cells that did not resolve cleanly against Evo-M1 (value_differs,
## matched_outside_metadata, provenance_gap), search the original messy source
## material for the value and report the file/row where it appears. This is
## where the compiler's number actually came from (e.g. Bauernfeind converted
## files, Visual_volumes.xlsx). Set HUNT <- FALSE to skip (it is slow: it reads
## every csv/xls/xlsx under the directory).
HUNT <- TRUE
hunt_dir <- "../datasets brain and regions species/Stephan and Frahm"
if (HUNT && dir.exists(hunt_dir)) {
  message("source hunt: indexing ", hunt_dir, " ...")
  hfiles <- list.files(hunt_dir, pattern = "\\.(csv|xls|xlsx)$", full.names = TRUE, recursive = TRUE)
  hfiles <- hfiles[!grepl("~\\$", hfiles)]
  hidx <- list(); hi <- 0L
  read_any <- function(f) {
    tryCatch({
      if (grepl("\\.csv$", f, ignore.case = TRUE)) {
        # read as latin1: accepts any byte, never errors on old encodings.
        # lazy = FALSE so any parse error surfaces here (inside tryCatch), not later.
        list(`_` = suppressWarnings(suppressMessages(
          readr::read_csv(f, col_names = FALSE, col_types = cols(.default = "c"),
                          locale = readr::locale(encoding = "latin1"),
                          lazy = FALSE, show_col_types = FALSE, progress = FALSE))))
      } else {
        sh <- readxl::excel_sheets(f)
        setNames(lapply(sh, function(s)
          suppressMessages(readxl::read_excel(f, sheet = s, col_names = FALSE, col_types = "text"))), sh)
      }
    }, error = function(e) NULL)
  }
  index_file <- function(f) {
    rel <- sub(paste0("^", hunt_dir, "/?"), "", f)
    frames <- read_any(f); if (is.null(frames)) return(invisible())
    for (shn in names(frames)) {
      df <- frames[[shn]]; if (is.null(df) || !nrow(df)) next
      cells_all <- lapply(df, function(cc) iconv(as.character(cc), to = "UTF-8", sub = ""))
      for (rr in seq_len(nrow(df))) {
        cells <- vapply(cells_all, `[`, character(1), rr)
        vals <- suppressWarnings(as.numeric(cells))
        vals <- vals[!is.na(vals)]; if (!length(vals)) next
        rowtext <- paste(cells, collapse = " | "); if (is.na(rowtext)) rowtext <- ""
        hi <<- hi + 1L
        hidx[[hi]] <<- list(vals = vals, rel = rel, sheet = if (shn == "_") "" else shn,
                            rowtext = substr(rowtext, 1, 150))
      }
    }
  }
  for (f in hfiles) tryCatch(index_file(f), error = function(e)
    message("  (skipped unreadable file: ", basename(f), ")"))
  hunt_one <- function(val, token) {
    hits <- Filter(function(h) any(abs(h$vals - val) <= 0.005 * pmax(abs(h$vals), abs(val))), hidx)
    if (!length(hits)) return(NULL)
    tok <- tolower(token)
    spmatch <- vapply(hits, function(h) grepl(tok, tolower(h$rowtext), fixed = TRUE), logical(1))
    hits <- hits[order(!spmatch)]; spmatch <- sort(spmatch, decreasing = TRUE)
    lapply(seq_len(min(3, length(hits))), function(k)
      tibble(hunt_species_match = ifelse(spmatch[k], "yes", "no"),
             hunt_file = hits[[k]]$rel, hunt_sheet = hits[[k]]$sheet, hunt_rowtext = hits[[k]]$rowtext))
  }
  review <- long %>% filter(status %in% c("value_differs", "matched_outside_metadata", "provenance_gap"))
  hunt_rows <- list(); hk <- 0L
  for (r in seq_len(nrow(review))) {
    token <- strsplit(review$Species[r], "_")[[1]][1]
    hh <- hunt_one(review$Stephan_value[r], token)
    base <- tibble(Species = review$Species[r], Stephan_column = review$Stephan_column[r],
                   Stephan_value = review$Stephan_value[r], status = review$status[r])
    if (is.null(hh)) {
      hk <- hk + 1L
      hunt_rows[[hk]] <- bind_cols(base, tibble(hunt_rank = 1L, hunt_species_match = "",
        hunt_file = "(value not found in Stephan & Frahm dir)", hunt_sheet = "", hunt_rowtext = ""))
    } else for (k in seq_along(hh)) {
      hk <- hk + 1L; hunt_rows[[hk]] <- bind_cols(base, mutate(hh[[k]], hunt_rank = k))
    }
  }
  if (length(hunt_rows))
    write_csv(bind_rows(hunt_rows) %>%
                select(Species, Stephan_column, Stephan_value, status, hunt_rank,
                       hunt_species_match, hunt_file, hunt_sheet, hunt_rowtext),
              file.path(out_dir, "Stephan_primates_source_hunt.csv"))
  message("source hunt: ", length(hunt_rows), " rows written")
}

## --- 11. console summary ----------------------------------------------------
message("cells written: ", nrow(long))
print(long %>% count(status) %>% arrange(desc(n)))
message("cells flagged for review (mismatches file): ", nrow(mism))
