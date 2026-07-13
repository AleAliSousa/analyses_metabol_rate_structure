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
#   (volumes_unfiltered.csv: Species x Variable x Value x Source-table x Year).
#
# INPUTS  (all read relative to the project root set below)
#   data_raw/Stephan_primates.csv                     - the compiled data
#   metadata/Stephan_primates_metadata.xlsx           - per-variable references
#   data_intermediate/volumes_unfiltered.csv          - Evo-M1 per-source values
#       (mirror of  Evo-M1-Trait-Data/__merging_volumes/volumes_unfiltered.csv;
#        refresh this file when the Evo-M1 merge is re-run.)
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
#   * When several sources carry the same value, the "preferred" reference is
#     the one consistent with the metadata's stated source(s) for that column,
#     breaking ties by earliest publication (the original measurement).
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
## Edit if you run from elsewhere. All paths below are relative to this.
proj <- "~/Library/CloudStorage/Dropbox/COLLABORATIVE/Do expensive brain regions increase less in humans/analyses_metabol_rate_structure"
if (dir.exists(proj)) setwd(proj)

f_data <- "data_raw/Stephan_primates.csv"
f_meta <- "metadata/Stephan_primates_metadata.xlsx"
f_prov <- "data_intermediate/volumes_unfiltered.csv"
out_dir <- "metadata"

stopifnot(file.exists(f_data), file.exists(f_meta), file.exists(f_prov))

## --- 1. species crosswalk (Stephan space-form -> Evo-M1 name) ---------------
sp_xwalk <- c(
  "Lagothrix lagotricha" = "Lagothrix lagothricha",  # spelling
  "Pongo pygmaeus"       = "Pongo sp.",              # Evo-M1 pools orangutans
  "Pongo abelii"         = "Pongo sp."               #   as "Pongo sp."
)
norm_sp <- function(x) {
  x <- gsub("_", " ", trimws(x))
  hit <- x %in% names(sp_xwalk)
  x[hit] <- unname(sp_xwalk[x[hit]])
  x
}

## --- 2. column -> Evo-M1 variable crosswalk ---------------------------------
## Columns not listed here have no clean Evo-M1 volumetric equivalent and are
## carried through as "metadata_only" (reference taken from the xlsx).
col2var <- c(
  Body_weight = "Body_Mass.g", Brain_weight = "Brain_Mass.mg",
  Ventricles = "Ventricles_Vol.mm3",
  Total_brain_net_volume = "Total_brain_net_volume_Vol.mm3",
  Medulla_oblongata = "Medulla_oblongata_Vol.mm3", Cerebellum = "Cerebellum_Vol.mm3",
  Mesencephalon = "Mesencephalon_Vol.mm3", Diencephalon = "Diencephalon_Vol.mm3",
  Telencephalon = "Telencephalon_Vol.mm3", Bulbus_olfactorius = "Bulbus_olfactorius_Vol.mm3",
  Bulbus_olfactorius_accessorius = "Bulbus_olfactorius_accessorius_Vol.mm3",
  Lobus_piriformis = "Lobus_piriformis_Vol.mm3", Septum = "Septum_Vol.mm3",
  Striatum = "Striatum_Vol.mm3", Schizocortex = "Schizo_cortex_Vol.mm3",
  Hippocampus = "Hippocampus_Vol.mm3", NeoWG = "Neocortex_Vol.mm3",
  Epithalamus = "Epithalamus_Vol.mm3", Thalamus = "Thalamus_Vol.mm3",
  Hypothalamus = "Hypothalamus_Vol.mm3", Subthalamus = "Subthalamus_Vol.mm3",
  Pallidum = "Pallidum_Vol.mm3", Nucleus_subthalamicus = "Nucleus_subthalamicus_Vol.mm3",
  Capsula_interna = "Capsula_interna_Vol.mm3", Tractus_opticus = "Tractus_opticus_Vol.mm3",
  Palaeocortex = "Palaeocortex_Vol.mm3", Amygdala = "Amygdala_Vol.mm3",
  Complexus_centromedialis = "Complexus_centromedialis_Vol.mm3",
  Nucleus_tractus_olfactorius = "Nucleus_tractus_olfactorius_Vol.mm3",
  Complexus_cortico_basolateralis = "Complexus_corticobasolateralis_Vol.mm3",
  Nucleus_amygdalae_basalis_parsmagnocellularis = "Nucleus_amygdalae_basalis_pars_magnocellularis_Vol.mm3",
  NeoG_Frahm = "Neocortex_grey_matter_Vol.mm3", NeoW_Frahm = "Neocortex_white_matter_Vol.mm3",
  ASG_Sousa = "Area_striata_grey_matter_Vol.mm3", LGN_Sousa = "Corpus_geniculatum_laterale_Vol.mm3",
  Pons = "Pons_Vol.mm3", Body_mass = "Body_Mass.g", Brain_mass = "Brain_Mass.mg",
  Granular_volume_L = "Granular_insular_cortex_left_Vol.mm3",
  Dysgranular_volume_L = "Dysgranular_insular_cortex_left_Vol.mm3",
  Agranular_volume_L = "Agranular_insular_cortex_left_Vol.mm3",
  FI_volume_L = "fronto_insular_cortex_left_Vol.mm3",
  Total_insula_volume_L = "Insula_left_Vol.mm3",
  Complexus.vestibularis = "Complexus_vestibularis_unilateral_Vol.mm3",
  Nucleus.vestibularis.superior = "Nucleus_vestibularis_superior_unilateral_Vol.mm3",
  Nucleus.vestibularis.lateralis = "Nucleus_vestibularis_lateralis_unilateral_Vol.mm3",
  Nucleus.vestibularis.medialis = "Nucleus_vestibularis_medialis_unilateral_Vol.mm3",
  Nucleus.vestibularis.descendens = "Nucleus_vestibularis_descendens_unilateral_Vol.mm3",
  Nucleus.septalis.triangularis = "Nucleus_septalis_triangularis_Vol.mm3",
  Corpus.subfornicale = "Corpus_subfornicale_Vol.mm3",
  Nucleus.habenularis.medialis = "Nucleus_habenularis_medialis_Vol.mm3",
  Corpus.pineale = "Corpus_pineal_Vol.mm3",
  Corpus.subcommissurale = "Corpus_subcommissurale_Vol.mm3",
  Cerebellar_nuclei_total = "Cerebellar_nuclei_total_Vol.mm3",
  Interpositus_cerebellar_nuclei = "Interpositus_cerebellar_nuclei_Vol.mm3",
  Lateral_cerebellar_nuclei = "Lateral_cerebellar_nuclei_Vol.mm3",
  Medial_cerebellar_nuclei = "Medial_cerebellar_nuclei_Vol.mm3",
  Body_weight_1985a = "Body_Mass.g"
)
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
  "Sherwood_etal_2005", "Sherwood et al. 2005 (J Hum Evol 48:45-84)",
  "Zilles_Rehkämper_1988", "Zilles & Rehkamper 1988 (Table 12-2)"
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
P$Value <- suppressWarnings(as.numeric(P$Value))
P <- P[!is.na(P$Value), ]
prov <- split(P[, c("Source", "Value", "Year")], list(P$Species, P$Variable), drop = TRUE)
prov_key <- function(sp, var) paste(sp, var, sep = ".")

## --- fallback provenance: Stephan et al 1981 combined Tables I-VI -----------
## Consulted only when the primary merge (volumes_unfiltered) yields no match.
## The merge drops some 1981 "revised" rows (e.g. Gorilla, Homo), keeping the
## older 1970 values, so those cells otherwise look like differences. Matching
## is value-gated, so the species synonym map below need not be exhaustive:
## wrong guesses simply fail to match and cost nothing.
f_s81 <- "data_intermediate/Stephan_etal_1981_TablesI-VI.csv"
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
S81 <- if (file.exists(f_s81)) read.csv(f_s81, stringsAsFactors = FALSE, check.names = FALSE) else NULL
s81_value <- function(sp_underscore, col) {
  if (is.null(S81) || !col %in% names(s81_col)) return(NA_real_)
  printed <- if (sp_underscore %in% names(s81_sp)) s81_sp[[sp_underscore]] else gsub("_", " ", sp_underscore)
  r <- which(S81$Species == printed); if (!length(r)) return(NA_real_)
  suppressWarnings(as.numeric(S81[[ s81_col[[col]] ]][r[1]]))
}

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

## --- 6. walk every data cell ------------------------------------------------
rows <- list(); i <- 0L
for (r in seq_len(nrow(S))) {
  sp_raw <- S$Species[r]; sp_evo <- S$Species_evo[r]
  for (col in names(col2var)) {
    sstr <- as.character(S[[col]][r]); if (is.na(sstr) || trimws(sstr) == "") next
    sval <- suppressWarnings(as.numeric(sstr)); if (is.na(sval)) next
    var  <- col2var[[col]]
    meta_refs <- meta_lookup(col)
    cand <- prov[[prov_key(sp_evo, var)]]
    scales <- if (col %in% mass_cols) c(1, 1000, 0.001) else 1
    matched_src <- character(0); facs <- numeric(0)
    if (!is.null(cand)) {
      for (k in seq_len(nrow(cand))) {
        for (sc in scales) {
          if (is_match(sval, sstr, cand$Value[k], sc)) {
            matched_src <- c(matched_src, cand$Source[k]); facs <- c(facs, sc); break
          }
        }
      }
    }
    matched_src_u <- sort(unique(matched_src))

    # preferred reference: prefer sources consistent with the metadata phrase(s)
    pref <- NA_character_
    if (length(matched_src_u)) {
      allow <- na.omit(unique(vapply(strsplit(meta_refs %||% "", ";\\s*")[[1]],
                                     meta_phrase_to_key, character(1))))
      pick <- matched_src_u[Reduce(`|`, lapply(allow, function(a) startsWith(matched_src_u, a)),
                                    rep(FALSE, length(matched_src_u)))]
      if (!length(pick)) pick <- matched_src_u
      yrs <- suppressWarnings(as.integer(sub(".*?(\\d{4}).*", "\\1", pick)))
      pref <- pick[order(yrs, pick)][1]
    }

    # fallback to Stephan 1981 combined table if the merge produced no match
    fb_used <- FALSE
    if (!length(matched_src_u)) {
      fv <- s81_value(sp_raw, col)
      if (!is.na(fv) && is_match(sval, sstr, fv, 1)) {
        matched_src_u <- "Stephan_etal_1981_TablesI-VI"
        pref <- "Stephan_etal_1981_TablesI-VI"; fb_used <- TRUE
      }
    }

    # status
    if (length(matched_src_u)) {
      status <- if (length(matched_src_u) == 1) "resolved_unique" else "resolved_multiple"
      closest <- if (fb_used)
        "matched via Stephan 1981 Tables I-VI (revised values; absent from Evo-M1 merge)" else NA_character_
    } else if (is.null(cand)) {
      status <- "provenance_gap"; closest <- NA_character_
    } else {
      status <- "value_differs"
      j <- which.min(rel_pct(sval, cand$Value))
      closest <- sprintf("%s = %s (%.1f%%)", cand$Source[j], cand$Value[j],
                         rel_pct(sval, cand$Value[j]))
    }
    scale_flag <- if (length(facs) && any(facs != 1)) paste(sort(unique(facs[facs != 1])), collapse = ";") else ""

    i <- i + 1L
    rows[[i]] <- tibble(
      Species = sp_raw, Species_evo = sp_evo, Stephan_column = col,
      EvoM1_variable = var, Stephan_value = sval,
      status = status,
      preferred_reference = pref,
      preferred_citation = cite_of(pref),
      all_matching_sources = paste(matched_src_u, collapse = "; "),
      n_candidate_sources = if (is.null(cand)) 0L else nrow(cand),
      metadata_reference = meta_refs,
      unit_rescale_flag = scale_flag,
      closest_nonmatching = closest
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
                    Stephan_value = sval, status = "metadata_only",
                    preferred_reference = NA_character_, preferred_citation = NA_character_,
                    all_matching_sources = "", n_candidate_sources = 0L,
                    metadata_reference = meta_lookup(col), unit_rescale_flag = "",
                    closest_nonmatching = NA_character_)
}
if (length(mo)) long <- bind_rows(long, bind_rows(mo))

long <- long %>% arrange(Species, Stephan_column)

## --- 8. write outputs -------------------------------------------------------
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write_csv(long, file.path(out_dir, "Stephan_primates_references_long.csv"))

mism <- long %>%
  filter(status == "value_differs" | unit_rescale_flag != "") %>%
  select(Species, Stephan_column, EvoM1_variable, Stephan_value, status,
         unit_rescale_flag, closest_nonmatching, all_matching_sources,
         metadata_reference) %>%
  arrange(desc(status == "value_differs"), Species, Stephan_column)
write_csv(mism, file.path(out_dir, "Stephan_primates_reference_mismatches.csv"))

by_col <- long %>%
  group_by(Stephan_column, EvoM1_variable) %>%
  summarise(n_cells = n(),
            n_resolved = sum(status %in% c("resolved_unique", "resolved_multiple")),
            n_unique   = sum(status == "resolved_unique"),
            n_multiple = sum(status == "resolved_multiple"),
            n_value_differs = sum(status == "value_differs"),
            n_gap = sum(status == "provenance_gap"),
            n_metadata_only = sum(status == "metadata_only"),
            metadata_reference = dplyr::first(metadata_reference),
            .groups = "drop") %>%
  arrange(desc(n_value_differs), Stephan_column)
write_csv(by_col, file.path(out_dir, "Stephan_primates_references_by_column.csv"))

## --- 9. console summary -----------------------------------------------------
message("cells written: ", nrow(long))
print(long %>% count(status) %>% arrange(desc(n)))
message("mismatches / rescales flagged for review: ", nrow(mism))
