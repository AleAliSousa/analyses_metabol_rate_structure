#!/usr/bin/env Rscript
# ============================================================
# Study 3 - Prepare the 100-tree phylogenetic sample
#
# Purpose:
#   1. Read a species list produced upstream in the pipeline.
#   2. Resolve each species to a tip of the Upham et al. (2019) posterior
#      tree sample, using only synonym tables already recorded in the
#      Evo-M1-Trait-Data repo.
#   3. Optionally restrict the set to one or more taxonomic orders
#      (default: Primates only).
#   4. Prune all 100 trees to that species set and relabel the tips with
#      the names used by the analysis data.
#   5. Write the tree sample, the source-tree IDs, and a coverage report
#      saying how every input species was handled.
#
# Runs directly after s3_prepare_Stephan_primate_data.R, which writes the
# default input (data_intermediate/s3_Stephan_species_list.csv), and before
# the s3 PGLS scripts, which read the tree sample this script produces.
#
# This script performs data preparation only.
# It does not fit models, run analyses, or create plots.
#
# NOTHING IS GRAFTED OR IMPUTED LOCALLY. A species that no published tip
# carries is reported absent and left off the trees. The Completed100 set's
# no-DNA placements are Upham et al.'s own published imputation, made afresh
# in each posterior tree, so that uncertainty propagates into any fit run
# across the sample - which is the point of using 100 trees rather than one.
#
# ---- USAGE -------------------------------------------------------------
#   Rscript scripts/s3_prepare_tree_sample.R
#
#   Override any setting in the CONFIG block from the command line:
#     --species=data_intermediate/other_species_list.csv
#     --column=accepted_name        # which column holds the species names
#     --orders=Primates             # comma-separated, or "all" for no filter
#     --out-prefix=s3v2_tree        # prefix for every output file
#     --tree-repo=/path/to/Evo-M1-Trait-Data
#     --source=nexus|repo_sample    # see TREE SOURCE below
#
#   Examples:
#     # all mammals in a different list, outputs named s4_tree_*
#     Rscript scripts/s3_prepare_tree_sample.R \
#       --species=data_intermediate/s4_species_list.csv \
#       --orders=all --out-prefix=s4_tree
#
#     # primates and scandentians from a plain text list
#     Rscript scripts/s3_prepare_tree_sample.R \
#       --species=species_list_Stephan_primates.txt \
#       --orders=Primates,Scandentia
#
# Requires: ape, here
# ============================================================

suppressWarnings(suppressMessages(library(ape)))
suppressWarnings(suppressMessages(library(here)))

options(scipen = 999)

# ============================================================
# CONFIG - edit here, or override on the command line
# ============================================================
CONFIG <- list(
  # Species list to build the trees for. Relative paths resolve against the
  # project root. .csv / .tsv / .txt all work; one species per row.
  species_file = file.path("data_intermediate", "s3_Stephan_species_list.csv"),

  # Column holding the binomials. NULL = autodetect (see SPECIES_COLUMNS).
  # Headerless single-column text files are handled automatically.
  species_column = NULL,

  # Orders to keep. Character vector, or NULL / "all" to keep every species
  # that reaches a tip. This is the primates-only trim.
  keep_orders = "Primates",

  # Prefix for every output file written into data_intermediate/.
  out_prefix = "s3_tree",

  # Evo-M1-Trait-Data repo: holds the frozen Upham trees and the synonym
  # tables. First existing path wins; --tree-repo overrides.
  tree_repo_candidates = c(
    "~/Library/CloudStorage/OneDrive-AllenInstitute/Species/Evo-M1-Trait-Data",
    "~/OneDrive - Allen Institute/Species/Evo-M1-Trait-Data",
    "~/Evo-M1-Trait-Data"
  ),

  # TREE SOURCE
  #   "nexus"       Upham_etal_2019/Completed100_topoCons_NDexp/output.nex
  #                 -- the full 5,911-tip published sample. Use this. It is
  #                 the only source that can reach a species the Evo-M1
  #                 project itself does not carry, and 3 of the 59 Stephan
  #                 primates are exactly that case (see NOTE below).
  #   "repo_sample" _keys/mammal_trees_sample100.nwk -- the same 100 trees
  #                 already pruned to the 203 Evo-M1 project species and
  #                 relabelled. Small and quick, but anything outside that
  #                 203 is unreachable. Fallback only.
  tree_source = "nexus",

  # Minimum tips for a usable tree.
  min_tips = 4L
)

# Candidate species-column names, in priority order, for autodetection.
SPECIES_COLUMNS <- c("Species", "species", "accepted_name", "Accepted_name",
                     "binomial", "Binomial", "taxon", "Taxon",
                     "tip_label", "species_sci", "name", "Name")

# ============================================================
# COMMAND-LINE OVERRIDES
# ============================================================
cli_args <- commandArgs(trailingOnly = TRUE)
cli_opt <- function(flag) {
  hit <- grep(paste0("^", flag, "="), cli_args, value = TRUE)
  if (!length(hit)) return(NULL)
  sub(paste0("^", flag, "="), "", tail(hit, 1))
}

if (!is.null(v <- cli_opt("--species")))    CONFIG$species_file   <- v
if (!is.null(v <- cli_opt("--column")))     CONFIG$species_column <- v
if (!is.null(v <- cli_opt("--out-prefix"))) CONFIG$out_prefix     <- v
if (!is.null(v <- cli_opt("--source")))     CONFIG$tree_source    <- v
if (!is.null(v <- cli_opt("--tree-repo")))  CONFIG$tree_repo_candidates <- v
if (!is.null(v <- cli_opt("--orders"))) {
  CONFIG$keep_orders <- trimws(strsplit(v, ",")[[1]])
  CONFIG$keep_orders <- CONFIG$keep_orders[nzchar(CONFIG$keep_orders)]
}

filter_orders <- !(is.null(CONFIG$keep_orders) ||
                   length(CONFIG$keep_orders) == 0 ||
                   identical(tolower(CONFIG$keep_orders), "all"))

if (!CONFIG$tree_source %in% c("nexus", "repo_sample"))
  stop("tree_source must be \"nexus\" or \"repo_sample\", not \"",
       CONFIG$tree_source, "\"")

# ============================================================
# FILE PATHS
# ============================================================
resolve_in_project <- function(p) if (grepl("^([/~]|[A-Za-z]:)", p)) path.expand(p) else here(p)

input_species <- resolve_in_project(CONFIG$species_file)

repo <- NULL
for (cand in CONFIG$tree_repo_candidates) {
  if (dir.exists(path.expand(cand))) { repo <- normalizePath(path.expand(cand)); break }
}
if (is.null(repo))
  stop("Evo-M1-Trait-Data repo not found. Tried:\n  ",
       paste(CONFIG$tree_repo_candidates, collapse = "\n  "),
       "\nPass --tree-repo=/path/to/Evo-M1-Trait-Data")

input_nexus     <- file.path(repo, "Upham_etal_2019", "Completed100_topoCons_NDexp", "output.nex")
input_repo_nwk  <- file.path(repo, "_keys", "mammal_trees_sample100.nwk")
input_crosswalk <- file.path(repo, "__merging_trees", "tree_tip_crosswalk.csv")
input_aliases   <- file.path(repo, "_keys", "species_display_aliases.csv")
input_reference <- file.path(repo, "_keys", "species_reference.csv")
input_upham_tips <- file.path(repo, "Upham_etal_2019", "Upham_etal_2019_DNAonlyMCC.csv")
input_sample_ids <- file.path(repo, "__merging_trees", "tree_sample_ids.csv")

output_dir <- here("data_intermediate")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

pre <- CONFIG$out_prefix
output_trees_nwk   <- file.path(output_dir, paste0(pre, "_sample100.nwk"))
output_trees_rds   <- file.path(output_dir, paste0(pre, "_sample100.rds"))
output_source_nwk  <- file.path(output_dir, paste0(pre, "_sample100_sourcelabels.nwk"))
output_ids_csv     <- file.path(output_dir, paste0(pre, "_sample_ids.csv"))
output_report_csv  <- file.path(output_dir, paste0(pre, "_coverage_report.csv"))

# ============================================================
# LOAD AND CHECK INPUTS
# ============================================================
if (!file.exists(input_species)) stop("File not found: ", input_species)

tree_file <- if (CONFIG$tree_source == "nexus") input_nexus else input_repo_nwk
if (!file.exists(tree_file))
  stop("Tree source not found: ", tree_file,
       "\n  See ", file.path(repo, "__merging_trees", "__merging_trees.README.md"))
if (!file.exists(input_crosswalk)) stop("File not found: ", input_crosswalk)

read_species_file <- function(path, column) {
  sep <- if (grepl("\\.tsv$|\\.txt$", path, ignore.case = TRUE)) "\t" else ","
  tbl <- tryCatch(
    read.csv(path, sep = sep, stringsAsFactors = FALSE,
             check.names = FALSE, fileEncoding = "UTF-8-BOM"),
    error = function(e)
      read.csv(path, sep = sep, stringsAsFactors = FALSE, check.names = FALSE)
  )
  # Headerless single-column list: the header row is itself a species name.
  if (ncol(tbl) == 1L && grepl("^[A-Z][a-z]+[ _][a-z]", names(tbl)[1])) {
    tbl <- data.frame(species = c(names(tbl)[1], tbl[[1]]), stringsAsFactors = FALSE)
    column <- "species"
  }
  if (is.null(column)) {
    hit <- intersect(SPECIES_COLUMNS, names(tbl))
    column <- if (length(hit)) hit[1] else names(tbl)[1]
  }
  if (!column %in% names(tbl))
    stop("Column \"", column, "\" not found in ", basename(path),
         ". Columns present: ", paste(names(tbl), collapse = ", "),
         "\n  Pass --column=<name>")
  list(values = tbl[[column]], column = column)
}

sp_read <- read_species_file(input_species, CONFIG$species_column)

# Keep the label exactly as supplied (underscored), so the tree tips join
# to the analysis data with no further renaming. Match on the spaced form.
raw <- trimws(gsub("[\r\n]", "", as.character(sp_read$values)))
raw <- raw[!is.na(raw) & nzchar(raw)]
label_in <- gsub("[[:space:]]+", "_", raw)
binom_in <- gsub("_+", " ", label_in)
keep_first <- !duplicated(tolower(binom_in))
label_in <- label_in[keep_first]
binom_in <- binom_in[keep_first]

if (!length(binom_in)) stop("No species names read from ", input_species)

cat("species list  :", input_species, "\n  column      :", sp_read$column,
    "\n  species     :", length(binom_in), "\n")

# ---- trees -------------------------------------------------------------
trees <- if (CONFIG$tree_source == "nexus") read.nexus(tree_file) else read.tree(tree_file)
if (!inherits(trees, "multiPhylo")) trees <- c(trees)

tipset1 <- sort(trees[[1]]$tip.label)
same <- vapply(trees, function(t) identical(sort(t$tip.label), tipset1), logical(1))
if (!all(same))
  stop("Tree ", which(!same)[1], " has a different tip set from tree 1 -- ",
       "the sample is not one consistent taxon set.")

cat("tree source   :", tree_file, "\n  trees       :", length(trees),
    "\n  tips/tree   :", length(trees[[1]]$tip.label), "\n")

# tip "Genus_species[_FAMILY_ORDER]" -> "Genus species"
tip_binom <- function(x) {
  p <- strsplit(gsub("[[:space:]]+", "_", trimws(x)), "_+")
  vapply(p, function(v) {
    v <- v[nzchar(v)]
    if (length(v) >= 2) paste(v[1], v[2]) else paste(v, collapse = " ")
  }, character(1))
}
tip_of <- setNames(trees[[1]]$tip.label, tolower(tip_binom(trees[[1]]$tip.label)))
tip_of <- tip_of[!duplicated(names(tip_of))]

# ============================================================
# SYNONYM TABLES
#
# Every candidate spelling comes from a table already in the repo; this
# script invents none. Sources, in the order they are consulted:
#   1. tree_tip_crosswalk.csv   accepted_name -> candidate, gated on
#                               auto_match == TRUE. auto_match == FALSE is a
#                               DIFFERENT taxon concept and is never used.
#   2. species_display_aliases.csv    variant <-> canonical, both directions
#   3. species_reference.csv          accepted_name <-> mdd_accepted_name
# ============================================================
syn <- list()
add_syn <- function(from, to, source) {
  from <- tolower(trimws(from)); to <- trimws(to)
  ok <- nzchar(from) & nzchar(to) & !is.na(from) & !is.na(to) & to != "NA"
  for (i in which(ok))
    syn[[from[i]]] <<- rbind(syn[[from[i]]],
                             data.frame(candidate = to[i], source = source,
                                        stringsAsFactors = FALSE))
}

xw <- read.csv(input_crosswalk, stringsAsFactors = FALSE, encoding = "UTF-8")
xw_auto <- xw[toupper(trimws(xw$auto_match)) == "TRUE" & nzchar(xw$candidate_name), ]
xw_auto <- xw_auto[order(suppressWarnings(as.integer(xw_auto$candidate_rank))), ]
add_syn(xw_auto$accepted_name, xw_auto$candidate_name, "tree_tip_crosswalk")

n_review <- sum(toupper(trimws(xw$auto_match)) == "FALSE" & nzchar(xw$candidate_name))

if (file.exists(input_aliases)) {
  al <- read.csv(input_aliases, stringsAsFactors = FALSE, encoding = "UTF-8")
  add_syn(al$canonical, al$variant, "species_display_aliases")
  add_syn(al$variant, al$canonical, "species_display_aliases")
} else warning("Alias table not found, skipping: ", input_aliases)

if (file.exists(input_reference)) {
  rf <- read.csv(input_reference, stringsAsFactors = FALSE, encoding = "UTF-8")
  add_syn(rf$mdd_accepted_name, rf$accepted_name, "species_reference")
  add_syn(rf$accepted_name, rf$mdd_accepted_name, "species_reference")
} else {
  rf <- NULL
  warning("Species reference not found, skipping: ", input_reference)
}

# ============================================================
# ORDER LOOKUP (for the taxonomic trim)
#
# Consulted in this order, for both the input name and the tip it resolved
# to, then as a last resort the genus. order_source records which answered.
# ============================================================
ord_maps <- list()
addmap <- function(name, keys, vals) {
  keys <- tolower(trimws(keys)); vals <- trimws(vals)
  ok <- nzchar(keys) & nzchar(vals) & !is.na(vals) & vals != "NA"
  m <- setNames(vals[ok], keys[ok])
  ord_maps[[name]] <<- m[!duplicated(names(m))]
}
addmap("tree_tip_crosswalk", xw$accepted_name, xw$order_expected)
if (!is.null(rf)) addmap("species_reference", rf$accepted_name, rf$Order_resolved)

genus_order <- character(0)
if (file.exists(input_upham_tips)) {
  ut <- read.csv(input_upham_tips, stringsAsFactors = FALSE, encoding = "UTF-8")
  addmap("Upham_tip_table", ut$species_sci, ut$Order_tip)
  g <- setNames(trimws(ut$Order_tip), tolower(trimws(ut$Genus)))
  g <- g[nzchar(names(g)) & nzchar(g) & g != "NA"]
  genus_order <- g[!duplicated(names(g))]
} else warning("Upham tip table not found, order lookup will be thinner: ",
               input_upham_tips)

lookup_order <- function(...) {
  keys <- tolower(trimws(unlist(list(...))))
  keys <- keys[nzchar(keys) & !is.na(keys)]
  for (nm in names(ord_maps)) {
    for (k in keys) {
      v <- ord_maps[[nm]][k]
      if (!is.na(v)) return(list(order = unname(v), source = nm))
    }
  }
  for (k in keys) {
    gk <- sub(" .*$", "", k)
    v <- genus_order[gk]
    if (!is.na(v)) return(list(order = unname(v), source = "Upham_genus"))
  }
  list(order = NA_character_, source = NA_character_)
}

# ============================================================
# RESOLVE EACH SPECIES TO ONE TIP
# ============================================================
resolve_one <- function(binom) {
  direct <- tip_of[tolower(binom)]
  if (!is.na(direct))
    return(list(tip = unname(direct), via = binom, source = "accepted_name",
                status = "matched_direct"))
  cands <- syn[[tolower(binom)]]
  if (!is.null(cands)) {
    cands <- cands[!duplicated(tolower(cands$candidate)), , drop = FALSE]
    for (i in seq_len(nrow(cands))) {
      hit <- tip_of[tolower(cands$candidate[i])]
      if (!is.na(hit))
        return(list(tip = unname(hit), via = cands$candidate[i],
                    source = cands$source[i], status = "matched_synonym"))
    }
  }
  list(tip = NA_character_, via = "", source = "", status = "absent_from_source_tree")
}

resolved <- lapply(binom_in, resolve_one)

# Direct matches claim a shared tip ahead of synonym matches, so the outcome
# does not depend on row order in the species file.
prio <- vapply(resolved, function(r)
  if (is.na(r$tip)) 2L else if (r$status == "matched_direct") 0L else 1L, integer(1))
claim_order <- order(prio, tolower(binom_in))

rep_rows <- vector("list", length(binom_in))
keep_tip <- character(0)
keep_lab <- character(0)
used <- character(0)

for (i in claim_order) {
  r   <- resolved[[i]]
  ord <- lookup_order(binom_in[i], r$tip, if (nzchar(r$via)) r$via else NULL)
  status <- r$status
  note <- ""
  tip_used <- ""

  if (!is.na(r$tip)) {
    if (r$tip %in% names(used)) {
      status <- "tip_conflict_not_used"
      note <- paste0("tip ", r$tip, " already assigned to ", used[[r$tip]],
                     "; a species-level tree cannot separate them. Resolve in ",
                     "tree_tip_crosswalk.csv / species_display_aliases.csv ",
                     "before relying on this run.")
    } else if (filter_orders &&
               (is.na(ord$order) || !(ord$order %in% CONFIG$keep_orders))) {
      status <- if (is.na(ord$order)) "excluded_order_unknown" else "excluded_other_order"
      note <- if (is.na(ord$order))
        "no order found in any repo table; excluded because an order filter is active"
      else paste0("order ", ord$order, " is not in keep_orders (",
                  paste(CONFIG$keep_orders, collapse = ", "), ")")
    } else {
      used[r$tip] <- binom_in[i]
      keep_tip <- c(keep_tip, r$tip)
      keep_lab <- c(keep_lab, label_in[i])
      tip_used <- r$tip
      if (status == "matched_synonym")
        note <- paste0("reached via \"", r$via, "\" (", r$source, ")")
    }
  } else {
    note <- paste0("no spelling recorded in the repo synonym tables occurs in ",
                   "the source tree; reported absent, never grafted",
                   if (n_review > 0)
                     paste0(". ", n_review, " review-only crosswalk candidates ",
                            "exist repo-wide and are deliberately not used")
                   else "")
  }

  rep_rows[[i]] <- data.frame(
    input_species    = label_in[i],
    binomial         = binom_in[i],
    status           = status,
    on_trees         = nzchar(tip_used),
    matched_tip      = if (is.na(r$tip)) "" else r$tip,
    matched_via      = r$via,
    match_source     = r$source,
    order_assigned   = if (is.na(ord$order)) "" else ord$order,
    order_source     = if (is.na(ord$source)) "" else ord$source,
    note             = note,
    stringsAsFactors = FALSE
  )
}

report <- do.call(rbind, rep_rows)

conflicts <- report$input_species[report$status == "tip_conflict_not_used"]
if (length(conflicts))
  warning("Tip conflict - two input species claim one tip: ",
          paste(conflicts, collapse = ", "),
          ". Fix the crosswalk; these rows were excluded.")

absent <- report$input_species[report$status == "absent_from_source_tree"]
if (length(absent))
  warning(length(absent), " species are not in the source tree and were left off: ",
          paste(absent, collapse = ", "))

unknown_order <- report$input_species[report$status == "excluded_order_unknown"]
if (length(unknown_order))
  warning(length(unknown_order), " species reached a tip but have no order in any ",
          "repo table, so the order filter excluded them: ",
          paste(unknown_order, collapse = ", "))

if (length(keep_tip) < CONFIG$min_tips)
  stop("Only ", length(keep_tip), " species reached the trees; need at least ",
       CONFIG$min_tips, ". Check the species column and the order filter.")

# ============================================================
# PRUNE EVERY TREE AND RELABEL
# ============================================================
relab <- setNames(keep_lab, keep_tip)

pruned_src <- lapply(unclass(trees), function(t) {
  p <- keep.tip(t, keep_tip)
  p$node.label <- NULL
  p
})
class(pruned_src) <- "multiPhylo"

pruned_out <- lapply(pruned_src, function(p) {
  p$tip.label <- unname(relab[p$tip.label]); p
})
class(pruned_out) <- "multiPhylo"

stopifnot(all(vapply(pruned_out, function(p)
  setequal(p$tip.label, keep_lab), logical(1))))

# ---- source-tree IDs ---------------------------------------------------
tree_ids <- if (!is.null(names(trees))) names(trees) else rep(NA_character_, length(trees))
if (all(is.na(tree_ids)) && file.exists(input_sample_ids)) {
  ids_tbl <- read.csv(input_sample_ids, stringsAsFactors = FALSE)
  if (nrow(ids_tbl) == length(trees)) tree_ids <- ids_tbl$tree_id
}
if (all(is.na(tree_ids))) tree_ids <- paste0("tree_", seq_along(trees))

ids <- data.frame(line = seq_along(trees), tree_id = tree_ids,
                  stringsAsFactors = FALSE)

# ---- ultrametricity check ---------------------------------------------
spreads <- vapply(pruned_out, function(p) {
  d <- node.depth.edgelength(p)[seq_len(Ntip(p))]
  max(d) - min(d)
}, numeric(1))
maxd <- vapply(pruned_out, function(p)
  max(node.depth.edgelength(p)[seq_len(Ntip(p))]), numeric(1))
w <- which.max(spreads)
ultrametric <- spreads[w] <= 1e-3 * maxd[w]
if (!ultrametric)
  warning("Pruned trees are NOT ultrametric (worst root-to-tip spread ",
          signif(spreads[w], 4), " on line ", w, ") -- check the source.")

# ============================================================
# SAVE OUTPUTS
# ============================================================
write.tree(pruned_out, output_trees_nwk)
saveRDS(pruned_out, output_trees_rds)
write.tree(pruned_src, output_source_nwk)
write.csv(ids, output_ids_csv, row.names = FALSE, quote = FALSE,
          fileEncoding = "UTF-8")
write.csv(report, output_report_csv, row.names = FALSE, na = "",
          fileEncoding = "UTF-8")

# ============================================================
# SUMMARY
# ============================================================
cat("\norder filter  :", if (filter_orders) paste(CONFIG$keep_orders, collapse = ", ") else "none (all orders)", "\n")
cat("on trees      :", length(keep_tip), "/", nrow(report), "input species\n")
tb <- sort(table(report$status), decreasing = TRUE)
for (i in seq_along(tb)) cat(sprintf("  %4d  %s\n", tb[i], names(tb)[i]))
cat(sprintf("\nroot age (tree 1): %.4g Ma   worst root-to-tip spread: %.3g  (%s)\n",
            maxd[1], spreads[w], if (ultrametric) "ultrametric" else "NOT ultrametric"))

message("Saved tree sample: ", output_trees_nwk)
message("Saved tree sample (RDS): ", output_trees_rds)
message("Saved published-label trees: ", output_source_nwk)
message("Saved source-tree IDs: ", output_ids_csv)
message("Saved coverage report: ", output_report_csv)

# ---- how the s3 PGLS scripts read this --------------------------------
#   trees <- ape::read.tree("data_intermediate/s3_tree_sample100.nwk")
#   # or  readRDS("data_intermediate/s3_tree_sample100.rds")
#   have  <- gsub(" ", "_", my_data$Species)
#   sub   <- lapply(trees, ape::keep.tip,
#                   tip = intersect(trees[[1]]$tip.label, have))
#   class(sub) <- "multiPhylo"
#   fits  <- lapply(sub, function(tr) <PGLS fit against tr>)
#
# Summarise across the 100 fits (median + spread); never cherry-pick a line.
# Cite which posterior trees they are via the *_sample_ids.csv table.
