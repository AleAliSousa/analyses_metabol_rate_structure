library(jsonlite)
library(tidyverse)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

flatten_allen_node <- function(node,
                               parent_id = NA_integer_,
                               parent_name = NA_character_,
                               parent_acronym = NA_character_,
                               name_path = character(),
                               acronym_path = character(),
                               depth = 0L) {
  
  this_id <- node$id %||% NA_integer_
  this_name <- node$name %||% NA_character_
  this_acronym <- node$acronym %||% NA_character_
  
  current_name_path <- c(name_path, this_name)
  current_acronym_path <- c(acronym_path, this_acronym)
  
  children <- node$children %||% list()
  
  this_row <- tibble(
    structure_id = this_id,
    acronym = this_acronym,
    name = this_name,
    parent_structure_id = parent_id,
    parent_acronym = parent_acronym,
    parent_name = parent_name,
    depth = depth,
    n_children = length(children),
    name_path = paste(current_name_path, collapse = " > "),
    acronym_path = paste(current_acronym_path, collapse = " > ")
  )
  
  child_rows <- map_dfr(
    children,
    flatten_allen_node,
    parent_id = this_id,
    parent_name = this_name,
    parent_acronym = this_acronym,
    name_path = current_name_path,
    acronym_path = current_acronym_path,
    depth = depth + 1L
  )
  
  bind_rows(this_row, child_rows)
}

allen_json <- fromJSON(
  "https://api.brain-map.org/api/v2/structure_graph_download/10.json",
  simplifyVector = FALSE
)

allen_regions <- map_dfr(allen_json$msg, flatten_allen_node) %>%
  arrange(name_path)

View(allen_regions)