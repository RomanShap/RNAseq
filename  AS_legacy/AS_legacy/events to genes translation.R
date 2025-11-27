suppressMessages({ library(dplyr); library(tibble); library(purrr) })

# Build EVENT -> GENE lookup from any number of splicing tables
# conflict = "first" keeps the first gene per EVENT; 
# conflict = "concat" concatenates multiple gene symbols with ';'
build_ev2gene <- function(..., conflict = c("first","concat")) {
  conflict <- match.arg(conflict)
  df <- bind_rows(...)
  if (!all(c("EVENT","GENE") %in% names(df))) {
    stop("Input tables must contain columns: EVENT, GENE")
  }
  df2 <- df %>%
    distinct(EVENT, GENE)
  if (conflict == "first") {
    df2 <- df2 %>% group_by(EVENT) %>% summarise(GENE = dplyr::first(GENE), .groups = "drop")
  } else {
    df2 <- df2 %>% group_by(EVENT) %>% summarise(GENE = paste(unique(GENE), collapse = ";"), .groups = "drop")
  }
  tibble::deframe(df2)  # named char vector: names = EVENT, values = GENE
}

# Translate one vector of EVENT IDs using an EVENT->GENE dictionary
translate_vec <- function(v, ev2gene, mode = c("label","gene","data")) {
  mode <- match.arg(mode)
  g <- unname(ev2gene[v])  # NA if not found
  if (mode == "gene") {
    return(unique(ifelse(is.na(g), v, g)))
  } else if (mode == "label") {
    return(unique(ifelse(is.na(g), v, paste0(v, " - ", g))))
  } else {
    return(tibble(EVENT = v, GENE = g))
  }
}

# Translate the entire partitioned_events structure across ALL event types
# Keeps the same nested structure and group names.
translate_partitioned_all <- function(partitioned_events, ev2gene, mode = c("label","gene","data")) {
  mode <- match.arg(mode)
  # For each event type (e.g., "INT", "ALTA", "SE", ...)
  out <- lapply(partitioned_events, function(by_group) {
    # For each group (e.g., "FUS_EWSR1", "EWSR1_TAF15", etc.)
    lapply(by_group, function(vec_ids) translate_vec(vec_ids, ev2gene, mode = mode))
  })
  out
}

# ---- Usage ----
# 1) Build the dictionary from ALL three (or more) tables
ev2gene <- build_ev2gene(EWSR1_splicing_table, FUS_splicing_table, TAF15_splicing_table, conflict = "first")

# 2a) Replace IDs by "EVENT - GENE" labels everywhere
partitioned_events_labeled <- translate_partitioned_all(partitioned_events, ev2gene, mode = "label")

# 2b) Replace IDs by GENE symbols only (fallback to EVENT when missing)
partitioned_events_genes <- translate_partitioned_all(partitioned_events, ev2gene, mode = "gene")

# 2c) Get data frames with EVENT/GENE per group (structure preserved)
partitioned_events_df <- translate_partitioned_all(partitioned_events, ev2gene, mode = "data")
