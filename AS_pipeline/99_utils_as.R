# 99_utils_as.R
#
# Splicing-specific helper functions for the alternative splicing pipeline.
# These functions complement the general utilities provided in shared/99_utils.R.
# They provide routines for reading VAST‑TOOLS output, filtering splicing events,
# splitting by event type and direction, computing summary metrics and building
# membership matrices for UpSet diagrams.  All functions are written in
# modular form and use tidyverse packages for ease of manipulation.

## Load required packages ------------------------------------------------------
suppressed_pkgs <- c("dplyr", "readr", "stringr", "purrr", "tidyr", "ggplot2")
suppressMessages(lapply(suppressed_pkgs, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed", pkg))
  }
  library(pkg, character.only = TRUE)
}))

## read_vast_tools_tables ------------------------------------------------------
#' Read VAST‑TOOLS splicing tables for each condition.
#'
#' This function expects that `vast_tools_dir` contains one subdirectory per
#' condition (matching the names of `sets` in the configuration).  Each
#' subdirectory should include at least one TSV file containing splicing
#' metrics (e.g. `inclusion_levels.tsv` or `MVs_scores.tsv`).  If multiple TSV
#' files are present the one with `inclusion` in the name is selected.  The
#' function reads all tables into a named list of data frames.
#'
#' @param vast_dir Character.  Path to the directory with VAST‑TOOLS results.
#' @param sets A named list where names are conditions and values are sample IDs.
#' @return A named list of data frames keyed by condition.
read_vast_tools_tables <- function(vast_dir, sets) {
  cond_names <- names(sets)
  tables <- map(cond_names, function(cond) {
    cond_dir <- file.path(vast_dir, cond)
    if (!dir.exists(cond_dir)) {
      stop(sprintf("Condition directory not found: %s", cond_dir))
    }
    tsv_files <- list.files(cond_dir, pattern = "\\.tsv$", full.names = TRUE)
    if (length(tsv_files) == 0) {
      stop(sprintf("No TSV files found in %s", cond_dir))
    }
    # Choose a file containing 'inclusion' if available; otherwise take the first
    file <- tsv_files[1]
    incl_files <- grep("inclusion", basename(tsv_files), value = TRUE, ignore.case = TRUE)
    if (length(incl_files) > 0) {
      file <- file.path(cond_dir, incl_files[1])
    }
    message(sprintf("Reading VAST‑TOOLS table for %s: %s", cond, file))
    df <- suppressMessages(readr::read_tsv(file, col_types = cols()))
    return(df)
  })
  names(tables) <- cond_names
  return(tables)
}

## filter_splicing_events ------------------------------------------------------
#' Filter splicing events by MV and ΔPSI thresholds and expression status.
#'
#' Each event table is expected to contain at least the following columns:
#' * `MV.dPsi._at_0.95` – the multi‑variate ΔPSI threshold at 95% confidence.
#' * `E.dPsi.` or `dPsi` – the observed ΔPSI.
#' * `EVENT` – the event identifier (e.g. HsaEX001234).
#' * Gene identifier (`GeneName` or `GENEID`).
#' Additional columns holding per‑sample PSI values are allowed and will be
#' retained.  Events are kept if they meet all of the following criteria:
#' 1. `MV.dPsi._at_0.95 >= mv_threshold` (if the column exists).
#' 2. `abs(E.dPsi. or dPsi) >= dpsi_threshold` (if the column exists).
#' 3. At least `min_valid_calls` replicate PSI values are non‑NA.
#' 4. If `expressed_only` is TRUE and `expressed_genes` is supplied, the event's
#'    gene must be present in `expressed_genes`.
#'
#' @param df Data frame of events for a single condition.
#' @param mv_threshold Numeric.  Minimum MV.dPsi threshold.
#' @param dpsi_threshold Numeric.  Minimum absolute ΔPSI.
#' @param expressed_only Logical.  Whether to restrict to events in expressed genes.
#' @param expressed_genes Character vector of gene names that are reliably expressed.
#' @param min_valid_calls Integer.  Minimum number of non‑NA PSI calls across replicates.
#' @return Filtered data frame of splicing events.
filter_splicing_events <- function(df,
                                   mv_threshold,
                                   dpsi_threshold,
                                   expressed_only = FALSE,
                                   expressed_genes = NULL,
                                   min_valid_calls = 2) {
  # Ensure event identifier exists
  if (!"EVENT" %in% names(df)) {
    stop("Input data frame lacks column 'EVENT'")
  }
  # Determine ΔPSI column
  dpsi_col <- intersect(c("E.dPsi.", "dPsi", "DeltaPSI"), names(df))
  mv_col <- intersect(c("MV.dPsi._at_0.95", "MV"), names(df))
  gene_col <- intersect(c("GeneName", "GENEID", "GENE"), names(df))
  if (length(dpsi_col) == 0) {
    warning("No dPSI column detected; skipping ΔPSI filtering")
    dpsi_col <- NULL
  }
  if (length(mv_col) == 0) {
    warning("No MV column detected; skipping MV filtering")
    mv_col <- NULL
  }
  # Count non‑NA PSI calls across replicates
  psi_cols <- grep("_PSI$|_psi$|PSI_", names(df), value = TRUE)
  # If no PSI columns are found try columns starting with sample IDs
  if (length(psi_cols) == 0) {
    psi_cols <- names(df)[grepl("^" , names(df))]
  }
  df <- df %>% mutate(non_na_calls = if (length(psi_cols) > 0) {
    rowSums(!is.na(select(., all_of(psi_cols))))
  } else {
    NA_integer_
  })
  # Apply filters
  filtered <- df %>%
    filter(
      (is.null(mv_col) || .data[[mv_col[1]]] >= mv_threshold),
      (is.null(dpsi_col) || abs(.data[[dpsi_col[1]]]) >= dpsi_threshold),
      is.na(non_na_calls) | non_na_calls >= min_valid_calls,
      (!expressed_only || is.null(expressed_genes) || .data[[gene_col[1]]] %in% expressed_genes)
    ) %>%
    select(-non_na_calls)
  return(filtered)
}

## split_event_types -----------------------------------------------------------
#' Split a splicing table into categories by event type.
#'
#' Event identifiers typically begin with a prefix indicating the type (e.g.
#' "HsaEX" for exon skipping, "HsaINT" for intron retention, "HsaALTA" for
#' alternative acceptor and "HsaALTD" for alternative donor).  This function
#' uses regular expressions to classify events into EX, INT, ALTA, ALTD.  Any
#' event not matching these patterns is assigned to the "OTHER" category.
#'
#' @param df Data frame of splicing events.
#' @return A named list of data frames, one per event type.
split_event_types <- function(df) {
  df <- df %>% mutate(
    event_type = case_when(
      str_detect(EVENT, "EX") ~ "EX",
      str_detect(EVENT, "INT") ~ "INT",
      str_detect(EVENT, "ALTA") ~ "ALTA",
      str_detect(EVENT, "ALTD") ~ "ALTD",
      TRUE ~ "OTHER"
    )
  )
  split(df, df$event_type)
}

## split_by_direction ----------------------------------------------------------
#' Divide events into all, positive and negative ΔPSI subsets.
#'
#' Uses the first available ΔPSI column to determine the sign.  Events with
#' missing ΔPSI are included only in the "all" subset.
#'
#' @param df Data frame of splicing events.
#' @return A list with elements `all`, `pos` and `neg`.
split_by_direction <- function(df) {
  dpsi_col <- intersect(c("E.dPsi.", "dPsi", "DeltaPSI"), names(df))
  if (length(dpsi_col) == 0) {
    warning("No dPSI column detected; returning only 'all' subset")
    return(list(all = df, pos = df[0, ], neg = df[0, ]))
  }
  dcol <- dpsi_col[1]
  list(
    all = df,
    pos = df %>% filter(.data[[dcol]] > 0),
    neg = df %>% filter(.data[[dcol]] < 0)
  )
}

## compute_splicing_metrics -----------------------------------------------------
#' Compute summary metrics for a filtered splicing table.
#'
#' Produces a tibble with counts of events, counts of positive/negative
#' ΔPSI events, mean and median MV and |ΔPSI| values and any other simple
#' statistics that may be useful for QC.
#'
#' @param df Data frame of filtered splicing events.
#' @return A one‑row tibble of summary metrics.
compute_splicing_metrics <- function(df) {
  dpsi_col <- intersect(c("E.dPsi.", "dPsi", "DeltaPSI"), names(df))
  mv_col <- intersect(c("MV.dPsi._at_0.95", "MV"), names(df))
  n_total <- nrow(df)
  n_pos <- if (length(dpsi_col) > 0) sum(df[[dpsi_col[1]]] > 0, na.rm = TRUE) else NA_integer_
  n_neg <- if (length(dpsi_col) > 0) sum(df[[dpsi_col[1]]] < 0, na.rm = TRUE) else NA_integer_
  mean_mv <- if (length(mv_col) > 0) mean(df[[mv_col[1]]], na.rm = TRUE) else NA_real_
  median_mv <- if (length(mv_col) > 0) median(df[[mv_col[1]]], na.rm = TRUE) else NA_real_
  mean_abs_dpsi <- if (length(dpsi_col) > 0) mean(abs(df[[dpsi_col[1]]]), na.rm = TRUE) else NA_real_
  median_abs_dpsi <- if (length(dpsi_col) > 0) median(abs(df[[dpsi_col[1]]]), na.rm = TRUE) else NA_real_
  tibble::tibble(
    total_events = n_total,
    positive_events = n_pos,
    negative_events = n_neg,
    mean_mv = mean_mv,
    median_mv = median_mv,
    mean_abs_dpsi = mean_abs_dpsi,
    median_abs_dpsi = median_abs_dpsi
  )
}

## build_event_membership_matrix ------------------------------------------------
#' Build a binary membership matrix from a list of event identifiers.
#'
#' Given a named list of event ID vectors (typically one per condition),
#' constructs a data frame with rows corresponding to all unique event IDs and
#' columns corresponding to conditions.  Values are 1 if the event is present
#' in that condition and 0 otherwise.  The resulting matrix can be passed to
#' plotting functions such as `plot_upset_complex()` from shared/99_utils.R.
#'
#' @param event_sets Named list of character vectors (event IDs).
#' @return A data frame of 0/1 membership (events × conditions).
build_event_membership_matrix <- function(event_sets) {
  all_events <- unique(unlist(event_sets))
  mat <- lapply(event_sets, function(ev) as.integer(all_events %in% ev)) %>%
    bind_cols() %>%
    setNames(names(event_sets))
  mat <- cbind(EVENT = all_events, mat)
  return(as.data.frame(mat))
}

## plot_event_upset -------------------------------------------------------------
#' Create an UpSet plot for splicing events.
#'
#' This is a convenience wrapper around `plot_upset_complex()` defined in
#' shared/99_utils.R.  It accepts a binary membership matrix (as produced by
#' `build_event_membership_matrix()`) and constructs an UpSet plot with
#' default options suitable for splicing events.
#'
#' @param membership_matrix Data frame with event IDs and binary membership.
#' @param title Character.  Title for the plot.
#' @param output_path Character.  File path to save the plot (PNG format).
#' @param colors Named vector of colours for conditions.
#' @return Invisibly returns the ggplot object.
plot_event_upset <- function(membership_matrix,
                             title = "Splicing event overlaps",
                             output_path = NULL,
                             colors = NULL) {
  # Exclude the EVENT column for plotting
  mat <- membership_matrix[, -1, drop = FALSE]
  # Use colours if provided
  if (!is.null(colors)) {
    plot <- plot_upset_complex(mat, title = title, set_colors = colors)
  } else {
    plot <- plot_upset_complex(mat, title = title)
  }
  if (!is.null(output_path)) {
    ggsave(output_path, plot = plot, width = 8, height = 6)
  }
  invisible(plot)
}

## Additional splicing utilities can be added below -----------------------------
# For example, functions to convert events to gene sets, to summarise GO terms
# across directions or to map event identifiers to gene symbols may be added
# here as needed.
