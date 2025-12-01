# 01_read_filter_splicing.R
#
# Step 1: Read and filter VAST‑TOOLS splicing tables.
#
# This script reads splicing tables for each condition defined in the
# configuration, applies MV and ΔPSI filtering, optionally restricts to
# events in reliably expressed genes and retains events with sufficient
# replicate calls.  It then writes the filtered tables and a combined PSI
# matrix to the output directory `01_filter/` under the project root.

message("Step 01: Reading and filtering splicing tables")

## Create output directory for this step --------------------------------------
step_output <- file.path(PROJECT_ROOT, "output", "01_filter")
if (!dir.exists(step_output)) dir.create(step_output, recursive = TRUE)

## Read VAST‑TOOLS tables ------------------------------------------------------
vast_dir <- cfg$vast_tools_dir
sample_sets <- cfg$sets
tables <- read_vast_tools_tables(vast_dir, sample_sets)

## Load expressed gene list (optional) ----------------------------------------
expressed_genes <- NULL
if (isTRUE(cfg$splicing_filter$expressed_only)) {
  expr_file <- cfg$expression_counts
  if (!is.null(expr_file) && file.exists(expr_file)) {
    message(sprintf("Loading expressed genes from %s", expr_file))
    expr_obj <- readRDS(expr_file)
    if (is.matrix(expr_obj) || is.data.frame(expr_obj)) {
      expressed_genes <- rownames(expr_obj)
    } else if (is.list(expr_obj) && !is.null(names(expr_obj))) {
      expressed_genes <- names(expr_obj)
    } else {
      stop("expression_counts file must be a matrix/data.frame with row names or a named object")
    }
  } else {
    warning("expressed_only is TRUE but expression_counts file is missing; skipping expression filter")
  }
}

## Apply filtering to each condition ------------------------------------------
filtered_tables <- list()
for (cond in names(tables)) {
  df <- tables[[cond]]
  message(sprintf("Filtering splicing events for %s", cond))
  filtered <- filter_splicing_events(
    df,
    mv_threshold = cfg$splicing_filter$mv_threshold,
    dpsi_threshold = cfg$splicing_filter$dpsi_threshold,
    expressed_only = isTRUE(cfg$splicing_filter$expressed_only),
    expressed_genes = expressed_genes,
    min_valid_calls = cfg$splicing_filter$min_valid_calls
  )
  filtered_tables[[cond]] <- filtered
  # Save individual condition table as TSV
  out_tsv <- file.path(step_output, sprintf("%s_filtered_events.tsv", cond))
  readr::write_tsv(filtered, out_tsv)
  message(sprintf("Written filtered events for %s to %s", cond, out_tsv))
}

## Save combined filtered tables as RDS ---------------------------------------
filtered_rds_path <- file.path(step_output, "filtered_splicing_tables.rds")
saveRDS(filtered_tables, filtered_rds_path)
message(sprintf("Saved filtered tables to %s", filtered_rds_path))

## Build combined PSI matrix ---------------------------------------------------
# Extract PSI columns for each condition and sample.  We attempt to match
# columns to sample IDs defined in cfg$sets.  If multiple columns match a
# sample ID the first match is used.  Columns are renamed to the sample ID
# for consistency.  Missing samples result in all NA values.
psi_mats <- list()
for (cond in names(filtered_tables)) {
  df <- filtered_tables[[cond]]
  cond_samples <- sample_sets[[cond]]
  # Determine candidate PSI columns (matching sample IDs)
  selected_cols <- list()
  for (sid in cond_samples) {
    matches <- grep(paste0("^", sid), names(df), value = TRUE)
    if (length(matches) == 0) {
      # Try to find columns ending with sample ID
      matches <- grep(paste0(sid, "$") , names(df), value = TRUE)
    }
    if (length(matches) > 0) {
      selected_cols[[sid]] <- matches[1]
    } else {
      warning(sprintf("No PSI column found for sample %s in condition %s", sid, cond))
    }
  }
  # Construct a data frame with EVENT and PSI values for matched samples
  sub_df <- df %>%
    select(EVENT, all_of(unlist(selected_cols)))
  # Rename columns to sample IDs
  colnames(sub_df) <- c("EVENT", names(selected_cols))
  psi_mats[[cond]] <- sub_df
}

# Merge PSI matrices by EVENT
psi_join <- purrr::reduce(psi_mats, function(x, y) {
  full_join(x, y, by = "EVENT")
})

# Create matrix with EVENTS as row names and samples as columns
psi_matrix <- psi_join %>% arrange(EVENT)
event_ids <- psi_matrix$EVENT
psi_matrix <- psi_matrix %>% select(-EVENT)
psi_mat <- as.matrix(psi_matrix)
rownames(psi_mat) <- event_ids

# Save PSI matrix
psi_rds_path <- file.path(step_output, "filtered_psi_matrix.rds")
saveRDS(psi_mat, psi_rds_path)
message(sprintf("Saved filtered PSI matrix to %s", psi_rds_path))

## Return objects to global environment for subsequent steps --------------------
assign("FILTERED_EVENTS", filtered_tables, envir = .GlobalEnv)
assign("PSI_MATRIX", psi_mat, envir = .GlobalEnv)
message("Step 01 completed")