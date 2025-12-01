# 02_summary_metrics.R
#
# Step 2: Compute general summary metrics for filtered splicing events.
#
# This script consumes the list of filtered event tables produced in step 01
# (available as the global object `FILTERED_EVENTS`) and computes basic
# statistics for each condition using `compute_splicing_metrics()`.  It also
# produces simple histograms of ΔPSI and MV.dPsi distributions for quick
# inspection.  Results are saved to the output directory `02_metrics/`.

message("Step 02: Computing summary metrics")

## Create output directory -----------------------------------------------------
step_output <- file.path(PROJECT_ROOT, "output", "02_metrics")
if (!dir.exists(step_output)) dir.create(step_output, recursive = TRUE)

# Check that filtered events are available
if (!exists("FILTERED_EVENTS", envir = .GlobalEnv)) {
  stop("FILTERED_EVENTS object is missing; please run step 01 first")
}
filtered_tables <- get("FILTERED_EVENTS", envir = .GlobalEnv)

## Compute metrics per condition ----------------------------------------------
metrics_list <- purrr::map(filtered_tables, compute_splicing_metrics)
metrics_df <- dplyr::bind_rows(metrics_list, .id = "condition")

## Save metrics to disk --------------------------------------------------------
metrics_tsv <- file.path(step_output, "summary_metrics.tsv")
readr::write_tsv(metrics_df, metrics_tsv)
metrics_rds <- file.path(step_output, "summary_metrics.rds")
saveRDS(metrics_df, metrics_rds)
message(sprintf("Summary metrics saved to %s and %s", metrics_tsv, metrics_rds))

## Plot distributions ----------------------------------------------------------
hist_output_mv <- file.path(step_output, "mv_distribution.png")
hist_output_dpsi <- file.path(step_output, "dpsi_distribution.png")

# Collect MV and dPSI values across all conditions
mv_values <- c()
dpsi_values <- c()
for (cond in names(filtered_tables)) {
  df <- filtered_tables[[cond]]
  mv_col <- intersect(c("MV.dPsi._at_0.95", "MV"), names(df))
  dpsi_col <- intersect(c("E.dPsi.", "dPsi", "DeltaPSI"), names(df))
  if (length(mv_col) > 0) mv_values <- c(mv_values, df[[mv_col[1]]])
  if (length(dpsi_col) > 0) dpsi_values <- c(dpsi_values, df[[dpsi_col[1]]])
}

library(ggplot2)

if (length(mv_values) > 0) {
  p_mv <- ggplot(data.frame(MV = mv_values), aes(x = MV)) +
    geom_histogram(binwidth = 0.02, fill = "steelblue", colour = "white") +
    labs(title = "Distribution of MV.dPsi", x = "MV.dPsi", y = "Count") +
    theme_minimal()
  ggsave(hist_output_mv, p_mv, width = 6, height = 4)
  message(sprintf("MV distribution plot saved to %s", hist_output_mv))
} else {
  warning("No MV values found; MV distribution plot not generated")
}

if (length(dpsi_values) > 0) {
  p_dpsi <- ggplot(data.frame(dPSI = dpsi_values), aes(x = dPSI)) +
    geom_histogram(binwidth = 0.01, fill = "firebrick", colour = "white") +
    labs(title = "Distribution of ΔPSI", x = "ΔPSI", y = "Count") +
    theme_minimal()
  ggsave(hist_output_dpsi, p_dpsi, width = 6, height = 4)
  message(sprintf("ΔPSI distribution plot saved to %s", hist_output_dpsi))
} else {
  warning("No ΔPSI values found; ΔPSI distribution plot not generated")
}

## Return metrics to global environment ---------------------------------------
assign("SUMMARY_METRICS", metrics_df, envir = .GlobalEnv)
message("Step 02 completed")