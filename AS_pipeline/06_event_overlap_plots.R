# 06_event_overlap_plots.R
#
# Step 6: Compute and visualise overlaps of splicing events across conditions.
#
# For each event type (EX, INT, ALTA, ALTD) and direction (all/pos/neg) this
# script builds a list of event identifiers per condition and constructs a
# binary membership matrix.  UpSet plots are generated using the helper
# function `plot_event_upset()` (which wraps around `plot_upset_complex()`
# from shared/99_utils.R).  The membership matrices are saved as RDS files
# for downstream steps.

message("Step 06: Building event overlap plots")

## Check prerequisites ---------------------------------------------------------
if (!exists("EVENT_SPLITS", envir = .GlobalEnv)) {
  stop("EVENT_SPLITS object is missing; please run step 03 first")
}
event_splits <- get("EVENT_SPLITS", envir = .GlobalEnv)

## Create output directory -----------------------------------------------------
step_output <- file.path(PROJECT_ROOT, "output", "06_event_overlap")
if (!dir.exists(step_output)) dir.create(step_output, recursive = TRUE)

## Build membership matrices and plots ----------------------------------------
conditions <- names(event_splits)
event_types <- unique(unlist(lapply(event_splits, names)))
directions <- c("all", "pos", "neg")
membership_matrices <- list()

for (etype in event_types) {
  membership_matrices[[etype]] <- list()
  for (dirn in directions) {
    # Build list of events per condition
    event_sets <- list()
    for (cond in conditions) {
      sub <- event_splits[[cond]][[etype]][[dirn]]
      if (!is.null(sub)) {
        event_sets[[cond]] <- sub$EVENT
      } else {
        event_sets[[cond]] <- character(0)
      }
    }
    # Build membership matrix
    mat <- build_event_membership_matrix(event_sets)
    membership_matrices[[etype]][[dirn]] <- mat
    # Plot UpSet
    plot_title <- sprintf("Splicing event overlaps (%s, %s)", etype, dirn)
    plot_path <- file.path(step_output, sprintf("UpSet_%s_%s.png", etype, dirn))
    plot_event_upset(mat, title = plot_title, output_path = plot_path, colors = cfg$colors)
  }
}

## Save membership matrices ----------------------------------------------------
membership_rds <- file.path(step_output, "event_membership_matrices.rds")
saveRDS(membership_matrices, membership_rds)
message(sprintf("Event membership matrices saved to %s", membership_rds))

## Return membership to global environment ------------------------------------
assign("EVENT_MEMBERSHIP", membership_matrices, envir = .GlobalEnv)
message("Step 06 completed")