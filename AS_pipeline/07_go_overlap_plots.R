# 07_go_overlap_plots.R
#
# Step 7: Analyse overlaps of GO terms across conditions and directions.
#
# This script reads the GO enrichment results from step 05 and for each
# event type and direction builds sets of GO term descriptions per condition.
# It then constructs binary membership matrices and produces UpSet plots to
# visualise shared and unique GO terms between conditions.  Membership
# matrices are saved as RDS files for further analysis.

message("Step 07: Building GO term overlap plots")

## Check prerequisites ---------------------------------------------------------
if (!exists("GO_RESULTS", envir = .GlobalEnv)) {
  stop("GO_RESULTS object is missing; please run step 05 first")
}
go_results <- get("GO_RESULTS", envir = .GlobalEnv)

## Create output directory -----------------------------------------------------
step_output <- file.path(PROJECT_ROOT, "output", "07_go_overlap")
if (!dir.exists(step_output)) dir.create(step_output, recursive = TRUE)

## Determine event types and directions ---------------------------------------
conditions <- names(go_results)
event_types <- unique(unlist(lapply(go_results, names)))
directions <- c("all", "pos", "neg")
go_membership <- list()

for (etype in event_types) {
  go_membership[[etype]] <- list()
  for (dirn in directions) {
    # Build GO term sets for each condition
    term_sets <- list()
    for (cond in conditions) {
      go_res <- go_results[[cond]][[etype]][[dirn]]
      if (!is.null(go_res) && nrow(go_res) > 0) {
        term_sets[[cond]] <- go_res$Description
      } else {
        term_sets[[cond]] <- character(0)
      }
    }
    # Skip if all term sets are empty
    if (all(lengths(term_sets) == 0)) {
      next
    }
    # Build binary membership matrix
    mat <- build_event_membership_matrix(term_sets)
    go_membership[[etype]][[dirn]] <- mat
    # Plot UpSet for GO terms
    plot_title <- sprintf("GO term overlaps (%s, %s)", etype, dirn)
    plot_path <- file.path(step_output, sprintf("GO_UpSet_%s_%s.png", etype, dirn))
    # Use wrapper directly for terms; pass colours if available
    plot_event_upset(mat, title = plot_title, output_path = plot_path, colors = cfg$colors)
  }
}

## Save GO membership matrices -------------------------------------------------
go_membership_rds <- file.path(step_output, "go_term_membership_matrices.rds")
saveRDS(go_membership, go_membership_rds)
message(sprintf("GO term membership matrices saved to %s", go_membership_rds))

## Return membership to global environment ------------------------------------
assign("GO_MEMBERSHIP", go_membership, envir = .GlobalEnv)
message("Step 07 completed")