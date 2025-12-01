# New helper functions for shared/99_utils.R
#
# The following functions are general utilities that can be added to
# `shared/99_utils.R`.  They provide reusable functionality for constructing
# membership matrices from sets and for plotting UpSet diagrams.  No existing
# functions in `99_utils.R` are modified; these definitions can simply be
# appended to the end of that file.

## build_event_membership_matrix ---------------------------------------------
#' Build a binary membership matrix from a list of sets.
#'
#' Given a named list of vectors (e.g. splicing event IDs or GO terms),
#' constructs a data frame where rows correspond to the union of all unique
#' identifiers and columns correspond to the list names.  Entries are 1 if
#' the identifier is present in the set and 0 otherwise.
#'
#' @param sets Named list of character vectors.
#' @return A data frame with an identifier column (named `id`) and one
#' column per set containing 0/1 indicators.
build_event_membership_matrix <- function(sets) {
  stopifnot(is.list(sets), length(sets) > 0, !is.null(names(sets)))
  ids <- unique(unlist(sets))
  membership <- lapply(sets, function(x) as.integer(ids %in% x))
  mat <- as.data.frame(membership)
  mat <- cbind(id = ids, mat)
  return(mat)
}

## plot_event_upset ------------------------------------------------------------
#' Plot an UpSet diagram for a binary membership matrix.
#'
#' This convenience wrapper around `plot_upset_complex()` (defined in
#' `shared/99_utils.R`) accepts a membership matrix created by
#' `build_event_membership_matrix()` and produces an UpSet diagram.  The first
#' column of the matrix is treated as the identifier column and is removed
#' before plotting.  Colours can be supplied to customise the set bars.
#'
#' @param membership_matrix Data frame with first column containing IDs and
#' subsequent columns containing 0/1 membership indicators.
#' @param title Character string for the plot title.
#' @param output_path Optional file name for saving the plot (PNG).
#' @param colors Optional named vector of colours for the sets.
#' @return Invisibly returns the ggplot object produced by
#' `plot_upset_complex()`.
plot_event_upset <- function(membership_matrix,
                             title = "UpSet diagram",
                             output_path = NULL,
                             colors = NULL) {
  stopifnot(is.data.frame(membership_matrix), ncol(membership_matrix) >= 2)
  mat <- membership_matrix[, -1, drop = FALSE]
  if (!is.null(colors)) {
    p <- plot_upset_complex(mat, title = title, set_colors = colors)
  } else {
    p <- plot_upset_complex(mat, title = title)
  }
  if (!is.null(output_path)) {
    ggplot2::ggsave(filename = output_path, plot = p, width = 8, height = 6)
  }
  invisible(p)
}
