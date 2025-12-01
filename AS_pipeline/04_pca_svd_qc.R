# 04_pca_svd_qc.R
#
# Step 4: PCA and SVD quality control.
#
# This script performs principal component analysis (PCA) and singular
# value decomposition (SVD) on the filtered PSI matrix produced in step 01.
# PCA and SVD are run separately for each major event type (EX, INT, ALTA,
# ALTD).  If fewer than two events or two samples remain for a given type
# the analysis is skipped for that type.  Results include PC1/PC2 scatter
# plots, SVD scree plots and, optionally, SVD on the residual matrix with
# the first component removed.  All results are written to the output
# directory `04_pca/`.

message("Step 04: Performing PCA and SVD QC")

## Check if PCA/SVD is enabled -----------------------------------------------
if (!isTRUE(cfg$pca$enabled)) {
  message("PCA is disabled in configuration; skipping step 04")
  return(invisible(NULL))
}

## Create output directory -----------------------------------------------------
step_output <- file.path(PROJECT_ROOT, "output", "04_pca")
if (!dir.exists(step_output)) dir.create(step_output, recursive = TRUE)

## Retrieve PSI matrix and event splits ---------------------------------------
if (!exists("PSI_MATRIX", envir = .GlobalEnv)) {
  stop("PSI_MATRIX object is missing; please run step 01 first")
}
psi_mat <- get("PSI_MATRIX", envir = .GlobalEnv)

## Determine sample information ------------------------------------------------
sample_info <- data.frame(
  sample = unlist(cfg$sets),
  condition = rep(names(cfg$sets), lengths(cfg$sets)),
  stringsAsFactors = FALSE
)
rownames(sample_info) <- sample_info$sample

## Helper to compute PCA -------------------------------------------------------
compute_pca <- function(mat) {
  # Replace NA values with row medians
  mat_clean <- mat
  na_idx <- which(is.na(mat_clean), arr.ind = TRUE)
  if (length(na_idx) > 0) {
    for (i in seq_len(nrow(mat_clean))) {
      row_vals <- mat_clean[i, ]
      row_med <- median(row_vals, na.rm = TRUE)
      row_vals[is.na(row_vals)] <- row_med
      mat_clean[i, ] <- row_vals
    }
  }
  # Transpose so that samples are rows and events are columns
  mat_t <- t(mat_clean)
  pca_res <- prcomp(mat_t, center = TRUE, scale. = TRUE)
  return(pca_res)
}

## Helper to compute SVD -------------------------------------------------------
compute_svd <- function(mat, subtract_pc1 = FALSE, max_components = 6) {
  # Replace NAs as for PCA
  mat_clean <- mat
  if (any(is.na(mat_clean))) {
    for (i in seq_len(nrow(mat_clean))) {
      row_vals <- mat_clean[i, ]
      row_med <- median(row_vals, na.rm = TRUE)
      row_vals[is.na(row_vals)] <- row_med
      mat_clean[i, ] <- row_vals
    }
  }
  # Transpose: samples as rows
  mat_t <- t(mat_clean)
  # Optionally subtract PC1 component
  if (isTRUE(subtract_pc1)) {
    # Compute PCA to get PC1 loadings
    pca_res <- prcomp(mat_t, center = TRUE, scale. = TRUE)
    pc1_scores <- pca_res$x[, 1]
    pc1_loadings <- pca_res$rotation[, 1]
    # Subtract PC1 contribution
    mat_resid <- mat_t - outer(pc1_scores, pc1_loadings)
    mat_centered <- scale(mat_resid, center = TRUE, scale = TRUE)
    svd_res <- svd(mat_centered)
  } else {
    mat_centered <- scale(mat_t, center = TRUE, scale = TRUE)
    svd_res <- svd(mat_centered)
  }
  # Compute variance explained
  svals <- svd_res$d
  var_explained <- svals^2 / sum(svals^2)
  # Keep up to max_components
  comps <- seq_len(min(length(svals), max_components))
  result <- list(
    d = svals[comps],
    u = svd_res$u[, comps, drop = FALSE],
    v = svd_res$v[, comps, drop = FALSE],
    var_explained = var_explained[comps]
  )
  return(result)
}

## Loop over event types -------------------------------------------------------
event_types <- c("EX", "INT", "ALTA", "ALTD")
for (etype in event_types) {
  message(sprintf("Processing event type: %s", etype))
  # Select events belonging to this type based on EVENT prefix
  event_idx <- grepl(etype, rownames(psi_mat))
  mat <- psi_mat[event_idx, , drop = FALSE]
  # Skip if not enough data
  if (nrow(mat) < 2 || ncol(mat) < 2) {
    warning(sprintf("Skipping PCA/SVD for %s: too few events or samples", etype))
    next
  }
  # Compute PCA
  pca_res <- compute_pca(mat)
  # PC1 vs PC2 plot
  scores <- data.frame(pca_res$x[, 1:2])
  scores$sample <- rownames(pca_res$x)
  scores$condition <- sample_info[scores$sample, "condition"]
  p <- ggplot(scores, aes(x = PC1, y = PC2, colour = condition)) +
    geom_point(size = 3) +
    labs(title = sprintf("PCA of PSI values (%s)", etype), x = "PC1", y = "PC2") +
    scale_colour_manual(values = cfg$colors) +
    theme_minimal()
  pca_plot_path <- file.path(step_output, sprintf("PCA_%s_PC1_PC2.png", etype))
  ggsave(pca_plot_path, p, width = 6, height = 4)
  # Save PCA object
  pca_rds_path <- file.path(step_output, sprintf("PCA_%s.rds", etype))
  saveRDS(pca_res, pca_rds_path)
  # Compute SVD if enabled
  if (isTRUE(cfg$pca$svd$enabled)) {
    svd_res <- compute_svd(mat,
                           subtract_pc1 = isTRUE(cfg$pca$svd$subtract_pc1),
                           max_components = cfg$pca$svd$max_components)
    # Scree plot
    svd_df <- data.frame(Component = seq_along(svd_res$var_explained),
                         VarianceExplained = svd_res$var_explained)
    svd_plot <- ggplot(svd_df, aes(x = Component, y = VarianceExplained)) +
      geom_col(fill = "darkorange") +
      geom_line(group = 1) +
      geom_point() +
      labs(title = sprintf("SVD Scree (%s)", etype), x = "Component", y = "Variance Explained") +
      theme_minimal()
    svd_plot_path <- file.path(step_output, sprintf("SVD_scree_%s.png", etype))
    ggsave(svd_plot_path, svd_plot, width = 6, height = 4)
    # Save SVD object
    svd_rds_path <- file.path(step_output, sprintf("SVD_%s.rds", etype))
    saveRDS(svd_res, svd_rds_path)
  }
}

message("Step 04 completed")