## scripts/03_pca_filtered_counts.R
## PCA on bimodal-filtered gene-level counts + PC1 vs PC2 plot.

if (!exists("cfg"))        stop("'cfg' not found (configuration list).")
if (!exists("output_dir")) stop("'output_dir' not found.")

out_dir_pca <- ensure_dir(output_dir, "03_pca")

## ----------------------------------------------------------------------
## 1) Check whether PCA step is enabled
## ----------------------------------------------------------------------

if (isFALSE(cfg$pca$enabled)) {
  message("PCA on filtered counts is disabled in cfg$pca$enabled.")
} else {

  ## --------------------------------------------------------------------
  ## 2) Load filtered counts matrix
  ## --------------------------------------------------------------------
  
  if (cfg$pca$counts_source == "kallisto_bimodal_filtered") {
    message("PCA will be performed on filtered counts from DESeq2 + bimodal filter.")
    title = "PCA_on_Bimodal-filtered_Kallisto_counts"
    counts_path <- file.path(output_dir, "02_deseq", "filtered_counts_matrix.rds")
    if (!file.exists(counts_path)) {
      stop("Filtered counts file not found: ", counts_path)
    }
    counts <- readRDS(counts_path)
  } else if (cfg$pca$counts_source == "kallisto_raw") {
    message("PCA will be performed on raw Kallisto counts.")
    title = "PCA_on_Raw_Kallisto_counts"
    counts <- as.matrix(COUNT_DATA_GENE)
    
  } else {
    stop("Unsupported cfg$pca$counts_source: ", cfg$pca$counts_source)
  }

  ## Convert raw counts to CPM
  lib_sizes <- colSums(counts)
  if (any(lib_sizes <= 0)) {
    stop("Some samples have non-positive library sizes; cannot compute CPM.")
  }

  counts <- t(t(counts) / lib_sizes * 1e6)

  ## --------------------------------------------------------------------
  ## 3) Obtain sample information
  ## --------------------------------------------------------------------
  ## Align sample_info to columns of counts
  
  if (!all(colnames(counts) %in% rownames(sample_info))) {
    stop("All count matrix columns must be present in rownames(sample_info).")
  }
  sample_info <- sample_info[colnames(counts), , drop = FALSE]

  ## Align sample_info to columns of counts
  if (!all(colnames(counts) %in% rownames(sample_info))) {
    stop("All count matrix columns must be present in rownames(sample_info).")
  }
  sample_info <- sample_info[colnames(counts), , drop = FALSE]

  ## Number of top variable features (if used by pca_from_counts)
  n_top_var <- cfg$pca$n_top_var
  if (is.null(n_top_var) || !is.finite(n_top_var) || n_top_var <= 0) {
    n_top_var <- length(rownames(counts))
  }

  group_colors <- unlist(COLOR_KO)

  # Keep a copy of raw counts: SVD uses counts (log2 + quantile), not CPM
  counts_for_svd <- counts
  counts_for_clust <- counts

  ## --------------------------------------------------------------------
  ## 4) Run PCA
  ## --------------------------------------------------------------------
  pca_res <- pca_from_counts(
    counts        = counts,
    sample_info   = sample_info,
    log_transform = TRUE,
    pseudocount   = 1,
    n_top_var     = n_top_var
  )

  saveRDS(pca_res, file.path(out_dir_pca, "pca_filtered_counts.rds"))

  ## --------------------------------------------------------------------
  ## 5) Plot PC1 vs PC2
  ## --------------------------------------------------------------------
  p <- plot_pca_scores(
    pca_result = pca_res,
    title      = title,
    group_colors = group_colors
  )

  ggplot2::ggsave(
    filename = file.path(out_dir_pca, paste0(title,".png")),
    plot     = p,
    width    = 6,
    height   = 5,
    dpi      = 300
  )

    ## --------------------------------------------------------------------
  ## 3) Hierarchical clustering of samples (replicates) from expression
  ## --------------------------------------------------------------------
  
  hc_cfg <- cfg$pca$hclust
  
  if (!is.null(hc_cfg) && isTRUE(hc_cfg$enabled)) {
    message("Sample hierarchical clustering from expression is enabled (cfg$pca$hclust).")
    
    # Use specific n_top_var for clustering if provided;
    # otherwise, reuse the same n_top_var used for PCA.
    n_top_var_hc <- hc_cfg$n_top_var
    if (is.null(n_top_var_hc)) {
      n_top_var_hc <- n_top_var
    }
    
    dist_type     <- if (is.null(hc_cfg$dist_type))     "pearson"   else hc_cfg$dist_type
    dist_method   <- if (is.null(hc_cfg$dist_method))   "euclidean" else hc_cfg$dist_method
    hclust_method <- if (is.null(hc_cfg$hclust_method)) "complete"  else hc_cfg$hclust_method
    
    # Build clustering directly from expression matrix (pre-PCA),
    # using the *aligned* sample_info, not the global SAMPLE_INFO.
    hc <- sample_hclust_from_matrix(
      counts        = counts_for_clust,
      sample_info   = sample_info,
      log_transform = TRUE,
      pseudocount   = 1,
      n_top_var     = n_top_var_hc,
      dist_type     = dist_type,
      dist_method   = dist_method,
      hclust_method = hclust_method
    )
    
    ## ------------------------------------------------------------------
    ## Derive labels (condition + replicate index) in the same order
    ## as hc$order, and pass them explicitly to plot().
    ## ------------------------------------------------------------------
    
    # Sample IDs in the hclust object (these are column names of counts)
    samples <- hc$labels
    
    # Use aligned sample_info (same object as used above)
    cond <- as.character(sample_info[samples, "condition"])
    
    # Extract replicate index from sample name (e.g. "SHSY5Y_WT_1" -> "1")
    rep_lab <- sub(".*_(\\d+)$", "\\1", samples)
    no_match <- rep_lab == samples
    rep_lab[no_match] <- samples  # fall back to full sample name
    
    labs <- paste(cond, rep_lab, sep = "_")
    
    # Reorder labels to match plotting order (hc$order)
    plot_labs <- labs[hc$order]
    
    hc_file <- file.path(out_dir_pca, "Sample_hclust_from_expression.png")
    
    grDevices::png(hc_file, width = 800, height = 600, res = 120)
    plot(
      hc,
      labels = plot_labs,
      main   = "Hierarchical clustering of samples\n(log2 expression, Pearson)",
      xlab   = "",
      sub    = ""
    )
    grDevices::dev.off()
    
    message("Saved sample clustering dendrogram to: ", hc_file)
  } else {
    message("Sample hierarchical clustering is disabled (cfg$pca$hclust$enabled is FALSE or missing).")
  }




  ## --------------------------------------------------------------------
  ## SVD on the same counts
  ## --------------------------------------------------------------------
  
  svd_cfg <- cfg$pca$svd
  
  if (!is.null(svd_cfg) && isTRUE(svd_cfg$enabled)) {
    message("Running SVD on filtered counts.")
    
    subtract_pc1 <- if (is.null(svd_cfg$subtract_pc1)) {
      TRUE
    } else {
      isTRUE(svd_cfg$subtract_pc1)
    }
    
    # SVD on counts_for_svd using SAMPLE_INFO for conditions
    svd_res <- run_svd(
      counts            = counts_for_svd,
      sample_info       = SAMPLE_INFO,
      log_pseudocount   = 1,
      quantile_normalise = TRUE,
      subtract_pc1      = subtract_pc1
    )
    
    # Save full SVD results
    saveRDS(
      svd_res,
      file.path(out_dir_pca, "svd_results.rds")
    )
    message("Saved SVD results to 03_pca/svd_results.rds")
    
    # Scree plot (after removing component 1)
    svd_scree <- plot_svd_scree(
      svd_res$variance_residual,
      title = "SVD (after removing component 1)"
    )
    
    ggplot2::ggsave(
      filename = file.path(out_dir_pca, "SVD_scree.png"),
      plot     = svd_scree,
      width    = 6,
      height   = 4,
      dpi      = 300
    )
    
    # Mean ± SD of component scores by condition
    max_comp <- svd_cfg$max_components
    if (is.null(max_comp) || !is.finite(max_comp) || max_comp <= 0L) {
      max_comp <- length(svd_res$svd_residual$d)
    }
    
    # Number of available components from residual SVD
    n_comp <- ncol(svd_res$svd_residual$v)
    max_comp <- min(max_comp, n_comp)
    
    for (j in seq_len(max_comp)) {
      p_cmp <- plot_svd_component_means(
        svd_result = svd_res,
        component  = j,
        palette    = group_colors
      )
      
      ggplot2::ggsave(
        filename = file.path(
          out_dir_pca,
          sprintf("SVD_PC%02d_mean_by_condition.png", j)
        ),
        plot     = p_cmp,
        width    = 5,
        height   = 4,
        dpi      = 300
      )
    }
  } else {
    message("SVD step is disabled (cfg$pca$svd$enabled is FALSE or missing).")
  }

}

# stop("PCA and SVD analysis completed successfully.")