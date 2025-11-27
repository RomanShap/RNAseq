## scripts/99_utils.R
## Generic helpers used by the whole quantification pipeline.

## ----------------------------------------------------------------------
## Directory helpers
## ----------------------------------------------------------------------

ensure_dir <- function(...) {
  path <- file.path(...)
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    message("Created directory: ", normalizePath(path, mustWork = FALSE))
  }
  path
}
## ----------------------------------------------------------------------
## Split title by multiple rows
## ----------------------------------------------------------------------
wrap_title_once <- function(x, width = 45L) {
  # Handle NA or non-character input defensively
  if (is.null(x) || is.na(x)) return(x)
  x <- as.character(x)
  
  # If short enough, return unchanged
  if (nchar(x) <= width) return(x)
  
  # Find all spaces
  space_pos <- gregexpr(" ", x, fixed = TRUE)[[1]]
  # Keep only spaces before or at 'width'
  space_pos <- space_pos[space_pos <= width]
  
  if (length(space_pos) == 0) {
    # No space before 'width': hard split
    return(
      paste0(
        substr(x, 1L, width),
        "\n",
        substr(x, width + 1L, nchar(x))
      )
    )
  }
  
  # Split at the last space before 'width'
  split_pos <- max(space_pos)
  part1 <- substr(x, 1L, split_pos - 1L)
  part2 <- substr(x, split_pos + 1L, nchar(x))
  
  paste0(part1, "\n", part2)
}


## ----------------------------------------------------------------------
## tx2gene and MSigDB caching
## ----------------------------------------------------------------------

load_or_build_tx2gene <- function(tx2gene_rds, transcript_ids,
                                  dataset = "hsapiens_gene_ensembl") {
  transcript_ids <- unique(transcript_ids)
  
  if (!is.null(tx2gene_rds) && file.exists(tx2gene_rds)) {
    message("Loading tx2gene mapping from RDS: ", tx2gene_rds)
    tx2gene <- readRDS(tx2gene_rds)
    return(tx2gene)
  }
  
  message("Building tx2gene mapping via biomaRt...")
  mart <- biomaRt::useMart("ensembl", dataset = dataset)
  
  tx2gene <- biomaRt::getBM(
    attributes = c("ensembl_transcript_id", "ensembl_gene_id"),
    filters    = "ensembl_transcript_id",
    values     = transcript_ids,
    mart       = mart
  )
  
  tx2gene <- tx2gene[!duplicated(tx2gene$ensembl_transcript_id), ]
  
  if (!is.null(tx2gene_rds)) {
    ensure_dir(dirname(tx2gene_rds))
    saveRDS(tx2gene, tx2gene_rds)
    message("Saved tx2gene mapping to: ", tx2gene_rds)
  }
  
  tx2gene
}

load_or_build_msigdb_c5 <- function(msigdb_c5_rds,
                                    species = "Homo sapiens",
                                    collection = "C5") {
  if (!is.null(msigdb_c5_rds) && file.exists(msigdb_c5_rds)) {
    message("Loading MSigDB C5 from RDS: ", msigdb_c5_rds)
    return(readRDS(msigdb_c5_rds))
  }
  
  message("Downloading MSigDB C5 gene sets via msigdbr...")
  msig <- msigdbr::msigdbr(species = species, collection = collection)
  ensure_dir(dirname(msigdb_c5_rds))
  saveRDS(msig, msigdb_c5_rds)
  msig
}

## ----------------------------------------------------------------------
## ID / annotation helpers
## ----------------------------------------------------------------------

add_gene_symbols <- function(res_df,
                             keycol = "ensembl_id",
                             orgdb = org.Hs.eg.db,
                             new_col = "symbol") {
  if (!keycol %in% names(res_df)) {
    stop("Column '", keycol, "' not found in result data frame.")
  }
  keys <- res_df[[keycol]]
  symbols <- AnnotationDbi::mapIds(
    orgdb,
    keys      = keys,
    column    = "SYMBOL",
    keytype   = "ENSEMBL",
    multiVals = "first"
  )
  res_df[[new_col]] <- unname(symbols)
  res_df
}

gene_label_from_res <- function(res_df,
                                id_col = "ensembl_id",
                                sym_col = "symbol") {
  id  <- res_df[[id_col]]
  sym <- if (sym_col %in% names(res_df)) res_df[[sym_col]] else NA_character_
  ifelse(is.na(sym) | sym == "", id, sym)
}

## ----------------------------------------------------------------------
## Bimodal expression filter
## ----------------------------------------------------------------------

bimodal_filter_counts <- function(counts,
                                  groups,
                                  max_sd   = 5,
                                  plot     = FALSE,
                                  plot_dir = NULL) {
  # Basic checks ------------------------------------------------------------
  if (!is.matrix(counts) && !is.data.frame(counts)) {
    stop("'counts' must be a matrix or data.frame.")
  }
  counts <- as.matrix(counts)
  
  if (length(groups) != ncol(counts)) {
    stop("Length of 'groups' must match number of columns in 'counts'.")
  }
  groups <- factor(groups)
  
  # Log-transform counts ----------------------------------------------------
  counts_log <- log2(1 + counts)
  
  # Fit 2-component Gaussian mixtures per sample ----------------------------
  bimdens <- lapply(seq_len(ncol(counts_log)), function(ix) {
    mclust::densityMclust(data = counts_log[, ix], G = 2, plot = FALSE)
  })
  
  # Derive per-sample thresholds (99% quantile of lower-mean component) -----
  lims <- vapply(bimdens, function(x) {
    lower_index      <- which.min(x$parameters$mean)
    background_mean  <- x$parameters$mean[lower_index]
    background_sd    <- sqrt(x$parameters$variance$sigmasq[lower_index])
    
    # Clamp extremely large SDs to avoid absurd thresholds
    if (background_sd > max_sd) {
    message("Warning: background SD (", background_sd,
                ") exceeds max_sd (", max_sd, "). Using max_sd.")
      background_sd <- max_sd
    }
    
    stats::qnorm(0.99, mean = background_mean, sd = background_sd)
  }, numeric(1))
  
  # Plot distributions and thresholds -----------------------------
  if (plot) {
    if (!is.null(plot_dir)) {
      # Use your helper if available; otherwise dir.create
      if (exists("ensure_dir", mode = "function")) {
        ensure_dir(plot_dir)
      } else {
        if (!dir.exists(plot_dir)) {
          dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
        }
      }
    }
    
    for (ix in seq_len(ncol(counts_log))) {
      sample_name <- colnames(counts_log)[ix]
      this_log    <- counts_log[, ix]
      thr         <- lims[ix]
      
      if (is.null(sample_name) || sample_name == "") {
        sample_name <- paste0("sample_", ix)
      }
      
      # Decide device: PNG if plot_dir set, otherwise just plot to current device
      if (!is.null(plot_dir)) {
        file_out <- file.path(
          plot_dir,
          paste0("bimodal_", sample_name, ".png")
        )
        png(file_out, width = 1200, height = 900, res = 150)
      }
      
      hist(
        this_log,
        breaks = 100,
        main   = paste0("Bimodal fit: ", sample_name),
        xlab   = "log2(1 + counts)",
        ylab   = "Frequency"
      )
      abline(v = thr, lwd = 2, lty = 2)
      legend(
        "topright",
        legend = sprintf("Threshold = %.2f", thr),
        bty    = "n"
      )
      
      if (!is.null(plot_dir)) {
        dev.off()
      }
    }
  }
  
  # Per-sample expression calls ---------------------------------------------
  is_expr_samples <- sapply(seq_len(ncol(counts_log)), function(ix) {
    counts_log[, ix] > lims[ix]
  })
  is_expr_samples <- as.matrix(is_expr_samples)
  colnames(is_expr_samples) <- colnames(counts)
  rownames(is_expr_samples) <- rownames(counts)
  
  # Per-group expression: must be TRUE in all samples of the group ----------
  is_expr_groups <- t(apply(is_expr_samples, 1L, function(z) {
    tapply(z, INDEX = groups, FUN = function(w) sum(w) == length(w))
  }))
  
  # Global mask: expressed in at least one group ----------------------------
  is_expr_global <- rowSums(is_expr_groups) >= 1
  
  # Return -------------------------------------------------------------------
  list(
    mask       = is_expr_global,
    per_sample = is_expr_samples,
    per_group  = is_expr_groups,
    thresholds = lims
  )
}

## ----------------------------------------------------------------------
## DESeq2 with bimodal filter and multiple contrasts
## ----------------------------------------------------------------------

run_deseq_bimodal <- function(counts,
                              sample_info,
                              ref_level,
                              padj_threshold = 1) {
  if (!"condition" %in% colnames(sample_info)) {
    stop("sample_info must contain a 'condition' column.")
  }
  
  counts <- as.matrix(counts)
  
  if (!all(colnames(counts) %in% rownames(sample_info))) {
    stop("All count matrix columns must be present in rownames(sample_info).")
  }
  
  sample_info <- sample_info[colnames(counts), , drop = FALSE]
  groups <- factor(sample_info$condition)
  
  filt <- bimodal_filter_counts(
    counts   = counts,       
    groups   = groups,
    plot     = TRUE,
    plot_dir = file.path(output_dir, "02_deseq", "bimodal_plots")
    )

  counts_filt <- counts[filt$mask, , drop = FALSE]
  counts_filt <- round(counts_filt)
  coldata <- data.frame(
    sample_info,
    row.names = rownames(sample_info),
    stringsAsFactors = FALSE
    )
  coldata$condition <- factor(coldata$condition)

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts_filt,
    colData   = coldata,
    design    = ~ condition
  )
  
  dds$condition <- relevel(dds$condition, ref = ref_level)
  dds <- DESeq2::DESeq(dds)
  
  cond_levels <- setdiff(levels(dds$condition), ref_level)
  
  res_list <- vector("list", length(cond_levels))
  names(res_list) <- cond_levels
  
  for (cond in cond_levels) {
    res <- DESeq2::results(dds, contrast = c("condition", cond, ref_level))
    res_df <- as.data.frame(res)
    res_df$ensembl_id <- gsub("\\..*", "", rownames(res_df))
    res_df <- add_gene_symbols(res_df, keycol = "ensembl_id")
    
    res_list[[cond]] <- res_df
  }
  
  list(
    dds      = dds,
    results  = res_list,
    mask     = filt$mask,
    groups   = groups,
    filters  = filt,
    filtered_counts = counts_filt
  )
}

## ============================================================================
## PCA HELPERS
## ============================================================================

# Prepare a matrix for PCA: log-transform, remove zero-variance genes,
# and optionally keep only the most variable genes.
prepare_pca_matrix <- function(counts,
                               log_transform = TRUE,
                               pseudocount = 1,
                               n_top_var = NULL) {
  counts <- as.matrix(counts)
  storage.mode(counts) <- "double"
  
  if (log_transform) {
    mat <- log2(counts + pseudocount)
  } else {
    mat <- counts
  }
  
  # Remove genes with zero variance (cannot contribute to PCA)
  var_vec <- apply(mat, 1L, var, na.rm = TRUE)
  keep <- var_vec > 0
  mat <- mat[keep, , drop = FALSE]
  var_vec <- var_vec[keep]
  
  # Optionally keep only the most variable genes
  if (!is.null(n_top_var) && is.finite(n_top_var) && n_top_var < nrow(mat)) {
    o <- order(var_vec, decreasing = TRUE)
    sel <- o[seq_len(n_top_var)]
    mat <- mat[sel, , drop = FALSE]
  }
  
  if (ncol(mat) < 2L || nrow(mat) < 2L) {
    stop("Not enough samples or features for PCA after preprocessing.")
  }
  
  mat
}

# Compute PCA on a count matrix using sample_info to define groups.
pca_from_counts <- function(counts,
                            sample_info,
                            log_transform = TRUE,
                            pseudocount = 1,
                            n_top_var = 2000) {
  counts <- as.matrix(counts)
  
  # Align sample_info rows to count matrix columns
  if (!all(colnames(counts) %in% rownames(sample_info))) {
    stop("All columns in 'counts' must be present in rownames(sample_info).")
  }
  sample_info <- sample_info[colnames(counts), , drop = FALSE]
  
  if (!"condition" %in% colnames(sample_info)) {
    stop("Column 'condition' not found in sample_info.")
  }
  
  mat <- prepare_pca_matrix(
    counts       = counts,
    log_transform = log_transform,
    pseudocount   = pseudocount,
    n_top_var     = n_top_var
  )
  
  # PCA expects rows = observations (samples), columns = variables (genes)
  pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  var_expl <- (pca$sdev^2) / sum(pca$sdev^2)
  
  groups <- factor(sample_info$condition)
  
  scores <- data.frame(
    sample = rownames(pca$x),
    PC1    = pca$x[, 1L],
    PC2    = pca$x[, 2L],
    group  = groups,
    stringsAsFactors = FALSE
  )
  
  list(
    pca        = pca,
    scores     = scores,
    var_expl   = var_expl,
    sample_info = sample_info
  )
}

# Simple PC1 vs PC2 plot, with optional replicate labels.
plot_pca_scores <- function(pca_result,
                            title = "PCA on filtered counts",
                            group_colors = NULL,
                            label_replicates = TRUE) {
  scores   <- pca_result$scores
  var_expl <- pca_result$var_expl
  
  # Derive a concise replicate label from sample names.
  # Example: "SHSY5Y_WT_1" -> "1".
  # If the pattern is not found, use the full sample name.
  rep_lab <- sub(".*_(\\d+)$", "\\1", scores$sample)
  # If nothing changed (no underscore/number at the end), fall back to sample name.
  no_match <- rep_lab == scores$sample
  rep_lab[no_match] <- scores$sample
  scores$replicate_label <- rep_lab
  
  pc1_lab <- paste0("PC1 (", scales::percent(var_expl[1L]), ")")
  pc2_lab <- paste0("PC2 (", scales::percent(var_expl[2L]), ")")
  
  p <- ggplot(scores, aes(x = PC1, y = PC2, color = group)) +
    geom_point(size = 3)
  
  if (isTRUE(label_replicates)) {
    p <- p +
      geom_text(
        aes(label = replicate_label),
        vjust = -1.0,
        size  = 3
      )
  }
  
  p <- p +
    labs(
      title = title,
      x     = pc1_lab,
      y     = pc2_lab,
      color = "Group"
    ) +
    theme_bw()
  
  if (!is.null(group_colors)) {
    # Keep only colors for groups that actually appear
    used_groups <- levels(scores$group)
    vals <- group_colors[used_groups]
    p <- p + scale_color_manual(values = vals)
  }
  
  p
}

## ============================================================================
## Sample hierarchical clustering from expression matrix (pre-PCA)
##  - counts:       gene x sample matrix (filtered counts)
##  - sample_info:  data.frame with rownames = sample IDs, column 'condition'
##  - uses prepare_pca_matrix() to log-transform and (optionally) select
##    top-variable genes, then clusters samples on that matrix.
##  - dist_type:
##      "pearson"   -> distance = 1 - Pearson correlation between samples
##      "euclidean" -> Euclidean distance on log2 expression
## ============================================================================

sample_hclust_from_matrix <- function(counts,
                                      sample_info,
                                      log_transform = TRUE,
                                      pseudocount   = 1,
                                      n_top_var     = NULL,
                                      dist_type     = c("pearson", "euclidean"),
                                      dist_method   = "euclidean",
                                      hclust_method = "complete") {
  counts <- as.matrix(counts)
  
  # Align sample_info rows to count matrix columns (samples)
  if (!all(colnames(counts) %in% rownames(sample_info))) {
    stop("All columns in 'counts' must be present in rownames(sample_info).")
  }
  sample_info <- sample_info[colnames(counts), , drop = FALSE]
  
  # Preprocess counts: log2, remove zero-variance genes, optional top-var selection
  mat <- prepare_pca_matrix(
    counts        = counts,        # genes x samples
    log_transform = log_transform,
    pseudocount   = pseudocount,
    n_top_var     = n_top_var
  )
  # mat: genes x samples
  
  dist_type <- match.arg(dist_type)
  
  if (dist_type == "pearson") {
    # Pearson correlation between SAMPLES (columns of mat)
    cor_mat <- stats::cor(
      mat,
      method = "pearson",
      use    = "pairwise.complete.obs"
    )
    d <- stats::as.dist(1 - cor_mat)
    
    # Labels are sample IDs
    labels <- colnames(mat)
  } else {
    # Euclidean (or other) distance on expression: rows must be samples
    X <- t(mat)  # samples x genes
    d <- stats::dist(X, method = dist_method)
    
    labels <- rownames(X)  # sample IDs
  }
  
  hc <- stats::hclust(d, method = hclust_method)
  hc$labels <- labels
  
  hc
}




## ============================================================================
## Sample hierarchical clustering based on PCA scores
## ============================================================================

# Build an hclust object using the first n_pcs principal components.
# This clusters samples (replicates) in the same space used by PCA.
sample_hclust_from_pca <- function(pca_result,
                                   n_pcs = 5L,
                                   dist_method = "euclidean",
                                   hclust_method = "complete") {
  if (is.null(pca_result$pca) || is.null(pca_result$pca$x)) {
    stop("pca_result$pca$x is required to compute sample clustering.")
  }
  
  X <- pca_result$pca$x
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  
  # Limit the number of PCs to a sensible value
  n_pcs <- as.integer(n_pcs)
  if (!is.finite(n_pcs) || n_pcs <= 0L) {
    n_pcs <- 5L
  }
  n_pcs <- min(n_pcs, ncol(X))
  
  X_sub <- X[, seq_len(n_pcs), drop = FALSE]
  
  d  <- dist(X_sub, method = dist_method)
  hc <- stats::hclust(d, method = hclust_method)
  
  hc
}


## ============================================================================
## Raphaëlle-style SVD helpers
##  - input: gene x sample counts
##  - steps: log2(1 + counts) -> quantile normalisation -> SVD
##           optionally remove PC1 and recompute SVD on residual matrix
##  - sample coordinates: columns of V
##  - gene projections: U %*% diag(D)
## ============================================================================

run_svd <- function(counts,
                              sample_info,
                              log_pseudocount      = 1,
                              quantile_normalise   = TRUE,
                              subtract_pc1         = TRUE) {
  # Basic checks -------------------------------------------------------------
  counts <- as.matrix(counts)
  storage.mode(counts) <- "double"
  
  if (!all(colnames(counts) %in% rownames(sample_info))) {
    stop("All columns of 'counts' must be present in rownames(sample_info).")
  }
  sample_info <- sample_info[colnames(counts), , drop = FALSE]
  
  if (!"condition" %in% colnames(sample_info)) {
    stop("sample_info must contain a 'condition' column.")
  }
  
  if (nrow(counts) < 2L || ncol(counts) < 2L) {
    stop("Not enough genes or samples for SVD.")
  }
  
  # Step 1: log2(1 + counts) -----------------------------------------------
  mat_log <- log2(counts + log_pseudocount)
  
  # Step 2: quantile normalisation across samples ---------------------------
  if (isTRUE(quantile_normalise)) {
    if (!requireNamespace("limma", quietly = TRUE)) {
      stop("Package 'limma' is required for quantile normalisation. ",
           "Please install it (install.packages('limma')).")
    }
    mat_norm <- limma::normalizeQuantiles(mat_log)
    rownames(mat_norm) <- rownames(mat_log)
    colnames(mat_norm) <- colnames(mat_log)
  } else {
    mat_norm <- mat_log
  }
  
  dat1 <- mat_norm
  
  # Step 3: first SVD on full matrix ---------------------------------------
  svd_full <- svd(dat1)
  pi       <- (svd_full$d ^ 2) / sum(svd_full$d ^ 2)
  
  variance_full_df <- data.frame(
    Component          = seq_along(pi),
    ExplainedVariance  = pi,
    CumulativeVariance = cumsum(pi)
  )
  
  # Step 4: optionally remove PC1 and recompute SVD -------------------------
  if (isTRUE(subtract_pc1)) {
    u1 <- svd_full$u[, 1L, drop = FALSE]
    v1 <- svd_full$v[, 1L, drop = FALSE]
    d1 <- svd_full$d[1L]
    
    dat_res <- dat1 - (u1 %*% t(v1)) * d1
    
    svd_res <- svd(dat_res)
    pip     <- (svd_res$d ^ 2) / sum(svd_res$d ^ 2)
    
    variance_res_df <- data.frame(
      Component          = seq_along(pip),
      ExplainedVariance  = pip,
      CumulativeVariance = cumsum(pip)
    )
  } else {
    dat_res          <- dat1
    svd_res          <- svd_full
    variance_res_df  <- variance_full_df
  }
  
  # Step 5: gene projections (U * D) ----------------------------------------
  gene_proj <- svd_res$u %*% diag(svd_res$d)
  rownames(gene_proj) <- rownames(dat1)
  colnames(gene_proj) <- paste0("PC", seq_len(ncol(gene_proj)))
  
  # Step 6: sample coordinates (V) ------------------------------------------
  sample_mat <- svd_res$v
  rownames(sample_mat) <- colnames(dat1)
  colnames(sample_mat) <- paste0("PC", seq_len(ncol(sample_mat)))
  
  sample_df <- data.frame(
    Sample    = rownames(sample_mat),
    sample_mat,
    Condition = sample_info$condition,
    stringsAsFactors = FALSE
  )
  # Keep condition as factor in original order
  sample_df$Condition <- factor(sample_df$Condition,
                                levels = unique(sample_info$condition))
  
  # Step 7: per-condition summary (mean ± SD per PC) ------------------------
  comp_cols <- grep("^PC", names(sample_df), value = TRUE)
  
  sample_summary <- sample_df |>
    dplyr::group_by(Condition) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(comp_cols),
        list(Mean = mean, SD = sd),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    )
  
  list(
    mat_log           = mat_log,
    mat_norm          = mat_norm,
    mat_residual      = dat_res,
    svd_full          = svd_full,
    svd_residual      = svd_res,
    variance_full     = variance_full_df,
    variance_residual = variance_res_df,
    gene_projections  = gene_proj,
    sample_scores     = sample_df,
    sample_summary    = sample_summary
  )
}

# Scree plot for SVD variance (works on the data.frame returned above)
plot_svd_scree <- function(variance_df,
                           title = "Variance explained by SVD components") {
  if (!all(c("Component", "ExplainedVariance") %in% names(variance_df))) {
    stop("variance_df must contain 'Component' and 'ExplainedVariance' columns.")
  }
  
  ggplot2::ggplot(
    variance_df,
    ggplot2::aes(x = Component, y = ExplainedVariance)
  ) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::labs(
      title = title,
      x     = "Component",
      y     = "Variance explained"
    ) +
    ggplot2::theme_minimal()
}

# Mean ± SD of a given PC by condition
plot_svd_component_means <- function(svd_result,
                                     component = 1L,
                                     palette   = NULL,
                                     title     = NULL) {
  scores <- svd_result$sample_scores
  
  comp_idx  <- as.integer(component)
  comp_name <- paste0("PC", comp_idx)
  
  if (!comp_name %in% colnames(scores)) {
    stop("Component column not found in sample_scores: ", comp_name)
  }
  if (!"Condition" %in% colnames(scores)) {
    stop("'sample_scores' in svd_result must contain a 'Condition' column.")
  }
  
  df <- scores[, c("Condition", comp_name)]
  colnames(df) <- c("Condition", "Score")
  
  df$Condition <- factor(df$Condition, levels = unique(df$Condition))
  
  summary_df <- df |>
    dplyr::group_by(Condition) |>
    dplyr::summarise(
      mean_score = mean(Score, na.rm = TRUE),
      sd_score   = stats::sd(Score, na.rm = TRUE),
      .groups    = "drop"
    )
  
  if (is.null(title)) {
    title <- sprintf("Mean of %s by Condition", comp_name)
  }
  y_lab <- sprintf("Component %s Value", comp_name)
  
  p <- ggplot2::ggplot(
    summary_df,
    ggplot2::aes(x = Condition, y = mean_score, color = Condition)
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = mean_score - sd_score,
        ymax = mean_score + sd_score
      ),
      width = 0.15
    ) +
    ggplot2::labs(
      title = title,
      x     = "Condition",
      y     = y_lab,
      color = "Condition"
    ) +
    ggplot2::theme_bw()
  
  if (!is.null(palette)) {
    p <- p + ggplot2::scale_color_manual(values = palette)
  }
  
  p
}

## ----------------------------------------------------------------------
## UpSet plot using ComplexUpset (colored intersections)
## ----------------------------------------------------------------------

plot_upset_complex <- function(bm,
                               outfile,
                               title      = NULL,
                               set_colors = NULL) {

  if (!requireNamespace("ComplexUpset", quietly = TRUE)) {
    stop("Package 'ComplexUpset' is required.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  df   <- as.data.frame(bm, stringsAsFactors = FALSE)
  sets_all <- colnames(df)

  ## 1) Fix set (row) order: F, E, T if possible
  if (all(c("FUS", "EWSR1", "TAF") %in% sets_all)) {
    sets <- c("FUS", "EWSR1", "TAF")
  } else {
    sets <- sets_all
  }

  ## Reorder columns so membership columns are in 'sets' order
  df <- df[, c(sets, setdiff(colnames(df), sets)), drop = FALSE]

  ## Convert membership columns to logical
  df[sets] <- lapply(df[sets], function(x) {
    if (is.logical(x)) {
      x
    } else if (is.numeric(x)) {
      x != 0
    } else if (is.factor(x)) {
      as.integer(x) != 0
    } else {
      x != 0
    }
  })

  ## 2) Single-set vs multi-set classification for bar coloring
  df$single_set <- apply(df[sets], 1L, function(row) {
    on <- which(row)
    if (length(on) == 1L) sets[on] else "multi"
  })
  df$single_set <- factor(df$single_set, levels = c(sets, "multi"))

  ## Colors: single-set intersections use KO colors, multi = black
  col_vec <- rep("black", length(sets))
  names(col_vec) <- sets

  if (!is.null(set_colors)) {
    sc <- set_colors
    if (is.list(sc)) sc <- unlist(sc, use.names = TRUE)
    if (!is.null(names(sc))) {
      common <- intersect(names(sc), sets)
      col_vec[common] <- sc[common]
    }
  }

  scale_vals <- c(col_vec, multi = "black")

  ## 3) Explicit intersection order: FET, FE, ET, FT, F, E, T
  custom_intersections <- NULL
  if (length(sets) == 3L) {
    custom_intersections <- list(
      sets,                     # F E T   (triple)
      sets[c(1, 2)],            # F E
      sets[c(2, 3)],            # E T
      sets[c(1, 3)],            # F T
      sets[1],                  # F
      sets[2],                  # E
      sets[3]                   # T
    )
  }

  ## 4) Build the plot
  p <- ComplexUpset::upset(
    df,
    intersect            = sets,
    intersections        = custom_intersections,
    sort_intersections   = FALSE,      # respect custom_intersections order
    mode                 = "exclusive_intersection",
    set_sizes            = FALSE,
    encode_sets          = FALSE,
    base_annotations     = list(
      "Intersection size" =
        ComplexUpset::intersection_size(
          mode                = "exclusive_intersection",
          mapping             = ggplot2::aes(fill = single_set),
          counts              = TRUE,
          bar_number_threshold = 2,    # always place numbers above bars
          text                = list(size = 8)
        ) +
        ggplot2::scale_fill_manual(values = scale_vals, guide = "none")
    )
  )

  if (!is.null(title) && nzchar(title)) {
    p <- p + ggplot2::ggtitle(title)
  }

  ## Larger fonts everywhere (~3x default)
  p <- p +
    ggplot2::theme(
      text        = ggplot2::element_text(size = 15),
      axis.text.x = ggplot2::element_text(size = 15),
      axis.text.y = ggplot2::element_text(size = 17),
      axis.title  = ggplot2::element_text(size = 14)
    )

  ggplot2::ggsave(
    filename = outfile,
    plot     = p,
    width    = 12,
    height   = 10,
    units    = "in",
    dpi      = 300
  )

  invisible(outfile)
}

## ----------------------------------------------------------------------
## GO over-representation convenience
## ----------------------------------------------------------------------

run_enrich_go_for_result <- function(res_df,
                                     ontology,
                                     padj_threshold,
                                     q_value_threshold,
                                     simplify_cutoff,
                                     orgdb = org.Hs.eg.db,
                                     count_min = 0,
                                     count_max = Inf,
                                     id_col = "ensembl_id",
                                     id_keytype = c("auto", "ENSEMBL", "SYMBOL")) {
  # id_col: column in res_df containing gene identifiers (symbols or Ensembl)
  # id_keytype:
  #   "auto"   -> try to infer whether the IDs are ENSEMBL or SYMBOL
  #   "ENSEMBL" -> force treat as Ensembl IDs
  #   "SYMBOL"  -> force treat as gene symbols

  id_keytype <- match.arg(id_keytype)

  if (!is.data.frame(res_df)) {
    stop("'res_df' must be a data.frame.")
  }
  if (!id_col %in% names(res_df)) {
    stop("Column '", id_col, "' not found in 'res_df'.")
  }

  # Extract unique gene IDs -----------------------------------------------
  genes <- as.character(res_df[[id_col]])
  genes <- unique(na.omit(genes))
  if (!length(genes)) {
    return(NULL)
  }

  # Decide keytype for AnnotationDbi::mapIds ------------------------------
  decide_keytype <- function(ids) {
    # Very simple heuristic: if majority look like ENSEMBL, use ENSEMBL,
    # otherwise use SYMBOL.
    is_ens <- grepl("^ENSG[0-9]+", ids)
    if (sum(is_ens) >= 1L && sum(is_ens) / length(ids) >= 0.5) {
      "ENSEMBL"
    } else {
      "SYMBOL"
    }
  }

  kt <- switch(
    id_keytype,
    auto    = decide_keytype(genes),
    ENSEMBL = "ENSEMBL",
    SYMBOL  = "SYMBOL"
  )

  map_with_keytype <- function(keytype) {
    AnnotationDbi::mapIds(
      orgdb,
      keys      = genes,
      column    = "ENTREZID",
      keytype   = keytype,
      multiVals = "first"
    )
  }

  # Try mapping; if auto and first keytype fails completely, fall back ----
  entrez_raw <- tryCatch(
    map_with_keytype(kt),
    error = function(e) {
      if (id_keytype == "auto") {
        alt_kt <- if (kt == "ENSEMBL") "SYMBOL" else "ENSEMBL"
        tryCatch(
          map_with_keytype(alt_kt),
          error = function(e2) {
            # If both fail, return NAs and handle downstream
            rep(NA_character_, length(genes))
          }
        )
      } else {
        stop(e)
      }
    }
  )

  entrez <- unique(na.omit(unname(entrez_raw)))
  if (!length(entrez)) {
    # Nothing mappable; just return NULL gracefully
    return(NULL)
  }

  # Run enrichGO -----------------------------------------------------------
  ego <- clusterProfiler::enrichGO(
    gene          = entrez,
    OrgDb         = orgdb,
    keyType       = "ENTREZID",
    ont           = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff  = padj_threshold,
    qvalueCutoff  = q_value_threshold,
    readable      = TRUE
  )
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0L) {
    return(NULL)
  }

  # Simplify and filter by Count ------------------------------------------
  ego_simpl <- clusterProfiler::simplify(
    ego,
    cutoff     = simplify_cutoff,
    by         = "p.adjust",
    select_fun = min
  )

  if (!is.null(ego_simpl) && nrow(as.data.frame(ego_simpl)) > 0L) {
    ego_df <- as.data.frame(ego_simpl)

    if ("Count" %in% colnames(ego_df)) {
      ego_df <- ego_df[
        ego_df$Count >= count_min & ego_df$Count <= count_max,
        ,
        drop = FALSE
      ]
    }

    if (!nrow(ego_df)) {
      return(NULL)
    }

    ego_simpl@result <- ego_df
  }

  ego_simpl
}



## Helper: barplot of top N GO terms for one direction and one contrast
plot_go_top_terms <- function(go_df,
                              n_top,
                              title,
                              fill_color,
                              outfile) {
  # Nothing to plot
  if (is.null(go_df) || !nrow(go_df)) {
    return(invisible(NULL))
  }

  # Order by adjusted p-value (most significant first)
  go_df <- go_df[order(go_df$p.adjust), , drop = FALSE]

  n_top <- min(as.integer(n_top), nrow(go_df))
  if (!is.finite(n_top) || n_top < 1L) {
    return(invisible(NULL))
  }

  # Keep only top N
  go_df <- go_df[seq_len(n_top), , drop = FALSE]

  # -log10(adjusted p-value) for x-axis
  go_df$neg_log10_padj <- -log10(pmax(go_df$p.adjust, .Machine$double.xmin))

  # Factor for y-axis ordering (most significant on top)
  go_df$Description <- factor(
    go_df$Description,
    levels = rev(go_df$Description)
  )

  ## --------------------------------------------------------------------
  ## Compute n and r for labels: n = Count, r = Count / total_genes
  ## total_genes is taken from the denominator of GeneRatio (Count / N)
  ## --------------------------------------------------------------------
  total_genes <- NA_real_

  if ("GeneRatio" %in% colnames(go_df)) {
    gr <- as.character(go_df$GeneRatio)
    gr <- gr[!is.na(gr)]
    if (length(gr) >= 1L && grepl("/", gr[1L], fixed = TRUE)) {
      parts <- strsplit(gr[1L], "/", fixed = TRUE)[[1L]]
      if (length(parts) == 2L) {
        denom <- suppressWarnings(as.numeric(parts[2L]))
        if (!is.na(denom) && denom > 0) {
          total_genes <- denom
        }
      }
    }
  }

  # Fallback: if GeneRatio is missing or unparsable, use max Count
  if (is.na(total_genes) && "Count" %in% colnames(go_df)) {
    total_genes <- max(go_df$Count, na.rm = TRUE)
  }

  if (!is.na(total_genes) && "Count" %in% colnames(go_df)) {
    go_df$n_term <- go_df$Count
    go_df$r_term <- go_df$n_term / total_genes
    go_df$label  <- sprintf(
      "n=%d  r=%.1f%%",
      go_df$n_term,
      100 * go_df$r_term
    )
  } else {
    go_df$label <- NA_character_
  }

  # Range for x-axis, extended a bit to make space for labels
  max_x <- max(go_df$neg_log10_padj, na.rm = TRUE)
  if (!is.finite(max_x) || max_x <= 0) {
    max_x <- 1
  }

  # p <- ggplot2::ggplot(
  #   go_df,
  #   ggplot2::aes(x = neg_log10_padj, y = Description)
  # ) +
  #   ggplot2::geom_col(fill = fill_color) +
  #   ggplot2::geom_text(
  #     ggplot2::aes(label = label),
  #     hjust = -0.1,    # push label slightly to the right of bar end
  #     size  = 3,
  #     na.rm = TRUE
  #   ) +
  #   ggplot2::coord_cartesian(xlim = c(0, max_x * 1.15)) +
  #   ggplot2::labs(
  #     x     = expression(-log[10](adjusted~p)),
  #     y     = "GO term",
  #     title = title
  #   ) +
  #   ggplot2::theme_minimal(base_size = 12) +
  #   ggplot2::theme(
  #     panel.grid  = ggplot2::element_blank(),
  #     axis.text.y = ggplot2::element_text(size = 8)
  #   )

  wrapped_title <- wrap_title_once(title, width = 35L)

  p <- ggplot2::ggplot(
      go_df,
      ggplot2::aes(x = neg_log10_padj, y = Description)
    ) +
      ggplot2::geom_col(fill = fill_color) +
      ggplot2::geom_text(
        ggplot2::aes(label = label),
        hjust = -0.1,
        size  = 3,
        na.rm = TRUE
      ) +
      ggplot2::coord_cartesian(xlim = c(0, max_x * 1.15)) +
      ggplot2::labs(
        x     = expression(-log[10](adjusted~p)),
        y     = "GO term",
        title = wrapped_title
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        panel.grid  = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(size = 8),
        plot.title  = ggplot2::element_text(
          size = ggplot2::rel(1.6),
          face = "bold"
        )
      )

  ggplot2::ggsave(
    filename = outfile,
    plot     = p,
    width    = 10,
    height   = max(3, 0.3 * n_top),
    dpi      = 300
  )
  invisible(p)
}



## ----------------------------------------------------------------------
## FGSEA convenience
## ----------------------------------------------------------------------

prepare_fgsea_ranks <- function(res_df,
                                lfc_col = "log2FoldChange",
                                symbol_col = "symbol") {
  if (!all(c(lfc_col, symbol_col) %in% names(res_df))) {
    stop("res_df must contain columns '", lfc_col, "' and '", symbol_col, "'.")
  }
  df <- res_df[!is.na(res_df[[lfc_col]]) & !is.na(res_df[[symbol_col]]), ]
  df <- df[!duplicated(df[[symbol_col]]), ]
  ranks <- df[[lfc_col]]
  names(ranks) <- df[[symbol_col]]
  sort(ranks, decreasing = TRUE)
}

## ----------------------------------------------------------------------
## Set / matrix helpers for UpSet and Venn
## ----------------------------------------------------------------------

build_binary_matrix_from_sets <- function(sets_list) {
  sets_list <- sets_list[!vapply(sets_list, is.null, logical(1))]
  if (!length(sets_list)) {
    stop("No sets provided.")
  }
  items <- sort(unique(unlist(sets_list, use.names = FALSE)))
  mat <- matrix(
    0L,
    nrow = length(items),
    ncol = length(sets_list),
    dimnames = list(items, names(sets_list))
  )
  for (nm in names(sets_list)) {
    mat[sets_list[[nm]], nm] <- 1L
  }
  as.data.frame(mat)
}

save_venn <- function(sets,
                      outfile,
                      fills       = NULL,
                      outline_col = "white",
                      alpha       = 0.6,
                      size        = 1200,
                      res         = 200) {
  n_sets <- length(sets)
  if (n_sets < 2L || n_sets > 4L) {
    message("Venn diagram only drawn for 2-4 sets. Got: ", n_sets)
    return(invisible(NULL))
  }
  if (is.null(fills)) {
    fills <- grDevices::rainbow(n_sets)
  }

  ## Keep VennDiagram from writing a log file to disk
  if (requireNamespace("futile.logger", quietly = TRUE)) {
    futile.logger::flog.threshold(
      futile.logger::ERROR,
      name = "VennDiagramLogger"
    )
    futile.logger::flog.appender(
      futile.logger::appender.console(),
      name = "VennDiagramLogger"
    )
  }

  vg <- VennDiagram::venn.diagram(
    x              = sets,
    category.names = names(sets),
    filename       = NULL,   # return a grob
    output         = FALSE,  # do not write image file directly
    fill           = fills,
    alpha          = alpha,
    lwd            = 2,
    col            = outline_col,
    cex            = 2.2,
    fontfamily     = "sans",
    cat.cex        = 2.2,
    cat.col        = fills
  )

  grDevices::png(outfile, width = size, height = size, res = res)
  grid::grid.newpage()
  grid::grid.draw(vg)
  grDevices::dev.off()

  invisible(outfile)
}

## ----------------------------------------------------------------------
## FGSEA plotting helpers
## ----------------------------------------------------------------------

# Dotplot of top enriched pathways for a single contrast.
plot_fgsea_top_pathways <- function(fgsea_res,
                                    top_n   = 20L,
                                    padj_max = 0.05,
                                    title   = "Top enriched pathways (fgsea)") {
  if (is.null(fgsea_res) || !nrow(fgsea_res)) {
    message("fgsea_res is empty; no dotplot will be produced.")
    return(NULL)
  }
  
  required_cols <- c("pathway", "NES", "padj", "size")
  if (!all(required_cols %in% colnames(fgsea_res))) {
    stop(
      "fgsea_res must contain columns: ",
      paste(required_cols, collapse = ", ")
    )
  }
  
  df <- fgsea_res
  df <- df[order(df$padj, -abs(df$NES)), , drop = FALSE]
  
  if (!is.null(padj_max) && is.finite(padj_max)) {
    df <- df[df$padj <= padj_max, , drop = FALSE]
  }
  
  if (!nrow(df)) {
    message("No pathways pass padj_max = ", padj_max, " for fgsea dotplot.")
    return(NULL)
  }
  
  if (!is.null(top_n) && top_n > 0L && nrow(df) > top_n) {
    df <- df[seq_len(top_n), , drop = FALSE]
  }
  
  df$pathway <- factor(df$pathway, levels = rev(df$pathway))
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = NES, y = pathway)) +
    ggplot2::geom_point(ggplot2::aes(size = size, color = padj)) +
    ggplot2::scale_color_continuous(trans = "reverse") +
    ggplot2::labs(
      title = title,
      x     = "Normalized Enrichment Score (NES)",
      y     = NULL,
      size  = "Genes in set",
      color = "padj"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 8)
    )
  
  p
}

# Build NES matrix (pathways x contrasts) from FGSEA_RESULTS.
build_fgsea_nes_matrix <- function(fgsea_results,
                                   padj_max     = 0.05,
                                   abs_nes_min  = 0,
                                   max_pathways = 50L) {
  if (is.null(fgsea_results) || !length(fgsea_results)) {
    message("fgsea_results is empty; NES matrix cannot be built.")
    return(NULL)
  }
  
  rows <- list()
  for (cond in names(fgsea_results)) {
    res <- fgsea_results[[cond]]
    if (is.null(res) || !nrow(res)) {
      next
    }
    required_cols <- c("pathway", "NES", "padj")
    if (!all(required_cols %in% colnames(res))) {
      stop(
        "Each fgsea result must contain columns: ",
        paste(required_cols, collapse = ", ")
      )
    }
    res$contrast <- cond
    rows[[cond]] <- res
  }
  
  if (!length(rows)) {
    message("No non-empty fgsea results to build NES matrix.")
    return(NULL)
  }
  
  df <- do.call(rbind, rows)
  df <- df[!is.na(df$padj) & !is.na(df$NES), , drop = FALSE]
  
  if (!is.null(padj_max) && is.finite(padj_max)) {
    df <- df[df$padj <= padj_max, , drop = FALSE]
  }
  
  if (!is.null(abs_nes_min) && abs_nes_min > 0) {
    df <- df[abs(df$NES) >= abs_nes_min, , drop = FALSE]
  }
  
  if (!nrow(df)) {
    message("No pathways pass fgsea filtering for NES heatmap.")
    return(NULL)
  }
  
  if (!is.null(max_pathways) && max_pathways > 0L) {
    min_p <- tapply(df$padj, df$pathway, min, na.rm = TRUE)
    min_p <- sort(min_p)
    sel_paths <- names(min_p)[seq_len(min(length(min_p), max_pathways))]
    df <- df[df$pathway %in% sel_paths, , drop = FALSE]
  }
  
  pathways  <- sort(unique(df$pathway))
  contrasts <- sort(unique(df$contrast))
  
  mat <- matrix(
    NA_real_,
    nrow = length(pathways),
    ncol = length(contrasts),
    dimnames = list(pathways, contrasts)
  )
  
  for (i in seq_len(nrow(df))) {
    pw <- df$pathway[i]
    ct <- df$contrast[i]
    mat[pw, ct] <- df$NES[i]
  }
  
  mat
}

# Create a ComplexHeatmap::Heatmap object for NES values.
make_fgsea_nes_heatmap <- function(nes_mat,
                                   name         = "NES",
                                   column_title = "FGSEA normalized enrichment scores") {
  if (is.null(nes_mat)) {
    return(NULL)
  }
  if (!is.matrix(nes_mat) || !is.numeric(nes_mat)) {
    stop("nes_mat must be a numeric matrix.")
  }
  if (!nrow(nes_mat) || !ncol(nes_mat)) {
    message("NES matrix is empty; nothing to plot.")
    return(NULL)
  }
  
  nes_lim <- max(abs(nes_mat), na.rm = TRUE)
  if (!is.finite(nes_lim) || nes_lim <= 0) {
    nes_lim <- 1
  }
  
  col_fun <- circlize::colorRamp2(
    c(-nes_lim, 0, nes_lim),
    c("#1E88E5", "white", "#FFC107")  # consistent with DE heatmap
  )
  
  ComplexHeatmap::Heatmap(
    nes_mat,
    name            = name,
    col             = col_fun,
    cluster_rows    = TRUE,
    cluster_columns = TRUE,
    column_title    = column_title,
    na_col          = "grey90"
  )
}
