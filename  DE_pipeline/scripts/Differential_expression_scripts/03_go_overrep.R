# ## scripts/04_go_overrep.R
# ## GO over-representation (enrichGO) for all / pos / neg DEGs.

# if (!exists("DE_RESULTS"))         stop("'DE_RESULTS' not found.")
# if (!exists("DE_LFC"))             stop("'DE_LFC' not found.")
# if (!exists("DE_PADJ"))            stop("'DE_PADJ' not found.")
# if (!exists("DE_QVAL"))            stop("'DE_QVAL' not found.")
# if (!exists("GO_ONTOLOGY"))        stop("'GO_ONTOLOGY' not found.")
# if (!exists("GO_SIMPLIFY_CUTOFF")) stop("'GO_SIMPLIFY_CUTOFF' not found.")
# if (!exists("GO_COUNT_MIN"))       stop("'GO_COUNT_MIN' not found.")
# if (!exists("GO_COUNT_MAX"))       stop("'GO_COUNT_MAX' not found.")
# if (!exists("output_dir"))         stop("'output_dir' not found.")
# if (!exists("GO_TOP_N")) {
#   GO_TOP_N <- 5L
# }

# out_dir_go <- ensure_dir(output_dir, "05_go_overrep")

# conds <- names(DE_RESULTS)

# GO_ALL_LIST <- list()
# GO_POS_LIST <- list()
# GO_NEG_LIST <- list()

# GO_ALL_DF_LIST <- list()
# GO_POS_DF_LIST <- list()
# GO_NEG_DF_LIST <- list()

# for (cond in conds) {
#   res_df <- DE_RESULTS[[cond]]
  
#   ego_all <- run_enrich_go_for_result(
#     res_df          = res_df,
#     direction       = "all",
#     ontology        = GO_ONTOLOGY,
#     lfc_threshold   = DE_LFC,
#     padj_threshold  = DE_PADJ,
#     q_value_threshold = DE_QVAL,
#     simplify_cutoff = GO_SIMPLIFY_CUTOFF,
#     orgdb           = org.Hs.eg.db,
#     count_min       = GO_COUNT_MIN,
#     count_max       = GO_COUNT_MAX
#   )
  
#   ego_pos <- run_enrich_go_for_result(
#     res_df          = res_df,
#     direction       = "pos",
#     ontology        = GO_ONTOLOGY,
#     lfc_threshold   = DE_LFC,
#     padj_threshold  = DE_PADJ,
#     q_value_threshold = DE_QVAL,
#     simplify_cutoff = GO_SIMPLIFY_CUTOFF,
#     orgdb           = org.Hs.eg.db,
#     count_min       = GO_COUNT_MIN,
#     count_max       = GO_COUNT_MAX
#   )
  
#   ego_neg <- run_enrich_go_for_result(
#     res_df          = res_df,
#     direction       = "neg",
#     ontology        = GO_ONTOLOGY,
#     lfc_threshold   = DE_LFC,
#     padj_threshold  = DE_PADJ,
#     q_value_threshold = DE_QVAL,
#     simplify_cutoff = GO_SIMPLIFY_CUTOFF,
#     orgdb           = org.Hs.eg.db,
#     count_min       = GO_COUNT_MIN,
#     count_max       = GO_COUNT_MAX
#   )
  
#   GO_ALL_LIST[[cond]] <- ego_all
#   GO_POS_LIST[[cond]] <- ego_pos
#   GO_NEG_LIST[[cond]] <- ego_neg

#   if (!is.null(ego_all)) {
#     df_all <- as.data.frame(ego_all)
#     GO_ALL_DF_LIST[[cond]] <- df_all
#     openxlsx::write.xlsx(
#       df_all,
#       file.path(out_dir_go, paste0("GO_", cond, "_ALL.xlsx")),
#       rowNames = FALSE
#     )

#     plot_go_top_terms(
#     go_df     = df_all,
#     n_top     = GO_TOP_N,
#     title     = paste0("GO over-representation (ALL) - ", cond),
#     fill_color = "#4CAF50",  # green
#     outfile   = file.path(
#       out_dir_go,
#       paste0("GO_", cond, "_ALL_top", GO_TOP_N, ".png")
#     )
#   )
#   }
#   if (!is.null(ego_pos)) {
#     df_pos <- as.data.frame(ego_pos)
#     GO_POS_DF_LIST[[cond]] <- df_pos
#     openxlsx::write.xlsx(
#       df_pos,
#       file.path(out_dir_go, paste0("GO_", cond, "_POS.xlsx")),
#       rowNames = FALSE
#     )
#     plot_go_top_terms(
#     go_df     = df_pos,
#     n_top     = GO_TOP_N,
#     title     = paste0("GO over-representation (POS) - ", cond),
#     fill_color = "#FFC107",  # yellow
#     outfile   = file.path(
#       out_dir_go,
#       paste0("GO_", cond, "_POS_top", GO_TOP_N, ".png")
#     )
#   )
#   }
#   if (!is.null(ego_neg)) {
#     df_neg <- as.data.frame(ego_neg)
#     GO_NEG_DF_LIST[[cond]] <- df_neg
#     openxlsx::write.xlsx(
#       df_neg,
#       file.path(out_dir_go, paste0("GO_", cond, "_NEG.xlsx")),
#       rowNames = FALSE
#     )
#     plot_go_top_terms(
#     go_df     = df_neg,
#     n_top     = GO_TOP_N,
#     title     = paste0("GO over-representation (NEG) - ", cond),
#     fill_color = "#1E88E5",  # blue
#     outfile   = file.path(
#       out_dir_go,
#       paste0("GO_", cond, "_NEG_top", GO_TOP_N, ".png")
#     )
#   )
#   }
# }

# ## Expose for downstream scripts (GO overlap, z-score heatmaps)
# ego_pos_list <- GO_POS_DF_LIST
# ego_neg_list <- GO_NEG_DF_LIST

# saveRDS(GO_ALL_LIST, file.path(out_dir_go, "GO_ALL_list.rds"))
# saveRDS(GO_POS_LIST, file.path(out_dir_go, "GO_POS_list.rds"))
# saveRDS(GO_NEG_LIST, file.path(out_dir_go, "GO_NEG_list.rds"))

## scripts/04_go_overrep.R
## GO over-representation (enrichGO) for all / pos / neg DEGs.

if (!exists("DE_RESULTS"))         stop("'DE_RESULTS' not found.")
if (!exists("DE_LFC"))             stop("'DE_LFC' not found.")
if (!exists("DE_PADJ"))            stop("'DE_PADJ' not found.")
if (!exists("DE_QVAL"))            stop("'DE_QVAL' not found.")
if (!exists("GO_ONTOLOGY"))        stop("'GO_ONTOLOGY' not found.")
if (!exists("GO_SIMPLIFY_CUTOFF")) stop("'GO_SIMPLIFY_CUTOFF' not found.")
if (!exists("GO_COUNT_MIN"))       stop("'GO_COUNT_MIN' not found.")
if (!exists("GO_COUNT_MAX"))       stop("'GO_COUNT_MAX' not found.")
if (!exists("output_dir"))         stop("'output_dir' not found.")
if (!exists("GO_TOP_N")) {
  GO_TOP_N <- 5L
}

out_dir_go <- ensure_dir(output_dir, "05_go_overrep")

conds <- names(DE_RESULTS)

GO_ALL_LIST    <- list()
GO_POS_LIST    <- list()
GO_NEG_LIST    <- list()

GO_ALL_DF_LIST <- list()
GO_POS_DF_LIST <- list()
GO_NEG_DF_LIST <- list()

for (cond in conds) {
  res_df <- DE_RESULTS[[cond]]

  # Basic sanity check on required columns ---------------------------------
  required_cols <- c("log2FoldChange", "padj", "ensembl_id")
  if (!all(required_cols %in% names(res_df))) {
    warning(
      "Skipping condition '", cond,
      "': required columns ", paste(required_cols, collapse = ", "),
      " not found."
    )
    next
  }

  ## ----------------------------------------------------------------------
  ## 1) Filter DEGs per direction (ALL / POS / NEG)
  ##    This defines which genes enter the GO test set.
  ## ----------------------------------------------------------------------

  idx_all <- which(
    res_df$padj < DE_PADJ &
      abs(res_df$log2FoldChange) >= DE_LFC
  )
  idx_pos <- which(
    res_df$padj < DE_PADJ &
      res_df$log2FoldChange >= DE_LFC
  )
  idx_neg <- which(
    res_df$padj < DE_PADJ &
      res_df$log2FoldChange <= -DE_LFC
  )

  res_all <- res_df[idx_all, , drop = FALSE]
  res_pos <- res_df[idx_pos, , drop = FALSE]
  res_neg <- res_df[idx_neg, , drop = FALSE]

  ## ----------------------------------------------------------------------
  ## 2) Run GO for each filtered set using the generic helper
  ## ----------------------------------------------------------------------

  ego_all <- if (nrow(res_all)) {
    run_enrich_go_for_result(
      res_df           = res_all,
      ontology         = GO_ONTOLOGY,
      padj_threshold   = DE_PADJ,
      q_value_threshold = DE_QVAL,
      simplify_cutoff  = GO_SIMPLIFY_CUTOFF,
      orgdb            = org.Hs.eg.db,
      count_min        = GO_COUNT_MIN,
      count_max        = GO_COUNT_MAX
    )
  } else {
    NULL
  }

  ego_pos <- if (nrow(res_pos)) {
    run_enrich_go_for_result(
      res_df           = res_pos,
      ontology         = GO_ONTOLOGY,
      padj_threshold   = DE_PADJ,
      q_value_threshold = DE_QVAL,
      simplify_cutoff  = GO_SIMPLIFY_CUTOFF,
      orgdb            = org.Hs.eg.db,
      count_min        = GO_COUNT_MIN,
      count_max        = GO_COUNT_MAX
    )
  } else {
    NULL
  }

  ego_neg <- if (nrow(res_neg)) {
    run_enrich_go_for_result(
      res_df           = res_neg,
      ontology         = GO_ONTOLOGY,
      padj_threshold   = DE_PADJ,
      q_value_threshold = DE_QVAL,
      simplify_cutoff  = GO_SIMPLIFY_CUTOFF,
      orgdb            = org.Hs.eg.db,
      count_min        = GO_COUNT_MIN,
      count_max        = GO_COUNT_MAX
    )
  } else {
    NULL
  }

  GO_ALL_LIST[[cond]] <- ego_all
  GO_POS_LIST[[cond]] <- ego_pos
  GO_NEG_LIST[[cond]] <- ego_neg

  ## ----------------------------------------------------------------------
  ## 3) Save tables + plots (same as before)
  ## ----------------------------------------------------------------------

  if (!is.null(ego_all)) {
    df_all <- as.data.frame(ego_all)
    GO_ALL_DF_LIST[[cond]] <- df_all

    openxlsx::write.xlsx(
      df_all,
      file.path(out_dir_go, paste0("GO_", cond, "_ALL.xlsx")),
      rowNames = FALSE
    )

    plot_go_top_terms(
      go_df      = df_all,
      n_top      = GO_TOP_N,
      title      = paste0("GO over-representation (ALL) - ", cond),
      fill_color = "#4CAF50",  # green
      outfile    = file.path(
        out_dir_go,
        paste0("GO_", cond, "_ALL_top", GO_TOP_N, ".png")
      )
    )
  }

  if (!is.null(ego_pos)) {
    df_pos <- as.data.frame(ego_pos)
    GO_POS_DF_LIST[[cond]] <- df_pos

    openxlsx::write.xlsx(
      df_pos,
      file.path(out_dir_go, paste0("GO_", cond, "_POS.xlsx")),
      rowNames = FALSE
    )

    plot_go_top_terms(
      go_df      = df_pos,
      n_top      = GO_TOP_N,
      title      = paste0("GO over-representation (POS) - ", cond),
      fill_color = "#FFC107",  # amber
      outfile    = file.path(
        out_dir_go,
        paste0("GO_", cond, "_POS_top", GO_TOP_N, ".png")
      )
    )
  }

  if (!is.null(ego_neg)) {
    df_neg <- as.data.frame(ego_neg)
    GO_NEG_DF_LIST[[cond]] <- df_neg

    openxlsx::write.xlsx(
      df_neg,
      file.path(out_dir_go, paste0("GO_", cond, "_NEG.xlsx")),
      rowNames = FALSE
    )

    plot_go_top_terms(
      go_df      = df_neg,
      n_top      = GO_TOP_N,
      title      = paste0("GO over-representation (NEG) - ", cond),
      fill_color = "#1E88E5",  # blue
      outfile    = file.path(
        out_dir_go,
        paste0("GO_", cond, "_NEG_top", GO_TOP_N, ".png")
      )
    )
  }
}

## Expose for downstream scripts (GO overlap, z-score heatmaps)
ego_pos_list <- GO_POS_DF_LIST
ego_neg_list <- GO_NEG_DF_LIST

saveRDS(GO_ALL_LIST, file.path(out_dir_go, "GO_ALL_list.rds"))
saveRDS(GO_POS_LIST, file.path(out_dir_go, "GO_POS_list.rds"))
saveRDS(GO_NEG_LIST, file.path(out_dir_go, "GO_NEG_list.rds"))
# stop("GO over-representation analysis completed successfully.")