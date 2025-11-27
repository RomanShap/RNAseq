## scripts/07_go_overlap_plots.R
## GO-level overlaps (positive / negative) + z-score heatmaps.

if (!exists("ego_pos_list"))   stop("'ego_pos_list' not found (run 04_go_overrep.R first).")
if (!exists("ego_neg_list"))   stop("'ego_neg_list' not found (run 04_go_overrep.R first).")
if (!exists("GO_COUNT_MIN"))   stop("'GO_COUNT_MIN' not found.")
if (!exists("GO_COUNT_MAX"))   stop("'GO_COUNT_MAX' not found.")
if (!exists("output_dir"))     stop("'output_dir' not found.")

out_dir_go_ov <- ensure_dir(output_dir, "08_go_overlap")

## Helper: build sets of GO descriptions (filtered by Count again, just in case)
go_sets_from_list <- function(ego_df_list) {
  sets <- list()
  for (nm in names(ego_df_list)) {
    df <- ego_df_list[[nm]]
    if (is.null(df) || !nrow(df)) next
    if (!"Description" %in% names(df) || !"Count" %in% names(df)) next
    df <- df[df$Count >= GO_COUNT_MIN & df$Count <= GO_COUNT_MAX, , drop = FALSE]
    if (!nrow(df)) next
    sets[[nm]] <- unique(df$Description)
  }
  sets
}

pos_sets <- go_sets_from_list(ego_pos_list)
neg_sets <- go_sets_from_list(ego_neg_list)

# if (length(pos_sets)) {
#   bm_pos <- build_binary_matrix_from_sets(pos_sets)
#   openxlsx::write.xlsx(
#     data.frame(
#       GO_term   = rownames(bm_pos),
#       bm_pos,
#       stringsAsFactors = FALSE
#     ),
#     file.path(out_dir_go_ov, "GO_overlap_POS.xlsx"),
#     rowNames = FALSE
#   )
  
#   png(file.path(out_dir_go_ov, "UpSet_GO_pos.png"), width = 12, height = 8, units = "in", res = 300)
#   grid::grid.newpage()
#   UpSetR::upset(
#     bm_pos,
#     sets             = colnames(bm_pos),
#     keep.order       = TRUE,
#     order.by         = "degree",
#     decreasing       = c(TRUE, TRUE),
#     mainbar.y.label  = "Intersection size",
#     sets.x.label     = "Set size",
#     matrix.color     = "black",
#     main.bar.color   = "black",
#     point.size       = 4.0,
#     line.size        = 1.0,
#     mb.ratio         = c(0.70, 0.30),
#     text.scale       = 3
#   )
#   grid::grid.text(label = "GO terms (upregulated DEGs)", y = 0.98, gp = grid::gpar(cex = 1.4))
#   grDevices::dev.off()
  
#   if (ncol(bm_pos) <= 4L) {
#     save_venn(
#       sets    = pos_sets,
#       outfile = file.path(out_dir_go_ov, "Venn_GO_pos.png"),
#       fills   = grDevices::rainbow(length(pos_sets))
#     )
#   }
# }
if (length(pos_sets)) {
  bm_pos <- build_binary_matrix_from_sets(pos_sets)
  openxlsx::write.xlsx(
    data.frame(
      GO_term   = rownames(bm_pos),
      bm_pos,
      stringsAsFactors = FALSE
    ),
    file.path(out_dir_go_ov, "GO_overlap_POS.xlsx"),
    rowNames = FALSE
  )

  ## KO colors (EWSR1/FUS/TAF) reused from DE overlaps if available
  ko_cols <- NULL
  if (exists("COLOR_KO")) {
    ko_cols <- COLOR_KO
  }

  plot_upset_complex(
    bm         = bm_pos,
    outfile    = file.path(out_dir_go_ov, "UpSet_GO_pos.png"),
    title      = "GO terms (upregulated DEGs)",
    set_colors = ko_cols
  )
  
  if (ncol(bm_pos) <= 4L) {
    save_venn(
      sets    = pos_sets,
      outfile = file.path(out_dir_go_ov, "Venn_GO_pos.png"),
      fills   = grDevices::rainbow(length(pos_sets))
    )
  }
}


# if (length(neg_sets)) {
#   bm_neg <- build_binary_matrix_from_sets(neg_sets)
#   openxlsx::write.xlsx(
#     data.frame(
#       GO_term   = rownames(bm_neg),
#       bm_neg,
#       stringsAsFactors = FALSE
#     ),
#     file.path(out_dir_go_ov, "GO_overlap_NEG.xlsx"),
#     rowNames = FALSE
#   )
  
#   png(file.path(out_dir_go_ov, "UpSet_GO_neg.png"), width = 12, height = 8, units = "in", res = 300)
#   grid::grid.newpage()
#   UpSetR::upset(
#     bm_neg,
#     sets             = colnames(bm_neg),
#     keep.order       = TRUE,
#     order.by         = "degree",
#     decreasing       = c(TRUE, TRUE),
#     mainbar.y.label  = "Intersection size",
#     sets.x.label     = "Set size",
#     matrix.color     = "black",
#     main.bar.color   = "black",
#     point.size       = 4.0,
#     line.size        = 1.0,
#     mb.ratio         = c(0.70, 0.30),
#     text.scale       = 3
#   )
#   grid::grid.text(label = "GO terms (downregulated DEGs)", y = 0.98, gp = grid::gpar(cex = 1.4))
#   grDevices::dev.off()
  
#   if (ncol(bm_neg) <= 4L) {
#     save_venn(
#       sets    = neg_sets,
#       outfile = file.path(out_dir_go_ov, "Venn_GO_neg.png"),
#       fills   = grDevices::rainbow(length(neg_sets))
#     )
#   }
# }
if (length(neg_sets)) {
  bm_neg <- build_binary_matrix_from_sets(neg_sets)
  openxlsx::write.xlsx(
    data.frame(
      GO_term   = rownames(bm_neg),
      bm_neg,
      stringsAsFactors = FALSE
    ),
    file.path(out_dir_go_ov, "GO_overlap_NEG.xlsx"),
    rowNames = FALSE
  )

  ko_cols <- NULL
  if (exists("COLOR_KO")) {
    ko_cols <- COLOR_KO
  }

  plot_upset_complex(
    bm         = bm_neg,
    outfile    = file.path(out_dir_go_ov, "UpSet_GO_neg.png"),
    title      = "GO terms (downregulated DEGs)",
    set_colors = ko_cols
  )
  
  if (ncol(bm_neg) <= 4L) {
    save_venn(
      sets    = neg_sets,
      outfile = file.path(out_dir_go_ov, "Venn_GO_neg.png"),
      fills   = grDevices::rainbow(length(neg_sets))
    )
  }
}


## Optional: z-score heatmaps (if zScore column exists)

plot_GO_zscore_heatmap <- function(bm, ego_df_list, title, fill_colors) {
  if (!length(ego_df_list)) return(invisible(NULL))
  
  bm$Description <- rownames(bm)
  bm_long <- tidyr::pivot_longer(
    bm,
    cols      = -Description,
    names_to  = "Condition",
    values_to = "Present"
  )
  
  z_long <- dplyr::bind_rows(lapply(names(ego_df_list), function(cond) {
    df <- ego_df_list[[cond]]
    if (is.null(df) || !nrow(df) || !"zScore" %in% names(df)) {
      return(data.frame(Description = character(0), zScore = numeric(0),
                        Condition = character(0)))
    }
    data.frame(
      Description = df$Description,
      zScore      = df$zScore,
      Condition   = cond,
      stringsAsFactors = FALSE
    )
  }))
  
  if (!nrow(z_long)) return(invisible(NULL))
  
  hm_df <- dplyr::left_join(bm_long, z_long, by = c("Description", "Condition"))
  hm_df$zScore[hm_df$Present == 0L] <- NA_real_
  
  hm_df <- hm_df |>
    dplyr::group_by(Description) |>
    dplyr::filter(any(!is.na(zScore))) |>
    dplyr::ungroup()
  
  hm_df$Description <- forcats::fct_rev(
    factor(hm_df$Description, levels = unique(hm_df$Description))
  )
  
  ggplot2::ggplot(hm_df, ggplot2::aes(x = Condition, y = Description, fill = zScore)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(
      ggplot2::aes(label = ifelse(!is.na(zScore), sprintf("%.2f", zScore), "")),
      size = 3
    ) +
    ggplot2::scale_fill_gradient2(
      low  = fill_colors[1],
      mid  = fill_colors[2],
      high = fill_colors[3],
      midpoint = 0,
      na.value = "grey90",
      name = "zScore"
    ) +
    ggplot2::labs(
      title = title,
      x     = "Condition",
      y     = "GO term"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid  = ggplot2::element_blank()
    )
}

if (exists("GO_POS_DF_LIST") && length(GO_POS_DF_LIST) && length(pos_sets)) {
  bm_pos <- build_binary_matrix_from_sets(pos_sets)
  p_pos <- plot_GO_zscore_heatmap(
    bm          = bm_pos,
    ego_df_list = GO_POS_DF_LIST,
    title       = "GO zScore heatmap (upregulated DEGs)",
    fill_colors = c("#d73027", "white", "#FFC107")
  )
  if (!is.null(p_pos)) {
    ggplot2::ggsave(
      file.path(out_dir_go_ov, "GO_zScore_heatmap_upregulated.png"),
      plot = p_pos,
      width = 10,
      height = 14,
      dpi = 300
    )
  }
}

if (exists("GO_NEG_DF_LIST") && length(GO_NEG_DF_LIST) && length(neg_sets)) {
  bm_neg <- build_binary_matrix_from_sets(neg_sets)
  p_neg <- plot_GO_zscore_heatmap(
    bm          = bm_neg,
    ego_df_list = GO_NEG_DF_LIST,
    title       = "GO zScore heatmap (downregulated DEGs)",
    fill_colors = c("#d73027", "white", "#1E88E5")
  )
  if (!is.null(p_neg)) {
    ggplot2::ggsave(
      file.path(out_dir_go_ov, "GO_zScore_heatmap_downregulated.png"),
      plot = p_neg,
      width = 10,
      height = 14,
      dpi = 300
    )
  }
}
