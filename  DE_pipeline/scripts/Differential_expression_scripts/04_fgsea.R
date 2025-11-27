## scripts/05_fgsea.R
## FGSEA for each contrast using MSigDB C5.

if (!exists("DE_RESULTS"))    stop("'DE_RESULTS' not found.")
if (!exists("msigdb_c5_rds")) stop("'msigdb_c5_rds' not found.")
if (!exists("output_dir"))    stop("'output_dir' not found.")

out_dir_fgsea <- ensure_dir(output_dir, "06_fgsea")


## FGSEA parameters (fall back to defaults if globals not set)
FGSEA_NPERM        <- if (exists("FGSEA_NPERM"))        FGSEA_NPERM        else 10000L
FGSEA_TOP_N        <- if (exists("FGSEA_TOP_N"))        FGSEA_TOP_N        else 20L
FGSEA_PADJ_MAX     <- if (exists("FGSEA_PADJ_MAX"))     FGSEA_PADJ_MAX     else 0.05
FGSEA_ABS_NES_MIN  <- if (exists("FGSEA_ABS_NES_MIN"))  FGSEA_ABS_NES_MIN  else 0
FGSEA_MAX_PATHWAYS <- if (exists("FGSEA_MAX_PATHWAYS")) FGSEA_MAX_PATHWAYS else 50L

## Load or build MSigDB C5
msig_c5 <- load_or_build_msigdb_c5(msigdb_c5_rds)

## Convert to list of gene sets (gs_name -> vector of symbols)
pathways_list <- split(msig_c5$gene_symbol, msig_c5$gs_name)

FGSEA_RESULTS <- list()

for (cond in names(DE_RESULTS)) {
  res_df <- DE_RESULTS[[cond]]
  ranks  <- prepare_fgsea_ranks(res_df, lfc_col = "log2FoldChange", symbol_col = "symbol")
  
  message("Running fgsea for contrast: ", cond)
  fg <- fgsea::fgsea(
    pathways = pathways_list,
    stats    = ranks,
    nperm    = FGSEA_NPERM
  )
  fg <- fg[order(fg$padj), ]
  
  FGSEA_RESULTS[[cond]] <- fg
  
  out_xlsx <- file.path(out_dir_fgsea, paste0("fgsea_", cond, ".xlsx"))
  out_rds  <- file.path(out_dir_fgsea, paste0("fgsea_", cond, ".rds"))
  openxlsx::write.xlsx(fg, out_xlsx, rowNames = FALSE)
  saveRDS(fg, out_rds)

  ## Per-contrast dotplot of top pathways
  p_dot <- plot_fgsea_top_pathways(
    fgsea_res = fg,
    top_n     = FGSEA_TOP_N,
    padj_max  = FGSEA_PADJ_MAX,
    title     = paste0("fgsea top pathways: ", cond)
  )
  
  if (!is.null(p_dot)) {
    ggplot2::ggsave(
      filename = file.path(
        out_dir_fgsea,
        paste0("fgsea_", cond, "_top", FGSEA_TOP_N, ".png")
      ),
      plot   = p_dot,
      width  = 8,
      height = 10,
      dpi    = 300
    )
  }
}

saveRDS(FGSEA_RESULTS, file.path(out_dir_fgsea, "fgsea_results_list.rds"))

## NES heatmap across contrasts
nes_mat <- build_fgsea_nes_matrix(
  fgsea_results = FGSEA_RESULTS,
  padj_max      = FGSEA_PADJ_MAX,
  abs_nes_min   = FGSEA_ABS_NES_MIN,
  max_pathways  = FGSEA_MAX_PATHWAYS
)

if (!is.null(nes_mat)) {
  ht <- make_fgsea_nes_heatmap(
    nes_mat,
    name         = "NES",
    column_title = "FGSEA normalized enrichment scores"
  )
  
  png(
    file.path(out_dir_fgsea, "fgsea_NES_heatmap.png"),
    width  = 1600,
    height = 1800,
    res    = 300
  )
  ComplexHeatmap::draw(ht)
  grDevices::dev.off()
} else {
  message("NES heatmap not produced (no pathways passed fgsea filters).")
}
# stop("FGSEA analysis completed successfully.")