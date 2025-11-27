## scripts/03_heatmap_lfc.R
## Shared DE genes heatmap across all contrasts.

if (!exists("DE_RESULTS")) stop("'DE_RESULTS' not found.")
if (!exists("DE_LFC"))     stop("'DE_LFC' not found.")
if (!exists("DE_PADJ"))    stop("'DE_PADJ' not found.")
if (!exists("COLOR_KO"))   stop("'COLOR_KO' not found.")
if (!exists("output_dir")) stop("'output_dir' not found.")

out_dir_hm <- ensure_dir(output_dir, "04_heatmap")

THRESH_LFC <- DE_LFC
ROW_CLUSTER_CUTOFF <- 800L

conds <- names(DE_RESULTS)
if (length(conds) < 1L) stop("DE_RESULTS has no contrasts.")

sel_cols <- c("ensembl_id", "symbol", "log2FoldChange", "padj")

## Rename columns per contrast and inner-join on ensembl_id
res_named <- lapply(conds, function(cn) {
  df <- DE_RESULTS[[cn]]
  df <- df[, sel_cols]
  names(df)[names(df) == "log2FoldChange"] <- paste0("lfc_", cn)
  names(df)[names(df) == "padj"]           <- paste0("padj_", cn)
  df
})

# res_all <- Reduce(function(x, y) dplyr::inner_join(x, y, by = "ensembl_id"), res_named)
res_all <- Reduce(
    function(x, y) dplyr::full_join(x, y, by = "ensembl_id"),
    res_named
   )

## Filter: significant in at least one contrast and LFC threshold
lfc_cols  <- paste0("lfc_",  conds)
padj_cols <- paste0("padj_", conds)

lfc_mat  <- as.matrix(res_all[, lfc_cols,  drop = FALSE])
padj_mat <- as.matrix(res_all[, padj_cols, drop = FALSE])

pass_mat <- (padj_mat < DE_PADJ) & (abs(lfc_mat) >= THRESH_LFC)
keep     <- rowSums(pass_mat, na.rm = TRUE) > 0

res_filt <- res_all[keep, ]

if (!nrow(res_filt)) {
  stop("No genes pass padj < ", DE_PADJ, " and |log2FC| > ", THRESH_LFC,
       " in any contrast.")
}

gene_label <- make.unique(paste0(
  ifelse(is.na(res_filt$symbol.x) | res_filt$symbol.x == "",
         res_filt$ensembl_id,
         res_filt$symbol.x),
  "|", res_filt$ensembl_id
))

lfc_mat <- res_filt[, lfc_cols, drop = FALSE]

rownames(lfc_mat) <- gene_label
lfc_mat <- as.matrix(lfc_mat)

padj_mat <- res_filt[, padj_cols, drop = FALSE]
rownames(padj_mat) <- gene_label
padj_mat <- as.matrix(padj_mat)

lfc_for_clust <- as.matrix(lfc_mat)

# Replace non-finite values with 0 just for distance calculation
lfc_for_clust[!is.finite(lfc_for_clust)] <- 0

many_rows <- nrow(lfc_for_clust) > ROW_CLUSTER_CUTOFF

row_clust <- if (!many_rows) {
  stats::hclust(stats::dist(lfc_for_clust))
} else {
  FALSE  # or NULL if you prefer no clustering when too many rows
}

## Colors
lim <- max(2, stats::quantile(abs(lfc_mat), 0.95, na.rm = TRUE))
col_fun <- circlize::colorRamp2(
  c(-lim, 0, lim),
  c("#1E88E5", "white", "#FFC107")
)

## Annotation colors: fall back to defaults if some conditions not in COLOR_KO
cond_colors <- vapply(
  conds,
  function(cn) {
    val <- COLOR_KO[[cn]]
    if (length(val) != 1L) {
      stop("COLOR_KO entry for ", cn, " must be a single color.")
    }
    val
  },
  FUN.VALUE = character(1L)
)

names(cond_colors) <- conds

ko_ann <- ComplexHeatmap::HeatmapAnnotation(
  Condition = conds,
  col = list(Condition = cond_colors),
  annotation_legend_param = list(title = "Condition")
)

lfc_mat <- as.matrix(lfc_mat)
padj_mat <- as.matrix(padj_mat)
rownames(lfc_mat)  <- gene_label
rownames(padj_mat) <- gene_label

ht <- ComplexHeatmap::Heatmap(
  lfc_mat,
  name              = "log2FC",
  col               = col_fun,
  na_col            = "#D3D3D3",  # only used if there are real NA values
  cluster_rows      = row_clust,
  cluster_columns   = TRUE,
  top_annotation    = ko_ann,
  show_row_names    = nrow(lfc_mat) <= 100,
  border            = FALSE,
  use_raster        = TRUE,
  heatmap_legend_param = list(title = "log2FC"),
  
  # Custom drawing for each cell
  cell_fun = function(j, i, x, y, width, height, fill) {
    # i = row index, j = column index
    padj_val <- padj_mat[i, j]
    lfc_val  <- lfc_mat[i, j]

    if (!is.na(padj_val) && padj_val >= DE_PADJ) {
      # Non-significant: draw grey cell over the default
      grid::grid.rect(
        x      = x,
        y      = y,
        width  = width,
        height = height,
        gp = grid::gpar(fill = "#D3D3D3", col = NA)
      )
    }

    # 2) Add text (log2FC) on top of the cell
    if (!is.na(lfc_val) && nrow(lfc_mat) <= 50) {
      grid::grid.text(
        label = sprintf("%.2f", lfc_val),  # 2 decimal places
        x     = x,
        y     = y,
        gp    = grid::gpar(fontsize = 10)   # adjust as needed
      )
    }
    # If padj < DE_PADJ, do nothing -> the default colored cell stays visible
  }
)

png(file.path(out_dir_hm, "DE_heatmap.png"), width = 2000, height = 3500, res = 300)
ComplexHeatmap::draw(ht)
grDevices::dev.off()
# stop("Heatmap_done")