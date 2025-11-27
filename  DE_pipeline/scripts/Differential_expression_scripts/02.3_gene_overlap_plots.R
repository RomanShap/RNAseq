## scripts/06_gene_overlap_plots.R
## Gene-level overlaps: UpSet + Venn + export of membership.

if (!exists("DE_RESULTS"))  stop("'DE_RESULTS' not found.")
if (!exists("DE_LFC"))      stop("'DE_LFC' not found.")
if (!exists("DE_PADJ"))     stop("'DE_PADJ' not found.")
if (!exists("output_dir"))  stop("'output_dir' not found.")

out_dir_gene_ov <- ensure_dir(output_dir, "07_gene_overlap")

conds <- names(DE_RESULTS)

## Build gene sets (labels = symbol if present, else ensembl_id)
gene_sets_all <- list()
gene_sets_pos <- list()
gene_sets_neg <- list()

for (cond in conds) {
  res_df <- DE_RESULTS[[cond]]
  label  <- gene_label_from_res(res_df, id_col = "ensembl_id", sym_col = "symbol")
  
  idx_all <- which(res_df$padj < DE_PADJ &
                     abs(res_df$log2FoldChange) > DE_LFC)
  idx_pos <- which(res_df$padj < DE_PADJ &
                     res_df$log2FoldChange > DE_LFC)
  idx_neg <- which(res_df$padj < DE_PADJ &
                     res_df$log2FoldChange < -DE_LFC)
  
  gene_sets_all[[cond]] <- unique(label[idx_all])
  gene_sets_pos[[cond]] <- unique(label[idx_pos])
  gene_sets_neg[[cond]] <- unique(label[idx_neg])
}

## Binary matrices
bm_all <- build_binary_matrix_from_sets(gene_sets_all)
bm_pos <- build_binary_matrix_from_sets(gene_sets_pos)
bm_neg <- build_binary_matrix_from_sets(gene_sets_neg)

## Export membership table (item + membership pattern)
membership_table <- function(bm) {
  pat <- apply(bm, 1L, function(x) {
    conds_present <- colnames(bm)[as.logical(x)]
    if (!length(conds_present)) "none" else paste(conds_present, collapse = " & ")
  })
  data.frame(
    gene      = rownames(bm),
    membership = pat,
    bm,
    stringsAsFactors = FALSE
  )
}

openxlsx::write.xlsx(
  membership_table(bm_all),
  file.path(out_dir_gene_ov, "Gene_overlap_ALL.xlsx"),
  rowNames = FALSE
)
openxlsx::write.xlsx(
  membership_table(bm_pos),
  file.path(out_dir_gene_ov, "Gene_overlap_POS.xlsx"),
  rowNames = FALSE
)
openxlsx::write.xlsx(
  membership_table(bm_neg),
  file.path(out_dir_gene_ov, "Gene_overlap_NEG.xlsx"),
  rowNames = FALSE
)

# ## UpSetR plots (all / pos / neg)

ko_cols <- NULL
if (exists("COLOR_KO")) {
  ko_cols <- COLOR_KO   # list or named vector
}

plot_upset_complex(
  bm         = bm_all,
  outfile    = file.path(out_dir_gene_ov, "UpSet_genes_all.png"),
  title      = "DEGs (all)",
  set_colors = ko_cols
)

plot_upset_complex(
  bm         = bm_pos,
  outfile    = file.path(out_dir_gene_ov, "UpSet_genes_pos.png"),
  title      = "DEGs (upregulated)",
  set_colors = ko_cols
)

plot_upset_complex(
  bm         = bm_neg,
  outfile    = file.path(out_dir_gene_ov, "UpSet_genes_neg.png"),
  title      = "DEGs (downregulated)",
  set_colors = ko_cols
)



## Venn diagrams (only if <= 4 contrasts)
if (length(conds) <= 4L) {
  fills <- c("#d53031", "#009E73", "#0072B2", "#FFB000")[seq_along(conds)]
  names(fills) <- conds
  
  save_venn(
    sets    = gene_sets_all,
    outfile = file.path(out_dir_gene_ov, "Venn_genes_all.png"),
    fills   = fills
  )
  save_venn(
    sets    = gene_sets_pos,
    outfile = file.path(out_dir_gene_ov, "Venn_genes_pos.png"),
    fills   = fills
  )
  save_venn(
    sets    = gene_sets_neg,
    outfile = file.path(out_dir_gene_ov, "Venn_genes_neg.png"),
    fills   = fills
  )
}

## ----------------------------------------------------------------------
## Save per-condition gene sets (for reuse by other scripts)
## ----------------------------------------------------------------------

saveRDS(
  gene_sets_all,
  file.path(out_dir_gene_ov, "gene_sets_all_by_condition.rds")
)
saveRDS(
  gene_sets_pos,
  file.path(out_dir_gene_ov, "gene_sets_pos_by_condition.rds")
)
saveRDS(
  gene_sets_neg,
  file.path(out_dir_gene_ov, "gene_sets_neg_by_condition.rds")
)

# stop("Gene overlap analysis completed successfully.")