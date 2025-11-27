## scripts/02_deseq_bimodal_filter.R
## DESeq2 + bimodal expression filter for all conditions vs reference.

if (!exists("COUNT_DATA_GENE")) stop("'COUNT_DATA_GENE' not found.")
if (!exists("SAMPLE_INFO"))   stop("'SAMPLE_INFO' not found.")
if (!exists("DE_LFC"))        stop("'DE_LFC' not found.")
if (!exists("DE_PADJ"))       stop("'DE_PADJ' not found.")
if (!exists("output_dir"))    stop("'output_dir' not found.")

out_dir_de <- ensure_dir(output_dir, "02_deseq")

## Choose reference condition
all_conditions <- unique(SAMPLE_INFO$condition)

ref_level <- if (!is.null(cfg$de$reference)) {
  cfg$de$reference
} else if ("NT" %in% all_conditions) {
  "NT"
} else if ("WT" %in% all_conditions) {
  "WT"
} else {
  all_conditions[1]
}

if (!ref_level %in% all_conditions) {
  stop("Reference level '", ref_level, "' not found in SAMPLE_INFO$condition.")
}

message("DESeq2 reference condition: ", ref_level)
## Run DESeq2 with bimodal filter
de_obj <- run_deseq_bimodal(
  counts         = COUNT_DATA_GENE,
  sample_info    = SAMPLE_INFO,
  ref_level      = ref_level,
  padj_threshold = DE_PADJ
)
# print(head(de_obj))
DESEQ_DDS <- de_obj$dds
DE_RESULTS <- de_obj$results
FILTERED_COUNTS <- de_obj$filtered_counts

## Save per-contrast results
for (cond in names(DE_RESULTS)) {
  res_df <- DE_RESULTS[[cond]]
  fname_xlsx <- file.path(out_dir_de, paste0("DE_", cond, "_vs_", ref_level, ".xlsx"))
  fname_rds  <- file.path(out_dir_de, paste0("DE_", cond, "_vs_", ref_level, ".rds"))
  openxlsx::write.xlsx(res_df, fname_xlsx, rowNames = FALSE)
  saveRDS(res_df, fname_rds)
  message("Saved DE result for ", cond, " vs ", ref_level, " to: ", fname_xlsx)
}

## Save combined list for downstream scripts
saveRDS(DE_RESULTS, file.path(out_dir_de, "DE_results_list.rds"))

## Save filtered counts matrix
saveRDS(FILTERED_COUNTS, file.path(out_dir_de, "filtered_counts_matrix.rds"))
fname_xlsx <- file.path(out_dir_de, "filtered_counts_matrix.xlsx")
openxlsx::write.xlsx(FILTERED_COUNTS, fname_xlsx, rowNames = TRUE)
message("Saved filtered counts matrix to: ", fname_xlsx)
# stop("Finished 02_deseq")