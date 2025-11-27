## scripts/01_kallisto_to_counts.R
## Read Kallisto 'abundance.tsv' files and aggregate to gene-level counts.

if (!exists("kallisto_paths")) {
  stop("Object 'kallisto_paths' not found. It should be created in Quantification_main.R.")
}
if (!exists("SAMPLE_INFO")) {
  stop("Object 'SAMPLE_INFO' not found.")
}
if (!exists("output_dir")) {
  stop("Object 'output_dir' not found.")
}
if (!exists("tx2gene_rds")) {
  stop("Object 'tx2gene_rds' not found.")
}

out_dir_counts <- ensure_dir(output_dir, "01_counts")

read_kallisto_counts <- function(path, sample_name) {
  df <- readr::read_tsv(path, col_types = readr::cols())
  if (!all(c("target_id", "est_counts") %in% names(df))) {
    stop("File ", path, " does not contain 'target_id' and 'est_counts' columns.")
  }
  df <- df[, c("target_id", "est_counts")]
  colnames(df)[2] <- sample_name
  df
}

## Read transcript-level counts for all samples
message("Reading Kallisto abundance.tsv files...")
count_list <- mapply(
  read_kallisto_counts,
  path        = kallisto_paths,
  sample_name = names(kallisto_paths),
  SIMPLIFY    = FALSE
)

count_data_ENST <- Reduce(function(x, y) merge(x, y, by = "target_id"), count_list)
message("Number of transcripts read: ", nrow(count_data_ENST))

## Clean transcript IDs (remove version suffix)
count_data_ENST$target_id <- sub("\\..*$", "", count_data_ENST$target_id)
transcript_ids <- unique(count_data_ENST$target_id)

## Load or build tx2gene mapping
tx2gene <- load_or_build_tx2gene(tx2gene_rds, transcript_ids)

## Map transcripts to genes
count_data_ENST$ENSG <- tx2gene$ensembl_gene_id[
  match(count_data_ENST$target_id, tx2gene$ensembl_transcript_id)
]

count_data_ENST <- count_data_ENST[!is.na(count_data_ENST$ENSG), ]

## Aggregate transcript counts to gene-level
count_data_gene <- count_data_ENST |>
  dplyr::select(-target_id) |>
  dplyr::group_by(ENSG) |>
  dplyr::summarise(dplyr::across(where(is.numeric), sum), .groups = "drop")

count_data_gene <- as.data.frame(count_data_gene)
rownames(count_data_gene) <- count_data_gene$ENSG
count_data_gene$ENSG <- NULL

## Reorder columns to follow SAMPLE_INFO
common_samples <- intersect(colnames(count_data_gene), rownames(SAMPLE_INFO))
count_data_gene <- count_data_gene[, common_samples, drop = FALSE]

message("Number of genes after aggregation: ", nrow(count_data_gene))
## Save to disk and expose for downstream scripts
count_rds_path <- file.path(out_dir_counts, "gene_counts.rds")
saveRDS(count_data_gene, count_rds_path)
message("Saved gene-level counts to: ", count_rds_path)
COUNT_DATA_GENE <- count_data_gene
rm(count_data_gene, count_list, count_data_ENST, tx2gene)
