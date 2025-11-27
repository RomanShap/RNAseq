make_gene_table <- function(partitioned_events_df,
                            event_type   = "EX",
                            ko_order     = c("EWSR1","FUS","TAF15"),
                            desired_cols = c("EWSR1","EWSR1_FUS","EWSR1_FUS_TAF15",
                                             "FUS","EWSR1_TAF15","FUS_TAF15","TAF15"),
                            out_file     = NULL,                # TSV (optional)
                            out_xlsx     = "output/genes_by_combo_EX.xlsx",  # XLSX (optional)
                            sort_genes   = TRUE) {
  stopifnot(event_type %in% names(partitioned_events_df))
  lst <- partitioned_events_df[[event_type]]
  
  norm_one <- function(nm) {
    parts <- unlist(strsplit(nm, "_", fixed = TRUE))
    parts <- intersect(ko_order, parts)
    if (length(parts) == 0L) return(nm)
    paste(parts, collapse = "_")
  }
  names(lst) <- vapply(names(lst), norm_one, character(1))
  
  gene_lists <- lapply(lst, function(tb) {
    v <- unique(tb$GENE)
    v <- v[!is.na(v) & nzchar(v)]
    if (sort_genes) sort(v) else v
  })
  
  missing_cols <- setdiff(desired_cols, names(gene_lists))
  for (mc in missing_cols) gene_lists[[mc]] <- character(0)
  gene_lists <- gene_lists[desired_cols]
  
  max_len <- max(vapply(gene_lists, length, integer(1), USE.NAMES = FALSE), 0L)
  padded  <- lapply(gene_lists, function(v) { length(v) <- max_len; v })
  out     <- as.data.frame(padded, stringsAsFactors = FALSE, check.names = FALSE)
  
  if (!is.null(out_file)) {
    dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
    write.table(out, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)
  }
  if (!is.null(out_xlsx)) {
    dir.create(dirname(out_xlsx), showWarnings = FALSE, recursive = TRUE)
    if (requireNamespace("openxlsx", quietly = TRUE)) {
      openxlsx::write.xlsx(out, file = out_xlsx, asTable = TRUE, overwrite = TRUE)
    } else if (requireNamespace("writexl", quietly = TRUE)) {
      writexl::write_xlsx(out, path = out_xlsx)
    } else {
      stop("Install either 'openxlsx' or 'writexl' to write .xlsx files.")
    }
  }
  
  out
}

# Example
out_ex <- make_gene_table(
  partitioned_events_df,
  event_type = "ALTD",
  out_xlsx   = "output/genes_by_combo_ALTD.xlsx"
)
