############################################################
## KEGG enrichment + pathway plotting + exports (R script) ##
############################################################

## -----------------------------
## 0) Required libraries
## -----------------------------
suppressPackageStartupMessages({
  library(dplyr)            # Data wrangling
  library(clusterProfiler)  # Enrichment analysis
  library(org.Hs.eg.db)     # Human gene annotations (Entrez/SYMBOL)
  library(pathview)         # KEGG pathway visualization
  library(KEGGREST)         # Access KEGG content (e.g., genes per pathway)
})

## ---------------------------------------------------------
## 1) Core function: run KEGG enrichment on a data frame
## ---------------------------------------------------------
#' Perform KEGG pathway enrichment using clusterProfiler.
#'
#' @param df Data frame containing a gene identifier column.
#' @param gene_col Column name in `df` holding gene symbols (default) or other IDs.
#'                 For this script, enrichment is done on Entrez IDs, so symbols
#'                 are mapped to Entrez via org.Hs.eg.db.
#' @param organism KEGG organism code (e.g., "hsa" for human).
#' @param pvalueCutoff Nominal p-value cutoff used by enrichKEGG.
#' @param qvalueCutoff FDR cutoff used by enrichKEGG.
#' @return A named list with:
#'         - enrich: enrichResult object (clusterProfiler)
#'         - mapping: data.frame with SYMBOL and ENTREZID used
#'         - params: list of parameters used for reproducibility
perform_kegg_enrichment <- function(df,
                                    gene_col = "gene",
                                    organism = "hsa",
                                    pvalueCutoff = 0.05,
                                    qvalueCutoff = 0.25) {
  # Basic checks
  if (!gene_col %in% colnames(df)) {
    stop(sprintf("Column '%s' not found in the input data frame.", gene_col))
  }
  
  # Extract symbols as character vector (remove missing)
  gene_symbols <- as.character(df[[gene_col]])
  gene_symbols <- gene_symbols[!is.na(gene_symbols) & nzchar(gene_symbols)]
  
  if (length(gene_symbols) == 0L) {
    stop("No non-missing gene symbols found to map.")
  }
  
  # Map SYMBOL -> ENTREZ using org.Hs.eg.db
  # Duplicated Entrez IDs are pruned to one entry (required by enrichKEGG)
  mapping <- bitr(gene_symbols,
                  fromType = "SYMBOL",
                  toType   = "ENTREZID",
                  OrgDb    = org.Hs.eg.db)
  mapping <- mapping[!duplicated(mapping$ENTREZID), , drop = FALSE]
  
  if (!nrow(mapping)) {
    stop("No SYMBOLs could be mapped to ENTREZ IDs. Check 'gene_col' and species.")
  }
  
  # Run KEGG enrichment
  enrich_res <- enrichKEGG(
    gene          = mapping$ENTREZID,
    organism      = organism,
    pvalueCutoff  = pvalueCutoff,
    qvalueCutoff  = qvalueCutoff
  )
  
  # Return a compact bundle
  return(list(
    enrich  = enrich_res,
    mapping = mapping,
    params  = list(
      gene_col      = gene_col,
      organism      = organism,
      pvalueCutoff  = pvalueCutoff,
      qvalueCutoff  = qvalueCutoff
    )
  ))
}

## ----------------------------------------------------------------
## 2) Plot a KEGG pathway with values mapped to Entrez gene IDs
## ----------------------------------------------------------------
#' Visualize a KEGG pathway with numeric values overlaid (via pathview).
#'
#' @param enrich_result Output list from perform_kegg_enrichment().
#' @param df Original data frame that contains the quantitative column.
#' @param value_col Column in `df` containing numeric values to color the pathway.
#' @param gene_col Column in `df` containing gene symbols to join with mapping.
#' @param pathway_id KEGG pathway identifier, e.g., "hsa04115".
#' @param out_suffix Suffix for output files produced by pathview.
#' @details This function:
#'   - joins `df` with the SYMBOL->ENTREZ mapping,
#'   - builds a named numeric vector keyed by Entrez IDs,
#'   - calls pathview to generate KEGG-native PNGs.
visualize_kegg_pathway <- function(enrich_result,
                                   df,
                                   value_col = "value",
                                   gene_col  = "gene",
                                   pathway_id,
                                   out_suffix = NULL) {
  # Validate enrichment container
  if (!all(c("mapping", "enrich") %in% names(enrich_result))) {
    stop("enrich_result must contain 'mapping' and 'enrich'.")
  }
  
  mapping <- enrich_result$mapping
  
  # Join user df (symbols) to mapping (SYMBOL, ENTREZID)
  if (!gene_col %in% colnames(df)) {
    stop(sprintf("Column '%s' not found in df.", gene_col))
  }
  
  df_join <- df %>%
    inner_join(mapping, by = setNames("SYMBOL", gene_col))
  
  if (!value_col %in% colnames(df_join)) {
    stop(sprintf("Column '%s' not found after join.", value_col))
  }
  if (!is.numeric(df_join[[value_col]])) {
    stop(sprintf("Column '%s' must be numeric.", value_col))
  }
  
  # Build named vector: names are ENTREZ IDs, values are numeric scores
  gene_data <- df_join[[value_col]]
  names(gene_data) <- df_join$ENTREZID
  
  # Scale limit: symmetric range around zero if mixed signs; otherwise max absolute
  lim <- max(abs(gene_data), na.rm = TRUE)
  
  # pathview call: species inferred from the first three characters of pathway_id
  pathview(
    gene.data    = gene_data,
    pathway.id   = pathway_id,
    species      = substr(pathway_id, 1, 3),
    limit        = list(gene = lim),
    kegg.native  = TRUE,
    out.suffix   = ifelse(is.null(out_suffix),
                          paste0("colored_", pathway_id),
                          out_suffix)
  )
  
  invisible(NULL)
}

## -------------------------------------------------------------------
## 3) Safety wrapper: only plot if there is actual gene-pathway overlap
## -------------------------------------------------------------------
#' Plot a pathway only if any mapped Entrez IDs are present in the pathway.
#'
#' @param res_object Output of perform_kegg_enrichment().
#' @param df Data frame with gene and value columns.
#' @param value_col Column with numeric values.
#' @param gene_col Column with gene symbols to join with mapping.
#' @param pathway_id KEGG pathway ID (e.g., "hsa04724").
#' @param out_suffix Output suffix for produced files.
safe_kegg_plot <- function(res_object,
                           df,
                           value_col,
                           gene_col,
                           pathway_id,
                           out_suffix) {
  # Extract mapped Entrez set
  mapped_entrez <- unique(res_object$mapping$ENTREZID)
  
  # Query KEGG for pathway gene members
  kdat <- KEGGREST::keggGet(pathway_id)
  if (length(kdat) == 0L || is.null(kdat[[1]]$GENE)) {
    message("Skipping ", out_suffix, ": no gene entries found in KEGG for ", pathway_id)
    return(invisible(NULL))
  }
  
  # KEGGREST returns the $GENE element for a pathway as a character vector
  # alternating between:
  #   1) Entrez Gene IDs (e.g., "7157")
  #   2) Gene symbol and description (e.g., "TP53, tumor protein p53")
  # To extract only the Entrez IDs, take the odd-numbered elements:
  #   positions 1, 3, 5, ... (seq(1, length(...), 2)).
  pathway_entrez <- kdat[[1]]$GENE[seq(1, length(kdat[[1]]$GENE), 2)]
  
  # Intersect with mapped Entrez IDs
  common <- intersect(mapped_entrez, pathway_entrez)
  
  if (length(common) == 0L) {
    message("Skipping ", out_suffix, ": no overlapping genes with ", pathway_id)
    return(invisible(NULL))
  }
  
  # Proceed to draw
  visualize_kegg_pathway(
    enrich_result = res_object,
    df            = df,
    value_col     = value_col,
    gene_col      = gene_col,
    pathway_id    = pathway_id,
    out_suffix    = out_suffix
  )
}

## --------------------------------------------------------------------
## 4) Export utilities: save enrichment results for later comparison
## --------------------------------------------------------------------
#' Save KEGG enrichment outputs to disk (CSV, filtered CSV, RDS, mapping, metadata).
#'
#' This utility serializes the results of perform_kegg_enrichment() so they can be
#' compared across samples and reloaded reproducibly. It writes:
#' - a full enrichment table (*_KEGG_enrichment.csv);
#' - an optional filtered table by FDR (*_KEGG_enrichment_padj_le_<thr>.csv);
#' - the raw enrichResult object (*.rds) for lossless reuse;
#' - the SYMBOL↔ENTREZ mapping used (*_SYMBOL_to_ENTREZ_mapping.csv);
#' - a small metadata/provenance file with session info (*_metadata.txt).
#'
#' Filenames are made filesystem-safe (no '<' or '>') by encoding the FDR threshold
#' as padj_le_<value> with dots replaced by underscores.
#'
#' @param result_list Output list from perform_kegg_enrichment(); must contain
#' elements enrich (an enrichResult) and mapping (data.frame with SYMBOL/ENTREZID).
#' @param sample_name Character label for the sample/contrast (e.g., "EWSR1").
#' @param out_dir Directory where files are written; created if it does not exist.
#' @param p_adjust_filter Numeric FDR (adjusted p-value) threshold used to produce
#' a concise, human-readable filtered CSV; set NULL to skip the filtered export.
#'
#' @return (invisible) named list of written file paths.
save_kegg_enrichment <- function(result_list,
                                 sample_name,
                                 out_dir = "kegg_enrichment_exports",
                                 p_adjust_filter = 0.05) {
  if (is.null(result_list) || !is.list(result_list) ||
      !all(c("enrich", "mapping") %in% names(result_list))) {
    stop("result_list must be the output of perform_kegg_enrichment().")
  }
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # Convert enrichResult to plain data.frame for CSV
  enrich_df <- as.data.frame(result_list$enrich)
  
  # Filtered table for quick human comparison (optional)
  if (!is.null(p_adjust_filter) && "p.adjust" %in% colnames(enrich_df)) {
    enrich_df_filtered <- enrich_df[enrich_df$p.adjust <= p_adjust_filter, , drop = FALSE]
  } else {
    enrich_df_filtered <- enrich_df
  }
  
  # Create a filesystem-safe tag for the p.adjust threshold (e.g., "padj_le_0_05")
  padj_tag <- if (!is.null(p_adjust_filter)) {
    paste0("padj_le_", gsub("\\.", "_", as.character(p_adjust_filter)))
  } else {
    NULL
  }
  
  # Paths (avoid '<' and '>')
  csv_full_path  <- file.path(out_dir, paste0(sample_name, "_KEGG_enrichment.csv"))
  csv_filt_path  <- if (!is.null(padj_tag)) {
    file.path(out_dir, paste0(sample_name, "_KEGG_enrichment_", padj_tag, ".csv"))
  } else NA_character_
  rds_path       <- file.path(out_dir, paste0(sample_name, "_KEGG_enrichment.rds"))
  map_path       <- file.path(out_dir, paste0(sample_name, "_SYMBOL_to_ENTREZ_mapping.csv"))
  meta_path      <- file.path(out_dir, paste0(sample_name, "_KEGG_enrichment_metadata.txt"))
  
  # Write files
  utils::write.csv(enrich_df, file = csv_full_path, row.names = FALSE)
  if (!is.na(csv_filt_path)) {
    utils::write.csv(enrich_df_filtered, file = csv_filt_path, row.names = FALSE)
  }
  saveRDS(result_list$enrich, file = rds_path)
  utils::write.csv(result_list$mapping, file = map_path, row.names = FALSE)
  
  # Metadata (basic provenance + session info)
  org   <- if (!is.null(result_list$params$organism))     result_list$params$organism    else "unknown"
  pcut  <- if (!is.null(result_list$params$pvalueCutoff)) result_list$params$pvalueCutoff else NA
  qcut  <- if (!is.null(result_list$params$qvalueCutoff)) result_list$params$qvalueCutoff else NA
  
  meta_lines <- c(
    paste0("sample_name: ", sample_name),
    paste0("export_time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("organism: ", org),
    paste0("pvalueCutoff: ", pcut),
    paste0("qvalueCutoff: ", qcut),
    paste0("n_rows_enrichment: ", nrow(enrich_df)),
    paste0("n_rows_enrichment_padj_filtered: ",
           if (is.null(p_adjust_filter)) NA_integer_ else nrow(enrich_df_filtered)),
    "sessionInfo:",
    capture.output(utils::sessionInfo())
  )
  writeLines(meta_lines, con = meta_path)
  
  invisible(list(
    csv_full     = csv_full_path,
    csv_filtered = csv_filt_path,
    rds          = rds_path,
    mapping_csv  = map_path,
    metadata_txt = meta_path
  ))
}

# Helper: infix to make metadata writing cleaner
`%||%` <- function(a, b) if (!is.null(a)) a else b

## ----------------------------------------------------------------------------
## 5) Optional: build a single wide comparison table across multiple samples
## ----------------------------------------------------------------------------
#' Combine per-sample filtered KEGG CSVs into one wide comparison table.
#' The table places p.adjust and Count side by side for each sample.
#'
#' @param sample_names Character vector of sample labels used in save_kegg_enrichment().
#' @param in_dir Directory where per-sample CSVs were written.
#' @param p_adjust_filter The same threshold used when exporting filtered CSVs.
#' @param outfile Output CSV path for the combined table.
#' @return (invisible) path to the written combined CSV.
combine_kegg_for_comparison <- function(sample_names,
                                        in_dir = "kegg_enrichment_exports",
                                        p_adjust_filter = 0.05,
                                        outfile = file.path(in_dir, "KEGG_comparison_combined.csv")) {
  # Same filesystem-safe tag used in save_kegg_enrichment()
  padj_tag <- if (!is.null(p_adjust_filter)) {
    paste0("padj_le_", gsub("\\.", "_", as.character(p_adjust_filter)))
  } else {
    stop("p_adjust_filter must be specified to locate filtered CSVs.")
  }
  
  # Read filtered CSVs and tag with Sample column
  dfs <- lapply(sample_names, function(s) {
    path <- file.path(in_dir, paste0(s, "_KEGG_enrichment_", padj_tag, ".csv"))
    if (!file.exists(path)) {
      warning("File not found: ", path)
      return(NULL)
    }
    df <- utils::read.csv(path, stringsAsFactors = FALSE)
    if (nrow(df) == 0L) return(NULL)
    df$Sample <- s
    df
  })
  dfs <- Filter(Negate(is.null), dfs)
  if (!length(dfs)) stop("No input filtered tables found to combine.")
  
  keep_cols <- Reduce(intersect, lapply(dfs, colnames))
  wanted <- intersect(
    c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Sample"),
    keep_cols
  )
  df_all <- do.call(rbind, lapply(dfs, function(x) x[, wanted, drop = FALSE]))
  
  suppressPackageStartupMessages(library(reshape2))
  padj_wide  <- reshape2::dcast(df_all, ID + Description ~ Sample, value.var = "p.adjust")
  count_wide <- reshape2::dcast(df_all, ID + Description ~ Sample, value.var = "Count",
                                fun.aggregate = sum, fill = 0)
  
  merged <- merge(padj_wide, count_wide, by = c("ID","Description"),
                  suffixes = c("_padj","_count"))
  
  utils::write.csv(merged, file = outfile, row.names = FALSE)
  invisible(outfile)
}


## -------------------------------------------------------------
## 6) Example pipeline with your specific objects and columns
## -------------------------------------------------------------
## Assumptions (rename if your objects differ):
## - EWSR1_splicing_filtered, FUS_splicing_filtered, TAF15_splicing_filtered
## - Column "MV.dPsi._at_0.95" exists for filtering by |ΔPSI| support
## - Column "E.dPsi." holds the ΔPSI values (numeric)
## - Column "GENE" holds HGNC gene symbols

## 6.1) Prepare per-sample filtered data frames
EWSR1_splicing_filtered_mv <- EWSR1_splicing_filtered %>%
  filter(MV.dPsi._at_0.95 >= 0.1) %>%
  mutate(E.dPsi. = round(E.dPsi., 2))

FUS_splicing_filtered_mv <- FUS_splicing_filtered %>%
  filter(MV.dPsi._at_0.95 >= 0.1) %>%
  mutate(E.dPsi. = round(E.dPsi., 2))

TAF15_splicing_filtered_mv <- TAF15_splicing_filtered %>%
  filter(MV.dPsi._at_0.95 >= 0.1) %>%
  mutate(E.dPsi. = round(E.dPsi., 2))

## 6.2) Perform KEGG enrichment per sample (on SYMBOLs in column "GENE")
##      Note: organism "hsa" for human; adjust p/q cutoffs if required.
res_EWSR1 <- perform_kegg_enrichment(EWSR1_splicing_filtered_mv, gene_col = "GENE", organism = "hsa")
res_FUS   <- perform_kegg_enrichment(FUS_splicing_filtered_mv,   gene_col = "GENE", organism = "hsa")
res_TAF15 <- perform_kegg_enrichment(TAF15_splicing_filtered_mv, gene_col = "GENE", organism = "hsa")

## 6.3) Define KEGG pathways to draw (IDs must be valid KEGG pathway IDs)
kegg_pathways <- list(
  Cholinergic   = "hsa04725",
  GABAergic     = "hsa04727",
  Serotonergic  = "hsa04726",
  Dopaminergic  = "hsa04728",
  Glutamatergic = "hsa04724",
  p53           = "hsa04115",
  # Newly detected/shared pathways (labels kept as provided)
  Regulation_of_actin_cytoskeleton           = "hsa04820",
  Protein_digestion_and_absorption           = "hsa04974",
  ECM_receptor_interaction                   = "hsa04512",
  Regulation_of_circadian_rhythm             = "hsa04810",
  MAPK_signaling_pathway                     = "hsa04082",
  Hypertrophic_cardiomyopathy                = "hsa05410",
  PI3K_Akt_signaling_pathway                 = "hsa04151",
  Adrenergic_signaling_in_cardiomyocytes     = "hsa04261",
  Focal_adhesion                             = "hsa04510",
  Arrhythmogenic_right_ventricular_cardiomyopathy = "hsa05412",
  Gap_junction                               = "hsa04024",
  Olfactory_transduction                     = "hsa04713",
  Neuroactive_ligand_receptor_interaction    = "hsa04080",
  Dilated_cardiomyopathy                     = "hsa05414",
  Calcium_signaling_pathway                  = "hsa04020",
  Gap_junction_adhesion                      = "hsa04540"
)

## 6.4) Set up a named list for looping over samples
samples <- list(
  EWSR1 = list(df = EWSR1_splicing_filtered_mv, res = res_EWSR1),
  FUS   = list(df = FUS_splicing_filtered_mv,   res = res_FUS),
  TAF15 = list(df = TAF15_splicing_filtered_mv, res = res_TAF15)
)

## 6.5) For each sample × pathway, draw pathway only if overlapping genes exist
for (sample_name in names(samples)) {
  sdat <- samples[[sample_name]]
  for (pathway_name in names(kegg_pathways)) {
    pid <- kegg_pathways[[pathway_name]]
    out_suffix <- paste0(sample_name, "_", pathway_name)
    safe_kegg_plot(
      res_object = sdat$res,
      df         = sdat$df,
      value_col  = "E.dPsi.",
      gene_col   = "GENE",
      pathway_id = pid,
      out_suffix = out_suffix
    )
  }
}

## 6.6) Persist enrichment outputs for later comparisons
dir_out <- "kegg_enrichment_exports"
save_kegg_enrichment(res_EWSR1, sample_name = "EWSR1", out_dir = dir_out, p_adjust_filter = 0.05)
save_kegg_enrichment(res_FUS,   sample_name = "FUS",   out_dir = dir_out, p_adjust_filter = 0.05)
save_kegg_enrichment(res_TAF15, sample_name = "TAF15", out_dir = dir_out, p_adjust_filter = 0.05)

## 6.7) Build a single comparison table across the three samples
combine_kegg_for_comparison(
  sample_names    = c("EWSR1", "FUS", "TAF15"),
  in_dir          = dir_out,
  p_adjust_filter = 0.05,
  outfile         = file.path(dir_out, "KEGG_comparison_combined.csv")
)

## -----------------------------
## End of script
## -----------------------------
