# Load required libraries
library(clusterProfiler)    # For enrichment analysis
library(org.Hs.eg.db)       # Human gene annotations
library(pathview)           # For KEGG pathway visualization
library(dplyr)              # For data manipulation

# ===============================
# Function 1: Perform KEGG enrichment
# ===============================

perform_kegg_enrichment <- function(df, gene_col = "gene", organism = "hsa") {
  # Check that gene column exists
  if (!gene_col %in% colnames(df)) {
    stop(paste("Column", gene_col, "not found in dataframe"))
  }
  
  # Convert gene symbols to Entrez IDs
  gene_symbols <- df[[gene_col]]
  gene_mapping <- bitr(gene_symbols,
                       fromType = "SYMBOL",
                       toType = "ENTREZID",
                       OrgDb = org.Hs.eg.db)
  
  # Remove duplicated mappings
  gene_mapping <- gene_mapping[!duplicated(gene_mapping$ENTREZID), ]
  
  # Run KEGG enrichment
  enrich_res <- enrichKEGG(gene = gene_mapping$ENTREZID,
                           organism = organism,
                           pvalueCutoff = 0.05)
  
  return(list(enrich = enrich_res, mapping = gene_mapping))
}

# ===============================
# Function 2: Visualize KEGG pathway
# ===============================

visualize_kegg_pathway <- function(enrich_result,
                                   df,
                                   value_col = "value",
                                   gene_col = "gene",
                                   pathway_id,
                                   out_suffix = NULL) {
  # Retrieve gene mapping
  gene_mapping <- enrich_result$mapping
  
  # Join user dataframe with gene mapping
  df <- df %>% inner_join(gene_mapping, by = setNames("SYMBOL", gene_col))
  
  # Validate value column
  if (!value_col %in% colnames(df)) {
    stop(paste("Column", value_col, "not found in dataframe"))
  }
  if (!is.numeric(df[[value_col]])) {
    stop(paste("Column", value_col, "must be numeric"))
  }
  
  # Create named vector of values
  gene_data <- df[[value_col]]
  names(gene_data) <- df$ENTREZID
  
  # Generate KEGG pathway image
  pathview(gene.data = gene_data,
           pathway.id = pathway_id,
           species = substr(pathway_id, 1, 3),
           limit = list(gene = max(abs(gene_data))),
           kegg.native = TRUE,
           out.suffix = ifelse(is.null(out_suffix),
                               paste0("colored_", pathway_id),
                               out_suffix))
}

# ======================================
# Helper: Run pathway plot only if gene overlap exists
# ======================================

safe_kegg_plot <- function(res_object, df, value_col, gene_col, pathway_id, out_suffix) {
  # Extract mapped Entrez IDs
  mapped_entrez <- res_object$mapping$ENTREZID
  
  # Get genes in the pathway
  kegg_data <- KEGGREST::keggGet(pathway_id)
  if (length(kegg_data[[1]]$GENE) == 0) {
    message("Skipping ", out_suffix, ": no gene entries in KEGG pathway ", pathway_id)
    return(invisible(NULL))
  }
  pathway_entrez <- kegg_data[[1]]$GENE[seq(1, length(kegg_data[[1]]$GENE), 2)]
  
  # Intersect
  common_genes <- intersect(mapped_entrez, pathway_entrez)
  
  if (length(common_genes) == 0) {
    message("Skipping ", out_suffix, ": no overlapping genes with ", pathway_id)
    return(invisible(NULL))
  }
  
  # Proceed with plot
  visualize_kegg_pathway(res_object,
                         df = df,
                         value_col = value_col,
                         gene_col = gene_col,
                         pathway_id = pathway_id,
                         out_suffix = out_suffix)
}


# ===============================
# Alternative splicing
#
# Synaptic signaling pathways:
#   Cholinergic: hsa04725
#   GABAergic: hsa04727
#   Serotonergic: hsa04726
#   Dopaminergic: hsa04728
#   Glutamatergic: hsa04724
# 
# Signaling pathway:
#   p53: hsa04115
# ===============================


# Prepare df for analysis
EWSR1_splicing_filtered_mv <- EWSR1_splicing_filtered %>%
  filter(MV.dPsi._at_0.95 >= 0.1)
#convert dPSI column to percetnage
EWSR1_splicing_filtered_mv$E.dPsi. <- round(EWSR1_splicing_filtered_mv$E.dPsi.,2)
FUS_splicing_filtered_mv <- FUS_splicing_filtered %>%
  filter(MV.dPsi._at_0.95 >= 0.1)
#convert dPSI column to percetnage
FUS_splicing_filtered_mv$E.dPsi. <- round(FUS_splicing_filtered_mv$E.dPsi.,2)
TAF15_splicing_filtered_mv <- TAF15_splicing_filtered %>%
  filter(MV.dPsi._at_0.95 >= 0.1)
#convert dPSI column to percetnage
TAF15_splicing_filtered_mv$E.dPsi. <- round(TAF15_splicing_filtered_mv$E.dPsi.,2)


# Named list of pathway names and KEGG identifiers
kegg_pathways <- list(
  Cholinergic   = "hsa04725",
  GABAergic     = "hsa04727",
  Serotonergic  = "hsa04726",
  Dopaminergic  = "hsa04728",
  Glutamatergic = "hsa04724",
  p53           = "hsa04115"
)

# Prepare input data for looping
samples <- list(
  EWSR1 = list(df = EWSR1_splicing_filtered_mv, res = res_EWSR1),
  FUS   = list(df = FUS_splicing_filtered_mv,   res = res_FUS),
  TAF15 = list(df = TAF15_splicing_filtered_mv, res = res_TAF15)
)

# Loop through each sample and pathway
for (sample_name in names(samples)) {
  sample_data <- samples[[sample_name]]
  
  for (pathway_name in names(kegg_pathways)) {
    pathway_id <- kegg_pathways[[pathway_name]]
    
    # Create output name like "FUS_GABAergic"
    out_suffix <- paste0(sample_name, "_", pathway_name)
    
    # Call safe plotting function
    safe_kegg_plot(res_object = sample_data$res,
                   df = sample_data$df,
                   value_col = "E.dPsi.",
                   gene_col = "GENE",
                   pathway_id = pathway_id,
                   out_suffix = out_suffix)
  }
}




