# KEGG Pathway Analysis and Visualization
#
# This script provides a set of functions to download and visualise KEGG pathway
# diagrams with user‑supplied gene‑level statistics.  It relies on Bioconductor
# packages such as `pathview` for generating the native KEGG diagrams and
# `ggkegg` combined with `ggraph` for producing more flexible graphics that can
# subsequently be made interactive via `plotly`.  The helper functions are
# designed to be general so they can be reused for multiple differential
# analyses (e.g. differential gene expression, alternative splicing, etc.).
# Where possible, gene identifiers are mapped to KEGG’s internal format using
# annotation databases (e.g. `org.Hs.eg.db`) so the functions work across
# species.
#
# Before running this script the required packages should be installed.  A
# convenient way to install Bioconductor packages is to use BiocManager.  For
# example:
#
  # if (!requireNamespace("BiocManager", quietly = TRUE)) {
  #   install.packages("BiocManager")
  # }
  # BiocManager::install(c("pathview", "ggkegg", "AnnotationDbi", "org.Hs.eg.db"))
  # install.packages(c("tidyverse", "plotly"))
#
# Load required libraries.  Suppress package startup messages to keep the
# console output clean when sourcing this script.
suppressPackageStartupMessages({
  library(pathview)      # Render native KEGG pathway diagrams【467489955321676†L38-L84】
  library(ggkegg)        # Parse KEGG pathways into tidy graph structures
  library(AnnotationDbi) # Map between gene identifier types
  library(org.Hs.eg.db)  # Annotation database for Homo sapiens; change as needed
  library(dplyr)         # Data manipulation
  library(tidyr)         # Data tidying
  library(ggraph)        # Grammar of graphics for network graphs
  library(ggplot2)       # General plotting
  library(plotly)        # Convert ggplot2 objects to interactive plots
})

################################################################################
# Gene identifier conversion
################################################################################

#' Convert gene identifiers to Entrez IDs and build a named numeric vector
#'
#' @param df A data frame containing at least two columns: one with gene
#'   identifiers (symbols, Ensembl IDs, etc.) and one with numeric values
#'   (e.g. log2 fold change).  Additional columns are ignored.
#' @param gene_col Character string naming the column in `df` containing gene
#'   identifiers.
#' @param value_col Character string naming the column in `df` containing
#'   numeric values associated with each gene.
#' @param keytype Character string specifying the type of gene identifier in
#'   `gene_col`.  Valid key types can be found via `AnnotationDbi::keytypes()`.
#' @param annotation_pkg The OrgDb package used for mapping.  For human data
#'   this defaults to `org.Hs.eg.db`, but other organisms can be specified.
#'
#' @return A named numeric vector where names are Entrez Gene IDs and values
#'   correspond to the data in `value_col`.  Genes without a mapping are
#'   silently dropped.  If multiple input identifiers map to the same Entrez
#'   identifier, their values are averaged.
convert_to_entrez_vector <- function(df,
                                     gene_col,
                                     value_col,
                                     keytype = "SYMBOL",
                                     annotation_pkg = org.Hs.eg.db) {
  # Extract the relevant columns and ensure the value column is numeric
  df_sub <- df %>%
    dplyr::select(gene = all_of(gene_col), value = all_of(value_col)) %>%
    dplyr::mutate(value = as.numeric(value)) %>%
    dplyr::filter(!is.na(value))

  # Map identifiers to Entrez Gene IDs using the specified OrgDb
  mapped <- AnnotationDbi::select(
    annotation_pkg,
    keys = unique(df_sub$gene),
    columns = c("ENTREZID"),
    keytype = keytype
  )

  # Join the mapping back to the original values
  df_mapped <- df_sub %>%
    dplyr::inner_join(mapped, by = c("gene" = keytype)) %>%
    dplyr::filter(!is.na(ENTREZID))

  # Average values if multiple identifiers map to the same Entrez ID
  summarised <- df_mapped %>%
    dplyr::group_by(ENTREZID) %>%
    dplyr::summarise(avg_value = mean(value, na.rm = TRUE), .groups = "drop")

  # Build the named numeric vector
  vec <- summarised$avg_value
  names(vec) <- summarised$ENTREZID
  return(vec)
}

################################################################################
# Pathview wrapper
################################################################################

#' Generate native KEGG diagrams coloured by gene statistics
#'
#' This function wraps the `pathview()` call to overlay user supplied gene
#' statistics onto KEGG pathway diagrams.  It supports either a single named
#' vector or a matrix with genes as rows and multiple conditions (e.g.
#' different analyses) as columns.  When multiple columns are supplied the
#' function can render all states in one figure using the multi‑state feature
#' introduced in pathview ≥1.1.6【467489955321676†L72-L74】.
#'
#' @param gene_data A named numeric vector or matrix.  When a matrix is
#'   supplied the row names must be Entrez IDs and columns correspond to
#'   different conditions or analyses.
#' @param pathway_id Character string giving the KEGG pathway identifier (e.g.
#'   "hsa04110" for human cell cycle).  IDs are case insensitive.
#' @param species Character string specifying the KEGG species code (e.g. "hsa" for
#'   human, "mmu" for mouse).  See pathview documentation for valid codes【467489955321676†L78-L105】.
#' @param out_suffix Character string appended to the output file names.
#' @param kegg_native Logical; if TRUE the native KEGG PNG is rendered.  If
#'   FALSE a Graphviz view is generated which can provide more flexibility for
#'   topology analysis【467489955321676†L83-L90】.
#' @param multi_state Logical; if TRUE and `gene_data` is a matrix, a single
#'   diagram will be created showing slices corresponding to each column
#'   (condition/analysis).  Otherwise separate diagrams will be produced for
#'   each column.
#' @param same_layer Logical; controls whether multiple states are drawn on the
#'   same layer (TRUE) or split into separate layers (FALSE).  Refer to
#'   `?pathview` for details.
#' @param low,mid,high Character vectors specifying colours for low, mid and
#'   high values respectively.  Each can be length one (applied to both genes
#'   and compounds) or named list with entries `gene` and `cpd`.
#' @param limit Numeric value or list specifying the fold change limits used to
#'   cap the colour scale.  See pathview documentation for details【467489955321676†L65-L71】.
#'
#' @return Invisibly returns the result of `pathview()`.  Files are saved to
#'   the current working directory with names incorporating `out_suffix`.
generate_kegg_map <- function(gene_data,
                              pathway_id,
                              species = "hsa",
                              out_suffix = "analysis",
                              kegg_native = TRUE,
                              multi_state = FALSE,
                              same_layer = TRUE,
                              low = list(gene = "green", cpd = "blue"),
                              mid = list(gene = "gray", cpd = "gray"),
                              high = list(gene = "red", cpd = "yellow"),
                              limit = list(gene = 2, cpd = 2)) {
  # Ensure row names are available when gene_data is a matrix
  if (is.matrix(gene_data) && is.null(rownames(gene_data))) {
    stop("gene_data matrix must have row names corresponding to Entrez IDs")
  }

  # Call pathview; wrap in tryCatch to gracefully handle errors
  pv_out <- tryCatch({
    pathview::pathview(
      gene.data   = gene_data,
      pathway.id  = pathway_id,
      species     = species,
      out.suffix  = out_suffix,
      kegg.native = kegg_native,
      multi.state = multi_state,
      same.layer  = same_layer,
      low         = low,
      mid         = mid,
      high        = high,
      limit       = limit
    )
  }, error = function(e) {
    warning(sprintf("pathview failed for pathway %s: %s", pathway_id, e$message))
    return(NULL)
  })
  invisible(pv_out)
}

################################################################################
# ggkegg interactive visualisation
################################################################################

#' Prepare an interactive KEGG pathway plot with multiple gene lists
#'
#' This function parses a KEGG pathway into a `tbl_graph` object using the
#' `ggkegg` package, merges user supplied gene statistics onto the node data and
#' produces a `ggplot2` object via `ggraph`.  The resulting plot can be turned
#' interactive with `plotly::ggplotly()`.  By default the underlying KEGG image
#' is overlaid using `overlay_raw_map()`; this can be toggled via an argument.
#'
#' @param pathway_id Character string giving the KEGG pathway identifier (e.g.
#'   "hsa04110").
#' @param gene_lists A list where each element is a data frame with at least two
#'   columns: one containing gene symbols (or the identifier specified by
#'   `keytype`) and one containing numeric scores.  Each list element should be
#'   named; these names become the category labels used to colour nodes.
#' @param gene_col, value_col Character strings naming the gene and score
#'   columns within each data frame.
#' @param keytype Character string specifying the identifier type used in
#'   `gene_col`.  For example "SYMBOL" or "ENSEMBL".  Must be a valid key type
#'   for the annotation database.
#' @param annotation_pkg OrgDb object used to map identifiers to Entrez IDs.
#' @param species Character string specifying the KEGG species code.  This is
#'   used by `ggkegg::pathway()` to fetch the pathway (defaults to "hsa").
#' @param overlay_raw Logical; whether to overlay the original KEGG PNG on the
#'   network plot.  Uses `ggkegg::overlay_raw_map()`【430613313600756†L67-L84】.
#'
#' @return A list containing two elements: `plot` (a ggplot2 object) and
#'   `interactive` (an interactive plotly object).  The interactive object is
#'   optional; it is only created if `plotly` is installed.
prepare_interactive_kegg_plot <- function(pathway_id,
                                          gene_lists,
                                          gene_col = "gene",
                                          value_col = "value",
                                          keytype = "SYMBOL",
                                          annotation_pkg = org.Hs.eg.db,
                                          species = "hsa",
                                          overlay_raw = TRUE) {
  # Validate input
  if (is.null(names(gene_lists)) || any(names(gene_lists) == "")) {
    stop("gene_lists must be a named list; each element name will be used as a category")
  }

  # Parse the KEGG pathway into a tidy graph.  use_cache=TRUE saves the KGML
  # locally for faster repeated access.
  graph <- ggkegg::pathway(pathway_id, use_cache = TRUE)

  # Extract node data and prepare a tibble for joining
  node_tbl <- graph %>%
    tidygraph::activate(nodes) %>%
    as_tibble() %>%
    # The "name" column contains KEGG identifiers, e.g. "hsa:7157"
    dplyr::mutate(entrez_id = stringr::str_replace(name, ".*:", ""))

  # Prepare a combined annotation of all gene lists
  annotation <- purrr::imap_dfr(gene_lists, function(df, category) {
    # Convert to Entrez vector
    vec <- convert_to_entrez_vector(df, gene_col, value_col, keytype, annotation_pkg)
    tibble(
      entrez_id = names(vec),
      score     = as.numeric(vec),
      category  = category
    )
  })

  # In case a gene appears in multiple categories, keep the record with the
  # largest absolute score.  This resolves overlaps between analyses by
  # favouring the most strongly affected dataset.
  annotation_unique <- annotation %>%
    dplyr::arrange(desc(abs(score))) %>%
    dplyr::distinct(entrez_id, .keep_all = TRUE)

  # Join the annotation onto the node table
  node_tbl_annot <- node_tbl %>%
    dplyr::left_join(annotation_unique, by = "entrez_id") %>%
    # Map Entrez IDs back to gene symbols for labels when available
    dplyr::mutate(
      gene_symbol = AnnotationDbi::mapIds(
        annotation_pkg, keys = entrez_id, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"
      )
    )

  # Define a colour palette for categories.  Use distinct hues for each category.
  categories <- unique(annotation_unique$category)
  pal <- scales::hue_pal()(length(categories))
  names(pal) <- categories

  # Build the ggraph plot
  p <- ggraph(graph, layout = "manual", x = node_tbl$x, y = node_tbl$y) +
    # Draw edges; arrowheads indicate direction where available
    geom_edge_parallel(aes(linetype = subtype_name),
                       arrow = arrow(length = unit(1, "mm"), type = "closed"),
                       end_cap = circle(1, "mm"), start_cap = circle(1, "mm"),
                       colour = "grey50") +
    # Draw nodes; use rectangles for genes (as per KEGG diagrams)
    geom_node_rect(
      aes(
        filter = type == "gene",
        fill = category
      ),
      colour = "black"
    ) +
    # Draw other nodes (compounds, etc.) as light grey
    geom_node_rect(
      aes(
        filter = type != "gene"
      ),
      fill = scales::alpha("lightgrey", 0.3),
      colour = "black"
    ) +
    # Add text labels for genes; show gene symbol when available
    geom_node_text(
      aes(
        label = ifelse(!is.na(gene_symbol), gene_symbol, entrez_id),
        filter = type == "gene"
      ),
      size = 2.5, vjust = 0.5
    ) +
    # Set category colours; use NA for nodes without category (transparent)
    scale_fill_manual(values = pal, na.value = NA) +
    guides(fill = guide_legend(title = "Analysis")) +
    theme_void()

  # Optionally overlay the raw KEGG map.  This uses a separate layer provided
  # by ggkegg and preserves the original diagram context【430613313600756†L67-L84】.
  if (overlay_raw) {
    p <- p + ggkegg::overlay_raw_map()
  }

  # Make interactive version using plotly.  Tooltips display gene symbol,
  # Entrez ID, category and score.  If plotly is not installed the
  # interactive component will be NULL.
  interactive_plot <- NULL
  if (requireNamespace("plotly", quietly = TRUE)) {
    # Add tooltip aesthetics; note that ggraph/ggplot objects created above
    # cannot directly specify tooltip text.  plotly picks up labels from
    # underlying data.  Here we customise it via ggplotly.
    interactive_plot <- plotly::ggplotly(
      p,
      tooltip = c("label", "fill")
    ) %>%
      plotly::style(
        hoverinfo = "text",
        text = ~paste(
          "Gene:", ifelse(is.na(node_tbl_annot$gene_symbol), node_tbl_annot$entrez_id, node_tbl_annot$gene_symbol),
          "<br>Entrez ID:", node_tbl_annot$entrez_id,
          "<br>Category:", ifelse(is.na(node_tbl_annot$category), "None", node_tbl_annot$category),
          "<br>Score:", ifelse(is.na(node_tbl_annot$score), "NA", sprintf("%.3f", node_tbl_annot$score))
        )
      )
  }

  return(list(plot = p, interactive = interactive_plot))
}

################################################################################
# Example usage
################################################################################

# The following commented example illustrates how to use the functions defined
# above.  It assumes you have multiple differential analysis results stored as
# data frames with gene identifiers and scores (e.g. log2 fold change).
# Replace `df_de`, `df_splicing`, etc. with your own data.  For demonstration
# purposes the code below will not run unless you uncomment it and supply
# suitable data.

# Example: prepare gene lists for multiple analyses
df_de       <- data.frame(gene = c("TP53", "CDKN2A", "CDC45"), value = c(2.1, -1.5, 0.8))
df_splicing <- data.frame(gene = c("TP53", "BCL2", "BRCA1"), value = c(1.2, -2.0, 0.5))
df_tf       <- data.frame(gene = c("NFkB1", "TP53"), value = c(-0.9, 1.0))

# Assemble a named list of gene lists.  Names define the category labels.
gene_lists <- list(
  DifferentialExpression = df_de,
  AlternativeSplicing   = df_splicing,
  TranscriptionFactors  = df_tf
)

# Convert gene lists into a gene × category matrix for pathview
vectors <- lapply(gene_lists, function(df) {
  convert_to_entrez_vector(df, gene_col = "gene", value_col = "value", keytype = "SYMBOL")
})

# Find the union of all Entrez IDs across analyses
all_entrez <- unique(unlist(lapply(vectors, names)))

# Build a matrix where rows are Entrez IDs and columns are analyses; missing
# values are filled with NA
gene_matrix <- sapply(vectors, function(vec) {
  v <- rep(NA_real_, length(all_entrez))
  names(v) <- all_entrez
  common <- intersect(names(vec), all_entrez)
  v[common] <- vec[common]
  return(v)
})
# Set row names for pathview
rownames(gene_matrix) <- all_entrez

# Generate a KEGG pathway map showing all analyses in one figure【467489955321676†L72-L74】
generate_kegg_map(
  gene_data   = gene_matrix,
  pathway_id  = "hsa04110",   # Replace with your KEGG pathway ID
  species     = "hsa",
  out_suffix  = "multi_analysis",
  kegg_native = TRUE,
  multi_state = TRUE,
  same_layer  = TRUE
)

# Create an interactive network plot overlaying the gene statistics
interactive_res <- prepare_interactive_kegg_plot(
  pathway_id  = "hsa04110",    # Same as above
  gene_lists  = gene_lists,
  gene_col    = "gene",
  value_col   = "value",
  keytype     = "SYMBOL",
  annotation_pkg = org.Hs.eg.db,
  species     = "hsa",
  overlay_raw = TRUE
)

# Display the static ggraph plot
print(interactive_res$plot)

# Display the interactive plot in an RMarkdown or Shiny context
if (!is.null(interactive_res$interactive)) {
  interactive_res$interactive
}

################################################################################
# End of script