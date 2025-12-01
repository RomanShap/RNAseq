# 08_go_shared_events.R
#
# Step 8: GO analysis on shared and condition‑specific splicing event groups.
#
# This script takes the event membership matrices from step 06 and defines
# exact membership groups using the helper function `build_exact_membership_groups()`
# from shared/99_utils.R.  For each group of events (e.g. events present
# exclusively in EWSR1, shared between FUS and TAF15, shared among all three,
# etc.) it extracts the corresponding genes and runs GO over‑representation
# analysis.  Results are written to the output directory `08_go_shared_events/`.

message("Step 08: GO analysis of shared and specific event groups")

## Check prerequisites ---------------------------------------------------------
if (!exists("EVENT_MEMBERSHIP", envir = .GlobalEnv)) {
  stop("EVENT_MEMBERSHIP object is missing; please run step 06 first")
}
if (!exists("FILTERED_EVENTS", envir = .GlobalEnv)) {
  stop("FILTERED_EVENTS object is missing; please run step 01 first")
}
event_membership <- get("EVENT_MEMBERSHIP", envir = .GlobalEnv)
filtered_tables <- get("FILTERED_EVENTS", envir = .GlobalEnv)

## Create output directory -----------------------------------------------------
step_output <- file.path(PROJECT_ROOT, "output", "08_go_shared_events")
if (!dir.exists(step_output)) dir.create(step_output, recursive = TRUE)

## Load required packages ------------------------------------------------------
suppressMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
})

## Function to map event IDs to gene symbols ----------------------------------
get_genes_for_events <- function(event_ids) {
  genes <- character()
  for (cond in names(filtered_tables)) {
    df <- filtered_tables[[cond]]
    gene_col <- intersect(c("GeneName", "GENEID", "GENE"), names(df))
    if (length(gene_col) > 0) {
      idx <- match(event_ids, df$EVENT)
      sel <- !is.na(idx)
      genes <- c(genes, df[[gene_col[1]]][idx[sel]])
    }
  }
  unique(genes)
}

## Helper for GO enrichment ----------------------------------------------------
run_go_enrichment <- function(gene_symbols) {
  entrez <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = gene_symbols,
                                  keytype = "SYMBOL", column = "ENTREZID",
                                  multiVals = "first")
  entrez <- entrez[!is.na(entrez)]
  if (length(entrez) == 0) return(NULL)
  res <- tryCatch({
    enrichGO(gene = entrez,
             OrgDb = org.Hs.eg.db,
             ont = cfg$go$ontology,
             keyType = "ENTREZID",
             pAdjustMethod = "BH",
             pvalueCutoff = 0.05,
             qvalueCutoff = 0.2)
  }, error = function(e) {
    warning(sprintf("enrichGO failed: %s", e$message))
    NULL
  })
  if (!is.null(res) && nrow(res) > 0) {
    res <- clusterProfiler::simplify(res, cutoff = cfg$go$simplify_cutoff)
    res <- res[res$Count >= cfg$go$count_min & res$Count <= cfg$go$count_max, ]
  }
  return(res)
}

## Iterate over event types and directions ------------------------------------
go_shared_results <- list()
for (etype in names(event_membership)) {
  go_shared_results[[etype]] <- list()
  for (dirn in names(event_membership[[etype]])) {
    mat <- event_membership[[etype]][[dirn]]
    # Skip empty matrices
    if (is.null(mat) || nrow(mat) == 0) next
    # Remove EVENT column for membership computation
    membership_matrix <- mat[, -1, drop = FALSE]
    rownames(membership_matrix) <- mat$EVENT
    groups <- build_exact_membership_groups(membership_matrix)
    # groups is a list mapping group names (e.g. "EWSR1&FUS") to event IDs
    group_res_list <- list()
    for (grp_name in names(groups)) {
      event_ids <- groups[[grp_name]]
      genes <- get_genes_for_events(event_ids)
      if (length(genes) == 0) next
      go_res <- run_go_enrichment(genes)
      group_res_list[[grp_name]] <- go_res
      # Save to file
      if (!is.null(go_res) && nrow(go_res) > 0) {
        timestamp <- format(Sys.time(), "%Y%m%d%H%M%S")
        out_file <- file.path(step_output, sprintf("GO_shared_%s_%s_%s_%s.tsv", etype, dirn, grp_name, timestamp))
        readr::write_tsv(as.data.frame(go_res), out_file)
      }
    }
    go_shared_results[[etype]][[dirn]] <- group_res_list
  }
}

## Save GO shared results ------------------------------------------------------
go_shared_rds <- file.path(step_output, "go_shared_results.rds")
saveRDS(go_shared_results, go_shared_rds)
message(sprintf("GO results for shared/specific event groups saved to %s", go_shared_rds))

## Return shared results to global environment --------------------------------
assign("GO_SHARED_RESULTS", go_shared_results, envir = .GlobalEnv)
message("Step 08 completed")