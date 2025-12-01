# 05_go_fgsea.R
#
# Step 5: GO over‑representation and FGSEA for splicing events.
#
# For each condition, event type and direction (all/pos/neg) this script
# extracts the unique gene set and performs Gene Ontology (GO) enrichment
# analysis using the clusterProfiler package.  It then computes a ranked
# vector of gene statistics based on the absolute mean ΔPSI across events and
# runs FGSEA to identify enriched pathways from the MSigDB C5 collection.
# Results are saved to the output directory `05_go_fgsea/` with separate
# subdirectories for GO and FGSEA.

message("Step 05: Running GO enrichment and FGSEA")

## Check prerequisites ---------------------------------------------------------
if (!exists("EVENT_SPLITS", envir = .GlobalEnv)) {
  stop("EVENT_SPLITS object is missing; please run steps 01–03 first")
}
event_splits <- get("EVENT_SPLITS", envir = .GlobalEnv)

## Create output directories ---------------------------------------------------
step_output <- file.path(PROJECT_ROOT, "output", "05_go_fgsea")
go_dir <- file.path(step_output, "go")
fgsea_dir <- file.path(step_output, "fgsea")
dir.create(go_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fgsea_dir, recursive = TRUE, showWarnings = FALSE)

## Load required libraries -----------------------------------------------------
suppressMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(fgsea)
  library(msigdbr)
  library(AnnotationDbi)
  library(ggplot2)
})

## Load MSigDB C5 gene sets for FGSEA ----------------------------------------
msigdb_list <- NULL
try({
  msigdb_df <- msigdbr(species = "Homo sapiens", category = "C5")
  msigdb_list <- split(msigdb_df$gene_symbol, msigdb_df$gs_name)
}, silent = TRUE)
if (is.null(msigdb_list)) {
  warning("Failed to load MSigDB C5 gene sets; FGSEA will be skipped")
}

## Helper function: convert gene symbols to Entrez IDs ------------------------
symbol_to_entrez <- function(symbols) {
  unique_symbols <- unique(symbols)
  entrez <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = unique_symbols,
                                  keytype = "SYMBOL", column = "ENTREZID",
                                  multiVals = "first")
  entrez <- entrez[!is.na(entrez)]
  return(entrez)
}

## Helper: prepare ranking vector from ΔPSI -----------------------------------
prepare_rank_vector <- function(df) {
  # Identify dPSI column
  dpsi_col <- intersect(c("E.dPsi.", "dPsi", "DeltaPSI"), names(df))
  if (length(dpsi_col) == 0) {
    warning("No dPSI column found; cannot compute ranking vector")
    return(NULL)
  }
  gene_col <- intersect(c("GeneName", "GENEID", "GENE"), names(df))
  if (length(gene_col) == 0) {
    warning("No gene column found; cannot compute ranking vector")
    return(NULL)
  }
  stats <- df %>%
    group_by(.data[[gene_col[1]]]) %>%
    summarise(stat = mean(abs(.data[[dpsi_col[1]]]), na.rm = TRUE), .groups = "drop")
  ranking <- stats$stat
  names(ranking) <- stats[[gene_col[1]]]
  # Remove duplicates
  ranking <- ranking[!is.na(ranking)]
  return(ranking)
}

## Iterate through conditions, event types and directions ----------------------
conditions <- names(event_splits)
event_types <- unique(unlist(lapply(event_splits, names)))
directions <- c("all", "pos", "neg")

go_results <- list()
fgsea_results <- list()

for (cond in conditions) {
  go_results[[cond]] <- list()
  fgsea_results[[cond]] <- list()
  for (etype in event_types) {
    go_results[[cond]][[etype]] <- list()
    fgsea_results[[cond]][[etype]] <- list()
    for (dirn in directions) {
      # Access events for this condition/type/direction
      if (is.null(event_splits[[cond]][[etype]][[dirn]])) next
      df <- event_splits[[cond]][[etype]][[dirn]]
      # Extract gene symbols
      gene_col <- intersect(c("GeneName", "GENEID", "GENE"), names(df))
      if (length(gene_col) == 0) {
        warning(sprintf("No gene column found in %s %s %s; skipping", cond, etype, dirn))
        next
      }
      genes <- unique(df[[gene_col[1]]])
      # GO over‑representation
      entrez <- symbol_to_entrez(genes)
      if (length(entrez) > 0) {
        go_res <- tryCatch({
          enrichGO(gene = entrez,
                   OrgDb = org.Hs.eg.db,
                   ont = cfg$go$ontology,
                   keyType = "ENTREZID",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2)
        }, error = function(e) {
          warning(sprintf("enrichGO failed for %s %s %s: %s", cond, etype, dirn, e$message))
          NULL
        })
        # Simplify terms if requested
        if (!is.null(go_res) && nrow(go_res) > 0) {
          go_res <- clusterProfiler::simplify(go_res, cutoff = cfg$go$simplify_cutoff)
          # Filter by size
          go_res <- go_res[go_res$Count >= cfg$go$count_min & go_res$Count <= cfg$go$count_max, ]
        }
        go_results[[cond]][[etype]][[dirn]] <- go_res
        # Save GO table
        if (!is.null(go_res) && nrow(go_res) > 0) {
          go_file <- file.path(go_dir, sprintf("GO_%s_%s_%s.tsv", cond, etype, dirn))
          readr::write_tsv(as.data.frame(go_res), go_file)
        }
      } else {
        go_results[[cond]][[etype]][[dirn]] <- NULL
      }
      # FGSEA
      if (!is.null(msigdb_list)) {
        ranking <- prepare_rank_vector(df)
        if (!is.null(ranking) && length(ranking) > 0) {
          # Keep only genes present in ranking and gene sets
          fg_res <- tryCatch({
            fgsea(pathways = msigdb_list,
                  stats = ranking,
                  nperm = cfg$fgsea$nperm,
                  minSize = 10,
                  maxSize = 500)
          }, error = function(e) {
            warning(sprintf("FGSEA failed for %s %s %s: %s", cond, etype, dirn, e$message))
            NULL
          })
          if (!is.null(fg_res)) {
            fg_res <- fg_res[order(fg_res$pval), ]
            # Filter by adjusted p-value and NES threshold
            fg_res <- fg_res[fg_res$padj <= cfg$fgsea$padj_max & abs(fg_res$NES) >= cfg$fgsea$abs_nes_min, ]
            fgsea_results[[cond]][[etype]][[dirn]] <- fg_res
            # Save FGSEA table
            fg_file <- file.path(fgsea_dir, sprintf("FGSEA_%s_%s_%s.tsv", cond, etype, dirn))
            readr::write_tsv(as.data.frame(fg_res), fg_file)
            # Plot top N pathways
            top_n <- min(cfg$fgsea$top_n, nrow(fg_res))
            if (top_n > 0) {
              plot_df <- fg_res[1:top_n, ]
              plot_df$gs_name <- factor(plot_df$pathway, levels = rev(plot_df$pathway))
              plt <- ggplot(plot_df, aes(x = gs_name, y = NES)) +
                geom_col(aes(fill = NES > 0)) +
                coord_flip() +
                labs(title = sprintf("FGSEA Top %d (%s %s %s)", top_n, cond, etype, dirn),
                     x = "Pathway", y = "Normalized Enrichment Score") +
                scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02"), guide = FALSE) +
                theme_minimal()
              plot_path <- file.path(fgsea_dir, sprintf("FGSEA_top_%s_%s_%s.png", cond, etype, dirn))
              ggsave(plot_path, plt, width = 7, height = 4 + 0.2 * top_n)
            }
          } else {
            fgsea_results[[cond]][[etype]][[dirn]] <- NULL
          }
        }
      }
    }
  }
}

## Save results as RDS ---------------------------------------------------------
go_rds <- file.path(step_output, "go_results.rds")
fgsea_rds <- file.path(step_output, "fgsea_results.rds")
saveRDS(go_results, go_rds)
saveRDS(fgsea_results, fgsea_rds)
message(sprintf("GO and FGSEA results saved to %s and %s", go_rds, fgsea_rds))

## Return results to global environment ---------------------------------------
assign("GO_RESULTS", go_results, envir = .GlobalEnv)
assign("FGSEA_RESULTS", fgsea_results, envir = .GlobalEnv)
message("Step 05 completed")