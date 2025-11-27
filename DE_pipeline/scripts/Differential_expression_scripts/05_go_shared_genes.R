## scripts/08_go_shared_genes.R
## GO over-representation for shared / non-shared gene groups
## derived from gene overlap sets (per-condition DEG sets).
##
## This script is generic and does not assume a specific set of conditions.
## It:
##  - loads gene_sets_all/pos/neg_by_condition.rds from 07_gene_overlap
##  - builds exact membership groups (e.g. A, A&B, A&B&C, ...)
##  - runs GO for each group using run_enrich_go_for_result()
##  - exports xlsx tables and barplots (top N terms) in 09_go_shared_genes

if (!exists("output_dir"))         stop("'output_dir' not found.")
if (!exists("DE_PADJ"))            stop("'DE_PADJ' not found.")
if (!exists("DE_QVAL"))            stop("'DE_QVAL' not found.")
if (!exists("GO_ONTOLOGY"))        stop("'GO_ONTOLOGY' not found.")
if (!exists("GO_SIMPLIFY_CUTOFF")) stop("'GO_SIMPLIFY_CUTOFF' not found.")
if (!exists("GO_COUNT_MIN"))       stop("'GO_COUNT_MIN' not found.")
if (!exists("GO_COUNT_MAX"))       stop("'GO_COUNT_MAX' not found.")
if (!exists("GO_TOP_N")) {
  GO_TOP_N <- 5L
}

## ----------------------------------------------------------------------
## Paths to input RDS files produced by 06_gene_overlap_plots.R
## ----------------------------------------------------------------------

gene_ov_dir <- file.path(output_dir, "07_gene_overlap")
if (!dir.exists(gene_ov_dir)) {
  stop("Gene overlap directory not found: ", gene_ov_dir)
}

path_all <- file.path(gene_ov_dir, "gene_sets_all_by_condition.rds")
path_pos <- file.path(gene_ov_dir, "gene_sets_pos_by_condition.rds")
path_neg <- file.path(gene_ov_dir, "gene_sets_neg_by_condition.rds")

if (!file.exists(path_all)) stop("File not found: ", path_all)
if (!file.exists(path_pos)) stop("File not found: ", path_pos)
if (!file.exists(path_neg)) stop("File not found: ", path_neg)

gene_sets_all <- readRDS(path_all)
gene_sets_pos <- readRDS(path_pos)
gene_sets_neg <- readRDS(path_neg)

## Consistency check on names
cond_names <- intersect(
  intersect(names(gene_sets_all), names(gene_sets_pos)),
  names(gene_sets_neg)
)
if (!length(cond_names)) {
  stop("No common condition names across gene_sets_all/pos/neg.")
}

gene_sets_all <- gene_sets_all[cond_names]
gene_sets_pos <- gene_sets_pos[cond_names]
gene_sets_neg <- gene_sets_neg[cond_names]

## ----------------------------------------------------------------------
## Output directory for shared-groups GO
## ----------------------------------------------------------------------

out_dir_go_shared <- ensure_dir(output_dir, "09_go_shared_genes")

## ----------------------------------------------------------------------
## Helpers: build exact membership groups and ordering
## ----------------------------------------------------------------------

# Build groups of genes defined by exact membership patterns across conditions.
# Input: named list of character vectors (Ensembl IDs) per condition.
# Output: named list of character vectors, where each name is a pattern
#         such as "CondA&CondB", "CondA", etc.
build_exact_membership_groups <- function(sets_list) {
  if (!length(sets_list)) {
    return(list())
  }

  # Clean sets: ensure unique, non-NA character values
  sets_list <- lapply(sets_list, function(v) {
    v <- as.character(v)
    unique(v[!is.na(v) & nzchar(v)])
  })

  conds <- names(sets_list)
  genes <- sort(unique(unlist(sets_list, use.names = FALSE)))

  if (!length(genes)) {
    return(list())
  }

  # Build membership matrix: rows = genes, columns = conditions
  membership_mat <- matrix(
    FALSE,
    nrow = length(genes),
    ncol = length(conds),
    dimnames = list(genes, conds)
  )

  for (nm in conds) {
    v <- sets_list[[nm]]
    idx <- match(v, genes, nomatch = 0L)
    idx <- idx[idx > 0L]
    membership_mat[idx, nm] <- TRUE
  }

  # Pattern is the set of conditions a gene belongs to
  patterns <- apply(membership_mat, 1L, function(x) {
    if (!any(x)) {
      return("")
    }
    paste(conds[which(x)], collapse = "&")
  })

  keep <- nzchar(patterns)
  genes_kept <- genes[keep]
  patterns_kept <- patterns[keep]

  split(genes_kept, patterns_kept)
}

# Order groups by size (largest first)
order_groups_by_size <- function(groups) {
  if (!length(groups)) return(groups)
  lens <- vapply(groups, length, integer(1L))
  groups[order(lens, decreasing = TRUE)]
}

# Sanitize group name for file system
sanitize_group_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

## ----------------------------------------------------------------------
## Build groups for ALL / POS / NEG
## ----------------------------------------------------------------------

groups_all <- build_exact_membership_groups(gene_sets_all)
groups_pos <- build_exact_membership_groups(gene_sets_pos)
groups_neg <- build_exact_membership_groups(gene_sets_neg)

groups_all <- order_groups_by_size(groups_all)
groups_pos <- order_groups_by_size(groups_pos)
groups_neg <- order_groups_by_size(groups_neg)

# Save group definitions for inspection / reuse
saveRDS(groups_all, file.path(out_dir_go_shared, "shared_groups_ALL_ensembl.rds"))
saveRDS(groups_pos, file.path(out_dir_go_shared, "shared_groups_POS_ensembl.rds"))
saveRDS(groups_neg, file.path(out_dir_go_shared, "shared_groups_NEG_ensembl.rds"))

## ----------------------------------------------------------------------
## Run GO for each group and direction using run_enrich_go_for_result()
## ----------------------------------------------------------------------

GO_SHARED_ALL     <- list()
GO_SHARED_POS     <- list()
GO_SHARED_NEG     <- list()

GO_SHARED_ALL_DF  <- list()
GO_SHARED_POS_DF  <- list()
GO_SHARED_NEG_DF  <- list()

run_go_for_groups <- function(groups_list,
                              direction_label,
                              fill_color,
                              go_obj_list,
                              go_df_list) {
  if (!length(groups_list)) {
    return(list(go_obj_list = go_obj_list, go_df_list = go_df_list))
  }

  for (grp_name in names(groups_list)) {
    genes <- groups_list[[grp_name]]
    if (!length(genes)) {
      next
    }

    # Build minimal df for run_enrich_go_for_result (expects id_col "ensembl_id")
    res_df <- data.frame(
      ensembl_id = as.character(genes),
      stringsAsFactors = FALSE
    )

    ego <- run_enrich_go_for_result(
      res_df            = res_df,
      ontology          = GO_ONTOLOGY,
      padj_threshold    = DE_PADJ,
      q_value_threshold = DE_QVAL,
      simplify_cutoff   = GO_SIMPLIFY_CUTOFF,
      orgdb             = org.Hs.eg.db,
      count_min         = GO_COUNT_MIN,
      count_max         = GO_COUNT_MAX
    )

    go_obj_list[[grp_name]] <- ego

    if (!is.null(ego)) {
      df <- as.data.frame(ego)
      go_df_list[[grp_name]] <- df

      grp_safe <- sanitize_group_name(grp_name)

      # Export table
      outfile_xlsx <- file.path(
        out_dir_go_shared,
        paste0("GO_shared_", grp_safe, "_", direction_label, ".xlsx")
      )
      openxlsx::write.xlsx(df, outfile_xlsx, rowNames = FALSE)

      # Plot top N terms (same style as go_overrep)
      outfile_plot <- file.path(
        out_dir_go_shared,
        paste0("GO_shared_", grp_safe, "_", direction_label,
               "_top", GO_TOP_N, ".png")
      )
      plot_go_top_terms(
        go_df      = df,
        n_top      = GO_TOP_N,
        title      = paste0(
          "GO over-representation (", direction_label,
          ") - shared group: ", grp_name
        ),
        fill_color = fill_color,
        outfile    = outfile_plot
      )
    }
  }

  list(go_obj_list = go_obj_list, go_df_list = go_df_list)
}

## ALL (all DE genes in each shared group)
res_all <- run_go_for_groups(
  groups_list  = groups_all,
  direction_label = "ALL",
  fill_color      = "#4CAF50",  # green, as in go_overrep
  go_obj_list     = GO_SHARED_ALL,
  go_df_list      = GO_SHARED_ALL_DF
)
GO_SHARED_ALL    <- res_all$go_obj_list
GO_SHARED_ALL_DF <- res_all$go_df_list

## POS (upregulated genes only)
res_pos <- run_go_for_groups(
  groups_list  = groups_pos,
  direction_label = "POS",
  fill_color      = "#FFC107",  # amber, as in go_overrep
  go_obj_list     = GO_SHARED_POS,
  go_df_list      = GO_SHARED_POS_DF
)
GO_SHARED_POS    <- res_pos$go_obj_list
GO_SHARED_POS_DF <- res_pos$go_df_list

## NEG (downregulated genes only)
res_neg <- run_go_for_groups(
  groups_list  = groups_neg,
  direction_label = "NEG",
  fill_color      = "#1E88E5",  # blue, as in go_overrep
  go_obj_list     = GO_SHARED_NEG,
  go_df_list      = GO_SHARED_NEG_DF
)
GO_SHARED_NEG    <- res_neg$go_obj_list
GO_SHARED_NEG_DF <- res_neg$go_df_list

## ----------------------------------------------------------------------
## Save GO objects for downstream use
## ----------------------------------------------------------------------

saveRDS(GO_SHARED_ALL,
        file.path(out_dir_go_shared, "GO_shared_ALL_list.rds"))
saveRDS(GO_SHARED_POS,
        file.path(out_dir_go_shared, "GO_shared_POS_list.rds"))
saveRDS(GO_SHARED_NEG,
        file.path(out_dir_go_shared, "GO_shared_NEG_list.rds"))

saveRDS(GO_SHARED_ALL_DF,
        file.path(out_dir_go_shared, "GO_shared_ALL_df_list.rds"))
saveRDS(GO_SHARED_POS_DF,
        file.path(out_dir_go_shared, "GO_shared_POS_df_list.rds"))
saveRDS(GO_SHARED_NEG_DF,
        file.path(out_dir_go_shared, "GO_shared_NEG_df_list.rds"))
