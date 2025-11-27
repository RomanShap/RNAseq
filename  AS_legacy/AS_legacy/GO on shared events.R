# --- Minimal GO enrichment from a gene list (human) ---------------------------
# Install once if needed:
# install.packages(c("writexl","readr"))
# BiocManager::install(c("clusterProfiler","org.Hs.eg.db"))
# 
# suppressPackageStartupMessages({
#   library(clusterProfiler)
#   library(org.Hs.eg.db)
#   library(writexl)
#   library(readr)
#   library(dplyr)
# })
# 
# # Read genes from a file (txt/csv) or accept a character vector
# read_gene_list <- function(x) {
#   if (is.character(x) && length(x) == 1 && file.exists(x)) {
#     # path
#     ext <- tolower(tools::file_ext(x))
#     if (ext %in% c("txt","tsv")) {
#       unique(na.omit(trimws(read_lines(x))))
#     } else {
#       df <- readr::read_delim(x, delim = ifelse(ext=="csv", ",", "\t"), show_col_types = FALSE)
#       col <- intersect(tolower(names(df)), c("symbol","gene","genes","hgnc_symbol"))[1]
#       if (is.na(col)) stop("File must have a 'symbol' column or be 1 gene per line.")
#       unique(na.omit(trimws(df[[col]])))
#     }
#   } else if (is.character(x)) {
#     unique(na.omit(trimws(x)))
#   } else stop("Provide a character vector of symbols or a path to a txt/csv.")
# }
# 
# # Map SYMBOL -> ENTREZ
# symbols_to_entrez <- function(sym_vec) {
#   if (!length(sym_vec)) return(character(0))
#   m <- suppressMessages(AnnotationDbi::select(
#     org.Hs.eg.db, keys = sym_vec, keytype = "SYMBOL", columns = "ENTREZID"
#   ))
#   unique(na.omit(m$ENTREZID))
# }
# 
# 
# go_enrich <- function(
#     genes,
#     out_xlsx = "GO_results.xlsx",
#     p_cut = 0.05, q_cut = 0.05,
#     min_size = 10, max_size = 5000,
#     universe = NULL,
#     # --- NEW: simplify options ---
#     simplify = TRUE,                  # apply clusterProfiler::simplify()
#     simplify_cutoff = 0.6,            # 0–1; higher = more aggressive merging
#     simplify_measure = c("Wang","Rel","Resnik","Jiang","Lin")[1],
#     simplify_by = c("p.adjust","pvalue","qvalue")[1],
#     simplify_select_fun = min,        # which term to keep within a cluster
#     keep_raw_sheet = FALSE            # also write a *_RAW sheet per ontology
# ) {
#   sym <- read_gene_list(genes)
#   entrez <- symbols_to_entrez(sym)
#   if (length(entrez) < min_size) {
#     message("Not enough mapped genes: n=", length(entrez), " < minGSSize=", min_size)
#     return(invisible(NULL))
#   }
# 
#   # Universe: if provided as symbols, map; if numeric, assume ENTREZ
#   if (!is.null(universe)) {
#     if (is.character(universe)) universe <- symbols_to_entrez(universe)
#     if (!length(universe)) universe <- NULL
#   }
# 
#   onts <- c("BP","MF","CC")
#   sheets <- list()
#   results <- list()
# 
#   for (ont in onts) {
#     eg <- tryCatch(
#       enrichGO(
#         gene          = entrez,
#         OrgDb         = org.Hs.eg.db,
#         keyType       = "ENTREZID",
#         ont           = ont,
#         universe      = universe,
#         pAdjustMethod = "BH",
#         pvalueCutoff  = p_cut,
#         qvalueCutoff  = q_cut,
#         minGSSize     = min_size,
#         maxGSSize     = max_size,
#         readable      = TRUE
#       ),
#       error = function(e) NULL
#     )
# 
#     if (is.null(eg) || nrow(as.data.frame(eg)) == 0) {
#       sheets[[ont]] <- tibble::tibble(
#         note = sprintf("No significant GO terms (p<=%.3f, FDR/q<=%.3f).", p_cut, q_cut),
#         Input_n = length(entrez),
#         Universe_n = ifelse(is.null(universe), NA_integer_, length(universe))
#       )
#       results[[ont]] <- sheets[[ont]]
#       next
#     }
# 
#     # (Optional) keep the raw table
#     raw_df <- as.data.frame(eg)
# 
#     # --- NEW: simplify to remove redundant terms ----------------------------
#     if (isTRUE(simplify)) {
#       eg_simpl <- tryCatch(
#         clusterProfiler::simplify(
#           eg,
#           cutoff     = simplify_cutoff,
#           by         = simplify_by,
#           select_fun = simplify_select_fun,
#           measure    = simplify_measure
#         ),
#         error = function(e) {
#           message("simplify() failed for ", ont, ": ", conditionMessage(e))
#           NULL
#         }
#       )
#       if (!is.null(eg_simpl) && nrow(as.data.frame(eg_simpl)) > 0) {
#         eg <- eg_simpl
#       }
#     }
# 
#     df <- as.data.frame(eg) %>%
#       dplyr::mutate(
#         Ontology  = ont,
#         Input_n   = length(entrez),
#         Universe_n= ifelse(is.null(universe), NA_integer_, length(universe))
#       ) %>%
#       dplyr::relocate(Ontology, .before = 1)
# 
#     # write main (simplified) sheet
#     sheets[[ont]]  <- df
#     results[[ont]] <- df
# 
#     # Optionally also write the raw (unsimplified) sheet for comparison
#     if (isTRUE(keep_raw_sheet)) {
#       sheets[[paste0(ont, "_RAW")]] <- raw_df %>%
#         dplyr::mutate(
#           Ontology  = ont,
#           Input_n   = length(entrez),
#           Universe_n= ifelse(is.null(universe), NA_integer_, length(universe))
#         ) %>%
#         dplyr::relocate(Ontology, .before = 1)
#     }
#   }
# 
#   writexl::write_xlsx(sheets, path = out_xlsx)
#   message("Saved: ", out_xlsx)
#   invisible(results)
# }
# 
# 
# # --------------------- EXAMPLES ---------------------
# # 1) From a character vector:
# my_genes <- final_upset_export$FUS
# go_res <- go_enrich(my_genes, out_xlsx = "GO_results.xlsx")

# 2) From a file (txt: one gene per line; or CSV with a 'symbol' column):
# go_res <- go_enrich("my_gene_list.txt", out_xlsx = "GO_results.xlsx")

# 3) With a custom background (universe) defined as all detected genes (symbols):
# detected <- read_lines("all_detected_symbols.txt")
# go_res <- go_enrich("my_gene_list.txt", out_xlsx = "GO_results.xlsx", universe = detected)


# ---- Assume go_enrich() from the previous message is already defined ----
# (If not, paste that function first.)

suppressPackageStartupMessages({
  library(readxl)
  library(fs)
  library(stringr)
})

# # Batch GO from "shared" Excel(s): one column = one gene list
# go_batch_from_shared_excels <- function(
#     xlsx_paths,
#     input_dir  = NULL, 
#     out_dir = "GO_by_column",
#     universe = NULL,
#     p_cut = 0.05, q_cut = 0.05,
#     min_size = 10, max_size = 5000
# ) {
#   dir_create(out_dir)
#   
#   xlsx_paths <- as.character(xlsx_paths)
#   for (path in xlsx_paths) {
#     if (!file.exists(path)) {
#       message("Skipping missing file: ", path)
#       next
#     }
#     message("\n=== Reading: ", path, " ===")
#     df <- readxl::read_xlsx(path)
#     
#     # Keep only non-empty columns
#     keep_cols <- names(df)[colSums(!is.na(df)) > 0]
#     df <- df[, keep_cols, drop = FALSE]
#     if (!ncol(df)) {
#       message("No non-empty columns in: ", path)
#       next
#     }
#     
#     base <- tools::file_path_sans_ext(basename(path))
#     for (col in names(df)) {
#       genes <- unique(na.omit(trimws(df[[col]])))
#       if (!length(genes)) {
#         message("  - ", col, ": empty, skipping.")
#         next
#       }
#       
#       # Safe filename: <file>__<col>_GO.xlsx
#       safe_col <- str_replace_all(col, "[^A-Za-z0-9._-]+", "_")
#       out_file <- file.path(out_dir, sprintf("%s__%s_GO.xlsx", base, safe_col))
#       
#       message("  - ", col, ": n=", length(genes), " → ", out_file)
#       go_enrich(
#         genes     = genes,
#         out_xlsx  = out_file,
#         universe  = universe,
#         p_cut     = p_cut,
#         q_cut     = q_cut,
#         min_size  = min_size,
#         max_size  = max_size
#       )
#     }
#   }
#   invisible(TRUE)
# }
go_batch_from_shared_excels <- function(
    xlsx_paths,
    input_dir  = NULL, 
    out_dir = "GO_by_column",
    universe = NULL,
    p_cut = 0.05, q_cut = 0.05,
    min_size = 10, max_size = 5000
) {
  dir_create(out_dir)
  
  # --- minimal addition: pull files from input_dir, if given ---
  xlsx_paths <- as.character(xlsx_paths)
  if (!is.null(input_dir)) {
    more <- list.files(input_dir, pattern = "\\.(xlsx|xlsm)$",
                       full.names = TRUE, ignore.case = TRUE)
    xlsx_paths <- unique(c(xlsx_paths, more))
  }
  xlsx_paths <- xlsx_paths[file.exists(xlsx_paths)]
  if (!length(xlsx_paths)) stop("No input Excel files found.")
  
  for (path in xlsx_paths) {
    if (!file.exists(path)) {
      message("Skipping missing file: ", path)
      next
    }
    message("\n=== Reading: ", path, " ===")
    df <- readxl::read_xlsx(path)
    
    # Keep only non-empty columns
    keep_cols <- names(df)[colSums(!is.na(df)) > 0]
    df <- df[, keep_cols, drop = FALSE]
    if (!ncol(df)) {
      message("No non-empty columns in: ", path)
      next
    }
    
    base <- tools::file_path_sans_ext(basename(path))
    for (col in names(df)) {
      genes <- unique(na.omit(trimws(df[[col]])))
      if (!length(genes)) {
        message("  - ", col, ": empty, skipping.")
        next
      }
      
      # Safe filename: <file>__<col>_GO.xlsx
      safe_col <- stringr::str_replace_all(col, "[^A-Za-z0-9._-]+", "_")
      out_file <- file.path(out_dir, sprintf("%s__%s_GO.xlsx", base, safe_col))
      
      message("  - ", col, ": n=", length(genes), " → ", out_file)
      go_enrich(
        genes     = genes,
        out_xlsx  = out_file,
        universe  = universe,
        p_cut     = p_cut,
        q_cut     = q_cut,
        min_size  = min_size,
        max_size  = max_size
      )
    }
  }
  invisible(TRUE)
}

# ----------------- EXAMPLES -----------------
# If your “shared” files are like:
#   /path/UpSet_gene_list.xlsx
#   /path/pos_UpSet_gene_list.xlsx
#   /path/neg_UpSet_gene_list.xlsx

# 1) Run all three in one shot
go_batch_from_shared_excels(
  xlsx_paths = c("output/genes_by_combo_EX.xlsx",
                 "output/genes_by_combo_INT.xlsx",
                 "output/genes_by_combo_ALTA.xlsx",
                 "output/genes_by_combo_ALTD.xlsx"),
  out_dir = "output/GO2"  # will be created if missing
)



