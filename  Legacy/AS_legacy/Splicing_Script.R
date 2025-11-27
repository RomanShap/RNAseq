# Clear the Environment ----
rm(list = ls())            # Remove all objects
graphics.off()             # Close all open plots
cat("\014")                # Clear console (works in RStudio)
# import necessary libraries ----
library(clusterProfiler)
library(org.Hs.eg.db)
library(writexl)
library(openxlsx)
library(dplyr)
library(ggplot2)
# import datasets directly from files ----
# Starting from Vast-tools tab files
base_dir <- "D:/FETseq/vast-tools_results"
sample_id <- dir(base_dir)
sample_id <- sample_id[!grepl("INCLUSION_LEVELS", sample_id)]
# create a list with kallisto tsv directories
splicing_tab_paths <- file.path(base_dir, sample_id)
names(splicing_tab_paths) <- sample_id
# Read each tab-separated file into a data frame
splicing_data_list <- lapply(splicing_tab_paths, function(path) {
  read.delim(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
})
EWSR1_splicing_table <- splicing_data_list[["diff_EWSR1KO_vs_WT.tab"]]
FUS_splicing_table <- splicing_data_list[["diff_FUSKO_vs_WT.tab"]]
TAF15_splicing_table <- splicing_data_list[["diff_TAF15KO_vs_WT.tab"]]

# Helper function to filter splicing table
filter_splicing <- function(tbl, mv_thresh = 0.05) {
  list(
    all = tbl[tbl$MV.dPsi._at_0.95 > mv_thresh, ],
    pos = tbl[tbl$MV.dPsi._at_0.95 > mv_thresh & tbl$E.dPsi. > 0, ],
    neg = tbl[tbl$MV.dPsi._at_0.95 > mv_thresh & tbl$E.dPsi. < 0, ]
  )
}

# Apply to all datasets
EWSR1_filtered <- filter_splicing(EWSR1_splicing_table)
FUS_filtered <- filter_splicing(FUS_splicing_table)
TAF15_filtered <- filter_splicing(TAF15_splicing_table)

# Access results
EWSR1_splicing_filtered      <- EWSR1_filtered$all
EWSR1_splicing_filtered_pos  <- EWSR1_filtered$pos
EWSR1_splicing_filtered_neg  <- EWSR1_filtered$neg

FUS_splicing_filtered        <- FUS_filtered$all
FUS_splicing_filtered_pos    <- FUS_filtered$pos
FUS_splicing_filtered_neg    <- FUS_filtered$neg

TAF15_splicing_filtered      <- TAF15_filtered$all
TAF15_splicing_filtered_pos  <- TAF15_filtered$pos
TAF15_splicing_filtered_neg  <- TAF15_filtered$neg

# Define samples, conditions, and event types
samples <- c("EWSR1", "FUS", "TAF15")
conditions <- c("splicing_filtered", "splicing_filtered_pos", "splicing_filtered_neg")
event_types <- c("EX", "INT", "ALTA", "ALTD")

# Helper function to get regex from event type
event_regex <- function(event_type) paste0("^Hsa", event_type)

# Create filtered lists and gene sets
filtered_lists <- list()
gene_sets <- list()

for (sample in samples) {
  for (condition in conditions) {
    df_name <- paste0(sample, "_", condition)
    df <- get(df_name)
    
    for (event in event_types) {
      pattern <- event_regex(event)
      filtered <- df[grepl(pattern, df$EVENT), ]
      
      # Assign to filtered list
      filtered_name <- paste0(sample, "_filtered_", event, ifelse(condition == "splicing_filtered", "", sub("splicing_filtered", "", condition)))
      assign(filtered_name, filtered)
      filtered_lists[[filtered_name]] <- filtered
      
      # Assign to gene list
      gene_name <- paste0("genes_", sample, "_", event, ifelse(condition == "splicing_filtered", "", sub("splicing_filtered", "", condition)))
      assign(gene_name, unique(filtered$GENE))
      gene_sets[[gene_name]] <- unique(filtered$GENE)
    }
  }
}
# Access any filtered table via filtered_lists[["EWSR1_filtered_EX"]], etc.
# Access gene sets via gene_sets[["genes_TAF15_ALTA_pos"]], etc.


# Run GO on the genes found by vast-tools and filtered by MV and dPSI ----

run_go <- function(gene_symbols, name = "GO_result") {
  # Convert gene symbols to Entrez IDs
  gene_ids <- bitr(gene_symbols,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
  
  # If empty after conversion, return NULL
  if (nrow(gene_ids) == 0) return(NULL)
  
  # Run GO enrichment
  enrichGO(gene         = gene_ids$ENTREZID,
           OrgDb        = org.Hs.eg.db,
           keyType      = "ENTREZID",
           ont          = "ALL",  # Biological Process
           pAdjustMethod = "BH",
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.2,
           readable     = TRUE)
}
# Gene lists with the whole range of dPSI ----
samples <- c("EWSR1", "FUS", "TAF15")
event_types <- c("EX", "INT", "ALTA", "ALTD")
conditions <- c("", "_pos", "_neg")

# Create gene_lists grouped by condition
gene_lists_all <- list()
gene_lists_pos <- list()
gene_lists_neg <- list()

for (sample in samples) {
  for (event in event_types) {
    for (cond in conditions) {
      var_name <- paste0("genes_", sample, "_", event, cond)
      list_name <- paste0(sample, "_", event, cond)
      value <- get(var_name)
      
      if (cond == "") {
        gene_lists_all[[list_name]] <- value
      } else if (cond == "_pos") {
        gene_lists_pos[[list_name]] <- value
      } else if (cond == "_neg") {
        gene_lists_neg[[list_name]] <- value
      }
    }
  }
}
# 
# # Whole range of dPSI
# go_results_all <- lapply(names(gene_lists_all), function(name) {
#   cat("Running GO for:", name, "\n")
#   run_go(gene_lists_all[[name]], name = name)
# })
# names(go_results_all) <- names(gene_lists_all)
# 
# for (name in names(go_results_all)) {
#   res <- go_results_all[[name]]
#   if (!is.null(res) && nrow(res) > 0) {
#     write_xlsx(as.data.frame(res), paste0("GO_", name, ".xlsx"))
#   }
# }
# # Positive
# go_results_pos <- lapply(names(gene_lists_pos), function(name) {
#   cat("Running GO for:", name, "\n")
#   run_go(gene_lists_pos[[name]], name = name)
# })
# names(go_results_pos) <- names(gene_lists_pos)
# 
# for (name in names(go_results_pos)) {
#   res <- go_results_pos[[name]]
#   if (!is.null(res) && nrow(res) > 0) {
#     write_xlsx(as.data.frame(res), paste0("GO_pos_", name, ".xlsx"))
#   }
# }
# # Negative
# go_results_neg <- lapply(names(gene_lists_neg), function(name) {
#   cat("Running GO for:", name, "\n")
#   run_go(gene_lists_neg[[name]], name = name)
# })
# names(go_results_neg) <- names(gene_lists_neg)
# 
# for (name in names(go_results_neg)) {
#   res <- go_results_neg[[name]]
#   if (!is.null(res) && nrow(res) > 0) {
#     write_xlsx(as.data.frame(res), paste0("GO_neg_", name, ".xlsx"))
#   }
# }
# for (name in names(go_results)) {
#   res <- go_results[[name]]
#   if (!is.null(res) && nrow(res) > 0) {
#     write.csv(as.data.frame(res), paste0("GO_", name, ".csv"), row.names = FALSE)
#   }
# }
# write_xlsx(go_results[["EWSR1_EX"]]@result, "GO_EWSR1_EX.xlsx")

run_go <- function(gene_symbols, name = "GO_result") {
  # Convert gene symbols to Entrez IDs
  gene_ids <- bitr(gene_symbols,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
  
  if (nrow(gene_ids) == 0) return(NULL)
  
  # Run GO enrichment
  GO_result <- enrichGO(
    gene         = gene_ids$ENTREZID,
    OrgDb        = org.Hs.eg.db,
    keyType      = "ENTREZID",
    ont          = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable     = TRUE
  )
  
  # If no enrichment, return NULL
  if (is.null(GO_result) || nrow(as.data.frame(GO_result)) == 0) return(NULL)
  
  # Try simplifying, catch errors
  GO_result_simple <- tryCatch({
    simplify(GO_result, cutoff = 0.7, by = "p.adjust", select_fun = min)
  }, error = function(e) {
    warning("Simplify failed: ", e$message)
    return(GO_result)
  })
  
  GO_result_simple_df <- as.data.frame(GO_result_simple)
  
  # If no rows after simplify, return just the enrichment result
  if (nrow(GO_result_simple_df) == 0) {
    return(list(GO = GO_result_simple, plot = NULL))
  }
  
  # Proceed to plotting
  plot_GO_result <- GO_result_simple_df %>%
    arrange(p.adjust) %>%
    slice_head(n = 20) %>%
    mutate(
      negLogFDR = -log10(p.adjust),
      Description = factor(Description, levels = rev(Description))
    )
  
  go_plot <- ggplot(plot_GO_result, aes(
    x = negLogFDR,
    y = Description,
    size = Count,
    color = RichFactor,
    shape = ONTOLOGY
  )) +
    geom_point() +
    scale_shape_manual(values = c(CC = 16, MF = 15, BP = 17)) +
    labs(
      title = paste("Top Enriched GO Terms for", name),
      x = expression(-log[10]("FDR")),
      y = "GO Term",
      size = "Mapped Genes",
      color = "RichFactor",
      shape = "Ontology"
    ) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      legend.text = element_text(size = 14, face = "bold"),
      legend.title = element_text(size = 14, face = "bold")
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 6))
  
  return(list(GO = GO_result_simple, plot = go_plot))
}

#Whole range of dPSI
go_results_all <- setNames(
  lapply(names(gene_lists_all), function(name) {
    cat("Running GO for:", name, "\n")
    run_go(gene_lists_all[[name]], name = name)
  }),
  names(gene_lists_all)
)

# Clean up failed runs
go_results_all <- go_results_all[!sapply(go_results_all, is.null)]

#Positive
go_results_pos <- setNames(
  lapply(names(gene_lists_pos), function(name) {
    cat("Running GO for:", name, "\n")
    run_go(gene_lists_pos[[name]], name = name)
  }),
  names(gene_lists_pos)
)

# Clean up failed runs
go_results_pos <- go_results_pos[!sapply(go_results_pos, is.null)]

#Negative
go_results_neg <- setNames(
  lapply(names(gene_lists_neg), function(name) {
    cat("Running GO for:", name, "\n")
    run_go(gene_lists_neg[[name]], name = name)
  }),
  names(gene_lists_neg)
)

# Clean up failed runs
go_results_neg <- go_results_neg[!sapply(go_results_neg, is.null)]

#saving the last part(no redundancy)

# Create a new workbook
wb <- createWorkbook()

# Loop through the results and add each table as a new sheet
for (i in seq_along(go_results_all)) {
  go_table_all <- go_results_all[[i]]$GO
  sheet_name <- names(go_results_all)[i]
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, go_table_all)
}

# Save workbook
saveWorkbook(wb, "GO_enrichment_results.xlsx", overwrite = TRUE)

# Save each plot 
for (i in seq_along(go_results_all)) {
  png_filename <- paste0(names(go_results_all)[i], "_GO.png") 
  png(png_filename, width = 1200, height = 800, res = 150)
  print(go_results_all[[i]]$plot)
  dev.off()
}

# Create a new workbook
wb <- createWorkbook()

# Loop through the results and add each table as a new sheet
for (i in seq_along(go_results_pos)) {
  go_table_pos <- go_results_pos[[i]]$GO
  sheet_name <- names(go_results_pos)[i]
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, go_table_pos)
}

# Save workbook
saveWorkbook(wb, "GO_pos_enrichment_results.xlsx", overwrite = TRUE)

# Save each plot
for (i in seq_along(go_results_pos)) {
  png_filename <- paste0(names(go_results_pos)[i], "_GO_pos.png")
  png(png_filename, width = 1200, height = 800, res = 150)
  print(go_results_pos[[i]]$plot)
  dev.off()
}


# Create a new workbook
wb <- createWorkbook()

# Loop through the results and add each table as a new sheet
for (i in seq_along(go_results_neg)) {
  go_table_neg <- go_results_neg[[i]]$GO
  sheet_name <- names(go_results_neg)[i]
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, go_table_neg)
}

# Save workbook
saveWorkbook(wb, "GO_neg_enrichment_results.xlsx", overwrite = TRUE)

# Save each plot
for (i in seq_along(go_results_neg)) {
  png_filename <- paste0(names(go_results_neg)[i], "_GO_neg.png")
  png(png_filename, width = 1200, height = 800, res = 150)
  print(go_results_neg[[i]]$plot)
  dev.off()
}








plot_go_heatmap <- function(go_result_list, direction = c("positive", "negative")) {
  direction <- match.arg(direction)
  
  # Combine zScores and condition labels
  z_long <- bind_rows(
    lapply(names(go_result_list), function(condition) {
      df <- as.data.frame(go_result_list[[condition]]$GO)
      if (!"Description" %in% names(df) || !"zScore" %in% names(df)) {
        return(NULL)
      }
      df %>%
        select(Description, zScore) %>%
        mutate(Condition = condition)
    }),
    .id = "source"
  ) %>% filter(!is.na(zScore)) %>% distinct()
  
  if (nrow(z_long) == 0) {
    warning("No zScore data to plot.")
    return(NULL)
  }
  
  # Create binary presence matrix
  all_terms <- unique(z_long$Description)
  conditions <- unique(z_long$Condition)
  
  binary_matrix <- expand.grid(Description = all_terms, Condition = conditions, stringsAsFactors = FALSE) %>%
    left_join(z_long, by = c("Description", "Condition")) %>%
    mutate(Present = !is.na(zScore))
  
  heatmap_df <- binary_matrix %>%
    mutate(zScore = ifelse(Present, zScore, NA)) %>%
    mutate(Description = fct_rev(factor(Description, levels = unique(Description))))
  
  # Plot
  p <- ggplot(heatmap_df, aes(x = Condition, y = Description, fill = zScore)) +
    geom_tile(color = "white") +
    geom_text(aes(label = ifelse(!is.na(zScore), sprintf("%.2f", zScore), "")), size = 3.5) +
    scale_fill_gradient2(
      low = "#D55E00", mid = "white", high = "#009E73", midpoint = 0, na.value = "grey90",
      name = "zScore"
    ) +
    labs(
      title = paste("GO Enrichment zScores by Condition -", direction),
      x = "Condition", y = "GO Term"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
  
  return(p)
}

# Plot and save heatmaps
go_heatmap_pos <- plot_go_heatmap(go_results_pos, direction = "positive")
go_heatmap_neg <- plot_go_heatmap(go_results_neg, direction = "negative")

# Save to file
ggsave("GO_heatmap_pos.png", go_heatmap_pos, width = 17, height = 30, dpi = 300)
ggsave("GO_heatmap_neg.png", go_heatmap_neg, width = 17, height = 30, dpi = 300)
library(ComplexUpset)
library(dplyr)
library(ggplot2)

# Function to restructure flat GO result list into nested structure by splicing type and KO factor
nest_go_results <- function(flat_list, types = c("EX", "INT", "ALTA", "ALTD"), factors = c("EWSR1", "FUS", "TAF15"), suffix = "_pos") {
  nested <- list()
  for (type in types) {
    nested[[type]] <- list()
    for (f in factors) {
      entry_name <- paste0(f, "_", type, suffix)
      if (entry_name %in% names(flat_list)) {
        nested[[type]][[f]] <- flat_list[[entry_name]]
      }
    }
  }
  return(nested)
}

# Function to generate a ComplexUpset plot from nested GO enrichment results
plot_upset_go <- function(go_result_per_type, type_label, direction_label) {
  # Extract GO terms for each KO
  EWSR1_terms <- go_result_per_type[["EWSR1"]]$GO$Description
  FUS_terms   <- go_result_per_type[["FUS"]]$GO$Description
  TAF15_terms <- go_result_per_type[["TAF15"]]$GO$Description
  
  all_terms <- unique(c(EWSR1_terms, FUS_terms, TAF15_terms))
  if (length(all_terms) == 0) {
    message("No GO terms to plot for ", type_label, " (", direction_label, ")")
    return(NULL)
  }
  
  # Construct binary matrix for intersections
  binary_matrix <- data.frame(
    Description = all_terms,
    EWSR1KO = all_terms %in% EWSR1_terms,
    FUSKO   = all_terms %in% FUS_terms,
    TAF15KO = all_terms %in% TAF15_terms
  )
  
  rownames(binary_matrix) <- binary_matrix$Description
  binary_matrix <- binary_matrix[, -1]  # remove Description column
  binary_matrix <- binary_matrix %>% mutate(across(everything(), as.integer))
  
  # Generate the plot
  p <- ComplexUpset::upset(
    binary_matrix,
    intersect = c("EWSR1KO", "FUSKO", "TAF15KO"),
    name = "Group",
    width_ratio = 0.15,
    wrap = TRUE,
    sort_intersections_by = 'degree',
    queries = list(
      upset_query(set = 'EWSR1KO', fill = '#EB5757'),
      upset_query(set = 'FUSKO', fill = '#27AE60'),
      upset_query(set = 'TAF15KO', fill = '#4A90E2')
    ),
    themes = upset_default_themes(text = element_text(size = 20, face = "bold")),
    base_annotations = list(
      'Intersection size' = intersection_size(
        text = element_text(size = 6, vjust = -1),
        fill = if (direction_label == "positive") "#009E73" else "#D55E00"
      )
    )
  ) + ggtitle(paste("GO UpSet -", type_label, "-", direction_label)) +
    theme(plot.title = element_text(size = 24, face = "bold"))
  
  return(p)
}

# Build nested lists for both directions
go_results_pos_nested <- nest_go_results(go_results_pos, suffix = "_pos")
go_results_neg_nested <- nest_go_results(go_results_neg, suffix = "_neg")

# Define target splicing types and directions
splicing_types <- c("EX", "INT", "ALTA", "ALTD")
directions <- c("positive", "negative")

# Save all plots to PDF
pdf("GO_ComplexUpset_All.pdf", width = 14, height = 10)
for (type in splicing_types) {
  for (dir in directions) {
    result_list <- if (dir == "positive") go_results_pos_nested[[type]] else go_results_neg_nested[[type]]
    plot <- plot_upset_go(result_list, type, dir)
    if (!is.null(plot)) print(plot)
  }
}
dev.off()

# Save plots individually to PNG
for (type in splicing_types) {
  for (dir in directions) {
    result_list <- if (dir == "positive") go_results_pos_nested[[type]] else go_results_neg_nested[[type]]
    plot <- plot_upset_go(result_list, type, dir)
    if (!is.null(plot)) {
      ggsave(
        filename = paste0("GO_UpSet_", type, "_", dir, ".png"),
        plot = plot, width = 12, height = 9, dpi = 300
      )
    }
  }
}

