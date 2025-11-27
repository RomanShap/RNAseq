# ---- Load required libraries ----
library(dplyr)
library(UpSetR)
library(ComplexUpset)
library(openxlsx)
library(scales)
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

# #only positive events
# EWSR1_filtered_EX <- EWSR1_filtered_EX %>%
# filter(`E.dPsi.` > 0)
# EWSR1_filtered_INT <- EWSR1_filtered_INT %>%
#   filter(`E.dPsi.` > 0)
# EWSR1_filtered_ALTA <- EWSR1_filtered_ALTA %>%
#   filter(`E.dPsi.` > 0)
# EWSR1_filtered_ALTD <- EWSR1_filtered_ALTD %>%
#   filter(`E.dPsi.` > 0)
# 
# FUS_filtered_EX <- FUS_filtered_EX %>%
#   filter(`E.dPsi.` > 0)
# FUS_filtered_INT <- FUS_filtered_INT %>%
#   filter(`E.dPsi.` > 0)
# FUS_filtered_ALTA <- FUS_filtered_ALTA %>%
#   filter(`E.dPsi.` > 0)
# FUS_filtered_ALTD <- FUS_filtered_ALTD %>%
#   filter(`E.dPsi.` > 0)
# 
# TAF15_filtered_EX <- TAF15_filtered_EX %>%
#   filter(`E.dPsi.` > 0)
# TAF15_filtered_INT <- TAF15_filtered_INT %>%
#   filter(`E.dPsi.` > 0)
# TAF15_filtered_ALTA <- TAF15_filtered_ALTA %>%
#   filter(`E.dPsi.` > 0)
# TAF15_filtered_ALTD <- TAF15_filtered_ALTD %>%
#   filter(`E.dPsi.` > 0)

#only negative events
EWSR1_filtered_EX <- EWSR1_filtered_EX %>%
  filter(`E.dPsi.` < 0)
EWSR1_filtered_INT <- EWSR1_filtered_INT %>%
  filter(`E.dPsi.` < 0)
EWSR1_filtered_ALTA <- EWSR1_filtered_ALTA %>%
  filter(`E.dPsi.` < 0)
EWSR1_filtered_ALTD <- EWSR1_filtered_ALTD %>%
  filter(`E.dPsi.` < 0)

FUS_filtered_EX <- FUS_filtered_EX %>%
  filter(`E.dPsi.` < 0)
FUS_filtered_INT <- FUS_filtered_INT %>%
  filter(`E.dPsi.` < 0)
FUS_filtered_ALTA <- FUS_filtered_ALTA %>%
  filter(`E.dPsi.` < 0)
FUS_filtered_ALTD <- FUS_filtered_ALTD %>%
  filter(`E.dPsi.` < 0)

TAF15_filtered_EX <- TAF15_filtered_EX %>%
  filter(`E.dPsi.` < 0)
TAF15_filtered_INT <- TAF15_filtered_INT %>%
  filter(`E.dPsi.` < 0)
TAF15_filtered_ALTA <- TAF15_filtered_ALTA %>%
  filter(`E.dPsi.` < 0)
TAF15_filtered_ALTD <- TAF15_filtered_ALTD %>%
  filter(`E.dPsi.` < 0)

# ---- Define your data ----

event_types <- list(
  EX = list(
    EWSR1 = EWSR1_filtered_EX$EVENT,
    FUS = FUS_filtered_EX$EVENT,
    TAF15 = TAF15_filtered_EX$EVENT
  ),
  INT = list(
    EWSR1 = EWSR1_filtered_INT$EVENT,
    FUS = FUS_filtered_INT$EVENT,
    TAF15 = TAF15_filtered_INT$EVENT
  ),
  ALTA = list(
    EWSR1 = EWSR1_filtered_ALTA$EVENT,
    FUS = FUS_filtered_ALTA$EVENT,
    TAF15 = TAF15_filtered_ALTA$EVENT
  ),
  ALTD = list(
    EWSR1 = EWSR1_filtered_ALTD$EVENT,
    FUS = FUS_filtered_ALTD$EVENT,
    TAF15 = TAF15_filtered_ALTD$EVENT
  )
)

# ---- Function to create an UpSet plot ----

plot_upset <- function(events, type_label) {
  # # #TEST
  # events <- event_types[["INT"]]
  # type_label <- "INT"

  all_events <- unique(unlist(events))
  
  if (length(all_events) == 0) {
    message("Skipping ", type_label, " - no events found.")
    return()
  }
  
  df <- data.frame(
    EVENT = all_events,
    EWSR1 = all_events %in% events[["EWSR1"]],
    FUS   = all_events %in% events[["FUS"]],
    TAF15 = all_events %in% events[["TAF15"]]
  )
  
  binary_df <- df[, c("EWSR1", "FUS", "TAF15")]
  rownames(binary_df) <- df[,1]
  
  pattern_counts <- binary_df %>%
    group_by(EWSR1, FUS, TAF15) %>%
    summarise(Freq = n(), .groups = "drop") %>%
    filter(Freq > 0)
  
  combo_names <- apply(pattern_counts[, 1:3], 1, function(x) {
    selected <- c("EWSR1", "FUS", "TAF15")[as.logical(x)]
    if (length(selected) == 0) return("")
    paste(selected, collapse = "&")
  })
  
  freq_vector <- setNames(pattern_counts$Freq, combo_names)
  freq_vector <- freq_vector[names(freq_vector) != ""]
  
  if (length(freq_vector) == 0) {
    message("No valid intersections for UpSet plot: ", type_label)
    return()
  }
  
  # --- Plot ---
  UpSetR::upset(fromExpression(freq_vector),
        #BLUE
        # main.bar.color = "#0072B2",
        # sets.bar.color = "#56B4E9",
        #GREEN
        # main.bar.color = "#009E73",
        # sets.bar.color = "#69D1AC",
        #RED
        main.bar.color = "#D55E00",
        sets.bar.color = "#F4A582",
        mainbar.y.label = paste("Intersection size -", type_label),
        sets.x.label = "Event count")
  
  ComplexUpset::upset(
    fromExpression(freq_vector),
    intersect = c("EWSR1", "FUS", "TAF15"),
    name = "Group",
    width_ratio = 0.15,
    wrap = TRUE,
    sort_intersections_by = 'degree',
    queries = list(
      ComplexUpset::upset_query(set = 'EWSR1', fill = '#EB5757'),
      ComplexUpset::upset_query(set = 'FUS', fill = '#27AE60'),
      ComplexUpset::upset_query(set = 'TAF15', fill = '#4A90E2')
    ),
    themes = upset_default_themes(text = element_text(size = 20,face = "bold")),
    base_annotations = list(
      'Intersection size' = intersection_size(
        text = element_text(size = 6, vjust = -1),
        # fill = "#0072B2" #BLUE
        # fill = "#009E73" #GREEN
        fill = "#D55E00" #RED
      )
    )
  # ) + ggtitle (paste('UpSet Plot of Alternative Splicing',type_label, 'Intersections ALL' )) +
  # ) + ggtitle (paste('UpSet Plot of Alternative Splicing',type_label, 'Intersections Positive' )) +
  ) + ggtitle (paste('UpSet Plot of Alternative Splicing',type_label, 'Intersections Negative' )) +
    theme(plot.title = element_text(size = 24, face = "bold"))



}
# ---- Create a library of splicing events to gene symbols ----
# Extract relevant columns
ews_event_gene <- dplyr::select(EWSR1_splicing_table, EVENT, GENE)
fus_event_gene <- dplyr::select(FUS_splicing_table, EVENT, GENE)
taf15_event_gene <- dplyr::select(TAF15_splicing_table, EVENT, GENE)

# Combine all mappings
all_event_gene <- bind_rows(ews_event_gene, fus_event_gene, taf15_event_gene)

# Remove duplicates (if an event points to the same gene) 
event_to_gene_library <- all_event_gene %>%
  distinct(EVENT, GENE)

# Optional: if you want to prioritize non-NA or non-empty GENE names
# event_to_gene_library <- event_to_gene_library %>%
#   filter(!is.na(GENE), GENE != "")

# View the library 
head(event_to_gene_library)

# ---- Create a library of splicing events to MV_95  ----
# Extract relevant columns
ews_event_MV <- dplyr::select(EWSR1_splicing_table, EVENT, MV.dPsi._at_0.95)
fus_event_MV <- dplyr::select(FUS_splicing_table, EVENT, MV.dPsi._at_0.95)
taf15_event_MV <- dplyr::select(TAF15_splicing_table, EVENT, MV.dPsi._at_0.95)

# Combine all mappings
event_to_MV_library <- bind_rows(ews_event_MV, fus_event_MV, taf15_event_MV)

event_to_MV_library <- event_to_MV_library %>%
  filter(MV.dPsi._at_0.95 > 0.0)

# View the library 
head(event_to_MV_library)

# ---- Create a library of spilicng events to chr coordinates ----
INCLUSION_LEVELS_FULL <- read.delim("D:/FETseq/vast-tools_results/INCLUSION_LEVELS_FULL-hg38-12-v251.tab", header=FALSE)
colnames(INCLUSION_LEVELS_FULL) <- INCLUSION_LEVELS_FULL[1,]
INCLUSION_LEVELS_FULL <- INCLUSION_LEVELS_FULL[-1,]
event_to_chr_library <- dplyr::select(INCLUSION_LEVELS_FULL, EVENT, COORD)
head(event_to_chr_library)
# ---- Function to create an UpSet info Excel export ----

export_upset_info <- function(events, type_label) {
  # TESTS
  # events <- event_types[["EX"]]
  # type_label <- "EX"
  
  all_events <- unique(unlist(events))
  
  if (length(all_events) == 0) {
    message("Skipping ", type_label, " - no events found.")
    return()
  }
  
  df <- data.frame(
    EVENT = all_events,
    EWSR1 = all_events %in% events[["EWSR1"]],
    FUS   = all_events %in% events[["FUS"]],
    TAF15 = all_events %in% events[["TAF15"]]
  )
  
  # Left join to add the gene symbols
  df <- df %>%
    left_join(event_to_gene_library, by = "EVENT")
  
  # Left join to add the MV
  df <- df %>%
    left_join(event_to_MV_library, by = "EVENT")
  
  # Left join to add the chr
  df <- df %>%
    left_join(event_to_chr_library, by = "EVENT")
  
  # Make a copy (for safety)
  df_numeric <- df
  
  # Identify only the columns you want to convert (not EVENT!)
  cols_to_convert <- c("EWSR1", "FUS", "TAF15")
  
  # Convert TRUE/FALSE → 1/0
  df_numeric[cols_to_convert] <- lapply(df_numeric[cols_to_convert], as.integer)
  
  head(df_numeric)
  
  # Preserve numeric binary columns
  e <- df_numeric$EWSR1
  f <- df_numeric$FUS
  t <- df_numeric$TAF15
  event_ids <- df_numeric$EVENT
  gene_names <- df_numeric$GENE
  MV_values <- df_numeric$MV.dPsi._at_0.95
  COORD <- df_numeric$COORD
  
  # EVENTS
  binary_matrix_groups_event <- data.frame(
    EWSR1 = ifelse(e == 1 & f == 0 & t == 0, event_ids, NA),
    EWSR1_FUS = ifelse(e == 1 & f == 1 & t == 0, event_ids, NA),
    EWSR1_FUS_TAF15 = ifelse(e == 1 & f == 1 & t == 1, event_ids, NA),
    FUS = ifelse(e == 0 & f == 1 & t == 0, event_ids, NA),
    EWSR1_TAF15 = ifelse(e == 1 & f == 0 & t == 1, event_ids, NA),
    TAF15_FUS = ifelse(e == 0 & f == 1 & t == 1, event_ids, NA),
    TAF15 = ifelse(e == 0 & f == 0 & t == 1, event_ids, NA)
  )
  
  # GENES
  binary_matrix_groups_genes <- data.frame(
    EWSR1 = ifelse(e == 1 & f == 0 & t == 0, gene_names, NA),
    EWSR1_FUS = ifelse(e == 1 & f == 1 & t == 0, gene_names, NA),
    EWSR1_FUS_TAF15 = ifelse(e == 1 & f == 1 & t == 1, gene_names, NA),
    FUS = ifelse(e == 0 & f == 1 & t == 0, gene_names, NA),
    EWSR1_TAF15 = ifelse(e == 1 & f == 0 & t == 1, gene_names, NA),
    TAF15_FUS = ifelse(e == 0 & f == 1 & t == 1, gene_names, NA),
    TAF15 = ifelse(e == 0 & f == 0 & t == 1, gene_names, NA)
  )
  
  # MV
  binary_matrix_groups_MV <- data.frame(
    EWSR1 = ifelse(e == 1 & f == 0 & t == 0, MV_values, NA),
    EWSR1_FUS = ifelse(e == 1 & f == 1 & t == 0, MV_values, NA),
    EWSR1_FUS_TAF15 = ifelse(e == 1 & f == 1 & t == 1, MV_values, NA),
    FUS = ifelse(e == 0 & f == 1 & t == 0, MV_values, NA),
    EWSR1_TAF15 = ifelse(e == 1 & f == 0 & t == 1, MV_values, NA),
    TAF15_FUS = ifelse(e == 0 & f == 1 & t == 1, MV_values, NA),
    TAF15 = ifelse(e == 0 & f == 0 & t == 1, MV_values, NA)
  )
  
  # COORD
  binary_matrix_groups_COORD <- data.frame(
    EWSR1 = ifelse(e == 1 & f == 0 & t == 0, COORD, NA),
    EWSR1_FUS = ifelse(e == 1 & f == 1 & t == 0, COORD, NA),
    EWSR1_FUS_TAF15 = ifelse(e == 1 & f == 1 & t == 1, COORD, NA),
    FUS = ifelse(e == 0 & f == 1 & t == 0, COORD, NA),
    EWSR1_TAF15 = ifelse(e == 1 & f == 0 & t == 1, COORD, NA),
    TAF15_FUS = ifelse(e == 0 & f == 1 & t == 1, COORD, NA),
    TAF15 = ifelse(e == 0 & f == 0 & t == 1, COORD, NA)
  )
  
  # Clean up NA and pad
  long_df_event <- lapply(binary_matrix_groups_event, function(x) x[!is.na(x)])
  long_df_genes <- lapply(binary_matrix_groups_genes, function(x) x[!is.na(x)])
  long_df_MV <- lapply(binary_matrix_groups_MV, function(x) x[!is.na(x)])
  long_df_COORD <- lapply(binary_matrix_groups_COORD, function(x) x[!is.na(x)])
  max_len <- max(sapply(c(long_df_event, long_df_genes), length))
  long_df_event <- lapply(long_df_event, function(x) { length(x) <- max_len; return(x) })
  long_df_genes <- lapply(long_df_genes, function(x) { length(x) <- max_len; return(x) })
  long_df_MV <- lapply(long_df_MV, function(x) { length(x) <- max_len; return(x) })
  long_df_COORD <- lapply(long_df_COORD, function(x) { length(x) <- max_len; return(x) })
  
  # Combine and export
  final_upset_export_event <- as.data.frame(long_df_event)
  final_upset_export_genes <- as.data.frame(long_df_genes)
  final_upset_export_MV <- as.data.frame(long_df_MV)
  final_upset_export_COORD <- as.data.frame(long_df_COORD)
  
  wb <- createWorkbook()
  addWorksheet(wb, "EVENTS")
  addWorksheet(wb, "GENES")
  addWorksheet(wb, "MV_values")
  addWorksheet(wb, "COORD")
  writeData(wb, "EVENTS", final_upset_export_event)
  writeData(wb, "GENES", final_upset_export_genes)
  writeData(wb, "MV_values", final_upset_export_MV)
  writeData(wb, "COORD", final_upset_export_COORD)
  saveWorkbook(wb, paste0("UpSet_", type_label, "_event_list.xlsx"), overwrite = TRUE)
  
  # Create workbook
  wb <- createWorkbook()
  
  # Get group names from columns (assuming they are the same across all 4 tables)
  group_names <- colnames(final_upset_export_event)
  
  # Loop through groups and create sheets
  for (group in group_names) {
    events <- final_upset_export_event[[group]]
    genes <- final_upset_export_genes[[group]]
    mv_vals <- final_upset_export_MV[[group]]
    coords <- final_upset_export_COORD[[group]]
    
    # Filter out NA rows (if needed)
    valid_rows <- !is.na(events)
    
    group_df <- data.frame(
      EVENT = events[valid_rows],
      GENE = genes[valid_rows],
      MV_value = mv_vals[valid_rows],
      COORD = coords[valid_rows],
      stringsAsFactors = FALSE
    )
    
    # Add sheet and write data
    addWorksheet(wb, group)
    writeData(wb, group, group_df)
  }
  
  # Save workbook
  saveWorkbook(wb, paste0("UpSet_", type_label, "_by_group.xlsx"), overwrite = TRUE)
}

# ---- Run for each splicing type to save shared and non shared events ----

export_upset_info(event_types[["EX"]], "EX")
export_upset_info(event_types[["INT"]], "INT")
export_upset_info(event_types[["ALTA"]], "ALTA")
export_upset_info(event_types[["ALTD"]], "ALTD")

# ---- Run for each splicing type to create UpSet plots----

plot_upset(event_types[["EX"]], "EX")
plot_upset(event_types[["INT"]], "INT")
plot_upset(event_types[["ALTA"]], "ALTA")
plot_upset(event_types[["ALTD"]], "ALTD")








library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

plot_splicing_heatmap <- function(splicing_tables, condition_name, event_to_gene, direction = c("positive", "negative")) {
  direction <- match.arg(direction)
  
  # Define filtering function
  dpsi_filter <- if (direction == "positive") {
    function(df) df %>% filter(`E.dPsi.` > 0)
  } else {
    function(df) df %>% filter(`E.dPsi.` < 0)
  }
  
  # Filter each condition-specific splicing table
  filtered_list <- lapply(splicing_tables, dpsi_filter)
  
  # Extract unique events across all 3 conditions
  all_events <- unique(unlist(lapply(filtered_list, function(df) df$EVENT)))
  
  # Build binary matrix
  binary_matrix <- data.frame(
    EVENT = all_events,
    EWSR1 = all_events %in% filtered_list[["EWSR1"]]$EVENT,
    FUS   = all_events %in% filtered_list[["FUS"]]$EVENT,
    TAF15 = all_events %in% filtered_list[["TAF15"]]$EVENT
  )
  
  binary_long <- binary_matrix %>%
    pivot_longer(cols = -EVENT, names_to = "Condition", values_to = "Present")
  
  # Create zScore table in long format
  z_long <- bind_rows(
    filtered_list[["EWSR1"]] %>% select(EVENT, `E.dPsi.`) %>% mutate(Condition = "EWSR1"),
    filtered_list[["FUS"]]   %>% select(EVENT, `E.dPsi.`) %>% mutate(Condition = "FUS"),
    filtered_list[["TAF15"]] %>% select(EVENT, `E.dPsi.`) %>% mutate(Condition = "TAF15")
  )
  colnames(z_long)[2] <- "zScore"
  
  # Merge zScores and gene names
  heatmap_df <- binary_long %>%
    left_join(z_long, by = c("EVENT", "Condition")) %>%
    mutate(zScore = ifelse(Present == 1, zScore, NA)) %>%
    left_join(event_to_gene, by = "EVENT") %>%
    mutate(Label = ifelse(!is.na(GENE), paste0(GENE, " (", EVENT, ")"), EVENT))
  
  heatmap_df$Label <- fct_rev(factor(heatmap_df$Label, levels = unique(heatmap_df$Label)))
  
  # Plot
  p <- ggplot(heatmap_df, aes(x = Condition, y = Label, fill = zScore)) +
    geom_tile(color = "white") +
    geom_text(aes(label = ifelse(!is.na(zScore), sprintf("%.2f", zScore), "")), size = 3.5) +
    scale_fill_gradient2(
      low = "#D55E00", mid = "white", high = "#009E73", midpoint = 0, na.value = "grey90",
      name = "ΔPSI"
    ) +
    labs(
      title = paste("ΔPSI for", condition_name, "-", direction),
      x = "Condition", y = "Splicing Event"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
  
  return(p)
}

# ---- Create a library of splicing events to gene symbols ----
ews_event_gene <- dplyr::select(EWSR1_splicing_table, EVENT, GENE)
fus_event_gene <- dplyr::select(FUS_splicing_table, EVENT, GENE)
taf15_event_gene <- dplyr::select(TAF15_splicing_table, EVENT, GENE)

# Combine all mappings
all_event_gene <- bind_rows(ews_event_gene, fus_event_gene, taf15_event_gene)

# Remove duplicates (if an event points to the same gene) 
event_to_gene_library <- all_event_gene %>%
  distinct(EVENT, GENE)
splicing_lists <- list(
  EX = list(
    EWSR1 = EWSR1_filtered_EX,
    FUS   = FUS_filtered_EX,
    TAF15 = TAF15_filtered_EX
  ),
  INT = list(
    EWSR1 = EWSR1_filtered_INT,
    FUS   = FUS_filtered_INT,
    TAF15 = TAF15_filtered_INT
  ),
  ALTA = list(
    EWSR1 = EWSR1_filtered_ALTA,
    FUS   = FUS_filtered_ALTA,
    TAF15 = TAF15_filtered_ALTA
  ),
  ALTD = list(
    EWSR1 = EWSR1_filtered_ALTD,
    FUS   = FUS_filtered_ALTD,
    TAF15 = TAF15_filtered_ALTD
  )
)

library(ggplot2)

pdf("All_Splicing_Heatmaps.pdf", width = 12, height = 10)

for (type_name in names(splicing_lists)) {
  for (direction in c("positive", "negative")) {
    p <- plot_splicing_heatmap(
      splicing_tables = splicing_lists[[type_name]],
      condition_name = type_name,
      event_to_gene = event_to_gene_library,
      direction = direction
    )
    print(p)
  }
}

dev.off()

