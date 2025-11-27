# Clear the Environment ----
rm(list = ls())            # Remove all objects
graphics.off()             # Close all open plots
cat("\014")                # Clear console (works in RStudio)
# ---- Load required libraries ----
suppressMessages({library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(UpSetR)
library(ComplexUpset)
library(openxlsx)
library(scales)})

# ---- Set base directory and read all splicing tables ----
base_dir <- "D:/FETseq/vast-tools_results"
sample_id <- dir(base_dir)
sample_id <- sample_id[!grepl("INCLUSION_LEVELS", sample_id)]
splicing_tab_paths <- file.path(base_dir, sample_id)
names(splicing_tab_paths) <- sample_id

splicing_data_list <- lapply(splicing_tab_paths, function(path) {
  read.delim(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
})

inclusion_levels_table <- read.table(file.path(base_dir, "INCLUSION_LEVELS_FULL-hg38-12-v251.tab"),
                   header = TRUE,       
                   sep = "\t",          
                   stringsAsFactors = FALSE,
                   quote = "",          
                   comment.char = "")   

EWSR1_splicing_table <- splicing_data_list[["diff_EWSR1KO_vs_WT.tab"]]
FUS_splicing_table    <- splicing_data_list[["diff_FUSKO_vs_WT.tab"]]
TAF15_splicing_table  <- splicing_data_list[["diff_TAF15KO_vs_WT.tab"]]

# Keep only the columns needed from inclusion levels
coords <- inclusion_levels_table[, c("EVENT", "COORD")]

# Left-join by EVENT only (preserves your GENE)
EWSR1_splicing_table <- merge(EWSR1_splicing_table, coords, by = "EVENT", all.x = TRUE)
FUS_splicing_table   <- merge(FUS_splicing_table,   coords, by = "EVENT", all.x = TRUE)
TAF15_splicing_table <- merge(TAF15_splicing_table, coords, by = "EVENT", all.x = TRUE)


# ---- Extract splicing TYPE from EVENT ID ----
extract_type <- function(event_id) {
  sub("^Hsa([A-Z]+).*", "\\1", event_id)
}

EWSR1_splicing_table$TYPE <- extract_type(EWSR1_splicing_table$EVENT)
FUS_splicing_table$TYPE   <- extract_type(FUS_splicing_table$EVENT)
TAF15_splicing_table$TYPE <- extract_type(TAF15_splicing_table$EVENT)

# ---- Split by splicing type ----
EWSR1_split <- split(EWSR1_splicing_table, EWSR1_splicing_table$TYPE)
FUS_split   <- split(FUS_splicing_table,   FUS_splicing_table$TYPE)
TAF15_split <- split(TAF15_splicing_table, TAF15_splicing_table$TYPE)

list_of_events <- list(
  EWSR1 = EWSR1_split,
  FUS = FUS_split,
  TAF15 = TAF15_split
)

# Function to compute all pairwise and triple intersections, and unique events
get_event_partitions <- function(event_df_list,column) {
  # Extract EVENT vectors per sample
  e <- event_df_list$EWSR1[[column]]
  f <- event_df_list$FUS[[column]]
  t <- event_df_list$TAF15[[column]]
  
  # Return named list of partitions
  list(
    FUS_EWSR1_TAF15 = intersect(intersect(f, e), t),
    FUS_EWSR1       = setdiff(intersect(f, e), t),
    FUS_TAF15       = setdiff(intersect(f, t), e),
    EWSR1_TAF15     = setdiff(intersect(e, t), f),
    FUS             = setdiff(f, union(e, t)),
    EWSR1           = setdiff(e, union(f, t)),
    TAF15           = setdiff(t, union(f, e))
  )
}

# Wrapper across all event types, with column as parameter
get_all_event_partitions <- function(list_of_events, column = "EVENT") {
  event_types <- intersect(
    intersect(names(list_of_events$EWSR1), names(list_of_events$FUS)),
    names(list_of_events$TAF15)
  )
  
  result <- lapply(event_types, function(type) {
    per_type_list <- lapply(list_of_events, function(x) x[[type]])
    get_event_partitions(per_type_list, column)
  })
  
  names(result) <- event_types
  return(result)
}

# Run the partitioning
partitioned_events <- get_all_event_partitions(list_of_events,column = "EVENT")

library(openxlsx)
library(dplyr)

# # Helper function to extract desired columns from a source table for a set of events
get_event_data <- function(df, events, prefix) {
  df %>%
    filter(EVENT %in% events) %>%
    select(EVENT, GENE, `E.dPsi.`, `MV.dPsi._at_0.95`, `COORD`, `TYPE`) %>%
    rename_with(~ paste0(prefix, "_", .), .cols = c(`E.dPsi.`, `MV.dPsi._at_0.95`))
}

# get_event_data <- function(df, events, prefix) {
#   # Ensure required columns exist
#   req <- c("EVENT", "GENE", "COORD", "TYPE", "E.dPsi.", "MV.dPsi._at_0.95")
#   miss <- setdiff(req, names(df))
#   if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))
#   
#   # Deduplicate and coerce once (important for speed)
#   events <- unique(as.character(events))
#   
#   # Fast membership via match()
#   keep <- match(df$EVENT, events, nomatch = 0L) > 0L
#   
#   out <- df[keep, req]
#   
#   # Prefix only the two dPsi columns
#   rename_map <- setNames(
#     paste0(prefix, "_", c("E.dPsi.", "MV.dPsi._at_0.95")),
#     c("E.dPsi.", "MV.dPsi._at_0.95")
#   )
#   names(out)[match(names(rename_map), names(out))] <- unname(rename_map)
#   
#   out
# }

# Main function to write Excel files for all event types and partitions
write_partitioned_event_tables <- function(partitioned_events, out_dir = ".", 
                                           EWSR1_df, FUS_df, TAF15_df) {
  
  for (event_type in names(partitioned_events)) {
    # Create a new workbook
    wb <- createWorkbook()
    
    # Loop over all 7 partitions
    for (group in names(partitioned_events[[event_type]])) {
      events <- partitioned_events[[event_type]][[group]]
      
      # Skip if no events
      if (length(events) == 0) next
      
      # Determine which samples are relevant for this group
      samples <- unlist(strsplit(group, "_"))
      
      # Initialize an empty list to hold partial data frames
      group_dfs <- list()
      
      # For each sample in this group, get its table and filtered rows
      if ("FUS" %in% samples) {
        group_dfs[["FUS"]] <- get_event_data(FUS_df, events, "FUS")
      }
      if ("EWSR1" %in% samples) {
        group_dfs[["EWSR1"]] <- get_event_data(EWSR1_df, events, "EWSR1")
      }
      if ("TAF15" %in% samples) {
        group_dfs[["TAF15"]] <- get_event_data(TAF15_df, events, "TAF15")
      }
      
      # Join all data frames by EVENT and GENE (same for all samples)
      merged_df <- Reduce(function(x, y) full_join(x, y, by = c("EVENT", "GENE")), group_dfs)
      
      # Add sheet to workbook
      addWorksheet(wb, sheetName = group)
      writeData(wb, sheet = group, merged_df)
    }
    
    # Save workbook
    saveWorkbook(wb, file = file.path(out_dir, paste0(event_type, ".xlsx")), overwrite = TRUE)
  }
}

write_partitioned_event_tables(
  partitioned_events = partitioned_events,
  out_dir = "Shared_Events",  
  EWSR1_df = EWSR1_splicing_table,
  FUS_df   = FUS_splicing_table,
  TAF15_df = TAF15_splicing_table
)
# ---- filterd tables ----

filter <- 0.05
EWSR1_splicing_table_filtered01 <- EWSR1_splicing_table %>%
  filter(`MV.dPsi._at_0.95` > filter)
FUS_splicing_table_filtered01 <- FUS_splicing_table %>%
  filter(`MV.dPsi._at_0.95` > filter)
TAF15_splicing_table_filtered01 <- TAF15_splicing_table %>%
  filter(`MV.dPsi._at_0.95` > filter)
# ---- Extract splicing TYPE from EVENT ID ----
extract_type <- function(event_id) {
  sub("^Hsa([A-Z]+).*", "\\1", event_id)
}

EWSR1_splicing_table_filtered01$TYPE <- extract_type(EWSR1_splicing_table_filtered01$EVENT)
FUS_splicing_table_filtered01$TYPE   <- extract_type(FUS_splicing_table_filtered01$EVENT)
TAF15_splicing_table_filtered01$TYPE <- extract_type(TAF15_splicing_table_filtered01$EVENT)

# ---- Split by splicing type ----
EWSR1_split <- split(EWSR1_splicing_table_filtered01, EWSR1_splicing_table_filtered01$TYPE)
FUS_split   <- split(FUS_splicing_table_filtered01,   FUS_splicing_table_filtered01$TYPE)
TAF15_split <- split(TAF15_splicing_table_filtered01, TAF15_splicing_table_filtered01$TYPE)

list_of_events <- list(
  EWSR1 = EWSR1_split,
  FUS = FUS_split,
  TAF15 = TAF15_split
)

# Function to compute all pairwise and triple intersections, and unique events
get_event_partitions <- function(event_df_list,column) {
  # Extract EVENT vectors per sample
  e <- event_df_list$EWSR1[[column]]
  f <- event_df_list$FUS[[column]]
  t <- event_df_list$TAF15[[column]]
  
  # Return named list of partitions
  list(
    FUS_EWSR1_TAF15 = intersect(intersect(f, e), t),
    FUS_EWSR1       = setdiff(intersect(f, e), t),
    FUS_TAF15       = setdiff(intersect(f, t), e),
    EWSR1_TAF15     = setdiff(intersect(e, t), f),
    FUS             = setdiff(f, union(e, t)),
    EWSR1           = setdiff(e, union(f, t)),
    TAF15           = setdiff(t, union(f, e))
  )
}

# Wrapper across all event types, with column as parameter
get_all_event_partitions <- function(list_of_events, column = "EVENT") {
  event_types <- intersect(
    intersect(names(list_of_events$EWSR1), names(list_of_events$FUS)),
    names(list_of_events$TAF15)
  )
  
  result <- lapply(event_types, function(type) {
    per_type_list <- lapply(list_of_events, function(x) x[[type]])
    get_event_partitions(per_type_list, column)
  })
  
  names(result) <- event_types
  return(result)
}

# Run the partitioning
partitioned_events <- get_all_event_partitions(list_of_events,column = "EVENT")

# Helper function to extract desired columns from a source table for a set of events
get_event_data <- function(df, events, prefix) {
  df %>%
    filter(EVENT %in% events) %>%
    select(EVENT, GENE, `E.dPsi.`, `MV.dPsi._at_0.95`) %>%
    rename_with(~ paste0(prefix, "_", .), .cols = c(`E.dPsi.`, `MV.dPsi._at_0.95`))
}

# Main function to write Excel files for all event types and partitions
write_partitioned_event_tables <- function(partitioned_events, out_dir = ".", 
                                           EWSR1_df, FUS_df, TAF15_df) {
  
  for (event_type in names(partitioned_events)) {
    # Create a new workbook
    wb <- createWorkbook()
    
    # Loop over all 7 partitions
    for (group in names(partitioned_events[[event_type]])) {
      events <- partitioned_events[[event_type]][[group]]
      
      # Skip if no events
      if (length(events) == 0) next
      
      # Determine which samples are relevant for this group
      samples <- unlist(strsplit(group, "_"))
      
      # Initialize an empty list to hold partial data frames
      group_dfs <- list()
      
      # For each sample in this group, get its table and filtered rows
      if ("FUS" %in% samples) {
        group_dfs[["FUS"]] <- get_event_data(FUS_df, events, "FUS")
      }
      if ("EWSR1" %in% samples) {
        group_dfs[["EWSR1"]] <- get_event_data(EWSR1_df, events, "EWSR1")
      }
      if ("TAF15" %in% samples) {
        group_dfs[["TAF15"]] <- get_event_data(TAF15_df, events, "TAF15")
      }
      
      # Join all data frames by EVENT and GENE (same for all samples)
      merged_df <- Reduce(function(x, y) full_join(x, y, by = c("EVENT", "GENE")), group_dfs)
      
      # Add sheet to workbook
      addWorksheet(wb, sheetName = group)
      writeData(wb, sheet = group, merged_df)
    }
    
    # Save workbook
    saveWorkbook(wb, file = file.path(out_dir, paste0(event_type, ".xlsx")), overwrite = TRUE)
  }
}

write_partitioned_event_tables(
  partitioned_events = partitioned_events,
  out_dir = "Shared_filtered_Events",  
  EWSR1_df = EWSR1_splicing_table,
  FUS_df   = FUS_splicing_table,
  TAF15_df = TAF15_splicing_table
)
