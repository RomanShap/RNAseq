# Here we are searching for the alternatively spliced TF taking ISMARA table,
# filtering it by abs(z)>1.5 and searching resulting genes in the AS table from vast tools. 
# ---- Clear the Environment ----
rm(list = ls())            # Remove all objects
graphics.off()             # Close all open plots
cat("\014")                # Clear console (works in RStudio)
# ---- Load required libraries ----
library(readxl)
library(tidyr)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)
library(writexl)

# ---- Import and filter ISMARA results ----
ISMARA_dir <- "D:/FETseq/FETseq-analysis/ISMARA"

# Read ISMARA differential result files for each sample
ISMARA_differential_result_list <- list(
  EWSR1 = read_excel(file.path(ISMARA_dir, "differential_results_EWSR1.xlsx")),
  FUS   = read_excel(file.path(ISMARA_dir, "differential_results_FUS.xlsx")),
  TAF15 = read_excel(file.path(ISMARA_dir, "differential_results_TAF15.xlsx"))
)

# Filter motifs with |z| ≥ 1.5
ISMARA_filtered <- lapply(ISMARA_differential_result_list, function(df) {
  filter(df, abs(z) >= 1.5)
})

# Separate multiple transcription factors (TFs) listed in the same row
split_tf_rows <- function(df) {
  separate_rows(df, Transcription_Factor, sep = "_")
}

ISMARA_separated <- lapply(ISMARA_filtered, split_tf_rows)
# ---- Import splicing tables ----
vast_tools_dir <- "D:/FETseq/vast-tools_results"
sample_id <- dir(vast_tools_dir)
sample_id <- sample_id[!grepl("INCLUSION_LEVELS", sample_id)]
splicing_tab_paths <- file.path(vast_tools_dir, sample_id)
names(splicing_tab_paths) <- sample_id

inclusion_levels_table <- read.table(file.path(vast_tools_dir, "INCLUSION_LEVELS_FULL-hg38-12-v251.tab"),
                                     header = TRUE,       
                                     sep = "\t",          
                                     stringsAsFactors = FALSE,
                                     quote = "",          
                                     comment.char = "")   

splicing_data_list <- lapply(splicing_tab_paths, function(path) {
  read.delim(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
})

EWSR1_splicing_table <- splicing_data_list[["diff_EWSR1KO_vs_WT.tab"]]
FUS_splicing_table    <- splicing_data_list[["diff_FUSKO_vs_WT.tab"]]
TAF15_splicing_table  <- splicing_data_list[["diff_TAF15KO_vs_WT.tab"]]

# --- helpers ---

# Extract event class from EVENT ("HsaEX...", "HsaINT...", "HsaALTA...", "HsaALTD...")
.event_class <- function(evt) {
  cls <- str_extract(evt, "(?<=^Hsa)[A-Z]+")
  ifelse(cls %in% c("EX","INT","ALTA","ALTD"), cls, NA_character_)
}

# Split a filtered DF into the four event classes (always return 4 tibbles)
.split_by_event <- function(df) {
  df <- df %>% mutate(.event = .event_class(EVENT))
  set_names(c("EX","INT","ALTA","ALTD")) |>
    map(\(k) df %>% filter(.event == k) %>% select(-.event))
}

# Add COORD to a splicing table using GENE + EVENT
.add_coords <- function(splicing_df, inclusion_levels_table) {
  coords <- inclusion_levels_table %>%
    select(GENE, EVENT, COORD) %>%
    distinct()   # protect against accidental duplicates
  
  as_tibble(splicing_df) %>%
    left_join(coords, by = c("GENE", "EVENT"))
}

# --- main ---

subset_splicing_by_ISMARA <- function(ISMARA_separated,
                                      EWSR1_splicing_table,
                                      FUS_splicing_table,
                                      TAF15_splicing_table,
                                      inclusion_levels_table) {
  # 1) Enrich each splicing table with COORD
  EWS_tbl <- .add_coords(EWSR1_splicing_table, inclusion_levels_table)
  FUS_tbl <- .add_coords(FUS_splicing_table,   inclusion_levels_table)
  TAF_tbl <- .add_coords(TAF15_splicing_table, inclusion_levels_table)
  
  # 2) ISMARA TF sets per condition
  ismara_genes <- list(
    EWS = ISMARA_separated$EWSR1 %>% pull(Transcription_Factor) %>% unique(),
    FUS = ISMARA_separated$FUS   %>% pull(Transcription_Factor) %>% unique(),
    TAF = ISMARA_separated$TAF15 %>% pull(Transcription_Factor) %>% unique()
  )
  
  # 3) Bundle for iteration
  splice_tbls <- list(EWS = EWS_tbl, FUS = FUS_tbl, TAF = TAF_tbl)
  
  # 4) For each condition: subset by ISMARA genes, then split by event class
  imap(splice_tbls, function(df, key) {
    df %>%
      filter(GENE %in% ismara_genes[[key]]) %>%
      .split_by_event()
  })
}

# --- Example run (objects already in memory) ---
result <- subset_splicing_by_ISMARA(
  ISMARA_separated,
  EWSR1_splicing_table,
  FUS_splicing_table,
  TAF15_splicing_table,
  inclusion_levels_table
)

write_splicing_result_to_xlsx <- function(result,
                                          folder = "splised_TF_results",
                                          file = "splicing_ISMARA_subsets.xlsx") {
  # Create folder if it doesn't exist
  if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
  
  # Full path to Excel file
  filepath <- file.path(folder, file)
  
  # Flatten nested list (EWS/FUS/TAF -> EX/INT/ALTA/ALTD)
  sheets_list <-
    imap(result, function(event_list, condition) {
      imap(event_list, function(df, evt) {
        df %>% arrange(GENE, EVENT)
      }) %>%
        set_names(~ paste(condition, .x, sep = "_"))
    }) %>%
    unlist(recursive = FALSE)
  
  # Guarantee non-empty sheets for Excel
  sheets_list <- imap(sheets_list, function(df, nm) {
    if (nrow(df) == 0) tibble(note = "No rows after filtering") else df
  })
  
  # Write Excel file
  write_xlsx(sheets_list, path = filepath)
  message("Wrote: ", normalizePath(filepath))
  invisible(filepath)
}

# ---- Run it ----
xlsx_path <- write_splicing_result_to_xlsx(result)

