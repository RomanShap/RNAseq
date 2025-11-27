# Clear the Environment ----
rm(list = ls())            # Remove all objects
graphics.off()             # Close all open plots
cat("\014")                # Clear console (works in RStudio)
# ---- Load required libraries ----
suppressMessages({library(dplyr)
  library(readxl)
  library(stringr)
  library(purrr)
  library(openxlsx)
  library(tidyr)
  })

# ---- Load data ----
SRRM3_targets <- read_excel("Juan-Mateu_NatMet2023_STable3.xlsx")
SRRM3_mouse_targets <- read_excel("Nakano_CellRep2019_TableS1.xlsx")
FET_EX_non_filtered <- read_excel("Shared_Events/EX.xlsx")
# FET_EX_filtered <- read_excel("Shared_filtered_Events/EX.xlsx")
FET_EX_non_filtered <- read_excel("Shared_Events/EX.xlsx", 
                 sheet = "FUS_EWSR1_TAF15")
FE_EX_non_filtered <- read_excel("Shared_Events/EX.xlsx", 
                                  sheet = "FUS_EWSR1")
FT_EX_non_filtered <- read_excel("Shared_Events/EX.xlsx", 
                                  sheet = "FUS_TAF15")
ET_EX_non_filtered <- read_excel("Shared_Events/EX.xlsx", 
                                  sheet = "EWSR1_TAF15")
F_EX_non_filtered <- read_excel("Shared_Events/EX.xlsx", 
                                  sheet = "FUS")
E_EX_non_filtered <- read_excel("Shared_Events/EX.xlsx", 
                                  sheet = "EWSR1")
T_EX_non_filtered <- read_excel("Shared_Events/EX.xlsx", 
                                  sheet = "TAF15")

# ---- Filter FET table for SRRM3 target genes ----
# 1) Extract SRRM3 gene list (second column; drop the header row "GENE")
srrm3_genes <- SRRM3_targets[[2]] %>%
  as.character() %>%
  discard(is.na) %>%
  str_trim() %>%
  setdiff("GENE") %>%     # remove the header text present as a row
  unique()

srrm3_mouse_symbols <- SRRM3_mouse_targets[[1]] %>%
  as.character() %>%
  replace_na("") %>%
  str_trim() %>%
  # drop header row(s)
  discard(~ .x == "" || str_to_upper(.x) == "GENE SYMBOL (COMMENT)") %>%
  # split "Fath, Fat1" -> c("Fath", "Fat1")
  str_split(",\\s*") %>%
  flatten_chr() %>%
  # remove any trailing comments like "GENE1 (alt)"
  str_replace("\\s*\\(.*\\)$", "") %>%
  str_trim() %>%
  str_to_upper() %>%
  unique()

# 2) Filter your FET table by those genes
# FET_EX_SRRM3 <- FET_EX_filtered %>%
  # filter(GENE %in% srrm3_genes)
FET_EX_SRRM3_nf <- FET_EX_non_filtered %>%
  filter(GENE %in% srrm3_genes)
FE_EX_SRRM3_nf <- FE_EX_non_filtered %>%
  filter(GENE %in% srrm3_genes)
FT_EX_SRRM3_nf <- FT_EX_non_filtered %>%
  filter(GENE %in% srrm3_genes)
ET_EX_SRRM3_nf <- ET_EX_non_filtered %>%
  filter(GENE %in% srrm3_genes)
F_EX_SRRM3_nf <- F_EX_non_filtered %>%
  filter(GENE %in% srrm3_genes)
E_EX_SRRM3_nf <- E_EX_non_filtered %>%
  filter(GENE %in% srrm3_genes)
T_EX_SRRM3_nf <- T_EX_non_filtered %>%
  filter(GENE %in% srrm3_genes)
# ---- Filter FET table for mouse SRRM3 target genes ----
FET_EX_SRRM3_mouse_nf <- FET_EX_non_filtered %>%
  filter(GENE %in% srrm3_mouse_symbols)
FE_EX_SRRM3_mouse_nf <- FE_EX_non_filtered %>%
  filter(GENE %in% srrm3_mouse_symbols)
FT_EX_SRRM3_mouse_nf <- FT_EX_non_filtered %>%
  filter(GENE %in% srrm3_mouse_symbols)
ET_EX_SRRM3_mouse_nf <- ET_EX_non_filtered %>%
  filter(GENE %in% srrm3_mouse_symbols)
F_EX_SRRM3_mouse_nf <- F_EX_non_filtered %>%
  filter(GENE %in% srrm3_mouse_symbols)
E_EX_SRRM3_mouse_nf <- E_EX_non_filtered %>%
  filter(GENE %in% srrm3_mouse_symbols)
T_EX_SRRM3_mouse_nf <- T_EX_non_filtered %>%
  filter(GENE %in% srrm3_mouse_symbols)

# ---- Save the filtered table ----
# ensure folder exists
dir.create("SRRM3_filtered_events", showWarnings = FALSE, recursive = TRUE)

# collect all tables into a named list (names = sheet names)
sheets <- list(
  FET_EX = FET_EX_SRRM3_nf,
  FE_EX  = FE_EX_SRRM3_nf,
  FT_EX  = FT_EX_SRRM3_nf,
  ET_EX  = ET_EX_SRRM3_nf,
  F_EX   = F_EX_SRRM3_nf,
  E_EX   = E_EX_SRRM3_nf,
  T_EX   = T_EX_SRRM3_nf
)

sheets_mouse <- list(
  FET_EX_mouse = FET_EX_SRRM3_mouse_nf,
  FE_EX_mouse  = FE_EX_SRRM3_mouse_nf,
  FT_EX_mouse  = FT_EX_SRRM3_mouse_nf,
  ET_EX_mouse  = ET_EX_SRRM3_mouse_nf,
  F_EX_mouse   = F_EX_SRRM3_mouse_nf,
  E_EX_mouse   = E_EX_SRRM3_mouse_nf,
  T_EX_mouse   = T_EX_SRRM3_mouse_nf
)



# (optional) drop empty or non-data-frame entries
sheets <- sheets[vapply(sheets, function(x) is.data.frame(x) && nrow(x) > 0, logical(1))]
sheets_mouse <- sheets_mouse[vapply(sheets_mouse, function(x) is.data.frame(x) && nrow(x) > 0, logical(1))]
# write one Excel with multiple sheets
write.xlsx(
  x = sheets,
  file = "SRRM3_filtered_events/SRRM3_filtered_events.xlsx",
  overwrite = TRUE
)
write.xlsx(
  x = sheets_mouse,
  file = "SRRM3_filtered_events/SRRM3_filtered_events_mouse_symbols.xlsx",
  overwrite = TRUE
)
