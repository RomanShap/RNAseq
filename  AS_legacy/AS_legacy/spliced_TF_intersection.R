
# OLD VERSION
# TAKES INTO ACCOUNT ONLY MV>0.06





# ---- Load required libraries ----
library(readxl)
library(dplyr)
library(tidyr)
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

# ---- Split splicing tables by event type ----
split_splicing_by_event_type <- function(splicing_tables) {
  lapply(splicing_tables, function(df) {
    if (!"EVENT" %in% names(df)) stop("Missing 'EVENT' column")
    df$Event_Type <- sub("^Hsa([A-Z]+).*", "\\1", df$EVENT)
    split(df, df$Event_Type)
  })
}

# ---- Prepare splicing input for all/pos/neg ----
# Assumes EWSR1_filtered, FUS_filtered, TAF15_filtered are pre-existing named lists
splicing_filtered <- list(
  all = list(EWSR1 = EWSR1_filtered$all, FUS = FUS_filtered$all, TAF15 = TAF15_filtered$all),
  pos = list(EWSR1 = EWSR1_filtered$pos, FUS = FUS_filtered$pos, TAF15 = TAF15_filtered$pos),
  neg = list(EWSR1 = EWSR1_filtered$neg, FUS = FUS_filtered$neg, TAF15 = TAF15_filtered$neg)
)

# Split splicing data into event types
splicing_by_event <- lapply(splicing_filtered, split_splicing_by_event_type)

# ---- Compare spliced genes to ISMARA TF targets ----

# List of filters, samples, and splicing event types
filter_types <- c("all", "pos", "neg")
samples <- c("EWSR1", "FUS", "TAF15")
event_types <- c("EX", "INT", "ALTA", "ALTD")

# Prepare output list for storing matched transcription factors
spliced_motifs <- list()

# Nested loop to compare splicing events with ISMARA TFs
for (filter_type in filter_types) {
  for (sample in samples) {
    for (event_type in event_types) {
      
      # Safely extract splicing dataframe for given filter/sample/event
      splicing_df <- tryCatch(
        splicing_by_event[[filter_type]][[sample]][[event_type]],
        error = function(e) NULL
      )
      if (is.null(splicing_df)) next
      
      # Get list of TFs from ISMARA for the current sample
      tf_list <- ISMARA_separated[[sample]][["Transcription_Factor"]]
      
      # Find overlap between spliced genes and TF targets
      matched <- intersect(splicing_df[["GENE"]], tf_list)
      
      # Store matched genes
      spliced_motifs[[filter_type]][[sample]][[event_type]] <- matched
    }
  }
}

# ---- Flatten results into a data frame ----
flat_spliced_motifs <- expand.grid(
  Filter = filter_types,
  Sample = samples,
  Event  = event_types,
  stringsAsFactors = FALSE
) %>%
  rowwise() %>%
  mutate(
    Genes = list(spliced_motifs[[Filter]][[Sample]][[Event]] %||% character(0))
  ) %>%
  unnest(cols = Genes)

# ---- Save results to Excel ----
write_xlsx(flat_spliced_motifs, path = "spliced_motifs_intersections.xlsx")

library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(purrr)

# Define your input objects
samples <- c("EWSR1", "FUS", "TAF15")
event_types <- c("EX", "INT", "ALTA", "ALTD")

get_long_dpsi <- function(sample_name, sample_filtered_list) {
  lapply(names(sample_filtered_list), function(filter_type) {
    df <- sample_filtered_list[[filter_type]]
    if (is.null(df) || !"GENE" %in% names(df) || !"E.dPsi." %in% names(df) || !"EVENT" %in% names(df)) {
      return(NULL)
    }
    
    df %>%
      mutate(
        Sample = sample_name,
        Filter = filter_type,
        Event = sub("^Hsa([A-Z]+).*", "\\1", EVENT)  # extract EX, ALTA, etc
      ) %>%
      dplyr::select(GENE, E.dPsi., Sample, Filter, Event)
  }) %>%
    bind_rows()
}


long_dpsi_df <- bind_rows(
  get_long_dpsi("EWSR1", EWSR1_filtered),
  get_long_dpsi("FUS", FUS_filtered),
  get_long_dpsi("TAF15", TAF15_filtered)
) %>%
  dplyr::rename(Gene = GENE, E.dPsi = `E.dPsi.`)

# 2. Filter flat_spliced_motifs to include only unique TF/Sample/Event combinations
heatmap_input <- flat_spliced_motifs %>%
  dplyr::rename(Gene = Genes) %>%
  dplyr::filter(Event %in% event_types, Sample %in% samples)

# # 3. Join with long-format E.dPsi values
# heatmap_data <- heatmap_input %>%
#   left_join(long_dpsi_df, by = c("Gene", "Sample", "Event"))
# 
# heatmap_data_collapsed <- heatmap_data %>%
#   group_by(Gene, Sample, Event) %>%
#   slice_max(order_by = abs(E.dPsi), n = 1, with_ties = FALSE) %>%
#   ungroup()
# 
# # 4. Loop to plot heatmaps by event type
# for (event_type in event_types) {
#   plot_df <- heatmap_data_collapsed %>%
#     filter(Event == event_type) %>%
#     mutate(Gene = fct_rev(factor(Gene, levels = unique(Gene))))
#   
#   p <- ggplot(plot_df, aes(x = Sample, y = Gene, fill = E.dPsi)) +
#     geom_tile(color = "white") +
#     geom_text(aes(label = ifelse(!is.na(E.dPsi), sprintf("%.2f", E.dPsi), "")), size = 3) +
#     scale_fill_gradient2(
#       low = "#d73027", mid = "white", high = "#1a9850", midpoint = 0, na.value = "grey90",
#       name = "E.dPsi"
#     ) +
#     labs(
#       title = paste("Spliced TFs with ISMARA Motifs — Event:", event_type),
#       x = "Sample", y = "Transcription Factor"
#     ) +
#     theme_minimal(base_size = 14) +
#     theme(
#       axis.text.x = element_text(angle = 45, hjust = 1),
#       panel.grid = element_blank()
#     )
#   
#   print(p)
# }

# 3. Join with long-format E.dPsi values
heatmap_data <- heatmap_input %>%
  left_join(long_dpsi_df, by = c("Gene", "Sample", "Event"))
heatmap_data_collapsed <- heatmap_data %>%
  group_by(Gene, Sample, Event) %>%
  summarise(E.dPsi = mean(E.dPsi, na.rm = TRUE), .groups = "drop") %>%
  ungroup()

# # 4. Loop to plot and save heatmaps by event type
# for (event_type in event_types) {
#   plot_df <- heatmap_data_collapsed %>%
#     filter(Event == event_type) %>%
#     mutate(Gene = fct_rev(factor(Gene, levels = unique(Gene))))
#   
#   p <- ggplot(plot_df, aes(x = Sample, y = Gene, fill = E.dPsi)) +
#     geom_tile(color = "white") +
#     geom_text(aes(label = ifelse(!is.na(E.dPsi), sprintf("%.2f", E.dPsi), "")), size = 3) +
#     scale_fill_gradient2(
#       low = "#1E88E5", mid = "white", high = "#FFC107", midpoint = 0, na.value = "grey90",
#       name = "E.dPsi"
#     ) +
#     labs(
#       title = paste("Spliced TFs with ISMARA Motifs — Event:", event_type),
#       x = "Sample", y = "Transcription Factor"
#     ) +
#     theme_minimal(base_size = 14) +
#     theme(
#       axis.text.x = element_text(angle = 45, hjust = 1),
#       panel.grid = element_blank()
#     )
#   
#   # Save plot to PNG
#   ggsave(
#     filename = paste0("heatmap_ISMARA_spliced_TFs_", event_type, ".png"),
#     plot = p,
#     width = 8,
#     height = max(4, nlevels(plot_df$Gene) * 0.3),  # adapt height to gene count
#     dpi = 300
#   )
#   
#   print(p)
# }
for (event_type in event_types) {
  plot_df <- heatmap_data_collapsed %>%
    filter(Event == event_type) %>%
    mutate(
      Gene = fct_rev(factor(Gene, levels = unique(Gene))),
      Sample = factor(Sample, levels = c("FUS", "EWSR1", "TAF15"))
    )
  
  p <- ggplot(plot_df, aes(x = Sample, y = Gene, fill = E.dPsi)) +
    geom_tile(color = "white") +
    geom_text(aes(label = ifelse(!is.na(E.dPsi), sprintf("%.2f", E.dPsi), "")), size = 3) +
    scale_fill_gradient2(
      low = "#1E88E5", mid = "white", high = "#FFC107", midpoint = 0, na.value = "grey90",
      name = "E.dPsi"
    ) +
    labs(
      title = paste("Spliced TFs with ISMARA Motifs — Event:", event_type),
      x = "Sample", y = "Transcription Factor"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
  
  ggsave(
    filename = paste0("heatmap_ISMARA_spliced_TFs_", event_type, ".png"),
    plot = p,
    width = 8,
    height = max(4, nlevels(plot_df$Gene) * 0.3),
    dpi = 300
  )
  
  print(p)
}

