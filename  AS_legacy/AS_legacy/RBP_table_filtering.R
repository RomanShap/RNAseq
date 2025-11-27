# ---- Load required libraries ----
library(readxl)
library(dplyr)
library(ComplexUpset)
library(ggplot2) 
library(writexl)
# ---- Import and shape RBP list ----
downloads_dir <- "C:/Users/Roman/Downloads"
RBP_Information_all_motifs <- read_excel(file.path(downloads_dir,"RBP_Information_all_motifs.xlsx"))
head(RBP_Information_all_motifs)
RBP_Information_human <- RBP_Information_all_motifs %>%
  filter(RBP_Species == "Homo_sapiens")
write_xlsx(RBP_Information_human, path = file.path(downloads_dir, "RBP_Information_human.xlsx"))

# ---- Split splicing tables by event type ----
split_splicing_by_event_type <- function(splicing_tables) {
  lapply(splicing_tables, function(df) {
    if (!"EVENT" %in% names(df)) stop("Missing 'EVENT' column")
    df$Event_Type <- sub("^Hsa([A-Z]+).*", "\\1", df$EVENT)
    split(df, df$Event_Type)
  })
}

make_upset_from_filtered <- function(EWSR1, FUS, TAF15, RBP_list, title) {
  splicing_list <- list(EWSR1 = EWSR1, FUS = FUS, TAF15 = TAF15)
  
  # Compute intersection for each condition
  rbp_sets <- lapply(splicing_list, function(df) intersect(RBP_list, df$GENE))
  all_genes <- unique(unlist(rbp_sets))
  
  # Build binary presence/absence data frame
  rbp_df <- data.frame(
    gene  = all_genes,
    EWSR1 = all_genes %in% rbp_sets$EWSR1,
    FUS   = all_genes %in% rbp_sets$FUS,
    TAF15 = all_genes %in% rbp_sets$TAF15
  )
  
  # Infer main color based on title
  main_color <- dplyr::case_when(
    grepl("pos", title, ignore.case = TRUE) ~ "#009E73",  # Green
    grepl("neg", title, ignore.case = TRUE) ~ "#D55E00",  # Red
    TRUE ~ "#0072B2"                                       # Blue (default)
  )
  
  # Define per-set color
  query_colors <- c(EWSR1 = "#EB5757", FUS = "#27AE60", TAF15 = "#4A90E2")
  
  plot <- ComplexUpset::upset(
    rbp_df,
    intersect = c("EWSR1", "FUS", "TAF15"),
    name = title,
    width_ratio = 0.15,
    sort_intersections_by = "degree",
    queries = lapply(names(query_colors), function(set) {
      ComplexUpset::upset_query(set = set, fill = query_colors[[set]])
    }),
    base_annotations = list(
      "Intersection size" = ComplexUpset::intersection_size(
        text = ggplot2::element_text(size = 8, vjust = -0.7),
        fill = main_color
      ) +
        ggplot2::theme(
          axis.text.y  = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(),
          axis.line.y  = ggplot2::element_blank()
        )
    ),
    themes = ComplexUpset::upset_default_themes(
      text = ggplot2::element_text(size = 20, face = "bold")
    ),
    set_sizes=(
      upset_set_size()
      + theme(axis.text.x=element_text(size = 16, face = "bold"))
      + ylab('Gene Number')
    )
  )
  
  
  list(plot = plot, data = rbp_df)
}


# ---- Prepare splicing input for all/pos/neg ----
splicing_filtered <- list(
  all = list(EWSR1 = EWSR1_filtered$all, FUS = FUS_filtered$all, TAF15 = TAF15_filtered$all),
  pos = list(EWSR1 = EWSR1_filtered$pos, FUS = FUS_filtered$pos, TAF15 = TAF15_filtered$pos),
  neg = list(EWSR1 = EWSR1_filtered$neg, FUS = FUS_filtered$neg, TAF15 = TAF15_filtered$neg)
)

# ---- Split each condition by event type ----
splicing_by_event <- lapply(splicing_filtered, split_splicing_by_event_type)

# ---- Generate all 12 plots ----
RBP_list <- RBP_Information_human$RBP_Name
filter_types <- names(splicing_by_event)
event_types <- c("EX", "INT", "ALTA", "ALTD")
# Store all RBP presence/absence matrices
rbp_presence_list <- list()

for (filter_type in filter_types) {
  for (etype in event_types) {
    EWSR1_df <- splicing_by_event[[filter_type]]$EWSR1[[etype]]
    FUS_df   <- splicing_by_event[[filter_type]]$FUS[[etype]]
    TAF15_df <- splicing_by_event[[filter_type]]$TAF15[[etype]]
    
    if (any(sapply(list(EWSR1_df, FUS_df, TAF15_df), is.null))) next
    
    label <- paste(toupper(filter_type), etype, sep = "_")
    plot_title <- paste("RBP –", label)
    cat("Generating plot:", plot_title, "\n")
    
    res <- make_upset_from_filtered(EWSR1_df, FUS_df, TAF15_df, RBP_list, title = plot_title)
    print(res$plot)
    
    # Save to PNG
    file_name <- paste0("RBP_upset_", label, ".png")
    ggplot2::ggsave(
      filename = file_name,
      plot = res$plot,
      width = 4000,
      height = 2800,
      units = "px",
      dpi = 300
    )
    
    rbp_presence_list[[label]] <- res$data
  }
}


# View shared RBPs in one comparison
head(rbp_presence_list$ALL_EX)

# Save all to Excel
write_xlsx(rbp_presence_list, path = "RBP_shared_presence.xlsx")





library(dplyr)
library(tidyr)
library(ggplot2)

# Loop over splicing event types
for (etype in c("EX", "INT", "ALTA", "ALTD")) {
  
  label <- paste("ALL", etype, sep = "_")
  cat("Processing", label, "\n")
  
  # Extract RBP logical matrix for this event type
  binary_matrix <- rbp_presence_list[[label]]
  if (is.null(binary_matrix)) next
  
  # Convert to long format
  binary_long <- binary_matrix %>%
    pivot_longer(cols = -gene, names_to = "Condition", values_to = "Present")
  
  # Extract the corresponding splicing data
  EWSR1_df <- splicing_by_event$all$EWSR1[[etype]]
  FUS_df   <- splicing_by_event$all$FUS[[etype]]
  TAF15_df <- splicing_by_event$all$TAF15[[etype]]
  
  # Check presence
  if (any(sapply(list(EWSR1_df, FUS_df, TAF15_df), is.null))) next
  
  # Gather E.dPsi. values into long format
  psi_long <- bind_rows(
    EWSR1_df %>% dplyr::select(gene = GENE, dPsi = E.dPsi.) %>% mutate(Condition = "EWSR1"),
    FUS_df   %>% dplyr::select(gene = GENE, dPsi = E.dPsi.) %>% mutate(Condition = "FUS"),
    TAF15_df %>% dplyr::select(gene = GENE, dPsi = E.dPsi.) %>% mutate(Condition = "TAF15")
  )
  
  # # Join presence + E.dPsi.
  # heatmap_df <- left_join(binary_long, psi_long, by = c("gene", "Condition")) %>%
  #   mutate(Present = ifelse(is.na(dPsi), FALSE, Present)) %>%
  #   group_by(gene) %>%
  #   filter(any(Present)) %>%
  #   ungroup()
  # Join presence + dPsi, then average duplicates
  heatmap_df <- left_join(binary_long, psi_long, by = c("gene", "Condition")) %>%
    group_by(gene, Condition) %>%
    summarise(
      dPsi = mean(dPsi, na.rm = TRUE),  # average across multiple splicing events
      Present = any(Present & !is.na(dPsi)),  # present only if at least one valid dPsi
      .groups = "drop"
    ) %>%
    filter(Present)
  
  # Set desired FET column order
  heatmap_df$Condition <- factor(heatmap_df$Condition, levels = c("FUS", "EWSR1", "TAF15"))
  # Remove self-gene events (i.e., FUS in FUSKO, etc.)
  heatmap_df <- heatmap_df %>%
    filter(!((gene == "FUS"   & Condition == "FUS") |
               (gene == "TAF15" & Condition == "TAF15") |
               (gene == "EWSR1" & Condition == "EWSR1")))
  
  
  # Plot heatmap
  p <- ggplot(heatmap_df %>% filter(Present),
              aes(x = Condition, y = factor(gene, levels = rev(unique(gene))), fill = dPsi)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", dPsi)), size = 3.5) +
    scale_fill_gradient2(
      low = "#1E88E5", mid = "white", high = "#FFC107",
      midpoint = 0, na.value = "grey90"
    ) +
    labs(
      title = paste("Alternative Splicing ΔPSI –", etype),
      x = "Condition", y = "Gene", fill = "E.dPsi."
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
  
  print(p)
  
  # Optional: save plot
  ggsave(filename = paste0("Splicing_dPsi_heatmap_", etype, ".png"),
         plot = p, width = 10, height = 8, units = "in", dpi = 300)
}

