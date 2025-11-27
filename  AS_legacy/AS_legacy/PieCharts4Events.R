# ---- Dependencies ----
suppressWarnings(suppressMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
}))

# ---- Function: pie chart of event-type counts (EX, INT, ALTA, ALTD) ----
plot_event_type_pie <- function(df,
                                event_col = "EVENT",
                                title = "Event type counts (EX / INT / ALTA / ALTD)",
                                colors = NULL,
                                save_path = NULL,
                                width = 5, height = 5, dpi = 300) {
  # Basic checks
  if (!is.data.frame(df)) stop("Input must be a data.frame or tibble.")
  if (!event_col %in% names(df)) stop(sprintf("Column '%s' not found.", event_col))
  
  # Extract the raw event strings
  events_raw <- as.character(df[[event_col]])
  
  # Try to extract the type from standard VAST-like IDs: e.g., "HsaEX...", "HsaINT...", "HsaALTA...", "HsaALTD..."
  # Fallback: detect presence of EX/INT/ALTA/ALTD anywhere in the string (case-insensitive).
  type_from_prefix <- str_match(events_raw, "^Hsa([A-Z]+)")[, 2]
  type_detected <- ifelse(!is.na(type_from_prefix), type_from_prefix, NA_character_)
  type_detected[is.na(type_detected) & str_detect(events_raw, regex("(^|[^A-Z])EX([^A-Z]|$)", ignore_case = TRUE))]   <- "EX"
  type_detected[is.na(type_detected) & str_detect(events_raw, regex("(^|[^A-Z])INT([^A-Z]|$)", ignore_case = TRUE))]  <- "INT"
  type_detected[is.na(type_detected) & str_detect(events_raw, regex("(^|[^A-Z])ALTA([^A-Z]|$)", ignore_case = TRUE))]                    <- "ALTA"
  type_detected[is.na(type_detected) & str_detect(events_raw, regex("(^|[^A-Z])ALTD([^A-Z]|$)", ignore_case = TRUE))]                    <- "ALTD"
  
  # Keep only the requested event classes, enforce order
  classes <- c("EX", "INT", "ALTA", "ALTD")
  type_factor <- factor(type_detected, levels = classes)
  
  # Build counts, including zeros for missing classes
  counts_tbl <- tibble(type = type_factor) |>
    filter(!is.na(type)) |>
    count(type, name = "n") |>
    tidyr::complete(type = factor(classes, levels = classes), fill = list(n = 0)) |>
    mutate(label = paste0(as.character(type), " (", n, ")"))
  
  # If all counts are zero, abort early with a blank plot
  total_n <- sum(counts_tbl$n)
  if (total_n == 0) {
    warning("No EX/INT/ALTA/ALTD events detected.")
    p <- ggplot() + theme_void() + ggtitle(title)
    return(list(counts = counts_tbl, plot = p))
  }
  
  # Build pie chart with counts; labels are counts placed on slices
  p <- ggplot(counts_tbl, aes(x = "", y = n, fill = type)) +
    geom_col(width = 1) +
    coord_polar(theta = "y") +
    # Show count labels inside slices; hide "0" labels automatically
    geom_text(aes(label = ifelse(n > 0, n, "")),
              position = position_stack(vjust = 0.5), size = 4) +
    labs(title = title, x = NULL, y = NULL, fill = "Event type") +
    theme_void(base_size = 12) +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5))
  
  # Optional manual colors if provided as a named vector c(EX="...", INT="...", ALTA="...", ALTD="...")
  if (!is.null(colors)) {
    # Only apply colors for provided names; others use defaults
    known_cols <- intersect(names(colors), classes)
    if (length(known_cols) > 0) {
      p <- p + scale_fill_manual(values = colors[known_cols], drop = FALSE)
    }
  }
  
  # Optional save to file
  if (!is.null(save_path)) {
    ggsave(filename = save_path, plot = p, width = width, height = height, dpi = dpi)
  }
  
  return(list(counts = counts_tbl, plot = p))
}

# ---- Example usage with your object ----
res <- plot_event_type_pie(EWSR1_filtered[["all"]],
                           title = "EWSR1: Event type counts")
res$counts   # table with counts
res$plot     # pie chart

res <- plot_event_type_pie(FUS_filtered[["all"]],
                           title = "FUS: Event type counts")
res$counts   # table with counts
res$plot     # pie chart

res <- plot_event_type_pie(TAF15_filtered[["all"]],
                           title = "TAF15: Event type counts")
res$counts   # table with counts
res$plot     # pie chart
