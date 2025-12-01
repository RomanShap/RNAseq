# 03_split_by_event_type.R
#
# Step 3: Split filtered splicing events by event type (EX, INT, ALTA, ALTD)
# and direction (all/pos/neg).
#
# This script uses the list of filtered splicing tables from step 01 to
# classify events by their type prefix and ΔPSI sign.  The output is a
# nested list structure `EVENT_SPLITS` where the top level keys are
# conditions, the second level keys are event types, and the third level keys
# are directions (all/pos/neg).  It also summarises the number of events
# belonging to each category and optionally plots the distribution of event
# types for each condition.

message("Step 03: Splitting events by type and direction")

## Create output directory -----------------------------------------------------
step_output <- file.path(PROJECT_ROOT, "output", "03_event_types")
if (!dir.exists(step_output)) dir.create(step_output, recursive = TRUE)

if (!exists("FILTERED_EVENTS", envir = .GlobalEnv)) {
  stop("FILTERED_EVENTS object is missing; please run step 01 first")
}
filtered_tables <- get("FILTERED_EVENTS", envir = .GlobalEnv)

## Split events ---------------------------------------------------------------
event_splits <- list()
summary_list <- list()
for (cond in names(filtered_tables)) {
  df <- filtered_tables[[cond]]
  type_list <- split_event_types(df)
  cond_split <- list()
  cond_summary <- data.frame()
  for (etype in names(type_list)) {
    sub_df <- type_list[[etype]]
    dir_list <- split_by_direction(sub_df)
    cond_split[[etype]] <- dir_list
    # Summarise counts
    counts <- sapply(dir_list, nrow)
    cond_summary <- rbind(cond_summary, data.frame(
      condition = cond,
      event_type = etype,
      all = counts["all"],
      pos = counts["pos"],
      neg = counts["neg"]
    ))
  }
  event_splits[[cond]] <- cond_split
  summary_list[[cond]] <- cond_summary
}

summary_df <- dplyr::bind_rows(summary_list)

## Save split tables ----------------------------------------------------------
splits_rds <- file.path(step_output, "event_splits.rds")
saveRDS(event_splits, splits_rds)
summary_tsv <- file.path(step_output, "event_type_summary.tsv")
readr::write_tsv(summary_df, summary_tsv)
message(sprintf("Saved event splits to %s and summary to %s", splits_rds, summary_tsv))

## Plot event type distribution per condition ---------------------------------
plot_path <- file.path(step_output, "event_type_distribution.png")
library(ggplot2)
dist_data <- summary_df %>%
  tidyr::pivot_longer(cols = c(all, pos, neg), names_to = "direction", values_to = "count")

p <- ggplot(dist_data, aes(x = event_type, y = count, fill = direction)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~ condition, scales = "free_y") +
  labs(title = "Distribution of splicing event types by direction",
       x = "Event type", y = "Number of events") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(plot_path, p, width = 8, height = 4 + 0.5 * length(unique(summary_df$condition)))
message(sprintf("Event type distribution plot saved to %s", plot_path))

## Return splits to global environment ----------------------------------------
assign("EVENT_SPLITS", event_splits, envir = .GlobalEnv)
message("Step 03 completed")