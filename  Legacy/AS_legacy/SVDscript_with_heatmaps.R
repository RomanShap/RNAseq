# Clear the Environment ----
rm(list = ls())            # Remove all objects
graphics.off()             # Close all open plots
cat("\014")                # Clear console (works in RStudio)
# import necessary libraries ----
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
# import datasets ----
# This data comes from kallisto filtering file in the same folder. With it we filter vast tools table 
# to leave only reliably expressed genes

split_splicing_data_EX <- t(read.delim("split_splicing_data_EX.tsv", check.names = FALSE))
split_splicing_data_INT <- t(read.delim("split_splicing_data_INT.tsv", check.names = FALSE))
split_splicing_data_ALTA <- t(read.delim("split_splicing_data_ALTA.tsv", check.names = FALSE))
split_splicing_data_ALTD <- t(read.delim("split_splicing_data_ALTD.tsv", check.names = FALSE))

# SVD analysis functions ----
noize_elimination_from_matrix <- function(matrix_with_noize){
  # This function removes noise from a given matrix using Singular Value Decomposition (SVD).
  # It first ensures that the matrix is properly formatted as numeric data.
  # The function then performs SVD on the transposed matrix to decompose it into three matrices: 
  # U (samples), D (singular values), and V (genes/features).
  # The first principal component, which often represents noise or unwanted variation, is reconstructed 
  # and subtracted from the original matrix to obtain a cleaned version.
  # This method is useful for reducing batch effects or unwanted technical variation in gene expression 
  # or motif-based data before further analysis.
  
  #TESTS
  # matrix_with_noize <- matrix_of_interest
  
  # Perform SVD      
  matrix_with_noize.svd <- svd(matrix_with_noize)
  
  # Extract U, D, and V matrices from the SVD results
  u_matrix <- matrix_with_noize.svd$u  # U matrix (samples)
  d_values <- matrix_with_noize.svd$d  # Singular values (variance contribution)
  v_matrix <- matrix_with_noize.svd$v  # V matrix (genes/features)
  # Select components that are considered noise and create a noise matrix
  noize <- (u_matrix[,1,drop=FALSE] %*% t(v_matrix[,1,drop=FALSE])) * d_values[1]
  # Subtract the noise from the original matrix
  matrix_without_noize <- matrix_with_noize - noize
  
  return(matrix_without_noize)
}
SVD_variance <- function(matrix_without_noize){
  # This function performs Singular Value Decomposition (SVD) on the input matrix and computes the projections 
  # of both features and samples onto principal components. It also extracts variance explained.
  #
  # The function a matrix to analyzed with SVD:
  # 1. Computes SVD on the transposed matrix, extracting the left singular vectors (U), singular values (D), 
  #    and right singular vectors (V).
  # 2. Calculates variance explained for each principal component and stores it in a data frame, which can be 
  #    used to visualize variance distribution.
  
  # Compute SVD
  matrix_without_noize.svd <- svd(matrix_without_noize)
  
  # Extract U, D, and V matrices from the SVD results
  u_matrix <- matrix_without_noize.svd$u  # U matrix (samples)
  d_values <- matrix_without_noize.svd$d  # Singular values (variance contribution)
  v_matrix <- matrix_without_noize.svd$v  # V matrix (genes/features)
  
  # Calculate the variance explained by each singular value
  explained_variance <- (d_values^2) / sum(d_values^2)
  cumulative_variance <- cumsum(explained_variance)
  
  # Create a data frame to plot the variance
  variance_df <- data.frame(
    Component = 1:length(explained_variance),
    ExplainedVariance = explained_variance,
    CumulativeVariance = cumulative_variance
  )
  
  return(list(
    u_matrix = u_matrix,
    d_values = d_values,
    v_matrix = v_matrix,
    variance_df = variance_df
  ))
}
plot_scree_bar <- function(variance_df) {
  # This function takes a data frame containing the variance explained by each principal component 
  # and generates a bar plot ("scree plot") showing the proportion of 
  # variance explained by each component. The input data frame must have at least two columns: 
  # 'Component' (the principal component index or name) and 'ExplainedVariance' (the proportion 
  # or percentage of variance explained by each component). Optionally, cumulative variance can 
  # be added by uncommenting the corresponding lines.
  ggplot(variance_df, aes(x = Component)) +
    geom_bar(aes(y = ExplainedVariance), stat = "identity", fill = "skyblue", alpha = 0.7) +
    # geom_line(aes(y = CumulativeVariance), color = "red", linewidth = 1) +
    # geom_point(aes(y = CumulativeVariance), color = "red", size = 2) +  # Add points to cumulative line
    labs(title = "Variance Explained by Each SVD Component", 
         x = "Principal Component", 
         y = "Variance Explained") +
    theme_minimal()
}
compute_features_projections_and_correlations <- function(matrix_of_interest, 
                                                          u_matrix, 
                                                          d_values, 
                                                          v_matrix) {
  # This function computes features projections, and correlations
  # for all components from Singular Value Decomposition (SVD). It processes
  # the left singular vectors (U matrix) and singular values (D matrix) to
  # analyze how samples align with principal components.
  
  num_features <- nrow(matrix_of_interest)
  num_components <- length(d_values)
  
  # Calculate feature projections (loadings)
  # features_projections is a num_features x num_components matrix.
  # For feature i and component j, features_projections[i, j] = u_matrix[i, j] * d_values[j].
  features_projections <- u_matrix %*% diag(d_values)
  # features_projections <- u_matrix
  colnames(features_projections) <- colnames(matrix_of_interest)
  rownames(features_projections) <- rownames(matrix_of_interest)
  
  
  # Initialize a matrix to store Pearson correlations for each feature with each component
  features_correlations <- matrix(NA, nrow = num_features, ncol = num_components)
  rownames(features_correlations) <- rownames(matrix_of_interest)
  colnames(features_correlations) <- paste0("Comp", 1:num_components)
  # Calculate features correltions according to https://academic.oup.com/nar/article/42/7/4180/2435916?login=false
  features_correlations <- features_projections / sqrt(rowSums(features_projections^2))
  
  # Create a df to be used for plotting of features projections
  scatter_data_projections <- matrix_to_component_df(features_projections)
  # Rename the column to clarify what it represents
  names(scatter_data_projections)[names(scatter_data_projections) == "value"] <- "Projections"
  # Create a df to be used for plotting of features correlations
  scatter_data_correlations <- matrix_to_component_df(features_correlations)
  # Rename the column to clarify what it represents
  names(scatter_data_correlations)[names(scatter_data_correlations) == "value"] <- "Correlations"
  scatter_data = scatter_data_projections
  scatter_data$Correlations <- scatter_data_correlations$Correlations
  # Add features names to the df
  rownames(matrix_of_interest) <- make.unique(rownames(matrix_of_interest))
  features_names <- rep(rownames(matrix_of_interest), times = ncol(matrix_of_interest))
  scatter_data$Features <- features_names 
  
  return(list(features_projections = features_projections,
              features_correlations = features_correlations,
              scatter_data = scatter_data))
}
matrix_to_component_df <- function(mat, 
                                   name_for_column) {
  # Function that takes a matrix and converts it into a single column recording to which component each of the values correspond
  # Check if the input is a matrix
  if (!is.matrix(mat)) {
    stop("Input must be a matrix.")
  }
  
  # Get the number of rows and columns
  nrows <- nrow(mat)
  ncols <- ncol(mat)
  
  # Convert the matrix to a vector (this stacks the columns by default)
  values <- as.vector(mat)
  
  # Create a vector for the component labels.
  # Each column's values are labeled "Component 1", "Component 2", etc.
  components <- rep(paste("Component", seq_len(ncols)), each = nrows)
  
  # Create a data frame with the component column as the first column
  df <- data.frame(Component = components, value = values)
  
  return(df)
}
classify_features_across_components <- function(scatter_data,
                                                matrix_of_interest,
                                                percent_threshold
                                                ) {
  # Takes scatter_data df with projections and correlations of a features and filters them to select only features of interest for each component

  # Set default classification to "Neutral"
  scatter_data$Selection <- "Neutral"

  # Initialize lists to store positive and negative features per component
  positive_tables <- list()
  negative_tables <- list()

  # Loop through all components dynamically
  for (i in 1:ncol(matrix_of_interest)) {

    component_name <- paste0("Component ", i)


    # Subset the data for the current component
    component_data <- scatter_data[scatter_data$Component == component_name, ]

    # Compute quantile thresholds for the current component
    # For positive labeling: the values must be in the top (1 - percent_threshold)
    positive_corr_threshold <- quantile(component_data$Correlation, probs = 1 - percent_threshold, na.rm = TRUE)
    positive_proj_threshold <- quantile(component_data$Projection, probs = 1 - percent_threshold, na.rm = TRUE)

    # For negative labeling: the values must be in the bottom percent_threshold
    negative_corr_threshold <- quantile(component_data$Correlation, probs = percent_threshold, na.rm = TRUE)
    negative_proj_threshold <- quantile(component_data$Projection, probs = percent_threshold, na.rm = TRUE)

    # Label the features in this component as "Positive" if both values are high
    component_data$Selection[
      component_data$Correlation >= positive_corr_threshold &
        component_data$Projection >= positive_proj_threshold
    ] <- "Positive"

    # Label as "Negative" if both values are low
    component_data$Selection[
      component_data$Correlation <= negative_corr_threshold &
        component_data$Projection <= negative_proj_threshold
    ] <- "Negative"

    # Now, write the updated selection back to scatter_data
    scatter_data[scatter_data$Component == component_name, "Selection"] <- component_data$Selection

    # Extract the "Features" column for positive and negative labels
    positive_features <- component_data[component_data$Selection == "Positive", "Features"]
    negative_features <- component_data[component_data$Selection == "Negative", "Features"]

    # Convert to data frames
    positive_table <- data.frame(Feature = positive_features)
    negative_table <- data.frame(Feature = negative_features)

    # Store results for each component
    positive_tables[[component_name]] <- positive_table
    negative_tables[[component_name]] <- negative_table
  }
  
  # Plots to evaluate the data distribuition
  component1_data <- scatter_data[scatter_data$Component == "Component 1", ]
  ggplot(component1_data, aes(x = Correlations)) +
    geom_density(fill = "skyblue", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
    theme_minimal() +
    labs(title = "Density of Correlations", x = "Correlation")
  ggplot(component1_data, aes(x = Projections)) +
    geom_density(fill = "skyblue", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
    theme_minimal() +
    labs(title = "Density of Projections", x = "Projections")
  ggplot(component1_data, aes(x = Projections, y = Correlations, color = Selection)) +
    geom_point(alpha = 0.6) +
    geom_vline(
      xintercept = quantile(component1_data$Projections, c(0.05, 0.95), na.rm = TRUE),
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = quantile(component1_data$Correlations, c(0.05, 0.95), na.rm = TRUE),
      linetype = "dashed"
    ) +
    theme_minimal()
  # Return results as a list containing classified motifs and gene lists
  return(list(
    positive_tables = positive_tables,
    negative_tables = negative_tables,
    scatter_data = scatter_data
  ))
}
# classify_features_across_components <- function(scatter_data, 
#                                                 matrix_of_interest, 
#                                                 percent_threshold) {
#   # Set default classification to "Neutral"
#   scatter_data$Selection <- "Neutral"
#   
#   # Initialize lists to store positive and negative features per component
#   positive_tables <- list()
#   negative_tables <- list()
#   
#   # Loop through all components dynamically
#   for (i in 1:ncol(matrix_of_interest)) {
#     
#     component_name <- paste0("Component ", i)
#     component_data <- scatter_data[scatter_data$Component == component_name, ]
#     
#     n_features <- nrow(component_data)
#     n_top <- ceiling(n_features * percent_threshold)
#     
#     # Top/bottom Projection
#     top_proj_features <- component_data$Features[order(component_data$Projection, decreasing = TRUE)[1:n_top]]
#     bottom_proj_features <- component_data$Features[order(component_data$Projection)[1:n_top]]
#     
#     # Top/bottom Correlation
#     top_corr_features <- component_data$Features[order(component_data$Correlation, decreasing = TRUE)[1:n_top]]
#     bottom_corr_features <- component_data$Features[order(component_data$Correlation)[1:n_top]]
#     
#     # Intersections
#     positive_features <- intersect(top_proj_features, top_corr_features)
#     negative_features <- intersect(bottom_proj_features, bottom_corr_features)
#     
#     # Label selections
#     component_data$Selection <- "Neutral"
#     component_data$Selection[component_data$Features %in% positive_features] <- "Positive"
#     component_data$Selection[component_data$Features %in% negative_features] <- "Negative"
#     
#     # Update main scatter_data
#     scatter_data[scatter_data$Component == component_name, "Selection"] <- component_data$Selection
#     
#     # Store results
#     positive_tables[[component_name]] <- data.frame(Feature = positive_features)
#     negative_tables[[component_name]] <- data.frame(Feature = negative_features)
#   }
#   
#   return(list(
#     positive_tables = positive_tables,
#     negative_tables = negative_tables,
#     scatter_data = scatter_data
#   ))
# }
compute_samples_projections <- function(v_matrix, 
                                        d_values,
                                        Sample_names,
                                        Condition_lables) {
  # This function computes the projections of samples onto all available principal components 
  # from Singular Value Decomposition (SVD). It takes as input the right singular vector matrix (V) 
  # and the singular values (D) obtained from SVD. The function scales each column of V by its 
  # corresponding singular value to obtain the projected coordinates in the principal component space.
  # The number of components is automatically determined based on the dimensions of the input data. 
  # The resulting projections are returned as a data frame, where each column represents a principal 
  # component and each row corresponds to a sample. This function is useful for analyzing the 
  # contributions of samples in reduced-dimensional representations, such as PCA or SVD-based analyses.
  
  # Ensure dimensions match
  num_components <- min(length(d_values), ncol(v_matrix))
  
  # Compute samples_projections for all components
  samples_projections <- diag(d_values) %*% t(v_matrix)
  # samples_projections <- t(v_matrix)
  
  # Convert to a data frame for easier handling and transposing it to add condition and sample_names to
  # the right orientation
  samples_projections_df <- as.data.frame(t(samples_projections))
  
  # Assign meaningful column names
  colnames(samples_projections_df) <- paste0("Component_", 1:num_components)
  
  # Add sample and condition columns to samples_projections_df
  samples_projections_df$Sample <- Sample_names
  samples_projections_df$Condition <- Condition_lables
  # Create a summary df from samples_projections_df
  samples_proj_summary <- samples_projections_df %>%
    group_by(Condition) %>%
    summarise(across(starts_with("Component"), list(Mean = mean, SD = sd), .names = "{.col}_{.fn}"))
  
  return(samples_proj_summary)
}
plot_All_components_by_condition <- function(samples_proj_summary) {
  # Identify all component mean and SD columns dynamically
  component_cols <- grep("_Mean$", colnames(samples_proj_summary), value = TRUE)
  
  # Define color mapping with WT first and gray
  unique_conditions <- unique(samples_proj_summary$Condition)
  unique_conditions <- c("SHSY5Y_WT", setdiff(unique_conditions, "SHSY5Y_WT"))  # Ensure WT is first
  
  color_mapping <- c("SHSY5Y_WT" = "gray", setNames(scales::hue_pal()(length(unique_conditions) - 1), setdiff(unique_conditions, "SHSY5Y_WT")))
  
  # Create a list to store all plots
  plot_list <- list()
  
  # Loop through all components and generate plots
  for (component in component_cols) {
    
    # Find corresponding SD column
    sd_col <- gsub("_Mean$", "_SD", component)
    component_number <- gsub("_Mean$", "", component)  # Extract component number
    
    # Plot function
    plot <- ggplot(samples_proj_summary, aes(x = factor(Condition, levels = unique_conditions), y = .data[[component]], group = 1)) +
      geom_point(aes(color = Condition), size = 3) +
      geom_errorbar(aes(ymin = .data[[component]] - .data[[sd_col]], ymax = .data[[component]] + .data[[sd_col]]), width = 0.2) +
      labs(title = paste("Mean of Component", component_number, "by Condition"),
           x = "Condition",
           y = paste("Component", component_number, "Value")) +
      theme_minimal() +
      scale_color_manual(values = color_mapping)
    
    # Store the plot in the list
    plot_list[[paste0("Component_", component_number)]] <- plot
  }
  
  return(plot_list)
}
# SVD master function ----
master_function_All_samples <- function(matrix_of_interest,
                                        noize_cancelling = FALSE,
                                        percent_threshold){
  # TESTS
  # percent_threshold = 0.05
  # matrix_of_interest <- split_splicing_data_EX

  
  # Sample names
  Sample_names <- matrix_of_interest[, 1]
  # Feature names
  Feature_names <- colnames(matrix_of_interest[,-1])
  
  # Create condition vector
  Condition_lables <- substr(Sample_names, 1, nchar(Sample_names) - 2)
  # Rename rows and shape the matrix in the right way
  rownames(matrix_of_interest) <- matrix_of_interest[,1]
  matrix_of_interest <- matrix_of_interest[, -1]
  matrix_of_interest <- apply(matrix_of_interest, 2, as.numeric)
  matrix_of_interest <- t(matrix_of_interest)
  rownames(matrix_of_interest) <- Feature_names
  # Check if the matrix is solid without nans
  check_no_missing <- function(mat) {
    if (any(is.na(mat))) {
      stop("there are NaN or NA in the matrix")
    }
  }
  check_no_missing(matrix_of_interest)
  
  # Eliminate the first component if necessary
  if (noize_cancelling){
    matrix_of_interest <- noize_elimination_from_matrix(matrix_of_interest)
  }
  
  # Compute SVD and variance
  svd_variance <- SVD_variance(matrix_of_interest)
  # Extract U, D, and V matrices from the SVD results
  u_matrix <- svd_variance$u_matrix  # U matrix (samples)
  d_values <- svd_variance$d_values  # Singular values (variance contribution)
  v_matrix <- svd_variance$v_matrix  # V matrix (genes/features)
  variance_df <- svd_variance$variance_df #Variance from SVD
  
  scree_variance <- plot_scree_bar(variance_df)
  print(scree_variance)
  
  samples_proj_summary <- compute_samples_projections(v_matrix, 
                                                      d_values,
                                                      Sample_names,
                                                      Condition_lables)
  
  #Check if samples are different between each other
  samples_projections_plots <- plot_All_components_by_condition(samples_proj_summary)
  print(samples_projections_plots)
  
  # Compute projections and correlations of FESTURES on components
  features_projections_correlations <- compute_features_projections_and_correlations(matrix_of_interest, 
                                                                                    u_matrix, 
                                                                                    d_values, 
                                                                                    v_matrix)
  features_projections <-  features_projections_correlations$features_projections
  features_correlations <-  features_projections_correlations$features_correlations
  # df that is used fore plotting
  scatter_data <- features_projections_correlations$scatter_data
  
  # Select features that are playing key roles in variance
  feature_tables <- classify_features_across_components(scatter_data, 
                                                        matrix_of_interest, 
                                                        percent_threshold)
  scatter_data <- feature_tables$scatter_data
  positive_tables <- feature_tables$positive_tables
  negative_tables <- feature_tables$negative_tables
  return(list(positive_tables = positive_tables,
              negative_tables = negative_tables))
}
# results ----

tables_ismara <- master_function_All_samples(split_splicing_data_EX,T,0.05)
positive_features <- tables_ismara[["positive_tables"]][["Component 1"]][["Feature"]]
negative_features <- tables_ismara[["negative_tables"]][["Component 1"]][["Feature"]]

selected_features <- unique(c(positive_features, negative_features))
heatmap_matrix <- t(split_splicing_data_EX)
colnames(heatmap_matrix) <- heatmap_matrix[1,]
heatmap_matrix = heatmap_matrix[-1,]
filtered_expression <- heatmap_matrix[rownames(heatmap_matrix) %in% selected_features, ]
filtered_expression <- apply(filtered_expression, 2, as.numeric)
# Scale the data (optional, for better visualization)
filtered_expression_scaled <- t(scale(t(filtered_expression)))  # Scale by gene (row)
# filtered_expression_scaled <- matrix_of_interest[!apply(filtered_expression_scaled, 1, function(x) any(is.nan(x))), ]

colnames(filtered_expression_scaled) <- colnames(heatmap_matrix)
heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
new_order <- c("SHSY5Y_WT_1", "SHSY5Y_WT_2", "SHSY5Y_WT_3",
               "SHSY5Y_EWSR1KO_1", "SHSY5Y_EWSR1KO_2", "SHSY5Y_EWSR1KO_3",
               "SHSY5Y_FUSKO_1", "SHSY5Y_FUSKO_2", "SHSY5Y_FUSKO_3",
               "SHSY5Y_TAF15KO_1", "SHSY5Y_TAF15KO_2", "SHSY5Y_TAF15KO_3")
filtered_expression_scaled <- filtered_expression_scaled[, new_order]

# Create the heatmap
pheatmap(filtered_expression_scaled,
         color = heatmap_colors,
         cluster_rows = nrow(filtered_expression_scaled) > 1,  # Cluster only if >1 gene
         # cluster_cols = ncol(filtered_expression_scaled) > 1,  # Cluster only if >1 sample
         cluster_cols = F,
         scale = "none",
         fontsize_row = 8,
         fontsize_col = 10,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "EX Heatmap of 1st Component Features")

positive_features <- tables_ismara[["positive_tables"]][["Component 2"]][["Feature"]]
negative_features <- tables_ismara[["negative_tables"]][["Component 2"]][["Feature"]]

selected_features <- unique(c(positive_features, negative_features))
heatmap_matrix <- t(split_splicing_data_EX)
colnames(heatmap_matrix) <- heatmap_matrix[1,]
heatmap_matrix = heatmap_matrix[-1,]
filtered_expression <- heatmap_matrix[rownames(heatmap_matrix) %in% selected_features, ]
filtered_expression <- apply(filtered_expression, 2, as.numeric)
# Scale the data (optional, for better visualization)
filtered_expression_scaled <- t(scale(t(filtered_expression)))  # Scale by gene (row)
# filtered_expression_scaled <- matrix_of_interest[!apply(filtered_expression_scaled, 1, function(x) any(is.nan(x))), ]

colnames(filtered_expression_scaled) <- colnames(heatmap_matrix)
heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
new_order <- c("SHSY5Y_WT_1", "SHSY5Y_WT_2", "SHSY5Y_WT_3",
               "SHSY5Y_EWSR1KO_1", "SHSY5Y_EWSR1KO_2", "SHSY5Y_EWSR1KO_3",
               "SHSY5Y_FUSKO_1", "SHSY5Y_FUSKO_2", "SHSY5Y_FUSKO_3",
               "SHSY5Y_TAF15KO_1", "SHSY5Y_TAF15KO_2", "SHSY5Y_TAF15KO_3")
filtered_expression_scaled <- filtered_expression_scaled[, new_order]

# Create the heatmap
pheatmap(filtered_expression_scaled,
         color = heatmap_colors,
         cluster_rows = nrow(filtered_expression_scaled) > 1,  # Cluster only if >1 gene
         # cluster_cols = ncol(filtered_expression_scaled) > 1,  # Cluster only if >1 sample
         cluster_cols = F,
         scale = "none",
         fontsize_row = 8,
         fontsize_col = 10,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "EX Heatmap of 2st Component Features")


tables_ismara <- master_function_All_samples(split_splicing_data_INT,T,0.05)
positive_features <- tables_ismara[["positive_tables"]][["Component 1"]][["Feature"]]
negative_features <- tables_ismara[["negative_tables"]][["Component 1"]][["Feature"]]

selected_features <- unique(c(positive_features, negative_features))
heatmap_matrix <- t(split_splicing_data_INT)
colnames(heatmap_matrix) <- heatmap_matrix[1,]
heatmap_matrix = heatmap_matrix[-1,]
filtered_expression <- heatmap_matrix[rownames(heatmap_matrix) %in% selected_features, ]
filtered_expression <- apply(filtered_expression, 2, as.numeric)
# Scale the data (optional, for better visualization)
filtered_expression_scaled <- t(scale(t(filtered_expression)))  # Scale by gene (row)
filtered_expression_scaled <- filtered_expression_scaled[complete.cases(filtered_expression_scaled), ]

colnames(filtered_expression_scaled) <- colnames(heatmap_matrix)
heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
new_order <- c("SHSY5Y_WT_1", "SHSY5Y_WT_2", "SHSY5Y_WT_3",
               "SHSY5Y_EWSR1KO_1", "SHSY5Y_EWSR1KO_2", "SHSY5Y_EWSR1KO_3",
               "SHSY5Y_FUSKO_1", "SHSY5Y_FUSKO_2", "SHSY5Y_FUSKO_3",
               "SHSY5Y_TAF15KO_1", "SHSY5Y_TAF15KO_2", "SHSY5Y_TAF15KO_3")
filtered_expression_scaled <- filtered_expression_scaled[, new_order]

# Create the heatmap
pheatmap(filtered_expression_scaled,
         color = heatmap_colors,
         cluster_rows = nrow(filtered_expression_scaled) > 1,  # Cluster only if >1 gene
         # cluster_cols = ncol(filtered_expression_scaled) > 1,  # Cluster only if >1 sample
         cluster_cols = F,
         scale = "none",
         fontsize_row = 8,
         fontsize_col = 10,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "INT Heatmap of 1st Component Features")

positive_features <- tables_ismara[["positive_tables"]][["Component 2"]][["Feature"]]
negative_features <- tables_ismara[["negative_tables"]][["Component 2"]][["Feature"]]

selected_features <- unique(c(positive_features, negative_features))
heatmap_matrix <- t(split_splicing_data_INT)
colnames(heatmap_matrix) <- heatmap_matrix[1,]
heatmap_matrix = heatmap_matrix[-1,]
filtered_expression <- heatmap_matrix[rownames(heatmap_matrix) %in% selected_features, ]
filtered_expression <- apply(filtered_expression, 2, as.numeric)
# Scale the data (optional, for better visualization)
filtered_expression_scaled <- t(scale(t(filtered_expression)))  # Scale by gene (row)
# filtered_expression_scaled <- matrix_of_interest[!apply(filtered_expression_scaled, 1, function(x) any(is.nan(x))), ]

colnames(filtered_expression_scaled) <- colnames(heatmap_matrix)
heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
new_order <- c("SHSY5Y_WT_1", "SHSY5Y_WT_2", "SHSY5Y_WT_3",
               "SHSY5Y_EWSR1KO_1", "SHSY5Y_EWSR1KO_2", "SHSY5Y_EWSR1KO_3",
               "SHSY5Y_FUSKO_1", "SHSY5Y_FUSKO_2", "SHSY5Y_FUSKO_3",
               "SHSY5Y_TAF15KO_1", "SHSY5Y_TAF15KO_2", "SHSY5Y_TAF15KO_3")
filtered_expression_scaled <- filtered_expression_scaled[, new_order]

# Create the heatmap
pheatmap(filtered_expression_scaled,
         color = heatmap_colors,
         cluster_rows = nrow(filtered_expression_scaled) > 1,  # Cluster only if >1 gene
         # cluster_cols = ncol(filtered_expression_scaled) > 1,  # Cluster only if >1 sample
         cluster_cols = F,
         scale = "none",
         fontsize_row = 8,
         fontsize_col = 10,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "INT Heatmap of 2st Component Features")

tables_ismara <- master_function_All_samples(split_splicing_data_ALTA,T,0.05)
positive_features <- tables_ismara[["positive_tables"]][["Component 1"]][["Feature"]]
negative_features <- tables_ismara[["negative_tables"]][["Component 1"]][["Feature"]]

selected_features <- unique(c(positive_features, negative_features))
heatmap_matrix <- t(split_splicing_data_ALTA)
colnames(heatmap_matrix) <- heatmap_matrix[1,]
heatmap_matrix = heatmap_matrix[-1,]
filtered_expression <- heatmap_matrix[rownames(heatmap_matrix) %in% selected_features, ]
filtered_expression <- apply(filtered_expression, 2, as.numeric)
# Scale the data (optional, for better visualization)
filtered_expression_scaled <- t(scale(t(filtered_expression)))  # Scale by gene (row)
# filtered_expression_scaled <- matrix_of_interest[!apply(filtered_expression_scaled, 1, function(x) any(is.nan(x))), ]

colnames(filtered_expression_scaled) <- colnames(heatmap_matrix)
heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
new_order <- c("SHSY5Y_WT_1", "SHSY5Y_WT_2", "SHSY5Y_WT_3",
               "SHSY5Y_EWSR1KO_1", "SHSY5Y_EWSR1KO_2", "SHSY5Y_EWSR1KO_3",
               "SHSY5Y_FUSKO_1", "SHSY5Y_FUSKO_2", "SHSY5Y_FUSKO_3",
               "SHSY5Y_TAF15KO_1", "SHSY5Y_TAF15KO_2", "SHSY5Y_TAF15KO_3")
filtered_expression_scaled <- filtered_expression_scaled[, new_order]

# Create the heatmap
pheatmap(filtered_expression_scaled,
         color = heatmap_colors,
         cluster_rows = nrow(filtered_expression_scaled) > 1,  # Cluster only if >1 gene
         # cluster_cols = ncol(filtered_expression_scaled) > 1,  # Cluster only if >1 sample
         cluster_cols = F,
         scale = "none",
         fontsize_row = 8,
         fontsize_col = 10,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "ALTA Heatmap of 1st Component Features")

tables_ismara <- master_function_All_samples(split_splicing_data_ALTD,T,0.05)
positive_features <- tables_ismara[["positive_tables"]][["Component 1"]][["Feature"]]
negative_features <- tables_ismara[["negative_tables"]][["Component 1"]][["Feature"]]

selected_features <- unique(c(positive_features, negative_features))
heatmap_matrix <- t(split_splicing_data_ALTD)
colnames(heatmap_matrix) <- heatmap_matrix[1,]
heatmap_matrix = heatmap_matrix[-1,]
filtered_expression <- heatmap_matrix[rownames(heatmap_matrix) %in% selected_features, ]
filtered_expression <- apply(filtered_expression, 2, as.numeric)
# Scale the data (optional, for better visualization)
filtered_expression_scaled <- t(scale(t(filtered_expression)))  # Scale by gene (row)
# filtered_expression_scaled <- matrix_of_interest[!apply(filtered_expression_scaled, 1, function(x) any(is.nan(x))), ]

colnames(filtered_expression_scaled) <- colnames(heatmap_matrix)
heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
new_order <- c("SHSY5Y_WT_1", "SHSY5Y_WT_2", "SHSY5Y_WT_3",
               "SHSY5Y_EWSR1KO_1", "SHSY5Y_EWSR1KO_2", "SHSY5Y_EWSR1KO_3",
               "SHSY5Y_FUSKO_1", "SHSY5Y_FUSKO_2", "SHSY5Y_FUSKO_3",
               "SHSY5Y_TAF15KO_1", "SHSY5Y_TAF15KO_2", "SHSY5Y_TAF15KO_3")
filtered_expression_scaled <- filtered_expression_scaled[, new_order]

# Create the heatmap
pheatmap(filtered_expression_scaled,
         color = heatmap_colors,
         cluster_rows = nrow(filtered_expression_scaled) > 1,  # Cluster only if >1 gene
         # cluster_cols = ncol(filtered_expression_scaled) > 1,  # Cluster only if >1 sample
         cluster_cols = F,
         scale = "none",
         fontsize_row = 8,
         fontsize_col = 10,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "ALTD Heatmap of 1st Component Features")

positive_features <- tables_ismara[["positive_tables"]][["Component 2"]][["Feature"]]
negative_features <- tables_ismara[["negative_tables"]][["Component 2"]][["Feature"]]

selected_features <- unique(c(positive_features, negative_features))
heatmap_matrix <- t(split_splicing_data_ALTD)
colnames(heatmap_matrix) <- heatmap_matrix[1,]
heatmap_matrix = heatmap_matrix[-1,]
filtered_expression <- heatmap_matrix[rownames(heatmap_matrix) %in% selected_features, ]
filtered_expression <- apply(filtered_expression, 2, as.numeric)
# Scale the data (optional, for better visualization)
filtered_expression_scaled <- t(scale(t(filtered_expression)))  # Scale by gene (row)
# filtered_expression_scaled <- matrix_of_interest[!apply(filtered_expression_scaled, 1, function(x) any(is.nan(x))), ]

colnames(filtered_expression_scaled) <- colnames(heatmap_matrix)
heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
new_order <- c("SHSY5Y_WT_1", "SHSY5Y_WT_2", "SHSY5Y_WT_3",
               "SHSY5Y_EWSR1KO_1", "SHSY5Y_EWSR1KO_2", "SHSY5Y_EWSR1KO_3",
               "SHSY5Y_FUSKO_1", "SHSY5Y_FUSKO_2", "SHSY5Y_FUSKO_3",
               "SHSY5Y_TAF15KO_1", "SHSY5Y_TAF15KO_2", "SHSY5Y_TAF15KO_3")
filtered_expression_scaled <- filtered_expression_scaled[, new_order]

# Create the heatmap
pheatmap(filtered_expression_scaled,
         color = heatmap_colors,
         cluster_rows = nrow(filtered_expression_scaled) > 1,  # Cluster only if >1 gene
         # cluster_cols = ncol(filtered_expression_scaled) > 1,  # Cluster only if >1 sample
         cluster_cols = F,
         scale = "none",
         fontsize_row = 8,
         fontsize_col = 10,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "ALTD Heatmap of 2st Component Features")