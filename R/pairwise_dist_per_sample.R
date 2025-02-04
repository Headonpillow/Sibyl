#' Calculate Average Pairwise Distance for Clusters Across Rarefaction Thresholds
#'
#' This function performs repeated rarefaction on a phyloseq object, computes 
#' ordination, and calculates the average pairwise distance for each sample 
#' across different rarefaction thresholds.
#'
#' @param input A `phyloseq` object.
#' @param repeats An integer. The number of rarefaction repeats. Default = 10.
#' @param t_min An integer. The minimum rarefaction threshold. Default = 50.
#' @param t_max An integer. The maximum rarefaction threshold. Default = 250.
#' @param t_step A numeric value. The step size for thresholds. Default = 10.
#' @param cores An integer. The number of cores for parallel processing. Default = 4.
#'
#' @return A list containing:
#'   - `plot`: A ggplot2 object showing average pairwise distances across thresholds.
#'   - `distances`: A data frame with calculated distances per sample.
#'
#' @importFrom phyloseq sample_data otu_table
#' @importFrom vegan vegdist
#' @importFrom dplyr mutate filter group_by summarise
#' @importFrom ggplot2 ggplot aes geom_line geom_point facet_wrap labs theme_minimal
#' @importFrom stats dist
#' @export
#' 
avg_pairwise_dist_per_sample <- function(input, repeats = 10, t_min = 50, t_max = 250, t_step = 10, cores = 4) {
  # Check if input is a Phyloseq object
  if (!inherits(input, "phyloseq")) {
    stop("Input must be a Phyloseq object, including a count table and sample data.")
  }
  
  # Extract sample IDs
  sample_ids <- rownames(phyloseq::sample_data(input)) 
  
  # Generate threshold sequence
  thresholds <- seq(t_min, t_max, by = t_step)
  
  # Initialize a data frame to store results
  avg_distances <- data.frame(Sample_ID = character(0), Threshold = numeric(0), Avg_Distance = numeric(0))
  
  # Loop through thresholds
  for (threshold in thresholds) {
    message(paste("Processing threshold:", threshold))
    
    # Perform repeated rarefaction
    rarefied_data <- rep_raref(data.frame(t(phyloseq::otu_table(input))), threshold, repeats)
    
    # Perform ordination
    ord_result <- ord_and_mean(rarefied_data$rarefied_matrix_list, repeats, cores)
    
    # Combine all aligned ordinations into one data frame
    all_points <- do.call(rbind, ord_result$aligned_ordinations)
    all_points_df <- as.data.frame(all_points)
    colnames(all_points_df) <- c("Dim1", "Dim2")
    all_points_df$sample_id <- rep(sample_ids, repeats)
    
    # Calculate average pairwise distance for each sample_id
    for (sample_id in sample_ids) {
      # Subset the points for this sample_id
      cluster_points <- all_points_df[all_points_df$sample_id == sample_id, 1:2]
      
      # Calculate pairwise distances
      distances <- as.matrix(dist(cluster_points))
      
      # Compute average pairwise distance
      avg_distance <- mean(distances[upper.tri(distances)])
      
      # Store the result
      avg_distances <- rbind(avg_distances, data.frame(Sample_ID = sample_id, Threshold = threshold, Avg_Distance = avg_distance))
    }
  }
  
  # Generate a single faceted plot
  plot <- ggplot(avg_distances, aes(x = Threshold, y = Avg_Distance, color = Sample_ID)) +
    geom_line() +
    geom_point() +
    facet_wrap(
      ~ Sample_ID, 
      scales = "free_x"
    ) +
    labs(
      title = "Average Pairwise Distance Across Rarefaction Thresholds",
      x = "Rarefaction Threshold",
      y = "Average Pairwise Distance"
    ) +
    theme_minimal()
  
  result <- list(plot = plot, distances = avg_distances)
  return(result)
}


