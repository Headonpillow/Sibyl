#' Calculate Average Pairwise Distance for Sample Clusters Across Rarefaction Thresholds
#'
#' This function performs repeated rarefaction on a `phyloseq` object, computes 
#' ordination, and calculates the average pairwise distance for each sample 
#' across different rarefaction thresholds.
#'
#' @param input A `phyloseq` object.
#' @param repeats An integer or a vector of integers. 
#' The number of times to repeat rarefaction. A value of 1 means no repeats. 
#' If using a vector, different rarefaction thresholds will be tested sequentially. 
#' Default = 50.
#' @param t_min An integer. The minimum value for the threshold testing range. 
#' Default = 50.
#' @param t_max An integer. The maximum value for the threshold testing range. 
#' Default = 250.
#' @param t_step A numeric value. The step size for the threshold testing range. 
#' A value between 0 and 1 will cause the same threshold to be tested multiple times. 
#' Default = 5.
#' @param cores An integer. The number of cores for parallel processing. Default = 4.
#'
#' @details
#' In the case the samples do not have enough reads to reach `t_max` value, their
#' APD are set to 0 after that threshold and a warning message is printed. 
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
avg_pairwise_dist_per_sample <- function(input, repeats = 50, t_min = 50, t_max = 250, t_step = 5, cores = 4) {
  # Check if input is a Phyloseq object
  if (!inherits(input, "phyloseq")) {
    stop("Input must be a Phyloseq object, including a count table and sample data.")
  }
  
  # Extract all sample IDs from the original phyloseq object
  sample_ids <- rownames(phyloseq::sample_data(input))
  
  # Generate threshold sequence
  thresholds <- seq(t_min, t_max, by = t_step)
  
  # Initialize a data frame to store results
  avg_distances <- data.frame(Sample_ID = character(0), Threshold = numeric(0), Avg_Distance = numeric(0))
  
  # Loop through thresholds
  for (threshold in thresholds) {
    message(paste("Processing threshold:", threshold))
    
    # Perform repeated rarefaction (note: rep_raref will remove samples that do not meet the threshold)
    rarefied_data <- rep_raref(data.frame(t(phyloseq::otu_table(input))), threshold, repeats)
    
    # Perform ordination and alignment across rarefaction replicates
    ord_result <- ord_and_mean(rarefied_data$rarefied_matrix_list, repeats, cores)
    
    # Determine which samples are present in the ordination results.
    # (Since the count table is filtered in rep_raref, samples with too low counts will be missing.)
    present_samples <- rownames(ord_result$aligned_ordinations[[1]])
    
    # Calculate average pairwise distance for each sample.
    for (sample_id in sample_ids) {
      if (!(sample_id %in% present_samples)) {
        # Sample did not meet threshold: assign distance = 0
        avg_distance <- 0
      } else {
        # Extract the coordinates for this sample from each replicate.
        # Here we loop over the aligned ordinations and grab the row corresponding to sample_id.
        cluster_points <- do.call(rbind, lapply(ord_result$aligned_ordinations, function(mat) {
          mat[sample_id, , drop = FALSE]
        }))
        
        # Compute pairwise distances among the replicate points
        distances <- as.matrix(dist(cluster_points))
        avg_distance <- mean(distances[upper.tri(distances)])
      }
      
      # Store the result
      avg_distances <- rbind(avg_distances,
                             data.frame(Sample_ID = sample_id,
                                        Threshold = threshold,
                                        Avg_Distance = avg_distance,
                                        stringsAsFactors = FALSE))
    }
  }
  
  # Generate a faceted plot (each sample in its own facet)
  plot <- ggplot(avg_distances, aes(x = Threshold, y = Avg_Distance, color = Sample_ID)) +
    geom_line() +
    geom_point() +
    facet_wrap(~ Sample_ID, scales = "free_x") +
    labs(
      title = "Average Pairwise Distance Across Rarefaction Thresholds",
      x = "Rarefaction Threshold",
      y = "Average Pairwise Distance"
    ) +
    theme_minimal()
  
  result <- list(plot = plot, distances = avg_distances)
  return(result)
}
