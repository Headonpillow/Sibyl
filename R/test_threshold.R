#' Test Different Rarefaction Thresholds Using Repeated Rarefaction
#'
#' This function is a wrapper around the steps of `repeated_rarefaction` that evaluates 
#' different rarefaction thresholds. It runs repeated rarefactions on a 
#' `phyloseq` object and calculates clustering performance using the Calinski-Harabasz index.
#' Moreover, the function also calculates the average pairwise distance for each sample 
#' across different rarefaction thresholds, to give a measure of sample concordance.
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
#' @param group A string. A string. Column name in `sample_data()` used 
#' for grouping the samples. It is also on this value that the Calinski-Harabasz
#' index is calculated.
#' @param cores An integer. Number of cores for parallel processing. Default = 4.
#' 
#' @details
#' Higher index values are better, and the plateauing of the index values when 
#' testing higher thresholds indicate that the clouds of repetitions are compact
#' enough to not being affected by increasing the threshold. 
#' The Calinski-Harabasz value takes into account both within cluster and between
#' cluster distances, and in the scope of this function is calculated on the
#' `group` parameter. 
#' A sudden drop in CH value means that samples might have been removed because 
#' they did not have enough reads to reach the threshold. In this case a warning
#' is printed listing the samples which have been removed.
#' In `avg_distances`, if samples do not have enough reads to reach `t_max` value, their
#' APD are set to 0 after that threshold and a warning message is printed. 
#'
#' @return A list containing:
#'   - `index_plot`: A `ggplot` object showing Calinski-Harabasz index vs. rarefaction threshold.
#'   - `ordination_plots`: A list containing the ordination plots for tested threshold values.
#'   The ordination plots have different levels corresponding to the number of repeats.
#'   - `avg_distances`: A data frame (or list of data frames if using multiple repeat 
#'   values) with average pairwise distances of the cloud generated from each 
#'   sample repetition attempts (for each threshold).
#'
#' @importFrom phyloseq sample_data otu_table sample_data<-
#' @importFrom vegan vegdist
#' @importFrom clusterSim index.G1
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot aes geom_point labs geom_smooth xlim ylim
#' @export
#' 
test_threshold <- function(input, repeats = 50, t_min = 50, t_max = 250, t_step = 5, group = "sample_id", cores = 4) {
  # Check if input is a Phyloseq object
  if (inherits(input, "phyloseq")) {
    physeq <- input
  } else {
    stop("Input must be a Phyloseq object, including a count table and sample data.")
  }
  
  # Make the rownames of the Phyloseq object a new "sample_id" variable for the sample data.
  # (this covers the case in which no sample_id column is present in the sample data)
  # Then set it to a separate variable, and we need one.
  ### TODO: when the sample_id does not exist tho
  sample_data(physeq)$sample_id <- rownames(sample_data(physeq))
  
  # Extract all sample IDs from the original phyloseq object (if samples are later
  # removed because of unmet rarefaction threshold)
  sample_ids <- rownames(phyloseq::sample_data(input))
  
  # ============ Checks and warnings
  if (!(group %in% names(sample_data(physeq)))) {
    stop(paste("'",group,"' is not a column name in the sample information in the inputed phyloseq object.
                  repeated_rarefaction needs an existing column to group samples by.", sep=""))
  }
  
  if (!(is.double(repeats))){
    stop(paste("Input for repeats: '", repeats, "' is not an integer.", sep=""))
  }
  
  # Determine the thresholds to loop over between min and max and specified step.
  thresholds <- seq.int(t_min, t_max, by = t_step)
  # Create the plot list which will hold the ordinations and the matrix which 
  # will hold index calculations results.
  plots <- list()
  index_data <- matrix(nrow = 0, ncol = 3)
  
  # Initialize a list to store avg_distances data frames for each repeat amount.
  avg_distances_list <- list()
  
  for (y in repeats) {
    # Create lists to store ordination points coordinates for a specific repeat 
    # amount and centroid coordinates.
    ordination_points <- list()
    consensus_coordinates <- list()
    # Initialize also a data frame to store pairwise distances of samples point clouds.
    avg_distances <- data.frame(Sample_ID = character(0), Threshold = numeric(0), Avg_Distance = numeric(0))
    
    for (x in thresholds) {
      message(paste("Running with", y, "repeats and", x, "threshold"))
      
      step1 <- rep_raref(data.frame(t(otu_table(physeq))), threshold = x, repeats = y)
      step2 <- ord_and_mean(step1$rarefied_matrix_list, repeats)
      
      # ============= AVERAGE PAIRWISE DISTANCE CALCULATION (THIS THRESHOLD)
      
      # Determine which samples are present in the ordination results.
      # (Since the count table is filtered in rep_raref, samples with too low counts will be missing.)
      present_samples <- rownames(step2$aligned_ordinations[[1]])
      # Calculate average pairwise distance for each sample.
      for (sample_id in sample_ids) {
        if (!(sample_id %in% present_samples)) {
          # Sample did not meet threshold: assign distance = 0
          avg_distance <- 0
        } else {
          # Extract the coordinates for this sample from each replicate.
          # Here we loop over the aligned ordinations and grab the row corresponding to sample_id.
          cluster_points <- do.call(rbind, lapply(step2$aligned_ordinations, function(mat) {
            mat[sample_id, , drop = FALSE]
          }))
          
          # Compute pairwise distances among the replicate points
          distances <- as.matrix(dist(cluster_points))
          avg_distance <- mean(distances[upper.tri(distances)])
        }
        
        # Store the result
        avg_distances <- rbind(avg_distances,
                               data.frame(Sample_ID = sample_id,
                                          Threshold = x,
                                          Avg_Distance = avg_distance,
                                          stringsAsFactors = FALSE))
      }
      
      # ============================== STORE COORDINATES (THIS THRESHOLD)
      
      # Store aligned ordination points for this threshold
      # Note: They may not actually be "aligned" in any external sense now,
      # but are just coordinates from each repetition's ordination.
      ordination_points[[paste0("threshold_", x)]] <- step2$aligned_ordinations
      
      # Store consensus coordinates for this threshold
      consensus_coords <- step2$consensus_coordinates
      colnames(consensus_coords) <- c("Dim1", "Dim2")
      consensus_coordinates[[paste0("threshold_", x)]] <- consensus_coords
    }
    
    # Save the avg_distances for the current repeat amount into the list.
    avg_distances_list[[paste0("repeat_number_", y)]] <- avg_distances
    
    # ================ FLATTEN REPEATS/THRESHOLD LIST (THIS REPEAT AMOUNT)
    
    # Flatten the list of lists containing individual repeat attempts for each threshold
    # After flattening, each element in flattened_list corresponds to a single threshold and
    # is a combined data frame (or matrix) of coordinates from all repeats for that threshold.
    flattened_list <- lapply(ordination_points, function(sublist) {
      do.call(rbind, sublist)
    })
    
    # ================================== NORMALIZATION OF THE COORDINATES
    
    # Normalize all coordinates to unit variance:
    
    # Combine all points across all thresholds to compute global scaling parameters
    all_points <- do.call(rbind, flattened_list)
    
    # Compute global means and sds for each axis
    axis_means <- colMeans(all_points, na.rm = TRUE)
    axis_sds   <- apply(all_points, 2, sd, na.rm = TRUE)
    
    # Define a normalization function
    normalize_coords <- function(mat, means, sds) {
      sweep(sweep(mat, 2, means, "-"), 2, sds, "/")
    }
    
    # Normalize each threshold's ordination points
    scaled_points <- lapply(flattened_list, function(x) normalize_coords(x, axis_means, axis_sds))
    
    # Determine global axis limits after normalization (with some tolerance)
    all_scaled_points <- do.call(rbind, scaled_points)
    x_limits <- range(all_scaled_points[, 1], na.rm = TRUE)
    y_limits <- range(all_scaled_points[, 2], na.rm = TRUE)
    x_limits <- c(x_limits[1] - 0.2, x_limits[2] + 0.2)
    y_limits <- c(y_limits[1] - 0.2, y_limits[2] + 0.2)
    
    # =================================== PLOTTING AND INDEX CALCULATION
    
    for (x in thresholds) {
      # Retrieve normalized data for this threshold
      # scaled_points[[paste0("threshold_", x)]] is a combined matrix with all repeats
      # We need to split by repeats again
      current_key <- paste0("threshold_", x)
      df <- as.data.frame(scaled_points[[current_key]])
      rownames(df) <- NULL
      df <- cbind("sample" = rownames(scaled_points[[current_key]]), df)
      
      # Split into repeats
      n_points <- nrow(df)
      # each threshold run had 'repeats' ordinations, so we split equally
      split_reps <- split(df, ceiling(seq_len(n_points) / (n_points/y)))
      
      # Plot
      # 'consensus_coordinates' should also be normalized.
      norm_consensus <- normalize_coords(consensus_coordinates[[current_key]], axis_means, axis_sds)
      
      step3 <- plot_rep_raref(split_reps, norm_consensus, sample_data(physeq), color = group, group = group, cloud = TRUE, ellipse = TRUE, title = paste0("Threshold: ", x))
      plot <- step3$plot
      
      plots[[paste0("repeat_number ", y)]][[current_key]] <- plot + xlim(x_limits) + ylim(y_limits)
      
      # =========================================== INDEX CALCULATION
      info <- step3$updated_info
      
      # Fix the normalized consensus coordinates adding group variable
      norm_consensus_df <- data.frame
      temp_df <- as.data.frame(norm_consensus)
      temp_df$sample_id <- split_reps[["1"]]$sample
      rownames(temp_df) <- split_reps[["1"]]$sample
      temp_df[[group]] <- info[[group]]
      norm_consensus_df <- rbind(temp_df)
      
      # Extract just position data for the individual points
      just_positions_consensus <- norm_consensus_df[,1:2]
      
      # Determine clusters according to the "group" variable
      clusters_consensus <- unlist(norm_consensus_df[[group]])
      existing_groups_consensus <- unique(clusters_consensus)
      
      # Convert cluster labels to numbers
      cluster_consensus_numbers <- match(clusters_consensus, existing_groups_consensus)
      
      # Calculate index
      index <- clusterSim::index.G1(just_positions_consensus, cluster_consensus_numbers)
      
      index_data <- rbind(index_data, c(toString(y), x, index))
    }
    
  }
  
  # Format the index_data matrix for plotting
  index_data <- as.data.frame(index_data)
  colnames(index_data) <- c("Repeat_Amount", "Threshold", "Index")
  index_data$Threshold <- as.numeric(index_data$Threshold)
  index_data$Index <- as.numeric(index_data$Index)
  # Plot the indexes
  index_plot <- ggplot(index_data, aes(x = Threshold, y = Index, color = Repeat_Amount)) +
    geom_point() +
    labs(x = "Rarefaction Threshold", y = "Calinski-Harabasz pseudo F-statistic") +
    geom_smooth(method = "loess", formula = 'y ~ x')
  
  output <- list(index_plot = index_plot, 
                 ordination_plots = plots,
                 avg_distances = avg_distances_list)
  class(output) <- "test_threshold"
  return(output)
}

#' Print Method for `test_threshold` Objects
#'
#' Prints the index plot from a `test_threshold` object.
#'
#' @param x A `test_threshold` object.
#' @param ... Additional arguments (not used).
#'
#' @export
#' 
print.test_threshold <- function(x, ...) {
  # Print only the index plot
  print(x[["index_plot"]])
  invisible(x)  # standard practice for print methods
}