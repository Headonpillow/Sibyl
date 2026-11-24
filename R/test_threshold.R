#' Test different rarefaction thresholds
#'
#' This function uses the same steps of `repeated_rarefaction` in a repeated 
#' fashion and summarizes the results across multiple thresholds.
#' 
#' Clustering performance is evaluated using the Calinski-Harabasz index.
#' The function also calculates the average pairwise distance for each sample 
#' cloud across different rarefaction thresholds, to give a measure of sample 
#' concordance (how stable the sample is at different thresholds).
#'
#' @param input A `phyloseq` object.
#' @param repeats An integer, or a vector of integers. The number of times to 
#' repeat rarefaction. A value of 1 means no repeats. If a vector, 
#' different rarefaction thresholds will be tested sequentially. 
#' @param t_min An integer. The minimum value for the threshold testing range. 
#' @param t_max An integer. The maximum value for the threshold testing range. 
#' @param t_step An integer. The step size for the threshold testing range. 
#' A value between 0 and 1 will cause the same threshold to be tested multiple 
#' times. 
#' @param group A string. Column name in `sample_data()`. Used to group the 
#' samples. The parameter is also used to draw an ellipse around the points. 
#' It the context of this function, is also the value that will be used for the 
#' calculation of the Calinski-Harabasz index.
#' @param cores An integer. Number of cores to use for parallel processing.
#' @param verbose A logical. If `TRUE`, prints messages during the execution.
#' @param ... Additional arguments are reserved to internal use.
#' @return A list containing:
#'   - `index_data`: Numeric CH results in native PCoA coordinates, by depth and repetition count.
#'   - `index_plot`: A `ggplot` object showing a scatterplot of 
#'   Calinski-Harabasz index vs. rarefaction threshold.
#'   - `ordination_plots`: A list of lists containing the ordination plots for 
#'   tested threshold values. If multiple repeat values were used, each one
#'   will be present in a different list.
#'   - `avg_distances`: A list of data frames, named by repetition count,
#'   containing each sample's average pairwise distance at each threshold.
#' @details
#' Higher CH values indicate stronger group separation in consensus coordinates.
#' A plateau alone does not establish low positional uncertainty; inspect APD
#' and retained-sample identities as well. CH uses native, unstandardized scores.
#' The Calinski-Harabasz value takes into account both within cluster and between
#' cluster distances, and in the scope of this function is calculated on the
#' `group` parameter. 
#' CH is missing when fewer than two groups remain or no within-group
#' degrees of freedom are available.
#' A sudden drop in CH value means that samples might have been removed because 
#' they did not have enough reads to reach the threshold. In this case a warning
#' is printed listing the samples which have been removed.
#' In `avg_distances`, if samples do not have enough reads to reach `t_max` value, their
#' APD are set to NA after that threshold and a warning message is printed.
#' @importFrom phyloseq sample_data otu_table sample_data<- subset_samples
#' @importFrom vegan vegdist
#' @importFrom clusterSim index.G1
#' @importFrom dplyr mutate group_by summarise ungroup
#' @importFrom ggplot2 ggplot aes geom_point labs geom_smooth xlim ylim
#' @importFrom utils capture.output
#' @export
#' @examples
#' library(Sibyl)
#' result <- test_threshold(adults, 
#'                          repeats = 10, 
#'                          t_min = 100, 
#'                          t_max = 500, 
#'                          t_step = 50, 
#'                          group = "location", 
#'                          verbose = FALSE, cores = 1)
#' result$index_plot
test_threshold <- function(input, repeats = 50, t_min = 50, t_max = 250, t_step = 5, group = "sample_id", cores = 2, verbose = TRUE, ...) {
  # additional args (reserved for internal use)
  
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
  
  if (!is.numeric(repeats) || any(!is.finite(repeats)) || any(repeats < 1) || any(repeats != floor(repeats))){
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

  # ============= Set up warning collection for samples which do not meet threshold.
  warningCollector <- new.env()
  warningCollector$table <- data.frame(
    Sample_ID = character(0),
    Threshold = numeric(0),
    stringsAsFactors = FALSE
  )
  
  for (y in repeats) {
    # Create lists to store ordination points coordinates for a specific repeat 
    # amount and centroid coordinates.
    ordination_points <- list()
    consensus_coordinates <- list()
    axis_variances <- list()
    # Initialize also a data frame to store pairwise distances of samples point clouds.
    avg_distances <- data.frame(Sample_ID = character(0), Threshold = numeric(0), Avg_Distance = numeric(0))
    
    for (x in thresholds) {
      if (verbose){
        message(paste("Running with", y, "repeats and", x, "threshold"))
      }
      # Perform rarefaction for this threshold and repeat amount
      step1 <- rep_raref(sample_count_matrix(physeq), threshold = x, repeats = y, cores = cores, warning_collector = warningCollector)
      step2 <- ord_and_mean(step1$rarefied_matrix_list, y, distance = "bray", cores = cores)
      axis_variances[[paste0("threshold_", x)]] <- step2$axis_variance
      # ============= AVERAGE PAIRWISE DISTANCE CALCULATION (THIS THRESHOLD)
      
      # Determine which samples are present in the ordination results.
      # (Since the count table is filtered in rep_raref, samples with too low counts will be missing.)
      present_samples <- rownames(step2$aligned_ordinations[[1]])
      # Calculate average pairwise distance for each sample.
      for (sample_id in sample_ids) {
        if (!(sample_id %in% present_samples)) {
          # An excluded sample has no cloud: missing, not zero uncertainty.
          avg_distance <- NA_real_
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
    
    # Keep native PCoA distance units for both plots and CH diagnostics.
    # Axis-wise standardization distorts the relative geometry of the axes.
    points_by_threshold <- flattened_list

    # =================================== PLOTTING AND INDEX CALCULATION
    
    for (x in thresholds) {
      # Retrieve native coordinates for this threshold
      # points_by_threshold[[paste0("threshold_", x)]] is a combined matrix with all repeats
      # We need to split by repeats again
      current_key <- paste0("threshold_", x)
      df <- as.data.frame(points_by_threshold[[current_key]])
      rownames(df) <- NULL
      df <- cbind("sample" = rownames(points_by_threshold[[current_key]]), df)
      
      # Split into repeats
      n_points <- nrow(df)
      # each threshold run had 'repeats' ordinations, so we split equally
      split_reps <- split(df, ceiling(seq_len(n_points) / (n_points/y)))
      
      # Plot
      # Consensus coordinates retain the same native units.
      norm_consensus <- consensus_coordinates[[current_key]]
      
      step3 <- plot_rep_raref(split_reps, norm_consensus, sample_data(physeq), color = group, group = group, cloud = TRUE, ellipse = TRUE, title = paste0("Threshold: ", x), axis_variance = axis_variances[[current_key]])
      plot <- step3$plot
      
      plots[[paste0("repeat_number ", y)]][[current_key]] <- plot
      
      # =========================================== INDEX CALCULATION
      info <- step3$updated_info
      
      # Add group metadata to the native consensus coordinates
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
      
      # CH requires at least two groups and residual within-group degrees of freedom.
      n_groups <- length(unique(cluster_consensus_numbers))
      index <- if (n_groups >= 2 && n_groups < nrow(just_positions_consensus)) {
        clusterSim::index.G1(just_positions_consensus, cluster_consensus_numbers)
      } else {
        NA_real_
      }
      
      index_data <- rbind(index_data, c(toString(y), x, index))
    }
    
  }
  
  # ============= WARNINGS
  if (nrow(warningCollector$table) > 0) {
    pooledTable <- warningCollector$table %>%
      group_by(Sample_ID) %>%
      summarise(Threshold = min(Threshold)) %>%
      ungroup()
    warning("The following samples have been removed during the analysis:\n",
            paste(capture.output(print(as.data.frame(pooledTable))), collapse = "\n"))
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
                 avg_distances = avg_distances_list,
                 index_data = index_data)
  class(output) <- "test_threshold"
  return(output)
}

#' Print Method for `test_threshold` Objects
#'
#' Prints the index plot from a `test_threshold` object.
#'
#' @param x A `test_threshold` object.
#' @param ... Additional arguments (not used).
#' @noRd
#' @export
print.test_threshold <- function(x, ...) {
  # Print only the index plot
  print(x[["index_plot"]])
  invisible(x)  # standard practice for print methods
}