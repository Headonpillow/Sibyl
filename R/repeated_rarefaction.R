#' Perform Repeated Rarefaction and Generate an Ordination Plot
#'
#' This function performs repeated rarefaction on a phyloseq object, 
#' computes ordination, and generates a PCA-based visualization.
#'
#' @param input A `phyloseq` object.
#' @param repeats An integer. The number of rarefaction repeats. Default = 50.
#' @param threshold An integer. The threshold value for rarefaction. Default = 250.
#' @param colorb A string. Column name in `sample_data()` used for point colors.
#' @param group A string. Column name in `sample_data()` used for grouping and ellipse calculation.
#' @param cloud A boolean. If `TRUE`, all repeated data points are shown. Default = FALSE.
#' @param ellipse A boolean. If `TRUE`, confidence ellipses are drawn. Default = TRUE.
#' @param cores An integer. Number of cores to use for parallel processing. Default = 4.
#'
#' @return A list containing:
#'   - `repeats`: Number of repeats.
#'   - `df_consensus_coordinates`: Median ordination positions.
#'   - `df_all`: All ordination positions.
#'   - `plot`: The generated ordination plot.
#'
#' @importFrom phyloseq sample_data otu_table
#' @importFrom dplyr mutate
#' @importFrom vegan rrarefy vegdist procrustes
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom ggplot2 ggplot aes geom_point stat_ellipse theme_minimal ggtitle xlab ylab
#' @export
#' 
repeated_rarefaction <- function(input, repeats = 50, threshold = 250, colorb="sample_id", group="sample_id", cloud = TRUE, ellipse = FALSE, cores = 4) {
  # Check if input is a Phyloseq object
  if (inherits(input, "phyloseq")) {
    physeq <- input
  } else {
    stop("Input must be a Phyloseq object, including a count table and sample data.")
  }
  
  # Make the rownames of the Phyloseq object a new "sample_id" variable for the sample data.
  # (this covers the case in which no sample_id column is present in the sample data)
  # Then set it to a separate variable.
  sample_data(physeq)$sample_id <- rownames(sample_data(physeq))

  # ============ Checks and warnings

  if (!(colorb %in% names(sample_data(physeq)))) {
    stop(paste("'",colorb,"' is not a column name in the sample information in the inputed phyloseq object.
                  repeated_rarefaction needs an existing column to color samples by.", sep=""))
  }
  if (!(group %in% names(sample_data(physeq)))) {
    stop(paste("'",group,"' is not a column name in the sample information in the inputed phyloseq object.
                  repeated_rarefaction needs an existing column to group samples by.", sep=""))
  }

  if (!(is.double(repeats))){
    stop(paste("Input for repeats: '", repeats, "' is not an integer.", sep=""))
  }

  if (!(repeats == round(repeats))){
    stop(paste("Input for repeats: '", repeats, "' is not an integer.", sep=""))
  }

  if (!(is.double(threshold))){
    stop(paste("Input for threshold: '" ,threshold, "' is not an integer.", sep=""))
  }

  if (repeats <=4 & ellipse == TRUE){
    warning("Too few repeats to draw confidence ellipses. Proceeding with the available data.")
    ellipse <- FALSE
  }

  # Perform the different steps of the repeated rarefaction algorithm
  step1 <- rep_raref(data.frame(t(otu_table(physeq))), threshold, repeats)
  step2 <- ord_and_mean(step1$rarefied_matrix_list, repeats, cores)
  step3 <- plot_rep_raref(step2$aligned_ordinations, step2$consensus_coordinates, sample_data(physeq), colorb, group, cloud, ellipse, "Aligned Ordinations with Consensus Overlaid")

  print(step3$plot)

  return(invisible(list("repeats" = repeats, "df_consensus_coordinates" = step3$consensus_df, "df_all" = step3$df_all, "plot" = step3$plot)))
}

#' Perform Repeated Rarefaction
#'
#' This function performs repeated rarefaction on a count table.
#'
#' @param count A matrix. OTU count table.
#' @param threshold An integer. The threshold for rarefaction.
#' @param repeats An integer. Number of repeats.
#'
#' @return A list containing:
#'   - `rarefied_matrix_list`: A list of rarefied matrices.
#'
#' @importFrom vegan rrarefy
#' @export
#' 
rep_raref <- function(count, threshold, repeats) {
  if (repeats == 0) {
    warning("repeats can't be 0. It needs to be a positive integer. Performs rarefaction without repetition.")
  }

  if (repeats < 0) {
    warning("repeats can't be negative. It needs to be a positive integer. Performs rarefaction without repetition.")
  }

  # Identify rows that do not meet the threshold
  row_totals <- rowSums(count)
  below_threshold <- which(row_totals < threshold)
  # Print a warning with the problematic sample indices
  if (length(below_threshold) > 0) {
    warning("The following samples have row sums less than ", threshold, " and have been removed: ", paste(names(below_threshold), collapse = ", "))
  }
  # Remove problematic samples
  count_filtered <- count[row_totals >= threshold, , drop = FALSE]
  count <- count_filtered

  # Set up working files
  rarefied_matrices <- list()

  # Perform repeated rarefaction and store the normalized results in a list
  for (i in 1:repeats) {
    rarefied_count <- rrarefy(count, sample = threshold)
    rarefied_matrices[[i]] <- rarefied_count
  }

  return(invisible(list("rarefied_matrix_list"=rarefied_matrices)))
}

#' Perform Ordination and Compute Consensus Coordinates
#'
#' Computes ordination using Bray-Curtis distance and aligns results across multiple rarefied datasets.
#'
#' @param rarefied_matrix_list A list of rarefied count tables.
#' @param repeats An integer. The number of rarefaction repeats.
#' @param cores An integer. The number of cores to use. Default = 4.
#'
#' @return A list containing:
#'   - `aligned_ordinations`: List of aligned ordinations.
#'   - `consensus_coordinates`: Consensus coordinates using Procrustes alignment.
#'
#' @importFrom vegan vegdist procrustes
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom geomorph mshape
#' @export
ord_and_mean <- function(rarefied_matrix_list, repeats, cores = 4) {

  #========================= ordinations and plots generation

  # Initialize a list to store ordinations
  ordinations <- list()

  # Set up parallel backend
  cl <- makeCluster(cores)
  registerDoParallel(cl)

  # Perform parallel computation using foreach
  results <- foreach(i = 1:length(rarefied_matrix_list), .packages = c('vegan', 'ggplot2')) %dopar% {
    # Calculate Bray-Curtis distance
    dist_matrix <- vegdist(rarefied_matrix_list[[i]], method = "bray")

    # Perform PCoA (Principal Coordinates Analysis)
    ordination <- cmdscale(dist_matrix, k = 2)

    # Convert ordination result to data frame
    ord_df <- as.data.frame(ordination)
    colnames(ord_df) <- c("PCoA1", "PCoA2")
    ord_df$Sample <- rownames(ordination)

    # Return a list containing the ordination and the plot
    list(ordination = ordination)

  }

  # Stop the cluster after computation
  stopCluster(cl)


  # Extract ordinations and plots from the results
  for (i in 1:length(results)) {
    ordinations[[i]] <- results[[i]]$ordination
  }

  #================================ procrustes

  # Perform Procrustes analysis to align all ordinations to the first one
  aligned_ordinations <- lapply(ordinations, function(x) procrustes(ordinations[[1]], x)$Yrot)

  # Convert list to array for consensus calculation
  aligned_array <- array(unlist(aligned_ordinations), dim = c(nrow(aligned_ordinations[[1]]), ncol(aligned_ordinations[[1]]), length(aligned_ordinations)))

  # Compute consensus using mean shape
  consensus_coords <- mshape(aligned_array)

  return(invisible(list("aligned_ordinations" = aligned_ordinations, "consensus_coordinates" = consensus_coords)))
}


#' Generate a Plot for Repeated Rarefaction
#'
#' Creates a PCA-based visualization of rarefied ordinations, highlighting consensus points.
#'
#' @param aligned_ordinations A list of aligned ordination matrices.
#' @param consensus_coordinates A matrix of consensus coordinates.
#' @param info A data frame containing sample metadata.
#' @param color A string. Column name in `info` for coloring.
#' @param group A string. Column name in `info` for grouping.
#' @param cloud A boolean. If `TRUE`, all repeated points are shown. Default = FALSE.
#' @param ellipse A boolean. If `TRUE`, confidence ellipses are drawn. Default = TRUE.
#' @param title A string. Title of the plot.
#'
#' @return A list containing:
#'   - `plot`: The generated ggplot2 plot.
#'   - `consensus_df`: A data frame of consensus coordinates.
#'   - `df_all`: A data frame of all ordination positions.
#'   - `updated_info`: Updated metadata.
#'
#' @importFrom ggplot2 ggplot aes geom_point stat_ellipse theme_minimal ggtitle xlab ylab
#' @export
#' 
plot_rep_raref <- function(aligned_ordinations, consensus_coordinates, info, color, group, cloud, ellipse, title) {
  
  # This code handles the occurrence of samples below the rarefaction threshold
  # which might have been removed from step1. 
  # Extract sample names from info
  data_sample_names <- rownames(info)
  # Extract row names from the first data frame in aligned_ordinations
  ordination_sample_names <- rownames(aligned_ordinations[[1]])
  # Identify samples to remove from info
  samples_to_remove <- setdiff(data_sample_names, ordination_sample_names)
  # Remove the extra samples from info
  if (length(samples_to_remove) > 0) {
    info <- info[!(rownames(info) %in% samples_to_remove), , drop = FALSE]
  }
  
  # Combine aligned ordinations into one data frame for plotting
  aligned_df <- data.frame()

  for (i in 1:length(aligned_ordinations)) {
    temp_df <- as.data.frame(aligned_ordinations[[i]])
    colnames(temp_df) <- c("Dim1", "Dim2")
    temp_df$sample_id <- rownames(aligned_ordinations[[i]])
    temp_df$ordination <- paste0("Ordination", i)
    temp_df[[color]] <- info[[color]]
    temp_df[[group]] <- info[[group]]
    aligned_df <- rbind(aligned_df, temp_df)
  }

  # Convert consensus_coordinates to a data frame
  consensus_df <- as.data.frame(consensus_coordinates)
  colnames(consensus_df) <- c("Dim1", "Dim2")
  consensus_df$sample_id <- rownames(aligned_ordinations[[1]])

  # Determine how many colors are needed to plot the levels of the color variable
  num_levels <- length(unique(consensus_df[[color]]))

  plot <- ggplot()

  if (num_levels <= 6) {
    # Use variable-based coloring if there are 6 or fewer categories
    plot <- plot +
      geom_point(
        data = aligned_df,
        aes(x = Dim1, y = Dim2, color = .data[[color]]),
        alpha = 0.3
      )
  } else {
    # Use a fixed grey color for all points if there are more than 6 categories
    plot <- plot +
      geom_point(
        data = aligned_df,
        aes(x = Dim1, y = Dim2),
        color = "grey70",
        alpha = 0.3
      )
  }

  if (!cloud) {
    plot$layers <- plot$layers[-1]
  }

  if (ellipse) {
    if (num_levels <= 6) {
      # Ellipses colored by the color variable
      plot <- plot +
        stat_ellipse(
          data = aligned_df,
          aes(x = Dim1, y = Dim2, color = .data[[color]], group = .data[[group]]),
          linetype = 1, lwd = 0.8
        )
    } else {
      # Ellipses in grey if more than 6 categories
      # Still grouping by the color variable to get multiple ellipses if there are multiple groups
      plot <- plot +
        stat_ellipse(
          data = aligned_df,
          aes(x = Dim1, y = Dim2, group = .data[[group]]),
          color = "grey70",
          linetype = 1, lwd = 0.8
        )
    }
  }

  # Plot all aligned ordinations (consensus points)
  plot <- plot +
    geom_point(data = consensus_df, aes(x = Dim1, y = Dim2), color = "red", size = 2) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(title) +
    xlab("Dimension 1") +
    ylab("Dimension 2")

  return(invisible(list("plot" = plot, "consensus_df" = consensus_df, "df_all" = aligned_df, "updated_info" = info)))
}

#' HLCYG_physeq_data
#'
#' This is a dataset containing stuff.
#'
#' @name HLCYG_physeq_data
#' 
#' @format A data frame with 53940 rows and 10 variables:
#' \describe{
#'   \item{location}{sample location}
"HLCYG_physeq_data"
