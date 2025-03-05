#' Generate a plot for repeated rarefaction
#'
#' Creates a PCoA-based visualization of rarefied ordinations, calculating median cluster points.
#' This is the third step of the algorithm.
#' @param aligned_ordinations A list of aligned ordination matrices.
#' @param consensus_coordinates A matrix of consensus coordinates.
#' @param info A data frame containing sample metadata.
#' @param color A string. Column name in `info` for coloring.
#' @param group A string. Column name in `info` for grouping.
#' @param cloud A boolean. If `TRUE`, all repeated points are shown.
#' @param ellipse A boolean. If `TRUE`, confidence ellipses are drawn.
#' @param title A string. Title of the plot.
#' @return A list containing:
#'   - `plot`: The generated ggplot2 plot.
#'   - `consensus_df`: A data frame of consensus coordinates.
#'   - `df_all`: A data frame of all ordination positions.
#'   - `updated_info`: Updated metadata (if samples have been removed according
#'   to threshold).
#' @importFrom ggplot2 ggplot aes geom_point stat_ellipse theme_minimal 
#' ggtitle xlab ylab theme element_text
#' @noRd
#' @keywords internal
plot_rep_raref <- function(aligned_ordinations, consensus_coordinates, info, color, group, cloud, ellipse, title) {
  
  # =========== Handle missing samples 
  # threshold which might have been removed from step1:
  # Extract sample names from info
  data_sample_names <- rownames(info)

  # Format the input, which differs if using `repeated_rarefaction` or `test_threshold`
  if (is.numeric(aligned_ordinations[[1]])) {
    aligned_ordinations <- lapply(aligned_ordinations, function(x) {
      sample_names <- rownames(x) 
      data.frame(sample = sample_names, V1 = x[,1], V2 = x[,2])
    })
  }
  
  # Extract row names from the first data frame in aligned_ordinations
  ordination_sample_names <- aligned_ordinations[[1]]$sample
  
  # Identify samples to remove using exact matching
  samples_to_remove <- data_sample_names[!data_sample_names %in% ordination_sample_names]
  
  # Remove the extra samples from info
  if (length(samples_to_remove) > 0) {
    info <- info[!(rownames(info) %in% samples_to_remove), , drop = FALSE]
  }
  
  # Combine aligned ordinations into one data frame for plotting
  aligned_df <- data.frame()
  
  for (i in 1:length(aligned_ordinations)) {
    temp_df <- as.data.frame(aligned_ordinations[[i]])
    colnames(temp_df) <- c("sample_id", "Dim1", "Dim2")
    temp_df$ordination <- paste0("Ordination", i)
    temp_df[[color]] <- info[[color]]
    temp_df[[group]] <- info[[group]]
    aligned_df <- rbind(aligned_df, temp_df)
  }
  
  # Remove rownames from the df since they became redundant
  rownames(aligned_df) <- NULL
  
  # Convert consensus_coordinates to a data frame
  consensus_df <- as.data.frame(consensus_coordinates)
  colnames(consensus_df) <- c("Dim1", "Dim2")
  consensus_df$sample_id <- aligned_ordinations[["1"]]$sample
  
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
          linetype = 1, lwd = 0.8,
          na.rm = TRUE
        )
    } else {
      # Ellipses in grey if more than 6 categories
      # Still grouping by the color variable to get multiple ellipses if there are multiple groups
      plot <- plot +
        stat_ellipse(
          data = aligned_df,
          aes(x = Dim1, y = Dim2, group = .data[[group]]),
          color = "grey70",
          linetype = 1, lwd = 0.8,
          na.rm = TRUE
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