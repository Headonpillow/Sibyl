#' Plot Average Pairwise Distance for Sample Clusters Across Rarefaction Thresholds
#'
#' This function draws a plot of `avg_distances` generated from the 
#' `threshold_testing` function. 
#'
#' @param avg_distances A dataframe of average distances obtained from `threshold_testing`.
#'
#' @return A ggplot2 object showing average pairwise distances across thresholds.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point facet_wrap labs theme_minimal
#' @importFrom stats dist
#' @export
#' 
avg_pairwise_dist_plot <- function(avg_distances) {
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
  
  return(plot)
}