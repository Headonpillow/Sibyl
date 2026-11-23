#' Plot average pairwise distance across thresholds
#'
#' This function it's a plotting function for `avg_distances` generated from 
#' \link{test_threshold}. 
#' @param avg_distances A dataframe of average distances obtained from `test_threshold`.
#' @param show_legend Logical. Show the redundant sample legend? Default FALSE.
#' @return A `ggplot` object showing average pairwise distances across thresholds.
#' @importFrom phyloseq subset_samples
#' @importFrom ggplot2 ggplot aes geom_line geom_point facet_wrap labs theme_minimal
#' @importFrom stats dist
#' @export
#' @examples
#' library(Sibyl)
#' # Creating a smaller subset of the data
#' adults_sub <- adults
#' result <- test_threshold(adults_sub, 
#'                          repeats = 10, 
#'                          t_min = 100, 
#'                          t_max = 500, 
#'                          t_step = 50, 
#'                          group = "location",
#'                          verbose = FALSE, cores = 1)
#' avg_pairwise_dist_plot(result$avg_distances$repeat_number_10)
avg_pairwise_dist_plot <- function(avg_distances, show_legend = FALSE) {
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
    theme_minimal() +
    ggplot2::scale_color_viridis_d(end = 0.85) +
    ggplot2::theme(legend.position = if (show_legend) "right" else "none")
  
  return(plot)
}