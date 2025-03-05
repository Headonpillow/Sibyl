#' Plot average pairwise distance across thresholds
#'
#' This function it's a plotting function for `avg_distances` generated from 
#' \link{test_threshold}. 
#' @param avg_distances A dataframe of average distances obtained from `test_threshold`.
#' @return A `ggplot` object showing average pairwise distances across thresholds.
#' @importFrom phyloseq subset_samples
#' @importFrom ggplot2 ggplot aes geom_line geom_point facet_wrap labs theme_minimal
#' @importFrom stats dist
#' @export
#' @examples
#' library(Sibyl)
#' # Creating a smaller subset of the data
#' adults_sub <- phyloseq::subset_samples(adults, location=="VK3")
#' result <- test_threshold(adults_sub, 
#'                          repeats = 10, 
#'                          t_min = 100, 
#'                          t_max = 500, 
#'                          t_step = 50, 
#'                          group = "location",
#'                          verbose = FALSE)
#' avg_pairwise_dist_plot(result$avg_distances$repeat_number_10)
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