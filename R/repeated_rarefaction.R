#' Perform Repeated Rarefaction and Generate an Ordination Plot using PCoA
#'
#' This function performs repeated rarefaction on a phyloseq object, 
#' computes ordination, and generates a PCoA-based visualization.
#' The same procedure is used from `threshold_testing` function when testing 
#' a range of thresholds. 
#'
#' @param input A `phyloseq` object.
#' @param repeats An integer. The number of times to repeat rarefaction. A value of 1 means no repeats. Default = 50.
#' @param threshold An integer. The threshold value to use for rarefaction. Default = 250.
#' @param colorb A string. Column name in the `sample_data()` bundled with the 
#' supplied `phyloseq` object, used to color sample points.
#' @param group A string. Column name in `sample_data()` used for grouping the samples. 
#' The parameter is also used to calculate an ellipse around the points. 
#' The supplied value should coincide with the `group` when using the function
#' `threshold_testing`, since it is on this parameter that also the clustering 
#' performance index is calculated.
#' @param cloud A boolean. If `TRUE`, all repeated data points are shown.
#' Otherwise, only the median points of each sample repetitions cloud are plotted.
#' Default = FALSE.
#' @param ellipse A boolean. If `TRUE`, confidence ellipses around samples are 
#' from the same group are drawn. Default = TRUE.
#' @param cores An integer. Number of cores to use for parallel processing. Default = 4.
#'
#' @return A list containing (While also showing the plot directly):
#'   - `repeats`: Number of repeats.
#'   - `df_consensus_coordinates`: Median ordination positions.
#'   - `df_all`: All ordination positions.
#'   - `plot`: a `ggplot` object.
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
  # Then set it to a separate variable because we need one.
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
