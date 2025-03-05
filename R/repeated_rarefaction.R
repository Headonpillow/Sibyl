#' Perform repeated rarefaction
#'
#' This function performs repeated rarefaction on a `phyloseq` object, 
#' computes ordination, and generates a PCoA-based visualization.
#' The same procedure is used from `threshold_testing` function when testing 
#' a range of thresholds. 
#'
#' @param input A `phyloseq` object.
#' @param repeats An integer. The number of times to repeat rarefaction. 
#' A value of 1 means no repeats. If too few repeats are selected it would be
#' not possible to draw an ellipse around the group.
#' @param threshold An integer. The threshold value to use for rarefaction. 
#' @param colorb A string. Column name in `sample_data()`. Used to color 
#' sample points.
#' @param group A string. Column name in `sample_data()`. Used to group the 
#' samples. The parameter is also used to draw an ellipse around the points. 
#' @param cloud A boolean. If `TRUE`, all the data points generated from
#' repetitions are shown. Otherwise, only the median points of each sample 
#' repetition cloud are plotted.
#' @param ellipse A boolean. If `TRUE`, confidence ellipses around sample
#' groups are drawn.
#' @param cores An integer. Number of cores to use for parallel processing.
#' @param ... Additional arguments are reserved to internal use.
#' @return A list containing (While also showing the plot directly):
#'   - `repeats`: Number of repeats.
#'   - `df_consensus_coordinates`: A data frame with coordinates of the median
#'   points of the sample clouds.
#'   - `df_all`: A data frame of coordinates ordered by ordination number,
#'   along with metatata.
#'   - `plot`: a `ggplot` object.
#' @importFrom phyloseq sample_data otu_table sample_data<-
#' @importFrom dplyr mutate
#' @importFrom vegan rrarefy vegdist procrustes
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom ggplot2 ggplot aes geom_point stat_ellipse theme_minimal ggtitle xlab ylab
#' @export
#' @examples
#' library(Sibyl)
#' # Running this with cloud = TRUE and ellipse = TRUE will generate a plot 
#' # where the samples belonging to the same group will be colored similarly 
#' # and an ellipse will be drawn around the group.
#' repeated_rarefaction(adults, 
#'                      repeats = 10, 
#'                      threshold = 250, 
#'                      group = "location", 
#'                      colorb = "location", 
#'                      cloud = TRUE, 
#'                      ellipse = TRUE)
#'                      
#' # We can run the function to highlight the spread of the single sample clouds
#' # too, setting the groupb parameter to the sample_id.
#' repeated_rarefaction(adults, 
#'                      repeats = 10, 
#'                      threshold = 250, 
#'                      group = "sample_id", 
#'                      colorb = "location", 
#'                      cloud = TRUE, 
#'                      ellipse = TRUE)
repeated_rarefaction <- function(input, repeats = 50, threshold = 250, colorb="sample_id", group="sample_id", cloud = TRUE, ellipse = FALSE, cores = 2, ...) {
  # list hidden arguments
  hidden_args <- list(...)
  
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
  
  # Extract all sample IDs from the original phyloseq object (if samples are later
  # removed because of unmet rarefaction threshold)
  sample_ids <- rownames(phyloseq::sample_data(input))

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

  if (!(is.double(threshold))){
    stop(paste("Input for threshold: '" ,threshold, "' is not an integer.", sep=""))
  }

  if (repeats <=4 & ellipse == TRUE){
    warning("Too few repeats to draw confidence ellipses. Proceeding with the available data.")
    ellipse <- FALSE
  }

  # Perform the different steps of the repeated rarefaction algorithm
  # Setting the seed might be done for testing purposes
  if(!is.null(hidden_args$seed)){
    step1 <- rep_raref(data.frame(t(otu_table(physeq))), threshold, repeats, cores = cores, seed = hidden_args$seed)
  } else {
      step1 <- rep_raref(data.frame(t(otu_table(physeq))), threshold, repeats, cores = cores)
  }
  step2 <- ord_and_mean(step1$rarefied_matrix_list, repeats, cores = cores)
  step3 <- plot_rep_raref(step2$aligned_ordinations, step2$consensus_coordinates, sample_data(physeq), colorb, group, cloud, ellipse, "Aligned Ordinations with Consensus Overlaid")

  print(step3$plot)

  return(invisible(list("repeats" = repeats, "df_consensus_coordinates" = step3$consensus_df, "df_all" = step3$df_all, "plot" = step3$plot)))
}
