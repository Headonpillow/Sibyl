#' Perform Dirichlet Monte Carlo ordination
#'
#' This function performs repeated Monte Carlo sampling from the Dirichlet
#' posterior distribution estimated by `ALDEx2::aldex.clr()` on a `phyloseq`
#' object, computes ordination using Aitchison distances, and generates a
#' PCoA-based visualization. It follows the same downstream pipeline as
#' `repeated_rarefaction()` (alignment by Procrustes rotation and consensus
#' coordinates), replacing repeated rarefaction with a probabilistic model
#' of compositional uncertainty. This enables a direct comparison of
#' uncertainty arising from stochastic subsampling versus probabilistic
#' compositional modeling, while preserving an identical downstream analysis
#' pipeline.
#'
#' @param input A `phyloseq` object.
#' @param draws An integer. The number of Monte Carlo draws to sample from
#' the Dirichlet posterior. If too few draws are selected it would not be
#' possible to draw an ellipse around the group.
#' @param colorb A string. Column name in `sample_data()`. Used to color
#' sample points.
#' @param group A string. Column name in `sample_data()`. Used to group the
#' samples, and to condition the Dirichlet sampling performed by
#' `ALDEx2::aldex.clr()`. The parameter is also used to draw an ellipse
#' around the points.
#' @param cloud A boolean. If `TRUE`, all the data points generated from the
#' Monte Carlo draws are shown. Otherwise, only the median points of each
#' sample draw cloud are plotted.
#' @param ellipse A boolean. If `TRUE`, confidence ellipses around sample
#' groups are drawn.
#' @param cores An integer. Number of cores to use for parallel processing.
#' @param ... Additional arguments are reserved to internal use.
#' @return A list containing (While also showing the plot directly):
#'   - `draws`: Number of Monte Carlo draws.
#'   - `df_consensus_coordinates`: A data frame with coordinates of the median
#'   points of the sample clouds.
#'   - `df_all`: A data frame of coordinates ordered by ordination number,
#'   along with metadata.
#'   - `plot`: a `ggplot` object.
#' @importFrom phyloseq sample_data otu_table sample_data<-
#' @importFrom vegan procrustes
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom ggplot2 ggplot aes geom_point stat_ellipse theme_minimal ggtitle xlab ylab
#' @export
#' @examples
#' library(Sibyl)
#' \donttest{
#' # Running this with cloud = TRUE and ellipse = TRUE will generate a plot
#' # where the samples belonging to the same group will be colored similarly
#' # and an ellipse will be drawn around the group.
#' dirichlet_ordination(adults,
#'                      draws = 10,
#'                      group = "location",
#'                      colorb = "location",
#'                      cloud = TRUE,
#'                      ellipse = TRUE)
#' }
dirichlet_ordination <- function(input, draws = 128, colorb="sample_id", group="sample_id", cloud = TRUE, ellipse = FALSE, cores = 2, ...) {
  # additional args (reserved for internal use)
  
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
                  dirichlet_ordination needs an existing column to color samples by.", sep=""))
  }
  if (!(group %in% names(sample_data(physeq)))) {
    stop(paste("'",group,"' is not a column name in the sample information in the inputed phyloseq object.
                  dirichlet_ordination needs an existing column to group samples by.", sep=""))
  }
  
  if (!(is.double(draws))){
    stop(paste("Input for draws: '", draws, "' is not an integer.", sep=""))
  }
  
  if (draws <=4 & ellipse == TRUE){
    warning("Too few MC draws to draw confidence ellipses. Proceeding with the available data.")
    ellipse <- FALSE
  }
  
  # Grab a vector of conditions by the group variable
  conds = as.character(sample_data(physeq)[[group]])
  
  # Perform the different steps of the dirichlet_ordination algorithm
  # Perform montecarlo draws from dirichlet distribution (each replicate handled inside rep_montecarlo_draws)
  step1 <- rep_mc_draws(data.frame(t(otu_table(physeq))), draws, conds = conds, cores = cores)
  step2 <- ord_and_mean(step1$dirichlet_matrix_list, draws, distance = "aitchison", cores = cores)
  step3 <- plot_rep_raref(step2$aligned_ordinations, step2$consensus_coordinates, sample_data(physeq), colorb, group, cloud, ellipse, "Aligned Ordinations with Consensus Overlaid")
  
  print(step3$plot)
  
  return(invisible(list("draws" = draws, "df_consensus_coordinates" = step3$consensus_df, "df_all" = step3$df_all, "plot" = step3$plot)))
}
