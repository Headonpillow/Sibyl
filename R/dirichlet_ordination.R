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
  step2 <- ord_and_mean2(step1$dirichlet_matrix_list, draws, distance = "aitchison", cores = cores)
  step3 <- plot_rep_raref(step2$aligned_ordinations, step2$consensus_coordinates, sample_data(physeq), colorb, group, cloud, ellipse, "Aligned Ordinations with Consensus Overlaid")
  
  print(step3$plot)
  
  return(invisible(list("draws" = draws, "df_consensus_coordinates" = step3$consensus_df, "df_all" = step3$df_all, "plot" = step3$plot)))
}
