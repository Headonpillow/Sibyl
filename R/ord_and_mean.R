#' Perform Ordination and Compute Consensus Coordinates
#'
#' Computes ordination (PCoA) using Bray-Curtis distance and aligns 
#' results across multiple rarefied datasets.
#' This is the second step of the algorithm.
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
#' @importFrom stats cmdscale
#' 
#' @keywords internal
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