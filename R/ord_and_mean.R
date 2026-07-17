#' Perform ordination and compute consensus coordinates
#'
#' Computes ordination (PCoA) using a configurable dissimilarity measure and
#' aligns results across multiple replicate compositional matrices (e.g.
#' rarefied count tables from `rep_raref()` or Dirichlet Monte Carlo draws
#' from `rep_mc_draws()`) via Procrustes rotation. This is the second step
#' shared by the repeated rarefaction and Dirichlet-based ordination
#' algorithms.
#' @param matrix_list A list of replicate count/compositional tables.
#' @param replicates An integer. The number of replicates (rarefaction repeats
#' or Monte Carlo draws).
#' @param distance A string. The dissimilarity measure passed to
#' `vegan::vegdist()`, e.g. `"bray"` for rarefied counts or `"aitchison"`
#' for CLR-transformed compositional data.
#' @param cores An integer. The number of cores to use.
#' @return A list containing:
#'   - `aligned_ordinations`: List of aligned ordinations.
#'   - `consensus_coordinates`: Consensus coordinates using Procrustes alignment.
#' @importFrom vegan vegdist procrustes
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom stats cmdscale
#' @noRd
#' @keywords internal
ord_and_mean <- function(matrix_list, replicates, distance = "bray", cores = 2) {
  
  #========================= ordinations and plots generation
  
  # Initialize a list to store ordinations
  ordinations <- list()
  
  # Set up parallel backend
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  # Perform parallel computation using foreach
  results <- foreach(
    i = 1:length(matrix_list), 
    .packages = c('vegan', 'stats')
    ) %dopar% {
    
    # Calculate distance
    dist_matrix <- vegdist(matrix_list[[i]], method = distance)
    # Perform PCoA (Principal Coordinates Analysis)
    ordination <- cmdscale(dist_matrix, k = 2)
    
    # # Convert ordination result to data frame
    # ord_df <- as.data.frame(ordination)
    # colnames(ord_df) <- c("PCoA1", "PCoA1")
    # ord_df$Sample <- rownames(ordination)
    
    # Ensure rownames from the original DF are transfered
    rownames(ordination) <- rownames(matrix_list[[i]])
    
    # Return a list containing the ordination and the plot
    list(ordination = ordination)
  }
  
  # Stop the cluster after computation
  suppressWarnings(stopCluster(cl))
  
  # Extract ordinations and plots from the results
  for (i in seq_along(results)) {
    ordinations[[i]] <- results[[i]]$ordination
  }
  
  #================================ procrustes
  
  # Perform Procrustes analysis to align all ordinations to the first one
  aligned_ordinations <- lapply(
    ordinations, 
    function(x) {
      fit <- procrustes(ordinations[[1]], x, scale = FALSE)
      y <- fit$Yrot
      rownames(y) <- rownames(ordinations[[1]])
      colnames(y) <- colnames(ordinations[[1]])
      y
    }
  )
  
  # Convert list to array for consensus calculation
  aligned_array <- array(
    unlist(aligned_ordinations), 
    dim = c(
      nrow(aligned_ordinations[[1]]), 
      ncol(aligned_ordinations[[1]]), 
      length(aligned_ordinations)
      )
    )
  
  # Compute consensus using mean shape (mean across the 3rd dimension: specimens)
  consensus_coords <- apply(aligned_array, c(1, 2), mean, na.rm = TRUE)
  rownames(consensus_coords) <- rownames(aligned_ordinations[[1]])
  colnames(consensus_coords) <- colnames(aligned_ordinations[[1]])
  
  return(
    invisible(list(
    "aligned_ordinations" = aligned_ordinations, 
    "consensus_coordinates" = consensus_coords
    ))
  )
}