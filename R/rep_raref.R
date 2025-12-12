#' Perform repeated rarefaction
#'
#' This function performs repeated rarefaction on a count table. 
#' @param count A matrix. OTU count table.
#' @param threshold An integer. The threshold for rarefaction.
#' @param repeats An integer. Number of repeats.
#' @return A list containing:
#'   - `rarefied_matrix_list`: A list of rarefied matrices.
#' @importFrom vegan rrarefy
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @noRd
#' @keywords internal
rep_raref <- function(count, threshold, repeats, cores = 2, ...) {
  hidden_args <- list(...)
  
  if (repeats == 0) {
    warning("repeats can't be 0. It needs to be a positive integer. Performs rarefaction without repetition.")
  }
  
  if (repeats < 0) {
    warning("repeats can't be negative. It needs to be a positive integer. Performs rarefaction without repetition.")
  }
  
  # Identify rows that do not meet the threshold
  row_totals <- rowSums(count)
  below_threshold <- which(row_totals < threshold)
  
  # Print a warning with the problematic sample indices
  if (length(below_threshold) > 0) {
    # Retrieve sample names below threshold
    sample_names <- names(below_threshold)
    warning_table <- data.frame(
      Sample_ID = sample_names,
      Threshold = threshold,
      stringsAsFactors = FALSE
      )
    # Check if a warning collector has been setup, add the sample to it.
    if (!is.null(hidden_args$warning_collector)) {
      hidden_args$warning_collector$table <- rbind(hidden_args$warning_collector$table, warning_table)
    }
    else{
      warning("The following samples have row sums less than ", threshold, " and have been removed: ", paste(sample_names, collapse = ", "))
    }
  }
  
  # Remove problematic samples
  count_filtered <- count[row_totals >= threshold, , drop = FALSE]
  count <- count_filtered
  
  # Set up working files
  rarefied_matrices <- list()
  
  # Set up parallel backend
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  # generate per-iteration seeds in the master process
  # This is randomly sampling integers up to the maximum integer size
  # This ensures different seeds across different runs (creates as many 
  # unique seeds as repeats)
  per_iter_seeds <- sample.int(.Machine$integer.max, size = repeats)

  # perform rarefactions: each worker sets its own seed from per_iter_seeds
  # to ensure reproducibility and avoid weird correlations
  rarefied_matrices <- foreach(i = 1:repeats, .packages = "vegan") %dopar% {
    set.seed(per_iter_seeds[i])
    rrarefy(count, sample = threshold)
  }
  
  # Stop the cluster
  suppressWarnings(stopCluster(cl))
  
  return(invisible(list("rarefied_matrix_list"=rarefied_matrices)))
}
