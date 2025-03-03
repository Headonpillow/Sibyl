#' Perform Repeated Rarefaction
#'
#' This function performs repeated rarefaction on a count table. 
#' This is the first step of the algorithm.
#'
#' @param count A matrix. OTU count table.
#' @param threshold An integer. The threshold for rarefaction.
#' @param repeats An integer. Number of repeats.
#'
#' @return A list containing:
#'   - `rarefied_matrix_list`: A list of rarefied matrices.
#'
#' @importFrom vegan rrarefy
#' @export
#' 
rep_raref <- function(count, threshold, repeats, ...) {
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
    # Check if a warning collector has been setup, add the sample to it.
    if (!is.null(hidden_args$warning_collector)) {
      hidden_args$warning_collector$warnings <- c(hidden_args$warning_collector$warnings, sample_names)
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
  
  # Perform repeated rarefaction and store the normalized results in a list
  for (i in 1:repeats) {
    rarefied_count <- rrarefy(count, sample = threshold)
    rarefied_matrices[[i]] <- rarefied_count
  }
  
  return(invisible(list("rarefied_matrix_list"=rarefied_matrices)))
}
