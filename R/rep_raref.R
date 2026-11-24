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
  count <- as.matrix(count_filtered)
  
  # Set up working files
  rarefied_matrices <- list()
  
  per_iter_seeds <- sample.int(.Machine$integer.max, size = repeats)
  rarefy_one <- function(seed, count, threshold) {
    set.seed(seed)
    vegan::rrarefy(count, sample = threshold)
  }
  # Workers need only the explicit arguments and vegan, not the package frame.
  environment(rarefy_one) <- baseenv()
  if (cores == 1) {
    rng_state <- .Random.seed
    on.exit(assign(".Random.seed", rng_state, envir = .GlobalEnv), add = TRUE)
    rarefied_matrices <- lapply(per_iter_seeds, rarefy_one,
                               count = count, threshold = threshold)
  } else {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    rarefied_matrices <- parallel::parLapply(
      cl, per_iter_seeds, rarefy_one, count = count, threshold = threshold)
  }

  return(invisible(list("rarefied_matrix_list"=rarefied_matrices)))
}
