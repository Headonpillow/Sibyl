#' Perform repeated Monte Carlo sampling from the Dirichlet posterior
#'
#' This function performs repeated Monte Carlo sampling from the Dirichlet
#' posterior distribution estimated by `ALDEx2::aldex.clr()` on a count
#' table. Each draw represents one plausible compositional realization of
#' the observed count data, providing an alternative to repeated rarefaction
#' for quantifying compositional uncertainty.
#' @param count A matrix. OTU count table.
#' @param draws An integer. Number of Monte Carlo draws to perform.
#' @param conds A character vector. Grouping variable used by
#' `ALDEx2::aldex.clr()`. If `NULL`, all samples are assigned to a single
#' group.
#' @param cores An integer. Number of cores to use. Values greater than 1
#' enable `ALDEx2::aldex.clr()`'s internal multicore processing (`useMC`).
#' @param ... Additional arguments are reserved to internal use.
#' @return A list containing:
#'   - `dirichlet_matrix_list`: A list of Monte Carlo draw matrices. Each
#'   element is a plausible compositional realization of the input count
#'   table, sampled from the Dirichlet posterior.
#' @importFrom ALDEx2 aldex.clr getDirichletSample
#' @noRd
#' @keywords internal
rep_mc_draws <- function(count, draws = 128, conds = NULL, cores = 2, ...) {
  hidden_args <- list(...)

  if (draws <= 0) {
    stop("draws can't be 0. It needs to be a positive integer")
  }

  # If no grouping variable is provided, assign all samples to one group
  if (is.null(conds)) {
    conds <- rep("all", nrow(count))
  }

  x <- ALDEx2::aldex.clr(
    t(count),
    conds = conds,
    mc.samples = draws,
    useMC = (cores > 1)
  )
  
  # Extract each Dirichlet Monte Carlo instance
  dirichlet_matrices <- lapply(seq_len(draws), function(i) {
    d.mc <- ALDEx2::getDirichletSample(x, i)
    t(d.mc)
  })
  
  # Return list of replicate matrices
  # Each element corresponds to one Monte Carlo draw
  return(invisible(list("dirichlet_matrix_list" = dirichlet_matrices)))
}
