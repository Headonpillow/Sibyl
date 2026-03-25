rep_mc_draws <- function(count, draws = 128, conds = NULL, cores = 2, ...) {
  hidden_args <- list(...)
  
  if (draws <= 0) {
    stop("repeats can't be 0. It needs to be a positive integer")
  }
  
  # If no grouping variable is provided, assign all samples to one group
  if (is.null(conds)) {
    conds <- rep("all", nrow(count))
  }
  
  x <- ALDEx2::aldex.clr(
    t(count),
    conds = conds,
    mc.samples = draws
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
