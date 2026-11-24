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
  ordinate <- function(counts, distance) {
    # Zero-only features contribute nothing to Bray-Curtis dissimilarities.
    if (identical(distance, "bray")) {
      counts <- counts[, colSums(counts) > 0, drop = FALSE]
    }
    d <- vegan::vegdist(counts, method = distance)
    fit <- stats::cmdscale(d, k = 2, eig = TRUE)
    rownames(fit$points) <- rownames(counts)
    list(ordination = fit$points, eigenvalues = fit$eig)
  }
  # Do not serialize the enclosing frame and its full replicate list to workers.
  environment(ordinate) <- baseenv()
  if (cores == 1) {
    results <- lapply(matrix_list, ordinate, distance = distance)
  } else {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    results <- parallel::parLapply(cl, matrix_list, ordinate, distance = distance)
  }
  ordinations <- lapply(results, `[[`, "ordination")
  aligned <- lapply(ordinations, function(x) {
    y <- vegan::procrustes(ordinations[[1]], x, scale = FALSE)$Yrot
    dimnames(y) <- dimnames(ordinations[[1]])
    y
  })
  consensus <- Reduce(`+`, aligned) / length(aligned)
  # Labels refer to the first (reference) configuration, not the rotated
  # replicate axes or a supposed eigendecomposition of the mean coordinates.
  eigenvalues <- results[[1]]$eigenvalues
  axis_variance <- 100 * pmax(eigenvalues[1:2], 0) / sum(abs(eigenvalues))
  invisible(list(aligned_ordinations = aligned,
                 consensus_coordinates = consensus,
                 axis_variance = axis_variance,
                 reference_eigenvalues = eigenvalues))
}
