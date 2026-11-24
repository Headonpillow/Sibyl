#' Generate a plot for repeated rarefaction
#'
#' Creates a PCoA-based visualization of rarefied ordinations, displaying consensus mean coordinates.
#' This is the third step of the algorithm.
#' @param aligned_ordinations A list of aligned ordination matrices.
#' @param consensus_coordinates A matrix of consensus coordinates.
#' @param info A data frame containing sample metadata.
#' @param color A string. Column name in `info` for coloring.
#' @param group A string. Column name in `info` for grouping.
#' @param cloud A boolean. If `TRUE`, all repeated points are shown.
#' @param ellipse A boolean. If `TRUE`, descriptive ellipses are drawn.
#' @param title A string. Title of the plot.
#' @param axis_variance Reference PCoA axis percentages, or NULL.
#' @return A list containing:
#'   - `plot`: The generated ggplot2 plot.
#'   - `consensus_df`: A data frame of consensus coordinates.
#'   - `df_all`: A data frame of all ordination positions.
#'   - `updated_info`: Updated metadata (if samples have been removed according
#'   to threshold).
#' @importFrom ggplot2 ggplot aes geom_point stat_ellipse theme_minimal 
#'   ggtitle xlab ylab theme element_text
#' @noRd
#' @keywords internal
plot_rep_raref <- function(aligned_ordinations, consensus_coordinates, info,
                           color, group, cloud, ellipse, title,
                           axis_variance = NULL) {
  as_coordinates <- function(x) {
    if (is.matrix(x)) {
      data.frame(sample = rownames(x), Dim1 = x[, 1], Dim2 = x[, 2])
    } else {
      names(x)[1:3] <- c("sample", "Dim1", "Dim2")
      x
    }
  }
  reps <- lapply(aligned_ordinations, as_coordinates)
  sample_ids <- reps[[1]]$sample
  # Match metadata by sample identity, rather than relying on incidental order.
  info <- methods::as(info, "data.frame")[sample_ids, , drop = FALSE]
  aligned_df <- do.call(rbind, lapply(seq_along(reps), function(i) {
    z <- reps[[i]]
    names(z)[1] <- "sample_id"
    z$ordination <- paste0("Ordination", i)
    z[[color]] <- info[z$sample_id, color]
    z[[group]] <- info[z$sample_id, group]
    z
  }))
  rownames(aligned_df) <- NULL
  consensus_df <- as.data.frame(consensus_coordinates)
  names(consensus_df) <- c("Dim1", "Dim2")
  consensus_df$sample_id <- rownames(consensus_coordinates)
  consensus_df[[color]] <- info[consensus_df$sample_id, color]
  consensus_df[[group]] <- info[consensus_df$sample_id, group]

  plot <- ggplot2::ggplot()
  if (cloud) {
    plot <- plot + ggplot2::geom_point(
      data = aligned_df,
      ggplot2::aes(x = Dim1, y = Dim2, color = .data[[color]]),
      alpha = 0.25, size = 1)
  }
  if (ellipse) {
    plot <- plot + ggplot2::stat_ellipse(
      data = aligned_df,
      ggplot2::aes(x = Dim1, y = Dim2, color = .data[[color]],
                   group = .data[[group]]),
      linewidth = 0.6, na.rm = TRUE)
  }
  # Hollow diamonds distinguish consensus points without a fixed red overlay.
  plot <- plot + ggplot2::geom_point(
    data = consensus_df,
    ggplot2::aes(x = Dim1, y = Dim2, color = .data[[color]]),
    shape = 23, fill = "white", size = 2.5, stroke = 0.8) +
    ggplot2::scale_color_viridis_d(end = 0.85) +
    ggplot2::coord_equal() + ggplot2::theme_minimal() +
    ggplot2::labs(title = title, color = color,
                  x = "PCoA 1", y = "PCoA 2")
  if (!is.null(axis_variance)) {
    plot <- plot + ggplot2::labs(
      x = sprintf("PCoA 1 (%.1f%%)", axis_variance[1]),
      y = sprintf("PCoA 2 (%.1f%%)", axis_variance[2]))
  }
  invisible(list(plot = plot, consensus_df = consensus_df,
                 df_all = aligned_df, updated_info = info))
}
