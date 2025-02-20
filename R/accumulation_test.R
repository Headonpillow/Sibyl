#' Accumulation Curve Analysis
#'
#' This function generates accumulation (rarefaction) curves for each sample in 
#' a given `phyloseq` object, fits a general accumulation model using the Abundance 
#' Coverage Estimator (ACE) as an asymptote, and identifies the sequencing 
#' depth at which 75% of the ACE value is reached. It also produces a density 
#' plot showing the distribution of these 75% completion thresholds.
#'
#' @param input A `phyloseq` object.
#' @param step  A numeric value. The step size for drawing the accumulation curve.
#' Default = 5.
#'
#' @return A list containing:
#'   - `accumulation_plot`: A ggplot2 object with all sites faceted.
#'   - `threshold_density`: A ggplot2 density/histogram plot of the 75% ACE thresholds.
#'   - `individual_plots`: A list of ggplot2 objects, one per site.
#'
#' @importFrom phyloseq otu_table
#' @importFrom vegan rarecurve estimateR
#' @importFrom dplyr left_join mutate group_by select rename tibble
#' @importFrom magrittr %>%
#' @importFrom tidyr nest unnest pivot_wider
#' @importFrom purrr map map2
#' @importFrom broom tidy augment
#' @importFrom stats setNames density
#' @importFrom ggplot2 ggplot aes geom_point geom_line scale_color_manual 
#'   geom_vline facet_wrap theme_minimal labs geom_histogram geom_density
#'   labeller after_stat
#'   
#' @export
#' 
accumulation_test <- function(input, step = 5) {
  # Transform the phyloseq to a df of counts
  counts <- as.data.frame(t(otu_table(input)))
  
  # Create a df out of rarefaction curves
  df <- rarecurve(counts, step = step, tidy = TRUE)
  
  # Estimate ACE from observed species
  res <- estimateR(as.matrix(counts))
  ace_df <- data.frame(
    Site = colnames(res),
    ACE  = res["S.ACE", ]   # row "S.ACE" contains the ACE estimate
  )
  df_joined <- dplyr::left_join(df, ace_df, by = "Site")
  
  # 1) Group by Site and nest the data
  df_nested <- df_joined %>%
    dplyr::group_by(Site) %>%
    tidyr::nest()
  
  # 2) Fit a generic accumulation model to each Site's data
  df_fitted <- df_nested %>%
    dplyr::mutate(
      fit = purrr::map(data, ~ nls(
        # Force asymptote = ACE
        # So the formula is: Species ~ (ACE * Sample) / (b + Sample)
        Species ~ (ACE * Sample) / (b + Sample),
        data = .x,
        start = list(b = median(.x$Sample, na.rm = TRUE))
      ))
    )
  
  # 3) Create a new data frame of Sample values for each site
  df_predicted <- df_fitted %>%
    dplyr::mutate(
      newdata = purrr::map(data, ~ tibble::tibble(
        Sample = seq(min(.x$Sample), max(.x$Sample), length.out = 100),
        ACE    = unique(.x$ACE)  # each site's ACE
      )),
      preds = purrr::map2(fit, newdata, ~ broom::augment(.x, newdata = .y))
    ) %>%
    dplyr::select(Site, preds) %>%
    tidyr::unnest(cols = preds)
  
  # Rename and keep needed columns
  df_final <- df_predicted %>%
    dplyr::rename(Species = .fitted) %>%
    dplyr::select(Site, Sample, Species)
  
  # Gather parameter estimates from nls
  params_by_site <- df_fitted %>%
    dplyr::mutate(
      params = purrr::map(fit, ~ broom::tidy(.x))
    ) %>%
    dplyr::select(Site, params) %>%
    tidyr::unnest(cols = params) %>%
    dplyr::left_join(ace_df, by = "Site")
  
  # Convert parameters into wide format
  params_wide <- params_by_site %>%
    dplyr::select(Site, term, estimate, ACE) %>%
    tidyr::pivot_wider(names_from = term, values_from = estimate) %>%
    dplyr::mutate(
      sample_75 = 3 * b # If the formula for 0.75 * ACE solves to 3*b
    )
  
  # Create a label vector for the faceted plot
  params_labeled <- params_wide %>%
    dplyr::mutate(
      label = paste0(
        Site, "\n",
        "ACE = ", round(ACE, 2), ", ",
        "75% = ", round(sample_75, 2)
      )
    )
  labels_vector <- setNames(params_labeled$label, params_labeled$Site)
  
  # We'll create a small data frame with just Site and sample_75 for vlines
  plateau_lines <- params_wide %>%
    dplyr::select(Site, sample_75)
  
  # ------------------------------------------------------------------------------
  # (1) Combined Facet Plot
  # ------------------------------------------------------------------------------
  accumulation_plot <- ggplot() +
    # Original data (Observed)
    geom_point(
      data = df,
      aes(x = Sample, y = Species, group = Site, color = "Observed"),
      alpha = 0.7
    ) +
    # Fitted lines
    geom_line(
      data = df_final,
      aes(x = Sample, y = Species, group = Site, color = "Fitted"),
      size = 1
    ) +
    # Vertical dotted lines at 75% threshold
    geom_vline(
      data = plateau_lines,
      aes(xintercept = sample_75),
      linetype = "dotted",
      color = "black",
      size = 1
    ) +
    # Color legend
    scale_color_manual(
      name = "Data Type",
      values = c("Observed" = "red", "Fitted" = "blue")
    ) +
    # Facet with custom labels
    facet_wrap(~ Site, scales = "free_x", labeller = labeller(Site = labels_vector)) +
    theme_minimal() +
    labs(
      title = "Species Accumulation Curves",
      x = "Sample Size",
      y = "Species Count"
    )
  
  # ------------------------------------------------------------------------------
  # (2) Threshold Density Plot
  # ------------------------------------------------------------------------------
  threshold_density <- ggplot(params_wide, aes(x = sample_75)) +
    geom_histogram(
      aes(y = after_stat(density)),
      bins = 30,
      fill = "skyblue",
      color = "skyblue",
      alpha = 0.6
    ) +
    geom_density(
      aes(y = after_stat(density)),
      color = "maroon",
      linewidth = 1
    ) +
    theme_minimal() +
    labs(
      x = "75% ACE Threshold",
      y = "Density",
      title = "Distribution of 75% ACE Thresholds"
    )
  
  # ------------------------------------------------------------------------------
  # (3) List of Individual Plots (one per Site)
  # ------------------------------------------------------------------------------
  unique_sites <- unique(df_final$Site)
  
  individual_plots <- lapply(unique_sites, function(st) {
    # Filter the data for this Site
    site_points <- df %>% dplyr::filter(Site == st)
    site_fitted <- df_final %>% dplyr::filter(Site == st)
    param_row   <- params_wide %>% dplyr::filter(Site == st)
    
    # Build a single-site plot
    p <- ggplot() +
      # Observed data
      geom_point(
        data = site_points,
        aes(x = Sample, y = Species),
        color = "red",
        alpha = 0.7
      ) +
      # Fitted line
      geom_line(
        data = site_fitted,
        aes(x = Sample, y = Species),
        color = "blue",
        size = 1
      ) +
      # Vertical line for 75% threshold
      geom_vline(
        aes(xintercept = param_row$sample_75),
        linetype = "dotted",
        color = "black",
        size = 1
      ) +
      theme_minimal() +
      labs(
        title = paste("Accumulation Curve for Sample:", st),
        subtitle = paste0(
          "ACE = ", round(param_row$ACE, 2),
          " | 75% = ", round(param_row$sample_75, 2)
        ),
        x = "Sample Size",
        y = "Species Count"
      )
    
    p
  })
  
  output <- list(accumulation_plot = accumulation_plot,
                 threshold_density = threshold_density,
                 individual_plots = individual_plots)
  class(output) <- "accumulation_test"
  return(output)
}

#' Print Method for `accumulation_test` Object
#'
#' @param x An `accumulation_test` object.
#' @param ... Additional arguments (not used).
#'
#' @export
#' 
print.accumulation_test <- function(x, ...) {
  # Print only the accumulation plot
  print(x[["accumulation_plot"]])
  invisible(x)  # standard practice for print methods
}
