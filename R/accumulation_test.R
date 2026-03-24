#' Accumulation curve analysis
#'
#' This function generates accumulation (rarefaction) curves for each sample in 
#' a given `phyloseq` object.
#' 
#' It fits a general accumulation model using the Abundance 
#' Coverage Estimator (ACE) as an asymptote, and identifies the sequencing 
#' depth at which 75% of the ACE value is reached. It also produces a density 
#' plot showing the distribution of these 75% completion thresholds.
#'
#' Sites that fail model fitting or quality-control criteria are excluded from
#' downstream analyses, and their identities are returned for inspection.
#' 
#' @param input A `phyloseq` object.
#' @param step  A numeric value. The step size for drawing the accumulation curve.
#' It influences the rarefaction curve calculation and the granularity of points plotted.
#' Default = 5.
#' @return A list containing:
#'   - `accumulation_plot`: A `ggplot` object with all sites faceted.
#'   - `threshold_density`: A `ggplot` density/histogram plot of the 75% ACE thresholds.
#'   - `individual_plots`: A list of `ggplot` objects, one per site.
#'   - `nls_failed_sites`: A character vector of site names for which the nonlinear
#'     model failed to converge.
#'   - `quality_failed_sites`: A character vector of site names for which the model
#'     converged but failed quality-control criteria (e.g. implausible threshold or
#'     poor fit).
#'   - `fitted_table`: A data frame containing per-site model results and diagnostics,
#'     including fitted parameters, threshold estimates, and quality-control flags.
#' @importFrom phyloseq otu_table subset_samples
#' @importFrom vegan rarecurve estimateR
#' @importFrom dplyr left_join mutate group_by select rename tibble
#' @importFrom magrittr %>%
#' @importFrom tidyr nest unnest pivot_wider
#' @importFrom purrr map map2 map_lgl
#' @importFrom broom tidy augment
#' @importFrom stats setNames density
#' @importFrom ggplot2 ggplot aes geom_point geom_line scale_color_manual 
#'   geom_vline facet_wrap theme_minimal labs geom_histogram geom_density
#'   labeller after_stat
#' @export
#' @examples
#' library(Sibyl)
#' # Creating a smaller subset of the data
#' adults_sub <- phyloseq::subset_samples(adults, location=="VK3")
#' # Running accumulation tests on a phyloseq object, higher step size reduces 
#' # execution time.
#' accumulation_test(adults_sub, step=50)
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
  safe_nls <- purrr::possibly(
    function(.x) {
      nls(
        # Force asymptote = ACE
        # So the formula is: Species ~ (ACE * Sample) / (b + Sample)
        Species ~ (ACE * Sample) / (b + Sample),
        data = .x,
        start = list(b = median(.x$Sample, na.rm = TRUE)),
        control = nls.control(maxiter = 50)
      )
    },
    # Unsuccessful fitting attempts are null and are removed
    otherwise = NULL
  )

  df_fitted <- df_nested %>%
    dplyr::mutate(
      fit = purrr::map(data, safe_nls),
      fit_ok = !purrr::map_lgl(fit, is.null)
    )
  
  # 2B) Quality filtering the fits to remove sub-optimal ones
  
  df_fitted <- df_fitted %>%
    dplyr::mutate(
      # Extract the fitted b parameter from each successful nls model.
      # If the fit failed (NULL), store NA.
      b_est = purrr::map_dbl(
        fit,
        ~ if (is.null(.x)) NA_real_ else stats::coef(.x)[["b"]]
      ),
      
      # Rule 1:
      # The sample size corresponding to 75% of the asymptote for the
      # Michaelis-Menten model is 3 * b.
      sample_75 = 3 * b_est,
      # Get the largest observed sample size for each Site.
      # This is used to check whether the estimated threshold is
      # unrealistically far beyond the observed range.
      max_sample = purrr::map_dbl(
        data,
        ~ max(.x$Sample, na.rm = TRUE)
      ),
      
      # Rule 2:
      # Compute RMSE between observed Species values and model predictions.
      # If the fit failed, store NA.
      rmse = purrr::map2_dbl(
        data, fit,
        ~ if (is.null(.y)) {
          NA_real_
        } else {
          pred <- stats::predict(.y, newdata = .x)
          sqrt(mean((.x$Species - pred)^2, na.rm = TRUE))
        }
      ),
      # Get the maximum observed Species count for each Site.
      # This is used to scale RMSE into a normalized RMSE (NRMSE),
      # so fit quality is judged relative to the size of the curve.
      max_species = purrr::map_dbl(
        data,
        ~ max(.x$Species, na.rm = TRUE)
      ),
      # Normalized RMSE:
      # smaller values indicate a better fit relative to the observed curve size.
      nrmse = rmse / max_species
    )
  
  # Global robust depth window based on the central part of site depths
  q1 <- stats::quantile(df_fitted$max_sample, 0.25, na.rm = TRUE)
  q3 <- stats::quantile(df_fitted$max_sample, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  global_upper <- q3 + 1.5 * iqr
  
  # Final fit acceptance and diagnostics
  # For each Site, we evaluate whether the fitted model is usable based on:
  #   (1) convergence of the nonlinear fit,
  #   (2) plausibility of the estimated threshold relative to the site's data,
  #   (3) plausibility relative to the global depth distribution,
  #   (4) overall goodness-of-fit (NRMSE).
  
  df_fitted <- df_fitted %>%
    dplyr::mutate(
      # Model convergence (safe_nls returns NULL if fitting fails)
      converged = !purrr::map_lgl(fit, is.null),
      # Rule 1 (local plausibility):
      # The estimated 75% threshold (3 * b) should not be far beyond the
      # observed range of the data for that site.
      pass_local_range = sample_75 <= 1.5 * max_sample,
      # Rule 2 (global plausibility):
      # The estimated threshold should also fall within a reasonable global
      # range of sequencing depths across all sites and not be an "outlier".
      # The upper bound is defined using the IQR of per-site maximum depths.
      pass_global_range = sample_75 <= global_upper,
      # Rule 3 (fit quality):
      # The normalized RMSE evaluates how well the model reproduces the observed
      # accumulation curve. Lower values indicate better fit.
      pass_nrmse = nrmse <= 0.1,
      
      # Final decision:
      # A fit is accepted only if:
      #   - the model converged, AND
      #   - all quality-control criteria are satisfied
      #
      # We also assign a readable label about failure reason.
      fit_ok = converged & pass_local_range & pass_global_range & pass_nrmse,
      failure_reason = dplyr::case_when(
        !converged ~ "nls_failed_to_converge",
        !pass_local_range ~ "threshold_above_local_range",
        !pass_global_range ~ "threshold_above_global_range",
        !pass_nrmse ~ "poor_fit_nrmse",
        TRUE ~ "passed"
      )
    )
  
  # Create a vector catching up the samples failing the NLS fitting and the 
  # Quality control for the other parameters
  nls_failed_sites <- df_fitted %>%
    dplyr::filter(purrr::map_lgl(fit, is.null)) %>%
    dplyr::pull(Site)
  
  quality_failed_sites <- df_fitted %>%
    dplyr::filter(!purrr::map_lgl(fit, is.null), !fit_ok) %>%
    dplyr::pull(Site)
  
  # Print warnings with the sample names
  if (length(nls_failed_sites) > 0) {
    warning(
      paste0(
        length(nls_failed_sites),
        " site(s) were skipped because the accumulation model failed to converge:\n",
        paste(nls_failed_sites, collapse = "\n")
      ),
      call. = FALSE
    )
  }
  
  if (length(quality_failed_sites) > 0) {
    warning(
      paste0(
        length(quality_failed_sites),
        " site(s) were skipped because the fitted accumulation curve did not pass quality-control criteria:\n",
        paste(quality_failed_sites, collapse = "\n"),
        "\nResults and suggested thresholds are based on the remaining samples."
      ),
      call. = FALSE
    )
  }
  
  if (!any(df_fitted$fit_ok)) {
    stop(
      "No sites passed model fitting and quality-control criteria. No accumulation thresholds could be estimated.",
      call. = FALSE
    )
  }
  
  # 3) Create a new data frame of Sample values for each site
  # We are considering only the fitted curves, plus only the 
  # Ones that survived quality control
  df_predicted <- df_fitted %>%
    dplyr::filter(fit_ok) %>%
    dplyr::mutate(
      newdata = purrr::map(data, ~ dplyr::tibble(
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
    dplyr::filter(fit_ok) %>%
    dplyr::mutate(
      params = purrr::map(fit, broom::tidy)
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
      linewidth = 1
    ) +
    # Vertical dotted lines at 75% threshold
    geom_vline(
      data = plateau_lines,
      aes(xintercept = sample_75),
      linetype = "dotted",
      color = "black",
      linewidth = 1
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
        linewidth = 1
      ) +
      # Vertical line for 75% threshold
      geom_vline(
        xintercept = param_row$sample_75,
        linetype = "dotted",
        color = "black",
        linewidth = 1
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
                 individual_plots = individual_plots,
                 nls_failed_sites = nls_failed_sites,
                 quality_failed_sites = quality_failed_sites,
                 fitted_table = df_fitted
                 )
  class(output) <- "accumulation_test"
  return(output)
}

#' Print Method for `accumulation_test` Object
#'
#' @param x An `accumulation_test` object.
#' @param ... Additional arguments (not used).
#' @noRd
#' @export
print.accumulation_test <- function(x, ...) {
  # Print only the accumulation plot
  print(x[["accumulation_plot"]])
  invisible(x)  # standard practice for print methods
}
