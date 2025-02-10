#' Accumulation Curve Analysis
#'
#' This function generates accumulation (rarefaction) curves for each sample in 
#' a given `phyloseq` object, fits a general accumulation model using the Abundance 
#' Coverage Estimator (ACE) as an asymptote, and identifies the sequencing 
#' depth at which 75% of the ACE value is reached. It also produces a density 
#' plot showing the distribution of these 75% completion thresholds.
#'
#' @param input A `phyloseq` object containing microbial community data.
#' @param step An integer specifying the step size used for rarefaction.
#'
#' @return A list containing:
#'   - `accumulation_plot`: A ggplot2 object visualizing the accumulation curves.
#'   - `threshold_density`: A ggplot2 density plot of the 75% ACE threshold.
#'
#' @importFrom phyloseq otu_table
#' @importFrom vegan rarecurve estimateR
#' @importFrom dplyr left_join mutate group_by select rename
#' @importFrom magrittr %>%
#' @importFrom tidyr nest unnest pivot_wider
#' @importFrom purrr map map2
#' @importFrom broom tidy augment
#' @importFrom stats setNames density
#' @importFrom ggplot2 ggplot aes geom_point geom_line scale_color_manual 
#'   geom_vline facet_wrap theme_minimal labs geom_histogram geom_density
#'   labeller
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
  df_joined <- df %>%
    left_join(ace_df, by = "Site")
  
  # 1) Group by SITE and nest the data
  df_nested <- df_joined %>%
    group_by(Site) %>%
    nest()
  
  # 2) Fit the Michaelis–Menten model to each SITE's data
  df_fitted <- df_nested %>%
    mutate(
      fit = map(data, ~nls(
        # Force asymptote = ACE
        # So the formula is: Species ~ (ACE * Sample) / (b + Sample)
        Species ~ (ACE * Sample) / (b + Sample),
        data = .x,
        start = list(b = median(.x$Sample, na.rm = TRUE))  # guess for b
      ))
    )
  
  # 3) Create a new data frame of SAMPLE values over which we want predictions (e.g., a smooth curve)
  df_predicted <- df_fitted %>%
    mutate(
      newdata = map(data, ~tibble(
        Sample = seq(min(.x$Sample), max(.x$Sample), length.out = 100),
        ACE    = unique(.x$ACE)  # each site's ACE
      )),
      preds = map2(fit, newdata, ~augment(.x, newdata = .y))
    ) %>%
    select(Site, preds) %>%
    unnest(cols = preds)
  
  df_final <- df_predicted %>%
    rename(Species = .fitted) %>%
    select(Site, Sample, Species)
  
  params_by_site <- df_fitted %>%
    mutate(
      params = map(fit, ~ broom::tidy(.x))  # tidy() gives a small tibble with term= "b", estimate= ...
    ) %>%
    select(Site, params) %>%
    unnest(cols = params) %>%              # unnest to get one row per parameter per Site
    left_join(ace_df, by = "Site")         # bring in the ACE column
  
  # params_by_site has columns: SITE, term, estimate (and maybe std.error, etc.)
  # We want to pivot so that 'term' (which is 'a' or 'b') becomes a column
  params_wide <- params_by_site %>%
    select(Site, term, estimate, ACE) %>%
    pivot_wider(names_from = term, values_from = estimate)
  
  params_labeled <- params_wide %>%
    mutate(
      # Format your labels as you like
      label = paste0(
        Site, "\n",
        "ACE = ", round(ACE, 2), ", ",
        "75% = ", round(b*3, 2)
      )
    )
  
  params_wide <- params_wide %>%
    mutate(
      sample_95=3*b  # from the math above
    )
  
  # Suppose SITE is a character (or factor). We map each SITE to its label
  labels_vector <- setNames(params_labeled$label, params_labeled$Site)
  
  # We'll create a small data frame with just SITE and sample_95
  plateau_lines <- params_wide %>%
    select(Site, sample_95)
  
  accumulation_plot <- ggplot() +
    # 1) Original data (points) in red
    geom_point(
      data = df,
      aes(
        x = Sample,
        y = Species,
        group = Site,            # group by SITE so each site's points connect properly if lines were used
        color = "Observed"       # a fixed label to distinguish in the legend
      ),
      alpha = 0.7
    ) +
    # 2) Fitted data (lines) in blue
    geom_line(
      data = df_final,
      aes(
        x = Sample,
        y = Species,
        group = Site,
        color = "Fitted"         # a fixed label to distinguish in the legend
      ),
      size = 1
    ) +
    # 3) Manually set colors for these two labels
    scale_color_manual(
      name = "Data Type",        # Legend title
      values = c("Observed" = "red", "Fitted" = "blue")
    ) +
    # 4) Vertical dotted line at sample_95 for each SITE
    geom_vline(
      data = plateau_lines,
      aes(xintercept = sample_95),
      linetype = "dotted",
      color = "black",
      size = 1
    ) +
    # 5) Facet if desired
    facet_wrap(
      ~ Site,
      scales = "free_x",
      labeller = labeller(Site = labels_vector)
    ) +
    theme_minimal() +
    labs(
      title = "Species Accumulation Curves",
      x = "Sample Size",
      y = "Species Count"
    )
  
  #####
  params_by_site <- params_by_site %>%
    mutate(three_b = 3 * estimate)
  
  # ------------------------------------------------------------------------------
  # 2) Plot a histogram of 'three_b' with a density curve, flipping coordinates
  #    so that 'three_b' is on the Y axis.
  # ------------------------------------------------------------------------------
  threshold_density <- ggplot(params_wide, aes(x = sample_95)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "skyblue", color = "skyblue", alpha = 0.6) +
    # Density line: also map y=..density.., same variable on x
    geom_density(aes(y = after_stat(density)), color = "maroon", linewidth = 1) +
    theme_minimal() +
    labs(
      x = "75%",
    )
  
  return(list(accumulation_plot = accumulation_plot, threshold_density = threshold_density))
}
