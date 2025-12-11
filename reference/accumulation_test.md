# Accumulation curve analysis

This function generates accumulation (rarefaction) curves for each
sample in a given `phyloseq` object.

## Usage

``` r
accumulation_test(input, step = 5)
```

## Arguments

- input:

  A `phyloseq` object.

- step:

  A numeric value. The step size for drawing the accumulation curve. It
  influences the rarefaction curve calculation and the granularity of
  points plotted. Default = 5.

## Value

A list containing:

- `accumulation_plot`: A `ggplot` object with all sites faceted.

- `threshold_density`: A `ggplot` density/histogram plot of the 75% ACE
  thresholds.

- `individual_plots`: A list of `ggplot` objects, one per site.

## Details

It fits a general accumulation model using the Abundance Coverage
Estimator (ACE) as an asymptote, and identifies the sequencing depth at
which 75% of the ACE value is reached. It also produces a density plot
showing the distribution of these 75% completion thresholds.

## Examples

``` r
library(Sibyl)
# Creating a smaller subset of the data
adults_sub <- phyloseq::subset_samples(adults, location=="VK3")
# Running accumulation tests on a phyloseq object, higher step size reduces 
# execution time.
accumulation_test(adults_sub, step=50)
```
