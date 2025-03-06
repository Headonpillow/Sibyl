# Sibyl <img src="man/figures/logo.svg" align="right" height="138" /></a>

<!-- badges: start -->
[![Codecov test coverage](https://codecov.io/gh/Headonpillow/Sibyl/graph/badge.svg)](https://app.codecov.io/gh/Headonpillow/Sibyl)
[![R-CMD-check](https://github.com/Headonpillow/Sibyl/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Headonpillow/Sibyl/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**Sibyl** is a package designed to test different rarefaction thresholds when 
normalizing microbial abundance data.  

When performing Principal Component Analysis (PCA) and other types of ordination, 
it is necessary to make sure of choosing a rarefaction threshold which does not 
impact the structure of the data, or our conclusions, during exploratory analysis.  

# Installation

### Quickest - Using Bioconductor Docker images (all OS)

We recommend getting a working bioconductor installation through
[docker images for bioconductor](https://bioconductor.org/help/docker/). We
find the use of containers very useful, allowing for isolation of 
R computing environments.

The bioconductor images are maintained and updated regularly, and already
include all the necessary system dependencies to make **Sibyl** and all its 
dependencies work.

Installation can be accomplished with:

``` r
BiocManager::install("Headonpillow/Sibyl", 
dependencies = TRUE)

```
Which installs the latest development version of **Sibyl**

### Slightly longer - Installing system dependencies manually (Linux only)

For Ubuntu/Debian run:

```
sudo apt update && sudo apt install -y \
    libcurl4-openssl-dev libssl-dev cmake libxml2-dev build-essential \
    liblapack-dev gfortran zlib1g-dev libgl1-mesa-dev libglu1-mesa-dev \
    freeglut3-dev libx11-dev libjpeg-dev libpng-dev
```

Some package dependencies are hosted on Bioconductor. In order to install them 
you will need to install `{BiocManager}`:

``` r
if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install(version = "3.20")
BiocManager::install("remotes")
```

Then it will be possible to install the latest development version 
of **Sibyl** with:

``` r
BiocManager::install("Headonpillow/Sibyl", 
dependencies = TRUE)

```

# Why using Sibyl

The **Sibyl** package was created with the aim of testing the lower limit of
rarefaction thresholds. 

While microbial abundance data is slowly moving away from methods like rarefaction
to account for differing library size, rarefaction has been extensively used 
and still in 2025 is present in much of the literature, and discussed.

When choosing a rarefaction threshold is common to need to operate a compromise 
between describing samples accurately (completeness), and including more samples,
which might sometimes not meet the selected threshold. 

**Sibyl** aims to solve that, allowing users to explore sample completeness, 
and the effect of rarefaction thresholds on ordinations, all with a single tool.
