# Sibyl

Sibyl is a package designed to test different rarefaction thresholds when 
normalizing microbial abundance data.  

When performing Principal Component Analysis (PCA) and other types of ordination, 
it is necessary to make sure of choosing a rarefaction threshold which does not 
impact the structure of the data, or our conclusions, during exploratory analysis.  

# Installation

Some package dependencies are hosted on Bioconductor. In order to install them 
you need to install `{BiocManager}`:

``` r
if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install(version = "3.20")
```

You can install the latest version of Sibyl from 
[Github](https://github.com/Headonpillow/Sibyl) with:

``` r
BiocManager::install("Headonpillow/Sibyl", 
dependencies = TRUE, build_vignettes = TRUE, force = TRUE)

```

# Examples



# Why this package was created