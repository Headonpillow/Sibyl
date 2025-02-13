# Sybil

Sybil is an R package that provides a convenient way of testing a range of 
thresholds for rarefaction, evaluating the effect of different values 
on ordination results.

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
