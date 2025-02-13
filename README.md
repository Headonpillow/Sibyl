# Sybil

Sybil is an R package that provides a convenient way of testing a range of thresholds for rarefaction, evaluating the effect of different values on ordination results.

Installation: 

if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install(version = "3.20")

BiocManager::install("Headonpillow/HLCYG", dependencies = TRUE, force = TRUE)