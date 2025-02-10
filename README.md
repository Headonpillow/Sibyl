# repeated_rarefaction

HLCYG (How low can you go) is an R package that provides a convenient way of testing a range of thresholds for rarefation, evaluating the effect of different values on ordination results, and integrating with commonly used packages in microbial ecology.

Installation guide and documentation can be found here: https://repeated-rarefaction.readthedocs.io/en/latest/#

Installation: 

if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install(version = "3.20")

if (!require("BiocManager")) {
  install.packages("BiocManager")
}

BiocManager::install("Headonpillow/HLCYG", dependencies = TRUE, force = TRUE)


On unix based systems loading the package will output a warning: 
Warning messages:
1: In rgl.init(initValue, onlyNULL) : RGL: unable to open X11 display
2: 'rgl.init' failed, running with 'rgl.useNULL = TRUE'. 

This does not affect the package behaviour. 