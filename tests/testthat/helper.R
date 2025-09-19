# Make rgl safe on headless machines
Sys.setenv(RGL_USE_NULL = "true")
options(rgl.useNULL = TRUE)
