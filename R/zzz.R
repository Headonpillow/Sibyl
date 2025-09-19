.onLoad <- function(libname, pkgname) {
  # Force rgl to use a headless device in non-interactive/headless environments
  if (identical(Sys.getenv("RGL_USE_NULL"), "true") ||
      !interactive() ||
      !capabilities("X11")) {
    options(rgl.useNULL = TRUE)
  }
}
