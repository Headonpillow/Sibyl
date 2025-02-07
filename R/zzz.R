.onLoad <- function(libname, pkgname) {
  suppressWarnings({
    # Code that triggers warnings
    # Disable "rgl" creation of graphical windows, since we do not require
    # this functionality from "geomorph".
    options(rgl.useNULL = TRUE)
  })
}
