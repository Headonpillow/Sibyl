.onLoad <- function(libname, pkgname) {
    # Disable "rgl" creation of graphical windows, since we do not require
    # this functionality from "geomorph".
    options(rgl.useNULL = TRUE)
}
