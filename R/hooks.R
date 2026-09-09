## For the transition release only: the package has two doors under names
## that swapped meaning, so the load says which is which rather than
## leaving a 0.9-x script to find out from a warning mid-fit. A standard
## packageStartupMessage, so suppressPackageStartupMessages silences it.
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "dbarts: 'bart' is the modern front door (formerly 'bart2'); 'bartBT' ",
    "is the BayesTree-style one, with 0.9-x's argument names and defaults. ",
    "A BayesTree-spelled 'bart' call is forwarded to 'bartBT' with a ",
    "warning for this release."
  )
}

.onUnload <- function(libpath) {
  ## gc is necessary to collect external pointers who have not yet been collected
  ## that have finalizers pointing to the soon-to-unloaded dll
  gc(FALSE)
  if (is.loaded("dbarts_finalize", PACKAGE = "dbarts")) {
    .Call(C_dbarts_finalize)
    library.dynam.unload("dbarts", libpath)
  }
}
