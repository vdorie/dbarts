# as_draws_array/as_draws_df are posterior's generics; posterior is
# Suggests-only, so the methods register dynamically whenever posterior's
# namespace loads - before or after dbarts, and without forcing the load
# here. envir locates the generic (posterior's namespace, not dbarts'):
# registerS3method stores the method in the generic's own S3 table.
.onLoad <- function(libname, pkgname) {
  registerPosteriorMethods <- function(...) {
    ns <- asNamespace("posterior")
    classes <- c(
      "bart",
      "bartMultinomial",
      "bartOrdinal",
      "bartNegbin",
      "bartHurdle"
    )
    for (class in classes) {
      registerS3method(
        "as_draws_array",
        class,
        get(paste0("as_draws_array.", class)),
        envir = ns
      )
      registerS3method(
        "as_draws_df",
        class,
        get(paste0("as_draws_df.", class)),
        envir = ns
      )
    }
  }
  setHook(packageEvent("posterior", "onLoad"), registerPosteriorMethods)
  if (isNamespaceLoaded("posterior")) registerPosteriorMethods()
}

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
