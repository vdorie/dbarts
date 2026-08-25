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
      "rbart",
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

.onUnload <- function(libpath) {
  ## gc is necessary to collect external pointers who have not yet been collected
  ## that have finalizers pointing to the soon-to-unloaded dll
  gc(FALSE)
  if (is.loaded("dbarts_finalize", PACKAGE = "dbarts")) {
    .Call(C_dbarts_finalize)
    library.dynam.unload("dbarts", libpath)
  }
}
