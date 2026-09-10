# Compiles inst/tinytest/capi/consumer.c against the INSTALLED headers with
# R CMD SHLIB and loads it, the way a LinkingTo package does. Shared by the
# test files that drive an entry point through that consumer - it is one
# source with several sets of entry points, and each caller uses its own.
#
# A non-NULL $skip in the returned list is a reason the toolchain produced no
# library, and the CALLER hands it to exit_file: tinytest masks exit_file in
# the test file's own environment, so a call from in here would reach the
# namespace version and silently do nothing. Under CI that same failure is an
# error, since there the toolchain is part of what is under test. Otherwise
# $CALL(name, ...) invokes a registered entry point, and $consumerSource and
# $includeDir serve a caller that compiles further variants of the source.
compileCapiConsumer <- function(prefix, label) {
  consumerSource <- system.file(
    "tinytest",
    "capi",
    "consumer.c",
    package = "dbarts"
  )
  if (consumerSource == "") {
    return(list(skip = "consumer source not installed"))
  }
  buildDir <- tempfile(prefix)
  dir.create(buildDir)
  file.copy(consumerSource, file.path(buildDir, "consumer.c"))

  includeDir <- system.file("include", package = "dbarts")
  headerPath <- file.path(includeDir, "dbarts", "dbarts.h")
  if (!nzchar(includeDir) || !file.exists(headerPath)) {
    msg <- paste0("dbarts.h not found under includeDir '", includeDir, "'")
    if (nzchar(Sys.getenv("CI", ""))) stop(msg)
    return(list(skip = msg))
  }

  # system2's env= is not reliably passed through to the child process on
  # Windows; a Makevars in the build dir is the portable channel for
  # PKG_CPPFLAGS across all platforms including Rtools.
  writeLines(
    sprintf('PKG_CPPFLAGS = -I"%s"', includeDir),
    file.path(buildDir, "Makevars")
  )
  owd <- setwd(buildDir)
  compileOutput <- tryCatch(
    suppressWarnings(system2(
      file.path(R.home("bin"), "R"),
      c("CMD", "SHLIB", "consumer.c"),
      stdout = TRUE,
      stderr = TRUE
    )),
    error = function(e) e
  )
  setwd(owd)

  sharedLib <- file.path(buildDir, paste0("consumer", .Platform$dynlib.ext))
  if (!file.exists(sharedLib)) {
    if (nzchar(Sys.getenv("CI", ""))) {
      stop(
        "could not compile ",
        label,
        " under CI:\n",
        paste(compileOutput, collapse = "\n")
      )
    }
    return(list(skip = paste0("could not compile ", label)))
  }

  dll <- dyn.load(sharedLib)
  list(
    skip = NULL,
    dll = dll,
    buildDir = buildDir,
    consumerSource = consumerSource,
    includeDir = includeDir,
    CALL = function(name, ...) .Call(getNativeSymbolInfo(name, dll), ...)
  )
}
