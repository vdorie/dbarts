# Compiles inst/tinytest/capi/consumer.c against the INSTALLED headers with
# R CMD SHLIB and loads it, the way a LinkingTo package does. Shared by the
# test files that drive an entry point through that consumer - it is one
# source with several sets of entry points, and each caller uses its own.
#
# A non-NULL $skip in the returned list means no toolchain was found to build
# with - no make, or no C compiler - and the CALLER hands it to exit_file:
# tinytest masks exit_file in the test file's own environment, so a call from
# in here would reach the namespace version and silently do nothing. That is
# the only skip. A missing source or header, and a consumer that does not
# compile, are failures of the shipped header or the package and stop, on
# CRAN as anywhere else. Otherwise $CALL(name, ...) invokes a registered
# entry point, and $consumerSource and $includeDir serve a caller that
# compiles further variants of the source (compileCapiSource).

# Whether R CMD SHLIB has a toolchain to run: make (MAKE, as R CMD SHLIB
# reads it) and the C compiler R was configured with, each found as a file or
# on the PATH. A compiler named with flags ("clang -arch arm64") or behind a
# launcher ("ccache gcc") is judged by its first word. NULL when both are
# there, otherwise the reason.
capiToolchainMissing <- function() {
  firstWord <- function(command) {
    strsplit(trimws(command), "[[:space:]]+")[[1L]][1L]
  }
  found <- function(program) {
    !is.na(program) &&
      nzchar(program) &&
      (file.exists(program) || nzchar(Sys.which(program)))
  }
  make <- firstWord(Sys.getenv("MAKE", "make"))
  if (!found(make)) {
    return(paste0("no make found ('", make, "')"))
  }
  compiler <- tryCatch(
    suppressWarnings(system2(
      file.path(R.home("bin"), "R"),
      c("CMD", "config", "CC"),
      stdout = TRUE,
      stderr = TRUE
    )),
    error = function(e) character()
  )
  compiler <- if (length(compiler) > 0L) firstWord(compiler[1L]) else NA
  if (!found(compiler)) {
    return(paste0("no C compiler found ('", compiler, "')"))
  }
  NULL
}

# Builds sourceName in dir with R CMD SHLIB and returns the library's path,
# stopping with the compiler's output when it does not build.
compileCapiSource <- function(dir, sourceName, label) {
  owd <- setwd(dir)
  on.exit(setwd(owd))
  output <- tryCatch(
    suppressWarnings(system2(
      file.path(R.home("bin"), "R"),
      c("CMD", "SHLIB", sourceName),
      stdout = TRUE,
      stderr = TRUE
    )),
    error = function(e) conditionMessage(e)
  )
  lib <- file.path(
    dir,
    paste0(tools::file_path_sans_ext(sourceName), .Platform$dynlib.ext)
  )
  if (!file.exists(lib)) {
    stop("could not compile ", label, ":\n", paste(output, collapse = "\n"))
  }
  lib
}

compileCapiConsumer <- function(prefix, label) {
  missingTool <- capiToolchainMissing()
  if (!is.null(missingTool)) {
    return(list(skip = missingTool))
  }
  consumerSource <- system.file(
    "tinytest",
    "capi",
    "consumer.c",
    package = "dbarts"
  )
  if (consumerSource == "") {
    stop("the C API consumer source is not installed")
  }
  includeDir <- system.file("include", package = "dbarts")
  headerPath <- file.path(includeDir, "dbarts", "dbarts.h")
  if (!nzchar(includeDir) || !file.exists(headerPath)) {
    stop("dbarts.h not found under includeDir '", includeDir, "'")
  }

  buildDir <- tempfile(prefix)
  dir.create(buildDir)
  file.copy(consumerSource, file.path(buildDir, "consumer.c"))
  # system2's env= is not reliably passed through to the child process on
  # Windows; a Makevars in the build dir is the portable channel for
  # PKG_CPPFLAGS across all platforms including Rtools.
  writeLines(
    sprintf('PKG_CPPFLAGS = -I"%s"', includeDir),
    file.path(buildDir, "Makevars")
  )
  sharedLib <- compileCapiSource(buildDir, "consumer.c", label)

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
