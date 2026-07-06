# Loads (building if necessary) the bartcore R shim and provides a
# bart()-compatible fit function for the equivalence harness.

bartcoreShimPath <- function() {
  root <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
  file.path(root, "tests", "cpp")
}

loadBartcoreShim <- function() {
  dir <- bartcoreShimPath()
  soFile <- file.path(dir, paste0("rshim", .Platform$dynlib.ext))
  sources <- c(
    file.path(dir, "rshim.cpp"),
    list.files(
      file.path(dirname(dir), "..", "src", "bartcore"),
      full.names = TRUE,
      pattern = "\\.hpp$"
    )
  )
  if (
    !file.exists(soFile) ||
      any(file.mtime(sources) > file.mtime(soFile))
  ) {
    message("building bartcore shim...")
    root <- dirname(dirname(dir))
    env <- c(
      paste0(
        "PKG_CPPFLAGS=",
        shQuote(paste0(
          "-I",
          file.path(root, "src", "include"),
          " -I",
          file.path(root, "src")
        ))
      ),
      "PKG_CXXFLAGS=-std=gnu++20",
      paste0(
        "PKG_LIBS=",
        shQuote(paste(
          file.path(root, "src", "misc.a"),
          file.path(root, "src", "external.a"),
          file.path(root, "src", "rc.a")
        ))
      )
    )
    status <- system2(
      file.path(R.home("bin"), "R"),
      c("CMD", "SHLIB", "-o", soFile, file.path(dir, "rshim.cpp")),
      env = env
    )
    if (status != 0L) stop("failed to build bartcore shim")
  }
  dyn.load(soFile)
  invisible(soFile)
}

# Matches bart()'s defaults and prior derivations; returns the fields the
# equivalence harness consumes (yhat.test, sigma, varcount as ndpost x .).
# dart: FALSE, or a list overriding any of alpha, update.alpha, a, b, rho.
bartcore_bart <- function(
  x.train,
  y.train,
  x.test = NULL,
  weights = NULL,
  ntree = 200L,
  ndpost = 1000L,
  nskip = 100L,
  k = 2.0,
  numcut = 100L,
  usequants = FALSE,
  sigdf = 3.0,
  sigquant = 0.90,
  splitprobs = NULL,
  dart = FALSE
) {
  binary <- length(unique(y.train)) == 2L

  sigest <- 1.0
  sigScale <- 1.0
  nodeScale <- if (binary) 3.0 else 0.5
  if (!binary) {
    sigest <- summary(lm(y.train ~ x.train, weights = weights))$sigma
    sigScale <- qchisq(1.0 - sigquant, sigdf) / sigdf
  }

  dartParams <- NULL
  if (!identical(dart, FALSE)) {
    spec <- if (is.list(dart)) dart else list()
    default <- function(name, value) {
      if (!is.null(spec[[name]])) spec[[name]] else value
    }
    dartParams <- as.double(c(
      default("alpha", 1.0),
      default("update.alpha", TRUE),
      default("a", 0.5),
      default("b", 1.0),
      default("rho", 0.0)
    ))
  }

  result <- .Call(
    "bartcore_fit",
    x.train,
    as.double(y.train),
    x.test,
    if (is.null(weights)) NULL else as.double(weights),
    NULL,
    binary,
    sigest,
    sigdf,
    sigScale,
    as.integer(ntree),
    k,
    nodeScale,
    as.integer(numcut),
    as.logical(usequants),
    as.integer(ndpost),
    as.integer(nskip),
    if (is.null(splitprobs)) NULL else as.double(splitprobs),
    dartParams
  )

  list(
    yhat.train = t(result$yhat.train),
    yhat.test = if (!is.null(result$yhat.test)) t(result$yhat.test) else NULL,
    sigma = result$sigma,
    varcount = t(result$varcount)
  )
}
