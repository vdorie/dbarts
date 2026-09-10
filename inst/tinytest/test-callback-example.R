# The vignette's running-mean recipe (vignette("dbarts-as-a-component"),
# recipe 7; docs/design/per-draw-callbacks.md section 5), in the plain-C form
# inst/tinytest/capi/consumer.c carries beside the counting consumer
# test-capi.R drives. Compiled the same way test-capi.R compiles its
# consumer - R CMD SHLIB against the installed headers, skipping wherever
# that fails - because the two are the SAME source file; only the entry
# points registered differ.
#
# This is the test that the RECIPE is correct, not merely that it compiles:
# a callback-accumulated running mean of the training channel must equal the
# same seeded fit's own yhat.train.mean, with keepFits = TRUE forced so that
# channel comes back for the comparison (a callback fit defaults it FALSE).

consumerSource <- system.file(
  "tinytest",
  "capi",
  "consumer.c",
  package = "dbarts"
)
if (consumerSource == "") {
  exit_file("consumer source not installed")
}

buildDir <- tempfile("capi-mean")
dir.create(buildDir)
file.copy(consumerSource, file.path(buildDir, "consumer.c"))

includeDir <- system.file("include", package = "dbarts")
headerPath <- file.path(includeDir, "dbarts", "dbarts.h")
if (!nzchar(includeDir) || !file.exists(headerPath)) {
  msg <- paste0("dbarts.h not found under includeDir '", includeDir, "'")
  if (nzchar(Sys.getenv("CI", ""))) stop(msg) else exit_file(msg)
}

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
      "could not compile the callback-example consumer under CI:\n",
      paste(compileOutput, collapse = "\n")
    )
  }
  exit_file("could not compile the callback-example consumer")
}

dll <- dyn.load(sharedLib)
CALL <- function(name, ...) .Call(getNativeSymbolInfo(name, dll), ...)

meanFn <- CALL("capi_mean_function")

n <- 40L
nChains <- 2L
set.seed(303, sample.kind = "Rejection")
x <- matrix(runif(n * 2L), n, 2L)
y <- x[, 1L] + rnorm(n, 0, 0.2)

# the ONLY per-observation allocation, per the vignette: REAL(acc) is what
# capi_mean_context_new takes, so acc must not be reassigned before the run
acc <- numeric(n * nChains)
ctx <- CALL("capi_mean_context_new", acc, n, nChains)

fit <- bart(
  x,
  y,
  n.chains = nChains,
  n.threads = nChains,
  n.trees = 20L,
  n.burn = 20L,
  n.samples = 40L,
  keepFits = TRUE, # explicit: overrides the callback's automatic FALSE
  callback = list(fn = meanFn, context = ctx),
  verbose = FALSE
)

expect_equal(CALL("capi_mean_status"), 0L)
means <- rowMeans(matrix(acc, n, nChains))
expect_equal(means, fit$yhat.train.mean, tolerance = 1e-8)
