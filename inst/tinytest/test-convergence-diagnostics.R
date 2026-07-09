source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# ---- shape round-trip: bart-convention arrays -> (iteration, chain,
# variable); this is dbarts:::bartDrawsArray, the mapping both
# as_draws_array/as_draws_df and summary() build on. No posterior
# dependency, so this runs unconditionally.

## combineChains = FALSE pinned deliberately: this fit exercises
## bartDrawsArray's reconstruction of the chain axis from an uncombined
## (n.chains x n.samples[ x n.vars]) stored shape, now that bart2's own
## stored-object default is combined (see the TRUE/FALSE pair below for the
## combined-shape side of the same reconstruction)
fit <- dbarts::bart2(
  testData$y ~ testData$x,
  n.chains = 3L,
  n.samples = 20L,
  n.burn = 10L,
  n.thin = 1L,
  n.trees = 5L,
  n.threads = 1L,
  verbose = FALSE,
  combineChains = FALSE
)

arr <- dbarts:::bartDrawsArray(fit, "sigma")
expect_equal(dim(arr), c(20L, 3L, 1L))
expect_equal(as.vector(arr[,, 1L]), as.vector(t(fit$sigma)))

varArr <- dbarts:::bartDrawsArray(fit, "varcount")
expect_equal(dim(varArr), c(20L, 3L, dim(fit$varcount)[3L]))
expect_equal(unname(varArr[5L, 2L, 3L]), unname(fit$varcount[2L, 5L, 3L]))

# combineChains at fit time flattens the chain axis of scalar fields; the
# conversion must reconstruct it from the object's n.chains
combinedFit <- dbarts::bart(
  testData$x,
  testData$y,
  ndpost = 20L,
  nskip = 5L,
  ntree = 5L,
  nchain = 3L,
  nthread = 1L,
  combinechains = TRUE,
  verbose = FALSE,
  seed = 7L
)
uncombinedFit <- dbarts::bart(
  testData$x,
  testData$y,
  ndpost = 20L,
  nskip = 5L,
  ntree = 5L,
  nchain = 3L,
  nthread = 1L,
  combinechains = FALSE,
  verbose = FALSE,
  seed = 7L
)
expect_null(dim(combinedFit$sigma))
expect_equal(
  dbarts:::bartDrawsArray(combinedFit, "sigma"),
  dbarts:::bartDrawsArray(uncombinedFit, "sigma")
)

# single chain: no chain axis on the fields to begin with
singleFit <- dbarts::bart2(
  testData$y ~ testData$x,
  n.chains = 1L,
  n.samples = 15L,
  n.burn = 5L,
  n.trees = 5L,
  n.threads = 1L,
  verbose = FALSE
)
expect_null(dim(singleFit$sigma))
expect_equal(dim(dbarts:::bartDrawsArray(singleFit, "sigma")), c(15L, 1L, 1L))

# ---- summary() degrades to a plain per-variable table regardless of
# posterior's availability

s <- summary(fit)
expect_equal(class(s), "summary.bart")
expect_true(is.data.frame(s$stats))
expect_equal(s$stats$variable, "sigma")

# no scalar parameters at all: a binary fit with a fixed (unmodeled) k
n.bin <- 60L
x.bin <- matrix(runif(n.bin * 2), n.bin, 2)
y.bin <- rbinom(n.bin, 1L, plogis(3 * x.bin[, 1] - 1.5))
noScalarFit <- dbarts::bart2(
  y.bin ~ x.bin,
  k = 2.0,
  n.samples = 10L,
  n.burn = 5L,
  n.chains = 2L,
  n.thin = 1L,
  n.trees = 5L,
  n.threads = 1L,
  verbose = FALSE
)
expect_null(noScalarFit[["sigma"]])
expect_null(noScalarFit[["k"]])
s0 <- summary(noScalarFit)
expect_null(s0$stats)
expect_error(dbarts:::bartDrawsArray(noScalarFit, c("sigma", "k", "tau")))

# rbart_vi contributes tau alongside sigma
# combineChains = FALSE pinned deliberately, same reason as the bart2 fit
# above: this checks the uncombined (n.chains x n.samples) tau shape
n.g <- 5L
g <- factor(sample(n.g, length(testData$y), replace = TRUE))
y.grouped <- testData$y + rnorm(n.g, 0, 1)[g]
rfit <- dbarts::rbart_vi(
  y.grouped ~ testData$x,
  group.by = g,
  n.samples = 15L,
  n.burn = 5L,
  n.thin = 1L,
  n.chains = 2L,
  n.trees = 5L,
  n.threads = 1L,
  verbose = FALSE,
  combineChains = FALSE
)
expect_equal(dim(rfit$tau), c(2L, 15L))
rArr <- dbarts:::bartDrawsArray(rfit, c("sigma", "tau"))
expect_equal(dim(rArr), c(15L, 2L, 2L))
expect_equal(sort(dimnames(rArr)[[3L]]), c("sigma", "tau"))
sr <- summary(rfit)
expect_equal(sort(sr$stats$variable), c("sigma", "tau"))

# print.summary.bart's notes, checked directly against handcrafted stats so
# the assertion does not depend on any particular fit crossing 1.01 by chance
poorFit <- structure(
  list(
    call = quote(f()),
    stats = data.frame(variable = "x", rhat = 1.5),
    posterior = TRUE
  ),
  class = "summary.bart"
)
expect_true(any(grepl("R-hat", capture.output(print(poorFit)), fixed = TRUE)))

goodFit <- structure(
  list(
    call = quote(f()),
    stats = data.frame(variable = "x", rhat = 1.0),
    posterior = TRUE
  ),
  class = "summary.bart"
)
expect_false(any(grepl("R-hat", capture.output(print(goodFit)), fixed = TRUE)))

# ---- posterior-backed diagnostics: draws accessors, and summary()'s
# rhat/ess against posterior's own reference computation on the same array

havePosterior <- requireNamespace("posterior", quietly = TRUE)
if (havePosterior) {
  da <- posterior::as_draws_array(fit)
  expect_equal(class(da)[1L], "draws_array")
  expect_equal(dim(da), c(20L, 3L, 1L))

  df <- posterior::as_draws_df(fit)
  expect_equal(class(df)[1L], "draws_df")

  mu <- posterior::extract_variable_matrix(da, "sigma")
  s <- summary(fit)
  expect_true(s$posterior)
  expect_true("rhat" %in% names(s$stats))
  expect_equal(s$stats$rhat[s$stats$variable == "sigma"], posterior::rhat(mu))
  expect_equal(
    s$stats$ess_bulk[s$stats$variable == "sigma"],
    posterior::ess_bulk(mu)
  )
  rm(da, df, mu)
}

# ---- graceful degrade without posterior, simulated via a namespace stub
# so the test does not depend on posterior actually being uninstalled

if (havePosterior) {
  oldPosteriorAvailable <- dbarts:::posteriorAvailable
  tryCatch(
    {
      assignInNamespace("posteriorAvailable", function() FALSE, ns = "dbarts")
      sDegraded <- summary(fit)
      expect_false(sDegraded$posterior)
      expect_true(is.data.frame(sDegraded$stats))
      expect_false("rhat" %in% names(sDegraded$stats))
      printed <- capture.output(print(sDegraded))
      expect_true(any(grepl("install 'posterior'", printed, fixed = TRUE)))
    },
    finally = assignInNamespace(
      "posteriorAvailable",
      oldPosteriorAvailable,
      ns = "dbarts"
    )
  )
  expect_true(summary(fit)$posterior)
  rm(oldPosteriorAvailable, sDegraded, printed)
}

rm(
  fit,
  arr,
  varArr,
  combinedFit,
  uncombinedFit,
  singleFit,
  s,
  n.bin,
  x.bin,
  y.bin,
  noScalarFit,
  s0,
  n.g,
  g,
  y.grouped,
  rfit,
  rArr,
  sr,
  poorFit,
  goodFit,
  havePosterior
)
