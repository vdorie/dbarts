source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# ---- shape round-trip: bart-convention arrays -> (iteration, chain,
# variable); this is dbarts:::bartDrawsArray, the mapping both draws() and
# summary() build on.

## combineChains = FALSE pinned deliberately: this fit exercises
## bartDrawsArray's reconstruction of the chain axis from an uncombined
## (n.chains x n.samples[ x n.vars]) stored shape, now that bart's own
## stored-object default is combined (see the TRUE/FALSE pair below for the
## combined-shape side of the same reconstruction)
fit <- dbarts::bart(
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
singleFit <- dbarts::bart(
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

# ---- summary() always produces a plain per-variable table

s <- summary(fit)
expect_equal(class(s), "summary.bart")
expect_true(is.data.frame(s$stats))
expect_equal(s$stats$variable, "sigma")

# no scalar parameters at all: a binary fit with a fixed (unmodeled) k
n.bin <- 60L
x.bin <- matrix(runif(n.bin * 2), n.bin, 2)
y.bin <- rbinom(n.bin, 1L, plogis(3 * x.bin[, 1] - 1.5))
noScalarFit <- dbarts::bart(
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
expect_error(
  dbarts:::bartDrawsArray(noScalarFit, c("sigma", "k", "tau")),
  "none of 'vars' \\(sigma, k, tau\\) are present on this fit"
)

# print.summary.bart's notes, checked directly against handcrafted stats so
# the assertion does not depend on any particular fit crossing 1.01 by chance
poorFit <- structure(
  list(
    call = quote(f()),
    stats = data.frame(variable = "x", rhat = 1.5)
  ),
  class = "summary.bart"
)
expect_true(any(grepl("R-hat", capture.output(print(poorFit)), fixed = TRUE)))

goodFit <- structure(
  list(
    call = quote(f()),
    stats = data.frame(variable = "x", rhat = 1.0)
  ),
  class = "summary.bart"
)
expect_false(any(grepl("R-hat", capture.output(print(goodFit)), fixed = TRUE)))

# ---- draws() is the exported accessor over the same array bartDrawsArray
# builds internally, and summary()'s rhat/ess_bulk/ess_tail columns are
# always present (no 'posterior' branch to degrade) ----

d <- draws(fit, "sigma")
expect_equal(dim(d), c(20L, 3L, 1L))
expect_equal(unclass(d), dbarts:::bartDrawsArray(fit, "sigma"))

s <- summary(fit)
expect_true(all(c("rhat", "ess_bulk", "ess_tail") %in% names(s$stats)))
rm(d)

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
  poorFit,
  goodFit
)

# summary()/as_draws_*() on a COMBINED (default) multi-chain fit must
# reconstruct a non-scalar field's chain axis from its combined
# (n.chains * n.samples) x n.vars layout, not the (n.chains, n.samples)
# layout an uncombined scalar field has - the two 2-D shapes are otherwise
# indistinguishable, and mistaking one for the other silently collapses
# every variable but the first. Ground truth: the same fit, uncombined,
# must produce an identical per-variable table.
combinedVc <- dbarts::bart(
  testData$x,
  testData$y,
  ndpost = 20L,
  nskip = 5L,
  ntree = 5L,
  nchain = 3L,
  nthread = 1L,
  combinechains = TRUE,
  verbose = FALSE,
  seed = 11L
)
uncombinedVc <- dbarts::bart(
  testData$x,
  testData$y,
  ndpost = 20L,
  nskip = 5L,
  ntree = 5L,
  nchain = 3L,
  nthread = 1L,
  combinechains = FALSE,
  verbose = FALSE,
  seed = 11L
)
sVcCombined <- summary(combinedVc, vars = "varcount")$stats
sVcUncombined <- summary(uncombinedVc, vars = "varcount")$stats
expect_equal(nrow(sVcCombined), ncol(testData$x))
expect_equal(sVcCombined, sVcUncombined)
rm(combinedVc, uncombinedVc, sVcCombined, sVcUncombined)

# scalarFields must list every scalar-per-draw field, not just sigma/k/tau:
# on a multi-chain, uncombined fit, first.sigma and resid.df are stored
# (n.chains, n) exactly as sigma is. Mistaking that shape for a
# per-variable field's (n.chains * n, n.vars) layout mis-splits the chain
# margin into one spurious variable per sample.
studentFit <- dbarts::bart(
  testData$y ~ testData$x,
  resid.dist = student(5),
  n.chains = 3L,
  n.samples = 10L,
  n.burn = 6L,
  n.thin = 1L,
  n.trees = 5L,
  n.threads = 1L,
  combineChains = FALSE,
  verbose = FALSE
)
expect_equal(dim(studentFit$first.sigma), c(3L, 6L))
expect_equal(dim(studentFit$resid.df), c(3L, 10L))

sFirstSigma <- summary(studentFit, vars = "first.sigma")$stats
expect_equal(nrow(sFirstSigma), 1L)
expect_equal(sFirstSigma$mean, mean(studentFit$first.sigma))

sResidDf <- summary(studentFit, vars = "resid.df")$stats
expect_equal(nrow(sResidDf), 1L)
expect_equal(sResidDf$mean, 5)

rm(studentFit, sFirstSigma, sResidDf)

# ---- known-limit checks on splitRhat/essBulk/essTail, dbarts's own
# rank-normalized split-Rhat and bulk/tail ESS (no 'posterior' involved
# anywhere in this file from here on) ----

set.seed(202)
n.iid <- 1000L
m.iid <- 4L
xIid <- matrix(rnorm(n.iid * m.iid), n.iid, m.iid)
expect_true(abs(dbarts:::splitRhat(xIid) - 1) < 0.01)
expect_true(dbarts:::essBulk(xIid) > 0.8 * n.iid * m.iid)
expect_true(dbarts:::essTail(xIid) > 0.8 * n.iid * m.iid)
rm(n.iid, m.iid, xIid)

set.seed(303)
n.bad <- 300L
half.bad <- n.bad %/% 2L
badChain <- function() c(rnorm(half.bad, 0), rnorm(n.bad - half.bad, 6))
xBad <- sapply(1:3, function(i) badChain())
expect_true(dbarts:::splitRhat(xBad) > 1.01)
rm(n.bad, half.bad, badChain, xBad)

# ---- pinned against posterior 1.7.0's own rhat()/ess_bulk()/ess_tail() on
# the same array, computed once by hand and recorded here as fixed literals
# - this file must run correctly with 'posterior' absent from the library,
# so nothing below calls it ----

set.seed(101)
n.pin <- 80L
m.pin <- 4L
phi.pin <- 0.5
makePinChain <- function() {
  e <- rnorm(n.pin)
  x <- numeric(n.pin)
  x[1L] <- e[1L]
  for (i in 2:n.pin) {
    x[i] <- phi.pin * x[i - 1L] + sqrt(1 - phi.pin^2) * e[i]
  }
  x
}
xPin <- sapply(seq_len(m.pin), function(j) makePinChain())

expect_equal(dbarts:::splitRhat(xPin), 1.0348352822317353, tolerance = 1e-12)
expect_equal(dbarts:::essBulk(xPin), 115.54086272289375, tolerance = 1e-12)
expect_equal(dbarts:::essTail(xPin), 176.89073381346304, tolerance = 1e-12)

arrPin <- array(xPin, c(n.pin, m.pin, 1L))
dimnames(arrPin) <- list(NULL, NULL, "v")
sPin <- dbarts:::summariseDraws(arrPin)
expect_equal(sPin$mean, -0.048772852699491616, tolerance = 1e-12)
expect_equal(sPin$median, -0.081883571354184684, tolerance = 1e-12)
expect_equal(sPin$sd, 0.93261434511023389, tolerance = 1e-12)
expect_equal(sPin$mad, 0.96465740153096868, tolerance = 1e-12)
expect_equal(sPin$q5, -1.5566449864298746, tolerance = 1e-12)
expect_equal(sPin$q95, 1.4312159689905077, tolerance = 1e-12)
expect_equal(sPin$rhat, 1.0348352822317353, tolerance = 1e-12)
expect_equal(sPin$ess_bulk, 115.54086272289375, tolerance = 1e-12)
expect_equal(sPin$ess_tail, 176.89073381346304, tolerance = 1e-12)

rm(n.pin, m.pin, phi.pin, makePinChain, xPin, arrPin, sPin)

# ---- a within-chain variance shift (same mean, larger scale in the second
# half) that the fold step exists to catch: bulk Rhat alone reads this as
# converged, and only the folded (tail) leg - fold, THEN split, THEN
# rank-normalize - reports the real non-convergence; pinned against
# posterior 1.7.0's rhat() the same way as above ----

set.seed(55)
n.var <- 400L
m.var <- 4L
makeVarShiftChain <- function() {
  c(rnorm(n.var / 2L, 0, 1), rnorm(n.var / 2L, 0, 3))
}
xVarShift <- sapply(seq_len(m.var), function(j) makeVarShiftChain())

expect_equal(
  dbarts:::splitRhat(xVarShift),
  1.1780351396473598,
  tolerance = 1e-12
)
expect_true(dbarts:::splitRhat(xVarShift) > 1.05)

rm(n.var, m.var, makeVarShiftChain, xVarShift)
