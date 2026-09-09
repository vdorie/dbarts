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

# summary()/draws() on a COMBINED (default) multi-chain fit must
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

# every literal pinned below was computed from a set.seed() draw taken
# under Mersenne-Twister/Inversion/Rejection; pin that kind explicitly
# (restored at the end of the file) so an earlier test file's leftover
# RNGkind() cannot change what these fixtures draw - both test_package's
# run order and which RNGkind() the process happens to carry in are
# outside this file's control
oldRNGkind <- RNGkind()
suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rejection"))

# every fixture below is rounded onto a 1/64 grid by dyadic(). splitRhat's
# tail leg folds the pooled draws around their median; with an even pooled
# count S (every fixture here has one) that median is the mean of the two
# middle order statistics, and the two folded values built from those two
# order statistics tie MATHEMATICALLY, exactly - abs(x - m) is the same for
# both by construction, for any input. Whether that tie survives as an
# exact floating-point equality is a separate question: R's mean() (which
# median() calls to average the two middle values) accumulates in a long
# double on x86_64 and in a plain double on arm64, so for an arbitrary
# non-dyadic median those two accumulations can round to different bit
# patterns, breaking the tie on one architecture and not the other;
# rank()'s ties.method = "average" then assigns the pair different ranks
# on one side, and the rank-normalized folded Rhat differs at the 1e-3
# level - far past the 1e-12 a pinned literal needs. Rounding every draw to
# a multiple of 1/64 makes the sum (and therefore the mean) of any two such
# draws exactly representable in double precision with no rounding at all,
# so a long-double accumulator and a plain-double accumulator produce the
# identical bit pattern: the tie (or non-tie) becomes an architecture-
# invariant fact about the data instead of an artifact of accumulator
# width. 'posterior' folds and medians the same way, so its own
# rhat()/ess_bulk()/ess_tail() share the sensitivity and the fix.
dyadic <- function(x) round(x * 64) / 64

set.seed(202)
n.iid <- 1000L
m.iid <- 4L
xIid <- dyadic(matrix(rnorm(n.iid * m.iid), n.iid, m.iid))
expect_true(abs(dbarts:::splitRhat(xIid) - 1) < 0.01)
expect_true(dbarts:::essBulk(xIid) > 0.8 * n.iid * m.iid)
expect_true(dbarts:::essTail(xIid) > 0.8 * n.iid * m.iid)
rm(n.iid, m.iid, xIid)

set.seed(303)
n.bad <- 300L
half.bad <- n.bad %/% 2L
badChain <- function() c(rnorm(half.bad, 0), rnorm(n.bad - half.bad, 6))
xBad <- dyadic(sapply(1:3, function(i) badChain()))
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
xPin <- dyadic(sapply(seq_len(m.pin), function(j) makePinChain()))

expect_equal(dbarts:::splitRhat(xPin), 1.0348827840167951, tolerance = 1e-12)
expect_equal(dbarts:::essBulk(xPin), 115.44342922345757, tolerance = 1e-12)
expect_equal(dbarts:::essTail(xPin), 176.89073381346304, tolerance = 1e-12)

arrPin <- array(xPin, c(n.pin, m.pin, 1L))
dimnames(arrPin) <- list(NULL, NULL, "v")
sPin <- dbarts:::summariseDraws(arrPin)
expect_equal(sPin$mean, -0.048974609375000151, tolerance = 1e-12)
expect_equal(sPin$median, -0.078125, tolerance = 1e-12)
expect_equal(sPin$sd, 0.93274625074016382, tolerance = 1e-12)
expect_equal(sPin$mad, 0.97295624999999997, tolerance = 1e-12)
expect_equal(sPin$q5, -1.5640624999999999, tolerance = 1e-12)
expect_equal(sPin$q95, 1.4242187500000005, tolerance = 1e-12)
expect_equal(sPin$rhat, 1.0348827840167951, tolerance = 1e-12)
expect_equal(sPin$ess_bulk, 115.44342922345757, tolerance = 1e-12)
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
xVarShift <- dyadic(sapply(seq_len(m.var), function(j) makeVarShiftChain()))

expect_equal(
  dbarts:::splitRhat(xVarShift),
  1.1783864715882801,
  tolerance = 1e-12
)
expect_true(dbarts:::splitRhat(xVarShift) > 1.05)

rm(n.var, m.var, makeVarShiftChain, xVarShift)

# ---- Geyer's initial monotone sequence indexes rho_hat_t[1:max_t], not
# seq_len(max_t): at max_t == 0 (short chains, few draws) 1:0 still selects
# element 1 in R, giving tau_hat = 2, exactly as posterior 1.7.0's own .ess
# does - not the empty sum seq_len(0) would give, which uncaps tau_hat and
# inflates ESS by roughly log10(ess). This exact shape (n.chains = 4L,
# n.samples = 10L) is reachable through summary(bart(...)) and previously
# reported ess_bulk 64.08 against posterior's 20 ----

set.seed(1)
xShort <- dyadic(matrix(rnorm(40), 10, 4))
expect_equal(dbarts:::splitRhat(xShort), 0.95575837010045062, tolerance = 1e-12)
expect_equal(dbarts:::essBulk(xShort), 20, tolerance = 1e-12)
expect_equal(dbarts:::essTail(xShort), 20, tolerance = 1e-12)
rm(xShort)

# same maxT == 0 shape, reached through a real fit's summary() rather than
# a hand-built matrix, pinned against posterior 1.7.0's own rhat()/
# ess_bulk() on the same draws array. Not made dyadic like the fixtures
# above: with n.samples = 10 split into half-chains of n = 5 rows, the
# essFromHalfChains() while loop condition (t < nrow(acov) - 5L) is false
# before it ever runs, so maxT stays 0 and tau_hat = 2 regardless of the
# draws' values - ess_bulk == 20 is a structural consequence of the
# (n.chains, n.samples) shape, not a statistic computed FROM the data, so
# it carries none of the median-tie sensitivity above.
set.seed(1)
xShort.x <- rnorm(50)
xShort.y <- rnorm(50)
fitShort <- dbarts::bart(
  xShort.x,
  xShort.y,
  n.chains = 4L,
  n.samples = 10L,
  n.burn = 5L,
  n.trees = 5L,
  n.threads = 1L,
  verbose = FALSE
)
sShort <- summary(fitShort)
expect_equal(sShort$stats$ess_bulk[1L], 20, tolerance = 1e-12)
rm(xShort.x, xShort.y, fitShort, sShort)

# ---- strongly negative lag-1 autocorrelation (alternating-sign draws),
# pinned against posterior 1.7.0's own ess_bulk()/ess_tail()/rhat() on the
# same array - a regime the fold/rank-normalize/Geyer-sequence machinery
# above is not otherwise exercised against ----

set.seed(2)
n.alt <- 60L
m.alt <- 4L
makeAltChain <- function() {
  e <- rnorm(n.alt)
  z <- numeric(n.alt)
  z[1L] <- e[1L]
  for (i in 2:n.alt) {
    z[i] <- -0.9 * z[i - 1L] + sqrt(1 - 0.81) * e[i]
  }
  z
}
xAlt <- dyadic(sapply(seq_len(m.alt), function(j) makeAltChain()))

expect_equal(dbarts:::splitRhat(xAlt), 1.0929339523688637, tolerance = 1e-12)
expect_equal(dbarts:::essBulk(xAlt), 571.25069801078541, tolerance = 1e-12)
expect_equal(dbarts:::essTail(xAlt), 63.823464509010577, tolerance = 1e-12)

rm(n.alt, m.alt, makeAltChain, xAlt)

# ---- an NA or Inf draw must propagate to NA rather than crash or report
# finite garbage: rank()'s na.last = TRUE default assigns NA a (last) rank
# unless explicitly undone, and stats::quantile() errors outright on an
# unguarded NA/NaN input. posterior 1.7.0's summarise_draws() reports NA in
# every column for an NA draw; an Inf draw leaves most columns finite (only
# the affected tail-ESS quantile goes NA) ----

set.seed(9)
xNA <- dyadic(matrix(rnorm(80), 20, 4))
xNA[3L, 2L] <- NA
arrNA <- array(xNA, c(20L, 4L, 1L))
dimnames(arrNA) <- list(NULL, NULL, "v")
sNA <- dbarts:::summariseDraws(arrNA)
expect_true(is.na(sNA$mean))
expect_true(is.na(sNA$median))
expect_true(is.na(sNA$sd))
expect_true(is.na(sNA$mad))
expect_true(is.na(sNA$q5))
expect_true(is.na(sNA$q95))
expect_true(is.na(sNA$rhat))
expect_true(is.na(sNA$ess_bulk))
expect_true(is.na(sNA$ess_tail))
rm(xNA, arrNA, sNA)

set.seed(9)
xInf <- dyadic(matrix(rnorm(80), 20, 4))
xInf[3L, 2L] <- Inf
arrInf <- array(xInf, c(20L, 4L, 1L))
dimnames(arrInf) <- list(NULL, NULL, "v")
sInf <- dbarts:::summariseDraws(arrInf)
expect_equal(sInf$rhat, 1.0360095978146742, tolerance = 1e-12)
expect_equal(sInf$ess_bulk, 72.841168535515024, tolerance = 1e-8)
expect_true(is.na(sInf$ess_tail))
rm(xInf, arrInf, sInf)

suppressWarnings(RNGkind(oldRNGkind[1L], oldRNGkind[2L], oldRNGkind[3L]))
rm(oldRNGkind, dyadic)
