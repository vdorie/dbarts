# The negative-binomial (count) R surface (docs/design/negative-binomial.md
# sections 4-5; docs/plans/negative-binomial.md C3): explicit family dispatch
# and refusals, response and dispersion validation, fit-object shapes, the
# dispersion draws and mean-count reporting, prediction, offset (log-exposure)
# semantics, state save/load through the R surface, setResponse mutation
# semantics (r kept), and a seeded statistical recovery smoke. The
# exact-posterior gate lives in benchmarks/R/negbin-exact.R.

source(system.file("common", "stateContinuation.R", package = "dbarts"))

set.seed(99)
n <- 200L
x <- matrix(runif(n * 3L), n, 3L)
etaTrue <- 0.9 * (x[, 1L] - 0.5) + 0.6 * (x[, 2L] - 0.5)
y <- rnbinom(n, size = 5L, mu = exp(etaTrue))
x.test <- matrix(runif(10L * 3L), 10L, 3L)

n.samples <- 60L
n.burn <- 30L
n.trees <- 15L

# --- explicit family = "nbinom" on a bart2 fit: no auto, a distinct class ---

fit <- bart2(
  x,
  y,
  test = x.test,
  family = "nbinom",
  n.samples = n.samples,
  n.burn = n.burn,
  n.trees = n.trees,
  n.chains = 1L,
  verbose = FALSE,
  keepTrees = TRUE
)
expect_inherits(fit, "bartNegbin")
expect_false(inherits(fit, "bart"))
expect_equal(fit$family, "nbinom")

# a count response under family = "auto" is a plain numeric (gaussian) fit -
# counts carry no unambiguous class, so nbinom is never inferred
suppressMessages(
  fitAuto <- bart2(
    x,
    y,
    n.samples = 20L,
    n.burn = 10L,
    n.trees = n.trees,
    n.chains = 1L,
    verbose = FALSE
  )
)
expect_equal(fitAuto$family, "gaussian")

# --- channel shapes: mean counts, latent psi, and the dispersion draws ---

expect_equal(dim(fit$yhat.train), c(n.samples, n))
expect_equal(dim(fit$yhat.test), c(n.samples, 10L))
expect_equal(dim(fit$latent.train), c(n.samples, n))
expect_equal(dim(fit$latent.test), c(n.samples, 10L))
expect_equal(length(fit$dispersion), n.samples)
expect_equal(dim(fit$varcount), c(n.samples, 3L))
expect_equal(length(fit$y), n)
# mean counts are positive; dispersion draws lie on the shipped positive grid
expect_true(all(fit$yhat.train > 0))
expect_true(all(
  fit$dispersion %in%
    c(1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 30, 50)
))
# mu = r exp(psi) holds draw by draw
expect_equal(fit$yhat.train, fit$dispersion * exp(fit$latent.train))

# --- fitted ---

meanCounts <- fitted(fit)
expect_equal(length(meanCounts), n)
expect_true(all(meanCounts > 0))
expect_equal(meanCounts, colMeans(fit$yhat.train))
latentMean <- fitted(fit, type = "bart")
expect_equal(length(latentMean), n)
expect_equal(latentMean, colMeans(fit$latent.train))

# --- extract ---

expect_equal(extract(fit, type = "ev"), fit$yhat.train)
expect_equal(extract(fit, type = "bart"), fit$latent.train)
expect_equal(extract(fit, type = "bart", sample = "test"), fit$latent.test)
ppd <- extract(fit, type = "ppd")
expect_equal(dim(ppd), c(n.samples, n))
expect_true(all(ppd >= 0 & ppd == round(ppd)))

# --- residuals: observed count minus posterior-mean count ---

res <- residuals(fit)
expect_equal(length(res), n)
expect_equal(res, y - fitted(fit))

# --- predict: replay reproduces the recorded test channels (draw-neutral) ---

expect_equal(predict(fit, x.test), fit$yhat.test)
expect_equal(predict(fit, x.test, type = "bart"), fit$latent.test)
ppdNew <- predict(fit, x.test, type = "ppd")
expect_equal(dim(ppdNew), c(n.samples, 10L))
expect_true(all(ppdNew >= 0 & ppdNew == round(ppdNew)))

# a log-exposure offset.test scales the mean counts multiplicatively
evPlain <- predict(fit, x.test, type = "ev")
evDoubled <- predict(fit, x.test, type = "ev", offset.test = rep(log(2), 10L))
expect_equal(evDoubled, 2 * evPlain)

# --- predict requires keepTrees ---

fitNoTrees <- bart2(
  x,
  y,
  family = "nbinom",
  n.samples = 20L,
  n.burn = 10L,
  n.trees = n.trees,
  n.chains = 1L,
  verbose = FALSE
)
expect_error(predict(fitNoTrees, x.test), pattern = "keepTrees")

# --- keepSampler retains $fit independent of keepTrees (D2) ---
fitKeepSampler <- bart2(
  x,
  y,
  family = "nbinom",
  n.samples = 20L,
  n.burn = 10L,
  n.trees = n.trees,
  n.chains = 1L,
  verbose = FALSE,
  keepSampler = TRUE,
  keepTrees = FALSE
)
expect_false(is.null(fitKeepSampler$fit))
expect_true(is.null(fitKeepSampler$bc))
rm(fitKeepSampler)

# --- the retained $fit is a host shell, as a multinomial fit's is ---
# bc holds the count model; the host carries only the design and priors it was
# built from, so a mutation through it used to be a silent no-op and
# $getDispersion() answered with the placeholder's own r
hostRefusal <- "host sampler of a bart2\\(family = \"nbinom\"\\) fit"
expect_error(fit$fit$setResponse(as.double(y)), hostRefusal)
expect_error(fit$fit$setData(dbartsData(x, as.double(y))), hostRefusal)
expect_error(fit$fit$run(0L, 1L), hostRefusal)
expect_error(fit$fit$getDispersion(), hostRefusal)
# the read surface predict() threads through is untouched
expect_equal(ncol(fit$fit$data@x), ncol(x))
expect_equal(predict(fit, x.test), fit$yhat.test)

# --- print names the family and reports the dispersion ---

printed <- capture.output(print(fit))
expect_true(any(grepl("negative binomial", printed)))
expect_true(any(grepl("dispersion", printed)))

# --- multi-chain shapes, uncombined ---

fit2c <- bart2(
  x,
  y,
  family = "nbinom",
  n.samples = 20L,
  n.burn = 10L,
  n.trees = n.trees,
  n.chains = 2L,
  combineChains = FALSE,
  verbose = FALSE,
  keepTrees = TRUE
)
expect_equal(dim(fit2c$yhat.train), c(2L, 20L, n))
expect_equal(dim(fit2c$dispersion), c(2L, 20L))
expect_equal(
  dim(predict(fit2c, x.test, combineChains = FALSE)),
  c(2L, 20L, 10L)
)

# --- fixed dispersion: r held at the supplied integer, no grid draw ---

fitFixed <- bart2(
  x,
  y,
  family = "nbinom",
  dispersion = 4,
  n.samples = 20L,
  n.burn = 10L,
  n.trees = n.trees,
  n.chains = 1L,
  verbose = FALSE
)
expect_true(all(fitFixed$dispersion == 4))

# --- validation refusals ---

expect_error(
  dbarts(x, y - 1, family = "nbinom"),
  pattern = "non-negative integer"
)
expect_error(
  dbarts(x, y + 0.5, family = "nbinom"),
  pattern = "non-negative integer"
)
expect_error(
  dbarts(x, y, family = "nbinom", dispersion = 2.5),
  pattern = "real dispersion is not yet supported"
)
expect_error(
  dbarts(x, y, family = "nbinom", dispersion = -3),
  pattern = "positive"
)
expect_error(
  dbarts(x, y, family = "nbinom", weights = runif(n, 0.5, 2)),
  pattern = "weights"
)
# weights identically 1 are treated as absent, the probit courtesy
samplerW1 <- dbarts(x, y, family = "nbinom", weights = rep(1, n))
expect_null(samplerW1$data@weights)

# rbart_vi and xbart refuse the count family (family-vector omission)
expect_error(
  rbart_vi(
    y ~ x,
    group.by = rep(1:4, each = 50L),
    family = "nbinom",
    n.samples = 10L,
    n.burn = 5L
  )
)
expect_error(xbart(x, y, family = "nbinom", n.samples = 10L, n.reps = 1L))

# --- dbarts()-direct sampler: coding, attribute, node.scale, and state ---

control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = n.trees,
  updateState = TRUE
)
sampler <- dbarts(x, y, family = "nbinom", control = control, verbose = FALSE)
expect_equal(sampler$model@family, "nbinom")
expect_equal(attr(sampler$control, "bartcore.dispersion"), -1)
expect_equal(sampler$model@node.scale, pi * sqrt(3))

invisible(sampler$run(20L, 5L))
state1 <- sampler$state
expect_true(state1[[1L]]$dispersion > 0)
expect_true(is.finite(state1[[1L]]$dispersion))
expect_equal(length(state1[[1L]]$latents), n)

# a zero (default) offset is kept for this fixed-unit-scale family, as probit's
offsetSampler <- dbarts(
  x,
  y,
  family = "nbinom",
  offset = rep(0, n),
  control = control,
  verbose = FALSE
)
expect_false(is.null(offsetSampler$data@offset))

# --- state save/load through the R surface: restore is deterministic ---

saved <- sampler$state
samplerB <- dbarts(x, y, family = "nbinom", control = control, verbose = FALSE)
samplerC <- dbarts(x, y, family = "nbinom", control = control, verbose = FALSE)
samplerB$setState(saved)
samplerC$setState(saved)
expect_true(statesAgree(samplerB$state, saved))
expect_identical(samplerB$state[[1L]]$dispersion, saved[[1L]]$dispersion)
rB <- samplerB$run(0L, 5L)
rC <- samplerC$run(0L, 5L)
expect_identical(rB$train, rC$train)

# an nbinom sampler refuses a state lacking its dispersion/latents block
gaussSampler <- dbarts(x, as.double(y), verbose = FALSE)
invisible(gaussSampler$run(5L, 2L))
expect_error(samplerB$setState(gaussSampler$state))

# --- mutation semantics: setResponse keeps r, redraws omega ---

before <- sampler$state
rBefore <- before[[1L]]$dispersion
latentsBefore <- before[[1L]]$latents

set.seed(101)
swap <- sample.int(n, 40L)
yNew <- y
yNew[swap] <- rnbinom(40L, size = 5L, mu = 3)
# updateState = TRUE re-stores the state so the post-mutation snapshot reflects
# the swap (the $state field is otherwise a stale stored copy)
sampler$setResponse(as.double(yNew), updateState = TRUE)

after <- sampler$state
# r is a slow-moving global the count swap keeps; the omega latents are redrawn
expect_identical(after[[1L]]$dispersion, rBefore)
expect_false(identical(after[[1L]]$latents, latentsBefore))
expect_true(all(after[[1L]]$latents > 0)) # omega are Polya-Gamma (positive)

# the support rule creation applies is applied at mutation too. A negative
# element is memory safety, not taste: it used to size the dispersion kernel's
# count histogram through static_cast<size_t>(lround(y)), underflowing into a
# ~1.8e19 allocation that took the process down uncatchably. Magnitude is the
# same allocation defect from the other side (see the cap below); a non-finite
# element is refused with the rest.
countRefusal <- "family \"nbinom\" requires a non-negative integer"
expect_error(
  sampler$setResponse(replace(as.double(yNew), 1L, -1)),
  countRefusal
)
expect_error(
  sampler$setResponse(replace(as.double(yNew), 1L, 2.5)),
  countRefusal
)
expect_error(
  sampler$setResponse(replace(as.double(yNew), 1L, Inf)),
  countRefusal
)
# a refused swap leaves the installed response alone, and a valid one still
# lands afterwards
expect_identical(sampler$data@y, as.double(yNew))
sampler$setResponse(as.double(replace(yNew, 1L, yNew[1L] + 1L)))
expect_identical(sampler$data@y, as.double(replace(yNew, 1L, yNew[1L] + 1L)))

# --- the magnitude cap: 1e6, at creation and on both mutation conduits ---
# NBDispersionPrior::computeKernel allocates maxCount + 1 doubles, 8 bytes per
# unit of the largest count, so an unbounded y is the same allocation defect a
# negative one is: y = 1e9 asks for 8 GB where no R error can be raised. The
# bound pins that at 8 MB, and creation and every conduit that swaps y state it
# alike (the flat C surface's half is in test-capi.R)
capRefusal <- "counts no larger than 1000000"
overBound <- as.double(replace(y, 1L, 1000001))
expect_error(dbarts(x, overBound, family = "nbinom"), capRefusal)

capSampler <- dbarts(
  x,
  as.double(y),
  family = "nbinom",
  control = control,
  verbose = FALSE
)
expect_error(capSampler$setResponse(overBound), capRefusal)
expect_error(capSampler$setData(dbartsData(x, overBound)), capRefusal)
# a refused swap leaves the installed response alone, and both conduits still
# take an in-bound one
expect_identical(capSampler$data@y, as.double(y))
capSampler$setResponse(as.double(replace(y, 1L, y[1L] + 1L)))
capSampler$setData(dbartsData(x, as.double(y)))
expect_identical(capSampler$data@y, as.double(y))
# the bound itself is IN - the comparison is > and not >= - so a response whose
# largest count is exactly 1e6 builds, at an 8 MB histogram
atBound <- dbarts(
  x,
  as.double(replace(y, 1L, 1000000)),
  family = "nbinom",
  control = control,
  verbose = FALSE
)
expect_equal(max(atBound$data@y), 1e6)

# --- recovery: mean-count calibration and r, statistically ---

set.seed(4242)
nRec <- 400L
xRec <- matrix(runif(nRec * 2L), nRec, 2L)
fRec <- 0.8 * sin(2 * pi * xRec[, 1L]) + 0.4 * (xRec[, 2L] - 0.5)
muRec <- exp(fRec)
rTrue <- 5L
yRec <- rnbinom(nRec, size = rTrue, mu = muRec)

fitRec <- bart2(
  xRec,
  yRec,
  family = "nbinom",
  n.samples = 250L,
  n.burn = 150L,
  n.trees = 50L,
  n.chains = 1L,
  verbose = FALSE
)
muHat <- fitted(fitRec)
# mean counts track the truth, and the estimated dispersion brackets rTrue = 5
expect_true(mean(abs(muHat - muRec)) < 0.6)
expect_true(abs(mean(fitRec$dispersion) - rTrue) < 3)
