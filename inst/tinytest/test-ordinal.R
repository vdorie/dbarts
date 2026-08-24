# The ordinal (cumulative-probit) R surface (docs/design/ordinal.md sections
# 5-6): family dispatch and
# refusals, fit-object shapes and the level round-trip, prediction, state
# save/load through the R surface, setResponse mutation semantics, and a
# seeded statistical recovery smoke. The exact-posterior gate lives in
# benchmarks/R/ordinal-exact.R.

source(
  system.file("common", "stateContinuation.R", package = "dbarts"),
  local = TRUE
)

set.seed(99)
n <- 200L
x <- matrix(runif(n * 3L), n, 3L)
etaTrue <- 2 * (x[, 1L] - 0.5)
z <- etaTrue + rnorm(n)
gammaTrue <- c(0, 0.8)
codes <- 1L + (z > gammaTrue[1L]) + (z > gammaTrue[2L])
lv <- c("lo", "mid", "hi")
y <- ordered(lv[codes], levels = lv)
x.test <- matrix(runif(10L * 3L), 10L, 3L)

n.samples <- 60L
n.burn <- 30L
n.trees <- 15L

# --- auto-dispatch on an ordered factor, with announcement ---

expect_message(
  fit <- bart2(
    x,
    y,
    test = x.test,
    n.samples = n.samples,
    n.burn = n.burn,
    n.trees = n.trees,
    n.chains = 1L,
    verbose = FALSE,
    keepTrees = TRUE
  ),
  pattern = "ordinal"
)
expect_inherits(fit, "bartOrdinal")
expect_false(inherits(fit, "bart"))
expect_equal(fit$family, "ordinal")
expect_equal(fit$K, 3L)

# --- level round-trip and channel shapes ---

expect_equal(fit$levels, lv)
expect_equal(dim(fit$yhat.train), c(n.samples, n, 3L))
expect_equal(dimnames(fit$yhat.train)[[3L]], lv)
expect_equal(dim(fit$yhat.test), c(n.samples, 10L, 3L))
expect_equal(dim(fit$cutpoints), c(n.samples, 2L))
expect_equal(dim(fit$latent.train), c(n.samples, n))
expect_equal(dim(fit$varcount), c(n.samples, 3L))
expect_true(is.ordered(fit$y))
expect_equal(levels(fit$y), lv)
# probabilities: in [0, 1] and summing to 1 across the K margin
expect_true(all(fit$yhat.train >= 0 & fit$yhat.train <= 1))
expect_true(
  max(abs(apply(fit$yhat.train, c(1L, 2L), sum) - 1)) < 1e-12
)
# cutpoints: gamma_1 pinned at 0, gamma_2 above it, every draw
expect_true(all(fit$cutpoints[, 1L] == 0))
expect_true(all(fit$cutpoints[, 2L] > 0))

# --- fitted ---

phat <- fitted(fit)
expect_equal(dim(phat), c(n, 3L))
expect_equal(colnames(phat), lv)
expect_equal(unname(rowSums(phat)), rep(1, n))
cl <- fitted(fit, type = "class")
expect_true(is.ordered(cl))
expect_equal(levels(cl), lv)
expect_equal(length(fitted(fit, type = "bart")), n)

# --- extract ---

expect_equal(extract(fit, type = "ev"), fit$yhat.train)
expect_equal(extract(fit, type = "bart"), fit$latent.train)
expect_equal(extract(fit, type = "bart", sample = "test"), fit$latent.test)
ppd <- extract(fit, type = "ppd")
expect_equal(dim(ppd), c(n.samples, n))
expect_true(all(ppd %in% 1:3))

# --- residuals ---

res <- residuals(fit)
expect_equal(dim(res), c(n, 3L))
expect_equal(unname(rowSums(res)), rep(0, n), tolerance = 1e-12)

# --- predict: replay matches the recorded test channel ---

expect_equal(predict(fit, x.test), fit$yhat.test)
expect_equal(predict(fit, x.test, type = "bart"), fit$latent.test)
ppdNew <- predict(fit, x.test, type = "ppd")
expect_equal(dim(ppdNew), c(n.samples, 10L))
expect_true(all(ppdNew %in% 1:3))

# --- print names the family and the ordered levels ---

printed <- capture.output(print(fit))
expect_true(any(grepl("ordinal", printed)))
expect_true(any(grepl("lo < mid < hi", printed)))

# --- predict requires keepTrees ---

suppressMessages(
  fitNoTrees <- bart2(
    x,
    y,
    n.samples = 20L,
    n.burn = 10L,
    n.trees = n.trees,
    n.chains = 1L,
    verbose = FALSE
  )
)
expect_error(predict(fitNoTrees, x.test), pattern = "keepTrees")

# --- keepSampler retains $fit independent of keepTrees ---
suppressMessages(
  fitKeepSampler <- bart2(
    x,
    y,
    n.samples = 20L,
    n.burn = 10L,
    n.trees = n.trees,
    n.chains = 1L,
    verbose = FALSE,
    keepSampler = TRUE,
    keepTrees = FALSE
  )
)
expect_false(is.null(fitKeepSampler$fit))
rm(fitKeepSampler)

# --- the retained $fit is the engine that ran, adopted from the abandoned
# first-created host: reads and mutations succeed ---
expect_equal(ncol(fit$fit$data@x), ncol(x))
expect_equal(predict(fit, x.test), fit$yhat.test)
expect_silent(fit$fit$setResponse(as.double(codes)))
expect_identical(fit$fit$data@y, as.double(codes))
expect_silent(fit$fit$setData(dbartsData(x, as.double(codes))))
expect_silent(invisible(fit$fit$run(0L, 1L)))

# --- multi-chain shapes, uncombined ---

suppressMessages(
  fit2c <- bart2(
    x,
    y,
    n.samples = 20L,
    n.burn = 10L,
    n.trees = n.trees,
    n.chains = 2L,
    combineChains = FALSE,
    verbose = FALSE,
    keepTrees = TRUE
  )
)
expect_equal(dim(fit2c$yhat.train), c(2L, 20L, n, 3L))
expect_equal(dim(fit2c$cutpoints), c(2L, 20L, 2L))
expect_equal(
  dim(predict(fit2c, x.test, combineChains = FALSE)),
  c(2L, 20L, 10L, 3L)
)

# --- explicit family = "ordinal" on an unordered factor: informational note ---

yUnordered <- factor(lv[codes], levels = lv)
expect_message(
  fitUn <- bart2(
    x,
    yUnordered,
    family = "ordinal",
    n.samples = 20L,
    n.burn = 10L,
    n.trees = n.trees,
    n.chains = 1L,
    verbose = FALSE
  ),
  pattern = "unordered"
)
expect_equal(fitUn$levels, lv)

# --- explicit family on a numeric response: sort(unique(y)) levels ---

fitNum <- bart2(
  x,
  codes,
  family = "ordinal",
  n.samples = 20L,
  n.burn = 10L,
  n.trees = n.trees,
  n.chains = 1L,
  verbose = FALSE
)
expect_equal(fitNum$K, 3L)
expect_equal(fitNum$levels, c("1", "2", "3"))

# --- a 2-level ordered factor stays binary (probit) under auto ---

y2 <- ordered(c("lo", "hi")[1L + (z > 0)], levels = c("lo", "hi"))
suppressMessages(
  fitBin <- bart2(
    x,
    y2,
    n.samples = 20L,
    n.burn = 10L,
    n.trees = n.trees,
    n.chains = 1L,
    verbose = FALSE
  )
)
expect_inherits(fitBin, "bart")
expect_equal(fitBin$family, "probit")

# --- refusals ---

expect_error(
  bart(x, y, ndpost = 10L, nskip = 5L, verbose = FALSE),
  pattern = "bart2"
)
expect_error(
  rbart_vi(
    y ~ x,
    group.by = rep(1:4, each = 50L),
    n.samples = 10L,
    n.burn = 5L
  ),
  pattern = "ordinal"
)
expect_error(
  xbart(x, y, n.samples = 10L, n.reps = 1L),
  pattern = "ordinal"
)
expect_error(
  dbarts(x, y, family = "ordinal", weights = runif(n, 0.5, 2)),
  pattern = "weights"
)
# weights identically 1 are treated as absent, the probit courtesy
samplerW1 <- dbarts(x, y, family = "ordinal", weights = rep(1, n))
expect_null(samplerW1$data@weights)

# --- dbarts()-direct sampler: coding, attribute, and state ---

control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = n.trees,
  updateState = TRUE
)
sampler <- dbarts(x, y, family = "ordinal", control = control, verbose = FALSE)
expect_equal(sampler$model@family, "ordinal")
expect_equal(attr(sampler$control, "bartcore.n.categories"), 3L)
expect_equal(range(sampler$data@y), c(1, 3))
expect_equal(sampler$data@response.levels, lv)
expect_equal(sampler$model@node.scale, 3.0)

r <- sampler$run(20L, 5L)
expect_equal(dim(r$train), c(n, 5L))

state1 <- sampler$state
expect_equal(length(state1[[1L]]$cutpoints), 2L)
expect_equal(state1[[1L]]$cutpoints[1L], 0)
expect_true(all(diff(state1[[1L]]$cutpoints) > 0))

# --- state save/load through the R surface: restore is deterministic ---

saved <- sampler$state
samplerB <- dbarts(x, y, family = "ordinal", control = control, verbose = FALSE)
samplerC <- dbarts(x, y, family = "ordinal", control = control, verbose = FALSE)
samplerB$setState(saved)
samplerC$setState(saved)
statesAgree(samplerB$state, saved)
expect_identical(samplerB$state[[1L]]$cutpoints, saved[[1L]]$cutpoints)
rB <- samplerB$run(0L, 5L)
rC <- samplerC$run(0L, 5L)
expect_identical(rB$train, rC$train)

# an ordinal sampler refuses a state lacking its cutpoint block
probitSampler <- dbarts(x, as.double(codes == 3L), verbose = FALSE)
invisible(probitSampler$run(5L, 2L))
expect_error(
  samplerB$setState(probitSampler$state),
  "'state' length must equal number of chains"
)

# --- mutation semantics: setResponse keeps cutpoints, redraws z ---

before <- sampler$state
cutBefore <- before[[1L]]$cutpoints
latentsBefore <- before[[1L]]$latents

set.seed(101)
swap <- sample.int(n, 40L)
codesNew <- codes
codesNew[swap] <- sample.int(3L, 40L, replace = TRUE)
# updateState = TRUE re-stores the state so the post-mutation snapshot below
# reflects the swap (the $state field is otherwise a stale stored copy)
sampler$setResponse(as.double(codesNew), updateState = TRUE)

after <- sampler$state
expect_identical(after[[1L]]$cutpoints, cutBefore)
expect_false(identical(after[[1L]]$latents, latentsBefore))
# the redrawn z respect the kept cutpoints' intervals for the NEW response:
# y = k implies gamma_{k-1} < z <= gamma_k (gamma_0 = -Inf, gamma_K = +Inf)
bounds <- c(-Inf, cutBefore, Inf)
zNew <- after[[1L]]$latents
expect_true(all(zNew > bounds[codesNew] & zNew <= bounds[codesNew + 1L]))

# the category coding is fixed at creation, so mutation accepts exactly what
# creation did: an out-of-support code silently drew latents against intervals
# it had no category for
ordinalRefusal <- "ordinal response must be an integer category index"
expect_error(sampler$setResponse(as.double(codesNew) - 1), ordinalRefusal)
expect_error(sampler$setResponse(rep(4, n)), ordinalRefusal)
expect_error(sampler$setResponse(as.double(codesNew) + 0.5), ordinalRefusal)
# a refused swap leaves the installed response alone, and a valid one still
# lands afterwards
expect_identical(sampler$data@y, as.double(codesNew))
sampler$setResponse(as.double(codes))
expect_identical(sampler$data@y, as.double(codes))

# --- recovery: probabilities and the free cutpoint, statistically ---

set.seed(4242)
nRec <- 300L
xRec <- matrix(runif(nRec * 2L), nRec, 2L)
fRec <- sin(2 * pi * xRec[, 1L]) * 0.9 + 0.5 * (xRec[, 2L] - 0.5)
zRec <- fRec + rnorm(nRec)
codesRec <- 1L + (zRec > 0) + (zRec > 0.8)
yRec <- ordered(lv[codesRec], levels = lv)
trueProbs <- cbind(
  pnorm(0 - fRec),
  pnorm(0.8 - fRec) - pnorm(0 - fRec),
  1 - pnorm(0.8 - fRec)
)

fitRec <- bart2(
  xRec,
  yRec,
  family = "ordinal",
  n.samples = 250L,
  n.burn = 150L,
  n.trees = 50L,
  n.chains = 1L,
  verbose = FALSE
)
phatRec <- fitted(fitRec)
expect_true(mean(abs(phatRec - trueProbs)) < 0.08)
expect_true(abs(mean(fitRec$cutpoints[, 2L]) - 0.8) < 0.35)
