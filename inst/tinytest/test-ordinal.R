# The ordinal (cumulative-probit) R surface: family dispatch and
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
# a single-category response has no cutpoint to place, and is refused by count
# rather than left to fail on an empty cutpoint vector downstream
expect_error(
  bart2(
    x,
    ordered(rep_len("only", n)),
    family = "ordinal",
    n.samples = 10L,
    n.burn = 5L,
    n.trees = n.trees,
    n.chains = 1L,
    verbose = FALSE
  ),
  pattern = "at least 2 categories"
)
# a continuous response has no plausible ordered-category reading; the
# distinct sorted values would otherwise silently become one category apiece
expect_error(
  bart2(
    x,
    rnorm(n),
    family = "ordinal",
    n.samples = 10L,
    n.burn = 5L,
    n.trees = n.trees,
    n.chains = 1L,
    verbose = FALSE
  ),
  "family \"ordinal\" requires a factor or an integer-valued response",
  fixed = TRUE
)
# an integer-valued numeric response is accepted, exactly as a factor is
expect_inherits(
  bart2(
    x,
    codes,
    family = "ordinal",
    n.samples = 10L,
    n.burn = 5L,
    n.trees = n.trees,
    n.chains = 1L,
    verbose = FALSE
  ),
  "bartOrdinal"
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
# and it MOVES: gamma_2 is redrawn every sweep, where a sampler that froze it
# at its cold start (spacing 1) would sit inside the band above with sd 0
expect_true(sd(fitRec$cutpoints[, 2L]) > 0.02)

# --- type synonyms: "response" and "link" are the predict.glm spellings of
# "ev" and "bart", accepted here exactly as on a "bart" fit ---

expect_identical(extract(fit, type = "response"), extract(fit, type = "ev"))
expect_identical(extract(fit, type = "link"), extract(fit, type = "bart"))
expect_identical(fitted(fit, type = "link"), fitted(fit, type = "bart"))
expect_identical(
  predict(fit, x.test, type = "response"),
  predict(fit, x.test, type = "ev")
)
expect_identical(
  predict(fit, x.test, type = "link"),
  predict(fit, x.test, type = "bart")
)
# a value this family does not offer refuses against the set it does
expect_error(predict(fit, x.test, type = "forest"), "type must be in 'ev'")

# --- extract's combineChains formal, honoured instead of returning the
# fit's own stored layout ---

combinedEv <- extract(fit2c, type = "ev")
expect_equal(dim(combinedEv), c(40L, n, 3L))
splitBart <- extract(fit2c, type = "bart", combineChains = FALSE)
expect_equal(dim(splitBart), c(2L, 20L, n))

# --- extract(type = "loglik"): log of the stored category probability at
# the observed level, an independently-coded oracle; extract-only,
# sample = "test" refused, shape drops the K margin ---

ll <- extract(fit, type = "loglik")
ev <- extract(fit, type = "ev")
k <- match(y, lv)
oracleLl <- log(vapply(seq_len(n), function(i) ev[1L, i, k[i]], numeric(1L)))
expect_equal(ll[1L, ], oracleLl, tolerance = 1e-12)
expect_equal(dim(ll), dim(ev)[-3L])
expect_error(
  extract(fit, type = "loglik", sample = "test"),
  "no test response exists"
)
expect_error(predict(fit, x.test, type = "loglik"), "type must be in")
expect_error(fitted(fit, type = "loglik"), "type must be in")

# --- fitted/predict's ci.level: (obs x K x 3) on type = "ev", a plain
# 3-column matrix on type = "bart" ---

fitCiEv <- fitted(fit, ci.level = 0.9)
expect_equal(dim(fitCiEv), c(n, 3L, 3L))
fitCiBart <- fitted(fit, type = "bart", ci.level = 0.9)
expect_equal(dim(fitCiBart), c(n, 3L))
expect_equal(colnames(fitCiBart), c("est", "ci.lower", "ci.upper"))
predCiEv <- predict(fit, x.test, ci.level = 0.9)
expect_equal(dim(predCiEv), c(10L, 3L, 3L))

# --- bart-family-only formals refused by name ---

expect_error(
  extract(fit, forest = 1L),
  "'forest' is not used by extract on a bartOrdinal fit",
  fixed = TRUE
)
expect_error(
  extract(fit, contribution = TRUE),
  "'contribution' is not used by extract on a bartOrdinal fit",
  fixed = TRUE
)
expect_error(
  predict(fit, x.test, forest = 1L),
  "'forest' is not used by predict on a bartOrdinal fit",
  fixed = TRUE
)
expect_error(
  fitted(fit, sample = "test"),
  "'sample' is not used by fitted on a bartOrdinal fit",
  fixed = TRUE
)

# --- own-class extract's 'sample' validation now shares bart/rbart's
# wording instead of a bare match.arg's "'arg' should be one of ..." ---

expect_error(extract(fit, sample = "bogus"), "sample must be in 'train'")

# --- plotTree/survivalProbabilities refused by name ---

expect_error(
  plotTree(fit),
  "plotTree is defined for bart, rbart_vi and dbartsSampler fits",
  fixed = TRUE
)
expect_error(
  survivalProbabilities(fit),
  "survivalProbabilities applies to a discrete-time hazard fit",
  fixed = TRUE
)

# --- plot(object): two panels (K > 2) with the cutpoint trace, one panel
# at K = 2 (no free cutpoint, so par is untouched); the K > 2 branch
# restores the caller's par afterward - assert restoration against a
# sentinel value distinct from both the plot's own layout and the device
# default ---

pdf(NULL)
par(mfrow = c(3L, 3L))
plot(fit)
mfrowK3 <- par("mfrow")

y2 <- ordered(ifelse(z > 0, "hi", "lo"), levels = c("lo", "hi"))
fit2 <- bart2(
  x,
  y2,
  family = "ordinal",
  n.samples = 20L,
  n.burn = 10L,
  n.trees = n.trees,
  n.chains = 1L,
  verbose = FALSE
)
dev.off()
pdf(NULL)
plot(fit2)
mfrowK2 <- par("mfrow")
dev.off()
expect_equal(mfrowK3, c(3L, 3L))
expect_equal(mfrowK2, c(1L, 1L))

# --- as_draws_array/df default to vars = c("cutpoints", "sigma", "k",
# "tau"); cutpoint[1] (pinned at 0) is kept ---

if (requireNamespace("posterior", quietly = TRUE)) {
  ad <- posterior::as_draws_array(fit)
  adNames <- dimnames(unclass(ad))[[3L]]
  expect_true(all(c("cutpoint[1]", "cutpoint[2]") %in% adNames))
  rm(ad, adNames)
}

rm(
  combinedEv,
  splitBart,
  ll,
  ev,
  k,
  oracleLl,
  fitCiEv,
  fitCiBart,
  predCiEv,
  mfrowK3,
  y2,
  fit2,
  mfrowK2
)
