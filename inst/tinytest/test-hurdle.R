# The hurdle.lognormal combine + retransformation generics
# (docs/design/hurdle.md sections 6, 13). The
# load-bearing gate is the ANALYTIC ORACLE: a bartHurdle stubbed from hand-set
# component draws whose combined natural-scale mean, prob/link channels, and a
# seeded bimodal ppd draw are checked against their closed forms to machine
# tolerance - the only genuinely new code. Surrounding gates: predict on new
# data (shape, finiteness, replay-equals-recorded), save/load round-trip, and a
# statistical recovery smoke.

# --- (1) analytic oracle: hand-set pi / f / sigma, closed-form combine --------

S <- 6L # posterior draws
m <- 4L # observations
set.seed(1L)
piMat <- matrix(runif(S * m, 0.1, 0.9), S, m) # occupancy probabilities
fMat <- matrix(rnorm(S * m, 1, 0.5), S, m) # positive-part log-scale means
sigmaDraws <- runif(S, 0.3, 0.8) # homoscedastic per-draw sigma

# a bart fit whose ev channel is exactly piMat (probit: ev = pnorm(latent))
occ <- structure(
  list(
    family = "probit",
    yhat.train = qnorm(piMat),
    n.chains = 1L,
    call = call("NULL")
  ),
  class = "bart"
)
# a gaussian fit whose full-n test channel (the log-scale fits) is fMat
pos <- structure(
  list(
    family = "gaussian",
    yhat.test = fMat,
    sigma = sigmaDraws,
    n.chains = 1L,
    call = call("NULL")
  ),
  class = "bart"
)
h <- structure(
  list(
    call = call("NULL"),
    family = "hurdle.lognormal",
    y = numeric(m),
    occupancy = occ,
    positive = pos
  ),
  class = "bartHurdle"
)

# E[y | y > 0, x]_s = exp(f_s + sigma_s^2 / 2); E[y | x]_s = pi_s * that, per
# draw (sigma_s^2 recycles down the draw margin of the S x m matrices)
expectedEv <- piMat * exp(fMat + 0.5 * sigmaDraws^2)

expect_equal(extract(h, type = "ev"), expectedEv)
expect_equal(extract(h, type = "response"), expectedEv) # alias
expect_equal(extract(h, type = "prob"), piMat)
expect_equal(extract(h, type = "bart"), fMat)
expect_equal(extract(h, type = "link"), fMat) # alias
expect_equal(extract(h, type = "log"), fMat) # alias

# fitted is the per-observation posterior mean over draws
expect_equal(fitted(h, type = "ev"), apply(expectedEv, 2L, mean))
expect_equal(fitted(h, type = "prob"), apply(piMat, 2L, mean))
expect_equal(fitted(h, type = "bart"), apply(fMat, 2L, mean))

# ci.level opts into est + credible band, est == the posterior mean
band <- fitted(h, type = "ev", ci.level = 0.9)
expect_equal(dim(band), c(m, 3L))
expect_equal(band[, "est"], apply(expectedEv, 2L, mean))

# residuals are natural-scale y - E[y | x]
h$y <- runif(m, 0, 3)
expect_equal(residuals(h), h$y - apply(expectedEv, 2L, mean))

# residuals is always against the training response; a caller-supplied
# 'sample' collided with the fixed sample = "train" this forwards to
# fitted.bartHurdle and raised a raw 'formal argument "sample" matched by
# multiple actual arguments' - refused by name instead
expect_error(
  residuals(h, sample = "train"),
  "'sample' is not used by residuals; residuals are always against the training response",
  fixed = TRUE
)

# ppd: per draw a Bernoulli(pi) spike at zero, else a lognormal - checked by
# replaying the exact RNG consumption order (rbinom over all draws, then rnorm)
seedVal <- 20250720L
n.total <- S * m
set.seed(seedVal)
ppd <- extract(h, type = "ppd")
set.seed(seedVal)
occDraw <- rbinom(n.total, 1L, as.vector(piMat))
z <- rnorm(n.total)
amount <- exp(as.vector(fMat) + rep_len(sigmaDraws, n.total) * z)
expect_equal(ppd, array(occDraw * amount, c(S, m)))
expect_true(all(ppd >= 0))
expect_true(any(ppd == 0)) # some unoccupied draws land exactly on the spike

# type validation and the absent test channel
expect_error(extract(h, type = "bogus"), pattern = "type must be")
expect_error(extract(h, sample = "test"), pattern = "test channel")

# extract(type = "loglik"): the composed natural-scale density from the same
# closed-form pi/f/sigma used above (NOT the sum of the two components' own
# loglik channels, no truncation - the positive part's lognormal support is
# already (0, Inf)), against an independently coded oracle covering both the
# y == 0 and y > 0 branches
hLL <- h
hLL$y <- c(0, 2.5, 0, 1.1)
llH <- extract(hLL, type = "loglik")
oracleLLH <- matrix(0, S, m)
for (s in seq_len(S)) {
  for (i in seq_len(m)) {
    oracleLLH[s, i] <- if (hLL$y[i] == 0) {
      log1p(-piMat[s, i])
    } else {
      log(piMat[s, i]) +
        dnorm(log(hLL$y[i]), fMat[s, i], sigmaDraws[s], log = TRUE) -
        log(hLL$y[i])
    }
  }
}
expect_equal(llH, oracleLLH, tolerance = 1e-12)
expect_equal(dim(llH), c(S, m))
rm(hLL, llH, oracleLLH, s, i)

# --- (2) predict on new data + replay reproduces the in-sample fitted ---------

set.seed(42L)
n <- 120L
p <- 3L
x <- matrix(runif(n * p), n, p)
pi.true <- pnorm(1.0 * x[, 1L] - 0.5)
mu.true <- 1.0 + 0.8 * x[, 2L]
occupied <- rbinom(n, 1L, pi.true) == 1L
y <- numeric(n)
y[occupied] <- exp(rnorm(sum(occupied), mu.true[occupied], 0.4))

fit <- bart2(
  x,
  y,
  family = "hurdle.lognormal",
  n.samples = 50L,
  n.burn = 30L,
  n.trees = 20L,
  n.chains = 1L,
  keepTrees = TRUE,
  verbose = FALSE,
  seed = 7L
)
expect_inherits(fit, "bartHurdle")

fittedTrain <- fitted(fit, type = "ev")
expect_equal(length(fittedTrain), n)
expect_true(all(is.finite(fittedTrain) & fittedTrain >= 0))

x.new <- matrix(runif(15L * p), 15L, p)
pr <- predict(fit, x.new, type = "ev")
expect_equal(dim(pr), c(50L, 15L)) # draws x n.new
expect_true(all(is.finite(pr) & pr >= 0))

prob <- predict(fit, x.new, type = "prob")
expect_true(all(prob >= 0 & prob <= 1))
link <- predict(fit, x.new, type = "bart")
expect_equal(dim(link), c(50L, 15L))
expect_true(all(is.finite(link)))

ppdNew <- predict(fit, x.new, type = "ppd")
expect_equal(dim(ppdNew), c(50L, 15L))
expect_true(all(ppdNew >= 0))

# replaying both forests at the training x reproduces the in-sample fitted
prTrain <- predict(fit, x, type = "ev")
expect_equal(apply(prTrain, 2L, mean), fittedTrain, tolerance = 1e-6)

# predict requires keepTrees
fitNoTrees <- bart2(
  x,
  y,
  family = "hurdle.lognormal",
  n.samples = 20L,
  n.burn = 10L,
  n.trees = 10L,
  n.chains = 1L,
  verbose = FALSE,
  seed = 7L
)
expect_error(predict(fitNoTrees, x.new), pattern = "keepTrees")

# samplerOnly stays refused: its $fit is a PAIR of samplers, a composition
# checkFamilyUnsupportedArgs's per-caller allow.samplerOnly flag does not
# unblock
expect_error(
  bart2(x, y, family = "hurdle.lognormal", samplerOnly = TRUE),
  pattern = "does not support 'samplerOnly'"
)

# --- keepSampler retains $fit on both component fits, independent of
# keepTrees ---
fitKeepSampler <- bart2(
  x,
  y,
  family = "hurdle.lognormal",
  n.samples = 20L,
  n.burn = 10L,
  n.trees = 10L,
  n.chains = 1L,
  verbose = FALSE,
  seed = 7L,
  keepSampler = TRUE,
  keepTrees = FALSE
)
expect_false(is.null(fitKeepSampler$occupancy$fit))
expect_false(is.null(fitKeepSampler$positive$fit))
rm(fitKeepSampler)

# multi-chain draw alignment (chain-interleaved sigma vs chain-blocked fits) and
# combineChains reshaping
fit2c <- bart2(
  x,
  y,
  family = "hurdle.lognormal",
  n.samples = 20L,
  n.burn = 10L,
  n.trees = 10L,
  n.chains = 2L,
  keepTrees = TRUE,
  verbose = FALSE,
  seed = 5L
)
expect_equal(
  dim(extract(fit2c, type = "ev", combineChains = FALSE)),
  c(2L, 20L, n)
)
expect_equal(dim(extract(fit2c, type = "ev", combineChains = TRUE)), c(40L, n))
expect_equal(dim(predict(fit2c, x.new, combineChains = FALSE)), c(2L, 20L, 15L))
llFit2c <- extract(fit2c, type = "loglik", combineChains = FALSE)
expect_equal(dim(llFit2c), c(2L, 20L, n))

# loglik is extract-only
expect_error(predict(fit, x.new, type = "loglik"), "type must be in")
expect_error(fitted(fit, type = "loglik"), "type must be in")

# bart-family-only formals refused by name (a hurdle fit's two components
# are each a single forest)
expect_error(
  extract(fit, forest = 1L),
  "'forest' is not used by extract on a bartHurdle fit",
  fixed = TRUE
)
expect_error(
  extract(fit, contribution = TRUE),
  "'contribution' is not used by extract on a bartHurdle fit",
  fixed = TRUE
)
expect_error(
  predict(fit, x.new, forest = 1L),
  "'forest' is not used by predict on a bartHurdle fit",
  fixed = TRUE
)

# plotTree/survivalProbabilities refused by name, naming the components
expect_error(
  plotTree(fit),
  "object$occupancy$fit",
  fixed = TRUE
)
expect_error(
  survivalProbabilities(fit),
  "survivalProbabilities applies to a discrete-time hazard fit",
  fixed = TRUE
)

# plot(object): 2x2 layout, all four panels render
pdf(NULL)
plot(fit)
usedMfrow <- par("mfrow")
dev.off()
expect_equal(usedMfrow, c(2L, 2L))

# as_draws_array/df: the union of both components' present scalar fields,
# dot-prefixed by component
if (requireNamespace("posterior", quietly = TRUE)) {
  ad <- posterior::as_draws_array(fit)
  adNames <- dimnames(unclass(ad))[[3L]]
  expect_true(any(startsWith(adNames, "occupancy.")))
  expect_true(any(startsWith(adNames, "positive.")))
  expect_true("positive.sigma" %in% adNames)
  rm(ad, adNames)
}

rm(llFit2c, usedMfrow)

# --- (3) save / load round-trip: fitted and predict are identical -------------

# store each component sampler's C++ state so its pointer survives the round
# trip (the general dbarts serialization requirement, ?`dbartsSampler-class`)
serialized <- tempfile(fileext = ".rds")
invisible(fit$occupancy$fit$state)
invisible(fit$positive$fit$state)
saveRDS(fit, serialized)
fitLoaded <- readRDS(serialized)
expect_equal(fitted(fitLoaded, type = "ev"), fitted(fit, type = "ev"))
expect_equal(
  predict(fitLoaded, x.new, type = "ev"),
  predict(fit, x.new, type = "ev")
)
unlink(serialized)

# --- (4) recovery smoke: pi(x), E[y | y > 0, x], and E[y | x] track truth -----

set.seed(2718L)
nR <- 500L
xR <- matrix(runif(nR * 2L), nR, 2L)
pi.trueR <- pnorm(1.5 * xR[, 1L] - 1.0 * xR[, 2L] - 0.2)
mu.trueR <- 0.5 + 1.2 * xR[, 1L] - 0.7 * xR[, 2L]
sigma.trueR <- 0.5
occupiedR <- rbinom(nR, 1L, pi.trueR) == 1L
yR <- numeric(nR)
yR[occupiedR] <- exp(rnorm(sum(occupiedR), mu.trueR[occupiedR], sigma.trueR))

fitR <- bart2(
  xR,
  yR,
  family = "hurdle.lognormal",
  n.samples = 120L,
  n.burn = 80L,
  n.trees = 40L,
  n.chains = 1L,
  keepTrees = TRUE,
  verbose = FALSE,
  seed = 11L
)

# recovered occupancy probability tracks truth
expect_true(cor(fitted(fitR, type = "prob"), pi.trueR) > 0.6)

# recovered E[y | y > 0, x] (retransformed positive part) tracks truth
fpos <- extract(fitR$positive, type = "bart", sample = "test")
eposHat <- apply(exp(fpos + 0.5 * fitR$positive$sigma^2), 2L, mean)
expect_true(cor(eposHat, exp(mu.trueR + 0.5 * sigma.trueR^2)) > 0.6)

# combined mean tracks the true combined mean
evTrue <- pi.trueR * exp(mu.trueR + 0.5 * sigma.trueR^2)
expect_true(cor(fitted(fitR, type = "ev"), evTrue) > 0.6)

printed <- capture.output(print(fitR))
expect_true(any(grepl("hurdle.lognormal", printed)))
