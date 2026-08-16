# The six compositions of vignette("dbarts-as-a-component"), asserted.
#
# A VIGNETTE ASSERTS NOTHING. A recipe that runs to completion and reports a
# WRONG answer - the defect class the embedding surface exists to prevent -
# builds cleanly and ships, so every recipe's claim is a cell here. The code
# below is the vignette's, compacted; where a tolerance is stated, the measured
# clean value and the value under the recipe's own mutation are recorded beside
# it, so a tolerance can be read as discriminating rather than merely wide.

n <- 200L
nBurn <- 100L
nDraws <- 100L

recipeControl <- function(seed, n.trees = 25L) {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = n.trees,
    updateState = FALSE,
    rngSeed = seed
  )
}

## 1. AUGMENTATION: a pure-R Albert-Chib probit over a gaussian sampler must
## reach the same posterior the engine's own probit family reaches. Both are
## given the SAME leaf prior, so this is an agreement test and not a prior
## comparison; they are two independent chains, so the tolerance is stated as a
## correlation and a root-mean-square difference between posterior means.
set.seed(101L)
x1 <- matrix(runif(n * 3L), n, 3L)
f1 <- 2 * (x1[, 1L] - 0.5) + sin(2 * pi * x1[, 2L])
o1 <- rep(0.25, n)
y1 <- rbinom(n, 1L, pnorm(f1 + o1))

host <- dbarts(
  x1,
  2 * y1 - 1 - o1,
  control = recipeControl(31L),
  resid.prior = fixed(1)
)
host$setCalibration(prior.scale = 2)
composed1 <- matrix(0, n, nDraws)
for (i in seq_len(nBurn + nDraws)) {
  z <- dbartsDrawLatents("probit", host$getFitsWithoutOffset(), y1, offset = o1)
  host$setResponse(dbartsWorkingResponse("probit", z, y1, offset = o1))
  host$run(0L, 1L)
  if (i > nBurn) composed1[, i - nBurn] <- host$getFitsWithoutOffset()
}
native <- dbarts(
  x1,
  y1,
  offset = o1,
  family = "probit",
  control = recipeControl(31L),
  node.prior = normal(k = 2, scale = 2)
)
nativeFit <- rowMeans(native$run(nBurn, nDraws)$train) - o1
composedFit <- rowMeans(composed1)

# measured: cor 0.9927, latent rms 0.1096, probability rms 0.0306
expect_true(cor(composedFit, nativeFit) > 0.97)
expect_true(sqrt(mean((composedFit - nativeFit)^2)) < 0.25)
expect_true(
  sqrt(mean((pnorm(composedFit + o1) - pnorm(nativeFit + o1))^2)) < 0.06
)
# the two pins the recipe rests on: no free residual scale, a stated leaf prior
expect_equal(as.numeric(host$getSigmas()), 1)
expect_equal(unname(host$getCalibration()[1L, "prior.scale"]), 2)

## 2. OFFSET BLOCK: BOTH blocks' truth. A partially linear model whose linear
## coefficient rides the offset channel must recover the coefficient AND the
## nonparametric part; a composition that recovered only one of them would be
## the failure this cell exists to see.
alphaTrue <- 1.25
tau <- 2
set.seed(11L)
x2 <- matrix(runif(n * 3L), n, 3L)
w2 <- rnorm(n)
f2 <- 3 * (x2[, 1L] - 0.5)^2 + sin(2 * pi * x2[, 2L])
y2 <- rnorm(n, f2 + alphaTrue * w2, 0.5)

alpha <- 0
sampler2 <- dbarts(x2, y2, offset = alpha * w2, control = recipeControl(8L))
keepAlpha <- rep(NA_real_, nDraws)
keepFit <- matrix(0, n, nDraws)
offsetIdentity <- rep(NA_real_, nDraws)
for (i in seq_len(nBurn + nDraws)) {
  installed <- alpha * w2
  sampler2$setOffset(installed)
  train <- as.numeric(sampler2$run(0L, 1L)$train)
  fit <- as.numeric(sampler2$getFitsWithoutOffset())
  sigma <- as.numeric(sampler2$getSigmas())
  v <- 1 / (sum(w2 * w2) / sigma^2 + 1 / tau^2)
  alpha <- rnorm(1L, v * sum(w2 * (y2 - fit)) / sigma^2, sqrt(v))
  if (i > nBurn) {
    keepAlpha[i - nBurn] <- alpha
    keepFit[, i - nBurn] <- fit
    offsetIdentity[i - nBurn] <- max(abs(train - fit - installed))
  }
}
# measured: alpha 1.2469 (sd 0.0358, 0.09 posterior sd from the truth),
# cor(f) 0.9690
expect_true(abs(mean(keepAlpha) - alphaTrue) < 3 * sd(keepAlpha))
expect_true(cor(rowMeans(keepFit), f2) > 0.93)
# the seeding contract is checkable: the offset in force during the sweep is
# the one installed before it, so the training channel is the fit plus exactly
# that vector on every kept draw
expect_true(max(offsetIdentity) < 1e-10)

## 2b. The same one-sweep step under simulation-based calibration, which is the
## end-to-end proof that the composition targets its posterior rather than
## merely recovering a truth once. The prior draw is BART's own, taken on the
## sampler itself so the chain starts AT theta0; the problem is small because
## every replication runs a chain.
nSbc <- 60L
sigma0 <- 1
tauSbc <- 1
set.seed(4L)
xSbc <- matrix(runif(nSbc * 2L), nSbc, 2L)
wSbc <- rnorm(nSbc)
sbcControl <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 10L,
  updateState = FALSE
)

newSbcSampler <- function() {
  s <- dbarts(
    xSbc,
    seq(-3, 3, length.out = nSbc),
    control = sbcControl,
    resid.prior = fixed(sigma0^2)
  )
  s$setCalibration(prior.scale = 2)
  s
}
drawSbcPrior <- function() {
  s <- newSbcSampler()
  s$sampleTreesFromPrior()
  s$sampleNodeParametersFromPrior()
  list(
    sampler = s,
    f = as.numeric(s$getFitsWithoutOffset()),
    alpha = rnorm(1L, 0, tauSbc)
  )
}
simulateSbc <- function(theta0) {
  rnorm(nSbc, theta0$f + theta0$alpha * wSbc, sigma0)
}
initSbc <- function(theta0, data) {
  theta0$sampler$setResponse(data)
  list(
    sampler = theta0$sampler,
    y = data,
    f = theta0$f,
    alpha = theta0$alpha
  )
}
# reportMean = TRUE is the defect the diagnostic exists to catch: the outer
# block's posterior MEAN where its draw belongs
sbcStep <- function(reportMean) {
  function(state) {
    state$sampler$setOffset(state$alpha * wSbc)
    state$sampler$run(0L, 1L)
    state$f <- as.numeric(state$sampler$getFitsWithoutOffset())
    v <- 1 / (sum(wSbc * wSbc) / sigma0^2 + 1 / tauSbc^2)
    m <- v * sum(wSbc * (state$y - state$f)) / sigma0^2
    state$alpha <- if (reportMean) m else rnorm(1L, m, sqrt(v))
    state
  }
}
sbcFunctionals <- function(state) {
  c(alpha = state$alpha, f.mean = mean(state$f))
}
validateRecipe2 <- function(reportMean) {
  dbartsValidateComposition(
    drawSbcPrior,
    simulateSbc,
    initSbc,
    sbcStep(reportMean),
    sbcFunctionals,
    n.replications = 40L,
    n.draws = 50L,
    n.thin = 5L,
    n.burn = 20L,
    seed = 12L
  )
}
sbcClean <- validateRecipe2(FALSE)
sbcMean <- validateRecipe2(TRUE)
# measured, band 0.2225: clean ecdf.diff 0.0863 (alpha) and 0.1309 (f.mean);
# under the mean-for-draw substitution alpha moves to 0.3824
expect_true(sbcClean$pass)
expect_identical(as.character(sbcClean$verdicts$verdict), c("PASS", "PASS"))
expect_false(sbcMean$pass)
expect_identical(sbcMean$verdicts$verdict[1L], "FLAG")
expect_true(sbcMean$verdicts$ecdf.diff[1L] > sbcMean$verdicts$band[1L])

## 3. K-FOREST: K samplers decomposing one additive fit. Two claims, and they
## are instrumented separately because only one of them moves when the rescale
## is dropped. The PRIOR claim is the rescale's own: each sampler was built on
## the whole response, so without setCalibration the sum carries sqrt(K) times
## the prior standard deviation a single sampler would. The FIT claim is that
## the composition reaches the same posterior mean a single sampler does.
kForest <- 2L
set.seed(21L)
x3 <- matrix(runif(n * 4L), n, 4L)
colnames(x3) <- paste0("x", 1:4)
f3 <- 2 * sin(pi * x3[, 1L]) + x3[, 2L] + 3 * (x3[, 3L] - 0.5)^2 - 2 * x3[, 4L]
y3 <- rnorm(n, f3, 0.5)

blocks <- list(x3[, 1:2, drop = FALSE], x3[, 3:4, drop = FALSE])
samplers <- lapply(blocks, function(xb) {
  dbarts(xb, y3, control = recipeControl(5L))
})
base <- samplers[[1L]]$getCalibration()[1L, "prior.scale"]
# MUTATION 1 lives here: drop this loop
for (s in samplers) {
  s$setCalibration(prior.scale = base / sqrt(kForest))
}

single <- dbarts(x3, y3, control = recipeControl(5L))
priorTotalSd <- function(group) {
  draws <- replicate(200L, {
    total <- 0
    for (s in group) {
      s$sampleTreesFromPrior()
      s$sampleNodeParametersFromPrior()
      total <- total + as.numeric(s$getFitsWithoutOffset())
    }
    total
  })
  mean(apply(draws, 1L, sd))
}
priorRatio <- priorTotalSd(samplers) / priorTotalSd(list(single))
# measured: 0.9706 with the rescale, 1.3726 without it (sqrt(2) = 1.414 less
# the finite-draw and cut-grid differences), so a tolerance of 0.15 clears the
# clean run by 5.1x and is exceeded by 2.5x under the mutation
expect_true(abs(priorRatio - 1) < 0.15)
expect_equal(
  unname(samplers[[1L]]$getCalibration()[1L, "prior.scale"]),
  unname(base) / sqrt(kForest)
)

fits <- rep(list(rep(0, n)), kForest)
total3 <- matrix(0, n, nDraws)
for (i in seq_len(nBurn + nDraws)) {
  for (k in seq_len(kForest)) {
    samplers[[k]]$setResponse(y3 - Reduce(`+`, fits[-k]))
    samplers[[k]]$run(0L, 1L)
    fits[[k]] <- as.numeric(samplers[[k]]$getFitsWithoutOffset())
  }
  if (i > nBurn) total3[, i - nBurn] <- Reduce(`+`, fits)
}
singleFit <- rowMeans(single$run(nBurn, nDraws)$train)
composedTotal <- rowMeans(total3)
# measured: cor 0.9868 against the single sampler, rms 0.1513. Dropping the
# rescale leaves both inside this band (0.9819, 0.1778), which is why the prior
# claim above carries the mutation and this one states the composition works
expect_true(cor(composedTotal, singleFit) > 0.97)
expect_true(sqrt(mean((composedTotal - singleFit)^2)) < 0.3)
# the response swap does NOT re-anchor the transform, which is what makes every
# sweep's partial residual readable on one scale
expect_equal(
  unname(samplers[[1L]]$getCalibration()[1L, "response.shift"]),
  unname(samplers[[2L]]$getCalibration()[1L, "response.shift"])
)

## 4. LATENT COVARIATE: the install mask read as an MH accept mask. The mask is
## non-vacuous by construction here (a deep tree prior on a small sample), which
## is the whole point: a host that ignores it disagrees with the samplers about
## the state on exactly the declined rows, silently and permanently.
set.seed(404L)
nJoint <- 80L
thetaTrue <- rnorm(nJoint)
d1 <- data.frame(theta = thetaTrue, u = runif(nJoint))
d1$y <- 1.5 * thetaTrue + rnorm(nJoint, 0, 0.3)
d2 <- data.frame(v = runif(nJoint), theta = thetaTrue)
d2$y <- rbinom(nJoint, 1L, pnorm(thetaTrue))
deepPrior <- dbartsPriors$cgm(power = 0.5)
sA <- dbarts(
  y ~ theta + u,
  d1,
  control = recipeControl(44L),
  tree.prior = deepPrior
)
sB <- dbarts(
  y ~ v + theta,
  d2,
  family = "probit",
  control = recipeControl(44L),
  tree.prior = deepPrior
)
invisible(sA$run(0L, 20L))
invisible(sB$run(0L, 20L))

jointLogLik <- function() {
  fA <- as.numeric(sA$getFitsWithoutOffset())
  fB <- as.numeric(sB$getFitsWithoutOffset())
  dnorm(d1$y, fA, as.numeric(sA$getSigmas()), log = TRUE) +
    ifelse(d2$y == 1, pnorm(fB, log.p = TRUE), pnorm(-fB, log.p = TRUE))
}
current <- sA$data@x[, "theta"]
llCurrent <- jointLogLik() + dnorm(current, 0, 1, log = TRUE)
proposal <- current + rnorm(nJoint, 0, 0.5)
u4 <- runif(nJoint)
installMask <- updatePredictorPerObservationJointly(
  list(sA, sB),
  proposal,
  column = "theta"
)
llProposal <- jointLogLik() + dnorm(proposal, 0, 1, log = TRUE)
accept <- installMask & (log(u4) < llProposal - llCurrent)
theta4 <- ifelse(accept, proposal, current)
reinstalled <- updatePredictorPerObservationJointly(
  list(sA, sB),
  theta4,
  column = "theta"
)

# measured: 77 of 80 installed, 37 accepted; the cell is vacuous if the mask is
# all TRUE, so its non-vacuity is asserted first
expect_true(sum(installMask) > 0L && sum(installMask) < nJoint)
expect_true(sum(accept) > 0L)
expect_true(all(theta4[!installMask] == current[!installMask]))
# the second call returns every row to the value the host settled on, and the
# samplers agree with the host and with each other
expect_true(all(reinstalled))
expect_equal(as.numeric(sA$data@x[, "theta"]), theta4)
expect_equal(as.numeric(sB$data@x[, "theta"]), theta4)
# the defect: dropping the mask from the conjunction moves the host's copy on
# the declined rows to values the samplers never installed
wrongTheta <- ifelse(log(u4) < llProposal - llCurrent, proposal, current)
expect_true(sum(wrongTheta != theta4) > 0L)
expect_true(all(which(wrongTheta != theta4) %in% which(!installMask)))

## 5. OUTER-OWNED SIGMA: the pin and the guard. resid.prior = fixed() suppresses
## the sampler's own draw, so what setSigma writes survives a sweep; without it
## the sampler redraws and the outer write is silently discarded.
set.seed(505L)
x5 <- matrix(runif(n * 3L), n, 3L)
f5 <- 2 * (x5[, 1L] - 0.5) + sin(2 * pi * x5[, 2L])
y5 <- rnorm(n, f5, 0.75)
sigma5 <- 1
sampler5 <- dbarts(
  x5,
  y5,
  control = recipeControl(55L),
  resid.prior = fixed(sigma5^2)
)
keepSigma <- rep(NA_real_, nDraws)
for (i in seq_len(nBurn + nDraws)) {
  sampler5$setSigma(sigma5)
  sampler5$run(0L, 1L)
  fit5 <- as.numeric(sampler5$getFitsWithoutOffset())
  sigma5 <- sqrt(sum((y5 - fit5)^2) / rchisq(1L, n))
  if (i > nBurn) keepSigma[i - nBurn] <- sigma5
}
# measured: posterior mean 0.678 against a truth of 0.75
expect_true(abs(mean(keepSigma) - 0.75) < 0.2)
# the pin: a write survives a sweep here and does not on a free-sigma sampler
sampler5$setSigma(2)
sampler5$run(0L, 1L)
expect_equal(as.numeric(sampler5$getSigmas()), 2)
freeSigma <- dbarts(x5, y5, control = recipeControl(55L))
freeSigma$setSigma(2)
freeSigma$run(0L, 1L)
expect_true(abs(as.numeric(freeSigma$getSigmas()) - 2) > 0.5)
# the guard, by message: a family with no residual scale of its own refuses the
# write rather than rescaling every leaf posterior precision for the run
probit5 <- dbarts(
  x5,
  rbinom(n, 1L, 0.5),
  family = "probit",
  control = recipeControl(55L)
)
expect_error(
  probit5$setSigma(1),
  "only gaussian and aft samplers carry a sigma"
)

## 6. INCREMENTAL USE: the slope bias. The outer block draws the coefficient of
## the covariate its own offset carries, from the latent residual. Reading
## run()$train where $getFitsWithoutOffset() belongs subtracts the offset twice,
## which reaches the outer block as a coefficient shrunk toward zero and as no
## error at all.
gammaTrue <- 1.5
set.seed(314L)
x6 <- matrix(runif(n * 3L), n, 3L)
v6 <- rnorm(n)
f6 <- 2 * (x6[, 1L] - 0.5) + sin(2 * pi * x6[, 2L])
y6 <- rbinom(n, 1L, pnorm(f6 + gammaTrue * v6))

gamma <- 0
sampler6 <- dbarts(
  x6,
  y6,
  offset = gamma * v6,
  family = "probit",
  control = recipeControl(20L)
)
keepGamma <- rep(NA_real_, nDraws)
shortfall <- rep(NA_real_, nDraws)
for (i in seq_len(nBurn + nDraws)) {
  offset6 <- gamma * v6
  sampler6$setOffset(offset6)
  train6 <- as.numeric(sampler6$run(0L, 1L)$train)
  z6 <- as.numeric(sampler6$getLatents())
  # MUTATION 2 lives here: read train6 in place of the accessor
  fit6 <- as.numeric(sampler6$getFitsWithoutOffset())
  vv <- 1 / (sum(v6 * v6) + 1 / 4)
  gamma <- rnorm(1L, vv * sum(v6 * (z6 - fit6)), sqrt(vv))
  if (i > nBurn) {
    keepGamma[i - nBurn] <- gamma
    shortfall[i - nBurn] <- max(abs((z6 - train6) - ((z6 - fit6) - offset6)))
  }
}
# the oracle standard error of a coefficient read off n latents with unit
# variance. Measured: the correct read sits 3.2 of them below the truth (Monte
# Carlo, and it walks to 2.1 at three times the draws), the run()$train read
# 16.78 below it (posterior mean 0.286 against a truth of 1.5) - a 13.6 SE
# separation across a threshold of 6
oracleSe <- 1 / sqrt(sum(v6 * v6))
expect_true(abs(mean(keepGamma) - gammaTrue) < 6 * oracleSe)
# and the arithmetic behind it: the two residuals differ by the installed
# offset exactly, which is why the wrong one is short by gamma * v
expect_true(max(shortfall) < 1e-10)
