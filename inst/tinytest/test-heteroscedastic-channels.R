# The exported channels that report a residual scale read a heteroscedastic
# fit's per-observation s(x), not the scalar sigma it also carries - which
# under this parameterization is a fixed unit residual times the response
# range, a constant with no posterior content. loglik scores at s(x),
# extract/predict draw their posterior predictive noise at it, and summary()
# reports it in place of that constant.

source(
  system.file("common", "captureWarnings.R", package = "dbarts"),
  local = TRUE
)

set.seed(21, sample.kind = "Rejection")

n <- 200L
x <- runif(n)
y <- 2 * x + ifelse(x < 0.5, 0.3, 1.5) * rnorm(n)
x.test <- matrix(seq.int(0.01, 0.99, length.out = 40L), 40L, 1L)

fit <- bart2(
  x,
  y,
  test = x.test,
  variance = varianceForest(n.trees = 20L),
  n.trees = 25L,
  n.samples = 200L,
  n.burn = 100L,
  n.chains = 1L,
  keepTrees = TRUE,
  verbose = FALSE
)

# the placeholder this whole file is about: sigma carries no posterior content
expect_equal(length(unique(fit$sigma)), 1L)
expect_equal(fit$sigma[1L], diff(range(y)))

# ---- loglik scores at s(x) ----
ev <- extract(fit, type = "ev")
warnings.loglik <- captureWarnings(loglik <- extract(fit, type = "loglik"))
expect_equal(length(warnings.loglik), 0L)

n.draws <- nrow(ev)
expected <- matrix(
  dnorm(
    rep(y, each = n.draws),
    as.vector(ev),
    as.vector(fit$s.train),
    log = TRUE
  ),
  n.draws,
  n
)
expect_equal(loglik, expected, tolerance = 1e-12)

# and is nowhere near what the scalar gave: the defect was a 3x error in the
# total, not a rounding difference
atScalar <- dnorm(
  rep(y, each = n.draws),
  as.vector(ev),
  fit$sigma[1L],
  log = TRUE
)
expect_true(sum(expected) - sum(atScalar) > 1000)

# ---- a case weight is a precision multiplier on s(x), not on sigma ----
set.seed(22, sample.kind = "Rejection")
w <- runif(n, 0.5, 2)
fitWeighted <- bart2(
  x,
  y,
  weights = w,
  variance = varianceForest(n.trees = 20L),
  n.trees = 25L,
  n.samples = 100L,
  n.burn = 100L,
  n.chains = 1L,
  verbose = FALSE
)
evWeighted <- extract(fitWeighted, type = "ev")
n.drawsWeighted <- nrow(evWeighted)
expect_equal(
  extract(fitWeighted, type = "loglik"),
  matrix(
    dnorm(
      rep(y, each = n.drawsWeighted),
      as.vector(evWeighted),
      as.vector(fitWeighted$s.train) / rep(sqrt(w), each = n.drawsWeighted),
      log = TRUE
    ),
    n.drawsWeighted,
    n
  ),
  tolerance = 1e-12
)

# ---- the draw and the scale pair correctly across chains ----
set.seed(23, sample.kind = "Rejection")
fitChains <- bart2(
  x,
  y,
  variance = varianceForest(n.trees = 20L),
  n.trees = 25L,
  n.samples = 50L,
  n.burn = 50L,
  n.chains = 2L,
  verbose = FALSE
)
evChains <- extract(fitChains, type = "ev")
expect_equal(
  extract(fitChains, type = "loglik"),
  matrix(
    dnorm(
      rep(y, each = nrow(evChains)),
      as.vector(evChains),
      as.vector(fitChains$s.train),
      log = TRUE
    ),
    nrow(evChains),
    n
  ),
  tolerance = 1e-12
)

# ---- the posterior predictive's noise sd tracks s(x), observation by
# observation ----
set.seed(31, sample.kind = "Rejection")
ppd <- extract(fit, type = "ppd")
set.seed(31, sample.kind = "Rejection")
expect_identical(extract(fit, type = "ppd"), ppd)

noiseSd <- apply(ppd - ev, 2L, sd)
sMean <- apply(fit$s.train, 2L, mean)
expect_true(cor(noiseSd, sMean) > 0.9)
# a per-row scale, not a global one: the ratio has no trend against s(x)
expect_true(abs(cor(noiseSd / sMean, sMean)) < 0.3)
# and it is nowhere near the constant it used to be drawn at
expect_true(mean(noiseSd) < 0.25 * fit$sigma[1L])

# ---- predict(type = "ppd") draws at s(x) too, at the new rows ----
predEv <- predict(fit, x.test, type = "ev")
sTest <- attr(predEv, "s")
expect_false(is.null(sTest))
set.seed(32, sample.kind = "Rejection")
predPpd <- predict(fit, x.test, type = "ppd")
set.seed(32, sample.kind = "Rejection")
expect_identical(predict(fit, x.test, type = "ppd"), predPpd)

noiseSd.test <- apply(predPpd - predEv, 2L, sd)
sMean.test <- apply(sTest, 2L, mean)
expect_true(cor(noiseSd.test, sMean.test) > 0.9)
expect_true(mean(noiseSd.test) < 0.25 * fit$sigma[1L])

# the stored test channel draws at the same scale as the replay
set.seed(33, sample.kind = "Rejection")
ppd.test <- extract(fit, type = "ppd", sample = "test")
noiseSd.stored <- apply(
  ppd.test - extract(fit, type = "ev", sample = "test"),
  2L,
  sd
)
expect_true(cor(noiseSd.stored, apply(fit$s.test, 2L, mean)) > 0.9)

# ---- summary() reports s(x), never the constant ----
fitSummary <- summary(fit)
expect_true("mean.s" %in% fitSummary$stats$variable)
expect_false("sigma" %in% fitSummary$stats$variable)
meanS <- fitSummary$stats[fitSummary$stats$variable == "mean.s", ]
expect_equal(meanS$mean, mean(fit$s.train), tolerance = 1e-8)
# the constant's giveaway was sd 0 and an undefined R-hat
expect_true(meanS$sd > 0)

# a homoscedastic fit keeps reporting sigma
set.seed(41, sample.kind = "Rejection")
fitHom <- bart2(
  x,
  y,
  n.trees = 25L,
  n.samples = 50L,
  n.burn = 50L,
  n.chains = 1L,
  keepTrees = TRUE,
  verbose = FALSE
)
expect_true("sigma" %in% summary(fitHom)$stats$variable)
expect_false("mean.s" %in% summary(fitHom)$stats$variable)

# ---- refusals: a scale the fit does not carry is named, not substituted ----
fitNoTestScale <- fit
fitNoTestScale$s.test <- NULL
noTestScaleRefusal <- tryCatch(
  extract(fitNoTestScale, type = "ppd", sample = "test"),
  error = function(e) conditionMessage(e)
)
expect_identical(
  noTestScaleRefusal,
  paste0(
    "posterior predictive sampling is not available at the test rows of a ",
    "heteroscedastic fit that stores no 's.test' draws"
  )
)

fitNoSurface <- fitHom
fitNoSurface$s.train <- fit$s.train
noSurfaceRefusal <- tryCatch(
  predict(fitNoSurface, x.test, type = "ppd"),
  error = function(e) conditionMessage(e)
)
expect_identical(
  noSurfaceRefusal,
  paste0(
    "posterior predictive sampling is not available on a heteroscedastic ",
    "fit whose sampler replays no variance surface"
  )
)
