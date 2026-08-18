# What a named prior.scale MEANS, measured against prior draws: it is the
# LEAF-PARAMETER scale of the forest total, equal to the prior sd of f(x) at
# every x for the constant leaf only. For the other three leaf models the prior
# of f(x) is x-dependent and prior.sd bounds it in a stated direction, so the
# rows below assert INEQUALITIES where the contract states a bound - an
# equality-only version would pass an engine doing the opposite.
#
# Relative tolerances are stated against the draw count: the standard error of
# an sd estimate is about sd / sqrt(2R), so 1500 draws support a 6% band.

set.seed(19)
n <- 200L
x <- cbind(x1 = runif(n), x2 = runif(n), x3 = runif(n))
y <- 4 * (x[, 1L] - x[, 2L]) + rnorm(n)

namedScale <- 1.5
fixedK <- 2
priorSd <- namedScale / fixedK
# prior.mean is the response transform's shift, which for a continuous
# response with no offset is the midpoint of the observed range
priorMean <- (max(y) + min(y)) / 2
numDraws <- 1500L

priorControl <- function() {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 20L,
    n.samples = 1L,
    updateState = FALSE,
    seed = 13L,
    keepTrees = FALSE
  )
}

# repeated draws of the whole forest from its prior, evaluated at chosen rows
priorDraws <- function(sampler, rows, draws = numDraws) {
  out <- matrix(0, draws, nrow(rows))
  for (draw in seq_len(draws)) {
    sampler$sampleTreesFromPrior(updateState = FALSE)
    sampler$sampleNodeParametersFromPrior(updateState = FALSE)
    out[draw, ] <- sampler$predict(rows)
  }
  out
}

# --- constant leaf: EXACT. prior.sd is the prior sd of f(x) at every x and
# prior.mean is its prior mean. ---
constantSampler <- dbarts(
  x,
  y,
  control = priorControl(),
  node.prior = normal(k = fixedK, scale = namedScale)
)
set.seed(3)
constantDraws <- priorDraws(constantSampler, x[1:4, , drop = FALSE])
expect_true(max(abs(apply(constantDraws, 2L, sd) / priorSd - 1)) < 0.06)
expect_true(max(abs(colMeans(constantDraws) - priorMean)) < 0.15 * priorSd)

# non-vacuity: the same fit without the name inherits the response range, which
# on this y is about three times the named scale (measured 3.02 over 4000
# draws; the 400 draws below carry a Monte Carlo se near 0.1 in these units,
# so the bar sits clear of the value rather than on it)
inheritedSampler <- dbarts(
  x,
  y,
  control = priorControl(),
  node.prior = normal(k = fixedK)
)
set.seed(3)
inheritedDraws <- priorDraws(inheritedSampler, x[1:4, , drop = FALSE], 400L)
expect_true(mean(apply(inheritedDraws, 2L, sd)) / priorSd > 2.5)

# --- linear leaf: prior.sd is a LOWER bound, attained at the standardized
# covariate origin, with sd(f(x)) = prior.sd * sqrt(1 + ||z(x)||^2). The rows
# exercise a constant column (which keeps sd 1 and contributes 0), training
# missingness, and a predict-time NA (which maps to z = 0). ---
xLinear <- x
xLinear[3L, 1L] <- NA_real_
xLinear <- cbind(xLinear, x4 = rep(1, n))
linearSampler <- dbarts(
  xLinear,
  y,
  control = priorControl(),
  node.prior = linear(c("x1", "x2", "x4"), k = fixedK, scale = namedScale)
)
standardize <- function(column, values) {
  observed <- xLinear[, column]
  center <- mean(observed, na.rm = TRUE)
  spread <- sd(observed, na.rm = TRUE)
  if (!(spread > 0)) {
    spread <- 1
  }
  ifelse(is.na(values), 0, (values - center) / spread)
}
linearRows <- rbind(
  # the standardized origin, where the bound is attained
  c(mean(xLinear[, 1L], na.rm = TRUE), mean(xLinear[, 2L]), 0.5, 1),
  c(0.9, 0.1, 0.5, 1),
  # a predict-time NA in a designated column
  c(NA_real_, mean(xLinear[, 2L]), 0.5, 1),
  # the row whose x1 is missing in training
  xLinear[3L, ]
)
colnames(linearRows) <- colnames(xLinear)
linearPredicted <- sqrt(
  1 +
    standardize("x1", linearRows[, 1L])^2 +
    standardize("x2", linearRows[, 2L])^2 +
    standardize("x4", linearRows[, 4L])^2
)
set.seed(3)
linearDraws <- priorDraws(linearSampler, linearRows)
linearMeasured <- apply(linearDraws, 2L, sd) / priorSd
expect_true(max(abs(linearMeasured / linearPredicted - 1)) < 0.07)
# the bound holds in the stated direction, and is not an equality
expect_true(all(linearMeasured > 1 - 0.06))
expect_true(max(linearMeasured) > 1.5)
# prior.mean is exact for the linear leaf
expect_true(max(abs(colMeans(linearDraws) - priorMean)) < 0.5 * priorSd)

# --- gp leaf: prior.sd is an UPPER bound over x, attained at rows reproducing
# a leaf member, decaying to 0 as x leaves the leaf's data cloud, where every
# prior draw equals prior.mean exactly. ---
gpSampler <- dbarts(
  x,
  y,
  control = priorControl(),
  node.prior = gp("x1", k = fixedK, scale = namedScale)
)
gpRows <- rbind(x[1L, ], x[2L, ], x[1L, ], x[1L, ], x[1L, ], x[1L, ])
gpRows[3L, 1L] <- 1.25
gpRows[4L, 1L] <- 1.5
gpRows[5L, 1L] <- 2
gpRows[6L, 1L] <- 20
set.seed(3)
gpDraws <- priorDraws(gpSampler, gpRows)
gpMeasured <- apply(gpDraws, 2L, sd) / priorSd
# member rows attain the bound
expect_true(max(abs(gpMeasured[1:2] - 1)) < 0.06)
# nowhere does it exceed it
expect_true(all(gpMeasured < 1 + 0.06))
# and it decays monotonically away from the training cloud, to nothing
expect_true(all(diff(gpMeasured[3:6]) < 0))
expect_true(gpMeasured[5L] < 0.2)
expect_true(gpMeasured[6L] < 1e-6)
# under full extrapolation the whole prior is the prior mean
expect_true(max(abs(gpDraws[, 6L] - priorMean)) < 1e-8)

# --- monotone leaf: prior.sd is a LOWER bound in the interior, and prior.mean
# is NOT the prior mean of f(x) - the constrained marginal is skew with an
# x-dependent mean tracking the constraint direction. Measured on THIS
# configuration: sd runs 8-14% above the bound and the standardized mean sweeps
# about 2 prior sds across the constrained axis. ---
monotoneSampler <- dbarts(
  x,
  y,
  control = priorControl(),
  monotone = c(x1 = 1),
  node.prior = normal(k = fixedK, scale = namedScale)
)
monotoneRows <- rbind(x[1L, ], x[1L, ], x[1L, ], x[1L, ], x[1L, ])
monotoneRows[, 1L] <- c(0.1, 0.35, 0.5, 0.65, 0.9)
set.seed(3)
monotoneDraws <- priorDraws(monotoneSampler, monotoneRows)
monotoneMeasured <- apply(monotoneDraws, 2L, sd) / priorSd
expect_true(all(monotoneMeasured > 1))
expect_true(max(monotoneMeasured) < 1.25)
standardizedMean <- (colMeans(monotoneDraws) - priorMean) / priorSd
# the constraint is increasing, so the mean tracks it
expect_true(all(diff(standardizedMean) > 0))
expect_true(diff(range(standardizedMean)) > 1)

# --- the anchor is honored on every family and decoration the single-forest
# conversion serves, each measured against its own response transform. The
# binary families carry latent units (transform 1), so the named scale IS the
# latent-scale prior sd; aft carries log-time units; the grouped decorator
# delegates the transform it wraps. ---
familyDraws <- 800L
familyBand <- 0.09
anchorRow <- x[1L, , drop = FALSE]

# one named sampler per family/decoration, each on its own response
anchorSampler <- function(response, ...) {
  dbarts(
    x,
    response,
    node.prior = normal(k = fixedK, scale = namedScale),
    ...
  )
}
measureAnchor <- function(sampler) {
  set.seed(3)
  sd(priorDraws(sampler, anchorRow, familyDraws)[, 1L])
}

yPositive <- exp(0.5 * (x[, 1L] - x[, 2L]) + rnorm(n, 0, 0.3))
yBinary <- rbinom(n, 1L, pnorm(x[, 1L] - x[, 2L]))
yCounts <- rnbinom(n, size = 4, mu = exp(0.5 * x[, 1L]) * 4)
yOrdered <- factor(
  cut(x[, 1L] - x[, 2L], breaks = c(-Inf, -0.3, 0.3, Inf), labels = FALSE),
  ordered = TRUE
)
groupControl <- priorControl()
attr(groupControl, "bartcore.groups") <- list(
  indices = rep_len(c(1L, 2L), n),
  n.groups = 2L,
  prior = "cauchy",
  rel.scale = 1,
  n.steps = 1L
)

anchorSamplers <- list(
  gaussian = anchorSampler(y, control = priorControl()),
  # resid.dist is NSE and cannot be forwarded through a helper's dots
  student = dbarts(
    x,
    y,
    control = priorControl(),
    resid.dist = student(5),
    node.prior = normal(k = fixedK, scale = namedScale)
  ),
  grouped = anchorSampler(y, control = groupControl),
  aft = anchorSampler(
    cbind(yPositive, rep(1L, n)),
    control = priorControl(),
    family = "aft"
  ),
  probit = anchorSampler(yBinary, control = priorControl()),
  logistic = anchorSampler(
    yBinary,
    control = priorControl(),
    family = "logistic"
  ),
  ordinal = anchorSampler(
    yOrdered,
    control = priorControl(),
    family = "ordinal"
  ),
  nbinom = anchorSampler(yCounts, control = priorControl(), family = "nbinom"),
  weightedLogistic = anchorSampler(
    yBinary,
    control = priorControl(),
    family = "logistic",
    weights = rep_len(c(1, 2), n)
  ),
  # the heteroscedastic path is ungated by construction rather than by a
  # guard - the conversion runs at the leaf-scale assignment, before the
  # variance-forest branch, and reads no family flag - so it needs a pin of
  # its own rather than an inherited claim. The variance forest is a
  # separate leaf model outside forests_ and carries no named calibration;
  # this measures the MEAN forest's.
  heteroscedastic = anchorSampler(
    y,
    control = priorControl(),
    variance = varianceForest(vars = ~x1, n.trees = 5L)
  )
)
for (familyName in names(anchorSamplers)) {
  measured <- measureAnchor(anchorSamplers[[familyName]])
  expect_true(
    abs(measured / priorSd - 1) < familyBand,
    info = paste0(familyName, ": ", measured, " vs ", priorSd)
  )
}
