# Frequency weights reach the engine through the model record's weights.type:
# a gaussian fit whose whole weights are read as counts is the fit to its rows
# written out that many times, draw for draw, once both start from one stated
# sigma estimate. Here the record is written by hand; the argument that
# writes it, and its refusals, belong to the fitting functions.

set.seed(20261106L)
n <- 30L
x <- matrix(runif(2L * n), n, 2L, dimnames = list(NULL, c("x1", "x2")))
y <- sin(6 * x[, 1L]) + 2 * x[, 2L] + rnorm(n, 0, 0.2)
w <- as.double(sample(1:3, n, replace = TRUE))
rows <- rep(seq_len(n), w)

control <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.samples = 60L,
  n.burn = 0L,
  updateState = FALSE,
  seed = 11L
)

sampler <- function(data, weights.type = NA_character_, ...) {
  spec <- dbarts::dbartsSpec(
    data,
    control = control,
    n.trees = 20L,
    sigest = 0.5,
    ...
  )
  spec$model@weights.type <- weights.type
  methods::validObject(spec$model)
  new("dbartsSampler", spec$control, spec$model, spec$data)
}

weighted <- sampler(dbarts::dbartsData(x, y, weights = w), "frequency")
expanded <- sampler(dbarts::dbartsData(x[rows, , drop = FALSE], y[rows]))
precision <- sampler(dbarts::dbartsData(x, y, weights = w))
drawsWeighted <- weighted$run()
drawsExpanded <- expanded$run()
drawsPrecision <- precision$run()

expect_equal(
  drawsWeighted$sigma,
  drawsExpanded$sigma,
  tolerance = 1e-12,
  info = "a frequency gaussian fit draws sigma as its expansion does"
)
expect_equal(
  drawsWeighted$train[rows, , drop = FALSE],
  drawsExpanded$train,
  tolerance = 1e-12,
  check.attributes = FALSE,
  info = "a frequency gaussian fit draws each row's fit as its copies do"
)
expect_true(
  max(abs(drawsPrecision$sigma - drawsWeighted$sigma)) > 1e-3,
  info = "the same weights read as precisions are another model"
)

# all-ones weights are one model under either reading, bitwise, Student-t's
# drawn df included
ones <- rep(1, n)
for (family in c("gaussian", "student")) {
  a <- sampler(
    dbarts::dbartsData(x, y, weights = ones),
    "frequency",
    family = family
  )$run()
  b <- sampler(
    dbarts::dbartsData(x, y, weights = ones),
    "precision",
    family = family
  )$run()
  expect_identical(a$sigma, b$sigma, info = family)
  expect_identical(a$train, b$train, info = family)
}

# a gp leaf's default lengthscale under frequency, rows at weight 0 included,
# is the median pairwise distance over the rows written out
wZero <- w
wZero[c(3L, 11L, 20L)] <- 0
rowsZero <- rep(seq_len(n), wZero)
expect_warning(
  dataZero <- dbarts::dbartsData(x, y, weights = wZero),
  "left out of the likelihood"
)
lengthscaleOf <- function(s) {
  s$storeState()
  s$state[[1L]]$forests[[1L]]$leaf.lengthscales
}
gpPrior <- dbarts::dbartsPriors$gp(1L)
expect_equal(
  lengthscaleOf(sampler(dataZero, "frequency", leaf.prior = gpPrior)),
  lengthscaleOf(sampler(
    dbarts::dbartsData(x[rowsZero, , drop = FALSE], y[rowsZero]),
    leaf.prior = gpPrior
  )),
  tolerance = 1e-12,
  info = "a frequency gp lengthscale is the expansion's"
)

# Student-t at whole weights of 2 and 3 draws the df its expansion draws, in
# distribution: a row's latent is the sum of its copies', and the df step's
# per-row constant for a count above 1 is what keeps the two together
set.seed(20261107L)
nT <- 150L
xT <- matrix(runif(2L * nT), nT, 2L, dimnames = list(NULL, c("x1", "x2")))
yT <- sin(6 * xT[, 1L]) + 2 * xT[, 2L] + 0.4 * rt(nT, 8)
wT <- as.double(sample(2:3, nT, replace = TRUE))
rowsT <- rep(seq_len(nT), wT)
studentDf <- function(data, weights.type) {
  spec <- dbarts::dbartsSpec(
    data,
    control = dbarts::dbartsControl(
      n.chains = 2L,
      n.threads = 1L,
      n.samples = 300L,
      n.burn = 100L,
      updateState = FALSE,
      seed = 3L
    ),
    n.trees = 10L,
    sigest = 0.5,
    family = "student"
  )
  spec$model@weights.type <- weights.type
  mean(new("dbartsSampler", spec$control, spec$model, spec$data)$run()$resid.df)
}
dfWeighted <- studentDf(dbarts::dbartsData(xT, yT, weights = wT), "frequency")
dfExpanded <- studentDf(
  dbarts::dbartsData(xT[rowsT, , drop = FALSE], yT[rowsT]),
  NA_character_
)
expect_true(
  abs(dfWeighted - dfExpanded) < 1.5,
  info = "a frequency Student-t fit at whole weights draws its expansion's df"
)

expect_error(
  {
    model <- weighted$model
    model@weights.type <- "count"
    methods::validObject(model)
  },
  "'weights.type' must be NA, \"precision\" or \"frequency\""
)

# a numeric basis that varies over rows whose counts total at most 1 has no
# scale to read a default from; the refusal and the kept-scale warning say so
# rather than calling the column constant
column <- runif(n)
basisSampler <- function(weights) {
  spec <- dbarts::dbartsSpec(
    dbarts::dbartsData(x, y, weights = weights),
    control = control,
    forests = list(forest(), do.call(forest, list(basis = column))),
    sigest = 0.5
  )
  spec$model@weights.type <- "frequency"
  new("dbartsSampler", spec$control, spec$model, spec$data)
}
expect_error(
  basisSampler(rep(1 / n, n)),
  paste0(
    "forest 2 states no sd and the rows its basis is read over carry a ",
    "total count of at most 1, so there is no standard deviation to take a ",
    "default from"
  ),
  fixed = TRUE,
  info = "a varying basis over counts totaling 1 is refused by its total"
)
fractional <- basisSampler(rep(1, n))
fractional$setWeights(rep(1 / n, n))
expect_warning(
  fractional$setForestBasis(2L, runif(n), updateBasisScale = TRUE),
  paste0(
    "the rows in the likelihood carry a total count of at most 1, so the ",
    "scale of forest 2's basis is kept as it was"
  ),
  fixed = TRUE,
  info = "a rescale over counts totaling 1 keeps the scale, saying why"
)
