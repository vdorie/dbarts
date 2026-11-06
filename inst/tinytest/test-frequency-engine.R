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

expect_error(
  {
    model <- weighted$model
    model@weights.type <- "count"
    methods::validObject(model)
  },
  "'weights.type' must be NA, \"precision\" or \"frequency\""
)
