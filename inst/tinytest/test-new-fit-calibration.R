# A new fit made from a data object, one that came out of a sampler included,
# calibrates its response-based priors from that data's response: the
# response's scale and sigma's starting estimate, unless the caller states
# sigest. A new fit is a new model; copy() is the same model again.

set.seed(20261111L)
n <- 80L
d <- data.frame(x = rnorm(n), z = rnorm(n))
d$y <- d$x + 0.3 * rnorm(n)
ctl <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.samples = 5L,
  n.burn = 3L,
  updateState = FALSE,
  verbose = FALSE
)
draws <- function(sampler) sampler$run(3L, 5L)$train
newFit <- function(...) dbarts(..., n.trees = 5L, control = ctl, seed = 3L)

s <- newFit(y ~ x, d, sigest = 0.37)
expect_identical(s$data@sigma, 0.37)
s$setResponse(10 * d$y)
# the sampler holds the scale it was made with and its own sigma record
expect_identical(s$data@sigma, 0.37)

fresh <- newFit(y ~ x, transform(d, y = 10 * y))
freshDraws <- draws(fresh)
expect_true(fresh$data@sigma != 0.37)

# dbarts and dbartsSpec start from the fresh fit's range and sigest
fromData <- newFit(s$data)
expect_identical(fromData$data@sigma, fresh$data@sigma)
expect_identical(
  dbarts:::dataResponseRange(fromData$data),
  dbarts:::dataResponseRange(fresh$data)
)
expect_identical(draws(fromData), freshDraws)
spec <- dbartsSpec(s$data, ctl, n.trees = 5L, seed = 3L)
expect_identical(spec$data@sigma, fresh$data@sigma)
expect_identical(
  draws(new("dbartsSampler", spec$control, spec$model, spec$data)),
  freshDraws
)
# bart starts from the same
bartFit <- function(data, ...) {
  bart(
    data,
    ...,
    n.trees = 5L,
    n.samples = 5L,
    n.burn = 3L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE,
    seed = 3L,
    keepTrees = TRUE
  )
}
expect_identical(
  bartFit(s$data)$sigest,
  bartFit(dbartsData(y ~ x, transform(d, y = 10 * y)))$sigest
)
expect_identical(
  bartFit(s$data)$yhat.train,
  bartFit(dbartsData(y ~ x, transform(d, y = 10 * y)))$yhat.train
)

# a stated sigest wins, a rule and "auto" estimate
expect_identical(newFit(s$data, sigest = 0.5)$data@sigma, 0.5)
expect_identical(bartFit(s$data, sigest = 0.5)$sigest, 0.5)
expect_identical(
  dbartsSpec(s$data, ctl, sigest = 0.5)$data@sigma,
  0.5
)
for (rule in list("auto", "dense", NULL)) {
  expect_identical(newFit(s$data, sigest = rule)$data@sigma, fresh$data@sigma)
  expect_identical(
    dbartsSpec(s$data, ctl, sigest = rule)$data@sigma,
    fresh$data@sigma
  )
}

# copy() is the same model again
copied <- s$copy()
expect_identical(copied$data@sigma, 0.37)

# a plain data object is estimated as before
plain <- dbartsData(y ~ x, d)
expect_identical(
  newFit(plain)$data@sigma,
  newFit(y ~ x, d)$data@sigma
)
