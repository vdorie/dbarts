# Inputs that never produced a valid fit are refused by name: an infinite case
# weight on every weighted entry, a per-column cut count below one on the data
# object's slot. Also pins pdbart's burn-in sigma.

set.seed(0)
n <- 20L
x <- matrix(rnorm(2L * n), n)
y <- x[, 1L] + rnorm(n)
yBinary <- as.integer(y > 0)
df <- data.frame(y = y, yBinary = yBinary, x)

control <- dbartsControl(
  n.samples = 1L,
  n.burn = 0L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)

## a logistic fit's weight is its Polya-Gamma count, so an infinite one never
## finishes a sweep; each family is refused at ingestion instead
for (bad in c(Inf, -Inf)) {
  w <- replace(rep(1, n), 3L, bad)
  rule <- if (bad > 0) "'weights' must all be finite" else "non-negative"
  for (family in c("gaussian", "logistic", "probit")) {
    response <- if (family == "gaussian") y else yBinary
    expect_error(
      bart(
        x,
        response,
        weights = w,
        family = family,
        n.samples = 1L,
        n.burn = 1L,
        n.chains = 1L,
        n.threads = 1L,
        verbose = FALSE
      ),
      rule,
      info = paste("bart", family, bad)
    )
    expect_error(
      dbarts(x, response, weights = w, family = family, control = control),
      rule,
      info = paste("dbarts", family, bad)
    )
  }
  for (response in list(y, yBinary)) {
    expect_error(
      xbart(
        x,
        response,
        weights = w,
        n.samples = 1L,
        n.burn = c(1L, 1L),
        n.reps = 1L,
        n.threads = 1L,
        verbose = FALSE
      ),
      rule
    )
    expect_error(
      bartBT(
        x,
        response,
        weights = w,
        ndpost = 1L,
        nskip = 1L,
        verbose = FALSE
      ),
      rule
    )
  }
  # the deprecation warning is the only one muffled
  expect_error(
    suppressWarnings(
      rbart_vi(
        y ~ X1 + X2,
        df,
        weights = w,
        group.by = rep(1:2, n / 2L),
        n.samples = 1L,
        n.burn = 1L,
        n.thin = 1L,
        n.chains = 1L,
        n.threads = 1L,
        verbose = FALSE
      ),
      classes = "deprecatedWarning"
    ),
    rule
  )
  expect_error(dbartsData(x, y, weights = w), rule)
  expect_error(dbartsData(y ~ X1 + X2, df, weights = w), rule)
}

## a live sampler's setter, and the bridge behind it for a logistic sampler,
## whose own R check is the same rule
gaussianSampler <- dbarts(x, y, control = control)
logisticSampler <- dbarts(x, yBinary, family = "logistic", control = control)
for (sampler in list(gaussianSampler, logisticSampler)) {
  expect_error(
    sampler$setWeights(replace(rep(1, n), 3L, Inf)),
    "'weights' must all be finite"
  )
  expect_error(
    sampler$setWeights(replace(rep(1, n), 3L, -Inf)),
    "non-negative"
  )
}
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_setWeights,
    logisticSampler$getPointer(),
    replace(rep(1, n), 3L, Inf)
  ),
  "logistic weights are observation counts"
)
rm(gaussianSampler, logisticSampler, sampler)

## a data object's slot edited past its validity check, handed to the two
## routes that take a data object, is refused before the starting sigma
## estimate rather than failing on it
for (bad in c(Inf, -Inf)) {
  rule <- if (bad > 0) "'weights' must all be finite" else "non-negative"
  data <- dbartsData(x, y)
  data@weights <- replace(rep(1, n), 3L, bad)
  expect_error(dbartsSpec(data, control), rule)
  expect_error(dbarts(data, control = control), rule)
}
rm(data)

## a continuous fit's posterior predictive weights are precision multipliers:
## an infinite one would draw the expected value and a negative one NaN
time <- exp(y)
status <- rep_len(c(1L, 1L, 0L), n)
for (family in c("gaussian", "student", "aft")) {
  fit <- bart(
    x,
    if (family == "aft") cbind(time, status) else y,
    family = if (family == "student") student(5) else family,
    n.samples = 2L,
    n.burn = 2L,
    n.trees = 5L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE,
    keepTrees = TRUE
  )
  for (bad in c(Inf, -1, -Inf)) {
    expect_error(
      predict(fit, x[1:3, ], type = "ppd", weights = c(bad, 1, 1)),
      paste0(
        "the posterior predictive 'weights' of a ",
        family,
        " fit are precision multipliers and must be finite and non-negative"
      ),
      fixed = TRUE
    )
  }
}
rm(fit, time, status)

## a cut count below one reaches the engine only through a data object's
## edited slot, which dbartsSpec keeps; quantile mode would divide by it, and
## uniform mode builds a column no stored state can hold. The control already
## refuses it at its own validity.
data <- dbartsData(x, y)
data@n.cuts <- c(0L, 100L)
expect_error(validObject(data), "'n.cuts' must contain only positive integers")
for (useQuantiles in c(FALSE, TRUE)) {
  specControl <- control
  specControl@useQuantiles <- useQuantiles
  spec <- dbartsSpec(data, specControl)
  expect_identical(spec$data@n.cuts, c(0L, 100L))
  expect_error(
    new("dbartsSampler", spec$control, spec$model, spec$data),
    "'n.cuts' of 0 for predictor 1 is below one"
  )
}
rm(data, specControl, spec)

## pdbart and pd2bart return the burn-in sigma draws a bart fit of the same
## call and seed returns
reference <- bart(
  x,
  y,
  n.samples = 3L,
  n.burn = 7L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 1L,
  verbose = FALSE
)$first.sigma
expect_equal(length(reference), 7L)
fit <- pdbart(
  x,
  y,
  xind = 1L,
  pl = FALSE,
  n.samples = 3L,
  n.burn = 7L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 1L,
  verbose = FALSE
)
expect_identical(fit$first.sigma, reference)
fit <- pd2bart(
  x,
  y,
  xind = 1:2,
  pl = FALSE,
  n.samples = 3L,
  n.burn = 7L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 1L,
  verbose = FALSE
)
expect_identical(fit$first.sigma, reference)
rm(fit, reference)
