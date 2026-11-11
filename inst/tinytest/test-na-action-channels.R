# The matrix interface hands na.action the weights and the offset beside the
# response and the predictors, as the formula interface does: na.omit and
# na.exclude drop a row whose weight or offset is missing (na.exclude pads it
# back); the default, na.pass and na.fail leave the missing value to the
# channel's own refusal, in its own words.

set.seed(20261110L)
n <- 40L
x <- cbind(x1 = rnorm(n), x2 = rnorm(n))
y <- x[, 1L] + rnorm(n)
yb <- as.integer(y > 0)
w <- runif(n, 1, 2)
o <- rnorm(n, 0, 0.1)
w[3L] <- NA_real_
o[7L] <- NA_real_

quick <- list(
  n.trees = 5L,
  n.samples = 5L,
  n.burn = 3L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)

# two channels, two actions, two families
for (family in c("gaussian", "probit")) {
  response <- if (family == "gaussian") y else yb
  weights <- if (family == "gaussian") w else ifelse(is.na(w), NA_real_, 1)
  for (action in list(na.omit, na.exclude)) {
    byWeight <- do.call(
      bart,
      c(
        list(
          x,
          response,
          weights = weights,
          family = family,
          na.action = action
        ),
        quick
      )
    )
    expect_equal(length(byWeight$yhat.train) / 5L, n - 1L)
    expect_equal(as.integer(byWeight$na.action), 3L)
    byOffset <- do.call(
      bart,
      c(
        list(x, response, offset = o, family = family, na.action = action),
        quick
      )
    )
    expect_equal(as.integer(byOffset$na.action), 7L)
    both <- do.call(
      bart,
      c(
        list(
          x,
          response,
          weights = weights,
          offset = o,
          family = family,
          na.action = action
        ),
        quick
      )
    )
    expect_equal(as.integer(both$na.action), c(3L, 7L))
  }
}

# the front door's sampler agrees
sampler <- dbarts(
  x,
  y,
  weights = w,
  offset = o,
  na.action = na.omit,
  n.trees = 5L,
  control = dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  )
)
expect_equal(length(sampler$data@y), n - 2L)

# na.exclude pads the dropped row back
padded <- do.call(
  bart,
  c(list(x, y, weights = w, na.action = na.exclude), quick)
)
expect_equal(length(fitted(padded)), n)
expect_true(is.na(fitted(padded)[3L]))
expect_false(anyNA(fitted(padded)[-3L]))

# every other action refuses in the channel's own words
for (action in list(na.pass, na.fail, na.keepPredictors)) {
  expect_error(
    do.call(bart, c(list(x, y, weights = w, na.action = action), quick)),
    "'weights' cannot be NA"
  )
  expect_error(
    do.call(bart, c(list(x, y, offset = o, na.action = action), quick)),
    "'offset' contains missing values"
  )
}
expect_error(
  do.call(bart, c(list(x, y, weights = w), quick)),
  "'weights' cannot be NA"
)

# a single missing weight or offset is not a row's: it reaches the channel's
# own refusal under every action, as on the formula interface
for (action in list(na.omit, na.exclude, na.pass, na.fail)) {
  expect_error(
    do.call(bart, c(list(x, y, weights = NA_real_, na.action = action), quick)),
    "'weights' cannot be NA"
  )
  expect_error(
    do.call(bart, c(list(x, y, offset = NA_real_, na.action = action), quick)),
    "'offset' contains missing values"
  )
}
