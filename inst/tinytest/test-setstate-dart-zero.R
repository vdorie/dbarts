# A state whose DART weight is 0 on a column the sampler may split on does not
# fit cleanly: setState returns FALSE with the sampler untouched, and a forced
# one installs it with the DART weights started fresh, equal over the columns
# the sampler may split on. Negative, missing and all-zero weights stay errors.

set.seed(2L)
n <- 80L
d <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
d$y <- d$x1 + d$x2 + rnorm(n, 0, 0.3)
control <- dbartsControl(n.chains = 1L, n.threads = 1L, updateState = FALSE)
make <- function(...) {
  dbarts(
    y ~ x1 + x2 + x3,
    d,
    split.probs = dart(),
    control = control,
    n.trees = 10L,
    ...
  )
}
weightsOf <- function(sampler) {
  sampler$storeState()
  sampler$state[[1L]]$dart.probabilities
}

restricted <- make(forests = forest(~ x1 + x2))
invisible(restricted$run(20L, 5L))
restricted$storeState()
state <- restricted$state
expect_identical(state[[1L]]$dart.probabilities[3L], 0)

# a sampler that may split on every column
full <- make()
invisible(full$run(5L, 5L))
before <- weightsOf(full)
full$storeState()
stateBefore <- full$state

expect_false(full$setState(state))
# untouched: the stored state is as it was
full$storeState()
expect_identical(full$state, stateBefore)
expect_identical(weightsOf(full), before)

# forced: installed (NULL, invisibly), the weights fresh and equal, the
# concentration back at its initial value and the delay counter restarted
shifted <- state
shifted[[1L]]$dart.alpha <- 7
shifted[[1L]]$dart.updates.skipped <- 0L
expect_null(full$setState(shifted, forceUpdate = TRUE))
expect_equal(weightsOf(full), rep(1 / 3, 3L))
full$storeState()
expect_equal(full$state[[1L]]$dart.alpha, state[[1L]]$dart.alpha)
expect_true(full$state[[1L]]$dart.updates.skipped > 0L)
full$setState(state, forceUpdate = TRUE)
# and it runs, the third column splittable again
samples <- full$run(5L, 5L)
expect_true(all(is.finite(samples$train)))

# the sampler the weights came from may leave the column: clean
twin <- make(forests = forest(~ x1 + x2))
expect_true(twin$setState(state))

# a state that fails the floors stays an error, forced or not
broken <- function(weights) {
  bad <- state
  bad[[1L]]$dart.probabilities <- weights
  bad
}
for (weights in list(c(0.5, 0.6, -0.1), c(0.5, NA, 0.5), c(0, 0, 0))) {
  for (force in c(FALSE, TRUE)) {
    expect_error(
      full$setState(broken(weights), forceUpdate = force),
      "not consistent",
      info = paste(weights, collapse = ",")
    )
  }
}
