# The bridge reads a binomial response off the data object's counts slot under
# family "logistic": successes in its first column and the trials its row
# sums, a row of k successes in n trials drawing n Polya-Gamma units per sweep.
# The fitting functions do not yet take such a response, so the entries are
# reached directly here.

set.seed(20261107L)
n <- 20L
x <- matrix(runif(2L * n), n, 2L, dimnames = list(NULL, c("x1", "x2")))
y <- as.double(runif(n) < plogis(2 * x[, 1L] - 1))
control <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.samples = 5L,
  n.burn = 0L,
  updateState = FALSE,
  seed = 3L
)
binary <- dbarts::dbarts(
  x,
  y,
  family = "logistic",
  control = control,
  n.trees = 5L
)
trials <- as.integer(1L + (seq_len(n) %% 4L))
successes <- as.integer(rbinom(n, trials, plogis(2 * x[, 1L] - 1)))
withCounts <- function(counts, weights = NULL) {
  data <- binary$data
  storage.mode(counts) <- "integer"
  data@counts <- counts
  data@y <- as.double(rowSums(counts))
  data@weights <- weights
  data
}
binomialData <- withCounts(cbind(successes, trials - successes))
create <- function(data, model = binary$model) {
  .Call(
    dbarts:::C_dbarts_bartcore_create,
    binary$control,
    model,
    data,
    "logistic",
    NULL,
    NULL
  )
}
run <- function(pointer, numSamples = 5L) {
  .Call(
    dbarts:::C_dbarts_bartcore_run,
    pointer,
    0L,
    numSamples,
    NULL,
    NULL,
    TRUE
  )
}
latents <- function(pointer) {
  .Call(dbarts:::C_dbarts_bartcore_getLatents, pointer, NULL)
}

binomial <- create(binomialData)
draws <- run(binomial, 20L)$train
expect_true(
  all(is.finite(draws)),
  info = "binomial: a sampler built from the counts slot runs"
)

# one trial a row is the binary response, draw for draw
expect_identical(
  run(create(withCounts(cbind(y, 1 - y))), 20L)$train,
  run(create(binary$data), 20L)$train,
  info = "binomial: one trial a row draws what the binary response draws"
)

# a fold over every row is the sampler over them, its trials carried with
# its successes
handle <- .Call(
  dbarts:::C_dbarts_bartcore_createDataHandle,
  binary$control,
  binomialData,
  NULL
)
fold <- .Call(
  dbarts:::C_dbarts_bartcore_createFromHandle,
  binary$control,
  binary$model,
  binomialData,
  handle,
  seq_len(n),
  NULL,
  "logistic",
  NULL
)
expect_identical(
  run(fold, 20L)$train,
  run(create(binomialData), 20L)$train,
  info = "binomial: a fold over every row draws what the sampler draws"
)

# a fold over a proper subset reads its own rows' trials: counts changed
# outside the fold leave its draws alone, counts changed inside move them
foldOf <- function(data, rows) {
  handle <- .Call(
    dbarts:::C_dbarts_bartcore_createDataHandle,
    binary$control,
    data,
    NULL
  )
  .Call(
    dbarts:::C_dbarts_bartcore_createFromHandle,
    binary$control,
    binary$model,
    data,
    handle,
    rows,
    NULL,
    "logistic",
    NULL
  )
}
foldRows <- 11:20
counts <- cbind(successes, trials - successes)
outerChanged <- counts
outerChanged[1:10, ] <- cbind(rep(7L, 10L), rep(5L, 10L))
innerChanged <- counts
innerChanged[foldRows, 2L] <- innerChanged[foldRows, 2L] + 3L
foldDraws <- run(foldOf(withCounts(counts), foldRows), 20L)$train
expect_identical(
  run(foldOf(withCounts(outerChanged), foldRows), 20L)$train,
  foldDraws,
  info = "binomial: a fold's draws ignore the counts of rows outside it"
)
expect_false(
  identical(
    run(foldOf(withCounts(innerChanged), foldRows), 20L)$train,
    foldDraws
  ),
  info = "binomial: a fold's draws follow the trials of its own rows"
)

# the counts' refusals
negative <- cbind(successes, trials - successes)
negative[c(2L, 5L), 2L] <- -1L
expect_error(
  create(withCounts(negative)),
  paste0(
    "a binomial response's successes and failures must be non-negative ",
    "whole numbers; row(s) 2, 5"
  ),
  fixed = TRUE,
  info = "binomial: a negative failure count is refused, naming its rows"
)
missing <- cbind(successes, trials - successes)
missing[3L, 1L] <- NA_integer_
expect_error(
  create(withCounts(missing)),
  "row(s) 3",
  fixed = TRUE,
  info = "binomial: a missing count is refused"
)
expect_error(
  create(withCounts(cbind(successes, trials - successes, 1L))),
  paste0(
    "the data carry an n x 3 count matrix: family \"logistic\" takes two ",
    "columns, successes and failures, and family \"multinomial\" fits 3 ",
    "categories"
  ),
  fixed = TRUE,
  info = "binomial: a count matrix of three columns is refused"
)
large <- cbind(successes, trials - successes)
large[1L, ] <- c(600000L, 400001L)
expect_error(
  create(withCounts(large)),
  "family \"logistic\" requires the trials no larger than 1000000",
  fixed = TRUE,
  info = "binomial: trials over one million are refused"
)
half <- cbind(successes, trials - successes)
half[1L, ] <- c(300000L, 300000L)
twice <- replace(rep(1, n), 1L, 2)
capText <- "'weights' times the trials no larger than 1000000"
expect_true(
  inherits(create(withCounts(half)), "externalptr"),
  info = "binomial: 600,000 trials at a weight of 1 are taken"
)
expect_error(
  create(withCounts(half, twice)),
  capText,
  fixed = TRUE,
  info = "binomial: 600,000 trials at a weight of 2 are refused"
)
expect_error(
  create(withCounts(
    cbind(successes, trials - successes),
    replace(twice, 2L, 0)
  )),
  "logistic weights are observation counts",
  fixed = TRUE,
  info = "binomial: a weight of 0 is refused, as on a binary response"
)

# $setCounts stores and draws nothing; the next sweep draws under the counts
setCounts <- function(pointer, counts) {
  storage.mode(counts) <- "integer"
  .Call(dbarts:::C_dbarts_bartcore_setCounts, pointer, counts)
}
before <- latents(binomial)
setCounts(binomial, cbind(trials, 0L))
expect_identical(
  latents(binomial),
  before,
  info = "binomial: $setCounts draws nothing"
)
expect_true(
  mean(run(binomial, 20L)$train) > mean(draws),
  info = "binomial: every trial a success moves the fit up"
)
expect_error(
  setCounts(binomial, cbind(trials, 0L, 0L)),
  "$setCounts: requires an integer matrix of 20 observations x 2",
  fixed = TRUE,
  info = "binomial: $setCounts to three columns is refused"
)
expect_error(
  setCounts(binomial, replace(cbind(trials, 0L), 4L, -2L)),
  "$setCounts: a binomial response's successes and failures must be",
  fixed = TRUE,
  info = "binomial: $setCounts to a negative count is refused"
)
expect_error(
  setCounts(binomial, large),
  "$setCounts: family \"logistic\" requires the trials no larger than 1000000",
  fixed = TRUE,
  info = "binomial: $setCounts to trials over one million is refused"
)
zeroTrials <- cbind(successes, trials - successes)
zeroTrials[c(6L, 9L), ] <- 0L
setCounts(binomial, zeroTrials)
expect_true(
  all(is.finite(run(binomial, 5L)$train)),
  info = "binomial: a sampler holding rows of no trials runs"
)

setWeights <- function(pointer, weights) {
  .Call(
    dbarts:::C_dbarts_bartcore_setWeights,
    pointer,
    as.double(weights),
    FALSE,
    FALSE,
    NULL
  )
}
weighted <- create(withCounts(half))
expect_error(
  setWeights(weighted, twice),
  paste0("$setWeights: family \"logistic\" requires ", capText),
  fixed = TRUE,
  info = "binomial: $setWeights to a weight of 2 at 600,000 trials is refused"
)
setWeights(weighted, replace(rep(1, n), 2L, 2))
halfTwice <- half
halfTwice[2L, ] <- c(300000L, 300000L)
expect_error(
  setCounts(weighted, halfTwice),
  capText,
  fixed = TRUE,
  info = "binomial: $setCounts to 600,000 trials at a weight of 2 is refused"
)

# a flat response or a whole data replacement cannot state the counts
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_setResponse,
    binomial,
    as.double(successes),
    FALSE,
    NULL
  ),
  paste0(
    "$setResponse: its response is the n x 2 count pair of successes and ",
    "failures, which a flat vector cannot express; replace it with $setCounts"
  ),
  fixed = TRUE,
  info = "binomial: $setResponse is refused"
)
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_setData,
    binomial,
    binomialData,
    "none",
    "keep",
    "keep",
    NULL,
    NULL
  ),
  "$setData: a binomial sampler's rows and their trials are fixed at creation",
  fixed = TRUE,
  info = "binomial: $setData is refused"
)
