# A count response is replaced as its formula writes it (dec-B628): on a
# binomial sampler $setResponse takes the n x 2 successes and failures, or a
# vector of successes against the trials held, and $setData takes new data
# whose counts state both, the number of rows free; on a multinomial sampler
# $setResponse takes the n x K counts and refuses a vector. A binomial sampler
# is built from a data object whose counts slot is written by hand, as the
# fitting functions do not yet take such a response.

set.seed(20261108L)
n <- 20L
x <- matrix(runif(2L * n), n, 2L, dimnames = list(NULL, c("x1", "x2")))
control <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.samples = 5L,
  n.burn = 0L,
  updateState = FALSE,
  seed = 3L
)
binaryOf <- function(x) {
  dbarts::dbarts(
    x,
    as.double(seq_len(nrow(x)) %% 2L),
    family = "logistic",
    control = control,
    n.trees = 5L
  )
}
binary <- binaryOf(x)
withCounts <- function(data, counts) {
  storage.mode(counts) <- "integer"
  data@counts <- counts
  data@y <- as.double(rowSums(counts))
  data
}
binomialOf <- function(data) {
  new("dbartsSampler", binary$control, binary$model, data)
}
drawsOf <- function(sampler) sampler$run(0L, 20L)$train

trials <- as.integer(1L + (seq_len(n) %% 4L))
successes <- as.integer(rbinom(n, trials, plogis(2 * x[, 1L] - 1)))
counts <- cbind(successes, trials - successes)
newSuccesses <- as.integer(rbinom(n, trials, 0.5))
newCounts <- cbind(newSuccesses, trials - newSuccesses)

# ---- $setResponse on a binomial sampler ---------------------------------------

# the two columns: the next sweep draws the latents under the new counts
# before the trees read them, so the sampler draws what one built on them
# draws, the latents in force at the call standing until then
fresh <- binomialOf(withCounts(binary$data, newCounts))
swapped <- binomialOf(withCounts(binary$data, counts))
before <- swapped$getLatents()
swapped$setResponse(newCounts)
expect_identical(swapped$getLatents(), before)
expect_identical(unname(swapped$data@counts), unname(fresh$data@counts))
expect_identical(swapped$data@y, as.double(trials))
expect_identical(drawsOf(swapped), drawsOf(fresh))

# new trials too
moreTrials <- cbind(successes + 2L, trials - successes + 1L)
fresh <- binomialOf(withCounts(binary$data, moreTrials))
swapped <- binomialOf(withCounts(binary$data, counts))
swapped$setResponse(moreTrials)
expect_identical(drawsOf(swapped), drawsOf(fresh))

# a vector gives the successes, each row's trials held
fresh <- binomialOf(withCounts(binary$data, newCounts))
swapped <- binomialOf(withCounts(binary$data, counts))
swapped$setResponse(newSuccesses)
expect_identical(unname(swapped$data@counts), unname(newCounts))
expect_identical(drawsOf(swapped), drawsOf(fresh))

# a swap after sweeps draws nothing and moves the fit
swapped <- binomialOf(withCounts(binary$data, counts))
startDraws <- drawsOf(swapped)
before <- swapped$getLatents()
swapped$setResponse(trials)
expect_identical(swapped$getLatents(), before)
expect_true(mean(drawsOf(swapped)) > mean(startDraws))

# the refusals
swapped <- binomialOf(withCounts(binary$data, counts))
expect_error(
  swapped$setResponse(replace(successes, c(2L, 7L), trials[c(2L, 7L)] + 1L)),
  paste0(
    "$setResponse: a vector gives a binomial sampler's successes, each row's ",
    "trials held; row(s) 2, 7 hold more successes than trials"
  ),
  fixed = TRUE
)
expect_error(
  swapped$setResponse(successes[-1L]),
  "$setResponse: a binomial sampler's successes must be numeric, one per row (20)",
  fixed = TRUE
)
expect_error(
  swapped$setResponse(cbind(counts, 0L)),
  "$setResponse: a binomial sampler's response is an n x 2 matrix of successes and failures, n = 20",
  fixed = TRUE
)
expect_error(
  swapped$setResponse(replace(counts, 4L, -1L)),
  paste0(
    "$setResponse: a binomial response's successes and failures must be ",
    "non-negative whole numbers; row(s) 4"
  ),
  fixed = TRUE
)
expect_error(
  swapped$setResponse(replace(successes, 3L, 0.5)),
  "non-negative whole numbers; row(s) 3",
  fixed = TRUE
)
expect_error(
  swapped$setResponse(replace(successes, 5L, NA)),
  "non-negative whole numbers; row(s) 5",
  fixed = TRUE
)
expect_error(
  swapped$setResponse(successes, status = rep(1, n)),
  "$setResponse: 'status' is an aft (survival) sampler's censoring structure; a count sampler takes none",
  fixed = TRUE
)
# a refusal leaves the response as it was
expect_identical(unname(swapped$data@counts), unname(counts))

# ---- $setData on a binomial sampler -------------------------------------------

# new rows and their counts: the sampler draws what one built on the new data
# draws, the number of rows free
m <- 30L
xNew <- matrix(runif(2L * m), m, 2L, dimnames = list(NULL, c("x1", "x2")))
trialsNew <- as.integer(1L + (seq_len(m) %% 3L))
successesNew <- as.integer(rbinom(m, trialsNew, plogis(xNew[, 2L])))
newData <- withCounts(
  binaryOf(xNew)$data,
  cbind(successesNew, trialsNew - successesNew)
)
fresh <- binomialOf(newData)
replaced <- binomialOf(withCounts(binary$data, counts))
replaced$setData(newData)
expect_identical(dim(replaced$data@counts), c(m, 2L))
expect_identical(drawsOf(replaced), drawsOf(fresh))

# the refusals: data that state no counts, and counts on a binary sampler
replaced <- binomialOf(withCounts(binary$data, counts))
expect_error(
  replaced$setData(binaryOf(xNew)$data),
  paste0(
    "$setData: a binomial sampler's data carry its response as the n x 2 ",
    "counts of successes and failures, which 'newData' does not"
  ),
  fixed = TRUE
)
negative <- newData
negative@counts[3L, 2L] <- -1L
expect_error(
  replaced$setData(negative),
  "$setData: a binomial response's successes and failures must be",
  fixed = TRUE
)
plain <- binaryOf(x)
expect_error(
  plain$setData(withCounts(binary$data, counts)),
  paste0(
    "$setData: 'newData' carries counts of successes and failures, which a ",
    "sampler created on a binary response does not take; make a new sampler"
  ),
  fixed = TRUE
)

# ---- $setResponse on a multinomial sampler ------------------------------------

# the n x K counts, against a sampler built on them (the parity arms of
# test-multinomial-counts-mutation.R go further); a vector is refused
categoryCounts <- matrix(rpois(3L * n, 2), n, 3L)
storage.mode(categoryCounts) <- "integer"
otherCounts <- matrix(rpois(3L * n, 2), n, 3L)
storage.mode(otherCounts) <- "integer"
multinomialOf <- function(counts) {
  dbarts::dbarts(
    dbarts::dbartsData(x, counts = counts),
    family = "multinomial",
    control = control,
    n.trees = 5L
  )
}
fresh <- multinomialOf(otherCounts)
multinomial <- multinomialOf(categoryCounts)
multinomial$setResponse(otherCounts)
expect_identical(multinomial$data@counts, fresh$data@counts)
expect_identical(multinomial$data@y, fresh$data@y)
expect_identical(drawsOf(multinomial), drawsOf(fresh))
expect_error(
  multinomial$setResponse(as.double(otherCounts[, 1L])),
  paste0(
    "$setResponse: a multinomial sampler's response is the n x K count ",
    "matrix; a vector cannot say how each row's trials split among the ",
    "categories"
  ),
  fixed = TRUE
)
expect_error(
  multinomial$setCounts(otherCounts),
  "setCounts"
)

rm(
  n,
  x,
  control,
  binaryOf,
  binary,
  withCounts,
  binomialOf,
  drawsOf,
  trials,
  successes,
  counts,
  newSuccesses,
  newCounts,
  fresh,
  swapped,
  before,
  moreTrials,
  startDraws,
  m,
  xNew,
  trialsNew,
  successesNew,
  newData,
  replaced,
  negative,
  plain,
  categoryCounts,
  otherCounts,
  multinomialOf,
  multinomial
)
