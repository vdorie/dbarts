# The bridge takes nbinom and multinomial weights as counts: a row of count w
# draws w times its units of Polya-Gamma draws per sweep, refused above one
# million for w > 1 at creation and at every weight or count change; at
# w = 1 today's caps stand. The fitting functions still refuse these weights,
# so the entries are reached directly here.

set.seed(20261107L)
n <- 20L
x <- matrix(runif(2L * n), n, 2L, dimnames = list(NULL, c("x1", "x2")))
y <- as.double(rpois(n, 3))
control <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.samples = 5L,
  n.burn = 0L,
  updateState = FALSE,
  seed = 3L
)
counts <- dbarts::dbarts(
  x,
  y,
  family = "nbinom",
  control = control,
  n.trees = 5L
)
create <- function(model, data, family) {
  .Call(
    dbarts:::C_dbarts_bartcore_create,
    counts$control,
    model,
    data,
    family,
    NULL,
    NULL
  )
}
# a refusal where a creation is taken fails its expectation rather than the
# file
created <- function(model, data, family) {
  tryCatch(
    inherits(create(model, data, family), "externalptr"),
    error = function(e) FALSE
  )
}
withCounts <- function(y, weights) {
  data <- counts$data
  data@y <- y
  data@weights <- weights
  data
}
capText <- "'weights' times (y + r) no larger than 1000000"

large <- replace(y, 1L, 999999)
expect_true(
  created(counts$model, withCounts(large, rep(1, n)), "nbinom"),
  info = "nbinom: y = 999,999 at a count of 1 is taken"
)
expect_error(
  create(counts$model, withCounts(large, replace(rep(1, n), 1L, 2)), "nbinom"),
  capText,
  fixed = TRUE,
  info = "nbinom: y = 999,999 at a count of 2 is refused"
)
# a drawn r is bounded by its grid's largest point, 50; a fixed one by itself
middle <- replace(y, 1L, 499975)
twice <- replace(rep(1, n), 1L, 2)
expect_error(
  create(counts$model, withCounts(middle, twice), "nbinom"),
  capText,
  fixed = TRUE,
  info = "nbinom: a drawn r is capped at its grid's largest point"
)
fixed <- counts$model
fixed@shape <- 3
expect_true(
  created(fixed, withCounts(middle, twice), "nbinom"),
  info = "nbinom: a fixed r of 3 takes 2 (499975 + 3)"
)
expect_error(
  create(counts$model, withCounts(y, replace(rep(1, n), 2L, 0.5)), "nbinom"),
  "nbinom weights are observation counts and must be non-negative whole numbers",
  fixed = TRUE,
  info = "nbinom: a fractional count is refused"
)

setWeights <- function(sampler, weights) {
  .Call(
    dbarts:::C_dbarts_bartcore_setWeights,
    sampler$getPointer(),
    as.double(weights),
    FALSE,
    FALSE,
    NULL
  )
}
counts$setResponse(large)
expect_error(
  setWeights(counts, twice),
  capText,
  fixed = TRUE,
  info = "nbinom: $setWeights to a count of 2 at y = 999,999 is refused"
)
counts$setResponse(y)
setWeights(counts, replace(twice, 3L, 0))
expect_error(
  counts$setResponse(large),
  capText,
  fixed = TRUE,
  info = "nbinom: $setResponse to y = 999,999 at a count of 2 is refused"
)
expect_true(
  all(is.finite(counts$run(0L, 5L)$train)),
  info = "nbinom: a sampler holding counts, one of them 0, runs"
)

# $setData swaps y and the weights together; the R method refuses counts
# until the fitting functions take them, so the entry is reached directly
setData <- function(sampler, data) {
  .Call(
    dbarts:::C_dbarts_bartcore_setData,
    sampler$getPointer(),
    data,
    "none",
    "keep",
    "keep",
    NULL,
    NULL
  )
}
swapped <- dbarts::dbarts(
  x,
  y,
  family = "nbinom",
  control = control,
  n.trees = 5L
)
expect_error(
  setData(swapped, withCounts(large, twice)),
  capText,
  fixed = TRUE,
  info = "nbinom: $setData to y = 999,999 at a count of 2 is refused"
)

# a $setModel shape install moves r under counts already installed: a fixed r
# of 3 at y = 499,975 and a count of 2 holds 999,956 units
specOf <- function(family) {
  fit <- dbarts::dbarts(
    x,
    y,
    family = family,
    control = control,
    n.trees = 5L
  )
  attr(fit$model, "family.spec", exact = TRUE)
}
modelWith <- function(sampler, family) {
  model <- sampler$model
  attr(model, "family.spec") <- specOf(family)
  model
}
shaped <- dbarts::dbarts(
  x,
  y,
  family = dbarts::dbartsFamilies$nbinom(shape = 3),
  control = control,
  n.trees = 5L
)
shaped$setResponse(middle)
setWeights(shaped, twice)
expect_error(
  shaped$setModel(modelWith(shaped, dbarts::dbartsFamilies$nbinom())),
  capText,
  fixed = TRUE,
  info = "nbinom: $setModel to a drawn r (bound 50) at y = 499,975 and a count of 2 is refused"
)
expect_error(
  shaped$setModel(
    modelWith(shaped, dbarts::dbartsFamilies$nbinom(shape = 600000))
  ),
  capText,
  fixed = TRUE,
  info = "nbinom: $setModel to a fixed r of 600,000 at a count of 2 is refused"
)
expect_true(
  tryCatch(
    {
      shaped$setModel(
        modelWith(shaped, dbarts::dbartsFamilies$nbinom(shape = 10))
      )
      TRUE
    },
    error = function(e) FALSE
  ),
  info = "nbinom: $setModel to a fixed r of 10 at a count of 2 (999,970) is taken"
)

# multinomial: the counts slot's row totals times the weights
# every row has a trial, so no zero-trial warning
y3 <- matrix(rpois(3L * n, 1), n, 3L) + rep(c(1L, 0L, 0L), each = n)
storage.mode(y3) <- "integer"
categories <- dbarts::dbarts(
  x,
  y3,
  family = "multinomial",
  control = control,
  n.trees = 5L
)
withRows <- function(weights, rowTotal = NULL) {
  data <- categories$data
  if (!is.null(rowTotal)) {
    data@counts[1L, ] <- c(rowTotal, 0L, 0L)
  }
  data@weights <- weights
  data
}
multinomialText <- "'weights' times the 'counts' row totals no larger than 1000000"
expect_true(
  created(categories$model, withRows(rep(1, n), 600000L), "multinomial"),
  info = "multinomial: a row total of 600,000 at a count of 1 is taken"
)
expect_error(
  create(categories$model, withRows(twice, 600000L), "multinomial"),
  multinomialText,
  fixed = TRUE,
  info = "multinomial: a row total of 600,000 at a count of 2 is refused"
)
expect_error(
  create(categories$model, withRows(replace(rep(1, n), 2L, -1)), "multinomial"),
  "multinomial weights are observation counts and must be non-negative whole numbers",
  fixed = TRUE,
  info = "multinomial: a negative count is refused"
)
setWeights(categories, replace(twice, 3L, 0))
expect_true(
  all(is.finite(categories$run(0L, 5L)$train)),
  info = "multinomial: a sampler holding counts, one of them 0, runs"
)
big <- y3
big[1L, ] <- c(600000L, 0L, 0L)
expect_error(
  .Call(dbarts:::C_dbarts_bartcore_setCounts, categories$getPointer(), big),
  multinomialText,
  fixed = TRUE,
  info = "multinomial: $setCounts to a row total of 600,000 at a count of 2 is refused"
)
# the other order: the row total first at a count of 1, then the weights
totals <- dbarts::dbarts(
  x,
  y3,
  family = "multinomial",
  control = control,
  n.trees = 5L
)
.Call(dbarts:::C_dbarts_bartcore_setCounts, totals$getPointer(), big)
expect_error(
  setWeights(totals, twice),
  multinomialText,
  fixed = TRUE,
  info = "multinomial: $setWeights to a count of 2 at a row total of 600,000 is refused"
)
