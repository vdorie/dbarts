# A cut grid changed with setCutPoints can be put back on any design, and
# leaves nothing behind: a later setData derives the grid n.cuts names
# whatever grid was set before it, as a copy and a reload do. The grid in
# force is read from the stored state.

cutPointsOf <- function(sampler) {
  sampler$storeState()
  attr(sampler$state, "cutPoints")
}
# the value of a call, or the message of the error it raised
outcomeOf <- function(expr) {
  tryCatch(expr, error = function(e) conditionMessage(e))
}
warmed <- function(...) {
  sampler <- dbarts(...)
  invisible(sampler$run(10L, 5L))
  sampler$storeState()
  sampler
}

set.seed(4127)
n <- 200L
x <- cbind(runif(n), rnorm(n))
y <- sin(4 * x[, 1L]) + x[, 2L] + rnorm(n, 0, 0.3)
xNew <- cbind(runif(n), rnorm(n))
yNew <- sin(4 * xNew[, 1L]) + xNew[, 2L] + rnorm(n, 0, 0.3)
newData <- dbartsData(xNew, yNew)
longer <- seq(0.01, 0.99, length.out = 50L)

for (useQuantiles in c(FALSE, TRUE)) {
  rule <- if (useQuantiles) "quantile" else "uniform"
  control <- dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 10L,
    n.cuts = 20L,
    useQuantiles = useQuantiles,
    updateState = FALSE,
    seed = 7L
  )

  # the undo leaves no residue: a longer grid set and put back, the state
  # restored, and the sampler derives at setData what an untouched twin, its
  # own copy and its reload derive, and draws what the twin draws
  sampler <- warmed(x, y, control = control)
  twin <- warmed(x, y, control = control)
  stored <- sampler$state
  original <- attr(stored, "cutPoints")
  expect_identical(lengths(original), c(20L, 20L), info = rule)
  sampler$setCutPoints(longer, 1L)
  expect_identical(lengths(cutPointsOf(sampler)), c(50L, 20L), info = rule)
  sampler$setCutPoints(original[[1L]], 1L)
  expect_true(sampler$setState(stored), info = rule)
  expect_true(twin$setState(twin$state), info = rule)
  path <- tempfile(fileext = ".rds")
  saveRDS(sampler, path)
  doors <- list(
    sampler = sampler,
    twin = twin,
    copy = sampler$copy(),
    reload = readRDS(path)
  )
  unlink(path)
  for (door in names(doors)) {
    doors[[door]]$setData(newData)
  }
  derived <- cutPointsOf(twin)
  expect_identical(lengths(derived), c(20L, 20L), info = rule)
  for (door in names(doors)) {
    expect_identical(
      cutPointsOf(doors[[door]]),
      derived,
      info = paste(rule, door)
    )
  }
  expect_identical(sampler$run(0L, 5L), twin$run(0L, 5L), info = rule)

  # a longer grid left in place does not reach the derivation either, and
  # the data object's count is what it was
  sampler <- warmed(x, y, control = control)
  sampler$setCutPoints(longer, 1L)
  sampler$setData(newData)
  expect_identical(cutPointsOf(sampler), derived, info = rule)
  expect_identical(sampler$data@n.cuts, c(20L, 20L), info = rule)

  # nor does one a state brought: a state stored on the longer grid installs
  # it, the grid is set back, and setData derives n.cuts points
  donor <- warmed(x, y, control = control)
  donor$setCutPoints(longer, 1L)
  donor$storeState()
  sampler <- warmed(x, y, control = control)
  sampler$setState(donor$state)
  expect_identical(lengths(cutPointsOf(sampler)), c(50L, 20L), info = rule)
  sampler$setCutPoints(original[[1L]], 1L)
  sampler$setData(newData)
  expect_identical(cutPointsOf(sampler), derived, info = rule)

  # the counts that stay: a refresh keeps the count the column holds, below
  # n.cuts and above it, and setData replaces a shorter grid by n.cuts points
  for (count in c(3L, 50L)) {
    info <- paste(rule, count)
    set <- seq(0.01, 0.99, length.out = count)
    sampler <- dbarts(x, y, control = control)
    sampler$setCutPoints(set, 1L)
    refreshed <- outcomeOf(
      sampler$setPredictor(xNew[, 1L], 1L, updateCutPoints = TRUE)
    )
    expect_identical(refreshed, TRUE, info = info)
    grid <- cutPointsOf(sampler)[[1L]]
    expect_identical(length(grid), count, info = info)
    expect_false(identical(grid, set), info = info)
  }
  sampler <- dbarts(x, y, control = control)
  sampler$setCutPoints(c(0.25, 0.5, 0.75), 1L)
  sampler$setData(newData)
  expect_identical(cutPointsOf(sampler), derived, info = rule)
}

# The uniform rule repeats a point over a column it cannot spread a grid
# over, a constant one or one narrower than doubles resolve; the sampler
# takes such a grid back, by column and as the whole list.
control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 10L,
  updateState = FALSE,
  seed = 7L
)
a <- runif(n)
const <- rep(1, n)
narrow <- 1 + (seq_len(n) %% 5L) * .Machine$double.eps
f <- factor(sample(letters[1L:4L], n, TRUE))
o <- factor(
  sample(c("lo", "mid", "hi"), n, TRUE),
  c("lo", "mid", "hi"),
  ordered = TRUE
)
y <- sin(4 * a) + (f == "b") + as.integer(o) / 2 + rnorm(n, 0, 0.3)

sampler <- warmed(cbind(a, const, narrow), y, control = control)
own <- attr(sampler$state, "cutPoints")
expect_identical(lengths(own), c(100L, 100L, 100L))
distinct <- lengths(lapply(own, unique))
expect_identical(distinct[1L:2L], c(100L, 1L))
expect_true(distinct[3L] > 1L && distinct[3L] < 100L)
expect_silent(sampler$setCutPoints(own[[2L]], 2L))
expect_silent(sampler$setCutPoints(own[[3L]], 3L))
expect_silent(sampler$setCutPoints(own))
expect_identical(cutPointsOf(sampler), own)

# a caller's grid with equal neighbours is taken too, and the sampler runs,
# stores and copies on it
expect_silent(sampler$setCutPoints(c(0.25, 0.5, 0.5, 0.75), 1L))
repeated <- cutPointsOf(sampler)
expect_identical(repeated[[1L]], c(0.25, 0.5, 0.5, 0.75))
expect_true(all(is.finite(sampler$run(0L, 5L)$train)))
expect_identical(cutPointsOf(sampler$copy()), repeated)
# a grid no sampler can hold is refused and leaves the grid as it was
for (cuts in list(c(0.6, 0.5), c(0.2, NaN), c(0.2, NA))) {
  expect_error(
    sampler$setCutPoints(cuts, 1L),
    pattern = "'cuts' must be sorted non-decreasingly and not contain NaN",
    fixed = TRUE
  )
}
expect_identical(cutPointsOf(sampler), repeated)

# The list a sampler reports has an entry per column, a factor's included.
# Given the whole list, the factor entries are skipped and not checked;
# naming a factor column stays refused, and so does a list of another length.
sampler <- warmed(y ~ a + f + o, data.frame(y, a, f, o), control = control)
own <- attr(sampler$state, "cutPoints")
expect_identical(lengths(own), c(100L, 0L, 2L))
expect_silent(sampler$setCutPoints(own))
expect_identical(cutPointsOf(sampler), own)
expect_silent(sampler$setCutPoints(list(c(0.25, 0.5), c(3, 1, NaN), NULL)))
expect_identical(cutPointsOf(sampler), c(list(c(0.25, 0.5)), own[-1L]))
expect_error(
  sampler$setCutPoints(0.5, 2L),
  pattern = "cannot set cut points for a categorical predictor"
)
expect_error(
  sampler$setCutPoints(0.5, "o"),
  pattern = "cannot set cut points for an ordered factor predictor"
)
expect_error(
  sampler$setCutPoints(own[-3L]),
  pattern = "requires one cut point vector per column"
)
expect_identical(cutPointsOf(sampler), c(list(c(0.25, 0.5)), own[-1L]))

# The whole undo on a design with a constant column and a factor: a grid is
# changed, the stored state's whole list handed back and the state restored,
# and the sampler continues as a twin that restored its own state.
design <- data.frame(y, a, const, f, o)
sampler <- warmed(y ~ a + const + f + o, design, control = control)
twin <- warmed(y ~ a + const + f + o, design, control = control)
stored <- sampler$state
sampler$setCutPoints(c(0.3, 0.6), 1L)
expect_identical(lengths(cutPointsOf(sampler)), c(2L, 100L, 0L, 2L))
expect_silent(sampler$setCutPoints(attr(stored, "cutPoints")))
expect_identical(cutPointsOf(sampler), attr(stored, "cutPoints"))
expect_true(sampler$setState(stored))
expect_true(twin$setState(twin$state))
expect_identical(sampler$run(0L, 5L), twin$run(0L, 5L))
