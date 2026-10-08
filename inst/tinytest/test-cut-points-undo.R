# A cut grid changed with setCutPoints can be put back on any design, and
# leaves nothing behind: a later setData derives the grid n.cuts names
# whatever grid was set before it, as a copy and a reload do. A caller's own
# grid is sorted and holds each point once, as every grid a sampler derives
# does, so the grids a sampler reports go back into it whole. The grid in
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
  controlWith <- function(n.trees) {
    dbartsControl(
      n.chains = 1L,
      n.threads = 1L,
      n.trees = n.trees,
      n.cuts = 20L,
      useQuantiles = useQuantiles,
      updateState = FALSE,
      seed = 7L
    )
  }
  control <- controlWith(10L)

  # the undo leaves no residue: a longer grid set and put back, the state
  # restored, and the sampler derives at setData what a twin that restored
  # its own state, its own copy and its reload derive, and draws what the
  # twin draws
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

  # no count stays: a refresh derives the grid setData derives for the
  # values, over a set grid below n.cuts and above it, and setData replaces a
  # shorter grid by n.cuts points
  for (count in c(3L, 50L)) {
    info <- paste(rule, count)
    set <- seq(0.01, 0.99, length.out = count)
    sampler <- dbarts(x, y, control = control)
    sampler$setCutPoints(set, 1L)
    refreshed <- outcomeOf(
      sampler$setPredictor(
        xNew[, 1L],
        1L,
        forceUpdate = TRUE,
        updateCutPoints = "position"
      )
    )
    expect_null(refreshed, info = info)
    expect_identical(cutPointsOf(sampler)[[1L]], derived[[1L]], info = info)
  }
  sampler <- dbarts(x, y, control = control)
  sampler$setCutPoints(c(0.25, 0.5, 0.75), 1L)
  sampler$setData(newData)
  expect_identical(cutPointsOf(sampler), derived, info = rule)

  # states that come and go leave a longer grid as it was: one on a shorter
  # grid installed and the longer grid set again, then one refused, and the
  # column still refreshes to the grid a derivation gives it
  sampler <- warmed(x, y, control = control)
  sampler$setCutPoints(longer, 1L)
  donor <- warmed(x, y, control = control)
  expect_true(sampler$setState(donor$state), info = rule)
  expect_identical(lengths(cutPointsOf(sampler)), c(20L, 20L), info = rule)
  sampler$setCutPoints(longer, 1L)
  donor <- warmed(x, y, control = controlWith(7L))
  expect_error(
    sampler$setState(donor$state),
    pattern = "state is not consistent with this sampler",
    info = rule
  )
  refreshed <- outcomeOf(
    sampler$setPredictor(
      xNew[, 1L],
      1L,
      forceUpdate = TRUE,
      updateCutPoints = "position"
    )
  )
  expect_false(is.character(refreshed), info = rule)
  expect_identical(cutPointsOf(sampler)[[1L]], derived[[1L]], info = rule)

  # a sparse column counts its refresh for itself, from n.cuts as a dense one
  # does
  if (requireNamespace("Matrix", quietly = TRUE)) {
    sparseColumn <- function() {
      Matrix::Matrix(
        cbind(s = ifelse(runif(n) < 0.5, 0, runif(n))),
        sparse = TRUE
      )
    }
    # the residual scale is stated so none is estimated on a sparse design
    sampler <- dbarts(
      cbind(
        Matrix::Matrix(x[, 1L, drop = FALSE], sparse = TRUE),
        sparseColumn()
      ),
      y,
      control = control,
      sigest = sd(y)
    )
    for (count in c(50L, 3L)) {
      info <- paste(rule, "sparse", count)
      set <- seq(0.01, 0.99, length.out = count)
      sampler$setCutPoints(set, 2L)
      refreshed <- outcomeOf(
        sampler$setPredictor(
          sparseColumn(),
          2L,
          forceUpdate = TRUE,
          updateCutPoints = "position"
        )
      )
      expect_null(refreshed, info = info)
      grid <- cutPointsOf(sampler)[[2L]]
      expect_identical(length(grid), 20L, info = info)
      expect_false(is.unsorted(grid, strictly = TRUE), info = info)
    }
  }
}

# The uniform rule holds one point over a constant column and five over one
# narrower than doubles resolve, where it once repeated them; the sampler
# takes the grid such a column holds back, by column and as the whole list.
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
stored <- sampler$state
own <- attr(stored, "cutPoints")
expect_identical(lengths(own), c(100L, 1L, 5L))
expect_identical(lengths(lapply(own, unique)), c(100L, 1L, 5L))
expect_identical(own[[2L]], 1)
expect_silent(sampler$setCutPoints(own[[2L]], 2L))
expect_silent(sampler$setCutPoints(own[[3L]], 3L))
expect_silent(sampler$setCutPoints(own))
expect_identical(cutPointsOf(sampler), own)

# A grid with equal neighbours is refused on any column, and so is a missing
# value, by column and as an entry of the whole list. The grid is left as it
# was.
refusal <- paste(
  "$setCutPoints: a cut point may appear only once in 'cuts'; for more",
  "splits near a value, give a denser grid around it"
)
nanRefusal <- "$setCutPoints: 'cuts' must not contain NA or NaN"
notHeld <- list(c(0.25, 0.5, 0.5, 0.75), c(0.6, 0.5, 0.6))
for (cuts in notHeld) {
  expect_error(sampler$setCutPoints(cuts, 1L), pattern = refusal, fixed = TRUE)
  expect_error(
    sampler$setCutPoints(list(cuts, own[[2L]], own[[3L]])),
    pattern = refusal,
    fixed = TRUE
  )
}
for (cuts in list(c(0.2, NaN), c(0.2, NA), c(0.6, NaN, 0.5))) {
  expect_error(
    sampler$setCutPoints(cuts, 1L),
    pattern = nanRefusal,
    fixed = TRUE
  )
  expect_error(
    sampler$setCutPoints(list(cuts, own[[2L]], own[[3L]])),
    pattern = nanRefusal,
    fixed = TRUE
  )
}
for (cuts in list(
  rep(2, 100L),
  rep(own[[2L]], 2L),
  rep(own[[3L]], each = 2L)
)) {
  expect_error(sampler$setCutPoints(cuts, 2L), pattern = refusal, fixed = TRUE)
}
expect_identical(cutPointsOf(sampler), own)
# over a column of zeros the grid is one zero; a zero beside a negative zero
# is one point given twice
zeros <- dbarts(cbind(a, zero = 0), y, control = control)
held <- cutPointsOf(zeros)[[2L]]
expect_identical(held, 0)
expect_silent(zeros$setCutPoints(held, 2L))
expect_error(zeros$setCutPoints(c(-0, 0), 2L), pattern = refusal, fixed = TRUE)

# A grid out of order is sorted and taken, by column, as an entry of the whole
# list and as a data frame: the sampler holds the sorted grid and draws what
# one handed the sorted grid draws. The two are twins and not a copy and its
# original: a copy holds its leaves' rows in another order and draws its
# original's draws to rounding only.
for (grid in list(c(0.75, 0.25, 0.5), c(0.5, 0.75, 0.25))) {
  sorted <- warmed(cbind(a, const, narrow), y, control = control)
  shuffled <- warmed(cbind(a, const, narrow), y, control = control)
  sorted$setCutPoints(c(0.25, 0.5, 0.75), 1L)
  shuffled$setCutPoints(grid, 1L)
  expect_identical(cutPointsOf(shuffled), cutPointsOf(sorted))
  expect_identical(shuffled$run(0L, 5L)$train, sorted$run(0L, 5L)$train)
}
# two changes in turn
sorted <- warmed(cbind(a, const, narrow), y, control = control)
shuffled <- warmed(cbind(a, const, narrow), y, control = control)
sorted$setCutPoints(c(0.25, 0.5, 0.75), 1L)
shuffled$setCutPoints(list(c(0.75, 0.25, 0.5)), 1L)
expect_identical(cutPointsOf(shuffled), cutPointsOf(sorted))
sorted$setCutPoints(own)
shuffled$setCutPoints(list(rev(own[[1L]]), own[[2L]], rev(own[[3L]])))
expect_identical(cutPointsOf(shuffled), own)
expect_identical(shuffled$run(0L, 5L)$train, sorted$run(0L, 5L)$train)
# a data frame is the list of its columns
sorted <- warmed(cbind(a, b = a^2), y, control = control)
shuffled <- warmed(cbind(a, b = a^2), y, control = control)
both <- cutPointsOf(sorted)
sorted$setCutPoints(both)
shuffled$setCutPoints(data.frame(rev(both[[1L]]), rev(both[[2L]])))
expect_identical(cutPointsOf(shuffled), both)
expect_identical(shuffled$run(0L, 5L)$train, sorted$run(0L, 5L)$train)
rm(grid, sorted, shuffled, both)

# A grid of distinct points is taken on either kind of column, the grids the
# sampler first held among them; with each point given twice they are
# refused, and the stored state brings grids and trees back.
expect_silent(sampler$setCutPoints(c(0.25, 0.5, 0.75), 1L))
expect_silent(sampler$setCutPoints(c(0.5, 1.5), 2L))
expect_identical(
  cutPointsOf(sampler),
  list(c(0.25, 0.5, 0.75), c(0.5, 1.5), own[[3L]])
)
expect_error(
  sampler$setCutPoints(lapply(own, rep, each = 2L)),
  pattern = refusal,
  fixed = TRUE
)
expect_silent(sampler$setCutPoints(own))
expect_identical(cutPointsOf(sampler), own)
expect_true(sampler$setState(stored))
expect_identical(cutPointsOf(sampler), own)

# A grid is a vector of numbers. What as.double would turn into one is
# refused by name with the rest, by column and as an entry of the whole
# list: a factor would go in as its codes and a Date as its day count.
notNumeric <- "$setCutPoints: 'cuts' must be numeric"
notGrids <- list(
  mean,
  c("0.2", "0.4"),
  c(FALSE, TRUE),
  factor(c("u", "v")),
  as.Date(c("2020-01-01", "2020-01-02")),
  NULL
)
for (cuts in notGrids) {
  expect_error(
    sampler$setCutPoints(list(cuts), 1L),
    pattern = notNumeric,
    fixed = TRUE
  )
  expect_error(
    sampler$setCutPoints(c(list(cuts), own[-1L])),
    pattern = notNumeric,
    fixed = TRUE
  )
}
expect_identical(cutPointsOf(sampler), own)
# whole numbers are numbers
expect_silent(sampler$setCutPoints(c(0L, 1L), 1L))
expect_identical(cutPointsOf(sampler)[[1L]], c(0, 1))

# The list a sampler reports has an entry per column, a factor's included.
# Given the whole list, what sits in a factor column's place is not read: it
# neither warns nor fails, whatever it is. The numeric entries are held to
# the rule; naming a factor column stays refused, and so does a list of
# another length.
sampler <- warmed(y ~ a + f + o, data.frame(y, a, f, o), control = control)
own <- attr(sampler$state, "cutPoints")
expect_identical(lengths(own), c(100L, 0L, 2L))
expect_silent(sampler$setCutPoints(own))
expect_identical(cutPointsOf(sampler), own)
expect_silent(sampler$setCutPoints(list(c(0.25, 0.5), levels(f), mean)))
set <- c(list(c(0.25, 0.5)), own[-1L])
expect_identical(cutPointsOf(sampler), set)
for (cuts in notHeld) {
  expect_error(
    sampler$setCutPoints(list(cuts, NULL, NULL)),
    pattern = refusal,
    fixed = TRUE
  )
}
expect_error(
  sampler$setCutPoints(list(mean, NULL, NULL)),
  pattern = notNumeric,
  fixed = TRUE
)
expect_error(
  sampler$setCutPoints(0.5, 2L),
  pattern = "cannot set cut points for a categorical predictor"
)
expect_error(
  sampler$setCutPoints(0.5, "o"),
  pattern = "cannot set cut points for an ordered factor predictor"
)
# a list of another length is refused, shorter or longer, before any entry
# is read: nothing in it is coerced and nothing warns
wrongLength <- "requires one cut point vector per column"
warnings <- character()
withCallingHandlers(
  for (cuts in list(own[-3L], c(own, list(1)), list(levels(f), mean))) {
    expect_error(sampler$setCutPoints(cuts), pattern = wrongLength)
  },
  warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
expect_identical(warnings, character())
expect_identical(cutPointsOf(sampler), set)
# entries are dropped by position only when no column is named: a column
# named three times on this design of three columns has each entry read
expect_silent(sampler$setCutPoints(list(0.2, 0.4, 0.6), c(1L, 1L, 1L)))
expect_identical(cutPointsOf(sampler), c(list(0.6), own[-1L]))
expect_error(
  sampler$setCutPoints(list(0.2, "0.4", 0.6), c(1L, 1L, 1L)),
  pattern = notNumeric,
  fixed = TRUE
)
# a data frame is a list of columns, and is taken as that list is
expect_silent(
  sampler$setCutPoints(data.frame(a = c(0.3, 0.6), f = c("x", "y"), o = 1:2))
)
expect_identical(cutPointsOf(sampler), c(list(c(0.3, 0.6)), own[-1L]))

# a design of factor columns alone leaves a whole list nothing to install:
# the call returns and the sampler is as it was
sampler <- warmed(y ~ f + o, data.frame(y, f, o), control = control)
stored <- sampler$state
expect_identical(lengths(attr(stored, "cutPoints")), c(0L, 2L))
expect_silent(sampler$setCutPoints(list(levels(f), mean)))
sampler$storeState()
expect_identical(sampler$state, stored)

# The whole undo on a design with a constant column and a factor: a grid is
# changed and the stored state's whole list handed back, which puts every
# grid back before the state is restored. The restore and the continuation
# after it are pins that hold without the hand-back too, a state bringing its
# own grid: the sampler continues as a twin that restored its own state.
design <- data.frame(y, a, const, f, o)
sampler <- warmed(y ~ a + const + f + o, design, control = control)
twin <- warmed(y ~ a + const + f + o, design, control = control)
stored <- sampler$state
sampler$setCutPoints(c(0.3, 0.6), 1L)
expect_identical(lengths(cutPointsOf(sampler)), c(2L, 1L, 0L, 2L))
expect_silent(sampler$setCutPoints(attr(stored, "cutPoints")))
expect_identical(cutPointsOf(sampler), attr(stored, "cutPoints"))
expect_true(sampler$setState(stored))
expect_true(twin$setState(twin$state))
expect_identical(sampler$run(0L, 5L), twin$run(0L, 5L))
