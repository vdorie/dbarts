# Whether a predictor column can hold a missing value is a property of the
# sampler: raised by the first one it is given, by any path, and never
# lowered. The first one draws the side of every rule already on the column,
# in the current trees and in the kept draws; the data object records the
# columns, so a copy and a reload hold them; predict refuses a missing value
# only in a column that has never held one. The draw itself - its coins,
# their order, and its undoing on a refusal - is held to the generator in
# tests/cpp (testMissingFirstSeen, testMonotoneMissingArrives).

set.seed(20261010L)
n <- 150L
complete <- data.frame(
  x1 = runif(n),
  x2 = runif(n),
  f = factor(sample(letters[1:4], n, replace = TRUE))
)
y <- 2 *
  (complete$x1 > 0.5) +
  complete$x2 +
  as.integer(complete$f) / 2 +
  rnorm(n, sd = 0.1)
controlOf <- function(n.chains = 2L, n.trees = 15L) {
  dbartsControl(
    n.chains = n.chains,
    n.threads = 1L,
    n.trees = n.trees,
    n.burn = 30L,
    n.samples = 20L,
    keepTrees = TRUE,
    updateState = FALSE,
    verbose = FALSE
  )
}
make <- function(data = complete, seed = 11L, control = controlOf()) {
  sampler <- dbarts(y ~ x1 + x2 + f, data, control = control, seed = seed)
  invisible(sampler$run())
  sampler
}
seen <- function(sampler) dbarts:::dataMissingSeen(sampler$data)
# the side each rule on x1 sends a missing value to, in the kept draws or
# the current trees; NULL while no column can hold one
sides <- function(sampler, current = FALSE, column = 1L) {
  trees <- sampler$getTrees(current = current)
  trees$missing[trees$var == column]
}
generators <- function(sampler) {
  sampler$storeState()
  lapply(sampler$state, `[[`, "rng.state")
}
gone <- c(3L, 40L, 77L)
holed <- complete
holed$x1[gone] <- NA
asked <- complete[1:2, ]
asked$x1[1L] <- NA
refusal <- "test predictors have missing values in 'x1', which carried none in training"

# ---- a column never missing: no record, and predict refuses (dec-B322) ----

sampler <- make()
expect_null(seen(sampler))
expect_null(sides(sampler))
expect_error(sampler$predict(asked), refusal, fixed = TRUE)

# ---- every path raises the column and draws, live and kept ----

# a numeric column and a factor, whose missing value is a missing label
for (column in c("x1", "f")) {
  j <- match(column, names(complete))
  brought <- complete
  brought[[column]][gone] <- NA
  paths <- list(
    "forced column" = function(s) {
      s$setPredictor(brought[[column]], column, forceUpdate = TRUE)
    },
    "unforced column" = function(s) s$setPredictor(brought[[column]], column),
    "whole frame" = function(s) s$setPredictor(brought, forceUpdate = FALSE),
    "row by row" = function(s) {
      s$setPredictor(brought[[column]], column, forceUpdate = "partial")
    },
    "jointly" = function(s) {
      updatePredictorPerObservationJointly(list(s), brought[[column]], column)
    },
    "setData" = function(s) s$setData(dbartsData(y ~ x1 + x2 + f, brought))
  )
  for (path in names(paths)) {
    info <- paste(column, path)
    sampler <- make()
    paths[[path]](sampler)
    expect_identical(seen(sampler), seq_len(3L) == j, info = info)
    # the kept draws were all recorded before the value arrived
    expect_true(all(c("L", "R") %in% sides(sampler, FALSE, j)), info = info)
    expect_true(all(c("L", "R") %in% sides(sampler, TRUE, j)), info = info)
    # predict answers for a missing value from every kept draw of each chain
    fits <- sampler$predict(brought[gone[1:2], ])
    expect_identical(dim(fits), c(2L, 20L, 2L), info = info)
    expect_true(all(is.finite(fits)), info = info)
    expect_true(all(is.finite(sampler$run(0L, 2L)$train)), info = info)
    # and the sampler installs the state it then stores
    sampler$storeState()
    expect_true(sampler$setState(sampler$state), info = info)
  }
}

# ---- two forests: the rules of each draw, live and kept ----

# The sides of the rules on x1 by forest, kept draws and then current trees;
# every rule sends a missing value left while no column can hold one.
z <- as.double(seq_len(n) %% 2L)
sidesByForest <- function(sampler) {
  unlist(
    lapply(c(FALSE, TRUE), function(current) {
      trees <- sampler$getTrees(current = current)
      onX1 <- trees$var == 1L
      side <- if (is.null(trees$missing)) "L" else trees$missing[onX1]
      split(rep_len(side, sum(onX1)), factor(trees$forest[onX1], 1:2))
    }),
    recursive = FALSE
  )
}
# the first seed whose forests each hold a rule on x1, kept and current,
# before the value arrives, and rules on both sides after it
found <- NA_integer_
for (seed in 11:40) {
  sampler <- dbarts(
    y ~ x1 + x2 + f,
    complete,
    forests = list(forest(), forest(basis = z)),
    control = controlOf(),
    seed = seed
  )
  invisible(sampler$run())
  before <- sidesByForest(sampler)
  sampler$setPredictor(holed$x1, "x1", forceUpdate = TRUE)
  after <- sidesByForest(sampler)
  held <- all(lengths(before) > 0L) && all(unlist(before) == "L")
  both <- vapply(after, function(side) all(c("L", "R") %in% side), NA)
  if (held && all(both)) {
    found <- seed
    break
  }
}
expect_false(is.na(found))
expect_identical(length(after), 4L)

# ---- the sides are fair coins: 40 seeds, one chain ----

right <- total <- c(kept = 0L, current = 0L)
for (seed in 1:40) {
  sampler <- make(seed = seed, control = controlOf(1L, 10L))
  sampler$setPredictor(holed$x1, "x1", forceUpdate = TRUE)
  drawn <- list(kept = sides(sampler), current = sides(sampler, TRUE))
  right <- right + vapply(drawn, function(side) sum(side == "R"), 0L)
  total <- total + lengths(drawn)
}
expect_true(all(total >= c(1000L, 40L)))
# a binomial band of 1e-3 around one half
expect_true(all(abs(right - total / 2) < qnorm(1 - 5e-4) * sqrt(total) / 2))

# ---- a second missing value draws nothing; filling keeps the sides ----

sampler <- make()
sampler$setPredictor(holed$x1, "x1", forceUpdate = TRUE)
drawn <- list(sides(sampler), generators(sampler))
more <- holed$x1
more[c(5L, 90L)] <- NA
expect_true(sampler$setPredictor(more, "x1"))
expect_identical(list(sides(sampler), generators(sampler)), drawn)
sampler$setPredictor(complete$x1, "x1", forceUpdate = TRUE)
expect_false(anyNA(as.matrix(sampler$data@x)))
expect_identical(seen(sampler), c(TRUE, FALSE, FALSE))
expect_identical(list(sides(sampler), generators(sampler)), drawn)
# predict's answer for the missing value, or a failed expectation and NULL
answered <- function(sampler) {
  fits <- NULL
  expect_silent(fits <- sampler$predict(asked))
  fits
}
expect_true(all(is.finite(answered(sampler))))
# setData to complete data keeps them too, and takes no record from its
# argument: the object's stale slot names x2, which stays refused
stale <- dbartsData(y ~ x1 + x2 + f, complete)
stale@missing.seen <- c(FALSE, TRUE, FALSE)
sampler$setData(stale)
expect_identical(seen(sampler), c(TRUE, FALSE, FALSE))
expect_identical(sides(sampler), drawn[[1L]])
expect_true(all(is.finite(answered(sampler))))
askedX2 <- complete[1:2, ]
askedX2$x2[1L] <- NA
expect_error(sampler$predict(askedX2), "missing values in 'x2'", fixed = TRUE)
never <- make()
never$setData(stale)
expect_null(seen(never))
expect_error(never$predict(askedX2), "missing values in 'x2'", fixed = TRUE)

# ---- a copy and a reload hold the columns, and continue alike ----

sampler$storeState()
duplicate <- sampler$copy()
path <- tempfile(fileext = ".rds")
saveRDS(sampler, path)
reloaded <- readRDS(path)
unlink(path)
for (route in list(duplicate, reloaded)) {
  expect_identical(seen(route), c(TRUE, FALSE, FALSE))
  expect_identical(answered(route), answered(sampler))
  # the engine holds the flag, not the record alone: a change writes it back
  route$setPredictor(complete$x1, "x1", forceUpdate = TRUE)
  expect_identical(seen(route), c(TRUE, FALSE, FALSE))
}
expect_identical(duplicate$run(0L, 5L), reloaded$run(0L, 5L))
# a sampler over a data handle, as xbart's folds are, takes the record with
# its data object, here complete: the rules on a recorded column draw their
# sides, and a view of some columns reads the record by the handle's
viewTrees <- function(columns = NULL) {
  view <- dbarts:::bartcoreSamplerFromHandle(
    dbarts:::bartcoreDataHandle(sampler$control, sampler$data),
    sampler$control,
    sampler$model,
    sampler$data,
    seq_len(n),
    columns = columns
  )
  invisible(dbarts:::bartcoreRun(view, 30L, 1L))
  getTrees <- dbarts:::C_dbarts_bartcore_getTrees
  .Call(getTrees, view$ptr, 1:2, NULL, 1:15, TRUE, NULL, NULL, 0L)
}
trees <- viewTrees()
expect_true(all(0:1 %in% trees$missing[trees$var == 1L]))
trees <- viewTrees(2:1)
expect_true(all(0:1 %in% trees$missing[trees$var == 2L]))
# a data object saved before the slot existed reads as no record, and a
# sampler created from it holds the columns its values name alone
old <- sampler$data
attr(old, "missing.seen") <- NULL
expect_false(methods::.hasSlot(old, "missing.seen"))
expect_null(dbarts:::dataMissingSeen(old))
expect_error(
  dbarts(old, control = controlOf(), seed = 3L)$predict(asked),
  refusal,
  fixed = TRUE
)
old@missing.seen <- NA
expect_error(validObject(old), "'missing.seen' must be NULL or a logical")

# ---- a refused first missing value leaves the sampler its twin ----

emptying <- holed
emptying$x2 <- 0.5
sampler <- make()
twin <- make()
expect_false(
  sampler$setPredictor(
    as.matrix(emptying[1:2]),
    c("x1", "x2"),
    forceUpdate = FALSE
  )
)
sampler$storeState()
twin$storeState()
expect_identical(sampler$state, twin$state)
expect_null(seen(sampler))
expect_identical(sampler$run(0L, 20L), twin$run(0L, 20L))
expect_false(sampler$setPredictor(emptying, forceUpdate = FALSE))
sampler$storeState()
twin$storeState()
expect_identical(sampler$state, twin$state)
expect_null(seen(sampler))
expect_equal(
  sampler$run(0L, 20L)$train,
  twin$run(0L, 20L)$train,
  tolerance = 1e-12
)

# ---- a rule sending every level one way (dec-B378) ----

# One tree, split once on f: no level right, a missing value right. On a
# sampler whose f has never held one the rule is malformed. Where f holds
# one it installs as stored; where f has held one and is complete again it
# is built, declined for its empty side and, forced, installed with that
# side merged.
oneWay <- function(sampler) {
  sampler$storeState()
  state <- sampler$state
  forest <- state[[1L]]$forests[[1L]]
  forest$tree.vars <- c(3L, -1L, -1L)
  forest$tree.values <- writeBin(c(0, -0.1, 0.1), raw())
  forest$tree.sizes <- 3L
  forest$tree.flags <- as.raw(c(5L, 0L, 0L))
  state[[1L]]$forests[[1L]] <- forest
  state
}
single <- controlOf(1L, 1L)
single@keepTrees <- FALSE
never <- dbarts(y ~ x1 + x2 + f, complete, control = single, seed = 5L)
expect_error(
  never$setState(oneWay(never)),
  "state is not consistent with this sampler",
  fixed = TRUE
)
holedF <- complete
holedF$f[gone] <- NA
held <- dbarts(y ~ x1 + x2 + f, holedF, control = single, seed = 5L)
expect_true(held$setState(oneWay(held)))
expect_identical(held$getTrees()$missing, c("R", NA, NA))
held$setPredictor(complete$f, "f", forceUpdate = TRUE)
expect_identical(seen(held), c(FALSE, FALSE, TRUE))
refilled <- oneWay(held)
expect_identical(held$setState(refilled), FALSE)
expect_null(held$setState(refilled, forceUpdate = TRUE))
expect_identical(held$state, refilled)
expect_identical(nrow(held$getTrees()), 1L)

# ---- row by row, a first missing value that would empty a leaf ----

# One tree, cut on x1 below every value but the first row's, which is alone
# in the left leaf. Where the draw sends a missing value right the row is
# declined and the draw taken back: no record, no side, and the generator and
# next draws of a twin given the column it holds. Where it sends it left the
# row stays where it is and is taken with its side. Seeds are tried until
# each has been met.
lone <- complete
lone$x1[1L] <- -5
loneLeft <- function(seed) {
  sampler <- dbarts(y ~ x1 + x2 + f, lone, control = single, seed = seed)
  sampler$storeState()
  state <- sampler$state
  forest <- state[[1L]]$forests[[1L]]
  forest$tree.vars <- c(1L, -1L, -1L)
  cut <- attr(state, "cutPoints")[[1L]][1L]
  forest$tree.values <- writeBin(c(cut, -0.1, 0.1), raw())
  forest$tree.sizes <- 3L
  forest$tree.flags <- as.raw(c(2L, 0L, 0L))
  state[[1L]]$forests[[1L]] <- forest
  stopifnot(isTRUE(sampler$setState(state)))
  sampler
}
brought <- replace(lone$x1, 1L, NA)
met <- c(declined = FALSE, taken = FALSE)
for (seed in 1:40) {
  sampler <- loneLeft(seed)
  installed <- sampler$setPredictor(brought, "x1", forceUpdate = "partial")
  arm <- if (installed[1L]) "taken" else "declined"
  if (met[[arm]]) {
    next
  }
  met[[arm]] <- TRUE
  expect_true(all(installed[-1L]), info = arm)
  if (arm == "taken") {
    expect_identical(seen(sampler), c(TRUE, FALSE, FALSE))
    expect_identical(sides(sampler, TRUE), "L")
  } else {
    twin <- loneLeft(seed)
    expect_true(all(twin$setPredictor(lone$x1, "x1", forceUpdate = "partial")))
    expect_null(seen(sampler))
    expect_null(sides(sampler, TRUE))
    expect_identical(generators(sampler), generators(twin))
    expect_identical(sampler$run(0L, 20L), twin$run(0L, 20L))
  }
  if (all(met)) {
    break
  }
}
expect_identical(met, c(declined = TRUE, taken = TRUE))
