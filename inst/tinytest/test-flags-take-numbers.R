# A TRUE/FALSE argument takes a single logical or number, 0 being FALSE, as
# as.logical reads it; NA, a string and another length are refused by name.
# updateState also takes NULL, and its NA keeps its warning (test-na-as-none.R).

set.seed(41L)
n <- 60L
x <- matrix(rnorm(n * 2L), n)
y <- x[, 1L] + rnorm(n, sd = 0.3)
w01 <- rep_len(c(1, 0, 1), n)
runsClean <- function(expr) {
  expect_identical(
    tryCatch(
      {
        expr
        "ok"
      },
      error = conditionMessage
    ),
    "ok"
  )
}
make <- function() {
  sampler <- dbarts::dbarts(
    x,
    y,
    n.trees = 10L,
    control = dbarts::dbartsControl(
      n.chains = 1L,
      n.threads = 1L,
      n.samples = 4L,
      keepTrees = TRUE,
      updateState = FALSE,
      verbose = FALSE
    )
  )
  invisible(sampler$run(5L, 4L))
  sampler$storeState()
  sampler
}
refused <- function(call, name, ...) {
  for (bad in list(NA, "TRUE", "yes", c(TRUE, FALSE), c(1, 0))) {
    expect_error(call(bad), paste0("'", name, "' must be TRUE"), fixed = TRUE)
  }
}

# updateState: a number is taken; NULL defers to the control
s <- make()
for (good in list(1L, 0L, 2, 0, TRUE, FALSE, NULL)) {
  runsClean(s$setSigma(2, updateState = good))
  runsClean(s$run(0L, 1L, updateState = good))
}
for (bad in list("TRUE", "yes", c(TRUE, FALSE), c(1, 0))) {
  expect_error(
    s$setSigma(2, updateState = bad),
    "'updateState' must be TRUE, FALSE or NULL",
    fixed = TRUE
  )
  expect_error(
    s$run(0L, 1L, updateState = bad),
    "'updateState' must be TRUE, FALSE or NULL",
    fixed = TRUE
  )
}
# the number is read as logical: 1 stores the state, 0 leaves it
s <- make()
before <- s$state
s$setSigma(2, updateState = 0)
expect_identical(s$state, before)
s$setSigma(2, updateState = 1)
expect_false(identical(s$state, before))

# setState's and installTrees' forceUpdate
s <- make()
donor <- make()
for (good in list(1, 0, TRUE, FALSE)) {
  runsClean(s$setState(s$state, forceUpdate = good))
  runsClean(s$installTrees(donor, forceUpdate = good))
}
expect_null(s$setState(s$state, forceUpdate = 1))
expect_true(s$setState(s$state, forceUpdate = 0))
refused(function(bad) s$setState(s$state, forceUpdate = bad), "forceUpdate")
refused(function(bad) s$installTrees(donor, forceUpdate = bad), "forceUpdate")
expect_null(s$setPredictor(x, forceUpdate = 1))
expect_true(s$setPredictor(x, forceUpdate = 0))
for (bad in list(NA, "yes", c(TRUE, FALSE))) {
  expect_error(
    s$setPredictor(x, forceUpdate = bad),
    "'forceUpdate' must be TRUE, FALSE or \"partial\"",
    fixed = TRUE
  )
}

# the response conduits' updateScale and updateStandardization
s <- make()
for (good in list(1, 0, TRUE, FALSE)) {
  runsClean(s$setResponse(y, updateScale = good))
  runsClean(s$setOffset(rep(0, n), updateScale = good))
  runsClean(s$setWeights(
    rep(1, n),
    updateScale = FALSE,
    updateStandardization = good
  ))
  runsClean(s$setData(s$data, updateScale = good, updateStandardization = good))
}
# a number rescales as TRUE does, and 0 as FALSE does
a <- make()
b <- make()
a$setResponse(y * 5, updateScale = 1)
b$setResponse(y * 5, updateScale = TRUE)
expect_identical(a$data@response.range, b$data@response.range)
refused(function(bad) s$setResponse(y, updateScale = bad), "updateScale")
refused(function(bad) s$setOffset(rep(0, n), updateScale = bad), "updateScale")
refused(function(bad) s$setWeights(w01, updateScale = bad), "updateScale")
refused(function(bad) s$setActiveRows(w01, updateScale = bad), "updateScale")
refused(
  function(bad) s$setPredictor(x, updateStandardization = bad),
  "updateStandardization"
)
refused(
  function(bad) s$setWeights(w01, updateStandardization = bad),
  "updateStandardization"
)
for (bad in list(NA, "TRUE", c(TRUE, FALSE))) {
  expect_error(
    s$setData(s$data, updateScale = bad),
    "'updateScale' must be TRUE, FALSE or NULL",
    fixed = TRUE
  )
}

# verbose at the creation doors
runsClean(dbarts::dbarts(x, y, n.trees = 5L, verbose = 0))
runsClean(dbarts::xbart(
  x,
  y,
  n.trees = 5L,
  n.reps = 1L,
  n.threads = 1L,
  k = 1,
  verbose = 0
))
for (bad in list(NA, "TRUE", c(TRUE, FALSE), NULL)) {
  expect_error(
    dbarts::dbarts(x, y, n.trees = 5L, verbose = bad),
    "'verbose' argument to dbarts must be TRUE or FALSE",
    fixed = TRUE
  )
}

# xbart's useQuantiles and drop
for (good in list(1, 0)) {
  runsClean(dbarts::xbart(
    x,
    y,
    n.trees = 5L,
    n.reps = 1L,
    n.threads = 1L,
    k = 1,
    useQuantiles = good
  ))
  runsClean(dbarts::xbart(
    x,
    y,
    n.trees = 5L,
    n.reps = 1L,
    n.threads = 1L,
    k = 1,
    drop = good
  ))
}
for (bad in list(NA, "TRUE", c(TRUE, FALSE))) {
  expect_error(
    dbarts::xbart(
      x,
      y,
      n.trees = 5L,
      n.reps = 1L,
      n.threads = 1L,
      k = 1,
      useQuantiles = bad
    ),
    "'useQuantiles' must be",
    fixed = TRUE
  )
  expect_error(
    dbarts::xbart(
      x,
      y,
      n.trees = 5L,
      n.reps = 1L,
      n.threads = 1L,
      k = 1,
      drop = bad
    ),
    "'drop' must be TRUE or FALSE",
    fixed = TRUE
  )
}

# DART's update.alpha and the partial-dependence doors' flags
runsClean(dbarts::dbartsPriors$dart(update.alpha = 1))
runsClean(dbarts::dbartsPriors$dart(update.alpha = 0))
for (bad in list(NA, "TRUE", c(TRUE, FALSE))) {
  expect_error(
    dbarts::dbartsPriors$dart(update.alpha = bad),
    "'update.alpha' must be TRUE or FALSE",
    fixed = TRUE
  )
  expect_error(
    dbarts:::pdbart.flag(bad, "keepSampler"),
    "'keepSampler' must be TRUE or FALSE",
    fixed = TRUE
  )
}
expect_identical(dbarts:::pdbart.flag(1, "keepSampler"), TRUE)
expect_identical(dbarts:::pdbart.flag(0L, "keepSampler"), FALSE)

# the flags of bartBT, rbart_vi, plot, the sampler's copy and getTrees, and
# the helper functions
flagBad <- list(NA, "yes", c(TRUE, FALSE))
smallBT <- function(...) {
  suppressWarnings(dbarts::bartBT(
    x,
    y,
    ndpost = 4L,
    nskip = 2L,
    ntree = 5L,
    verbose = FALSE,
    ...
  ))
}
for (arg in c(
  "combinechains",
  "keeptrees",
  "keepcall",
  "sampleronly",
  "keepsampler"
)) {
  for (bad in flagBad) {
    expect_error(
      do.call(smallBT, stats::setNames(list(bad), arg)),
      paste0("'", arg, "' must be TRUE or FALSE"),
      fixed = TRUE
    )
  }
}
runsClean(smallBT(combinechains = 1, keepcall = 0))

dd <- data.frame(y = y, x1 = x[, 1L], g = gl(5L, n / 5L))
rbartFlag <- function(...) {
  suppressWarnings(dbarts::rbart_vi(
    y ~ x1,
    data = dd,
    group.by = g,
    n.samples = 4L,
    n.burn = 2L,
    n.thin = 1L,
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 5L,
    verbose = FALSE,
    ...
  ))
}
for (arg in c("combineChains", "keepCall", "keepTestFits")) {
  for (bad in flagBad) {
    expect_error(
      do.call(rbartFlag, stats::setNames(list(bad), arg)),
      paste0("'", arg, "' must be TRUE or FALSE"),
      fixed = TRUE
    )
  }
}

sampler <- make()
for (bad in flagBad) {
  expect_error(
    sampler$copy(shallow = bad),
    "'shallow' must be TRUE or FALSE",
    fixed = TRUE
  )
  expect_error(
    sampler$getTrees(current = bad),
    "'current' must be TRUE or FALSE",
    fixed = TRUE
  )
  expect_error(
    sampler$getVariance(test = bad),
    "'test' must be TRUE or FALSE",
    fixed = TRUE
  )
  expect_error(
    dbarts::guessNumCores(bad),
    "'logical' must be TRUE or FALSE",
    fixed = TRUE
  )
  expect_error(
    dbarts::makeModelMatrixFromDataFrame(data.frame(a = gl(2L, 2L)), bad),
    "'drop' must be TRUE or FALSE",
    fixed = TRUE
  )
  expect_error(
    dbarts::makeind(data.frame(a = gl(2L, 2L)), bad),
    "'all' must be TRUE or FALSE",
    fixed = TRUE
  )
}
runsClean(sampler$copy(shallow = 1))
runsClean(sampler$getTrees(current = 1))
expect_identical(
  dim(dbarts::makeModelMatrixFromDataFrame(data.frame(a = gl(2L, 2L)), 1)),
  dim(dbarts::makeModelMatrixFromDataFrame(data.frame(a = gl(2L, 2L)), TRUE))
)
expect_identical(
  dim(dbarts::makeModelMatrixFromDataFrame(data.frame(a = gl(2L, 2L)), 0)),
  dim(dbarts::makeModelMatrixFromDataFrame(data.frame(a = gl(2L, 2L)), FALSE))
)
fitBT <- dbarts::bart(
  x,
  y,
  n.samples = 4L,
  n.burn = 2L,
  n.trees = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
for (bad in flagBad) {
  expect_error(
    extract(fitBT, type = "forest", contribution = bad),
    "'contribution' must be TRUE or FALSE",
    fixed = TRUE
  )
}
for (bad in flagBad) {
  expect_error(
    dbarts:::plot.pd2bart(list(), justmedian = bad),
    "'justmedian' must be TRUE or FALSE",
    fixed = TRUE
  )
}
