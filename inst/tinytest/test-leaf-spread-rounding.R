# A linear or gp leaf covariate whose spread about its center, over the rows
# in the likelihood, is under 64 times machine precision of its size
# (||z - center|| < 64 * .Machine$double.eps * ||z||) is read as one value:
# refused at creation, at a setData that derives and in an xbart fold, and
# kept, warning once per sampler, by the setters that re-derive. A column
# above it runs as built.

set.seed(20261107L, sample.kind = "Rejection")
n <- 200L
eps <- .Machine$double.eps
x1 <- runif(n)
y <- x1 + rnorm(n, 0, 0.3)
w <- rep(c(1, 0), each = 100L)
# 0.1 * 3 is 0.30000000000000004, one unit in the last place from 0.3
rounding <- rep(c(0.3, 0.1 * 3), length.out = n)
real <- runif(n, 0.2, 0.4)
# rounding on the rows of positive weight, real values on the rest
roundingIn <- ifelse(w == 1, rounding, real)
matrixOf <- function(z) cbind(x1 = x1, z = z)

controlWith <- function(...) {
  settings <- list(
    n.chains = 1L,
    n.threads = 1L,
    n.samples = 5L,
    n.burn = 0L,
    updateState = FALSE,
    seed = 7L
  )
  do.call(dbarts::dbartsControl, utils::modifyList(settings, list(...)))
}
make <- function(z, weights = NULL, leaf.prior = dbartsPriors$linear("z")) {
  dbarts::dbarts(
    matrixOf(z),
    y,
    weights = weights,
    leaf.prior = leaf.prior,
    control = controlWith(),
    n.trees = 10L
  )
}
recordOf <- function(s) lapply(s$data@standardization[[1L]], unname)
fallbacks <- function(expr) {
  caught <- character()
  value <- withCallingHandlers(
    expr,
    dbartsScaleFallbackWarning = function(w) {
      caught <<- c(caught, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = caught)
}
holds <- paste0(
  "the leaf covariate \"z\" holds values whose spread about their center is ",
  "under 64 times machine precision of their size (||z - center|| < 64 * ",
  ".Machine$double.eps * ||z||)"
)
cannot <- ", so its center and scale cannot be taken from them; "
# subtracting a constant from the stored values restores no lost digit, so
# the refusal advises recomputing from the source
recompute <- paste0(
  "or, if its values truly differ, recompute it from its source about an ",
  "origin near its values, "
)
refusal <- function(count, gp = FALSE) {
  paste0(
    holds,
    " over the rows in the likelihood (",
    count,
    " of 200 observations are in it)",
    cannot,
    if (gp) {
      "state its lengthscale in gp(), drop it from gp(), "
    } else {
      "drop it from linear(), "
    },
    recompute,
    "or change the rows in (the weights, the active-row mask, or subset)"
  )
}
keepText <- paste0(
  "fewer than two values of the leaf covariate \"z\" that differ by more ",
  "than rounding are in the likelihood, so its center and scale are kept ",
  "as they were (said once per sampler)"
)

# ---- creation ----------------------------------------------------------------

expect_true(length(unique(rounding)) == 2L)
expect_error(make(rounding), refusal(200L), fixed = TRUE)
expect_error(
  make(rounding, leaf.prior = dbartsPriors$gp("z")),
  refusal(200L, gp = TRUE),
  fixed = TRUE
)
# bart() reaches the same refusal
expect_error(
  bart(
    y ~ x1 + z,
    data.frame(y = y, x1 = x1, z = rounding),
    leaf.prior = linear("z"),
    n.samples = 2L,
    n.burn = 2L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  ),
  refusal(200L),
  fixed = TRUE
)
# a constant column keeps its own words
expect_error(make(rep(0.3, n)), "holds one value over the rows", fixed = TRUE)

# seconds since 1970 over 100 seconds (relative spread about 1.7e-8) and a
# column of sd 1e-6 at 0.3 have spread, and are scaled by their sd
times <- 1.7e9 + runif(n, 0, 100)
small <- 0.3 + rnorm(n, 0, 1e-6)
for (z in list(times, small)) {
  s <- make(z)
  expect_equal(recordOf(s), list(center = mean(z), scale = sd(z)))
}

# either side of the tolerance: half the rows 1 and half 1 + d, a spread of
# d / 2 about a center near 1
expect_error(
  make(rep(c(1, 1 + 126 * eps), length.out = n)),
  refusal(200L),
  fixed = TRUE
)
s <- make(rep(c(1, 1 + 130 * eps), length.out = n))
expect_true(recordOf(s)$scale > 0 && recordOf(s)$scale < 1e-13)

# the rows in decide it: rounding on the rows of positive weight is refused
# whatever the rows of weight 0 hold, and those rows' spread reads at weight 1
expect_error(make(roundingIn, weights = w), refusal(100L), fixed = TRUE)
s <- make(roundingIn)
expect_equal(
  recordOf(s),
  list(center = mean(roundingIn), scale = sd(roundingIn))
)

# a gp leaf that states its lengthscale takes the placeholder, divided by 1,
# and records it as it records a column of one value
s <- make(rounding, leaf.prior = dbartsPriors$gp("z", lengthscale = 0.5))
constant <- make(
  rep(0.3, n),
  leaf.prior = dbartsPriors$gp("z", lengthscale = 0.5)
)
expect_true(is.na(recordOf(s)$scale))
expect_identical(recordOf(s)$scale, recordOf(constant)$scale)
# centered at the corrected mean, within a unit in the last place of 0.3
expect_true(abs(recordOf(s)$center - 0.3) <= 1e-16)

# ---- setData -----------------------------------------------------------------

# deriving, at NULL and TRUE, refused and the sampler as it was; FALSE keeps
setDataRefusal <- paste0(
  "$setData: ",
  holds,
  " over the rows in the likelihood (200 of 200 observations are in it)",
  cannot,
  "drop it from linear(), recompute it from its source about an origin near ",
  "its values if its values truly differ, change the rows in (the weights, ",
  "the active-row mask, or subset), or keep the sampler's with ",
  "updateStandardization = FALSE"
)
s <- make(real)
before <- recordOf(s)
expect_error(
  s$setData(dbartsData(matrixOf(rounding), y)),
  setDataRefusal,
  fixed = TRUE
)
expect_error(
  s$setData(dbartsData(matrixOf(rounding), y), updateStandardization = TRUE),
  setDataRefusal,
  fixed = TRUE
)
expect_identical(recordOf(s), before)
s$setData(dbartsData(matrixOf(rounding), y), updateStandardization = FALSE)
expect_identical(recordOf(s), before)
# a column with spread derives
s$setData(dbartsData(matrixOf(times), y))
expect_equal(recordOf(s), list(center = mean(times), scale = sd(times)))

# ---- the setters that re-derive keep -----------------------------------------

# setPredictor, one column and the whole matrix
s <- make(real)
before <- recordOf(s)
kept <- fallbacks(
  s$setPredictor(rounding, "z", TRUE, updateStandardization = TRUE)
)
expect_identical(kept$warnings, keepText)
expect_identical(recordOf(s), before)
s <- make(real)
kept <- fallbacks(s$setPredictor(
  matrixOf(rounding),
  forceUpdate = TRUE,
  updateStandardization = TRUE
))
expect_identical(kept$warnings, keepText)
expect_identical(recordOf(s), before)
# a column with spread re-derives, in silence
s <- make(real)
kept <- fallbacks(
  s$setPredictor(small, "z", TRUE, updateStandardization = TRUE)
)
expect_identical(kept$warnings, character())
expect_equal(recordOf(s), list(center = mean(small), scale = sd(small)))

# setWeights and setActiveRows: the column has spread over every row, and
# differs by rounding only on the rows the call leaves in
s <- make(roundingIn)
before <- recordOf(s)
kept <- fallbacks(s$setWeights(w, updateStandardization = TRUE))
expect_identical(kept$warnings, keepText)
expect_identical(recordOf(s), before)
s <- make(roundingIn)
kept <- fallbacks(s$setActiveRows(w, updateStandardization = TRUE))
expect_identical(kept$warnings, keepText)
expect_identical(recordOf(s), before)
# rows in with spread re-derive over them
s <- make(roundingIn)
kept <- fallbacks(s$setWeights(1 - w, updateStandardization = TRUE))
expect_identical(kept$warnings, character())
expect_equal(
  recordOf(s),
  list(center = mean(real[w == 0]), scale = sd(real[w == 0]))
)
# a gp leaf keeps as the linear leaf does
s <- make(real, leaf.prior = dbartsPriors$gp("z"))
before <- recordOf(s)
kept <- fallbacks(
  s$setPredictor(rounding, "z", TRUE, updateStandardization = TRUE)
)
expect_identical(kept$warnings, keepText)
expect_identical(recordOf(s), before)

# ---- xbart -------------------------------------------------------------------

expect_error(
  xbart(
    matrixOf(rounding),
    y,
    leaf.prior = dbartsPriors$linear("z"),
    n.reps = 1L,
    n.samples = 2L,
    n.burn = c(2L, 2L),
    n.trees = 5L,
    n.threads = 1L,
    seed = 1L
  ),
  paste0(
    "in a fold of xbart, ",
    holds,
    " over the fold's training rows in the likelihood",
    cannot,
    "drop it from linear(), ",
    recompute,
    "or change the rows (the weights, subset or n.test)"
  ),
  fixed = TRUE
)

# ---- counts ------------------------------------------------------------------

# where the rows are counts (frequency weights, a binomial row's trials, an
# nbinom fit's whole-number weights) the test reads each row at its count, as
# the center and scale do, so it is the test of the rows written out: a
# spread of about 2.1 times the tolerance over the rows read once each is
# about 0.23 times it at these counts
spreadOnce <- rep(c(1, 1 + 1e-13), c(180L, 20L))
heavy <- ifelse(spreadOnce == 1, 100, 1)
countsFit <- function(
  response = y,
  z = spreadOnce,
  leaf.prior = dbartsPriors$linear("z"),
  ...
) {
  dbarts::dbarts(
    matrixOf(z),
    response,
    leaf.prior = leaf.prior,
    control = controlWith(),
    n.trees = 10L,
    ...
  )
}
# precision weights do not count rows, so the column has spread
s <- countsFit(weights = heavy)
expect_equal(
  recordOf(s),
  list(center = mean(spreadOnce), scale = sd(spreadOnce))
)
expect_error(
  countsFit(weights = heavy, weights.type = "frequency"),
  refusal(200L),
  fixed = TRUE
)
expect_error(
  countsFit(rpois(n, 3), weights = heavy, family = "nbinom"),
  refusal(200L),
  fixed = TRUE
)
successes <- rbinom(n, as.integer(heavy), 0.4)
expect_error(
  dbarts::dbarts(
    cbind(k, f) ~ x1 + z,
    data.frame(k = successes, f = heavy - successes, x1 = x1, z = spreadOnce),
    family = binomial,
    leaf.prior = dbartsPriors$linear("z"),
    control = controlWith(),
    n.trees = 5L
  ),
  refusal(200L),
  fixed = TRUE
)
# setWeights to those counts under frequency keeps, under precision derives
s <- countsFit(weights = rep(1, n), weights.type = "frequency")
before <- recordOf(s)
kept <- fallbacks(s$setWeights(heavy, updateStandardization = TRUE))
expect_identical(kept$warnings, keepText)
expect_identical(recordOf(s), before)
s <- countsFit(weights = rep(1, n))
kept <- fallbacks(s$setWeights(heavy, updateStandardization = TRUE))
expect_identical(kept$warnings, character())
# the center is corrected at the counts: half the rows at 1, count 1000, and
# half at 1 + 1e-13, count 1, are 0.22 times the tolerance written out
expect_error(
  countsFit(
    z = rep(c(1, 1 + 1e-13), each = 100L),
    weights = rep(c(1000, 1), each = 100L),
    weights.type = "frequency"
  ),
  refusal(200L),
  fixed = TRUE
)

# the doors that derive again read the counts: a binomial sampler's trials at
# setPredictor and setData, frequency weights at setActiveRows and in an xbart
# fold
realCounts <- runif(n, 0.5, 1.5)
binomialFrame <- function(z) {
  data.frame(k = successes, f = heavy - successes, x1 = x1, z = z)
}
binomialFit <- function() {
  dbarts::dbarts(
    cbind(k, f) ~ x1 + z,
    binomialFrame(realCounts),
    family = binomial,
    leaf.prior = dbartsPriors$linear("z"),
    control = controlWith(),
    n.trees = 5L
  )
}
s <- binomialFit()
before <- recordOf(s)
kept <- fallbacks(
  s$setPredictor(spreadOnce, "z", TRUE, updateStandardization = TRUE)
)
expect_identical(kept$warnings, keepText)
expect_identical(recordOf(s), before)
s <- binomialFit()
expect_error(
  s$setData(dbartsData(
    matrixOf(spreadOnce),
    counts = cbind(successes, heavy - successes)
  )),
  setDataRefusal,
  fixed = TRUE
)
expect_identical(recordOf(s), before)
# 20 more rows of real values give the column spread until setActiveRows
# leaves them out
extended <- function(...) {
  dbarts::dbarts(
    cbind(x1 = c(x1, runif(20L)), z = c(spreadOnce, runif(20L, 0.5, 1.5))),
    c(y, rnorm(20L)),
    weights = c(heavy, rep(1, 20L)),
    leaf.prior = dbartsPriors$linear("z"),
    control = controlWith(),
    n.trees = 5L,
    ...
  )
}
active <- rep(c(1, 0), c(n, 20L))
s <- extended(weights.type = "frequency")
before <- recordOf(s)
kept <- fallbacks(s$setActiveRows(active, updateStandardization = TRUE))
expect_identical(kept$warnings, keepText)
expect_identical(recordOf(s), before)
s <- extended()
kept <- fallbacks(s$setActiveRows(active, updateStandardization = TRUE))
expect_identical(kept$warnings, character())
expect_equal(
  recordOf(s),
  list(center = mean(spreadOnce), scale = sd(spreadOnce))
)
expect_error(
  xbart(
    matrixOf(spreadOnce),
    y,
    weights = heavy,
    weights.type = "frequency",
    leaf.prior = dbartsPriors$linear("z"),
    n.reps = 1L,
    n.samples = 2L,
    n.burn = c(2L, 2L),
    n.trees = 5L,
    n.threads = 1L,
    seed = 1L
  ),
  "in a fold of xbart, the leaf covariate \"z\" holds values whose spread",
  fixed = TRUE
)

# a gp leaf stating its lengthscale records the column as the rows written
# out record it: one value, its scale NA
written <- rep(seq_len(n), heavy)
gpStated <- dbartsPriors$gp("z", lengthscale = 0.5)
s <- countsFit(
  weights = heavy,
  weights.type = "frequency",
  leaf.prior = gpStated
)
writtenOut <- dbarts::dbarts(
  cbind(x1 = x1[written], z = spreadOnce[written]),
  y[written],
  leaf.prior = gpStated,
  control = controlWith(),
  n.trees = 5L
)
expect_true(is.na(recordOf(s)$scale))
expect_identical(recordOf(s)$scale, recordOf(writtenOut)$scale)
# where fewer than two rows have a positive count the mark reads every row, a
# count of 0 as 1 and a positive row at its count: one row of count 5 over a
# column within rounding of its center on every row records scale NA, as its
# rows written out, five copies of one row, do
caught <- list()
s <- withCallingHandlers(
  countsFit(
    z = rounding,
    weights = c(5, rep(0, n - 1L)),
    weights.type = "frequency",
    leaf.prior = gpStated
  ),
  warning = function(w) {
    caught[[length(caught) + 1L]] <<- w
    invokeRestart("muffleWarning")
  }
)
expect_true(is.na(recordOf(s)$scale))
expect_identical(length(caught), 2L)
expect_identical(
  conditionMessage(caught[[1L]]),
  paste0(
    "rows with 'weights' of 0 are left out of the likelihood but stay in ",
    "the fit; use 'subset' to drop them"
  )
)
expect_true(inherits(caught[[2L]], "dbartsScaleFallbackWarning"))

# a setData reads the old mark at the old counts, before they move to the new
# rows: a linear leaf at the placeholder scale over a column within rounding
# of its center at counts 1, 1 and 1000, and not read once each, set onto a
# record at another center whose rows hold one value keeps each block for the
# column, neither side having spread; read at the new counts the old side has
# spread, and the intercept moves by the slope times 1000
onRecord <- function(x, weights, center) {
  data <- dbarts::dbartsData(
    x,
    x[, 1L] + rnorm(nrow(x), 0, 0.1),
    weights = weights
  )
  data@standardization <- list(list(
    center = c(z = center),
    scale = c(z = NA_real_)
  ))
  data
}
s <- dbarts::dbarts(
  onRecord(cbind(x1 = runif(40L), z = 1000), rep(c(1, 1000), c(2L, 38L)), 1000),
  leaf.prior = dbartsPriors$linear("z"),
  weights.type = "frequency",
  control = controlWith(),
  n.trees = 5L
)
invisible(s$run())
expect_true(
  s$setPredictor(c(1000 + 1e-10, 1000 - 1e-10, rep(1000, 38L)), "z")
)
newRows <- cbind(x1 = runif(80L), z = 2000)
s$setData(onRecord(newRows, rep(1, 80L), 2000))
expect_true(max(abs(s$predict(newRows))) < 10)
