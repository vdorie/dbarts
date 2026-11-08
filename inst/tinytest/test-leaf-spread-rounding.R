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
