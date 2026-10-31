# $setData takes the cut count and the quantile rule its data object states, one
# carrying a grid bringing what made it, and the sampler's record follows, so a
# later refresh derives by them. An object stating neither keeps the sampler's.

set.seed(61)
n <- 150L
makeFrame <- function(rate) {
  x1 <- rexp(n, rate)
  x2 <- runif(n)
  data.frame(y = x1 + x2 + rnorm(n, sd = 0.3), x1 = x1, x2 = x2)
}
frame <- makeFrame(1)
other <- makeFrame(0.2)
run <- function() {
  dbartsControl(n.samples = 2L, n.burn = 0L, n.chains = 1L, n.threads = 1L)
}
fit <- function(data = frame, ...) {
  sampler <- dbarts(y ~ x1 + x2, data, ..., control = run())
  invisible(sampler$run(0L, 2L))
  sampler
}
# an even grid's points are equally spaced; quantiles on the skewed x1 are not
isEven <- function(grid) {
  isTRUE(all.equal(diff(grid), rep(mean(diff(grid)), length(grid) - 1L)))
}
grids <- function(sampler) sampler$data@cutPoints
expectRecord <- function(sampler, cuts, rule, info = "") {
  expect_identical(sampler$data@n.cuts, cuts, info = info)
  expect_identical(sampler$data@useQuantiles, rule, info = info)
  expect_identical(sampler$control@useQuantiles, rule, info = info)
  expect_identical(lengths(grids(sampler)), cuts, info = info)
  expect_identical(isEven(grids(sampler)[[1L]]), !rule, info = info)
}
refresh <- function(sampler) {
  sampler$setPredictor(
    rexp(n, 0.5),
    1L,
    updateCutPoints = "position",
    forceUpdate = TRUE
  )
  sampler
}

# both statements, with updateCutPoints left NULL and given to derive: the grid
# is derived at what the object states and a refresh goes on deriving by it
statesCount <- dbartsData(y ~ x1 + x2, other, n.cuts = 10L)
statesRule <- dbartsData(y ~ x1 + x2, other, useQuantiles = TRUE)
for (given in list(list(), list(updateCutPoints = "value"), list(updateCutPoints = "position"))) {
  info <- paste(unlist(given), collapse = "")
  sampler <- fit(n.cuts = 25L)
  do.call(sampler$setData, c(list(statesCount), given))
  expectRecord(sampler, c(10L, 10L), FALSE, info)
  expectRecord(refresh(sampler), c(10L, 10L), FALSE, info)

  sampler <- fit(n.cuts = 25L)
  do.call(sampler$setData, c(list(statesRule), given))
  expectRecord(sampler, c(25L, 25L), TRUE, info)
  expectRecord(refresh(sampler), c(25L, 25L), TRUE, info)
}

# a call that keeps the grid in force keeps what made it, whatever is stated
for (object in list(statesCount, statesRule)) {
  sampler <- fit(n.cuts = 25L)
  kept <- grids(sampler)
  sampler$setData(object, updateCutPoints = "none")
  expect_identical(grids(sampler), kept)
  expectRecord(sampler, c(25L, 25L), FALSE)
}

# another fit's data object brings its grid and what made it, so a refresh of
# one column gives that column 10 quantile points, not 25 even beside 10 quantile
tenQuantile <- fit(other, n.cuts = 10L, useQuantiles = TRUE)
sampler <- fit(n.cuts = 25L)
sampler$setData(tenQuantile$data)
expectRecord(sampler, c(10L, 10L), TRUE)
expect_identical(grids(sampler), grids(tenQuantile))
expectRecord(refresh(sampler), c(10L, 10L), TRUE)

# an object carrying a grid and no statement of a rule reads as even, on a
# sampler at quantiles
tenEven <- fit(other, n.cuts = 10L)
sampler <- fit(n.cuts = 25L, useQuantiles = TRUE)
sampler$setData(tenEven$data)
expectRecord(sampler, c(10L, 10L), FALSE)
expect_false(dbarts:::dataStatesRule(sampler$data))
expectRecord(refresh(sampler), c(10L, 10L), FALSE)
handGrid <- dbartsData(y ~ x1 + x2, other)
handGrid@cutPoints <- list(c(1, 2, 4), c(0.25, 0.75))
sampler <- fit(n.cuts = 25L, useQuantiles = TRUE)
sampler$setData(handGrid)
expect_identical(lengths(grids(sampler)), c(3L, 2L))
expect_identical(sampler$data@n.cuts, c(25L, 25L))
expect_false(sampler$data@useQuantiles)
expect_true(isEven(grids(refresh(sampler))[[1L]]))

# one given to derive brings the same, its grid not taken
sampler <- fit(n.cuts = 25L)
sampler$setData(tenQuantile$data, updateCutPoints = "value")
expectRecord(sampler, c(10L, 10L), TRUE)

# an object stating nothing keeps the sampler's, and whether it was stated
sampler <- fit(n.cuts = 25L, useQuantiles = TRUE)
sampler$setData(dbartsData(y ~ x1 + x2, other))
expectRecord(sampler, c(25L, 25L), TRUE)
expect_true(dbarts:::dataStatesRule(sampler$data))
expectRecord(refresh(sampler), c(25L, 25L), TRUE)
sampler <- fit(n.cuts = 25L)
sampler$setData(dbartsData(y ~ x1 + x2, other))
expectRecord(sampler, c(25L, 25L), FALSE)
expect_false(dbarts:::dataStatesRule(sampler$data))

# a stated FALSE is a statement, on a sampler at quantiles
sampler <- fit(n.cuts = 25L, useQuantiles = TRUE)
sampler$setData(dbartsData(y ~ x1 + x2, other, useQuantiles = FALSE))
expectRecord(sampler, c(25L, 25L), FALSE)
expect_true(dbarts:::dataStatesRule(sampler$data))

# a count moves with its columns
sampler <- fit(n.cuts = c(30L, 30L))
sampler$setData(dbartsData(y ~ x2 + x1, other, n.cuts = c(7L, 12L)))
expect_identical(sampler$data@n.cuts, c(12L, 7L))
expect_identical(lengths(grids(sampler)), c(12L, 7L))
expect_identical(colnames(sampler$data@x), c("x1", "x2"))

# what it takes is the record a copy and a reload are made on, with no message
sampler <- fit(n.cuts = 25L)
sampler$setData(tenQuantile$data)
expect_silent(copied <- sampler$copy())
expectRecord(copied, c(10L, 10L), TRUE)
saved <- tempfile(fileext = ".rds")
sampler$storeState()
saveRDS(sampler, saved)
expect_silent({
  reloaded <- readRDS(saved)
  reloaded$getPointer()
})
unlink(saved)
expectRecord(reloaded, c(10L, 10L), TRUE)
expectRecord(refresh(reloaded), c(10L, 10L), TRUE)

# $setControl keeps refusing both, against what the sampler now holds
written <- function(...) suppressWarnings(dbartsControl(n.chains = 1L, n.threads = 1L, ...))
expect_error(
  sampler$setControl(written(n.cuts = 25L)),
  "changing 'n.cuts' is not available on an existing sampler"
)
expect_error(
  sampler$setControl(written(useQuantiles = FALSE)),
  "changing 'useQuantiles' is not available on an existing sampler"
)
