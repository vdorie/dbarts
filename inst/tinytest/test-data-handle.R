# The shared data handle and its row-subset view samplers
# (public-surface.md section 5; internal). Views copy the handle's cut grid
# and gather their rows' codes, so folds bin identically to the full data;
# they refuse raw-predictor mutation.

set.seed(42)
n <- 150L
p <- 4L
x <- matrix(runif(n * p), n, p)
y <- rowSums(x[, 1:2]) + rnorm(n, 0, 0.5)

control <- dbartsControl(n.chains = 1L, n.threads = 1L, n.trees = 20L,
                         updateState = FALSE)
sampler <- dbarts(x, y, control = control)
handle <- dbarts:::bartcoreDataHandle(sampler$control, sampler$data)

# a view over every row reproduces the ordinary sampler bitwise: same
# store contents, same single-chain draws through R's generator
bc.full <- dbarts:::bartcoreSampler(sampler)
bc.view <- dbarts:::bartcoreSamplerFromHandle(handle, sampler$control,
                                              sampler$model, sampler$data,
                                              seq_len(n))
set.seed(7)
r.full <- dbarts:::bartcoreRun(bc.full, 20L, 30L)
set.seed(7)
r.view <- dbarts:::bartcoreRun(bc.view, 20L, 30L)
expect_identical(r.view$sigma, r.full$sigma)
expect_identical(r.view$train, r.full$train)

# a fold: test rows index the handle's observations
testRows  <- seq(1L, n, by = 5L)
trainRows <- setdiff(seq_len(n), testRows)
bc.fold <- dbarts:::bartcoreSamplerFromHandle(handle, sampler$control,
                                              sampler$model, sampler$data,
                                              trainRows, testRows)
r.fold <- dbarts:::bartcoreRun(bc.fold, 20L, 30L)
expect_equal(dim(r.fold$train), c(length(trainRows), 30L))
expect_equal(dim(r.fold$test), c(length(testRows), 30L))
expect_true(all(is.finite(r.fold$test)))

# raw-predictor mutation is refused on views...
x.fold <- x[trainRows, , drop = FALSE]
expect_error(dbarts:::bartcoreSetPredictor(bc.fold, x.fold),
             pattern = "owns its predictors")
expect_error(dbarts:::bartcoreUpdatePredictor(bc.fold, x.fold[, 1L], 1L),
             pattern = "owns its predictors")
expect_error(
  dbarts:::bartcoreUpdatePredictorPerObservation(bc.fold, x.fold[, 1L], 1L),
  pattern = "owns its predictors")
expect_error(dbarts:::bartcoreSetCutPoints(bc.fold, list(c(0.25, 0.5)), 1L),
             pattern = "owns its predictors")
expect_error(dbarts:::bartcoreSetData(bc.fold, sampler$data),
             pattern = "owns its predictors")
expect_error(
  dbarts:::bartcoreSetState(bc.fold, dbarts:::bartcoreStoreState(bc.fold)),
  pattern = "owns its predictors")

# ...while the response side and raw test data stay available (test
# quantization needs only the copied cut grid)
expect_silent(dbarts:::bartcoreSetResponse(bc.fold, y[trainRows] + 0.1))
expect_silent(
  dbarts:::bartcoreSetTestPredictor(bc.fold, x[testRows, , drop = FALSE]))
predictions <- dbarts:::bartcorePredict(bc.fold, x[1:3, , drop = FALSE])
expect_true(all(is.finite(predictions)))

# view samplers bin on the full data's grid: hold out the rows carrying
# column 1's extremes and check every split against the parent's cut set
x.wide <- x
x.wide[1L, 1L] <- 0
x.wide[2L, 1L] <- 1
control.trees <- dbartsControl(n.chains = 1L, n.threads = 1L, n.trees = 10L,
                               keepTrees = TRUE, n.samples = 10L,
                               updateState = FALSE)
sampler.wide <- dbarts(x.wide, y, control = control.trees)
handle.wide <- dbarts:::bartcoreDataHandle(sampler.wide$control,
                                           sampler.wide$data)
bc.interior <- dbarts:::bartcoreSamplerFromHandle(
  handle.wide, sampler.wide$control, sampler.wide$model, sampler.wide$data,
  3:n, 1:2)
invisible(dbarts:::bartcoreRun(bc.interior, 20L, 10L))
trees <- dbarts:::bartcoreGetTrees(bc.interior, chainNums = 1L,
                                   treeNums = 1:10, current = TRUE)
splits <- trees[trees$var > 0L, ]
fullGrid <- lapply(seq_len(p), function(j) {
  r <- range(x.wide[, j])
  r[1L] + seq_len(100L) * diff(r) / 101
})
onGrid <- vapply(seq_len(nrow(splits)), function(i) {
  min(abs(fullGrid[[splits$var[i]]] - splits$value[i])) < 1e-12
}, logical(1L))
expect_true(length(onGrid) > 0L && all(onGrid))

# shape mismatches between the data object and the handle are refused
data.short <- dbartsData(x[-1L, , drop = FALSE], y[-1L])
data.short@n.cuts <- rep_len(100L, p)
data.short@sigma <- sampler$data@sigma
expect_error(
  dbarts:::bartcoreSamplerFromHandle(handle, sampler$control, sampler$model,
                                     data.short, seq_len(n - 1L)),
  pattern = "shape")
expect_error(
  dbarts:::bartcoreSamplerFromHandle(handle, sampler$control, sampler$model,
                                     sampler$data, c(0L, 1L)),
  pattern = "out of range")
expect_error(
  dbarts:::bartcoreSamplerFromHandle(handle, sampler$control, sampler$model,
                                     sampler$data, seq_len(n), n + 1L),
  pattern = "out of range")

# offsets: the view slices the regular offset by row and takes its test
# offset from offset[testRows]; a constant offset must appear in test fits
offset <- rep(100, n)
sampler.off <- dbarts(x, y + offset, offset = offset, control = control)
handle.off <- dbarts:::bartcoreDataHandle(sampler.off$control,
                                          sampler.off$data)
bc.off <- dbarts:::bartcoreSamplerFromHandle(handle.off, sampler.off$control,
                                             sampler.off$model,
                                             sampler.off$data,
                                             trainRows, testRows)
r.off <- dbarts:::bartcoreRun(bc.off, 20L, 30L)
expect_true(abs(mean(r.off$test) - mean(y[testRows] + 100)) < 1)
