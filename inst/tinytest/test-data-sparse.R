# sparse predictor matrices: a Matrix::dgCMatrix enters through the x/y
# interface, columns below a density threshold take rank-bitmap storage in
# the engine, denser ones densify their codes, and the raw-x mutation
# surface is fixed at creation (docs/design/sparse-columns.md)

if (!requireNamespace("Matrix", quietly = TRUE)) exit_file("Matrix not available")

set.seed(42)
n <- 400L
p <- 6L
x.dense <- matrix(0, n, p)
for (j in seq_len(p)) {
  nonzero <- runif(n) < 0.1
  x.dense[nonzero, j] <- 0.5 + runif(sum(nonzero))
}
colnames(x.dense) <- paste0("x", seq_len(p))
x.sparse <- methods::as(x.dense, "CsparseMatrix")
y <- 2 * x.dense[, 1L] - 1.5 * x.dense[, 2L] + rnorm(n, 0, 0.3)

# ingestion: the data object holds the dgCMatrix itself, all columns ordinal
data.sparse <- dbartsData(x.sparse, y)
expect_inherits(data.sparse@x, "dgCMatrix")
expect_equal(data.sparse@varTypes, rep.int(0L, p))
expect_equal(length(data.sparse@y), n)

# subsetting composes
data.subset <- dbartsData(x.sparse, y, subset = 1:200)
expect_equal(length(data.subset@y), 200L)
expect_equal(nrow(data.subset@x), 200L)

# a sparse fit runs and recovers the signal a dense fit of the same values
# finds; draws are not bitwise comparable (rank columns partition through a
# different kernel), so the comparison is on fit quality
fit.sparse <- bart(x.sparse, y, ndpost = 300L, nskip = 100L, ntree = 50L,
                   verbose = FALSE)
fit.dense <- bart(x.dense, y, ndpost = 300L, nskip = 100L, ntree = 50L,
                  verbose = FALSE)
f <- 2 * x.dense[, 1L] - 1.5 * x.dense[, 2L]
sse.sparse <- sum((fit.sparse$yhat.train.mean - f)^2)
sse.dense <- sum((fit.dense$yhat.train.mean - f)^2)
sse.mean <- sum((mean(y) - f)^2)
expect_true(sse.sparse < 0.25 * sse.mean)
expect_true(sse.sparse < 2.0 * sse.dense + 1.0)
# the sparse sigma estimate is the marginal fallback, not the linear fit's
expect_true(is.finite(fit.sparse$sigest) && fit.sparse$sigest > 0)

# variable names ride through to varcount
expect_equal(colnames(fit.sparse$varcount), colnames(x.dense))

# missing values are stored NaN entries and route through MIA (bart pins
# missing = "error", so the NA fit goes through bart2)
x.na <- x.sparse
x.na[3L, 1L] <- NA_real_
fit.na <- bart2(x.na, y, n.samples = 20L, n.burn = 20L, n.trees = 25L,
                verbose = FALSE)
expect_equal(length(fitted(fit.na)), n)
expect_error(dbartsData(x.na, y, missing = "error"),
             pattern = "missing values")
expect_error(bart(x.na, y, ndpost = 20L, nskip = 20L, ntree = 25L,
                  verbose = FALSE),
             pattern = "missing values")

# test predictions work off a dense x.test
x.test <- x.dense[1:20, , drop = FALSE]
fit.test <- bart(x.sparse, y, x.test, ndpost = 50L, nskip = 50L,
                 ntree = 25L, verbose = FALSE)
expect_equal(ncol(fit.test$yhat.test), 20L)

# the sampler surface: response-side mutation stays open, the raw-x
# surface and whole-data replacement are fixed at creation, and grouped
# rbart_vi is reserved
control <- dbartsControl(n.samples = 10L, n.burn = 0L, n.trees = 25L,
                         updateState = FALSE)
sampler <- dbarts(x.sparse, y, control = control)
invisible(sampler$run())
expect_silent(sampler$setResponse(y))
expect_error(sampler$setPredictor(x.dense),
             pattern = "sparse predictors fix the design")
expect_error(sampler$setData(dbartsData(x.sparse, y)),
             pattern = "sparse")
expect_error(dbarts(x.sparse, y, node.prior = linear(c("x1", "x2"))),
             pattern = "sparse")
# the refusal names the wrappers that do accept sparse, rather than dead-ending
expect_error(rbart_vi(x.sparse, y, group.by = rep_len(1:4, n),
                      n.samples = 5L, n.burn = 5L, n.trees = 25L,
                      n.chains = 1L, n.threads = 1L),
             pattern = "dbarts.* and bart2.* do")

# save/load: sampler re-creation from the stored dgCMatrix restores state
control.state <- dbartsControl(n.samples = 10L, n.burn = 0L, n.trees = 25L,
                               n.chains = 2L, n.threads = 1L,
                               updateState = TRUE)
sampler.state <- dbarts(x.sparse, y, control = control.state)
set.seed(99)
run.before <- sampler.state$run(numBurnIn = 20L, numSamples = 10L)
sampler.state$storeState()
tempFile <- tempfile(fileext = ".rds")
saveRDS(sampler.state, tempFile)
sampler.loaded <- readRDS(tempFile)
file.remove(tempFile)
set.seed(101)
run.original <- sampler.state$run(numBurnIn = 0L, numSamples = 10L)
set.seed(101)
run.restored <- sampler.loaded$run(numBurnIn = 0L, numSamples = 10L)
expect_equal(run.original$sigma, run.restored$sigma)
expect_equal(run.original$train, run.restored$train)

# xbart goes through the data-handle path; folds densify internally
xval <- xbart(x.sparse, y, n.samples = 40L, n.burn = c(20L, 5L, 5L),
              n.trees = 25L, n.reps = 2L, n.test = 4L, n.threads = 1L)
expect_true(all(is.finite(xval)))
