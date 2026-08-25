# sparse predictor matrices: a Matrix::dgCMatrix enters through the x/y
# interface, columns below a density threshold take rank-bitmap storage in
# the engine, denser ones densify their codes, and the raw-x mutation
# surface is fixed at creation

if (!requireNamespace("Matrix", quietly = TRUE)) {
  exit_file("Matrix not available")
}

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
fit.sparse <- bart(
  x.sparse,
  y,
  ndpost = 300L,
  nskip = 100L,
  ntree = 50L,
  verbose = FALSE
)
fit.dense <- bart(
  x.dense,
  y,
  ndpost = 300L,
  nskip = 100L,
  ntree = 50L,
  verbose = FALSE
)
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
fit.na <- bart2(
  x.na,
  y,
  n.samples = 20L,
  n.burn = 20L,
  n.trees = 25L,
  verbose = FALSE
)
expect_equal(length(fitted(fit.na)), n)
expect_error(dbartsData(x.na, y, missing = "error"), pattern = "missing values")
expect_error(
  bart(x.na, y, ndpost = 20L, nskip = 20L, ntree = 25L, verbose = FALSE),
  pattern = "missing values"
)

# test predictions work off a dense x.test
x.test <- x.dense[1:20, , drop = FALSE]
fit.test <- bart(
  x.sparse,
  y,
  x.test,
  ndpost = 50L,
  nskip = 50L,
  ntree = 25L,
  verbose = FALSE
)
expect_equal(ncol(fit.test$yhat.test), 20L)

# the sampler surface: response-side mutation stays open; a sparse column
# accepts column-granular and whole-matrix between-sweep mutation;
# per-observation and whole-data replacement stay fixed at creation, and
# grouped rbart_vi is reserved
control <- dbartsControl(
  n.samples = 10L,
  n.burn = 0L,
  n.trees = 25L,
  updateState = FALSE
)
sampler <- dbarts(x.sparse, y, control = control)
invisible(sampler$run())
expect_silent(sampler$setResponse(y))

# column-granular mutation of a sparse column: pattern-preserving (rescale the
# stored nonzeros) then pattern-changing (plant a nonzero at a former zero),
# forced so it always installs; data@x is maintained and stays a dgCMatrix
newCol1 <- x.dense[, 1L]
newCol1[newCol1 != 0] <- newCol1[newCol1 != 0] * 1.5
expect_silent(
  sampler$setPredictor(newCol1, column = 1L, forceUpdate = TRUE)
)
expect_inherits(sampler$data@x, "dgCMatrix")
expect_equal(as.numeric(sampler$data@x[, 1L]), unname(newCol1))

newCol2 <- newCol1
newCol2[which(newCol2 == 0)[1:3]] <- c(0.6, 0.7, 0.8)
expect_silent(
  sampler$setPredictor(newCol2, column = 1L, forceUpdate = TRUE)
)
expect_equal(as.numeric(sampler$data@x[, 1L]), unname(newCol2))
# a transactional (non-forced) column update returns a logical accept flag and
# leaves the sampler runnable
accepted <- sampler$setPredictor(newCol1, column = 1L)
expect_true(is.logical(accepted) && length(accepted) == 1L)
expect_silent(invisible(sampler$run()))

# whole-matrix replacement of a sparse design: the container is spliced R-side
# and stays a dgCMatrix, with the replaced columns densifying their storage
# (every row now differs from the implicit zero)
x.whole <- x.dense
x.whole[, 2L] <- x.whole[, 2L] + 0.5
expect_silent(sampler$setPredictor(x.whole, forceUpdate = TRUE))
expect_inherits(sampler$data@x, "dgCMatrix")
expect_equal(unname(as.matrix(sampler$data@x)), unname(x.whole))
expect_equal(diff(sampler$data@x@p)[2L], n)
# a rejected transactional whole-matrix replacement leaves data@x untouched
x.before <- as.matrix(sampler$data@x)
accepted.whole <- sampler$setPredictor(matrix(0, n, p), forceUpdate = FALSE)
expect_false(isTRUE(accepted.whole))
expect_equal(as.matrix(sampler$data@x), x.before)

# out of scope, refused by name: per-observation replacement of a sparse
# column, and whole-data replacement
expect_error(
  sampler$setPredictor(newCol1, column = 1L, forceUpdate = "partial"),
  pattern = "per-observation updates require a dense-backed column"
)
expect_error(sampler$setData(dbartsData(x.sparse, y)), pattern = "sparse")
expect_error(
  dbarts(x.sparse, y, node.prior = linear(c("x1", "x2"))),
  pattern = "sparse"
)
# the refusal names the wrappers that do accept sparse, rather than dead-ending
expect_error(
  rbart_vi(
    x.sparse,
    y,
    group.by = rep_len(1:4, n),
    n.samples = 5L,
    n.burn = 5L,
    n.trees = 25L,
    n.chains = 1L,
    n.threads = 1L
  ),
  pattern = "dbarts.* and bart2.* do"
)

# save/load: sampler re-creation from the stored dgCMatrix restores state
control.state <- dbartsControl(
  n.samples = 10L,
  n.burn = 0L,
  n.trees = 25L,
  n.chains = 2L,
  n.threads = 1L,
  updateState = TRUE
)
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
xval <- xbart(
  x.sparse,
  y,
  n.samples = 40L,
  n.burn = c(20L, 5L, 5L),
  n.trees = 25L,
  n.reps = 2L,
  n.test = 4L,
  n.threads = 1L
)
expect_true(all(is.finite(xval)))

# a bare dgCMatrix 'test' set: accepted symmetrically with a bare dgCMatrix
# 'x' - it stays resident as an all-sparse mixed container instead of
# silently densifying, taking exactly the path a mixed container's own
# sparse test columns already do
x.test.dense <- x.dense[1:20, , drop = FALSE]
x.test.sparse <- x.sparse[1:20, , drop = FALSE]

# ingestion: train sparse + test sparse
data.sparseTest <- dbartsData(x.sparse, y, test = x.test.sparse)
expect_inherits(data.sparseTest@x.test, "dbartsMixedMatrix")
expect_true(dbarts:::predictorSourceIsSparse(data.sparseTest@x.test))

# ingestion: train dense (all-ordinal) + test sparse also stays resident -
# quantization only reads the training cut grid, which exists regardless of
# how the training store itself was built
data.denseTrainSparseTest <- dbartsData(x.dense, y, test = x.test.sparse)
expect_inherits(data.denseTrainSparseTest@x.test, "dbartsMixedMatrix")
expect_true(dbarts:::predictorSourceIsSparse(data.denseTrainSparseTest@x.test))

fitControl <- dbartsControl(
  n.samples = 15L,
  n.burn = 10L,
  n.trees = 15L,
  n.chains = 1L,
  n.threads = 1L,
  updateState = FALSE
)

# train sparse + test sparse fits bitwise-identically to train sparse + the
# same test values supplied dense
set.seed(201)
run.sparseTrain.sparseTest <- dbarts(
  x.sparse,
  y,
  test = x.test.sparse,
  control = fitControl
)$run()
set.seed(201)
run.sparseTrain.denseTest <- dbarts(
  x.sparse,
  y,
  test = x.test.dense,
  control = fitControl
)$run()
expect_identical(
  run.sparseTrain.sparseTest$test,
  run.sparseTrain.denseTest$test
)

# train dense + test sparse fits bitwise-identically to train dense + the
# same test values supplied dense
set.seed(202)
run.denseTrain.sparseTest <- dbarts(
  x.dense,
  y,
  test = x.test.sparse,
  control = fitControl
)$run()
set.seed(202)
run.denseTrain.denseTest <- dbarts(
  x.dense,
  y,
  test = x.test.dense,
  control = fitControl
)$run()
expect_identical(run.denseTrain.sparseTest$test, run.denseTrain.denseTest$test)

# column alignment goes through the same name-matching code as a dense test
# matrix: a reordered, named sparse test set matches by name, not position
permutation <- c(3L, 1L, 4L, 2L, 6L, 5L)
set.seed(203)
run.permuted.sparseTest <- dbarts(
  x.dense,
  y,
  test = x.test.sparse[, permutation, drop = FALSE],
  control = fitControl
)$run()
set.seed(203)
run.permuted.denseTest <- dbarts(
  x.dense,
  y,
  test = x.test.dense[, permutation, drop = FALSE],
  control = fitControl
)$run()
expect_identical(run.permuted.sparseTest$test, run.permuted.denseTest$test)

# predict() on a bart2 fit takes a bare dgCMatrix newdata the same path
fit.forPredict <- bart2(
  x.dense,
  y,
  n.samples = 10L,
  n.burn = 5L,
  n.trees = 10L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrees = TRUE
)
expect_identical(
  predict(fit.forPredict, x.test.sparse),
  predict(fit.forPredict, x.test.dense)
)

# train dense with a categorical (factor) column + a bare dgCMatrix test:
# a numeric sparse matrix carries no factor levels, so it cannot supply that
# column's values - refuse informatively rather than error out of the bridge
df.withFactor <- data.frame(
  x1 = x.dense[, 1L],
  x2 = x.dense[, 2L],
  f = factor(sample(c("a", "b"), n, replace = TRUE))
)
expect_error(
  dbartsData(
    df.withFactor,
    y,
    test = x.test.sparse[, 1:2, drop = FALSE]
  ),
  pattern = "categorical training column"
)

# the class validity union accepts a bare dgCMatrix 'x.test' directly - it
# is converted to the resident container form before storage (so a
# validated object never actually holds one), but the union itself stays
# symmetric with 'x' rather than rejecting it outright
data.rawSparseTestSlot <- data.sparseTest
data.rawSparseTestSlot@x.test <- x.test.sparse
expect_silent(methods::validObject(data.rawSparseTestSlot))
