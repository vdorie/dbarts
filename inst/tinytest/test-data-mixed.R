# mixed dense and sparse predictors: a data frame holding numeric and
# factor columns alongside Matrix::sparseVector or dgCMatrix columns enters
# through the x/y interface, dense-backed columns keep categorical splits
# and linear leaves, and CSC-backed columns tier exactly as an all-sparse
# design does (docs/design/sparse-columns.md)

if (!requireNamespace("Matrix", quietly = TRUE)) exit_file("Matrix not available")

set.seed(42)
n <- 400L
x1 <- rnorm(n)
f <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
sv <- Matrix::sparseVector(x = 0.5 + runif(60L), i = sort(sample.int(n, 60L)),
                           length = n)
sm.dense <- matrix(0, n, 2L)
for (j in 1:2) {
  nonzero <- runif(n) < 0.1
  sm.dense[nonzero, j] <- 0.5 + runif(sum(nonzero))
}
sm <- methods::as(sm.dense, "CsparseMatrix")
x.frame <- data.frame(x1 = x1, f = f)
x.frame$sv <- sv
x.frame$sm <- sm
f.true <- 2 * as.double(sv) - 1.5 * x1 + 0.5 * (f == "b")
y <- f.true + rnorm(n, 0, 0.3)

# the categorical builder passes sparse columns through unexpanded
mm <- dbarts:::makeCategoricalModelMatrix(x.frame)
expect_inherits(mm, "dbartsMixedMatrix")
expect_equal(dim(mm), c(n, 5L))
expect_equal(colnames(mm), c("x1", "f", "sv", "sm.1", "sm.2"))
expect_equal(mm$map, c(1L, 2L, -1L, -2L, -3L))
expect_equal(attr(mm, "varTypes"), c(0L, 1L, 0L, 0L, 0L))
expect_equal(attr(mm, "factor.levels")[[2L]], levels(f))
expect_inherits(mm$sparse, "dgCMatrix")
expect_equal(as.double(mm$sparse[, 1L]), as.double(sv))

# the indicators builder splices sparse columns among the dummy expansions
mm.ind <- makeModelMatrixFromDataFrame(x.frame)
expect_inherits(mm.ind, "dbartsMixedMatrix")
expect_equal(colnames(mm.ind), c("x1", "f.a", "f.b", "f.c", "sv", "sm.1", "sm.2"))
expect_equal(length(attr(mm.ind, "drop")), 4L)

# row subsetting keeps the container; column extraction densifies
mm.rows <- mm[1:100, , drop = FALSE]
expect_inherits(mm.rows, "dbartsMixedMatrix")
expect_equal(nrow(mm.rows), 100L)
expect_equal(mm[, 3L], as.double(sv))

# ingestion: the data object holds the container, the factor stays
# categorical
data.mixed <- dbartsData(x.frame, y)
expect_inherits(data.mixed@x, "dbartsMixedMatrix")
expect_equal(data.mixed@varTypes, c(0L, 1L, 0L, 0L, 0L))
data.subset <- dbartsData(x.frame, y, subset = 1:200)
expect_equal(length(data.subset@y), 200L)
expect_equal(nrow(data.subset@x), 200L)

# a mixed fit recovers the signal a fully dense fit of the same values
# finds; sigma estimates differ by construction (marginal fallback), so the
# comparison is on fit quality
fit.mixed <- bart2(x.frame, y, n.samples = 300L, n.burn = 100L,
                   n.trees = 50L, n.chains = 1L, n.threads = 1L,
                   verbose = FALSE)
x.dense.equiv <- as.matrix(data.mixed@x)
fit.dense <- bart2(x.dense.equiv, y, n.samples = 300L, n.burn = 100L,
                   n.trees = 50L, n.chains = 1L, n.threads = 1L,
                   verbose = FALSE)
sse.mixed <- sum((fitted(fit.mixed) - f.true)^2)
sse.dense <- sum((fitted(fit.dense) - f.true)^2)
sse.mean <- sum((mean(y) - f.true)^2)
expect_true(sse.mixed < 0.25 * sse.mean)
expect_true(sse.mixed < 2.0 * sse.dense + 1.0)

# variable names ride through to varcount
expect_equal(colnames(fit.mixed$varcount), colnames(mm))

# the indicators surface composes too
fit.bart <- bart(x.frame, y, ndpost = 50L, nskip = 50L, ntree = 25L,
                 verbose = FALSE)
expect_equal(colnames(fit.bart$varcount), colnames(mm.ind))

# missing values in either representation route through MIA
sv.na <- sv
sv.na@x[1L] <- NA_real_
x.frame.na <- data.frame(x1 = x1, f = f)
x.frame.na$x1[2L] <- NA_real_
x.frame.na$sv <- sv.na
fit.na <- bart2(x.frame.na, y, n.samples = 20L, n.burn = 20L, n.trees = 25L,
                n.chains = 1L, n.threads = 1L, verbose = FALSE)
expect_equal(length(fitted(fit.na)), n)
expect_error(dbartsData(x.frame.na, y, missing = "error"),
             pattern = "missing values")

# test predictions expand a test frame to a dense matrix over all columns
x.test <- data.frame(x1 = x1[1:20], f = f[1:20])
x.test$sv <- as.double(sv)[1:20]
x.test$sm <- sm.dense[1:20, , drop = FALSE]
fit.test <- bart2(x.frame, y, test = x.test, n.samples = 50L, n.burn = 50L,
                  n.trees = 25L, n.chains = 1L, n.threads = 1L,
                  verbose = FALSE)
expect_equal(length(fit.test$yhat.test.mean), 20L)

# linear leaves designate dense-backed columns; sparse-backed ones refuse
control <- dbartsControl(n.samples = 10L, n.burn = 0L, n.trees = 25L,
                         updateState = FALSE)
sampler.linear <- dbarts(x.frame, y, node.prior = linear("x1"),
                         control = control)
run.linear <- sampler.linear$run()
expect_true(all(is.finite(run.linear$sigma)))
expect_error(dbarts(x.frame, y, node.prior = linear("sv")),
             pattern = "sparse-backed")

# the raw-x surface and whole-data replacement are fixed at creation;
# grouped rbart_vi is reserved
sampler <- dbarts(x.frame, y, control = control)
invisible(sampler$run())
expect_silent(sampler$setResponse(y))
expect_error(sampler$setPredictor(x.dense.equiv),
             pattern = "sparse predictors fix the design")
expect_error(sampler$setData(dbartsData(x.frame, y)),
             pattern = "sparse")
expect_error(rbart_vi(x.frame, y, group.by = rep_len(1:4, n),
                      n.samples = 5L, n.burn = 5L, n.trees = 25L,
                      n.chains = 1L, n.threads = 1L),
             pattern = "sparse")

# save/load: sampler re-creation from the stored container restores state
control.state <- dbartsControl(n.samples = 10L, n.burn = 0L, n.trees = 25L,
                               n.chains = 2L, n.threads = 1L,
                               updateState = TRUE)
sampler.state <- dbarts(x.frame, y, control = control.state)
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
xval <- xbart(x.frame, y, n.samples = 40L, n.burn = c(20L, 5L, 5L),
              n.trees = 25L, n.reps = 2L, n.test = 4L, n.threads = 1L)
expect_true(all(is.finite(xval)))

# a frame of only sparse columns still builds (the container's dense part
# is empty)
x.all.sparse <- data.frame(row.names = seq_len(n))
x.all.sparse$sv <- sv
x.all.sparse$sm <- sm
data.all.sparse <- dbartsData(x.all.sparse, y)
expect_inherits(data.all.sparse@x, "dbartsMixedMatrix")
expect_true(is.null(data.all.sparse@x$dense))
fit.all.sparse <- bart2(x.all.sparse, y, n.samples = 20L, n.burn = 20L,
                        n.trees = 25L, n.chains = 1L, n.threads = 1L,
                        verbose = FALSE)
expect_equal(length(fitted(fit.all.sparse)), n)
