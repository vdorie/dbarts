# mixed dense and sparse predictors: a data frame holding numeric and
# factor columns alongside Matrix::sparseVector or dgCMatrix columns enters
# through the x/y interface, dense-backed columns keep categorical splits
# and linear leaves, and CSC-backed columns tier exactly as an all-sparse
# design does (docs/design/sparse-columns.md)

if (!requireNamespace("Matrix", quietly = TRUE)) {
  exit_file("Matrix not available")
}

set.seed(42)
n <- 400L
x1 <- rnorm(n)
f <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
sv <- Matrix::sparseVector(
  x = 0.5 + runif(60L),
  i = sort(sample.int(n, 60L)),
  length = n
)
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
expect_equal(
  colnames(mm.ind),
  c("x1", "f.a", "f.b", "f.c", "sv", "sm.1", "sm.2")
)
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
fit.mixed <- bart2(
  x.frame,
  y,
  n.samples = 300L,
  n.burn = 100L,
  n.trees = 50L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
x.dense.equiv <- as.matrix(data.mixed@x)
fit.dense <- bart2(
  x.dense.equiv,
  y,
  n.samples = 300L,
  n.burn = 100L,
  n.trees = 50L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
sse.mixed <- sum((fitted(fit.mixed) - f.true)^2)
sse.dense <- sum((fitted(fit.dense) - f.true)^2)
sse.mean <- sum((mean(y) - f.true)^2)
expect_true(sse.mixed < 0.25 * sse.mean)
expect_true(sse.mixed < 2.0 * sse.dense + 1.0)

# variable names ride through to varcount
expect_equal(colnames(fit.mixed$varcount), colnames(mm))

# the indicators surface composes too
fit.bart <- bart(
  x.frame,
  y,
  ndpost = 50L,
  nskip = 50L,
  ntree = 25L,
  verbose = FALSE
)
expect_equal(colnames(fit.bart$varcount), colnames(mm.ind))

# missing values in either representation route through MIA
sv.na <- sv
sv.na@x[1L] <- NA_real_
x.frame.na <- data.frame(x1 = x1, f = f)
x.frame.na$x1[2L] <- NA_real_
x.frame.na$sv <- sv.na
fit.na <- bart2(
  x.frame.na,
  y,
  n.samples = 20L,
  n.burn = 20L,
  n.trees = 25L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_equal(length(fitted(fit.na)), n)
expect_error(
  dbartsData(x.frame.na, y, missing = "error"),
  pattern = "missing values"
)

# test predictions expand a test frame to a dense matrix over all columns
x.test <- data.frame(x1 = x1[1:20], f = f[1:20])
x.test$sv <- as.double(sv)[1:20]
x.test$sm <- sm.dense[1:20, , drop = FALSE]
fit.test <- bart2(
  x.frame,
  y,
  test = x.test,
  n.samples = 50L,
  n.burn = 50L,
  n.trees = 25L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_equal(length(fit.test$yhat.test.mean), 20L)

# linear leaves designate dense-backed columns; sparse-backed ones refuse
control <- dbartsControl(
  n.samples = 10L,
  n.burn = 0L,
  n.trees = 25L,
  updateState = FALSE
)
sampler.linear <- dbarts(
  x.frame,
  y,
  node.prior = linear("x1"),
  control = control
)
run.linear <- sampler.linear$run()
expect_true(all(is.finite(run.linear$sigma)))
expect_error(
  dbarts(x.frame, y, node.prior = linear("sv")),
  pattern = "sparse-backed"
)

# a mixed dense/sparse container takes column-granular AND whole-matrix
# mutation, maintained R-side; whole-data replacement stays fixed at creation,
# and grouped rbart_vi is reserved
sampler <- dbarts(x.frame, y, control = control)
invisible(sampler$run())
expect_silent(sampler$setResponse(y))
expect_error(sampler$setData(dbartsData(x.frame, y)), pattern = "sparse")

# a dense-backed column, then a CSC-backed one (pattern-changing), then both
# kinds named in the SAME call; data@x stays a container throughout
x1.new <- x1 * 1.1 + 0.05
expect_silent(sampler$setPredictor(x1.new, column = 1L, forceUpdate = TRUE))
expect_inherits(sampler$data@x, "dbartsMixedMatrix")
expect_equal(as.matrix(sampler$data@x)[, 1L], x1.new)

sv.new <- as.double(sv)
sv.new[sv.new != 0] <- sv.new[sv.new != 0] * 1.5
sv.new[c(1L, 2L)] <- c(0.9, 1.1)
expect_silent(sampler$setPredictor(sv.new, column = 3L, forceUpdate = TRUE))
expect_equal(as.matrix(sampler$data@x)[, 3L], sv.new)
expect_equal(diff(sampler$data@x$sparse@p)[1L], sum(sv.new != 0))

sm1.new <- sm.dense[, 1L] * 2
expect_silent(
  sampler$setPredictor(
    cbind(x1, sm1.new),
    column = c(1L, 4L),
    forceUpdate = TRUE
  )
)
expect_equal(as.matrix(sampler$data@x)[, 1L], x1)
expect_equal(as.matrix(sampler$data@x)[, 4L], sm1.new)
expect_equal(diff(sampler$data@x$sparse@p)[2L], sum(sm1.new != 0))
# the untouched CSC column is bit-identical
expect_equal(as.matrix(sampler$data@x)[, 5L], sm.dense[, 2L])

# per-observation updates run on a DENSE-backed column of a mixed design (the
# IRT latent case) and stay refused on a CSC-backed one
installed <- sampler$setPredictor(x1.new, column = 1L, forceUpdate = "partial")
expect_true(is.logical(installed) && length(installed) == n)
expect_equal(as.matrix(sampler$data@x)[installed, 1L], x1.new[installed])
expect_error(
  sampler$setPredictor(sv.new, column = 3L, forceUpdate = "partial"),
  pattern = "per-observation updates require a dense-backed column"
)

# whole-matrix replacement of a mixed design: the container is spliced R-side
# and STAYS a container, dense-backed and CSC-backed columns alike, with a
# replaced sparse column densifying its storage
x.whole <- x.dense.equiv
x.whole[, 1L] <- x.whole[, 1L] * 0.9 + 0.02
x.whole[, 3L] <- x.whole[, 3L] + 0.25
expect_silent(sampler$setPredictor(x.whole, forceUpdate = TRUE))
expect_inherits(sampler$data@x, "dbartsMixedMatrix")
expect_equal(unname(as.matrix(sampler$data@x)), unname(x.whole))
expect_equal(diff(sampler$data@x$sparse@p)[1L], n)
# the categorical column keeps its declared levels through the splice
expect_equal(attr(sampler$data@x, "factor.levels")[[2L]], levels(f))
# a rejected transactional whole-matrix replacement leaves data@x untouched
x.before <- as.matrix(sampler$data@x)
accepted.whole <- sampler$setPredictor(
  matrix(0, n, ncol(x.dense.equiv)),
  forceUpdate = FALSE
)
expect_false(isTRUE(accepted.whole))
expect_equal(as.matrix(sampler$data@x), x.before)

# the mirror gate: re-creating the sampler from the mutated container (the
# transparent rebuild after save/load) continues the chains from the same
# state, so the R-side splice describes the design the engine now holds
sampler.mirror <- dbarts(
  x.frame,
  y,
  control = dbartsControl(
    n.samples = 10L,
    n.burn = 0L,
    n.trees = 25L,
    n.chains = 1L,
    n.threads = 1L,
    updateState = TRUE
  )
)
invisible(sampler.mirror$run(numBurnIn = 20L, numSamples = 5L))
expect_silent(sampler.mirror$setPredictor(x.whole, forceUpdate = TRUE))
sampler.mirror$storeState()
mirrorFile <- tempfile(fileext = ".rds")
saveRDS(sampler.mirror, mirrorFile)
sampler.recreated <- readRDS(mirrorFile)
file.remove(mirrorFile)
set.seed(313)
run.mutated <- sampler.mirror$run(numBurnIn = 0L, numSamples = 10L)
set.seed(313)
run.recreated <- sampler.recreated$run(numBurnIn = 0L, numSamples = 10L)
expect_equal(run.mutated$sigma, run.recreated$sigma)
expect_equal(run.mutated$train, run.recreated$train)

# ALIGNMENT: a container assembled from another data set carries that set's
# level order, so its DENSE-backed factor column is re-coded against the
# training table at the validateXTest funnel - it then predicts exactly as the
# equivalent data frame does, and a genuinely new level refuses by name
alignControl <- dbartsControl(
  n.samples = 10L,
  n.burn = 0L,
  n.trees = 25L,
  n.chains = 1L,
  n.threads = 1L,
  updateState = FALSE
)
set.seed(515L)
n.align <- 15L
f.align <- factor(
  sample(c("a", "b", "c"), n.align, replace = TRUE),
  levels = c("c", "a", "b")
)
sm.align.dense <- matrix(0, n.align, 2L)
sm.align.dense[c(1L, 5L), 1L] <- c(0.7, 0.3)
sm.align.dense[3L, 2L] <- 0.4
x1.align <- rnorm(n.align)
sv.align <- numeric(n.align)
sv.align[2L] <- 0.6

frame.align <- data.frame(x1 = x1.align, f = f.align)
frame.align$sv <- Matrix::sparseVector(x = 0.6, i = 2L, length = n.align)
frame.align$sm <- methods::as(sm.align.dense, "CsparseMatrix")
container.align <- dbartsData(frame.align, rnorm(n.align))@x
expect_equal(levels(container.align$dense[[2L]]), c("c", "a", "b"))

data.align <- dbartsData(x.frame, y, test = container.align)
expect_equal(
  as.matrix(data.align@x.test)[, 2L],
  as.double(match(as.character(f.align), levels(f)) - 1L)
)

frame.align.dense <- data.frame(x1 = x1.align, f = f.align)
frame.align.dense$sv <- sv.align
frame.align.dense$sm <- sm.align.dense
alignmentGate <- function(test, seed) {
  set.seed(seed)
  dbarts(x.frame, y, test = test, control = alignControl)$run()
}
expect_identical(
  alignmentGate(container.align, 616L)$test,
  alignmentGate(frame.align.dense, 616L)$test
)

frame.align.new <- frame.align
frame.align.new$f <- factor(c("z", rep("a", n.align - 1L)))
expect_error(
  dbartsData(
    x.frame,
    y,
    test = dbartsData(frame.align.new, rnorm(n.align))@x
  ),
  pattern = "has levels not present in the"
)
expect_error(
  rbart_vi(
    x.frame,
    y,
    group.by = rep_len(1:4, n),
    n.samples = 5L,
    n.burn = 5L,
    n.trees = 25L,
    n.chains = 1L,
    n.threads = 1L
  ),
  pattern = "sparse"
)

# save/load: sampler re-creation from the stored container restores state
control.state <- dbartsControl(
  n.samples = 10L,
  n.burn = 0L,
  n.trees = 25L,
  n.chains = 2L,
  n.threads = 1L,
  updateState = TRUE
)
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
xval <- xbart(
  x.frame,
  y,
  n.samples = 40L,
  n.burn = c(20L, 5L, 5L),
  n.trees = 25L,
  n.reps = 2L,
  n.test = 4L,
  n.threads = 1L
)
expect_true(all(is.finite(xval)))

# a frame of only sparse columns still builds (the container's dense part
# is empty)
x.all.sparse <- data.frame(row.names = seq_len(n))
x.all.sparse$sv <- sv
x.all.sparse$sm <- sm
data.all.sparse <- dbartsData(x.all.sparse, y)
expect_inherits(data.all.sparse@x, "dbartsMixedMatrix")
expect_true(is.null(data.all.sparse@x$dense))
fit.all.sparse <- bart2(
  x.all.sparse,
  y,
  n.samples = 20L,
  n.burn = 20L,
  n.trees = 25L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_equal(length(fitted(fit.all.sparse)), n)

# METADATA: a container whose per-sparse-column reference metadata does not
# describe its sparse block is refused where it arrives - creation validates
# the pair exactly as the mutation entrances do, rather than ignoring it and
# leaving the same object to be refused at the next entrance
data.short <- data.mixed
data.short@x$sparseReference <- data.short@x$sparseReference[1L]
expect_error(
  dbarts(data.short, control = alignControl),
  pattern = "malformed mixed predictor container"
)
data.no.count <- data.mixed
data.no.count@x$sparseCategoryCount <- NULL
expect_error(
  dbarts(data.no.count, control = alignControl),
  pattern = "malformed mixed predictor container"
)
# a map entry that names no block at all - NA included, which is the negative
# extreme of the integer type the sparse arm reads - is refused the same way
data.na.map <- data.mixed
data.na.map@x$map[1L] <- NA_integer_
expect_error(
  dbarts(data.na.map, control = alignControl),
  pattern = "malformed mixed predictor container"
)
data.wild.map <- data.mixed
data.wild.map@x$map[1L] <- -100L
expect_error(
  dbarts(data.wild.map, control = alignControl),
  pattern = "malformed mixed predictor container"
)

sampler.meta <- dbarts(x.frame, y, control = alignControl)
expect_error(
  sampler.meta$setPredictor(data.short@x, forceUpdate = TRUE),
  pattern = "malformed mixed predictor container"
)
# the well-formed container of the same design still installs and runs
expect_silent(sampler.meta$setPredictor(data.mixed@x, forceUpdate = TRUE))
expect_true(all(is.finite(sampler.meta$run()$sigma)))
