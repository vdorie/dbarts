# a sparse-valued mutation argument - any sparse Matrix class, a sparseVector,
# or a mixed container - lands on a sparse or mixed design without the caller
# densifying. The bridge hands its sparse columns to the engine as stored
# entries, read under the STORE's implicit rule, so every accepted shape must
# land BITWISE what its dense equivalent lands: fits, cut grid, predictions

source(
  system.file("common", "strictData.R", package = "dbarts"),
  local = TRUE
)
if (!requireNamespace("Matrix", quietly = TRUE)) {
  exit_file("Matrix not available")
}

control <- dbartsControl(
  n.samples = 5L,
  n.burn = 0L,
  n.trees = 20L,
  n.chains = 1L,
  n.threads = 1L,
  updateState = FALSE
)

## A sampler over `design` run to a fixed state, then mutated: the fits it
## reports afterwards and the source it maintains R-side are what the two
## spellings of one mutation must agree on.
mutatedState <- function(design, y, mutate, ctrl = control) {
  set.seed(91L)
  sampler <- dbarts(design, y, control = ctrl)
  invisible(sampler$run(30L, 5L))
  result <- mutate(sampler)
  sampler$storeState()
  after <- sampler$run(0L, 5L)
  list(
    result = result,
    ssr = sampler$getSumsOfSquaredResiduals(),
    train = after$train,
    sigma = after$sigma,
    cutPoints = attr(sampler$state, "cutPoints"),
    predicted = sampler$predict(design),
    x = as.matrix(sampler$data@x),
    class = class(sampler$data@x)
  )
}

## The sparse argument and its dense equivalent, against the same design.
expectTwinsAgree <- function(
  design,
  y,
  sparse.arg,
  dense.arg,
  column = NULL,
  forceUpdate = TRUE,
  updateCutPoints = FALSE,
  ctrl = control
) {
  mutate <- function(value) {
    function(sampler) {
      if (is.null(column)) {
        sampler$setPredictor(
          value,
          forceUpdate = forceUpdate,
          updateCutPoints = updateCutPoints
        )
      } else {
        sampler$setPredictor(
          value,
          column = column,
          forceUpdate = forceUpdate,
          updateCutPoints = updateCutPoints
        )
      }
    }
  }
  sparse <- mutatedState(design, y, mutate(sparse.arg), ctrl)
  # nolint next: object_usage_linter. tinytest attaches expect_* at run time.
  expect_identical(sparse, mutatedState(design, y, mutate(dense.arg), ctrl))
  invisible(sparse)
}

sparseBlock <- function(values) {
  methods::as(values, "CsparseMatrix")
}

set.seed(9081)
n <- 240L

# DESIGN A: a bare dgCMatrix, every column CSC-backed and ordinal
p.a <- 4L
a.values <- matrix(0, n, p.a)
for (j in seq_len(p.a)) {
  nonzero <- runif(n) < 0.2
  a.values[nonzero, j] <- 0.5 + runif(sum(nonzero))
}
x.a <- sparseBlock(a.values)
y.a <- as.double(a.values %*% c(1.5, -1, 0.5, 2)) + rnorm(n, 0, 0.3)

# DESIGN B: a mixed container - a dense numeric column, a dense factor, an
# ordinal sparse column, and a CSC-backed categorical column
levels.f <- c("s1", "s2", "s3")
b.frame <- data.frame(
  u = rnorm(n),
  g = factor(sample(c("a", "b", "c"), n, replace = TRUE))
)
b.frame$s <- Matrix::sparseVector(
  x = 0.5 + runif(48L),
  i = sort(sample.int(n, 48L)),
  length = n
)
b.frame$f <- sparseFactor(
  factor(sample(levels.f, n, replace = TRUE), levels = levels.f),
  reference = "s2"
)
y.b <- b.frame$u - as.double(b.frame$s) + rnorm(n, 0, 0.3)
expect_inherits(dbartsData(b.frame, y.b)@x, "dbartsMixedMatrix")

# the replacement values, drawn once so both spellings carry the same doubles
newColumn <- function(density) {
  column <- rep(0, n)
  nonzero <- runif(n) < density
  column[nonzero] <- 0.5 + runif(sum(nonzero))
  column
}
a.new <- matrix(0, n, p.a)
for (j in seq_len(p.a)) {
  a.new[, j] <- newColumn(0.25)
}

# --- DESIGN A x {whole, single, multi} x {dgCMatrix, sparseVector, container}

expectTwinsAgree(x.a, y.a, sparseBlock(a.new), a.new)

# the bridge takes the argument itself: the R layer hands the sparse object
# over rather than densifying it, so the agreement above is the engine's
# entry-wise install and not a coincidence of two dense calls
set.seed(91L)
sampler.bridge <- dbarts(x.a, y.a, control = control)
expect_true(.Call(
  dbarts:::C_dbarts_bartcore_setPredictor,
  sampler.bridge$getPointer(),
  sparseBlock(a.new),
  TRUE,
  FALSE
))

a.one <- sparseBlock(a.new[, 2L, drop = FALSE])
expectTwinsAgree(x.a, y.a, a.one, as.double(a.new[, 2L]), column = 2L)

a.two <- sparseBlock(a.new[, c(1L, 3L)])
expectTwinsAgree(x.a, y.a, a.two, a.new[, c(1L, 3L)], column = c(1L, 3L))

a.vector <- Matrix::sparseVector(
  x = a.new[a.new[, 4L] != 0, 4L],
  i = which(a.new[, 4L] != 0),
  length = n
)
expectTwinsAgree(x.a, y.a, a.vector, a.new[, 4L], column = 4L)

a.container <- dbarts:::wrapSparseTestMatrix(sparseBlock(a.new))
expectTwinsAgree(x.a, y.a, a.container, a.new)

# --- DESIGN B x {whole, single, multi} x {dgCMatrix, sparseVector, container}

b.new <- data.frame(
  u = rnorm(n),
  g = factor(
    sample(c("a", "b", "c"), n, replace = TRUE),
    levels = c("a", "b", "c")
  )
)
b.new$s <- Matrix::sparseVector(
  x = 0.5 + runif(60L),
  i = sort(sample.int(n, 60L)),
  length = n
)
b.new$f <- sparseFactor(
  factor(sample(levels.f, n, replace = TRUE), levels = levels.f),
  reference = "s2"
)
b.container <- dbarts:::makeCategoricalModelMatrix(b.new)
expectTwinsAgree(b.frame, y.b, b.container, as.matrix(b.container))

# a sparse argument onto the SPARSE-backed ordinal column (3), by index
b.sparse.col <- sparseBlock(matrix(a.new[, 1L], n, 1L))
expectTwinsAgree(b.frame, y.b, b.sparse.col, a.new[, 1L], column = 3L)

# a sparse argument onto a DENSE-backed column (1) of the same container: the
# dense kernels run, so the argument's storage is purely an ingestion detail
expectTwinsAgree(b.frame, y.b, b.sparse.col, a.new[, 1L], column = 1L)

# by name, spanning a dense-backed and a sparse-backed column at once
b.two <- sparseBlock(a.new[, c(2L, 3L)])
expectTwinsAgree(b.frame, y.b, b.two, a.new[, c(2L, 3L)], column = c("u", "s"))

b.vector <- Matrix::sparseVector(
  x = a.new[a.new[, 2L] != 0, 2L],
  i = which(a.new[, 2L] != 0),
  length = n
)
expectTwinsAgree(b.frame, y.b, b.vector, a.new[, 2L], column = 3L)

# --- explicit stored zeros: they materialize to the implicit, and the R-side
# splice must canonicalize them away or data@x's pattern diverges from the
# store's
zero.column <- a.new[, 1L]
zero.column[c(3L, 11L, 40L)] <- 0
explicit <- Matrix::sparseMatrix(
  i = seq_len(n),
  j = rep(1L, n),
  x = zero.column,
  dims = c(n, 1L),
  repr = "C"
)
expect_true(any(explicit@x == 0))
expectTwinsAgree(b.frame, y.b, explicit, zero.column, column = 3L)

set.seed(91L)
sampler.zeros <- dbarts(b.frame, y.b, control = control)
sampler.zeros$setPredictor(explicit, column = 3L, forceUpdate = TRUE)
stored.zeros <- sampler.zeros$data@x$sparse
kept <- seq.int(
  stored.zeros@p[1L] + 1L,
  length.out = stored.zeros@p[2L] - stored.zeros@p[1L]
)
expect_false(any(stored.zeros@x[kept] == 0))

# --- an all-implicit sparse column materializes constant, and a constant
# column is degenerate: the transaction rejects it exactly as it rejects the
# dense equivalent, and the refused update leaves data@x byte-untouched
empty <- Matrix::sparseMatrix(
  i = integer(0),
  j = integer(0),
  x = double(0),
  dims = c(n, 1L),
  repr = "C"
)
expectTwinsAgree(x.a, y.a, empty, rep(0, n), column = 1L)

set.seed(91L)
sampler.degenerate <- dbarts(x.a, y.a, control = control)
invisible(sampler.degenerate$run(30L, 5L))
before.degenerate <- sampler.degenerate$data@x
expect_false(
  sampler.degenerate$setPredictor(empty, column = 1L, forceUpdate = FALSE)
)
expect_identical(sampler.degenerate$data@x, before.degenerate)

# --- a container declaring a reference level against a store-ORDINAL column is
# malformed: the container reads its implicit rows as the reference, the engine
# reads them as zero
bad.frame <- data.frame(z = rep(0, n))
bad.frame$z <- sparseFactor(
  factor(sample(levels.f, n, replace = TRUE), levels = levels.f),
  reference = "s2"
)
bad.container <- dbarts:::makeCategoricalModelMatrix(bad.frame)
expect_false(is.na(bad.container$sparseReference[1L]))
set.seed(91L)
sampler.reference <- dbarts(b.frame, y.b, control = control)
expect_error(
  sampler.reference$setPredictor(
    bad.container,
    column = 3L,
    forceUpdate = TRUE
  ),
  pattern = "reference level"
)

# --- a categorical column takes its labels, so a sparse matrix of codes - here
# one past the store's category count - is refused by name
over.k <- Matrix::sparseMatrix(
  i = c(2L, 7L),
  j = c(1L, 1L),
  x = c(1, 9),
  dims = c(n, 1L),
  repr = "C"
)
expect_error(
  sampler.reference$setPredictor(over.k, column = 4L, forceUpdate = TRUE),
  pattern = "column 'f' is categorical"
)

# --- an NA in a sparse argument is a stored NaN, and trips a data object
# whose missing policy is "error"
na.column <- Matrix::sparseMatrix(
  i = c(4L, 9L),
  j = c(1L, 1L),
  x = c(NA_real_, 0.7),
  dims = c(n, 1L),
  repr = "C"
)
sampler.strict <- dbarts(strictData(x.a, y.a), control = control)
expect_error(
  sampler.strict$setPredictor(na.column, column = 1L, forceUpdate = TRUE),
  pattern = "missing values"
)
# with missingness admitted the NA lands where its dense equivalent lands
expectTwinsAgree(x.a, y.a, na.column, as.double(na.column), column = 1L)

# --- mis-shaped arguments are refused, never reinterpreted: a sparseVector
# carries no column layout, so it can only fill one column's worth of rows
expect_error(
  sampler.strict$setPredictor(a.vector, column = c(1L, 2L), forceUpdate = TRUE),
  pattern = "'x' must have length"
)
expect_error(
  sampler.strict$setPredictor(a.one, column = integer(0), forceUpdate = TRUE),
  pattern = "."
)

# --- the re-quantize surface after a sparse-valued mutation: the store keeps
# its retained slices, so setCutPoints and setState still run
set.seed(91L)
sampler.after <- dbarts(x.a, y.a, control = control)
invisible(sampler.after$run(30L, 5L))
sampler.after$setPredictor(sparseBlock(a.new), forceUpdate = TRUE)
expect_null(dbarts:::rawPredictorMatrix(sampler.after$data@x))
ssr.after <- sampler.after$getSumsOfSquaredResiduals()
sampler.after$storeState()
stored.state <- sampler.after$state
sampler.after$setCutPoints(list(c(0.1, 0.5, 0.9)), 1L)
sampler.after$setState(stored.state)
expect_equal(sampler.after$getSumsOfSquaredResiduals(), ssr.after)

# --- the mirror gate: a sampler re-created from the mutated source (the
# transparent rebuild after save/load) continues the chain from the same state,
# so the R-side splice describes the design the engine now holds
mirror.control <- dbartsControl(
  n.samples = 5L,
  n.burn = 0L,
  n.trees = 20L,
  n.chains = 1L,
  n.threads = 1L,
  updateState = TRUE
)
set.seed(91L)
sampler.mirror <- dbarts(b.frame, y.b, control = mirror.control)
invisible(sampler.mirror$run(20L, 5L))
sampler.mirror$setPredictor(b.container, forceUpdate = TRUE)
sampler.mirror$storeState()
mirrorFile <- tempfile(fileext = ".rds")
saveRDS(sampler.mirror, mirrorFile)
sampler.recreated <- readRDS(mirrorFile)
file.remove(mirrorFile)
set.seed(313)
run.mutated <- sampler.mirror$run(0L, 10L)
set.seed(313)
run.recreated <- sampler.recreated$run(0L, 10L)
expect_equal(run.mutated$train, run.recreated$train)
expect_equal(run.mutated$sigma, run.recreated$sigma)

# --- a cut refresh reads the replacement's stored entries: uniform and
# quantile grids, whole and by column, land the dense twin's grid
quantile.control <- control
quantile.control@useQuantiles <- TRUE
# denser than the design, so a quantile grid keeps its count
a.fine <- vapply(seq_len(p.a), function(j) newColumn(0.5), double(n))
b.cut <- cbind(rnorm(n), a.fine[, 3L])
for (ctrl in list(control, quantile.control)) {
  expectTwinsAgree(
    x.a,
    y.a,
    sparseBlock(a.fine),
    a.fine,
    updateCutPoints = TRUE,
    ctrl = ctrl
  )
  expectTwinsAgree(
    b.frame,
    y.b,
    sparseBlock(b.cut),
    b.cut,
    column = c(1L, 3L),
    updateCutPoints = TRUE,
    ctrl = ctrl
  )
}

# --- transactional: an accepted replacement and a rolled-back one
a.near <- a.new
a.near[a.near != 0] <- a.near[a.near != 0] * 1.01
accepted <- expectTwinsAgree(
  x.a,
  y.a,
  sparseBlock(a.near),
  a.near,
  forceUpdate = FALSE
)
expect_true(accepted$result)
rolled <- expectTwinsAgree(
  x.a,
  y.a,
  empty,
  rep(0, n),
  column = 1L,
  forceUpdate = FALSE
)
expect_false(rolled$result)

# --- the sparseFactor column (4), whose reference is "s2", takes its labels:
# as a sparseFactor, whose absent rows read its own reference, or as a plain
# character vector, installing the same store
f.labels <- levels.f[sample(1:3, n, replace = TRUE)]
expectTwinsAgree(
  b.frame,
  y.b,
  sparseFactor(f.labels, levels = levels.f, reference = "s1"),
  f.labels,
  column = 4L
)

# --- a refused middle column leaves data@x and the next sweeps those of an
# untouched twin: codes for a categorical column and a quantile precheck
# (errors), and an all-absent column under forceUpdate = FALSE (FALSE, rolled
# back)
expectRefusalLeavesTwin <- function(
  design,
  y,
  value,
  columns,
  updateCutPoints,
  pattern,
  ctrl = control
) {
  runAfter <- function(mutate) {
    set.seed(91L)
    sampler <- dbarts(design, y, control = ctrl)
    invisible(sampler$run(30L, 5L))
    before <- sampler$data@x
    outcome <- mutate(sampler)
    list(outcome, identical(sampler$data@x, before), sampler$run(0L, 5L))
  }
  twin <- runAfter(function(sampler) FALSE)
  # each entry of columns is one call's column argument, NULL for the whole
  for (col in columns) {
    refused <- runAfter(function(sampler) {
      call <- function() {
        if (is.null(col)) {
          sampler$setPredictor(
            value,
            forceUpdate = FALSE,
            updateCutPoints = updateCutPoints
          )
        } else {
          sampler$setPredictor(
            value,
            column = col,
            forceUpdate = FALSE,
            updateCutPoints = updateCutPoints
          )
        }
      }
      if (is.null(pattern)) {
        return(call())
      }
      expect_error(call(), pattern = pattern)
      FALSE
    })
    expect_identical(refused[1:2], twin[1:2])
    # a refusal writes nothing, so the sweeps are the twin's bit for bit; a
    # rollback restores the store and trees exactly but rebuilds the cached
    # fits, whose sums round in another order
    if (is.null(pattern)) {
      expect_equal(refused[[3L]], twin[[3L]])
    } else {
      expect_identical(refused[[3L]], twin[[3L]])
    }
  }
}
three <- a.new[, 1:3]
three.empty <- three
three.empty[, 2L] <- 0
x.three <- sparseBlock(a.values[, 1:3])
expectRefusalLeavesTwin(
  x.three,
  y.a,
  sparseBlock(three.empty),
  list(NULL, 1:3),
  TRUE,
  "induced cut points",
  quantile.control
)
expectRefusalLeavesTwin(
  x.three,
  y.a,
  sparseBlock(three.empty),
  list(NULL, 1:3),
  FALSE,
  NULL
)
off.table <- sparseBlock(cbind(a.new[, 3L], c(9, rep(0, n - 1L)), a.new[, 2L]))
expectRefusalLeavesTwin(
  b.frame,
  y.b,
  off.table,
  list(c(3L, 4L, 1L)),
  FALSE,
  "column 'f' is categorical"
)

# --- a coded container is no label vector, so it is refused on the
# categorical column as numbers are
wide.levels <- paste0("s", 1:6)
wide.frame <- data.frame(z = rep(0, n))
wide.frame$z <- sparseFactor(
  factor(
    sample(c("s1", "s2", "s6"), n, replace = TRUE, prob = c(1, 1, 8)),
    levels = wide.levels
  ),
  reference = "s6"
)
wide.container <- dbarts:::makeCategoricalModelMatrix(wide.frame)
expect_error(
  sampler.reference$setPredictor(
    wide.container,
    column = 4L,
    forceUpdate = TRUE
  ),
  pattern = "column 'f' is categorical"
)

# --- a column named twice and rolled back unwinds exactly, on a sparse and a
# dense design
for (design in list(x.a, a.values)) {
  expectRefusalLeavesTwin(
    design,
    y.a,
    sparseBlock(cbind(a.near[, 1L], 0)),
    list(c(1L, 1L)),
    FALSE,
    NULL
  )
}

# --- any sparse Matrix class is the dgCMatrix it holds, never densified
a.pattern <- a.new != 0
other.classes <- list(
  methods::as(sparseBlock(a.new), "TsparseMatrix"),
  methods::as(sparseBlock(a.new), "RsparseMatrix"),
  methods::as(sparseBlock(a.pattern), "lMatrix"),
  methods::as(methods::as(sparseBlock(a.pattern), "nMatrix"), "CsparseMatrix")
)
expected.classes <- c("dgTMatrix", "dgRMatrix", "lgCMatrix", "ngCMatrix")
for (k in seq_along(other.classes)) {
  arg <- other.classes[[k]]
  expect_inherits(arg, expected.classes[k])
  twin <- if (k <= 2L) a.new else a.pattern + 0
  whole <- expectTwinsAgree(x.a, y.a, arg, sparseBlock(twin))
  expect_true(whole$class == "dgCMatrix")
  expectTwinsAgree(
    b.frame,
    y.b,
    arg[, 2:3],
    sparseBlock(twin[, 2:3]),
    column = c(1L, 3L)
  )
}
# a test set of another class stays resident and predicts as its dgCMatrix
testPredictions <- function(x.test) {
  set.seed(91L)
  sampler <- dbarts(x.a, y.a, control = control)
  sampler$setTestPredictor(x.test)
  list(
    class = class(sampler$data@x.test),
    test = sampler$run(10L, 5L)$test,
    predicted = sampler$predict(x.test)
  )
}
triplet.test <- testPredictions(other.classes[[1L]])
expect_true("dbartsMixedMatrix" %in% triplet.test$class)
expect_identical(triplet.test, testPredictions(sparseBlock(a.new)))

# --- signed zero: a 0/1 column replaced by zeros with one stored -0 under a
# quantile refresh stores the dense twin's cut, and it is +0
x.binary <- sparseBlock(cbind(
  as.double(seq_len(n) %% 3L == 0L),
  a.values[, 2L]
))
negative <- Matrix::sparseMatrix(
  i = 1L,
  j = 1L,
  x = -0,
  dims = c(n, 1L),
  repr = "C"
)
expect_true(1 / negative@x < 0)
signedCuts <- function(value) {
  mutatedState(
    x.binary,
    y.a,
    function(sampler) {
      sampler$setPredictor(
        value,
        column = 1L,
        forceUpdate = TRUE,
        updateCutPoints = TRUE
      )
    },
    quantile.control
  )$cutPoints[[1L]]
}
sparse.cuts <- signedCuts(negative)
expect_true(identical(
  sparse.cuts,
  signedCuts(c(-0, rep(0, n - 1L))),
  num.eq = FALSE
))
expect_true(all(sparse.cuts == 0 & 1 / sparse.cuts > 0))
