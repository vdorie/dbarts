# sparseFactor: the sparse unordered-factor wrapper - a high-cardinality
# unordered factor stored CSC over level codes with the reference level
# implicit. The categorical x/y path ingests it as one categorical column
# that bins bitwise-identically to a dense factor of the same values; the
# formula path and indicator expansion still refuse it.

set.seed(41)

# dense construction: non-reference entries become the stored ones,
# positions 0-based, level codes 1-based
f <- factor(c("a", "b", "a", "c", "a", "b"))
sf <- sparseFactor(f, reference = "a")
expect_inherits(sf, "sparseFactor")
expect_equal(sf@levels, c("a", "b", "c"))
expect_equal(sf@reference, "a")
expect_equal(sf@length, 6L)
expect_equal(sf@i, c(1L, 3L, 5L))
expect_equal(sf@values, c(2L, 3L, 2L))
expect_equal(length(sf), 6L)

# character input takes factor()'s level convention and the first level as
# the default reference
sf.chr <- sparseFactor(c("b", "a", "c", "a"))
expect_equal(sf.chr@levels, c("a", "b", "c"))
expect_equal(sf.chr@reference, "a")
expect_equal(sf.chr@i, c(0L, 2L))
expect_equal(sf.chr@values, c(2L, 3L))

# sparse construction: 1-based i in (the Matrix convention), stored
# 0-based and ascending; explicit reference-level entries drop
sf.sp <- sparseFactor(
  c("c", "b", "a"),
  levels = c("a", "b", "c"),
  reference = "a",
  i = c(4L, 2L, 6L),
  length = 8L
)
expect_equal(sf.sp@i, c(1L, 3L))
expect_equal(sf.sp@values, c(2L, 3L))
expect_equal(sf.sp@length, 8L)

# integer codes construct against an explicit level table
sf.int <- sparseFactor(
  c(2L, 3L),
  levels = c("a", "b", "c"),
  i = c(1L, 5L),
  length = 5L
)
expect_equal(sf.int@values, c(2L, 3L))
expect_error(sparseFactor(c(1L, 2L)), pattern = "requires explicit 'levels'")

# validation failures, one per rule
expect_error(
  sparseFactor(f, reference = "z"),
  pattern = "'reference' must be a single element of 'levels'"
)
expect_error(
  sparseFactor(c("a", "zzz"), levels = c("a", "b")),
  pattern = "absent from 'levels'"
)
expect_error(
  sparseFactor(c(0L, 1L), levels = c("a", "b"), i = c(1L, 2L), length = 2L),
  pattern = "level codes"
)
expect_error(
  sparseFactor(factor(c("a", NA))),
  pattern = "missing values are not supported"
)
expect_error(
  sparseFactor(c("b", "b"), levels = c("a", "b"), i = c(1L, 9L), length = 5L),
  pattern = "1-based positions"
)
expect_error(
  sparseFactor(c("b", "b"), levels = c("a", "b"), i = c(3L, 3L), length = 5L),
  pattern = "duplicated positions"
)
expect_error(
  sparseFactor(c("b", "b"), levels = c("a", "b"), i = c(1L, 2L)),
  pattern = "'length' is required"
)
expect_error(
  sparseFactor(c("b", "b"), levels = c("a", "b"), i = c(1L, 2.5), length = 5L),
  pattern = "integer positions"
)

# print shows length, stored-entry count, levels, and the reference
expect_stdout(print(sf), pattern = "sparseFactor of length 6")
expect_stdout(print(sf), pattern = "3 stored entries")
expect_stdout(print(sf), pattern = "levels: a, b, c")
expect_stdout(print(sf), pattern = "reference \\(implicit\\): a")

# narrowed refusals (the surfaces a sparseFactor still cannot enter; none of
# these reach the mixed-matrix assembly, so they need no Matrix)
n <- 40L
df <- data.frame(x1 = rnorm(n))
df$sf <- sparseFactor(
  sample(c("b", "c"), 8L, replace = TRUE),
  levels = c("a", "b", "c"),
  reference = "a",
  i = sort(sample.int(n, 8L)),
  length = n
)
y <- rnorm(n)

# the formula interface still refuses a sparseFactor (a bare S4 term would
# die inside model.frame; the pre-scan refuses cleanly ahead of it)
df.formula <- df
df.formula$y <- y
expect_error(
  dbartsData(y ~ x1 + sf, df.formula),
  pattern = "sparse predictors must be specified through the x/y interface"
)

# indicator expansion cannot dummy-code a sparse factor without densifying it
expect_error(
  dbartsData(df, y, factors = "indicators"),
  pattern = "sparse categorical predictors require"
)

# a bare sparseFactor as 'x' is unrecognized, exactly as a bare sparseVector is
expect_error(
  dbartsData(df$sf, y),
  pattern = "unrecognized 'formula' type"
)

if (!requireNamespace("Matrix", quietly = TRUE)) {
  exit_file("Matrix not available")
}

# ingestion of a frame mixing dense numeric, dense factor, sparseVector, and
# sparseFactor: the container types the two factors categorical and carries
# per-sparse-column reference/K metadata (0-based reference code, level count)
n <- 300L
set.seed(5L)
x1 <- rnorm(n)
fdense <- factor(sample(c("p", "q", "r"), n, replace = TRUE))
sv <- Matrix::sparseVector(
  x = 0.5 + runif(40L),
  i = sort(sample.int(n, 40L)),
  length = n
)
sfCodes <- sample.int(5L, n, replace = TRUE, prob = c(1, 10, 1, 1, 1))
sfCodes[1L] <- 5L
sf <- sparseFactor(
  factor(paste0("c", sfCodes), levels = paste0("c", 1:5)),
  reference = "c2"
)
frame <- data.frame(x1 = x1, fdense = fdense)
frame$sv <- sv
frame$sf <- sf
y <- 0.4 * sfCodes + rnorm(n)

data <- dbartsData(frame, y)
expect_inherits(data@x, "dbartsMixedMatrix")
expect_equal(data@varTypes, c(0L, 1L, 0L, 1L))
expect_equal(data@x$sparseCategoryCount, c(0L, 5L))
expect_equal(data@x$sparseReference, c(NA_integer_, 1L))
expect_equal(attr(data@x, "factor.levels")[[4L]], paste0("c", 1:5))

# a short fit over the mixed frame runs and stays finite
fit <- bart2(
  frame,
  y,
  n.samples = 40L,
  n.burn = 20L,
  n.trees = 25L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_true(all(is.finite(fitted(fit))))

# THE BITWISE GATE: a sparseFactor column's draws are byte-identical to the
# same values given as a dense factor. Both storage kinds take their category
# count from the DECLARED level table, so the top level is deliberately left
# unobserved - a genuine gap, which only a declared count can bin. The
# reference is a middle level (never first alphabetically), whose own
# level-order code the implicit rows carry.
runGateFit <- function(frame, y) {
  set.seed(2718L)
  sampler <- dbarts(
    frame,
    y,
    sigma = 1.0,
    control = dbartsControl(
      n.trees = 25L,
      n.chains = 1L,
      n.threads = 1L,
      updateState = FALSE
    )
  )
  sampler$run(20L, 40L)
}
bitwiseGate <- function(codes, levels, reference, seed) {
  set.seed(seed)
  codes[codes == length(levels)] <- length(levels) - 1L # the top level is a gap
  g <- factor(levels[codes], levels = levels)
  sf <- sparseFactor(g, reference = reference)
  m <- length(codes)
  x1 <- rnorm(m)
  y <- 0.3 * codes + rnorm(m)
  frame.dense <- data.frame(x1 = x1, g = g)
  frame.sparse <- data.frame(x1 = x1)
  frame.sparse$g <- sf
  dense <- runGateFit(frame.dense, y)
  sparse <- runGateFit(frame.sparse, y)
  list(
    sf = sf,
    sigma = identical(dense$sigma, sparse$sigma),
    train = identical(dense$train, sparse$train)
  )
}

# rank tier: a dominant reference keeps the stored fraction below the 0.2
# densification threshold
set.seed(20L)
rankCodes <- sample.int(6L, 400L, replace = TRUE, prob = c(1, 1, 50, 1, 1, 1))
rank <- bitwiseGate(rankCodes, paste0("L", 1:6), "L3", 31L)
expect_true(length(rank$sf@i) < 0.2 * 400L)
expect_true(rank$sigma)
expect_true(rank$train)

# densified tier: a rare reference pushes the stored fraction above 0.2, so
# the column scatters to dense categorical codes and runs the dense kernel
set.seed(21L)
denseCodes <- sample.int(6L, 400L, replace = TRUE, prob = c(3, 3, 1, 3, 3, 3))
densified <- bitwiseGate(denseCodes, paste0("L", 1:6), "L3", 32L)
expect_true(length(densified$sf@i) > 0.2 * 400L)
expect_true(densified$sigma)
expect_true(densified$train)

# pooled path: more than 63 levels exercises the wide-mask membership kernel
set.seed(22L)
pooledProb <- rep(1, 80L)
pooledProb[40L] <- 400
pooledCodes <- sample.int(80L, 400L, replace = TRUE, prob = pooledProb)
pooled <- bitwiseGate(pooledCodes, sprintf("V%02d", 1:80), "V40", 33L)
expect_true(length(pooled$sf@i) < 0.2 * 400L)
expect_true(pooled$sigma)
expect_true(pooled$train)

# sparse ordinal columns are unaffected (the S-CAT boundary): a sparseVector
# frame still ingests as ordinal
df.ord <- data.frame(x1 = rnorm(40L))
df.ord$sv <- Matrix::sparseVector(
  x = 0.5 + runif(5L),
  i = sort(sample.int(40L, 5L)),
  length = 40L
)
expect_inherits(dbartsData(df.ord, rnorm(40L)), "dbartsData")

# THE MATERIALIZER: as.matrix on a mixed container fills a sparse-
# categorical column's implicit rows with its reference code, not zero - an
# explicit entry can legitimately be code 0, so the fill has to precede the
# scatter, not follow it (mixedMatrix.R as.matrix.dbartsMixedMatrix)
n <- 10L
frame.mat <- data.frame(x1 = rnorm(n))
frame.mat$sf <- sparseFactor(
  c("a", "c"),
  levels = c("a", "b", "c"),
  reference = "b",
  i = c(1L, 6L),
  length = n
)
expectedCodes <- rep(1, n) # "b" is level 2, the implicit reference code
expectedCodes[1L] <- 0L # explicit "a" - the first level, code 0, stored anyway
expectedCodes[6L] <- 2L # explicit "c"
mm.mat <- dbarts:::makeCategoricalModelMatrix(frame.mat)
expect_equal(as.matrix(mm.mat)[, 2L], as.double(expectedCodes))

sampler.mat <- dbarts(
  frame.mat,
  rnorm(n),
  sigma = 1.0,
  control = dbartsControl(
    n.trees = 5L,
    n.chains = 1L,
    n.threads = 1L,
    updateState = FALSE
  )
)
expect_equal(extract(sampler.mat, "predictors")[, 2L], as.double(expectedCodes))

# a resident sparse test set assembles a mixed container, which needs Matrix
if (!requireNamespace("Matrix", quietly = TRUE)) {
  exit_file("Matrix not available")
}

# RESIDENT-SPARSE x.test: a sparseFactor x.test column recodes over the
# training level table and predicts exactly like a dense factor of the same
# values, held sparse (a mixed container) into the engine
set.seed(50L)
n <- 200L
x1.tr <- rnorm(n)
levels.g <- paste0("L", 1:4)
codes.tr <- sample.int(4L, n, replace = TRUE)
codes.tr[1L] <- 1L
g.tr <- factor(levels.g[codes.tr], levels = levels.g)
y.tr <- 0.3 * codes.tr + rnorm(n)
train.dense <- data.frame(x1 = x1.tr, g = g.tr)

n.test <- 25L
x1.te <- rnorm(n.test)
codes.te <- sample.int(4L, n.test, replace = TRUE)
codes.te[1L] <- 1L # the first level appears explicitly in the sparse test column
g.te <- factor(levels.g[codes.te], levels = levels.g)
test.dense <- data.frame(x1 = x1.te, g = g.te)
test.sparse <- data.frame(x1 = x1.te)
test.sparse$g <- sparseFactor(g.te, reference = "L2")

runTestGate <- function(test, seed) {
  set.seed(seed)
  sampler <- dbarts(
    train.dense,
    y.tr,
    sigma = 1.0,
    test = test,
    control = dbartsControl(
      n.trees = 20L,
      n.chains = 1L,
      n.threads = 1L,
      updateState = FALSE
    )
  )
  sampler$run(10L, 20L)
}
result.dense <- runTestGate(test.dense, 1234L)
result.sparse <- runTestGate(test.sparse, 1234L)
expect_identical(result.dense$test, result.sparse$test)

# a test level absent from training errors the same way for a sparse test
# column as it does for a dense one
unseenLabels <- c("L5", as.character(g.te[-1L]))
expect_error(
  dbarts(
    train.dense,
    y.tr,
    test = data.frame(x1 = x1.te, g = factor(unseenLabels)),
    control = dbartsControl(n.trees = 5L, updateState = FALSE)
  ),
  pattern = "has levels not present in the"
)
test.unseen <- data.frame(x1 = x1.te)
test.unseen$g <- sparseFactor(
  unseenLabels,
  levels = c(levels.g, "L5"),
  reference = "L2"
)
expect_error(
  dbarts(
    train.dense,
    y.tr,
    test = test.unseen,
    control = dbartsControl(n.trees = 5L, updateState = FALSE)
  ),
  pattern = "has levels not present in the"
)

# RESIDENCY: the sparse test set is held as a mixed container, not densified
set.seed(1234L)
sampler.resident <- dbarts(
  train.dense,
  y.tr,
  sigma = 1.0,
  test = test.sparse,
  control = dbartsControl(
    n.trees = 20L,
    n.chains = 1L,
    n.threads = 1L,
    updateState = FALSE
  )
)
expect_inherits(sampler.resident$data@x.test, "dbartsMixedMatrix")
expect_true(dbarts:::predictorSourceIsSparse(sampler.resident$data@x.test))

# MIXED RESIDENT-SPARSE x.test: a dense, a sparseFactor, and a sparse ordinal
# column ingest as one container and fit bitwise-identically to the same data
# supplied as a dense matrix
set.seed(60L)
n.mix <- 200L
a.tr <- rnorm(n.mix)
codes.tr.mix <- sample.int(4L, n.mix, replace = TRUE)
g.tr.mix <- factor(levels.g[codes.tr.mix], levels = levels.g)
s.tr <- rnorm(n.mix)
y.mix <- 0.4 * codes.tr.mix - 0.6 * a.tr + 0.5 * s.tr + rnorm(n.mix)
train.mix <- data.frame(a = a.tr, g = g.tr.mix, s = s.tr)

n.mtest <- 30L
a.te.mix <- rnorm(n.mtest)
codes.te.mix <- sample.int(4L, n.mtest, replace = TRUE)
g.te.mix <- factor(levels.g[codes.te.mix], levels = levels.g)
s.rows <- sort(sample.int(n.mtest, 8L))
s.vals <- 0.5 + runif(8L)
s.dense <- numeric(n.mtest) # implicit rows are numeric zero, as CSC densifies
s.dense[s.rows] <- s.vals

test.mix.dense <- data.frame(a = a.te.mix, g = g.te.mix, s = s.dense)
test.mix.sparse <- data.frame(a = a.te.mix)
test.mix.sparse$g <- sparseFactor(g.te.mix, reference = "L2")
test.mix.sparse$s <- Matrix::sparseVector(
  x = s.vals,
  i = s.rows,
  length = n.mtest
)

makeMixedSampler <- function(test, seed) {
  set.seed(seed)
  dbarts(
    train.mix,
    y.mix,
    sigma = 1.0,
    test = test,
    control = dbartsControl(
      n.trees = 20L,
      n.chains = 1L,
      n.threads = 1L,
      updateState = FALSE
    )
  )
}
sampler.mix.sparse <- makeMixedSampler(test.mix.sparse, 77L)
expect_inherits(sampler.mix.sparse$data@x.test, "dbartsMixedMatrix")
result.mix.sparse <- sampler.mix.sparse$run(10L, 20L)

sampler.mix.dense <- makeMixedSampler(test.mix.dense, 77L)
expect_true(is.matrix(sampler.mix.dense$data@x.test))
result.mix.dense <- sampler.mix.dense$run(10L, 20L)

expect_identical(result.mix.dense$test, result.mix.sparse$test)

# LEAF-COVARIATE REFUSAL: a linear leaf covariate that would land on a sparse
# test column is refused - sparse storage serves no dense raw test covariate
expect_error(
  dbarts(
    train.mix,
    y.mix,
    test = test.mix.sparse,
    node.prior = linear("s"),
    control = dbartsControl(
      n.trees = 5L,
      n.chains = 1L,
      n.threads = 1L,
      updateState = FALSE
    )
  ),
  pattern = "leaf covariate"
)

# SET-TEST MUTATION: setTestPredictor with a mixed/sparse container replaces the
# whole test store and fits identically to the dense-equivalent replacement
sampler.set.sparse <- makeMixedSampler(test.mix.dense, 88L)
sampler.set.sparse$setTestPredictor(test.mix.sparse)
expect_inherits(sampler.set.sparse$data@x.test, "dbartsMixedMatrix")
result.set.sparse <- sampler.set.sparse$run(10L, 20L)

sampler.set.dense <- makeMixedSampler(test.mix.dense, 88L)
sampler.set.dense$setTestPredictor(test.mix.dense)
expect_true(is.matrix(sampler.set.dense$data@x.test))
result.set.dense <- sampler.set.dense$run(10L, 20L)

expect_identical(result.set.dense$test, result.set.sparse$test)

# DENSE MUTATION preserved: per-column update by index and by name, plus NULL
# removal, against a dense test matrix
sampler.dense.mut <- makeMixedSampler(test.mix.dense, 99L)
new.a <- rnorm(n.mtest)
sampler.dense.mut$setTestPredictor(new.a, column = 1L)
expect_equal(as.numeric(sampler.dense.mut$data@x.test[, 1L]), new.a)
sampler.dense.mut$setTestPredictor(2 * new.a, column = "a")
expect_equal(as.numeric(sampler.dense.mut$data@x.test[, "a"]), 2 * new.a)
sampler.dense.mut$setTestPredictor(NULL)
expect_null(sampler.dense.mut$data@x.test)

# CONTAINER PER-COLUMN UPDATE: a container's per-column storage decision
# (dense vs CSC-backed) is preserved across a single-column replacement, and
# the result matches a whole-object install of the same spliced container AND
# the dense equivalent, by index and by name - for the dense-backed column
# ('a') and the sparse-backed one ('s')
new.a.container <- rnorm(n.mtest)
sampler.col.index <- makeMixedSampler(test.mix.sparse, 101L)
sampler.col.index$setTestPredictor(new.a.container, column = 1L)
expect_inherits(sampler.col.index$data@x.test, "dbartsMixedMatrix")
result.col.index <- sampler.col.index$run(10L, 20L)

sampler.col.name <- makeMixedSampler(test.mix.sparse, 101L)
sampler.col.name$setTestPredictor(new.a.container, column = "a")
result.col.name <- sampler.col.name$run(10L, 20L)
expect_identical(result.col.index$test, result.col.name$test)

test.mix.sparse.a <- test.mix.sparse
test.mix.sparse.a$a <- new.a.container
result.col.whole <- makeMixedSampler(test.mix.sparse.a, 101L)$run(10L, 20L)
expect_identical(result.col.index$test, result.col.whole$test)

test.mix.dense.a <- test.mix.dense
test.mix.dense.a$a <- new.a.container
result.col.dense <- makeMixedSampler(test.mix.dense.a, 101L)$run(10L, 20L)
expect_identical(result.col.index$test, result.col.dense$test)

new.s.container <- Matrix::sparseVector(
  x = 0.4 + runif(10L),
  i = sort(sample.int(n.mtest, 10L)),
  length = n.mtest
)
sampler.sparse.col.index <- makeMixedSampler(test.mix.sparse, 111L)
sampler.sparse.col.index$setTestPredictor(new.s.container, column = 3L)
result.sparse.col.index <- sampler.sparse.col.index$run(10L, 20L)

sampler.sparse.col.name <- makeMixedSampler(test.mix.sparse, 111L)
sampler.sparse.col.name$setTestPredictor(new.s.container, column = "s")
result.sparse.col.name <- sampler.sparse.col.name$run(10L, 20L)
expect_identical(result.sparse.col.index$test, result.sparse.col.name$test)

test.mix.sparse.s <- test.mix.sparse
test.mix.sparse.s$s <- new.s.container
result.sparse.col.whole <- makeMixedSampler(test.mix.sparse.s, 111L)$run(
  10L,
  20L
)
expect_identical(result.sparse.col.index$test, result.sparse.col.whole$test)

test.mix.dense.s <- test.mix.dense
test.mix.dense.s$s <- as.double(new.s.container)
result.sparse.col.dense <- makeMixedSampler(test.mix.dense.s, 111L)$run(
  10L,
  20L
)
expect_identical(result.sparse.col.index$test, result.sparse.col.dense$test)

# BRIDGE REFUSAL ROLLS BACK: a per-column value that plants an out-of-range
# category code past the training K is refused at the bridge exactly as a
# whole-object container is (the codeMessage above), and the refusal leaves
# @x.test/@offset.test as the prior, accepted container
sampler.col.refuse <- makeMixedSampler(test.mix.sparse, 121L)
before.col.refuse.x.test <- sampler.col.refuse$data@x.test
before.col.refuse.offset.test <- sampler.col.refuse$data@offset.test
bad.g.column <- Matrix::sparseVector(x = 9, i = 1L, length = n.mtest)
expect_error(
  sampler.col.refuse$setTestPredictor(bad.g.column, column = "g"),
  pattern = "existing category codes"
)
expect_identical(sampler.col.refuse$data@x.test, before.col.refuse.x.test)
expect_identical(
  sampler.col.refuse$data@offset.test,
  before.col.refuse.offset.test
)

# setTestPredictorAndOffset installs a container together with an offset,
# matching the dense-equivalent bitwise; a mismatched offset length still errors
off.mtest <- rnorm(n.mtest)
sampler.pao.sparse <- makeMixedSampler(test.mix.dense, 202L)
sampler.pao.sparse$setTestPredictorAndOffset(test.mix.sparse, off.mtest)
result.pao.sparse <- sampler.pao.sparse$run(10L, 20L)

sampler.pao.dense <- makeMixedSampler(test.mix.dense, 202L)
sampler.pao.dense$setTestPredictorAndOffset(test.mix.dense, off.mtest)
result.pao.dense <- sampler.pao.dense$run(10L, 20L)
expect_identical(result.pao.dense$test, result.pao.sparse$test)

expect_error(
  sampler.pao.sparse$setTestPredictorAndOffset(test.mix.sparse, off.mtest[-1L]),
  pattern = "'offset.test' must have the same number of rows as 'x.test'"
)

# LEAF-COVARIATE REFUSAL ON MUTATION: swapping in a container whose leaf
# covariate would be CSC-backed is refused, the prior dense test store stays
# intact, and the sampler still fits
sampler.leaf.mut <- dbarts(
  train.mix,
  y.mix,
  sigma = 1.0,
  test = test.mix.dense,
  node.prior = linear("s"),
  control = dbartsControl(
    n.trees = 10L,
    n.chains = 1L,
    n.threads = 1L,
    updateState = FALSE
  )
)
expect_error(
  sampler.leaf.mut$setTestPredictor(test.mix.sparse),
  pattern = "leaf covariate"
)
expect_true(is.matrix(sampler.leaf.mut$data@x.test))
result.leaf.mut <- sampler.leaf.mut$run(10L, 20L)
expect_false(anyNA(result.leaf.mut$test))

# CODE BOUNDS: every categorical ingestion entrance bounds codes against the
# TRAINING-side category count - the K its host declares whatever backs it,
# and the inferred max + 1 only where nothing declares one - whatever view
# carries them: a dense
# x.test matrix, a container's dense or CSC slice, or the reference code the
# container's implicit rows read. A code at or past that count mis-bins inline,
# shifts past a tree's category mask, or over-reads the pooled bitmap, so it is
# refused rather than clamped. The refusal lands at the sampler constructor;
# dbartsData never sees the training grid.
boundControl <- dbartsControl(
  n.trees = 10L,
  n.chains = 1L,
  n.threads = 1L,
  updateState = FALSE
)
sparseFrame <- function(labels, levels, reference) {
  frame <- data.frame(x1 = rnorm(length(labels)))
  frame$f <- sparseFactor(labels, levels = levels, reference = reference)
  frame
}

set.seed(70L)
levels.small <- c("a", "b", "c")
levels.big <- c("a", "b", "c", "d", "e")
n.bound <- 60L
labels.bound <- sample(levels.small, n.bound, replace = TRUE)
y.bound <- rnorm(n.bound) + match(labels.bound, levels.small)
train.bound <- sparseFrame(labels.bound, levels.small, "a")
sampler.bound <- dbarts(train.bound, y.bound, control = boundControl)

# two foreign containers: one whose STORED codes run past the training K = 3,
# one whose stored codes are in range but whose REFERENCE code - what every
# implicit row reads - is not
over.codes <- dbartsData(
  sparseFrame(c("a", "d", "e", "b"), levels.big, "a"),
  rnorm(4L)
)@x
over.reference <- dbartsData(
  sparseFrame(c("a", "b", "a", "b"), levels.big, "e"),
  rnorm(4L)
)@x
expect_equal(sort(over.codes$sparse@x), c(1, 3, 4))
expect_equal(over.reference$sparseReference, 4L)
expect_true(all(over.reference$sparse@x < 3))

codeMessage <- "categorical test predictors must hold existing category codes"
# A CONTAINER carries its own level table, so a code past the training count is
# a level the training data never saw: the alignment funnel (validateXTest)
# re-codes a foreign container against the training levels and refuses by name
# before the code bound ever sees it. The code bound stays the backstop for
# every view that carries no level table - a dense x.test matrix, and a
# container whose levels DO line up but whose codes were planted past them
# (below) - and for the training side.
levelMessage <- "has levels not present in the training data"
expect_error(
  dbarts(
    dbartsData(train.bound, y.bound, test = over.codes),
    control = boundControl
  ),
  pattern = levelMessage
)
expect_error(
  dbarts(
    dbartsData(train.bound, y.bound, test = over.reference),
    control = boundControl
  ),
  pattern = levelMessage
)
expect_error(sampler.bound$setTestPredictor(over.codes), pattern = levelMessage)
expect_error(
  sampler.bound$setTestPredictor(over.reference),
  pattern = levelMessage
)
expect_error(
  sampler.bound$setTestPredictorAndOffset(over.codes, rep(0, 4L)),
  pattern = levelMessage
)
expect_error(
  sampler.bound$setTestPredictorAndOffset(over.reference, rep(0, 4L)),
  pattern = levelMessage
)
# the refusals are inert: no test store was installed
expect_null(sampler.bound$data@x.test)

# a DENSE x.test view of the same CSC-trained column takes the declared K too
test.dense.over <- cbind(rnorm(3L), c(0, 1, 7))
colnames(test.dense.over) <- c("x1", "f")
expect_error(
  dbarts(
    dbartsData(train.bound, y.bound, test = test.dense.over),
    control = boundControl
  ),
  pattern = codeMessage
)

# a DENSE-trained categorical column bounds a CSC test container the same way
train.bound.dense <- data.frame(
  x1 = train.bound$x1,
  f = factor(labels.bound, levels = levels.small)
)
sampler.bound.dense <- dbarts(
  train.bound.dense,
  y.bound,
  control = boundControl
)
expect_error(
  sampler.bound.dense$setTestPredictor(over.codes),
  pattern = levelMessage
)
expect_error(
  dbarts(
    dbartsData(train.bound.dense, y.bound, test = over.reference),
    control = boundControl
  ),
  pattern = levelMessage
)

# POOLED (more than 63 levels, the wide-mask kernel): the same bound, where an
# unbounded code over-reads the membership bitmap rather than mis-binning
set.seed(71L)
levels.pooled <- sprintf("P%03d", 1:70)
prob.pooled <- rep(1, 70L)
prob.pooled[1L] <- 400
labels.pooled <- sample(levels.pooled, 200L, replace = TRUE, prob = prob.pooled)
y.pooled <- rnorm(200L)
train.pooled <- sparseFrame(labels.pooled, levels.pooled, "P001")
sampler.pooled <- dbarts(train.pooled, y.pooled, control = boundControl)
over.pooled <- dbartsData(
  sparseFrame(c("P002", "P150", "P200"), sprintf("P%03d", 1:220), "P001"),
  rnorm(3L)
)@x
expect_error(
  sampler.pooled$setTestPredictor(over.pooled),
  pattern = levelMessage
)
expect_error(
  dbarts(
    dbartsData(train.pooled, y.pooled, test = over.pooled),
    control = boundControl
  ),
  pattern = levelMessage
)

# TRAINING side: a container whose stored codes reach its own declared K is
# refused - that K becomes the store's category count
data.over.k <- dbartsData(train.bound, y.bound)
data.over.k@x$sparseCategoryCount <- 2L
expect_error(
  dbarts(data.over.k, control = boundControl),
  pattern = "categorical predictors must hold integer category codes"
)

# LOCK-IN: the entrances that already bounded against the store keep doing so
sampler.saved <- dbarts(
  train.bound,
  y.bound,
  control = dbartsControl(
    n.trees = 10L,
    n.chains = 1L,
    n.threads = 1L,
    n.samples = 5L,
    n.burn = 0L,
    updateState = FALSE,
    keepTrees = TRUE
  )
)
invisible(sampler.saved$run(20L, 5L))
expect_error(sampler.saved$predict(over.codes), pattern = levelMessage)
expect_error(
  sampler.saved$getTrees(newdata = test.dense.over),
  pattern = "categorical predictor values must be existing category codes"
)

# the code bound is untouched where no level table can speak for the codes: a
# container whose levels match training exactly leaves the alignment nothing to
# do, so a planted code - stored or reference - goes straight to the bridge
aligned.container <- dbartsData(
  sparseFrame(c("b", "c", "b"), levels.small, "a"),
  rnorm(3L)
)@x
planted.codes <- aligned.container
planted.codes$sparse@x[1L] <- 7
expect_error(
  sampler.bound$setTestPredictor(planted.codes),
  pattern = codeMessage
)
planted.reference <- aligned.container
planted.reference$sparseReference[1L] <- 7L
planted.reference$sparseCategoryCount[1L] <- 8L
expect_error(
  sampler.bound$setTestPredictor(planted.reference),
  pattern = codeMessage
)

# NO FALSE REFUSALS: a container is judged on its codes, not on the K it
# declares, so both a larger-K and a smaller-K foreign table are accepted while
# every code stays inside the training count; NA is the reserved missing code
in.range.big <- dbartsData(
  sparseFrame(c("b", "c", "b", "c"), levels.big, "a"),
  rnorm(4L)
)@x
sampler.bound$setTestPredictor(in.range.big)
expect_true(all(is.finite(sampler.bound$run(10L, 10L)$test)))
in.range.small <- dbartsData(
  sparseFrame(c("a", "b", "b"), c("a", "b"), "a"),
  rnorm(3L)
)@x
sampler.bound$setTestPredictorAndOffset(in.range.small, rep(0, 3L))
expect_true(all(is.finite(sampler.bound$run(10L, 10L)$test)))

# MUTATION LIFT: a design carrying a sparse CATEGORICAL column takes
# column-granular mutation. The engine keys the replacement's nonzero
# pattern on the column's kind - {i : code != refCode} for a categorical
# column, not the ordinal {i : value != 0} - and the R side mirrors it into
# data@x, leaving the reference and the declared level count
# creation-pinned.
codes.bound <- as.matrix(sampler.bound$data@x)[, 2L]
new.codes <- (codes.bound + 1) %% 3
expect_silent(
  sampler.bound$setPredictor(new.codes, column = 2L, forceUpdate = TRUE)
)
expect_equal(as.matrix(sampler.bound$data@x)[, 2L], new.codes)
expect_equal(sampler.bound$data@x$sparseReference, 0L)
expect_equal(sampler.bound$data@x$sparseCategoryCount, 3L)
# the mirrored block stores exactly the non-reference cells
expect_equal(diff(sampler.bound$data@x$sparse@p), sum(new.codes != 0))
# a dense-backed column of the same design mutates alongside it
expect_silent(
  sampler.bound$setPredictor(
    train.bound$x1 * 2,
    column = 1L,
    forceUpdate = TRUE
  )
)
expect_equal(as.matrix(sampler.bound$data@x)[, 1L], train.bound$x1 * 2)

# the flat entries the bridge exposes accept it too (an identity re-install, so
# the data@x these bypass stays current), and the sampler still fits
predictors.bound <- as.matrix(sampler.bound$data@x)
storage.mode(predictors.bound) <- "double"
expect_true(.Call(
  dbarts:::C_dbarts_bartcore_setPredictor,
  sampler.bound$getPointer(),
  predictors.bound,
  TRUE,
  FALSE
))
expect_true(.Call(
  dbarts:::C_dbarts_bartcore_updatePredictor,
  sampler.bound$getPointer(),
  predictors.bound[, 2L],
  2L,
  TRUE,
  FALSE
))
expect_true(all(is.finite(sampler.bound$run(10L, 10L)$train)))

# the mutated sparse design draws exactly as the dense-equivalent design does
dense.bound <- data.frame(
  x1 = train.bound$x1,
  f = factor(labels.bound, levels = levels.small)
)
mutationGate <- function(frame, seed) {
  set.seed(seed)
  sampler <- dbarts(frame, y.bound, sigma = 1.0, control = boundControl)
  invisible(sampler$run(5L, 5L))
  invisible(sampler$setPredictor(new.codes, column = 2L, forceUpdate = TRUE))
  sampler$run(0L, 10L)
}
result.mutated.sparse <- mutationGate(train.bound, 4242L)
result.mutated.dense <- mutationGate(dense.bound, 4242L)
expect_identical(result.mutated.sparse$train, result.mutated.dense$train)
expect_identical(result.mutated.sparse$sigma, result.mutated.dense$sigma)

# EXPLICIT ZEROS: a reference level that is not levels[1] is the only shape
# whose sparse block stores code 0 explicitly, so it is the shape a matrix-wide
# drop0 would silently corrupt. Mutating a sparse ORDINAL column beside it must
# leave the categorical column's cells bit-identical.
set.seed(84L)
n.zero <- 120L
labels.zero <- sample(levels.small, n.zero, replace = TRUE)
frame.zero <- data.frame(x1 = rnorm(n.zero))
frame.zero$g <- sparseFactor(
  labels.zero,
  levels = levels.small,
  reference = "b"
)
frame.zero$s <- Matrix::sparseVector(
  x = 0.5 + runif(10L),
  i = sort(sample.int(n.zero, 10L)),
  length = n.zero
)
y.zero <- rnorm(n.zero) + match(labels.zero, levels.small)
sampler.zero <- dbarts(frame.zero, y.zero, control = boundControl)
invisible(sampler.zero$run(5L, 5L))
categorical.rank <- -sampler.zero$data@x$map[2L]
entriesFor <- function(container, rank) {
  pointers <- container$sparse@p
  seq.int(
    pointers[rank] + 1L,
    length.out = pointers[rank + 1L] - pointers[rank]
  )
}
zeros.before <- sum(
  sampler.zero$data@x$sparse@x[entriesFor(
    sampler.zero$data@x,
    categorical.rank
  )] ==
    0
)
expect_true(zeros.before > 0L)
categorical.before <- as.matrix(sampler.zero$data@x)[, 2L]

new.ordinal <- numeric(n.zero)
new.ordinal[c(3L, 40L, 99L)] <- c(1.5, 2.5, 3.5)
expect_silent(
  sampler.zero$setPredictor(new.ordinal, column = 3L, forceUpdate = TRUE)
)
expect_equal(as.matrix(sampler.zero$data@x)[, 3L], new.ordinal)
expect_equal(as.matrix(sampler.zero$data@x)[, 2L], categorical.before)
expect_equal(
  sum(
    sampler.zero$data@x$sparse@x[
      entriesFor(sampler.zero$data@x, categorical.rank)
    ] ==
      0
  ),
  zeros.before
)

# and mutating the categorical column keeps the ordinal one bit-identical
ordinal.before <- as.matrix(sampler.zero$data@x)[, 3L]
new.zero.codes <- (as.matrix(sampler.zero$data@x)[, 2L] + 2) %% 3
expect_silent(
  sampler.zero$setPredictor(new.zero.codes, column = 2L, forceUpdate = TRUE)
)
expect_equal(as.matrix(sampler.zero$data@x)[, 2L], new.zero.codes)
expect_equal(as.matrix(sampler.zero$data@x)[, 3L], ordinal.before)
# the reference level is code 1 here, so the stored entries are the other two
expect_equal(
  diff(sampler.zero$data@x$sparse@p)[categorical.rank],
  sum(new.zero.codes != 1)
)

# a rejected transactional update leaves data@x exactly as it was
before.transaction <- as.matrix(sampler.zero$data@x)
candidate <- rep(0.25, n.zero)
accepted <- sampler.zero$setPredictor(candidate, column = 1L)
after.transaction <- as.matrix(sampler.zero$data@x)
expect_equal(
  after.transaction[, 1L],
  if (isTRUE(accepted)) candidate else before.transaction[, 1L]
)
expect_equal(after.transaction[, -1L], before.transaction[, -1L])

# SAVE/RE-CREATE: the C sampler is rebuilt from data@x after a mutation, so a
# restored continuation has to reproduce the mutated one
stateControl <- dbartsControl(
  n.trees = 10L,
  n.chains = 1L,
  n.threads = 1L,
  n.samples = 5L,
  n.burn = 0L,
  updateState = TRUE
)
sampler.state <- dbarts(frame.zero, y.zero, control = stateControl)
invisible(sampler.state$run(5L, 5L))
invisible(sampler.state$setPredictor(
  new.ordinal,
  column = 3L,
  forceUpdate = TRUE
))
invisible(
  sampler.state$setPredictor(new.zero.codes, column = 2L, forceUpdate = TRUE)
)
sampler.state$storeState()
stateFile <- tempfile(fileext = ".rds")
saveRDS(sampler.state, stateFile)
sampler.restored <- readRDS(stateFile)
invisible(file.remove(stateFile))
set.seed(313L)
run.mutated <- sampler.state$run(0L, 5L)
set.seed(313L)
run.recreated <- sampler.restored$run(0L, 5L)
expect_equal(run.mutated$train, run.recreated$train)
expect_equal(run.mutated$sigma, run.recreated$sigma)

# ALIGNMENT: a foreign container - one assembled from another data set, so
# carrying its own level order - is re-coded against the training level table
# at the validateXTest funnel, and then predicts exactly as the equivalent
# data frame does
set.seed(85L)
n.align <- 20L
labels.align <- sample(levels.small, n.align, replace = TRUE)
frame.align <- data.frame(x1 = rnorm(n.align))
frame.align$f <- sparseFactor(
  labels.align,
  levels = rev(levels.small),
  reference = "c"
)
container.align <- dbartsData(frame.align, rnorm(n.align))@x
expect_equal(attr(container.align, "factor.levels")[[2L]], rev(levels.small))

sampler.align <- dbarts(train.bound, y.bound, control = boundControl)
sampler.align$setTestPredictor(container.align)
aligned <- sampler.align$data@x.test
expect_equal(
  as.matrix(aligned)[, 2L],
  as.double(match(labels.align, levels.small) - 1L)
)
expect_equal(aligned$sparseReference, 2L)
expect_equal(aligned$sparseCategoryCount, 3L)

frame.align.dense <- data.frame(
  x1 = frame.align$x1,
  f = factor(labels.align, levels = rev(levels.small))
)
alignmentGate <- function(test, seed) {
  set.seed(seed)
  dbarts(
    train.bound,
    y.bound,
    sigma = 1.0,
    test = test,
    control = boundControl
  )$run(5L, 10L)
}
expect_identical(
  alignmentGate(container.align, 606L)$test,
  alignmentGate(frame.align.dense, 606L)$test
)

# a level the training design never saw still refuses, by name
frame.align.new <- data.frame(x1 = rnorm(3L))
frame.align.new$f <- sparseFactor(
  c("a", "z", "b"),
  levels = c(levels.small, "z"),
  reference = "a"
)
expect_error(
  sampler.align$setTestPredictor(dbartsData(frame.align.new, rnorm(3L))@x),
  pattern = levelMessage
)

# the reference-level pin: a container's declared sparse reference level meant one thing to
# predict()'s R-side densification (as.matrix.dbartsMixedMatrix, gated only on
# is.na) and another to the engine (a reference is read only for a
# store-categorical column, 0 otherwise) - reachable by training with an
# ORDERED factor column (typed ordinal, but still carrying a factor.levels
# entry) and supplying that column as a sparseFactor test column;
# mapFactorColumnsToTrainingLevels remaps it rather than refusing, since its
# numeric-training guard keys on a NULL level table, not on varTypes. Both
# funnels now refuse it identically, at the validateXTest choke point they
# both share.
levels.ord <- c("lo", "mid", "hi")
f.ord <- factor(
  sample(levels.ord, 120L, replace = TRUE),
  levels = levels.ord,
  ordered = TRUE
)
y.ord <- c(-10, 0, 10)[as.integer(f.ord)] + rnorm(120L, 0, 0.5)
train.ord <- data.frame(f = f.ord)
expect_equal(
  attr(dbarts:::makeCategoricalModelMatrix(train.ord), "varTypes"),
  0L
)

sampler.ord <- dbarts(
  train.ord,
  y.ord,
  control = dbartsControl(
    n.trees = 1L,
    n.chains = 1L,
    n.threads = 1L,
    updateState = FALSE
  )
)
invisible(sampler.ord$run(50L, 1L))

test.ord.bad <- data.frame(
  f = factor(rep("hi", 5L), levels = levels.ord, ordered = TRUE)
)
test.ord.bad$f <- sparseFactor(test.ord.bad$f, reference = "hi")

referenceMessage <-
  "a sparse predictor column may declare a reference level only for a categorical predictor"
expect_error(sampler.ord$predict(test.ord.bad), pattern = referenceMessage)
expect_error(
  sampler.ord$setTestPredictor(test.ord.bad),
  pattern = referenceMessage
)

# the bridge's own check is independent of validateXTest: a container fed
# directly to the bridge entry point (bypassing the R wrapper, the way the
# equivalence harness does) is refused there too
bc.ord <- dbarts:::bartcoreSampler(sampler.ord)
bad.container.ord <- dbarts:::makeCategoricalModelMatrix(test.ord.bad)
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_setTestPredictor,
    bc.ord$ptr,
    bad.container.ord
  ),
  pattern = referenceMessage
)

# NO FALSE REFUSALS, continued: an NA test code must not be read as an
# out-of-range CATEGORY. That only holds on a column whose training values
# carried an NA too (D8's new-route rule), so this needs its own NA-bearing
# training design rather than mutating train.bound/sampler.bound, which every
# code-bound assertion above shares. sparseFactor() itself refuses NA
# (missing values are not supported in a sparseFactor), so the categorical
# column here is dense, the train.bound.dense idiom above. Placed at the
# file's end since a new rnorm() call here would shift every seeded draw that
# follows it.
labels.na <- labels.bound
labels.na[2L] <- NA
train.na <- data.frame(
  x1 = train.bound$x1,
  f = factor(labels.na, levels = levels.small)
)
test.dense.na <- cbind(rnorm(3L), c(0, NA, 2))
colnames(test.dense.na) <- c("x1", "f")
expect_inherits(
  dbarts(
    dbartsData(train.na, y.bound, test = test.dense.na),
    control = boundControl
  ),
  "dbartsSampler"
)
