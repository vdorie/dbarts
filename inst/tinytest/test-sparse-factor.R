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
  pattern = "sparse categorical predictors must be supplied through the x/y interface"
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
# same values given as a dense factor. Both storage kinds see K categories,
# so the top level is planted; the reference is a middle level (never first
# alphabetically), whose own level-order code the implicit rows carry.
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
  codes[1L] <- length(levels) # a dense factor's count is max(code) + 1
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
test.mix.sparse$s <- Matrix::sparseVector(x = s.vals, i = s.rows, length = n.mtest)

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

# a per-column update against a resident sparse test store is refused: the
# container takes whole-object replacement only, and the refusal is inert
sampler.container.mut <- makeMixedSampler(test.mix.sparse, 101L)
expect_inherits(sampler.container.mut$data@x.test, "dbartsMixedMatrix")
expect_error(
  sampler.container.mut$setTestPredictor(rnorm(n.mtest), column = 1L),
  pattern = "single column"
)
expect_inherits(sampler.container.mut$data@x.test, "dbartsMixedMatrix")

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
  pattern = "length of test offset"
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
