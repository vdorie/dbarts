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
  pattern = "sparse categorical predictors are not yet supported"
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
