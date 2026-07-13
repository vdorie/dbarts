# sparseFactor: the sparse unordered-factor wrapper (decision S-CAT:
# recognized by ingestion, refused at data construction until the
# CSC-categorical engine path lands)

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

# refusal at data construction, x/y interface: the column is recognized as
# sparse and refused before anything reaches the bridge
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
expect_error(
  dbartsData(df, y),
  pattern = "sparse categorical predictors are not yet supported"
)
expect_error(
  dbarts(df, y),
  pattern = "sparse categorical predictors are not yet supported"
)
expect_error(
  dbartsData(df, y, factors = "indicators"),
  pattern = "sparse categorical predictors are not yet supported"
)

# refusal on the formula path, ahead of model.frame's bare S4 type error
df.formula <- df
df.formula$y <- y
expect_error(
  dbartsData(y ~ x1 + sf, df.formula),
  pattern = "sparse categorical predictors are not yet supported"
)
expect_error(
  dbarts(y ~ x1 + sf, df.formula),
  pattern = "sparse categorical predictors are not yet supported"
)

# a bare sparseFactor as 'x' refuses too
expect_error(
  dbartsData(df$sf, y),
  pattern = "sparse categorical predictors are not yet supported"
)

# sparse ordinal columns are untouched by the refusal (S-CAT boundary)
if (requireNamespace("Matrix", quietly = TRUE)) {
  df.ord <- data.frame(x1 = rnorm(n))
  df.ord$sv <- Matrix::sparseVector(
    x = 0.5 + runif(5L),
    i = sort(sample.int(n, 5L)),
    length = n
  )
  expect_inherits(dbartsData(df.ord, y), "dbartsData")
}
