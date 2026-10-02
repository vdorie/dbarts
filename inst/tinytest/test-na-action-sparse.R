# na.action on the sparse containers: the rows a stored NA marks are found off
# the stored entries, and are the rows base R's na.omit drops from the same
# predictors densified. A binary response keeps the starting sigma (and its
# sparse-design fallback warning) out of it

if (!requireNamespace("Matrix", quietly = TRUE)) {
  exit_file("Matrix not available")
}

omittedRows <- function(y, x) {
  as.integer(attr(stats::na.omit(data.frame(y = y, x)), "na.action"))
}

set.seed(7801)
n <- 10L
y <- rep(c(0, 1), length.out = n)

# a bare dgCMatrix: an NA entry at row 3 (stored 0-based as row index 2)
x.sparse <- Matrix::sparseMatrix(
  i = c(2L, 3L, 5L, 7L),
  j = c(1L, 1L, 2L, 2L),
  x = c(1, NA, 2, 3),
  dims = c(n, 2L)
)
data.sparse <- dbartsData(x.sparse, y, na.action = na.omit)
expect_equal(
  as.integer(data.sparse@na.action),
  omittedRows(y, as.matrix(x.sparse))
)
expect_equal(data.sparse@y, y[-3L])

# a mixed container: a dense column's NA at row 4 beside a sparse block with
# none still drops row 4
x1 <- runif(n)
x1[4L] <- NA
sv <- Matrix::sparseVector(x = c(1, 2), i = c(2L, 9L), length = n)
x.mixed <- data.frame(x1 = x1)
x.mixed$sv <- sv
data.mixed <- dbartsData(x.mixed, y, na.action = na.omit)
expect_inherits(data.mixed@x, "dbartsMixedMatrix")
expect_equal(
  as.integer(data.mixed@na.action),
  omittedRows(y, data.frame(x1 = x1, sv = as.vector(sv)))
)
expect_equal(data.mixed@y, y[-4L])

rm(data.sparse, data.mixed, x.sparse, x.mixed, x1, sv, y, n, omittedRows)
