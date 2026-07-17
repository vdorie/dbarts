# C_dbarts_assignInPlace validates its R-supplied index and source before the
# per-iteration in-place write (the rbart/partialDependence hot path). The
# valid callers are exercised throughout the rbart and partial-dependence
# tests; here we pin the guards it adds on bad input.
assignInPlace <- dbarts:::C_dbarts_assignInPlace

target <- matrix(0, 3L, 2L)

# a column index past the trailing dimension is rejected, not written OOB
expect_error(
  .Call(assignInPlace, target, 5L, c(1, 2, 3)),
  "out of bounds"
)

# a source whose length does not fill the leading dimension is rejected
expect_error(
  .Call(assignInPlace, target, 1L, c(1, 2)),
  "leading dimension"
)

# a non-integer index is rejected before it is read as one
expect_error(
  .Call(assignInPlace, target, 1.5, c(1, 2, 3)),
  "must be integer"
)

# a scalar write past the end of an undimensioned target is rejected
flat <- numeric(4L)
expect_error(
  .Call(assignInPlace, flat, 9L, 1),
  "out of bounds"
)

rm(assignInPlace, target, flat)
