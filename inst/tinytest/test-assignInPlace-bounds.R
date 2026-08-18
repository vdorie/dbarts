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

# a non-scalar source is rejected once the index names every dimension
expect_error(
  .Call(assignInPlace, target, c(1L, 1L), c(1, 2)),
  "source must be a scalar"
)

# an index naming neither all dimensions nor all but the first is rejected
expect_error(
  .Call(assignInPlace, target, c(1L, 1L, 1L), c(1, 2, 3)),
  "all but the first array dimension"
)

# an integer source against a double target is rejected, not coerced
expect_error(
  .Call(assignInPlace, flat, 1L, 1L),
  "must be double"
)

rm(assignInPlace, target, flat)
