# Shape predicate for an xbart() result over the (n.trees, k, power, base)
# grid: array class, dim, no NA, and dimnames threaded from the grid values
# by name. Every source() of this file must pass local = TRUE for these
# assertions to reach the run's masked expect_* bindings.
checkXvalShape <- function(xval, n.reps, n.trees, k, power, base) {
  # nolint next: object_usage_linter. tinytest attaches expect_* at run time.
  expect_inherits(xval, "array")
  # nolint next: object_usage_linter. tinytest attaches expect_* at run time.
  expect_equal(
    dim(xval),
    c(n.reps, length(n.trees), length(k), length(power), length(base))
  )
  # nolint next: object_usage_linter. tinytest attaches expect_* at run time.
  expect_true(!anyNA(xval))
  # nolint next: object_usage_linter. tinytest attaches expect_* at run time.
  expect_equal(
    dimnames(xval),
    list(
      rep = NULL,
      n.trees = as.character(n.trees),
      k = as.character(k),
      power = as.character(power),
      base = as.character(base)
    )
  )
}
