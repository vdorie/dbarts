# tinytest's expect_warning(expr, pattern) passes as soon as any warning
# expr emits matches pattern; it does not check how many warnings fired, so
# an extra unrelated warning or a duplicate of the expected one still
# passes. captureWarnings() runs expr once, muffles and records every
# warning condition it raises (letting expr run to completion), and returns
# them in emission order so a call site can assert an exact count alongside
# the pattern/class match.
captureWarnings <- function(expr) {
  observed <- list()
  withCallingHandlers(
    expr,
    warning = function(w) {
      observed[[length(observed) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  observed
}
