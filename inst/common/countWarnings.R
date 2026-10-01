# Counts the warnings expr raises that inherit from class or, when pattern is
# given, whose message matches it, muffling each so expr runs to completion.
countWarnings <- function(expr, class = NULL, pattern = NULL) {
  count <- 0L
  withCallingHandlers(
    expr,
    warning = function(w) {
      hit <- if (is.null(pattern)) {
        inherits(w, class)
      } else {
        grepl(pattern, conditionMessage(w))
      }
      if (hit) {
        count <<- count + 1L
      }
      invokeRestart("muffleWarning")
    }
  )
  count
}
