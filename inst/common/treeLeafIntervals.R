# A getTrees() data frame walked in prefix order reconstructs each leaf's
# (lo, hi] split interval; leafOfObservations() then maps values to the
# index of the leaf interval that contains them.
parseLeaves <- function(rows) {
  i <- 0L
  leaves <- list()
  recurse <- function(lo, hi) {
    i <<- i + 1L
    row <- rows[i, ]
    if (row$var == -1) {
      leaves[[length(leaves) + 1L]] <<- list(lo = lo, hi = hi, n = row$n)
      return(invisible())
    }
    v <- row$value
    recurse(lo, v) # left child:  x <= value
    recurse(v, hi) # right child: x >  value
  }
  recurse(-Inf, Inf)
  leaves
}
leafOfObservations <- function(leaves, values) {
  vapply(
    values,
    function(xj) {
      which(vapply(leaves, function(l) l$lo < xj & xj <= l$hi, logical(1)))
    },
    integer(1)
  )
}
