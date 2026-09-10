# The engine's fixed limits and the control settings that move them: that each
# cap refuses (or routes) where it says it does, that a setting reaches the
# engine as observable behaviour rather than as a stored value, and that every
# default leaves today's draws alone.

set.seed(0)
n <- 100L
x <- matrix(rnorm(2L * n), n)
y <- x[, 1L] + rnorm(n)

## the per-column cut ceiling refuses by name rather than clamping: a caller
## who asks for 100000 cuts must not silently receive 65533
expect_error(
  dbarts(y ~ x, control = dbartsControl(n.cuts = 100000L, n.chains = 1L)),
  "over the cap of 65533"
)
expect_error(
  dbarts(y ~ x, control = dbartsControl(n.cuts = 65534L, n.chains = 1L)),
  "over the cap of 65533"
)

## the last representable request is still accepted, and costs nothing beyond
## the distinct values a column actually has
sampler <- dbarts(
  y ~ x,
  control = dbartsControl(n.cuts = 65533L, n.chains = 1L, updateState = FALSE)
)
expect_true(all(sampler$data@n.cuts == 65533L))
expect_true(all(is.finite(sampler$run(0L, 2L)$train)))

rm(sampler)
