# Stripping every column to zero cut points leaves a root-only sampler with no
# available split variable anywhere. Birth/death must treat each tree's move as
# a no-op rather than force a birth and draw a rule for an invalid variable.
# Regression for the degenerate-root guard (src/bartcore/moves.hpp): unfixed,
# sampler$run segfaults here.

set.seed(0)
n <- 60L
p <- 2L
x <- matrix(runif(n * p), n, p)
y <- rnorm(n)

control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 5L,
  updateState = FALSE
)
sampler <- dbarts(x, y, control = control)

sampler$setCutPoints(list(numeric(0), numeric(0)), 1:2)

res <- sampler$run(0L, 5L)

# the run completes with finite output instead of crashing
expect_true(all(is.finite(res$train)))
expect_true(all(is.finite(res$sigma)))

# every tree stays a lone root: one node per tree, all leaves (var == -1)
trees <- sampler$getTrees(current = TRUE)
expect_equal(nrow(trees), 5L)
expect_true(all(trees$var == -1L))

rm(sampler, res, trees, control, x, y, n, p)
