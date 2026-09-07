# The `proposal.probs` surface at three structural names. Swap ships at zero
# but stays a legal name: a caller who asks for it gets it, and the one-NA fill
# and the sum-to-one rule run over all three.

set.seed(31L)
n <- 60L
p <- 3L
x <- matrix(runif(n * p), n, p, dimnames = list(NULL, c("a", "b", "c")))
y <- x[, 1L] - x[, 3L] + rnorm(n, 0, 0.2)
control <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 5L,
  updateState = FALSE
)
fit <- function(probs, ...) {
  dbarts::dbarts(x, y, control = control, proposal.probs = probs, ...)
}

# ---- the shipped default ---------------------------------------------------

defaulted <- dbarts::dbarts(x, y, control = control)$model
expect_equal(defaulted@p.birth_death, 0.6)
expect_equal(defaulted@p.swap, 0)
expect_equal(defaulted@p.change, 0.4)

# spelling the default explicitly agrees with defaulting it
explicit <- fit(c(birth_death = 0.6, swap = 0, change = 0.4, birth = 0.5))$model
expect_equal(explicit@p.birth_death, 0.6)
expect_equal(explicit@p.swap, 0)
expect_equal(explicit@p.change, 0.4)

# ---- a caller-supplied three-move mixture ----------------------------------

# swap is the only move that rotates a child's rule up the tree, so a
# single-tree fit is the case that needs it positive
threeMove <- c(birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5)
oneTree <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 1L,
  n.burn = 0L,
  n.samples = 20L,
  updateState = FALSE
)
sampler <- dbarts::dbarts(x, y, control = oneTree, proposal.probs = threeMove)
expect_equal(sampler$model@p.birth_death, 0.5)
expect_equal(sampler$model@p.swap, 0.1)
expect_equal(sampler$model@p.change, 0.4)

set.seed(9L)
samples <- sampler$run()
expect_true(all(is.finite(samples$train)))
expect_true(all(is.finite(samples$sigma)))
# the mixture the sampler was created with survives the run
expect_equal(sampler$model@p.swap, 0.1)

# the creation printout names all three probabilities
printed <- capture.output(
  dbarts::dbarts(
    x,
    y,
    control = dbarts::dbartsControl(
      n.chains = 1L,
      n.threads = 1L,
      n.trees = 1L,
      updateState = FALSE,
      verbose = TRUE
    ),
    proposal.probs = threeMove
  )
)
expect_true(any(grepl(
  "birth/death 0.50, swap 0.10, change 0.40",
  printed,
  fixed = TRUE
)))

# ---- the fill and the sum --------------------------------------------------

# one missing name takes the residual
partial <- fit(c(birth_death = 0.7, swap = 0))$model
expect_equal(partial@p.birth_death, 0.7)
expect_equal(partial@p.swap, 0)
expect_equal(partial@p.change, 0.3)
expect_equal(fit(c(birth_death = 0.5, change = 0.4))$model@p.swap, 0.1)

# all three missing falls back to the default
expect_equal(fit(c(birth = 0.25))$model@p.birth_death, 0.6)
expect_equal(fit(c(birth = 0.25))$model@p.birth, 0.25)

# more than one missing name is not a fill, and the remainder must sum to one
expect_error(
  fit(c(birth_death = 0.7, swap = 0.1, change = 0.4)),
  "sum to 1"
)

# ---- the monotone rewrite --------------------------------------------------

# a defaulted mixture is rewritten birth/death-only; a non-default one errors
forced <- fit(
  c(birth_death = 0.6, swap = 0, change = 0.4, birth = 0.5),
  monotone = c(a = "+")
)$model
expect_equal(forced@p.birth_death, 1)
expect_equal(forced@p.swap, 0)
expect_equal(forced@p.change, 0)
expect_error(fit(threeMove, monotone = c(a = "+")), "proposal.probs")
