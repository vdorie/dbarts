# The `proposal.probs` surface at two structural names. The swap proposal is
# gone from the kernel, so a vector still naming it is refused BY NAME at every
# entry point rather than dropped, zeroed or folded into birth_death; the check
# runs ahead of both the monotone rewrite (which would discard a swap-carrying
# vector matching the new default) and the two-forest refusal (which would
# report only a sum). What survives is the one-NA fill over the two names and
# the sum-to-one rule.

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

removed <- "swap"

# ---- refused by name, with and without a monotone constraint --------------

# the old documented default
oldDefault <- c(birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5)
expect_error(fit(oldDefault), removed)
expect_error(fit(oldDefault, monotone = c(a = "+")), removed)

# a vector semantically identical to the new default, but naming the move
zeroed <- c(birth_death = 0.6, swap = 0, change = 0.4, birth = 0.5)
expect_error(fit(zeroed), removed)
expect_error(fit(zeroed, monotone = c(a = "+")), removed)

# the name alone
expect_error(fit(c(swap = 0.1)), removed)
expect_error(fit(c(swap = 0.1), monotone = c(a = "+")), removed)

# the same message reaches dbartsSpec() and a directly built model
expect_error(
  dbarts::dbartsSpec(
    dbarts::dbartsData(x, y),
    control = control,
    proposal.probs = oldDefault
  ),
  removed
)
expect_error(
  methods::new("dbartsModel", proposal.probs = oldDefault),
  removed
)

# ---- accepted --------------------------------------------------------------

# the new default proceeds, and under a constraint is rewritten silently
newDefault <- c(birth_death = 0.6, change = 0.4, birth = 0.5)
expect_equal(fit(newDefault)$model@p.birth_death, 0.6)
expect_equal(fit(newDefault)$model@p.change, 0.4)
forced <- fit(newDefault, monotone = c(a = "+"))$model
expect_equal(forced@p.birth_death, 1)
expect_equal(forced@p.change, 0)

# one name fills the other in from the residual
partial <- fit(c(birth_death = 0.7))$model
expect_equal(partial@p.birth_death, 0.7)
expect_equal(partial@p.change, 0.3)
expect_equal(fit(c(change = 0.25))$model@p.birth_death, 0.75)

# a partial vector is not the default vector, so it still trips the monotone
# stop rather than being rewritten
expect_error(
  fit(c(birth_death = 0.7), monotone = c(a = "+")),
  "proposal.probs"
)

# the defaulted call and the explicit new default agree
expect_equal(
  dbarts::dbarts(x, y, control = control)$model@p.birth_death,
  0.6
)
