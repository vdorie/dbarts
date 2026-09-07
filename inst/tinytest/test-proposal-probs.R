# The `proposal.probs` surface at four structural names. Swap and perturb ship
# at zero but stay legal names: a caller who asks for either gets it, and the
# one-NA fill and the sum-to-one rule run over all four. Perturb resolves ahead
# of the fill, so every resolution the three-name rule pinned is preserved
# element for element.

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
expect_equal(defaulted@p.perturb, 0)

# spelling the default explicitly agrees with defaulting it
explicit <- fit(
  c(birth_death = 0.6, swap = 0, change = 0.4, perturb = 0, birth = 0.5)
)$model
expect_equal(explicit@p.birth_death, 0.6)
expect_equal(explicit@p.swap, 0)
expect_equal(explicit@p.change, 0.4)
expect_equal(explicit@p.perturb, 0)

# and so does the four-name spelling that omits perturb, which is the one
# every consumer forwarding the documented default writes
omitted <- fit(c(birth_death = 0.6, swap = 0, change = 0.4, birth = 0.5))$model
expect_equal(omitted@p.perturb, 0)

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
  "birth/death 0.50, swap 0.10, change 0.40, perturb 0.00; birth 0.50",
  printed,
  fixed = TRUE
)))

# ---- the fill and the sum --------------------------------------------------

# one unnamed element takes the residual, whichever it is
partial <- fit(c(birth_death = 0.7, swap = 0))$model
expect_equal(partial@p.birth_death, 0.7)
expect_equal(partial@p.swap, 0)
expect_equal(partial@p.change, 0.3)
expect_equal(fit(c(birth_death = 0.5, change = 0.4))$model@p.swap, 0.1)

# two unnamed, one of them swap: swap takes its zero and the other the residual
bdOnly <- fit(c(birth_death = 0.7))$model
expect_equal(bdOnly@p.birth_death, 0.7)
expect_equal(bdOnly@p.swap, 0)
expect_equal(bdOnly@p.change, 0.3)
changeOnly <- fit(c(change = 0.25))$model
expect_equal(changeOnly@p.birth_death, 0.75)
expect_equal(changeOnly@p.swap, 0)
expect_equal(changeOnly@p.change, 0.25)

# a zero-default move named alone leaves the birth/death-versus-change split
# undetermined, whichever one it is
expect_error(fit(c(swap = 0.1)), "name at least one of")
expect_error(fit(c(perturb = 0.16)), "name at least one of")
expect_error(fit(c(swap = 0.1, perturb = 0.16)), "name at least one of")

# all three unnamed falls back to the default
expect_equal(fit(c(birth = 0.25))$model@p.birth_death, 0.6)
expect_equal(fit(c(birth = 0.25))$model@p.swap, 0)
expect_equal(fit(c(birth = 0.25))$model@p.change, 0.4)
expect_equal(fit(c(birth = 0.25))$model@p.perturb, 0)
expect_equal(fit(c(birth = 0.25))$model@p.birth, 0.25)

# perturb resolves BEFORE the fill: it takes its zero rather than the
# residual, so the residual is taken against 1 - perturb
withPerturb <- fit(c(birth_death = 0.5, change = 0.34, perturb = 0.16))$model
expect_equal(withPerturb@p.swap, 0)
expect_equal(withPerturb@p.perturb, 0.16)
twoUnnamed <- fit(c(birth_death = 0.84, perturb = 0.16))$model
expect_equal(twoUnnamed@p.swap, 0)
expect_equal(twoUnnamed@p.change, 0)
expect_equal(twoUnnamed@p.perturb, 0.16)

# all four named must sum to one
expect_error(
  fit(c(birth_death = 0.7, swap = 0.1, change = 0.4)),
  "sum to 1"
)
expect_error(
  fit(c(birth_death = 0.6, swap = 0, change = 0.4, perturb = 0.1)),
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
expect_equal(forced@p.perturb, 0)
expect_error(fit(threeMove, monotone = c(a = "+")), "proposal.probs")

# the refusal must not fire on a caller who spells the documented default and
# omits the move that ships at zero: the comparison fills it first
forcedFive <- fit(
  c(birth_death = 0.6, swap = 0, change = 0.4, perturb = 0, birth = 0.5),
  monotone = c(a = "+")
)$model
expect_equal(forcedFive@p.birth_death, 1)
expect_equal(forcedFive@p.perturb, 0)

# ---- a perturb-dominant run --------------------------------------------

# One tree, most of its proposals cut displacements. A birth or a death moves
# the node count, so consecutive draws holding it fixed took a displacement or
# nothing: their split variables and shape must agree exactly, and at most one
# cut may have moved.
perturbControl <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 1L,
  n.burn = 0L,
  n.samples = 60L,
  keepTrees = TRUE,
  updateState = FALSE
)
perturbSampler <- dbarts::dbarts(
  x,
  y,
  control = perturbControl,
  proposal.probs = c(
    birth_death = 0.2,
    swap = 0,
    change = 0,
    perturb = 0.8,
    birth = 0.5
  )
)
set.seed(41L)
perturbSamples <- perturbSampler$run()
expect_true(all(is.finite(perturbSamples$train)))
perturbTrees <- perturbSampler$getTrees()
byDraw <- split(perturbTrees[, c("var", "value")], perturbTrees$sample)
heldShape <- 0L
cutMoves <- 0L
for (i in seq_len(length(byDraw) - 1L)) {
  before <- byDraw[[i]]
  after <- byDraw[[i + 1L]]
  if (nrow(before) != nrow(after)) {
    next
  }
  heldShape <- heldShape + 1L
  expect_identical(before$var, after$var)
  interior <- before$var != -1L
  moved <- sum(before$value[interior] != after$value[interior])
  expect_true(moved <= 1L)
  cutMoves <- cutMoves + moved
}
# the run exercised both arms: the shape moved somewhere, and cuts moved
# without it
expect_true(length(unique(vapply(byDraw, nrow, 0L))) > 1L)
expect_true(heldShape > 10L)
expect_true(cutMoves > 0L)
