# The `proposal.probs` surface at five structural names. Swap, perturb and
# rule_gibbs ship at zero but stay legal names: a caller who asks for any of
# them gets it, and the one-NA fill and the sum-to-one rule run over all five.
# The two zero-default moves resolve ahead of the fill, so every resolution the
# three-name rule pinned is preserved element for element.

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
# the mixture is a control setting; a fit's resolved copy is read back off the
# sampler's own control
mixtureControl <- function(probs) {
  dbarts::dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 5L,
    updateState = FALSE,
    proposal.probs = probs
  )
}
fit <- function(probs, ...) {
  dbarts::dbarts(x, y, control = mixtureControl(probs), ...)
}

# ---- the shipped default ---------------------------------------------------

defaulted <- dbarts::dbarts(x, y, control = control)$control
expect_equal(defaulted@proposal.probs[["birth_death"]], 0.6)
expect_equal(defaulted@proposal.probs[["swap"]], 0)
expect_equal(defaulted@proposal.probs[["change"]], 0.4)
expect_equal(defaulted@proposal.probs[["perturb"]], 0)
expect_equal(defaulted@proposal.probs[["rule_gibbs"]], 0)

# spelling the default explicitly agrees with defaulting it
explicit <- fit(
  c(
    birth_death = 0.6,
    swap = 0,
    change = 0.4,
    perturb = 0,
    rule_gibbs = 0,
    birth = 0.5
  )
)$control
expect_equal(explicit@proposal.probs[["birth_death"]], 0.6)
expect_equal(explicit@proposal.probs[["swap"]], 0)
expect_equal(explicit@proposal.probs[["change"]], 0.4)
expect_equal(explicit@proposal.probs[["perturb"]], 0)
expect_equal(explicit@proposal.probs[["rule_gibbs"]], 0)

# and so does the four-name spelling that omits both moves shipping at zero,
# which is the one every consumer forwarding the documented default writes
omitted <- fit(c(
  birth_death = 0.6,
  swap = 0,
  change = 0.4,
  birth = 0.5
))$control
expect_equal(omitted@proposal.probs[["perturb"]], 0)
expect_equal(omitted@proposal.probs[["rule_gibbs"]], 0)

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
  updateState = FALSE,
  proposal.probs = threeMove
)
sampler <- dbarts::dbarts(x, y, control = oneTree)
expect_equal(sampler$control@proposal.probs[["birth_death"]], 0.5)
expect_equal(sampler$control@proposal.probs[["swap"]], 0.1)
expect_equal(sampler$control@proposal.probs[["change"]], 0.4)

set.seed(9L)
samples <- sampler$run()
expect_true(all(is.finite(samples$train)))
expect_true(all(is.finite(samples$sigma)))
# the mixture the sampler was created with survives the run
expect_equal(sampler$control@proposal.probs[["swap"]], 0.1)

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
      verbose = TRUE,
      proposal.probs = threeMove
    )
  )
)
expect_true(any(grepl(
  paste(
    "birth/death 0.50, swap 0.10, change 0.40, perturb 0.00,",
    "rule_gibbs 0.00; birth 0.50"
  ),
  printed,
  fixed = TRUE
)))

# ---- the fill and the sum --------------------------------------------------

# one unnamed element takes the residual, whichever it is
partial <- fit(c(birth_death = 0.7, swap = 0))$control
expect_equal(partial@proposal.probs[["birth_death"]], 0.7)
expect_equal(partial@proposal.probs[["swap"]], 0)
expect_equal(partial@proposal.probs[["change"]], 0.3)
expect_equal(
  fit(c(birth_death = 0.5, change = 0.4))$control@proposal.probs[["swap"]],
  0.1
)

# two unnamed, one of them swap: swap takes its zero and the other the residual
bdOnly <- fit(c(birth_death = 0.7))$control
expect_equal(bdOnly@proposal.probs[["birth_death"]], 0.7)
expect_equal(bdOnly@proposal.probs[["swap"]], 0)
expect_equal(bdOnly@proposal.probs[["change"]], 0.3)
changeOnly <- fit(c(change = 0.25))$control
expect_equal(changeOnly@proposal.probs[["birth_death"]], 0.75)
expect_equal(changeOnly@proposal.probs[["swap"]], 0)
expect_equal(changeOnly@proposal.probs[["change"]], 0.25)

# a zero-default move named alone leaves the birth/death-versus-change split
# undetermined, whichever one it is
expect_error(fit(c(swap = 0.1)), "name at least one of")
expect_error(fit(c(perturb = 0.16)), "name at least one of")
expect_error(fit(c(swap = 0.1, perturb = 0.16)), "name at least one of")
expect_error(fit(c(rule_gibbs = 0.16)), "name at least one of")
expect_error(
  fit(c(perturb = 0.16, rule_gibbs = 0.16)),
  "name at least one of"
)

# all three unnamed falls back to the default
expect_equal(fit(c(birth = 0.25))$control@proposal.probs[["birth_death"]], 0.6)
expect_equal(fit(c(birth = 0.25))$control@proposal.probs[["swap"]], 0)
expect_equal(fit(c(birth = 0.25))$control@proposal.probs[["change"]], 0.4)
expect_equal(fit(c(birth = 0.25))$control@proposal.probs[["perturb"]], 0)
expect_equal(fit(c(birth = 0.25))$control@proposal.probs[["rule_gibbs"]], 0)
expect_equal(fit(c(birth = 0.25))$control@proposal.probs[["birth"]], 0.25)

# perturb resolves BEFORE the fill: it takes its zero rather than the
# residual, so the residual is taken against 1 - perturb
withPerturb <- fit(c(birth_death = 0.5, change = 0.34, perturb = 0.16))$control
expect_equal(withPerturb@proposal.probs[["swap"]], 0)
expect_equal(withPerturb@proposal.probs[["perturb"]], 0.16)
twoUnnamed <- fit(c(birth_death = 0.84, perturb = 0.16))$control
expect_equal(twoUnnamed@proposal.probs[["swap"]], 0)
expect_equal(twoUnnamed@proposal.probs[["change"]], 0)
expect_equal(twoUnnamed@proposal.probs[["perturb"]], 0.16)

# and the same for rule_gibbs, alone and beside perturb: the residual is
# taken against 1 minus their sum
withGibbs <- fit(c(birth_death = 0.5, change = 0.34, rule_gibbs = 0.16))$control
expect_equal(withGibbs@proposal.probs[["swap"]], 0)
expect_equal(withGibbs@proposal.probs[["rule_gibbs"]], 0.16)
bothZeroDefault <- fit(
  c(birth_death = 0.5, change = 0.18, perturb = 0.16, rule_gibbs = 0.16)
)$control
expect_equal(bothZeroDefault@proposal.probs[["swap"]], 0)
expect_equal(bothZeroDefault@proposal.probs[["perturb"]], 0.16)
expect_equal(bothZeroDefault@proposal.probs[["rule_gibbs"]], 0.16)

# all four named must sum to one
expect_error(
  fit(c(birth_death = 0.7, swap = 0.1, change = 0.4)),
  "sum to 1"
)
expect_error(
  fit(c(birth_death = 0.6, swap = 0, change = 0.4, perturb = 0.1)),
  "sum to 1"
)
expect_error(
  fit(c(birth_death = 0.6, swap = 0, change = 0.4, rule_gibbs = 0.1)),
  "sum to 1"
)

# ---- the monotone rewrite --------------------------------------------------

# a defaulted mixture is rewritten birth/death-only; a non-default one errors
forced <- fit(
  c(birth_death = 0.6, swap = 0, change = 0.4, birth = 0.5),
  monotone = c(a = "+")
)$control
expect_equal(forced@proposal.probs[["birth_death"]], 1)
expect_equal(forced@proposal.probs[["swap"]], 0)
expect_equal(forced@proposal.probs[["change"]], 0)
expect_equal(forced@proposal.probs[["perturb"]], 0)
expect_equal(forced@proposal.probs[["rule_gibbs"]], 0)
expect_error(fit(threeMove, monotone = c(a = "+")), "proposal.probs")

# the refusal must not fire on a caller who spells the documented default and
# omits the move that ships at zero: the comparison fills it first
forcedFull <- fit(
  c(
    birth_death = 0.6,
    swap = 0,
    change = 0.4,
    perturb = 0,
    rule_gibbs = 0,
    birth = 0.5
  ),
  monotone = c(a = "+")
)$control
expect_equal(forcedFull@proposal.probs[["birth_death"]], 1)
expect_equal(forcedFull@proposal.probs[["perturb"]], 0)
expect_equal(forcedFull@proposal.probs[["rule_gibbs"]], 0)

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
  updateState = FALSE,
  proposal.probs = c(
    birth_death = 0.2,
    swap = 0,
    change = 0,
    perturb = 0.8,
    birth = 0.5
  )
)
perturbSampler <- dbarts::dbarts(x, y, control = perturbControl)
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

# ---- the frozen mixture ----------------------------------------------------

# All four structural probabilities exactly zero is legal and proposes no
# structure at all: the trees stand where they are while the leaf values and
# sigma keep being drawn, which is how a fitted forest is re-sampled as a
# fixed basis.
frozen <- c(birth_death = 0, swap = 0, change = 0, perturb = 0)
frozenModel <- fit(frozen)$control
expect_equal(frozenModel@proposal.probs[["birth_death"]], 0)
expect_equal(frozenModel@proposal.probs[["swap"]], 0)
expect_equal(frozenModel@proposal.probs[["change"]], 0)
expect_equal(frozenModel@proposal.probs[["perturb"]], 0)
# birth is unread when nothing is proposed and keeps its default
expect_equal(frozenModel@proposal.probs[["birth"]], 0.5)

# the residual is a share of structural mass, and an all-zero mixture has
# none: an unnamed swap keeps its zero rather than taking the whole of it
bothZero <- fit(c(birth_death = 0, change = 0))$control
expect_equal(bothZero@proposal.probs[["birth_death"]], 0)
expect_equal(bothZero@proposal.probs[["swap"]], 0)
expect_equal(bothZero@proposal.probs[["change"]], 0)
expect_equal(bothZero@proposal.probs[["perturb"]], 0)

# an unnamed birth/death or change still takes the residual, so a single
# named zero fills exactly as before
singleZero <- fit(c(change = 0))$control
expect_equal(singleZero@proposal.probs[["birth_death"]], 1)
expect_equal(singleZero@proposal.probs[["swap"]], 0)
expect_equal(singleZero@proposal.probs[["change"]], 0)
expect_equal(singleZero@proposal.probs[["perturb"]], 0)

# only exact zero freezes: a mixture that merely rounds to nothing is still a
# mixture, and must sum to one
expect_error(
  fit(c(birth_death = 0, swap = 0, change = 1e-8, perturb = 0)),
  "sum to 1"
)

frozenControl <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 3L,
  n.burn = 50L,
  n.samples = 25L,
  keepTrees = TRUE,
  updateState = FALSE
)
frozenSampler <- dbarts::dbarts(x, y, control = frozenControl)
frozenAfter <- frozenControl
frozenAfter@proposal.probs <- frozenModel@proposal.probs
set.seed(23L)
invisible(frozenSampler$run())
frozenSampler$setControl(frozenAfter)
frozenSamples <- frozenSampler$run(0L, 25L)

frozenTrees <- frozenSampler$getTrees()
frozenDraws <- split(
  frozenTrees[, c("tree", "n", "var", "value")],
  frozenTrees$sample
)
# variables, cut points and node counts, but not the leaf values
frozenStructure <- vapply(
  frozenDraws,
  function(d) {
    paste(
      d$tree,
      d$n,
      d$var,
      ifelse(d$var != -1L, d$value, NA_real_),
      collapse = "|"
    )
  },
  character(1L)
)
expect_equal(length(unique(frozenStructure)), 1L)
# and the structure that stood was not the trivial all-root one
expect_true(any(frozenTrees$var != -1L))

# the leaf values and sigma are still moving
frozenLeaves <- vapply(
  frozenDraws,
  function(d) paste(d$value[d$var == -1L], collapse = "|"),
  character(1L)
)
expect_equal(length(unique(frozenLeaves)), length(frozenDraws))
expect_true(length(unique(as.vector(frozenSamples$sigma))) > 1L)
expect_true(all(is.finite(frozenSamples$train)))

# and this is the mixture that switches the level-fibre Gibbs step on:
# dbartsControl's levelGibbs defaults to NA, which takes the step for a
# forest exactly where that forest's structures are frozen. The freeze
# arrives by setControl after fifty growing sweeps, so the decision is taken
# per sweep and not once at creation - the run above and the arm below stand
# at the same forest when it lands, and part only from the sweep after it.
frozenOff <- frozenControl
frozenOff@levelGibbs <- FALSE
frozenOffSampler <- dbarts::dbarts(x, y, control = frozenOff)
frozenOffAfter <- frozenOff
frozenOffAfter@proposal.probs <- frozenModel@proposal.probs
set.seed(23L)
invisible(frozenOffSampler$run())
frozenOffSampler$setControl(frozenOffAfter)
expect_false(
  identical(frozenOffSampler$run(0L, 25L)$train, frozenSamples$train)
)
rm(frozenOff, frozenOffAfter, frozenOffSampler)

# ---- the default replays bitwise ------------------------------------------

# the fill's new branch must not touch the shipped default: spelling it out
# reproduces the unset one draw for draw
replayControl <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 5L,
  n.burn = 10L,
  n.samples = 20L,
  updateState = FALSE
)
# a sampler takes its seed off the R stream at creation, so the seed is set
# ahead of each one rather than ahead of each run
set.seed(59L)
replayA <- dbarts::dbarts(x, y, control = replayControl)$run()
set.seed(59L)
replayControlSpelled <- replayControl
replayControlSpelled@proposal.probs <- dbarts:::resolveProposalProbs(c(
  birth_death = 0.6,
  swap = 0,
  change = 0.4,
  perturb = 0,
  birth = 0.5
))
replayB <- dbarts::dbarts(x, y, control = replayControlSpelled)$run()
expect_identical(replayA$train, replayB$train)
expect_identical(replayA$sigma, replayB$sigma)

# ---- a rule_gibbs-dominant run ---------------------------------------------

# One tree on a mixed design, most of its proposals exact rule draws. Only one
# proposal is made per tree per sweep, so consecutive draws holding the node
# count fixed took a rule draw or nothing: the shape must agree exactly, every
# rule that moved must sit at a nog node - an interior node whose two children
# are both leaves - and no categorical rule may move at all, a node carrying
# one being a fixed point of this move.
gibbsX <- data.frame(
  a = runif(n),
  b = runif(n),
  f = factor(sample(letters[1:3L], n, replace = TRUE))
)
gibbsY <- gibbsX$a -
  0.5 * gibbsX$b +
  as.integer(gibbsX$f) / 3 +
  rnorm(n, 0, 0.2)
gibbsControl <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 1L,
  n.burn = 0L,
  n.samples = 80L,
  keepTrees = TRUE,
  updateState = FALSE
)
gibbsControl@proposal.probs <- dbarts:::resolveProposalProbs(c(
  birth_death = 0.2,
  swap = 0,
  change = 0,
  perturb = 0,
  rule_gibbs = 0.8,
  birth = 0.5
))
gibbsSampler <- dbarts::dbarts(gibbsY ~ ., gibbsX, control = gibbsControl)
set.seed(43L)
gibbsSamples <- gibbsSampler$run()
expect_true(all(is.finite(gibbsSamples$train)))
expect_true(all(is.finite(gibbsSamples$sigma)))

# the nog nodes of a depth-first (parent, left subtree, right subtree) listing,
# read off the leaf pattern alone
nogNodes <- function(isLeaf) {
  nog <- logical(length(isLeaf))
  position <- 1L
  walk <- function() {
    here <- position
    position <<- position + 1L
    if (isLeaf[here]) {
      return(TRUE)
    }
    leftIsLeaf <- walk()
    rightIsLeaf <- walk()
    nog[here] <<- leftIsLeaf && rightIsLeaf
    FALSE
  }
  walk()
  nog
}

sameOrBothMissing <- function(before, after) {
  (is.na(before) & is.na(after)) |
    (!is.na(before) & !is.na(after) & before == after)
}

gibbsTrees <- gibbsSampler$getTrees()
byGibbsDraw <- split(gibbsTrees, gibbsTrees$sample)
heldShape <- 0L
ruleMoves <- 0L
categoricalRules <- 0L
for (i in seq_len(length(byGibbsDraw) - 1L)) {
  before <- byGibbsDraw[[i]]
  after <- byGibbsDraw[[i + 1L]]
  categoricalRules <- categoricalRules +
    sum(before$var != -1L & is.na(before$value))
  if (nrow(before) != nrow(after)) {
    next
  }
  heldShape <- heldShape + 1L
  # one proposal per sweep, so an unchanged node count is an unchanged shape
  expect_identical(before$var == -1L, after$var == -1L)
  interior <- before$var != -1L
  # a categorical rule reports its split in 'directions' and no cut value, an
  # ordinal one the other way round, so every column decides and each of them
  # is NA where the other kind of rule (or a leaf) sits
  sameRule <- sameOrBothMissing(before$var, after$var) &
    sameOrBothMissing(before$value, after$value)
  if (!is.null(before$directions)) {
    sameRule <- sameRule &
      sameOrBothMissing(before$directions, after$directions)
  }
  moved <- interior & !sameRule
  expect_true(sum(moved) <= 1L)
  expect_true(all(nogNodes(before$var == -1L)[moved]))
  expect_true(!any(is.na(before$value[moved])))
  ruleMoves <- ruleMoves + sum(moved)
}
# the run exercised both arms, and the design put categorical rules in front of
# the move for it to leave alone
expect_true(length(unique(vapply(byGibbsDraw, nrow, 0L))) > 1L)
expect_true(heldShape > 10L)
expect_true(ruleMoves > 0L)
expect_true(categoricalRules > 0L)
