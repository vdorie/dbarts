# The three hand-picked 1e-10 almost-equal checks on user-supplied
# probability vectors: the tree-prior rule proposal probabilities
# (R/A_class.R dbartsControl validity
# and the bartcore bridge's parseProposalProbs) and the split probabilities
# (the bridge only) move to the conventional almost-equal tolerance,
# sqrt(DBL_EPSILON) - R's own all.equal default. A vector mis-normalized by
# 1e-9 (refused before this change) is now ACCEPTED; one off by 1e-7 stays
# refused. Error messages are unchanged. The CGM prior's own splitProbabilities
# simplex check (R/A_class.R, dbartsCGMPrior validity) is a separate,
# untouched 1e-10 check outside this file's scope; the bridge sites below are
# reached by direct slot assignment, which bypasses R validity entirely and
# so cannot trip it.

set.seed(19)
n <- 60L
p <- 3L
x <- matrix(runif(n * p), n, p)
y <- rowSums(x) + rnorm(n, sd = 0.2)
control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 5L,
  updateState = FALSE
)
sampler <- dbarts(x, y, control = control)
create <- function(control = sampler$control, model = sampler$model) {
  .Call(dbarts:::C_dbarts_bartcore_create, control, model, sampler$data, "")
}

# --- R/A_class.R dbartsControl validity: rule proposal probabilities ------
makeControl <- function(delta) {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    proposal.probs = c(
      birth_death = 1.0 - delta,
      swap = 0.0,
      change = 0.0,
      birth = 0.5
    )
  )
}
expect_inherits(makeControl(1e-9), "dbartsControl")
expect_error(makeControl(1e-7), "rule proposal probabilities must sum to 1")

# the same pair with every structural name spelled, the two zero-default
# moves included: an unnamed perturb or rule_gibbs is resolved to zero ahead
# of the fill and never takes the residual, so naming them changes neither
# verdict
makeFullControl <- function(delta) {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    proposal.probs = c(
      birth_death = 1.0 - delta,
      swap = 0.0,
      change = 0.0,
      perturb = 0.0,
      rule_gibbs = 0.0,
      birth = 0.5
    )
  )
}
expect_inherits(makeFullControl(1e-9), "dbartsControl")
expect_error(makeFullControl(1e-7), "rule proposal probabilities must sum to 1")

# --- bartcore bridge, parseProposalProbs: rule proposal probabilities -----
# Reached directly via a mutated copy of a valid control, bypassing R
# validity (plain slot assignment does not call it) so the C-side check is
# pinned in isolation.
pushProposal <- function(delta) {
  control <- sampler$control
  control@proposal.probs[["birth_death"]] <-
    control@proposal.probs[["birth_death"]] - delta
  control
}
expect_silent(create(pushProposal(1e-9)))
expect_error(
  create(pushProposal(1e-7)),
  "rule proposal probabilities must sum to 1.0"
)

# --- bartcore bridge, parseModel: split probabilities ----------------------
splitPrior <- methods::new(
  "dbartsCGMPrior",
  power = 2,
  base = 0.95,
  splitProbabilities = rep(1 / p, p)
)
pushSplit <- function(delta) {
  model <- sampler$model
  model@tree.prior <- splitPrior
  model@tree.prior@splitProbabilities[1L] <-
    model@tree.prior@splitProbabilities[1L] + delta
  model
}
expect_silent(create(model = pushSplit(1e-9)))
expect_error(
  create(model = pushSplit(1e-7)),
  "split probabilities must sum to 1.0"
)
