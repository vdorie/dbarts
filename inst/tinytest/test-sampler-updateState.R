# dbartsSampler's mutators (setData, setResponse, setOffset, setWeights,
# setSigma, setPredictor, setCutPoints, and the rest) resolve an unset
# updateState against control@updateState, exactly as run() (and
# sampleTreesFromPrior/sampleLeafParametersFromPrior) do: NA (the default)
# stores when the control says to and skips when it says not to, and an
# explicit TRUE/FALSE overrides the control either way. This matters only
# once $state has already been forced (read or stored) at least once; an
# unforced state promise always materializes CURRENT (post-mutation) state
# on first access regardless (see man/dbartsSampler-class.Rd).

set.seed(0)
n <- 200L
x <- matrix(runif(n), n, 1)
y <- ifelse(x[, 1] > 0.5, 1, -1) + rnorm(n, 0, 0.1)

control <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 20L,
  n.samples = 5L,
  n.cuts = 50L,
  updateState = TRUE
)

# setCutPoints prunes leaves that end up empty under the new grid, so its
# effect on $state's content is not a no-op (contrast setWeights/setSigma/
# setResponse/setOffset, which mutate 'data' rather than tree/RNG structure
# and so leave $state's content identical either way, whether or not it is
# re-stored)
sampler <- dbarts::dbarts(y ~ x, control = control)
invisible(sampler$run(50L, 5L))
stateBefore <- sampler$state # forces the promise once

# control@updateState = TRUE: the default (NA) stores
sampler$setCutPoints(list(c(0.5)), 1L)
expect_false(identical(sampler$state, stateBefore))
stateAfterDefault <- sampler$state

# an explicit FALSE overrides a TRUE control and skips the store
sampler$setCutPoints(list(c(0.5, 0.6)), 1L, updateState = FALSE)
expect_identical(sampler$state, stateAfterDefault)

# an explicit TRUE stores regardless
sampler$setCutPoints(list(seq(0.1, 0.9, 0.1)), 1L, updateState = TRUE)
expect_false(identical(sampler$state, stateAfterDefault))

rm(sampler, stateBefore, stateAfterDefault)

# control@updateState = FALSE: the default (NA) stores nothing
controlNoUpdate <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 20L,
  n.samples = 5L,
  n.cuts = 50L,
  updateState = FALSE
)
sampler2 <- dbarts::dbarts(y ~ x, control = controlNoUpdate)
invisible(sampler2$run(50L, 5L))
stateBefore2 <- sampler2$state # forces the promise once

sampler2$setCutPoints(list(c(0.5)), 1L)
expect_identical(sampler2$state, stateBefore2)

# an explicit TRUE still forces a store under a FALSE control
sampler2$setCutPoints(list(seq(0.1, 0.9, 0.1)), 1L, updateState = TRUE)
expect_false(identical(sampler2$state, stateBefore2))

rm(sampler2, stateBefore2, controlNoUpdate)

# an unforced state promise is unaffected by the mutator's resolved default:
# a mutate-then-first-read sequence still reflects the mutation, since the
# delayedAssign fires (and captures CURRENT state) at first access
sampler3 <- dbarts::dbarts(y ~ x, control = control)
invisible(sampler3$run(50L, 5L))
sampler3$setCutPoints(list(c(0.5)), 1L) # default updateState = NA, control TRUE
stateAfterMutation <- sampler3$state # first access: forces to current state
sampler3$setCutPoints(list(seq(0.1, 0.9, 0.1)), 1L, updateState = TRUE)
expect_false(identical(sampler3$state, stateAfterMutation))

rm(sampler3, stateAfterMutation)

# internal setter calls store nothing of their own: a sampler built with 0/1
# probit weights (installed as the active-row mask at creation) keeps its
# state promise, so a run that skips the store still reads the post-run
# state; and setData's own updateState = FALSE is not overridden by the mask
# it reinstalls
yBin <- as.integer(y > 0)
w <- rep(1, n)
w[1:5] <- 0
sampler5 <- dbarts::dbarts(x, yBin, weights = w, control = control)
invisible(sampler5$run(20L, 5L, updateState = FALSE))
statePromised <- sampler5$state
sampler5$storeState()
expect_identical(statePromised, sampler5$state)
sampler5$setData(
  dbarts::dbartsData(x, yBin, weights = rev(w)),
  updateState = FALSE
)
expect_identical(sampler5$state, statePromised)

rm(sampler5, statePromised, yBin, w)

# all seven mutators accept updateState = TRUE without error (a wiring
# smoke test, independent of whether content actually changes for a given
# one)
sampler4 <- dbarts::dbarts(y ~ x, control = control)
invisible(sampler4$run(10L, 2L))
expect_silent(sampler4$setData(dbarts::dbartsData(y ~ x), updateState = TRUE))
expect_silent(sampler4$setResponse(y, updateState = TRUE))
expect_silent(sampler4$setOffset(0, updateState = TRUE))
expect_silent(sampler4$setWeights(rep(1, n), updateState = TRUE))
expect_silent(sampler4$setSigma(1, updateState = TRUE))
expect_silent(sampler4$setPredictor(x, updateState = TRUE))
expect_silent(
  sampler4$setCutPoints(list(seq(0.1, 0.9, 0.1)), 1L, updateState = TRUE)
)

rm(sampler4, control, x, y, n)
