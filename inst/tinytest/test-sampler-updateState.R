# dbartsSampler's seven mutators (setData, setResponse, setOffset,
# setWeights, setSigma, setPredictor, setCutPoints) wire updateState as an
# opt-in store: explicit TRUE calls storeState() after a successful
# mutation; NA (the default) and FALSE store nothing, even when the
# sampler's control@updateState is TRUE - contrast run()'s (and
# sampleTreesFromPrior's/sampleNodeParametersFromPrior's) NA ->
# control@updateState convention. This matters only once $state has
# already been forced (read or stored) at least once; an unforced state
# promise always materializes CURRENT (post-mutation) state on first
# access regardless (see man/dbartsSampler-class.Rd).

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

sampler$setCutPoints(list(c(0.5)), 1L)
expect_identical(sampler$state, stateBefore)

sampler$setCutPoints(list(c(0.5, 0.6)), 1L, updateState = FALSE)
expect_identical(sampler$state, stateBefore)

sampler$setCutPoints(list(seq(0.1, 0.9, 0.1)), 1L, updateState = TRUE)
expect_false(identical(sampler$state, stateBefore))

rm(sampler, stateBefore)

# an unforced state promise is unaffected by the mutator default: a
# mutate-then-first-read sequence still reflects the mutation, since the
# delayedAssign fires (and captures CURRENT state) at first access
sampler2 <- dbarts::dbarts(y ~ x, control = control)
invisible(sampler2$run(50L, 5L))
sampler2$setCutPoints(list(c(0.5)), 1L) # default updateState = NA, no store
stateAfterMutation <- sampler2$state # first access: forces to current state
sampler2$setCutPoints(list(seq(0.1, 0.9, 0.1)), 1L, updateState = TRUE)
expect_false(identical(sampler2$state, stateAfterMutation))

rm(sampler2, stateAfterMutation)

# all seven mutators accept updateState = TRUE without error (a wiring
# smoke test, independent of whether content actually changes for a given
# one)
sampler3 <- dbarts::dbarts(y ~ x, control = control)
invisible(sampler3$run(10L, 2L))
expect_silent(sampler3$setData(dbarts::dbartsData(y ~ x), updateState = TRUE))
expect_silent(sampler3$setResponse(y, updateState = TRUE))
expect_silent(sampler3$setOffset(0, updateState = TRUE))
expect_silent(sampler3$setWeights(rep(1, n), updateState = TRUE))
expect_silent(sampler3$setSigma(1, updateState = TRUE))
expect_silent(sampler3$setPredictor(x, updateState = TRUE))
expect_silent(
  sampler3$setCutPoints(list(seq(0.1, 0.9, 0.1)), 1L, updateState = TRUE)
)

rm(sampler3, control, x, y, n)
