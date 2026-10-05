# The count cap. A logistic count weight and a fixed negative-binomial shape
# cost one Polya-Gamma draw per unit per row per sweep, so both are refused
# above 1e6, the bound a negative-binomial response already carries, on every
# surface that takes them, as is a multinomial count row totalling more
# trials; the bound itself is in. The interrupt inside a sweep's draws cannot
# be simulated here without racing the run's 100 ms poll throttle, so
# tests/cpp drives it with an injected cancel.

set.seed(5309L, sample.kind = "Rejection")
n <- 20L
x <- matrix(runif(n), n)
y <- as.double(rbinom(n, 1L, 0.5))
counts <- as.double(rpois(n, 3))
control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 5L,
  updateState = FALSE,
  seed = 17L
)
atCap <- replace(rep(1, n), 1L, 1e6)
overCap <- replace(rep(1, n), 1L, 1e6 + 1)
# every cap refusal names its entry, the family, the argument and the cap
refusal <- function(caller, family, what) {
  sprintf(
    "%s: family \"%s\" requires %s no larger than 1000000",
    caller,
    family,
    what
  )
}
weightRefusal <- function(caller) refusal(caller, "logistic", "'weights'")

# --- logistic weights: creation, $setWeights, $setData ---
expect_error(
  dbarts(x, y, weights = overCap, family = "logistic", control = control),
  weightRefusal("sampler creation"),
  fixed = TRUE
)
sampler <- dbarts(
  x,
  y,
  weights = atCap,
  family = "logistic",
  control = control
)
expect_identical(sampler$data@weights, atCap)
expect_error(
  sampler$setWeights(overCap),
  weightRefusal("$setWeights"),
  fixed = TRUE
)
sampler$setWeights(atCap)
expect_error(
  sampler$setData(dbartsData(x, y, weights = overCap)),
  weightRefusal("$setData"),
  fixed = TRUE
)
sampler$setData(dbartsData(x, y, weights = atCap))
expect_identical(sampler$data@weights, atCap)

# --- a fixed nbinom shape at creation; a drawn one at a state install ---
expect_error(
  dbarts(x, counts, family = nbinom(shape = 1e6 + 1), control = control),
  refusal("sampler creation", "nbinom", "a fixed 'shape'"),
  fixed = TRUE
)
fixed <- dbarts(x, counts, family = nbinom(shape = 1e6), control = control)
expect_equal(fixed$getShape(), 1e6)

drawn <- dbarts(
  x,
  counts,
  family = "nbinom",
  control = dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 5L,
    updateState = TRUE,
    seed = 17L
  )
)
invisible(drawn$run(5L, 1L))
state <- drawn$state
state[[1L]]$shape <- 1e6
drawn$setState(state)
expect_equal(drawn$getShape(), 1e6)
state[[1L]]$shape <- 1e6 + 1
expect_error(drawn$setState(state), "not consistent with this sampler")

# --- the augmentation helpers draw the same variates, so take the same cap ---
expect_error(
  dbartsDrawLatents("logistic", 0, 1, weights = 1e6 + 1),
  refusal("dbartsDrawLatents", "logistic", "'weights'"),
  fixed = TRUE
)
expect_error(
  dbartsDrawLatents("nbinom", 0, 2, shape = 1e6 + 1),
  refusal("dbartsDrawLatents", "nbinom", "'shape'"),
  fixed = TRUE
)
expect_error(
  dbartsWorkingResponse("nbinom", 1, 2, shape = 1e6 + 1),
  refusal("dbartsWorkingResponse", "nbinom", "'shape'"),
  fixed = TRUE
)
expect_true(dbartsDrawLatents("logistic", 0, 1, weights = 1e6) > 0)
expect_true(dbartsDrawLatents("nbinom", 0, 2, shape = 1e6) > 0)

# --- multinomial: a count row's trial total, at creation and $setCounts ---
trialRefusal <- function(caller) {
  refusal(caller, "multinomial", "'counts' row totals")
}
mnCounts <- function(first) {
  counts <- matrix(1L, n, 3L)
  counts[1L, ] <- as.integer(first)
  counts
}
atTrialCap <- mnCounts(c(1e6 - 2, 1, 1))
overTrialCap <- mnCounts(c(1e6 - 1, 1, 1))
mnControl <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 5L,
  updateState = TRUE,
  seed = 17L
)
expect_error(
  dbarts(
    dbartsData(x, counts = overTrialCap),
    family = "multinomial",
    control = mnControl
  ),
  trialRefusal("sampler creation"),
  fixed = TRUE
)
expect_error(
  bart(x, overTrialCap, family = "multinomial", control = mnControl),
  trialRefusal("sampler creation"),
  fixed = TRUE
)
multinomial <- dbarts(
  dbartsData(x, counts = atTrialCap),
  family = "multinomial",
  control = mnControl
)
expect_error(
  multinomial$setCounts(overTrialCap),
  trialRefusal("$setCounts"),
  fixed = TRUE
)
expect_identical(multinomial$data@counts, atTrialCap)
multinomial$setCounts(atTrialCap)
