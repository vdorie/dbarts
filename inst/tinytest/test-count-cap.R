# The count cap. A logistic count weight and a fixed negative-binomial shape
# cost one Polya-Gamma draw per unit per row per sweep, so both are refused
# above 1e6, the bound a negative-binomial response already carries, on every
# surface that takes them; the bound itself is in. A user interrupt reaches
# inside a sweep's Polya-Gamma draws and leaves the sampler usable.

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
weightRefusal <- "logistic 'weights' are observation counts and must be no larger than 1000000"

# --- logistic weights: creation, $setWeights, $setData ---
expect_error(
  dbarts(x, y, weights = overCap, family = "logistic", control = control),
  weightRefusal
)
sampler <- dbarts(
  x,
  y,
  weights = atCap,
  family = "logistic",
  control = control
)
expect_identical(sampler$data@weights, atCap)
expect_error(sampler$setWeights(overCap), weightRefusal)
sampler$setWeights(atCap)
expect_error(
  sampler$setData(dbartsData(x, y, weights = overCap)),
  weightRefusal
)
sampler$setData(dbartsData(x, y, weights = atCap))
expect_identical(sampler$data@weights, atCap)

# --- a fixed nbinom shape at creation; a drawn one at a state install ---
shapeRefusal <- "nbinom 'shape' must be no larger than 1000000"
expect_error(
  dbarts(x, counts, family = nbinom(shape = 1e6 + 1), control = control),
  shapeRefusal
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
  "dbartsDrawLatents: 'weights' must be counts no larger than 1000000"
)
expect_error(
  dbartsDrawLatents("nbinom", 0, 2, shape = 1e6 + 1),
  "dbartsDrawLatents: 'shape' must be no larger than 1000000"
)
expect_error(
  dbartsWorkingResponse("nbinom", 1, 2, shape = 1e6 + 1),
  "dbartsWorkingResponse: 'shape' must be no larger than 1000000"
)
expect_true(dbartsDrawLatents("logistic", 0, 1, weights = 1e6) > 0)
expect_true(dbartsDrawLatents("nbinom", 0, 2, shape = 1e6) > 0)

# --- an interrupt inside the first sweep's Polya-Gamma draws ---
# The run polls at once as its first sweep starts and then no sooner than
# 100 ms later; at 2e5 draws per row the first sweep's refresh runs well past
# that, so the hook's second poll lands inside it. The rows it had not reached
# keep their cold start, w / 4, where an interrupt between sweeps would have
# found every row redrawn.
pollHook <- function(after) {
  invisible(.Call(
    dbarts:::C_dbarts_bartcore_setMonotoneCountHooks,
    NA_real_,
    FALSE,
    as.integer(after)
  ))
}
heavy <- rep(2e5, n)
slow <- dbarts(x, y, weights = heavy, family = "logistic", control = control)
pollHook(2L)
message <- tryCatch(
  {
    slow$run(0L, 1L)
    "not interrupted"
  },
  error = conditionMessage,
  finally = pollHook(0L)
)
expect_true(grepl("sampler run interrupted", message))
latents <- slow$getLatents()
coldRows <- latents == heavy / 4
expect_true(any(coldRows) && !all(coldRows))
expect_true(all(is.finite(latents) & latents > 0))
# the interrupted state stores and installs, and the chain runs on from it
slow$storeState()
slow$setState(slow$state)
expect_identical(slow$getLatents(), latents)
expect_true(is.list(slow$run(0L, 1L)))
expect_false(any(slow$getLatents() == heavy / 4))

# --- multinomial: a count row's trial total, at creation and $setCounts ---
trialRefusal <- "multinomial 'counts' rows must total no more than 1000000 trials"
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
  trialRefusal
)
expect_error(
  bart(x, overTrialCap, family = "multinomial", control = mnControl),
  trialRefusal
)
multinomial <- dbarts(
  dbartsData(x, counts = atTrialCap),
  family = "multinomial",
  control = mnControl
)
expect_error(multinomial$setCounts(overTrialCap), trialRefusal)
expect_identical(multinomial$data@counts, atTrialCap)
multinomial$setCounts(atTrialCap)
