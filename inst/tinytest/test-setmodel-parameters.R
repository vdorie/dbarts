# $setModel installs the parameters a model carries beyond the tree and leaf
# priors on every chain: a Student-t df, a negative-binomial shape, a gp
# leaf's lengthscale and DART's settings. Before any run a model draws what
# creating with it draws, except that one handed from held to drawn starts at
# the value in force; mid-run the live sampler, its copy and its reload draw
# the same.

student <- dbartsFamilies$student
nbinom <- dbartsFamilies$nbinom
dart <- dbartsPriors$dart
gp <- dbartsPriors$gp

set.seed(23L)
n <- 80L
x <- matrix(runif(n * 3L), n, 3L)
y <- 2 * x[, 1L] + x[, 2L] + rnorm(n, sd = 0.4)
yCount <- rpois(n, exp(0.5 + x[, 1L]))

makeControl <- function(keepTrees) {
  dbartsControl(
    n.chains = 2L,
    n.threads = 1L,
    n.samples = 10L,
    seed = 9L,
    updateState = FALSE,
    keepTrees = keepTrees
  )
}
build <- function(..., response = y, keepTrees = FALSE) {
  dbarts(x, response, ..., control = makeControl(keepTrees), n.trees = 10L)
}
maker <- function(...) function() build(...)
stateOf <- function(sampler, name) {
  sampler$storeState()
  vapply(sampler$state, function(chain) chain[[name]], 0)
}

# before any run, the model draws what creating with it draws
asCreated <- function(make, made) {
  sampler <- make()
  created <- made()
  sampler$setModel(created$model)
  expect_identical(sampler$run(0L, 15L)$train, created$run(0L, 15L)$train)
  sampler
}
# mid-run, the live sampler, its copy and its reload draw the same, to the
# order in which fits are summed
agreeAfter <- function(make, model) {
  sampler <- make()
  invisible(sampler$run(0L, 10L))
  sampler$setModel(model)
  sampler$storeState()
  copied <- sampler$copy()
  reloaded <- unserialize(serialize(sampler, NULL))
  live <- sampler$run(0L, 10L)$train
  expect_equal(copied$run(0L, 10L)$train, live, tolerance = 1e-12)
  expect_equal(reloaded$run(0L, 10L)$train, live, tolerance = 1e-12)
  sampler
}

# ---- Student-t degrees of freedom -------------------------------------------

held3 <- maker(family = student(df = 3))
drawnDf <- maker(family = student())
asCreated(held3, maker(family = student(df = 30)))
fixed <- asCreated(drawnDf, held3)
expect_true(all(fixed$run(0L, 10L)$resid.df == 3))
# held to drawn starts at the value in force, then draws
handed <- held3()
handed$setModel(drawnDf()$model)
expect_identical(stateOf(handed, "resid.df"), c(3, 3))
expect_true(length(unique(as.vector(handed$run(0L, 30L)$resid.df))) > 1L)
agreeAfter(held3, build(family = student(df = 30))$model)
agreeAfter(held3, drawnDf()$model)
agreeAfter(drawnDf, held3()$model)

# ---- negative-binomial shape ------------------------------------------------

held4 <- maker(family = nbinom(shape = 4), response = yCount)
drawnShape <- maker(family = nbinom(), response = yCount)
held10 <- maker(family = nbinom(shape = 10), response = yCount)
changed <- asCreated(held4, held10)
expect_identical(changed$model@shape, 10)
fixed <- asCreated(drawnShape, held4)
expect_true(all(fixed$run(0L, 10L)$shape == 4))
handed <- held4()
handed$setModel(drawnShape()$model)
expect_identical(handed$model@shape, -1)
expect_identical(stateOf(handed, "shape"), c(4, 4))
expect_true(length(unique(as.vector(handed$run(0L, 30L)$shape))) > 1L)
agreeAfter(held4, held10()$model)
agreeAfter(held4, drawnShape()$model)
agreeAfter(drawnShape, held4()$model)

# ---- gp lengthscale ---------------------------------------------------------

lengthscalesOf <- function(sampler) {
  sampler$storeState()
  sampler$state[[1L]]$forests[[1L]]$leaf.lengthscales
}
heuristic <- maker(leaf.prior = gp(1:2))
stated <- maker(leaf.prior = gp(1:2, lengthscale = 0.3))
asCreated(heuristic, stated)
asCreated(stated, maker(leaf.prior = gp(1:2, lengthscale = c(0.5, 2))))
# a stated lengthscale, then a model stating none: the stated one stays in
# force and on the record a copy re-creates from
statesNone <- agreeAfter(heuristic, stated()$model)
statesNone$setModel(heuristic()$model)
expect_identical(statesNone$model@leaf.prior@lengthscale, c(0.3, 0.3))
expect_identical(lengthscalesOf(statesNone), c(0.3, 0.3))
agreeAfter(
  function() {
    sampler <- stated()
    sampler$setModel(heuristic()$model)
    sampler
  },
  heuristic()$model
)
# a new lengthscale is refused while kept draws are held, before anything
# else in the model goes in: the sampler is left as a twin never asked
withDraws <- function() {
  sampler <- build(leaf.prior = gp(1:2), keepTrees = TRUE)
  invisible(sampler$run(0L, 5L))
  sampler
}
sampler <- withDraws()
twin <- withDraws()
before <- sampler$model
model <- stated()$model
model@tree.prior@power <- 3
expect_error(
  sampler$setModel(model),
  "cannot change a gp leaf's lengthscale on a sampler holding kept draws"
)
expect_identical(sampler$model, before)
expect_identical(lengthscalesOf(sampler), lengthscalesOf(twin))
expect_identical(sampler$run(0L, 5L)$train, twin$run(0L, 5L)$train)
# a model stating none is taken after the data move, kept draws and all, the
# lengthscale in force unchanged
inForce <- lengthscalesOf(sampler)
sampler$setPredictor(x[, 2L] * 1.01, 2L)
sampler$setModel(heuristic()$model)
sampler$setResponse(y + 0.1)
model <- heuristic()$model
model@tree.prior@power <- 3
sampler$setModel(model)
expect_identical(lengthscalesOf(sampler), inForce)
expect_null(sampler$model@leaf.prior@lengthscale)

# ---- DART's settings --------------------------------------------------------

# update.delay = 0 throughout: at the default (half of n.burn) a run this short
# never updates DART, and the pins would compare the initial state with itself
dart0 <- function(...) dart(..., update.delay = 0)
drawnAlpha <- maker(split.probs = dart0())
for (made in list(
  maker(split.probs = dart0(a = 0.7, b = 2)),
  maker(split.probs = dart0(rho = 2)),
  maker(split.probs = dart0(alpha = 3)),
  maker(split.probs = dart(update.delay = 1)),
  maker(split.probs = dart0(alpha = 2, update.alpha = FALSE))
)) {
  asCreated(drawnAlpha, made)
}
heldAlpha <- maker(split.probs = dart0(alpha = 2, update.alpha = FALSE))
asCreated(
  heldAlpha,
  maker(split.probs = dart0(alpha = 5, update.alpha = FALSE))
)
# held to drawn starts at the value in force; a drawn alpha restarts at a new
# stated start and keeps its value under the same one
handed <- heldAlpha()
handed$setModel(build(split.probs = dart0(alpha = 5))$model)
expect_identical(stateOf(handed, "dart.alpha"), c(2, 2))
sampler <- drawnAlpha()
invisible(sampler$run(0L, 30L))
expect_true(all(stateOf(sampler, "dart.alpha") != 1))
alpha <- stateOf(sampler, "dart.alpha")
sampler$setModel(build(split.probs = dart0(a = 0.7))$model)
expect_identical(stateOf(sampler, "dart.alpha"), alpha)
sampler$setModel(build(split.probs = dart0(a = 0.7, alpha = 3))$model)
expect_identical(stateOf(sampler, "dart.alpha"), c(3, 3))
# the split weights and the count of skipped updates are kept mid-run
sampler <- drawnAlpha()
invisible(sampler$run(0L, 30L))
sampler$storeState()
weights <- lapply(sampler$state, `[[`, "dart.probabilities")
expect_true(all(vapply(weights, function(w) any(abs(w - 1 / 3) > 1e-6), TRUE)))
skipped <- lapply(sampler$state, `[[`, "dart.updates.skipped")
sampler$setModel(build(split.probs = dart(b = 3, update.delay = 50))$model)
sampler$storeState()
expect_identical(lapply(sampler$state, `[[`, "dart.probabilities"), weights)
expect_identical(lapply(sampler$state, `[[`, "dart.updates.skipped"), skipped)
agreeAfter(drawnAlpha, build(split.probs = dart0(a = 0.7, rho = 2))$model)
agreeAfter(drawnAlpha, heldAlpha()$model)
agreeAfter(heldAlpha, drawnAlpha()$model)

# a non-finite alpha or rho is refused by the constructor, and by $setModel
# before anything is written (a record built by slot passes the constructor)
expect_error(
  dart(alpha = Inf),
  "'alpha' must be a single positive finite number"
)
expect_error(dart(rho = Inf), "'rho' must be a single positive finite number")
expect_error(dart(a = Inf), "'a' must be a single positive finite number")
expect_error(dart(b = Inf), "'b' must be a single positive finite number")
expect_error(dart(a = NaN), "'a' must be a single positive finite number")
expect_error(
  dart(alpha = NA_real_),
  "'alpha' must be a single positive finite number"
)
for (arm in c("alpha", "rho", "a", "b")) {
  sampler <- drawnAlpha()
  twin <- drawnAlpha()
  invisible(sampler$run(0L, 10L))
  invisible(twin$run(0L, 10L))
  model <- sampler$model
  if (arm == "alpha") {
    model@tree.prior@alpha <- Inf
  } else if (arm == "rho") {
    model@tree.prior@rho <- Inf
  } else if (arm == "a") {
    model@tree.prior@a <- Inf
  } else {
    model@tree.prior@b <- Inf
  }
  model@tree.prior@power <- 3
  expect_error(sampler$setModel(model), "must be finite")
  expect_identical(sampler$model, twin$model)
  expect_identical(sampler$run(0L, 10L)$train, twin$run(0L, 10L)$train)
}

# ---- negative-binomial shape: the unchanged model moves nothing --------------

for (made in list(held4, drawnShape)) {
  sampler <- made()
  twin <- made()
  invisible(sampler$run(0L, 10L))
  invisible(twin$run(0L, 10L))
  sampler$setModel(sampler$model)
  expect_identical(sampler$run(0L, 10L)$train, twin$run(0L, 10L)$train)
}
