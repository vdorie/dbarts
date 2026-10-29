# $setModel changes a model's parameters and refuses its structure, by name and
# before anything is stored: the sampler's model and its next draws are
# untouched by a refusal.

set.seed(11L)
n <- 90L
x <- matrix(runif(n * 3L), n, 3L)
colnames(x) <- c("a", "b", "c")
y <- 2 * x[, 1L] + x[, 2L] + rnorm(n, sd = 0.3)
yBinary <- as.integer(y > median(y))
yCount <- rpois(n, exp(0.5 + x[, 1L]))

chisq <- dbartsPriors$chisq
fixed <- dbartsPriors$fixed
chi <- dbartsPriors$chi
normal <- dbartsPriors$normal

makeControl <- function() {
  dbartsControl(
    n.chains = 2L,
    n.threads = 1L,
    n.samples = 12L,
    n.burn = 6L,
    seed = 5L,
    n.trees = 12L
  )
}
build <- function(response = y, ...) {
  dbarts(x, response, control = makeControl(), ...)
}

# the refusal names its subject, and leaves the stored model and the next
# draws as a twin that was never asked
refused <- function(make, model, pattern) {
  sampler <- make()
  twin <- make()
  before <- sampler$model
  expect_error(sampler$setModel(model), pattern, info = pattern)
  expect_identical(sampler$model, before, info = pattern)
  expect_identical(
    sampler$run(2L, 4L)$train,
    twin$run(2L, 4L)$train,
    info = pattern
  )
}

# the family, in each direction
refused(build, build(yBinary, family = probit())$model, "response family")
refused(
  function() build(yBinary, family = probit()),
  build(y)$model,
  "response family"
)
refused(
  build,
  build(family = student(df = 3))$model,
  "response family.*'gaussian'.*'student'"
)
refused(
  function() build(family = student(df = 3)),
  build()$model,
  "response family.*'student'.*'gaussian'"
)
# a model written "auto" is resolved against the sampler's response as a
# fitting function resolves it
local({
  auto <- build()$model
  auto@family <- "auto"
  expect_silent(build()$setModel(auto))
  auto <- build(yBinary, family = probit())$model
  auto@family <- "auto"
  expect_silent(build(yBinary, family = probit())$setModel(auto))
  # a gaussian fit of a 0/1 response is not what "auto" resolves it to
  gaussianBinary <- function() build(yBinary, family = gaussian())
  auto <- gaussianBinary()$model
  auto@family <- "auto"
  refused(gaussianBinary, auto, "response family.*'gaussian'.*'probit'")
})

# the leaf kind and its columns, in each direction
refused(build, build(leaf.prior = linear(1L))$model, "leaf model")
refused(build, build(leaf.prior = gp(1L))$model, "leaf model")
refused(
  function() build(leaf.prior = linear(1L)),
  build()$model,
  "leaf model"
)
refused(
  function() build(leaf.prior = gp(1L)),
  build(leaf.prior = linear(1L))$model,
  "leaf model"
)
refused(
  function() build(leaf.prior = linear(1L)),
  build(leaf.prior = linear(2L))$model,
  "leaf covariate columns"
)
refused(
  function() build(leaf.prior = gp(1L)),
  build(leaf.prior = gp(2L))$model,
  "leaf covariate columns"
)
refused(
  function() build(leaf.prior = gp(1L)),
  build(leaf.prior = gp(1L, max.leaf.size = 2L))$model,
  "max.leaf.size"
)

# interactions and blocks: added, removed, changed; an absent attribute is a
# statement of none
makeInteractions <- function() build(interactions = interactions(max.order = 2))
refused(build, makeInteractions()$model, "interaction limit")
refused(makeInteractions, build()$model, "interaction limit")
refused(
  makeInteractions,
  build(interactions = interactions(max.order = 3))$model,
  "interaction limit"
)
refused(
  makeInteractions,
  build(
    interactions = interactions(max.order = 2, forbid = list(c(1L, 2L)))
  )$model,
  "forbidden interactions"
)
makeBlocks <- function() build(blocks = blocks(list(c(1L, 2L), 3L)))
refused(build, makeBlocks()$model, "blocks")
refused(makeBlocks, build()$model, "blocks")
refused(
  makeBlocks,
  build(blocks = blocks(list(1L, c(2L, 3L))))$model,
  "blocks"
)
refused(
  makeBlocks,
  build(blocks = blocks(list(c(1L, 2L), 3L), c(8L, 4L)))$model,
  "blocks"
)

# DART against the ordinary prior, in each direction
refused(build, build(tree.prior = dart())$model, "tree prior's kind")
refused(
  function() build(tree.prior = dart()),
  build()$model,
  "tree prior's kind"
)

# the columns a forest may split on
refused(
  function() build(forests = list(forest(vars = c(1L, 2L)))),
  build()$model,
  "columns a forest may split on"
)

# monotone: the bridge's message stands
refused(
  function() build(monotone = c(1, 0, 0)),
  build()$model,
  "monotone"
)
refused(
  build,
  build(monotone = c(1, 0, 0))$model,
  "monotone"
)

# parameters the engine does not install yet are refused by name
refused(
  function() build(family = student(df = 3)),
  build(family = student(df = 30))$model,
  "Student-t degrees of freedom.*yet"
)
refused(
  function() build(family = student(df = 3)),
  build(family = student())$model,
  "Student-t degrees of freedom.*yet"
)
refused(
  function() build(family = student()),
  build(family = student(df = 3))$model,
  "Student-t degrees of freedom.*yet"
)
refused(
  function() build(yCount, family = nbinom(shape = 4)),
  build(yCount, family = nbinom(shape = 10))$model,
  "negative-binomial shape.*yet"
)
refused(
  function() build(yCount, family = nbinom(shape = 4)),
  build(yCount, family = nbinom())$model,
  "negative-binomial shape.*yet"
)
refused(
  function() build(leaf.prior = gp(1L, lengthscale = 1)),
  build(leaf.prior = gp(1L, lengthscale = 0.2))$model,
  "lengthscale.*yet"
)
dartSettings <- list(
  a = function() build(tree.prior = dart(a = 1)),
  b = function() build(tree.prior = dart(b = 2)),
  rho = function() build(tree.prior = dart(rho = 2)),
  alpha = function() build(tree.prior = dart(alpha = 2)),
  update.alpha = function() build(tree.prior = dart(update.alpha = FALSE)),
  update.delay = function() build(tree.prior = dart(update.delay = 1))
)
for (setting in names(dartSettings)) {
  refused(
    function() build(tree.prior = dart()),
    dartSettings[[setting]]()$model,
    paste0("DART's '", setting, "'.*yet")
  )
}

# a family with no residual scale takes no residual prior; the sampler's own
# stands, and a family with one takes any
local({
  make <- function() build(yBinary, family = probit())
  for (prior in list(fixed(4), chisq(10, 0.5))) {
    model <- make()$model
    model@resid.prior <- prior
    refused(make, model, "family 'probit' has no residual scale")
  }
  sampler <- make()
  expect_silent(sampler$setModel(sampler$model))
  sampler <- build()
  model <- sampler$model
  model@resid.prior <- chisq(10, 0.5)
  expect_silent(sampler$setModel(model))
  expect_identical(sampler$model@resid.prior, chisq(10, 0.5))
})

# a raw leaf specification out of step with its resolved slots, k or sd alone
local({
  model <- build()$model
  model@leaf.prior@k <- 9
  refused(build, model, "does not resolve.*setLeafPrior")
  model <- build()$model
  model@leaf.prior@prior.sd <- 0.3
  refused(build, model, "does not resolve.*setLeafPrior")
  model <- build(yBinary, family = probit())$model
  model@leaf.prior@k <- 9
  refused(
    function() build(yBinary, family = probit()),
    model,
    "does not resolve"
  )
})

# a consistent one is taken, in each spelling
for (leaf in list(normal(k = 3), normal(sd = 0.4), normal(k = chi(2, 3)))) {
  sampler <- build()
  expect_silent(sampler$setModel(build(leaf.prior = leaf)$model))
  expect_identical(sampler$model@leaf.prior, leaf)
}

# the sampler's own model is taken at every point of its life: after
# setResponse with and without updateScale, setOffset, setLeafPrior and the
# route setControl takes through setModel
for (leaf in list(normal(), normal(sd = 0.4), normal(k = chi(2, 3)))) {
  sampler <- build(leaf.prior = leaf)
  expect_silent(sampler$setModel(sampler$model))
  sampler$run(2L, 3L)
  sampler$setResponse(y + 1, updateScale = FALSE)
  expect_silent(sampler$setModel(sampler$model))
  sampler$setResponse(2 * y, updateScale = TRUE)
  expect_silent(sampler$setModel(sampler$model))
  sampler$setOffset(rep(0.5, n))
  expect_silent(sampler$setModel(sampler$model))
  sampler$setLeafPrior(normal(sd = 0.7))
  expect_silent(sampler$setModel(sampler$model))
  moves <- dbartsControl(
    n.chains = 2L,
    n.threads = 1L,
    n.samples = 12L,
    n.burn = 6L,
    seed = 5L,
    n.trees = 12L,
    proposal.probs = c(birth_death = 0.5, swap = 0.1, change = 0.4)
  )
  expect_silent(sampler$setControl(moves))
}

# a DART sampler takes a DART model of its own settings with power, base, k
# and the residual prior changed, drawing before any run what creating with it
# draws
local({
  sampler <- build(tree.prior = dart())
  created <- build(
    tree.prior = dart(power = 3, base = 0.8),
    leaf.prior = normal(k = 4),
    family = gaussian(sigma = chisq(5, 0.8))
  )
  sampler$setModel(created$model)
  expect_identical(sampler$model@tree.prior, created$model@tree.prior)
  expect_identical(sampler$model@resid.prior, chisq(5, 0.8))
  expect_identical(sampler$getK(), created$getK())
  expect_identical(sampler$run(3L, 5L)$train, created$run(3L, 5L)$train)
  # the same on the other side of a run, from the engine's current state
  again <- sampler$model
  again@tree.prior@power <- 1.5
  again@leaf.prior@k <- 2
  again@leaf.hyperprior <- dbarts:::newValidated("dbartsFixedHyperprior", k = 2)
  expect_silent(sampler$setModel(again))
  expect_identical(sampler$model@tree.prior@power, 1.5)
  expect_identical(nrow(sampler$run(1L, 2L)$train), n)
})

# xbart's cells vary power, base and k through the bridge and are unaffected
local({
  fit <- xbart(
    x,
    y,
    n.samples = 10L,
    n.test = 4L,
    n.reps = 2L,
    n.trees = c(5L, 10L),
    k = c(1, 2),
    power = 2,
    base = 0.95,
    n.threads = 1L,
    seed = 3L,
    verbose = FALSE
  )
  expect_true(all(is.finite(fit)))
  again <- xbart(
    x,
    y,
    n.samples = 10L,
    n.test = 4L,
    n.reps = 2L,
    n.trees = c(5L, 10L),
    k = c(1, 2),
    power = 2,
    base = 0.95,
    n.threads = 1L,
    seed = 3L,
    verbose = FALSE
  )
  expect_identical(fit, again)
})
