# Inside forest() a prior's name is dbarts's prior only where an argument takes
# a prior (amplitude =); in sd = and the other number arguments a caller's own
# function of that name is called. And every refusal that names a forest names
# it by its label, as the fit's results do.

forest <- dbartsForests$forest

set.seed(41L)
n <- 60L
dat <- data.frame(
  x1 = runif(n),
  x2 = runif(n),
  dose = runif(n),
  z = rbinom(n, 1L, 0.5)
)
dat$y <- dat$x1 + dat$dose * dat$x2 + rnorm(n, sd = 0.3)

# every fit is a literal call: a forests list handed through a wrapper's
# argument is read where it was written, outside the argument that takes it
ctl <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.burn = 0L,
  n.samples = 2L,
  updateState = FALSE,
  seed = 3L
)
sdsOf <- function(fit) vapply(fit$model@forests, function(f) f@sd, 0)
refusal <- function(expr) {
  tryCatch(
    {
      expr
      "no error"
    },
    error = function(e) conditionMessage(e)
  )
}

# ---- B708: a prior's name in a number argument is the caller's ----

calls <- 0L
normal <- function(v) {
  calls <<- calls + 1L
  2 * v
}
fixed <- function(v) 3 * v

fit <- dbarts(
  y ~ x1 + x2,
  dat,
  forests = list(forest(), forest(x2, basis = dose, sd = normal(0.5))),
  n.trees = 5L,
  control = ctl
)
expect_equal(sdsOf(fit)[2L], 1)
expect_equal(calls, 1L)
fit <- dbarts(
  y ~ x1 + x2,
  dat,
  forests = list(forest(), forest(x2, basis = dose, sd = fixed(0.5))),
  n.trees = 5L,
  control = ctl
)
expect_equal(sdsOf(fit)[2L], 1.5)

# the formula's forest() term reads sd the same way
fit <- dbarts(
  y ~ forest(x1) + forest(x2, basis = dose, sd = normal(0.5)),
  dat,
  n.trees = 5L,
  control = ctl
)
expect_equal(sdsOf(fit)[2L], 1)

# $setForest takes sd the same way
sampler <- dbarts(
  y ~ x1 + x2,
  dat,
  forests = list(forest(), forest(x2, basis = dose)),
  n.trees = 5L,
  control = ctl
)
sampler$setForest("dose", sd = normal(0.5))
expect_equal(sampler$getForest("dose")$sd, 1)

# where the argument takes a prior, the name is dbarts's and the caller's
# function is not called
calls <- 0L
fit <- dbarts(
  y ~ x1 + x2,
  dat,
  forests = list(forest(), forest(x2, basis = dose, amplitude = normal())),
  n.trees = 5L,
  control = ctl
)
expect_equal(calls, 0L)
expect_false(fit$model@forests[[2L]]@held.amplitude)

# a basis term that calls a prior's name is refused with a binding too
expect_match(
  refusal(dbarts(
    y ~ x1 + x2,
    dat,
    forests = list(forest(), forest(x2, basis = normal(dose))),
    n.trees = 5L,
    control = ctl
  )),
  "calls normal(), which is not a column",
  fixed = TRUE
)

# a caller's variable of the name that is no function still wins, as it did
rm(fixed)
fixed <- 0.7
fit <- dbarts(
  y ~ x1 + x2,
  dat,
  forests = list(forest(), forest(x2, basis = dose, sd = fixed)),
  n.trees = 5L,
  control = ctl
)
expect_equal(sdsOf(fit)[2L], 0.7)

# with no binding, a prior in sd = is refused as before
rm(normal, fixed, calls)
expect_equal(
  refusal(dbarts(
    y ~ x1 + x2,
    dat,
    forests = list(forest(), forest(x2, basis = dose, sd = normal(0.5))),
    n.trees = 5L,
    control = ctl
  )),
  "forest \"dose\"'s 'sd' must be a number, not a dbartsNormalPrior"
)
expect_equal(
  refusal(dbarts(
    y ~ x1 + x2,
    dat,
    forests = list(forest(), forest(x2, basis = dose, sd = fixed(0.5))),
    n.trees = 5L,
    control = ctl
  )),
  "forest \"dose\"'s 'sd' must be a number, not a dbartsFixedPrior"
)

# ---- B713: a refusal names a forest by its label ----

mustBe <- "'n.perturb.cuts' must be a single positive whole number or Inf"
one <- function(label) paste0("forest ", label, "'s ", mustBe)

# unnamed, the second forest carrying a basis: labeled "dose", at creation and
# at $setForest by position and by name
expect_equal(
  refusal(dbarts(
    y ~ x1 + x2,
    dat,
    forests = list(forest(), forest(basis = ~dose, n.perturb.cuts = 0)),
    n.trees = 5L,
    control = ctl
  )),
  one("\"dose\"")
)
sampler <- dbarts(
  y ~ x1 + x2,
  dat,
  forests = list(forest(), forest(basis = ~dose)),
  n.trees = 5L,
  control = ctl
)
expect_equal(
  refusal(sampler$setForest(2, n.perturb.cuts = 0)),
  one("\"dose\"")
)
expect_equal(
  refusal(sampler$setForest("dose", n.perturb.cuts = 0)),
  one("\"dose\"")
)

# named mu and tau: tau at creation and at both doors
expect_equal(
  refusal(dbarts(
    y ~ x1 + x2,
    dat,
    forests = list(
      mu = forest(),
      tau = forest(basis = ~dose, n.perturb.cuts = 0)
    ),
    n.trees = 5L,
    control = ctl
  )),
  one("\"tau\"")
)
sampler <- dbarts(
  y ~ x1 + x2,
  dat,
  forests = list(mu = forest(), tau = forest(basis = ~dose)),
  n.trees = 5L,
  control = ctl
)
expect_equal(
  refusal(sampler$setForest(2, n.perturb.cuts = 0)),
  one("\"tau\"")
)
expect_equal(
  refusal(sampler$setForest("tau", n.perturb.cuts = 0)),
  one("\"tau\"")
)

# named mu and dose, the name being the basis's text: dose at both doors
sampler <- dbarts(
  y ~ x1 + x2,
  dat,
  forests = list(mu = forest(), dose = forest(basis = ~dose)),
  n.trees = 5L,
  control = ctl
)
expect_equal(
  refusal(sampler$setForest(2, n.perturb.cuts = 0)),
  one("\"dose\"")
)
expect_equal(
  refusal(sampler$setForest("dose", n.perturb.cuts = 0)),
  one("\"dose\"")
)

# a forest with no name and no basis text stays "forest 2"
handed <- dat$dose
expect_equal(
  refusal(dbarts(
    y ~ x1 + x2,
    dat,
    forests = list(
      forest(),
      do.call(forest, list(basis = handed, n.perturb.cuts = 0))
    ),
    n.trees = 5L,
    control = ctl
  )),
  one(2L)
)
sampler <- dbarts(
  y ~ x1 + x2,
  dat,
  forests = list(forest(), do.call(forest, list(basis = handed))),
  n.trees = 5L,
  control = ctl
)
expect_equal(refusal(sampler$setForest(2, n.perturb.cuts = 0)), one(2L))

# ---- B708, the rest: other number arguments, an alias, a nested call ----

invchi <- function(v) v + 1
normal <- function(v) 2 * v
fit <- dbarts(
  y ~ x1 + x2,
  dat,
  forests = list(forest(), forest(x2, basis = dose, n.trees = invchi(4))),
  n.trees = 5L,
  control = ctl
)
expect_equal(fit$model@forests[[2L]]@n.trees, 5L)

# a call through an alias of forest() reads the same way
fc <- dbartsForests$forest
fit <- dbarts(
  y ~ x1 + x2,
  dat,
  forests = list(fc(), fc(x2, basis = dose, sd = normal(0.5))),
  n.trees = 5L,
  control = ctl
)
expect_equal(sdsOf(fit)[2L], 1)

# and through the exported list and the namespace's own spelling, so that
# forest(), an alias and the two spellings agree
fit <- dbarts(
  y ~ x1 + x2,
  dat,
  forests = list(
    forest(),
    dbartsForests$forest(x2, basis = dose, sd = normal(0.5))
  ),
  n.trees = 5L,
  control = ctl
)
expect_equal(sdsOf(fit)[2L], 1)
fit <- dbarts(
  y ~ x1 + x2,
  dat,
  forests = list(forest(), dbarts:::forest(x2, basis = dose, sd = normal(0.5))),
  n.trees = 5L,
  control = ctl
)
expect_equal(sdsOf(fit)[2L], 1)

# and a prior's name nested in an expression
fit <- dbarts(
  y ~ x1 + x2,
  dat,
  forests = list(forest(), forest(x2, basis = dose, sd = 1 + normal(0.5))),
  n.trees = 5L,
  control = ctl
)
expect_equal(sdsOf(fit)[2L], 2)
rm(invchi, normal, fc)

# ---- B713, getTrees names a labeled forest by its label ----

ctlKeep <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.burn = 0L,
  n.samples = 2L,
  updateState = FALSE,
  seed = 3L,
  keepTrees = TRUE
)
sampler <- dbarts(
  y ~ x1 + x2,
  dat,
  forests = list(mu = forest(), dose = forest(basis = ~dose)),
  n.trees = 5L,
  control = ctlKeep
)
expect_equal(
  refusal(sampler$getTrees(forest = "dose", treeNums = 99L)),
  "'treeNums' must be in [1, 50] for forest \"dose\""
)
expect_equal(
  refusal(sampler$getTrees(forest = 2L, treeNums = 99L)),
  "'treeNums' must be in [1, 50] for forest \"dose\""
)
