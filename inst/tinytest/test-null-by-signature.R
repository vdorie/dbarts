# A NULL written at creation is taken where the argument's default in the
# signature is NULL, and there it is the name left out; anywhere else it is
# refused, saying to leave the name out. It is judged by the value, so a
# wrapper's argument is read by what it holds.

set.seed(5L)
n <- 50L
d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
d$y <- d$x1 + rnorm(n)
fit <- function(..., n.trees = 5L) {
  bart(
    ...,
    n.samples = 5L,
    n.burn = 5L,
    n.chains = 1L,
    n.threads = 1L,
    n.trees = n.trees,
    seed = 1L,
    verbose = FALSE
  )$yhat.train
}
refusal <- function(name) {
  paste0("'", name, "' cannot be NULL; leave it out for the default")
}

# a NULL where the default is NULL is the name left out
base <- fit(y ~ x1 + x2, d)
expect_equal(fit(y ~ x1 + x2, d, k = NULL), base)
expect_equal(fit(y ~ x1 + x2, d, leaf.prior = NULL), base)
expect_equal(fit(y ~ x1 + x2, d, split.probs = NULL), base)
sized <- fit(y ~ forest(x1 + x2, sd = 2), d)
expect_equal(fit(y ~ forest(x1 + x2, sd = 2), d, k = NULL), sized)
expect_equal(fit(y ~ forest(x1 + x2, sd = 2), d, leaf.prior = NULL), sized)
expect_equal(
  fit(y ~ forest(x1 + x2, sd = NULL), d, k = 3),
  fit(y ~ x1 + x2, d, k = 3)
)
powered <- fit(y ~ forest(x1 + x2), d, power = 3)
expect_equal(fit(y ~ forest(x1 + x2, power = NULL), d, power = 3), powered)
expect_equal(
  fit(y ~ forest(x1 + x2, n.trees = NULL), d, n.trees = 10L),
  fit(y ~ forest(x1 + x2), d, n.trees = 10L)
)
# a real value is stated
expect_error(fit(y ~ forest(x1 + x2, sd = 2), d, k = 3), "size twice")
expect_error(
  fit(y ~ forest(x1 + x2, sd = 2), d, leaf.prior = normal(3)),
  "size twice"
)
# the default's own value is stated too
expect_error(fit(y ~ forest(x1 + x2, power = 2), d, power = 2), "given")

# a NULL where the default is not NULL is refused at every door, whatever
# else the call holds
nulled <- function(name) setNames(list(NULL), name)
refused <- function(door, name, ...) {
  expect_error(do.call(door, c(list(...), nulled(name))), refusal(name))
}
quiet <- list(n.chains = 1L, n.threads = 1L, verbose = FALSE)
for (name in c("n.trees", "power", "base")) {
  refused(bart, name, y ~ x1 + x2, d, n.samples = 5L, n.burn = 5L)
  refused(bart, name, y ~ forest(x1 + x2, power = 3), d, n.samples = 5L)
  refused(xbart, name, y ~ x1 + x2, d, n.samples = 5L, n.burn = c(5L, 5L))
  refused(rbart_vi, name, y ~ x1 + x2, d, group.by = rep(1:5, 10L))
}
expect_error(dbartsControl(n.trees = NULL), refusal("n.trees"))
for (name in c("power", "base", "leaf.prior")) {
  refused(dbarts, name, y ~ x1 + x2, d)
  refused(dbartsSpec, name, dbartsData(y ~ x1 + x2, d))
}
refused(rbart_vi, "k", y ~ x1 + x2, d, group.by = rep(1:5, 10L))
for (name in c("ntree", "power", "base", "k")) {
  refused(dbarts:::bartBT, name, d[, 1:2], d$y)
}
# a forest's own arguments default to NULL, so a NULL there is left out
expect_equal(
  fit(y ~ forest(x1 + x2, power = NULL, base = NULL), d),
  fit(y ~ x1 + x2, d)
)
expect_equal(
  fit(y ~ x1 + x2, d, variance = ~ varianceForest(x1, power = NULL)),
  fit(y ~ x1 + x2, d, variance = ~ varianceForest(x1))
)

# judged by the value: a wrapper's variable, or one passed on by '...'
wrap <- function(formula, data, k = NULL, ...) {
  fit(formula, data, k = k, ...)
}
expect_equal(wrap(y ~ forest(x1 + x2, sd = 2), d), sized)
expect_equal(wrap(y ~ forest(x1 + x2, sd = 2), d, k = NULL), sized)
hold <- function(formula, data, power = NULL) {
  fit(formula, data, power = power)
}
expect_error(hold(y ~ x1 + x2, d), refusal("power"))
expect_error(hold(y ~ forest(x1 + x2, power = 3), d), refusal("power"))
expect_equal(hold(y ~ x1 + x2, d, power = 2), base)
through <- function(formula, data, ...) fit(formula, data, ...)
expect_equal(through(y ~ forest(x1 + x2, sd = 2), d, k = NULL), sized)
expect_error(through(y ~ x1 + x2, d, power = NULL), refusal("power"))
p <- NULL
expect_equal(fit(y ~ forest(x1 + x2, power = p), d, power = 3), powered)

# a written call is evaluated once, not once to test it for NULL and again to
# use it: k = runif() draws what it draws without the test
set.seed(7L)
drawn <- fit(y ~ x1 + x2, d, k = runif(1L, 1, 3))
set.seed(7L)
held <- runif(1L, 1, 3)
expect_equal(
  drawn,
  fit(y ~ x1 + x2, d, k = held),
  info = "k = runif() draws once"
)
ticks <- 0L
tick <- function(value) {
  ticks <<- ticks + 1L
  value
}
once <- function(info, ...) {
  ticks <<- 0L
  suppressWarnings(suppressMessages(...))
  expect_identical(ticks, 1L, info = info)
}
once("bart k", fit(y ~ x1 + x2, d, k = tick(2)))
once(
  "xbart k",
  xbart(
    y ~ x1 + x2,
    d,
    n.samples = 5L,
    n.burn = c(5L, 5L),
    n.reps = 1L,
    n.threads = 1L,
    n.test = 5,
    verbose = FALSE,
    k = tick(2)
  )
)
once(
  "xbart sd",
  xbart(
    y ~ x1 + x2,
    d,
    n.samples = 5L,
    n.burn = c(5L, 5L),
    n.reps = 1L,
    n.threads = 1L,
    n.test = 5,
    verbose = FALSE,
    sd = tick(1)
  )
)
once("dbarts power", dbarts(y ~ x1 + x2, d, power = tick(2)))
once("dbarts base", dbarts(y ~ x1 + x2, d, base = tick(0.95)))
once(
  "dbartsSpec power",
  dbartsSpec(dbartsData(y ~ x1 + x2, d), power = tick(2))
)
once(
  "bartBT k",
  dbarts:::bartBT(
    d[, 1:2],
    d$y,
    ndpost = 5L,
    nskip = 5L,
    ntree = 5L,
    verbose = FALSE,
    k = tick(2)
  )
)
once(
  "rbart_vi k",
  rbart_vi(
    y ~ x1 + x2,
    d,
    group.by = rep(1:5, 10L),
    n.samples = 5L,
    n.burn = 5L,
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 5L,
    verbose = FALSE,
    k = tick(2)
  )
)

# settings that default to a value refuse a NULL as well
expect_error(dbartsControl(proposal.probs = NULL), refusal("proposal.probs"))
expect_error(fit(y ~ x1 + x2, d, samplerOnly = NULL), refusal("samplerOnly"))
expect_error(
  xbart(
    y ~ x1 + x2,
    d,
    n.samples = 5L,
    n.burn = c(5L, 5L),
    n.reps = 1L,
    parallel = NULL
  ),
  refusal("parallel")
)

# node.prior defaults to NULL in the signature, as tree.prior does: a written
# NULL is the name left out
expect_silent(nodeNull <- dbarts(y ~ x1 + x2, d, node.prior = NULL))
expect_silent(dbarts(y ~ x1 + x2, d, tree.prior = NULL))
expect_identical(
  nodeNull$model@leaf.prior,
  dbarts(y ~ x1 + x2, d)$model@leaf.prior
)

# a held formula selects the forest's predictors as the bare spelling does
f <- ~ x1 + x2
expect_equal(
  fit(y ~ forest(f), d),
  fit(y ~ forest(x1 + x2), d),
  info = "held formula door"
)
