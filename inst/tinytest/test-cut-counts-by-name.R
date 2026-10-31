# A per-column n.cuts reads as split.probs reads at every door that takes it: by
# column or term name with .default for the rest, or unnamed, one value or one
# per column. Anything else is refused.

set.seed(71)
n <- 60L
frame <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n))
frame$y <- frame$x1 + rnorm(n, sd = 0.3)
frame$f <- factor(sample(c("a", "b", "c"), n, TRUE))
X <- as.matrix(frame[c("x1", "x2", "x3")])
y <- frame$y
run <- function(...) {
  dbartsControl(n.samples = 2L, n.burn = 0L, n.chains = 1L, n.threads = 1L, ...)
}
written <- function(...) suppressWarnings(run(...))
form <- y ~ x1 + x2 + x3
doors <- list(
  dbarts = function(v) dbarts(form, frame, n.cuts = v, control = run())$data@n.cuts,
  matrix = function(v) dbarts(X, y, n.cuts = v, control = run())$data@n.cuts,
  data = function(v) dbartsData(form, frame, n.cuts = v)@n.cuts,
  control = function(v) {
    suppressWarnings(dbarts(form, frame, control = run(n.cuts = v)))$data@n.cuts
  },
  spec = function(v) {
    suppressWarnings(dbartsSpec(dbartsData(form, frame), run(n.cuts = v)))$data@n.cuts
  },
  bart = function(v) {
    bart(
      form,
      frame,
      n.cuts = v,
      n.samples = 2L,
      n.burn = 0L,
      n.chains = 1L,
      n.threads = 1L,
      verbose = FALSE,
      samplerOnly = TRUE
    )$data@n.cuts
  }
)
for (door in names(doors)) {
  read <- doors[[door]]
  expect_identical(read(c(x3 = 5L, x1 = 25L, x2 = 10L)), c(25L, 10L, 5L), info = door)
  expect_identical(read(c(x3 = 5, x1 = 25, x2 = 10)), c(25L, 10L, 5L), info = door)
  expect_identical(read(c(x2 = 10L, .default = 7L)), c(7L, 10L, 7L), info = door)
  expect_identical(read(c(.default = 7L)), rep(7L, 3L), info = door)
  expect_identical(read(7L), rep(7L, 3L), info = door)
  expect_identical(read(c(5L, 6L, 7L)), c(5L, 6L, 7L), info = door)
  expect_error(read(c(x2 = 10L)), "missing values for columns 'x1', 'x3'", info = door)
  expect_error(
    read(c(zz = 10L, x1 = 25L, x2 = 5L, x3 = 1L)),
    "unrecognized variable name 'zz'",
    info = door
  )
  expect_error(read(c(25L, 10L)), "'n.cuts' has 2 values", info = door)
  expect_error(read(1:4), "'n.cuts' has 4 values", info = door)
  expect_error(
    read(c(.default = 1L, .default = 2L)),
    "default specified multiple times",
    info = door
  )
  expect_error(read(c(x1 = 5L, 6L, 7L)), "unrecognized variable name ''", info = door)
}

# a single named value is named, taken where it names the only column
expect_identical(
  dbarts(y ~ x2, frame, n.cuts = c(x2 = 10L), control = run())$data@n.cuts,
  10L
)
# a term names each of its columns
expect_identical(
  dbarts(
    y ~ x1 + f,
    frame,
    factors = "indicators",
    n.cuts = c(f = 3L, x1 = 9L),
    control = run()
  )$data@n.cuts,
  c(9L, 3L, 3L, 3L)
)
# a matrix with no names has no columns to name
unnamedX <- unname(X)
expect_error(
  dbarts(unnamedX, y, n.cuts = c(x1 = 5L, .default = 6L), control = run()),
  "model matrix has no column names"
)
expect_identical(
  dbarts(unnamedX, y, n.cuts = c(.default = 6L), control = run())$data@n.cuts,
  rep(6L, 3L)
)

# the count reaches the column it names
sampler <- dbarts(form, frame, n.cuts = c(x3 = 5L, x1 = 25L, x2 = 10L), control = run())
invisible(sampler$run(0L, 2L))
expect_identical(lengths(sampler$data@cutPoints), c(25L, 10L, 5L))

# the other doors that read one: bartBT names its own argument
expect_error(
  suppressWarnings(bartBT(X, y, numcut = c(zz = 5L, .default = 6L), ndpost = 2L, nskip = 0L, verbose = FALSE)),
  "cannot assign 'numcut': unrecognized variable name 'zz'"
)
expect_error(
  suppressWarnings(rbart_vi(
    form,
    frame,
    group.by = rep(1:3, 20L),
    n.cuts = c(zz = 5L, .default = 6L),
    n.trees = 5L,
    n.thin = 1L,
    n.samples = 2L,
    n.burn = 0L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  )),
  "unrecognized variable name 'zz'"
)
xb <- function(...) {
  xbart(
    X,
    y,
    n.trees = 5L,
    ...,
    n.samples = 2L,
    method = "k-fold",
    n.test = 2L,
    n.reps = 1L,
    n.burn = c(2L, 1L),
    n.threads = 1L,
    seed = 3L
  )
}
expect_identical(xb(n.cuts = c(x3 = 5, x1 = 25, x2 = 10)), xb(n.cuts = c(25L, 10L, 5L)))
expect_error(xb(n.cuts = c(zz = 5L, .default = 6L)), "unrecognized variable name 'zz'")
expect_error(xb(n.cuts = c(25L, 10L)), "'n.cuts' has 2 values")

# a count written on a control is read the same beside a sampler's and a door's:
# agreeing by name is agreeing, and differing by position is differing
sampler <- dbarts(form, frame, n.cuts = c(25L, 10L, 5L), control = run())
expect_silent(sampler$setControl(written(n.cuts = c(x3 = 5L, x1 = 25L, x2 = 10L))))
expect_error(
  sampler$setControl(written(n.cuts = c(5L, 25L, 10L))),
  "changing 'n.cuts' is not available on an existing sampler"
)
expect_error(
  sampler$setControl(written(n.cuts = c(zz = 5L, .default = 6L))),
  "unrecognized variable name 'zz'"
)
notUsedAtDirect <- function(...) {
  told <- character()
  withCallingHandlers(
    new("dbartsSampler", written(...), sampler$model, sampler$data),
    dbartsControlNotUsedWarning = function(w) {
      told <<- c(told, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    warning = function(w) invokeRestart("muffleWarning")
  )
  told
}
expect_identical(notUsedAtDirect(n.cuts = c(x3 = 5L, x1 = 25L, x2 = 10L)), character())
expect_identical(length(notUsedAtDirect(n.cuts = c(5L, 25L, 10L))), 1L)
expect_error(notUsedAtDirect(n.cuts = c(25L, 10L)), "'n.cuts' has 2 values")

# split.probs reads as it did, a single value losing its name
probs <- function(v) {
  sampler <- dbarts(form, frame, split.probs = v, control = run())
  sampler$model@tree.prior@splitProbabilities
}
expect_equal(
  unname(probs(c(x3 = 0.7, x1 = 0.2, x2 = 0.1))),
  c(0.2, 0.1, 0.7)
)
expect_identical(length(probs(c(x2 = 0.9))), 0L)
expect_error(
  probs(c(x1 = 0.5, x2 = 0.5)),
  "cannot assign split probabilities: missing values for columns 'x3'"
)
expect_error(
  probs(c(zz = 0.5, x1 = 0.5, x2 = 0.5, x3 = 0.5)),
  "cannot assign split probabilities: unrecognized variable name 'zz'"
)
