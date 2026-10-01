# Every warning the package raises carries a class that inherits from
# dbartsWarning; the retired spellings' once-per-session warnings are
# dbartsDeprecatedWarning, the substitution warnings dbartsFallbackWarning.

resetOnce <- function() {
  env <- dbarts:::onceWarnState
  rm(list = ls(env, all.names = TRUE), envir = env)
}
expectClassed <- function(expr, class) {
  resetOnce()
  w <- tryCatch(
    withCallingHandlers(expr, message = function(m) {
      invokeRestart("muffleMessage")
    }),
    warning = function(w) w
  )
  expect_true(inherits(w, class))
  expect_true(inherits(w, "dbartsWarning"))
}

set.seed(1)
x <- matrix(rnorm(60), 30L, 2L, dimnames = list(NULL, c("a", "b")))
y <- rnorm(30)
d <- data.frame(y, x)
quick <- list(
  n.samples = 4L,
  n.burn = 2L,
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 3L,
  verbose = FALSE
)
control <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 3L,
  n.samples = 2L,
  updateState = FALSE
)

expectClassed(
  do.call(dbarts::bart2, c(list(y ~ ., d), quick)),
  "dbartsDeprecatedWarning"
)
expectClassed(
  dbarts::bart(
    x.train = x,
    y.train = y,
    ndpost = 4L,
    nskip = 2L,
    nchain = 1L,
    nthread = 1L,
    ntree = 3L,
    verbose = FALSE
  ),
  "dbartsDeprecatedWarning"
)
expectClassed(
  do.call(dbarts::bart, c(list(y ~ ., d, sigdf = 3), quick)),
  "dbartsDeprecatedWarning"
)
expectClassed(
  do.call(dbarts::bart, c(list(y ~ ., d, rngSeed = 1L), quick)),
  "dbartsDeprecatedWarning"
)
expectClassed(
  dbarts::dbarts(y ~ ., d, sigma = 1, control = control),
  "dbartsDeprecatedWarning"
)
expectClassed(
  do.call(dbarts::bart, c(list(y ~ ., d, seed = NA_integer_), quick)),
  "dbartsDeprecatedWarning"
)
expectClassed(
  dbarts::dbarts(
    y ~ .,
    d,
    node.prior = dbarts::dbartsPriors$normal(3),
    control = control
  ),
  "dbartsDeprecatedWarning"
)
expectClassed(
  dbarts::dbartsPriors$chi(degreesOfFreedom = 2),
  "dbartsDeprecatedWarning"
)

sampler <- dbarts::dbarts(y ~ ., d, control = control)
expectClassed(sampler$startThreads(), "dbartsDeprecatedWarning")
expectClassed(sampler$run(1L, 1L, n.threads = 1L), "dbartsDeprecatedWarning")
expectClassed(
  sampler$sampleNodeParametersFromPrior(updateState = FALSE),
  "dbartsDeprecatedWarning"
)
expectClassed(
  dbarts::rbart_vi(
    y ~ a + b,
    d,
    group.by = factor(rep(1:3, 10L)),
    n.samples = 4L,
    n.burn = 2L,
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 3L,
    verbose = FALSE
  ),
  "dbartsDeprecatedWarning"
)

sf <- dbarts::sparseFactor(c("u", "v", "u", "v"))
expectClassed(sf[1L] <- "w", "dbartsFallbackWarning")
expectClassed(sf < sf, "dbartsFallbackWarning")
