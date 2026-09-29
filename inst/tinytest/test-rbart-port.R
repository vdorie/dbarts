## rbart_vi as a deprecated one-release feature: the warning, each fix to the
## 0.9-x loop, the refusals, recovery, and the call shapes bartCause makes.

fitRbart <- function(...) {
  args <- list(
    n.samples = 5L,
    n.burn = 0L,
    n.thin = 1L,
    n.chains = 1L,
    n.trees = 5L,
    n.threads = 1L,
    verbose = FALSE
  )
  args[names(list(...))] <- list(...)
  suppressWarnings(do.call(dbarts::rbart_vi, args))
}

simulateGrouped <- function(n, n.g, tau, sigma = 1, seed = 1L) {
  set.seed(seed)
  x <- matrix(rnorm(n * 3L), n, 3L)
  g <- factor(sample(n.g, n, replace = TRUE), levels = seq_len(n.g))
  b <- rnorm(n.g, 0, tau)
  eta <- 2 * x[, 1L] + x[, 2L]^2 + b[g]
  list(x = x, g = g, b = b, eta = eta, sigma = sigma, n = n)
}

# --- the deprecation warning: once per session, from rbart_vi only ---
onceState <- dbarts:::onceWarnState
onceState[["tombstone.rbart_vi"]] <- NULL
sim <- simulateGrouped(60L, 4L, 1)
y <- sim$eta + rnorm(sim$n)
x <- sim$x
g <- sim$g

countWarnings <- function(expr) {
  messages <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, messages = messages)
}
rbartCall <- function() {
  dbarts::rbart_vi(
    y ~ x,
    group.by = g,
    n.samples = 5L,
    n.burn = 0L,
    n.thin = 1L,
    n.chains = 1L,
    n.trees = 5L,
    n.threads = 1L,
    verbose = FALSE
  )
}

first <- countWarnings(rbartCall())
expect_equal(length(first$messages), 1L)
expect_true(grepl("stan4bart", first$messages, fixed = TRUE))
expect_true(grepl("1.1-0", first$messages, fixed = TRUE))
second <- countWarnings(rbartCall())
expect_equal(length(second$messages), 0L)

# the methods never warn, even in a session that has not yet
onceState[["tombstone.rbart_vi"]] <- NULL
fit <- first$value
methodWarnings <- countWarnings({
  predict(fit, x, g)
  dbarts::extract(fit)
  fitted(fit)
  residuals(fit)
  capture.output(print(fit))
  pdf(NULL)
  plot(fit)
  dev.off()
})$messages
expect_equal(sum(grepl("deprecated", methodWarnings)), 0L)
onceState[["tombstone.rbart_vi"]] <- TRUE
rm(first, second, fit, methodWarnings, countWarnings, rbartCall, y, x, g, sim)

# --- D1: the first intercept draw starts from a sweep's fit ---
# On a skewed response the midpoint of the range is far from its mean, and
# starting from the midpoint put that gap into the intercepts for good
# (posterior mean tau near 100 for a truth of 1).
set.seed(1L)
n <- 1000L
x <- matrix(rnorm(n * 5L), n, 5L)
f <- 10 *
  sin(pi * x[, 1L] * x[, 2L]) +
  20 * (x[, 3L] - 0.5)^2 +
  10 * x[, 4L] +
  5 * x[, 5L]
g <- factor(sample(10L, n, TRUE))
b <- rnorm(10L)
y <- f + b[g] + rnorm(n)
fit <- fitRbart(
  y ~ x,
  group.by = g,
  n.samples = 500L,
  n.burn = 500L
)
expect_true(mean(fit$tau) < 5)
rm(n, x, f, g, b, y, fit)

# --- D2: a group.by that is not a column of data is not the first column ---
sim <- simulateGrouped(60L, 4L, 1)
df <- data.frame(x_1 = sim$x[, 1L], x_2 = sim$x[, 2L])
df$y <- sim$eta + rnorm(sim$n)
g <- sim$g
# the calls are written out: a symbol has to reach rbart_vi as a symbol
rbartSymbols <- function(...) {
  suppressWarnings(dbarts::rbart_vi(
    y ~ x_1 + x_2,
    df,
    n.samples = 5L,
    n.burn = 0L,
    n.thin = 1L,
    n.chains = 1L,
    n.trees = 5L,
    n.threads = 1L,
    verbose = FALSE,
    ...
  ))
}
fit <- rbartSymbols(group.by = g)
expect_equal(nlevels(fit$group.by), 4L)
expect_error(rbartSymbols(group.by = not_a_symbol), "'group.by' not found")
g.test <- g[1:10]
fit <- rbartSymbols(test = df[1:10, ], group.by = g, group.by.test = g.test)
expect_equal(as.character(fit$group.by.test), as.character(g.test))
expect_error(
  rbartSymbols(
    test = df[1:10, ],
    group.by = g,
    group.by.test = not_a_symbol
  ),
  "'group.by.test' not found"
)
rm(sim, df, g, g.test, fit, rbartSymbols)

# --- D3: a new level with several chains, combined draws ---
sim <- simulateGrouped(60L, 4L, 1)
x <- sim$x
y <- sim$eta + rnorm(sim$n)
g <- sim$g
g.new <- factor(c("1", "2", "9", "9"), levels = c("1", "2", "9"))
x.new <- x[1:4, ]
for (fitCombined in c(FALSE, TRUE)) {
  fit <- fitRbart(
    y ~ x,
    group.by = g,
    n.chains = 2L,
    n.samples = 6L,
    combineChains = fitCombined
  )
  set.seed(3L)
  combined <- suppressWarnings(predict(fit, x.new, g.new, combineChains = TRUE))
  expect_equal(dim(combined), c(12L, 4L))
  expect_true(all(is.finite(combined)))
  set.seed(3L)
  split <- suppressWarnings(predict(fit, x.new, g.new, combineChains = FALSE))
  expect_equal(dim(split), c(2L, 6L, 4L))
  # the two layouts hold the same measured-level draws, chain-major combined
  expect_equal(combined[1:6, 1:2], split[1, , 1:2])
  expect_equal(combined[7:12, 1:2], split[2, , 1:2])
}
# the new level's draws are scaled by the posterior tau of their own draw
set.seed(4L)
fit <- fitRbart(
  y ~ x,
  group.by = g,
  n.chains = 2L,
  n.samples = 6L,
  combineChains = FALSE
)
fit$tau[1L, ] <- 1e-8
fit$tau[2L, ] <- 1e3
ranefNew <- suppressWarnings(
  predict(fit, x.new, g.new, type = "ranef", combineChains = TRUE)
)
expect_true(all(abs(ranefNew[1:6, "9"]) < 1e-3))
expect_true(all(abs(ranefNew[7:12, "9"]) > 1e-3))
rm(sim, x, y, g, g.new, x.new, fit, combined, split, ranefNew, fitCombined)

# --- D4 and PSOCK: predict and trees from a fit that came back from a worker ---
sim <- simulateGrouped(60L, 4L, 1)
x <- sim$x
y <- sim$eta + rnorm(sim$n)
g <- sim$g
fit <- fitRbart(
  y ~ x,
  group.by = g,
  n.chains = 2L,
  n.threads = 2L,
  n.samples = 6L
)
expect_equal(
  dim(predict(fit, x, g, type = "bart", combineChains = FALSE)),
  c(2L, 6L, 60L)
)
trees <- dbarts::extract(fit, type = "trees")
expect_true(is.data.frame(trees))
expect_true(all(c("sample", "chain", "tree") %in% names(trees)))
expect_equal(sort(unique(trees$chain)), 1:2)

# a serial fit that is saved and read back predicts what it did before
fit <- fitRbart(y ~ x, group.by = g, n.chains = 2L, n.samples = 6L)
before <- predict(fit, x, g)
path <- tempfile(fileext = ".rds")
saveRDS(fit, path)
after <- predict(readRDS(path), x, g)
unlink(path)
expect_equal(after, before)
rm(sim, x, y, g, fit, trees, before, after, path)

# --- D5: the intercepts follow the weights as precisions ---
# Rescaling every weight rescales sigma with it, so the posterior of tau and
# of the intercepts must not move; a group's intercept was drawn as if the
# weights were all 1, which left tau moving with their scale.
n <- 300L
sim <- simulateGrouped(n, 10L, 1, seed = 5L)
x <- sim$x
g <- sim$g
y <- sim$eta + rnorm(n)
meanTau <- function(w) {
  set.seed(6L)
  fit <- fitRbart(
    y ~ x,
    group.by = g,
    weights = w,
    n.samples = 400L,
    n.burn = 200L
  )
  mean(fit$tau)
}
tauOne <- meanTau(rep(1, n))
tauSmall <- meanTau(rep(1 / n, n))
expect_equal(tauSmall, tauOne, tolerance = 0.05)
rm(n, sim, x, g, y, meanTau, tauOne, tauSmall)

# --- D6: test rows do not inherit the training intercepts ---
n <- 80L
sim <- simulateGrouped(n, 5L, 1, seed = 7L)
x <- sim$x
g <- sim$g
y <- sim$eta + rnorm(n)
fit <- fitRbart(
  y ~ x,
  group.by = g,
  test = x,
  group.by.test = g,
  offset = rep(0.5, n),
  n.samples = 6L
)
expect_equal(
  dbarts::extract(fit, type = "bart", sample = "test"),
  dbarts::extract(fit, type = "bart", sample = "train")
)
rm(n, sim, x, g, y, fit)

# --- D7: a seeded fit leaves the caller's random stream alone ---
sim <- simulateGrouped(40L, 4L, 1)
x <- sim$x
g <- sim$g
y <- sim$eta + rnorm(sim$n)
for (nThreads in 1:2) {
  set.seed(11L)
  before <- .Random.seed
  fitRbart(
    y ~ x,
    group.by = g,
    seed = 3L,
    n.chains = 2L,
    n.threads = nThreads
  )
  expect_identical(.Random.seed, before)

  rm(".Random.seed", envir = globalenv())
  fitRbart(
    y ~ x,
    group.by = g,
    seed = 3L,
    n.chains = 2L,
    n.threads = nThreads
  )
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
}
set.seed(12L)
rm(sim, x, g, y, nThreads, before)

# --- M9 and the refusals before any chain starts ---
n <- 60L
sim <- simulateGrouped(n, 4L, 1)
x <- sim$x
g <- sim$g
yBinary <- as.integer(sim$eta > median(sim$eta))
w <- runif(n, 0.5, 2)
for (nThreads in 1:2) {
  seen <- character()
  result <- withCallingHandlers(
    tryCatch(
      dbarts::rbart_vi(
        yBinary ~ x,
        group.by = g,
        weights = w,
        n.chains = 2L,
        n.threads = nThreads,
        n.samples = 5L,
        n.burn = 0L,
        n.thin = 1L,
        n.trees = 5L,
        verbose = FALSE
      ),
      error = function(e) e
    ),
    warning = function(cond) {
      seen <<- c(seen, conditionMessage(cond))
      invokeRestart("muffleWarning")
    }
  )
  expect_inherits(result, "error")
  expect_true(grepl("rbart_vi", conditionMessage(result), fixed = TRUE))
  expect_true(grepl("0 and 1 weights", conditionMessage(result), fixed = TRUE))
  expect_false(any(grepl("defaulting to single", seen, fixed = TRUE)))
}
# refusals a worker would otherwise swallow surface once, up front
seen <- character()
for (nThreads in 1:2) {
  three <- factor(rep(c("a", "b", "c"), length.out = n))
  for (refused in list(
    quote(dbarts::rbart_vi(three ~ x, group.by = g, n.threads = nThreads)),
    quote(dbarts::rbart_vi(
      sim$eta ~ x,
      group.by = g,
      k = -1,
      n.threads = nThreads
    ))
  )) {
    refused$n.chains <- 2L
    refused$n.samples <- 5L
    refused$n.burn <- 0L
    refused$n.thin <- 1L
    refused$n.trees <- 5L
    refused$verbose <- FALSE
    result <- withCallingHandlers(
      tryCatch(eval(refused), error = function(e) e),
      warning = function(cond) {
        seen <<- c(seen, conditionMessage(cond))
        invokeRestart("muffleWarning")
      }
    )
    expect_inherits(result, "error")
  }
}
expect_false(any(grepl("defaulting to single", seen, fixed = TRUE)))
# 0 and 1 weights are fine
expect_inherits(
  fitRbart(yBinary ~ x, group.by = g, weights = rep(c(0, 1), n / 2L)),
  "rbart"
)
# a survival response would resolve to aft
if (requireNamespace("survival", quietly = TRUE)) {
  status <- rep(c(0L, 1L), n / 2L)
  expect_error(
    dbarts::rbart_vi(
      survival::Surv(exp(sim$eta), status) ~ x,
      group.by = g,
      n.chains = 1L,
      n.threads = 1L,
      n.samples = 5L,
      n.burn = 0L,
      n.thin = 1L,
      verbose = FALSE
    ),
    "continuous or binary"
  )
  rm(status)
}

# keepFits = FALSE is overridden: the loop reads the fits every sweep
expect_inherits(
  fitRbart(sim$eta ~ x, group.by = g, keepFits = FALSE),
  "rbart"
)

# a group whose rows all have weight 0 informs nothing: its intercepts are
# draws from the prior, as spread as tau
set.seed(41L)
nz <- 400L
dz <- data.frame(x = rnorm(nz), g = factor(sample(6L, nz, TRUE)))
bz <- rnorm(6L)
dz$z <- rbinom(nz, 1L, pnorm(dz$x + bz[as.integer(dz$g)]))
fitZero <- fitRbart(
  z ~ x,
  dz,
  group.by = dz$g,
  weights = as.numeric(dz$g != "1"),
  n.samples = 500L,
  n.burn = 200L
)
ratio <- sd(fitZero$ranef[, "1"]) / mean(fitZero$tau)
expect_true(ratio > 0.6 && ratio < 1.6)
rm(nz, dz, bz, fitZero, ratio)

# a fit saved by 0.9-x has samplers this version cannot re-create
fit <- fitRbart(sim$eta ~ x, group.by = g)
fit$fit <- list(new.env())
expect_error(predict(fit, x, g), "saved by dbarts 0.9-x")
expect_error(dbarts::extract(fit, type = "trees"), "saved by dbarts 0.9-x")
rm(n, sim, x, g, yBinary, w, nThreads, seen, result, three, refused, fit)

# --- recovery on a simulated design ---
recover <- function(binary) {
  n <- if (binary) 2000L else 1000L
  sim <- simulateGrouped(n, 20L, 1, seed = 21L + binary)
  x <- sim$x
  g <- sim$g
  y <- if (binary) {
    as.integer(sim$eta + rnorm(n) > median(sim$eta))
  } else {
    sim$eta + rnorm(n)
  }
  set.seed(22L)
  fit <- fitRbart(
    y ~ x,
    group.by = g,
    n.chains = 2L,
    n.samples = 500L,
    n.burn = 500L,
    n.trees = 50L
  )
  list(
    tau = mean(fit$tau),
    correlation = cor(fit$ranef.mean, sim$b - mean(sim$b)),
    fit = fit
  )
}
gaussian <- recover(FALSE)
expect_true(gaussian$tau > 0.5 && gaussian$tau < 1.8)
expect_true(gaussian$correlation > 0.9)
probit <- recover(TRUE)
expect_true(probit$tau > 0.3 && probit$tau < 2.2)
expect_true(probit$correlation > 0.8)
rm(gaussian, probit, recover)

# --- the call shapes bartCause makes ---
sim <- simulateGrouped(80L, 4L, 1, seed = 31L)
df <- data.frame(x_1 = sim$x[, 1L], x_2 = sim$x[, 2L], x_3 = sim$x[, 3L])
df$y <- sim$eta + rnorm(sim$n)
df.test <- df[1:15, ]
g <- sim$g
data <- dbarts::dbartsData(y ~ x_1 + x_2 + x_3, df, test = df.test)
fit <- fitRbart(
  data,
  group.by = g,
  group.by.test = g[1:15],
  n.chains = 3L,
  n.samples = 6L
)
expect_inherits(fit, "rbart")
expect_equal(dim(fit$yhat.train), c(3L, 6L, 80L))
expect_true(!is.null(fit$fit))
extracted <- dbarts::extract(
  fit,
  type = "bart",
  sample = "test",
  combineChains = FALSE
)
predicted <- predict(
  fit,
  df.test[, c("x_1", "x_2", "x_3")],
  g[1:15],
  type = "bart",
  combineChains = FALSE
)
expect_equal(as.vector(extracted), as.vector(predicted))
rm(sim, df, df.test, g, data, fit, extracted, predicted)

rm(onceState, fitRbart, simulateGrouped)
