# A forest with no basis beside any other (dec-B509): held unless it states a
# prior, its leaf prior the fitting function's as one forest's is, whatever
# stands beside it.

forest <- dbartsForests$forest
normal <- dbartsPriors$normal
fixed <- dbartsPriors$fixed
student <- dbartsPriors$student
chi <- dbartsPriors$chi

set.seed(20261030)
n <- 90L
x1 <- runif(n)
x2 <- runif(n)
z <- factor(rbinom(n, 1L, 0.5))
y <- sin(3 * x1) + x2 + (z == "1") * 0.5 + rnorm(n, sd = 0.3)
yb <- as.double(y > stats::median(y))
d <- data.frame(y, yb, x1, x2, z)

ctl <- function(...) {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.samples = 4L,
    n.burn = 0L,
    updateState = FALSE,
    seed = 3L,
    ...
  )
}
families <- c("gaussian", "probit", "logistic")
fit <- function(family, rhs, ...) {
  formula <- as.formula(paste(
    if (family == "gaussian") "y" else "yb",
    "~",
    rhs
  ))
  args <- list(formula, d, control = ctl(), ...)
  if (family == "logistic") {
    args$family <- "logistic"
  }
  do.call(dbarts, args)
}
withBasis <- "forest(x1) + forest(x2, basis = z)"
withSd <- "forest(x1) + forest(x2, sd = 2)"
# the coefficient at every draw, and the k at every draw
coefficientDraws <- function(sampler, forest = 1L) {
  vapply(
    1:4,
    function(i) {
      sampler$run(0L, 1L)
      c(sampler$getForestAmplitudes(forest))[1L]
    },
    0
  )
}
kDraws <- function(sampler) {
  vapply(
    1:4,
    function(i) {
      sampler$run(0L, 1L)
      sampler$getK(forest = 1L)[[1L]]
    },
    0
  )
}
reader <- function(sampler) sampler$getLeafPrior(forest = 1L)
# one forest's own leaf prior under the family
scaleOf <- c(
  gaussian = diff(range(y)) / 2,
  probit = 3,
  logistic = 3 * pi / sqrt(3)
)
unitOf <- c(gaussian = sd(y), probit = 1, logistic = pi / sqrt(3))

## --- the default, beside a forest with a basis and beside one stating sd -----
for (family in families) {
  for (rhs in c(withBasis, withSd)) {
    info <- paste(family, rhs)
    sampler <- fit(family, rhs)
    prior <- reader(sampler)
    expect_equal(unname(prior$k.scale), scaleOf[[family]], info = info)
    expect_identical(prior$prior.sd.of, "forest", info = info)
    expect_true(is(prior$amplitude, "dbartsFixedPrior"), info = info)
    expect_true(all(coefficientDraws(sampler) == 1), info = info)
    if (family == "gaussian") {
      expect_equal(unname(prior$sd), scaleOf[[family]] / 2, info = info)
      expect_identical(prior$leaf.prior@k, 2, info = info)
      expect_identical(kDraws(sampler), rep(2, 4L), info = info)
    } else {
      expect_identical(unname(prior$sd), NA_real_, info = info)
      expect_true(is(prior$leaf.prior@k, "dbartsChiHyperprior"), info = info)
      expect_identical(prior$leaf.prior@k@degreesOfFreedom, 1.5, info = info)
      expect_identical(prior$leaf.prior@k@scale, 2, info = info)
      expect_true(length(unique(kDraws(sampler))) > 1L, info = info)
    }
  }
}

## --- the fitting function's leaf prior: taken, not refused ------------------
for (family in families) {
  scale <- scaleOf[[family]]
  atK <- fit(family, withBasis, leaf.prior = normal(k = 3))
  expect_equal(unname(reader(atK)$sd), scale / 3, info = family)
  expect_identical(atK$getK(forest = 1L)[[1L]], 3, info = family)
  expect_identical(atK$getK(forest = 2L)[[1L]], 1, info = family)
  drawn <- fit(family, withBasis, leaf.prior = normal(k = chi(1.5, 3)))
  expect_identical(reader(drawn)$leaf.prior@k@scale, 3, info = family)
  expect_true(length(unique(kDraws(drawn))) > 1L, info = family)
  named <- fit(family, withBasis, leaf.prior = normal(sd = 1.1))
  expect_equal(unname(reader(named)$sd), 1.1, info = family)
  expect_equal(unname(reader(named)$k.scale), scale, info = family)
  expect_equal(named$getK(forest = 1L)[[1L]], scale / 1.1, info = family)
  expect_true(all(coefficientDraws(named) == 1), info = family)
}
expect_true(is(
  bart(
    y ~ forest(x1) + forest(x2, basis = z),
    d,
    k = 3,
    n.samples = 3L,
    n.burn = 0L,
    verbose = FALSE
  ),
  "bart"
))
# the forests list and a data object read the same
frame <- cbind(x1 = x1, x2 = x2)
listed <- dbarts(
  frame,
  y,
  forests = list(forest(), forest(basis = z)),
  leaf.prior = normal(k = 3),
  control = ctl()
)
expect_equal(unname(reader(listed)$sd), scaleOf[["gaussian"]] / 3)
fromData <- dbarts(
  dbartsData(frame, y, bases = list(NULL, z)),
  leaf.prior = normal(k = 3),
  control = ctl()
)
expect_equal(unname(reader(fromData)$sd), scaleOf[["gaussian"]] / 3)

## --- fixed() restates the default; a stated sd is held at it ----------------
for (family in families) {
  info <- family
  plain <- fit(family, withBasis)
  restated <- fit(
    family,
    "forest(x1, amplitude = fixed()) + forest(x2, basis = z)"
  )
  expect_identical(restated$run()$train, plain$run()$train, info = info)
  stated <- fit(family, "forest(x1, sd = 0.7) + forest(x2, basis = z)")
  statedHeld <- fit(
    family,
    "forest(x1, sd = 0.7, amplitude = fixed()) + forest(x2, basis = z)"
  )
  expect_identical(unname(reader(stated)$sd), 0.7, info = info)
  expect_equivalent(reader(stated)$k.scale, 0.7, info = info)
  expect_true(is(reader(stated)$amplitude, "dbartsFixedPrior"), info = info)
  expect_identical(stated$getK(forest = 1L)[[1L]], 1, info = info)
  expect_identical(statedHeld$run()$train, stated$run()$train, info = info)
  expect_true(all(coefficientDraws(stated) == 1), info = info)
}

## --- a stated prior draws, its median the size the forest would hold --------
medianOf <- c(
  gaussian = diff(range(y)) / 4,
  probit = 3 / (2 * sqrt(stats::qchisq(0.5, 1.5))),
  logistic = 3 * pi / sqrt(3) / (2 * sqrt(stats::qchisq(0.5, 1.5)))
)
studentRhs <- "forest(x1, amplitude = student(1)) + forest(x2, basis = z)"
for (family in families) {
  info <- family
  drawn <- fit(family, studentRhs)
  prior <- reader(drawn)
  expect_equal(unname(prior$sd), medianOf[[family]], info = info)
  expect_equal(unname(prior$k.scale), unitOf[[family]], info = info)
  expect_identical(prior$prior.sd.of, "forest (prior median)", info = info)
  expect_identical(drawn$getK(forest = 1L)[[1L]], 1, info = info)
  expect_true(length(unique(coefficientDraws(drawn))) > 1L, info = info)
  # at a held k, a named sd, and a stated sd
  atK <- fit(family, studentRhs, leaf.prior = normal(k = 3))
  expect_equal(unname(reader(atK)$sd), scaleOf[[family]] / 3, info = info)
  named <- fit(family, studentRhs, leaf.prior = normal(sd = 1.1))
  expect_equal(unname(reader(named)$sd), 1.1, info = info)
  stated <- fit(
    family,
    "forest(x1, sd = 0.7, amplitude = student(1)) + forest(x2, basis = z)"
  )
  expect_identical(unname(reader(stated)$sd), 0.7, info = info)
  # a flat prior on k gives no size
  expect_error(
    fit(family, studentRhs, leaf.prior = normal(k = chi(1.5, Inf))),
    "forest 1: amplitude = student() with no 'sd' takes its median from the fitting function's leaf prior, which here gives no size",
    fixed = TRUE,
    info = info
  )
  expect_true(
    is(
      fit(family, withBasis, leaf.prior = normal(k = chi(1.5, Inf))),
      "dbartsSampler"
    ),
    info = info
  )
}
# one forest has no coefficient: fixed() restates it, student() is refused
expect_identical(
  dbarts(y ~ forest(x1, amplitude = fixed()), d, control = ctl())$run()$train,
  dbarts(y ~ forest(x1), d, control = ctl())$run()$train
)
expect_error(
  dbarts(y ~ forest(x1, amplitude = student(1)), d, control = ctl()),
  "this model has one forest, which has no coefficient to give a prior",
  fixed = TRUE
)

## --- with no forest to take it, the leaf prior is refused by name -----------
stated <- "forest(x1, sd = 1) + forest(x2, sd = 2)"
expect_true(is(
  fit("gaussian", stated, leaf.prior = normal(2)),
  "dbartsSampler"
))
for (prior in list(normal(k = 3), normal(k = chi(1.5, 2)))) {
  expect_error(
    fit("gaussian", stated, leaf.prior = prior),
    "no forest of this model takes a leaf prior's 'k'",
    fixed = TRUE
  )
}
expect_error(
  fit("gaussian", stated, leaf.prior = normal(sd = 1)),
  "no forest of this model takes a leaf prior's 'sd'",
  fixed = TRUE
)

## --- $setLeafPrior: a stated sd leaves the leaf prior ----------------------
for (family in families) {
  info <- family
  sampler <- fit(family, withBasis)
  sampler$run(0L, 1L)
  # the reader's own value goes back as nothing written
  before <- sampler$getLeafPrior()
  sampler$setLeafPrior(forests = lapply(before, `[[`, "leaf.prior"))
  expect_identical(sampler$getLeafPrior(), before, info = info)
  sampler$setLeafPrior(forests = list(forest(sd = 1.5), forest()))
  prior <- reader(sampler)
  expect_identical(unname(prior$sd), 1.5, info = info)
  expect_equivalent(prior$k.scale, 1.5, info = info)
  expect_identical(sampler$getK(forest = 1L)[[1L]], 1, info = info)
  expect_true(all(coefficientDraws(sampler) == 1), info = info)
  # no forest takes the leaf prior now, so a copy passes the bridge and runs on
  copied <- sampler$copy()
  expect_identical(unname(reader(copied)$sd), 1.5, info = info)
  expect_identical(copied$run(0L, 1L)$k, sampler$run(0L, 1L)$k, info = info)
}
# restating what is in force is taken; another k through this door is not
gaussian <- fit("gaussian", withBasis)
gaussian$setLeafPrior(normal(k = 2))
expect_equivalent(reader(gaussian)$k.scale, scaleOf[["gaussian"]])
expect_error(gaussian$setLeafPrior(normal(k = 3)), "forests =", fixed = TRUE)
probit <- fit("probit", withBasis)
expect_error(probit$setLeafPrior(normal(k = 2)), "forests =", fixed = TRUE)
probit$setLeafPrior(reader(probit)$leaf.prior)

## --- a copy and a reload continue as one, and as the live sampler -----------
for (family in families) {
  for (rhs in c(withBasis, withSd)) {
    info <- paste(family, rhs)
    args <- list(
      as.formula(paste(if (family == "gaussian") "y" else "yb", "~", rhs)),
      d,
      control = dbartsControl(
        n.chains = 1L,
        n.threads = 1L,
        n.samples = 3L,
        n.burn = 3L,
        seed = 13L,
        updateState = TRUE
      )
    )
    if (family == "logistic") {
      args$family <- "logistic"
    }
    kept <- do.call(dbarts, args)
    invisible(kept$run())
    kept$storeState()
    copied <- kept$copy()
    path <- tempfile(fileext = ".rds")
    saveRDS(kept, path)
    reloaded <- readRDS(path)
    unlink(path)
    liveNext <- kept$run(0L, 3L)
    copyNext <- copied$run(0L, 3L)
    expect_identical(reloaded$run(0L, 3L), copyNext, info = info)
    expect_equal(liveNext$train, copyNext$train, tolerance = 1e-10, info = info)
    expect_identical(copied$getK(), reloaded$getK(), info = info)
  }
}
