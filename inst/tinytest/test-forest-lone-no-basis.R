# A forest with no basis beside any other (dec-B509): held unless it states a
# prior, its leaf prior the fitting function's as one forest's is, whatever
# stands beside it.

forest <- dbartsForests$forest
normal <- dbartsPriors$normal
fixed <- dbartsPriors$fixed
student <- dbartsPriors$student
chi <- dbartsPriors$chi
invchi <- dbartsPriors$invchi

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
  expect_equal(unname(prior$k.scale), medianOf[[family]], info = info)
  expect_identical(
    prior$prior.sd.of,
    "forest (at a coefficient of 1)",
    info = info
  )
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
    "forest 1: amplitude = student() with no 'sd' takes its size from the fitting function's leaf prior, which here gives no size",
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
  do.call(sampler$setForest, c(list(1L), sampler$getForest(1L)))
  do.call(sampler$setForest, c(list(2L), sampler$getForest(2L)))
  expect_identical(sampler$getLeafPrior(), before, info = info)
  sampler$setForest(1, sd = 1.5)
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
expect_error(gaussian$setLeafPrior(normal(k = 3)), "$setForest", fixed = TRUE)
probit <- fit("probit", withBasis)
expect_error(probit$setLeafPrior(normal(k = 2)), "$setForest", fixed = TRUE)
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

## --- printed and read as the fit's: the leaf prior's size and its drawn k ----
fullBlock <- function(sampler) {
  dbarts:::samplerForestLines(sampler, full = TRUE)
}
expect_true(any(grepl(
  "no basis; sd [0-9.]+ \\(default\\)",
  fullBlock(fit("gaussian", withBasis))
)))
expect_true(any(grepl(
  "no basis; sd [0-9.]+ \\(leaf prior\\)",
  fullBlock(fit("gaussian", withBasis, leaf.prior = normal(k = 3)))
)))
drawnFit <- bart(
  yb ~ forest(x1) + forest(x2, basis = z),
  d,
  n.samples = 5L,
  n.burn = 0L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_identical(length(extract(drawnFit, type = "k", forest = 1L)), 5L)
expect_identical(length(unique(c(extract(drawnFit, "k", forest = 2L)))), 1L)
expect_true("k" %in% summary(drawnFit)$vars)

## --- the two spellings of one drawn prior give one median --------------------
kScale <- scaleOf[["gaussian"]]
chiMedian <- sqrt(qchisq(0.5, 3))
asK <- fit("gaussian", studentRhs, leaf.prior = normal(k = chi(3, 1)))
asSd <- fit(
  "gaussian",
  studentRhs,
  leaf.prior = normal(sd = invchi(3, kScale))
)
expect_equal(unname(reader(asK)$sd), kScale / chiMedian)
expect_equal(unname(reader(asSd)$sd), kScale / chiMedian)
expect_equal(asSd$run(0L, 4L)$train, asK$run(0L, 4L)$train, tolerance = 1e-8)
# the scale itself is not the median, which a scale of 1 shows plainly
expect_equal(
  unname(
    reader(fit(
      "gaussian",
      studentRhs,
      leaf.prior = normal(sd = invchi(3, 1))
    ))$sd
  ),
  1 / chiMedian
)
# beside a stated forest, and through a copy
keptCtl <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.samples = 4L,
  n.burn = 0L,
  updateState = TRUE,
  keepTrees = TRUE,
  seed = 3L
)
kept <- function(rhs, ...) {
  dbarts(as.formula(paste("y ~", rhs)), d, control = keptCtl, ...)
}
beside <- kept(
  "forest(x1, amplitude = student(1)) + forest(x2, sd = 2)",
  leaf.prior = normal(sd = invchi(3, kScale))
)
expect_equal(unname(reader(beside)$sd), kScale / chiMedian)
expect_identical(beside$copy()$getLeafPrior(1L)$sd, reader(beside)$sd)

## --- the sd in force written on a lone forest is the same prior --------------
for (family in families) {
  for (written in list(normal(k = 3), normal(sd = 1.1))) {
    info <- paste(family, deparse(written@k), deparse(written@prior.sd))
    build <- function() {
      args <- list(
        as.formula(paste(
          if (family == "gaussian") "y" else "yb",
          "~",
          withBasis
        )),
        d,
        control = keptCtl,
        leaf.prior = written
      )
      if (family == "logistic") {
        args$family <- "logistic"
      }
      do.call(dbarts, args)
    }
    wrote <- build()
    twin <- build()
    inForce <- unname(wrote$getLeafPrior(1L)$sd)
    wrote$setForest(1, sd = inForce)
    expect_identical(unname(wrote$getLeafPrior(1L)$sd), inForce, info = info)
    expect_identical(
      wrote$getK(forest = 1L),
      twin$getK(forest = 1L),
      info = info
    )
    expect_identical(
      wrote$run(0L, 6L)$train,
      twin$run(0L, 6L)$train,
      info = info
    )
  }
}

## --- a lone forest's written sd: what the sampler keeps of the leaf prior -----
for (written in list(normal(k = 3), normal(sd = 1.1))) {
  sampler <- kept(withBasis, leaf.prior = written)
  invisible(sampler$run(0L, 2L))
  sampler$setForest(1, sd = 0.4)
  # nothing of the leaf prior stays on the model, and a copy passes the bridge
  expect_null(sampler$model@leaf.prior@k)
  expect_null(sampler$model@leaf.prior@prior.sd)
  expect_identical(sampler$model@prior.scale, NA_real_)
  expect_equal(
    tryCatch(
      unname(sampler$copy()$getLeafPrior(1L)$sd),
      error = conditionMessage
    ),
    0.4
  )
}

## --- $setLeafPrior(leaf.prior) on a lone forest restates, and only that -------
sampler <- kept(withBasis, leaf.prior = normal(sd = 1.1))
expect_error(sampler$setLeafPrior(normal(sd = 9)), "state a forest's size")
expect_equal(unname(reader(sampler)$sd), 1.1)
sampler <- kept(withBasis)
expect_silent(sampler$setLeafPrior(normal(k = 2L)))

## --- a size the leaf prior cannot give is refused only where it is asked for --
expect_error(
  kept(studentRhs, leaf.prior = normal(k = Inf)),
  "student.*gives no size"
)
flat <- normal(k = chi(1.5, Inf))
flatTaken <- function(rhs) {
  tryCatch(
    {
      kept(rhs, leaf.prior = flat)
      TRUE
    },
    error = conditionMessage
  )
}
expect_true(flatTaken(
  "forest(x1) + forest(x2, basis = z, amplitude = student(1))"
))
expect_true(flatTaken(
  "forest(x1, sd = 0.4) + forest(x2) + forest(x2, basis = z)"
))
expect_true(flatTaken(
  "forest(x1, sd = 0.4, amplitude = student(1)) + forest(x2) + forest(x2, basis = z)"
))

## --- "(default)" is printed only where it is true ----------------------------
lineOf <- function(sampler) fullBlock(sampler)[2L]
expect_true(endsWith(lineOf(fit("gaussian", withBasis)), "(default)"))
expect_true(endsWith(
  lineOf(fit("gaussian", withBasis, leaf.prior = normal(k = 2))),
  "(default)"
))
for (family in c("probit", "logistic")) {
  expect_true(
    endsWith(
      lineOf(fit(family, withBasis, leaf.prior = normal(k = 2))),
      "(leaf prior)"
    ),
    info = family
  )
}
expect_true(endsWith(
  lineOf(fit("gaussian", studentRhs, leaf.prior = normal(k = 3))),
  "(leaf prior); coefficient student(1)"
))
expect_true(endsWith(
  lineOf(fit(
    "gaussian",
    "forest(x1, sd = 2, amplitude = student(1)) + forest(x2, basis = z)"
  )),
  "(stated); coefficient student(1)"
))
# the short line carries no source
expect_false(any(grepl(
  "(leaf prior)",
  dbarts:::samplerForestLines(fit("gaussian", studentRhs)),
  fixed = TRUE
)))
# a named drawn total prints as its prior, on the sd
expect_true(grepl(
  "no basis; sd invchi(1.5, 2)",
  lineOf(fit("gaussian", withBasis, leaf.prior = normal(sd = invchi(1.5, 2)))),
  fixed = TRUE
))

## --- verbose creation names a group only where there is one -------------------
verboseLines <- function(rhs) {
  capture.output(dbarts(
    as.formula(paste("y ~", rhs)),
    d,
    control = ctl(verbose = TRUE)
  ))
}
expect_false(any(grepl("sharing one prior total", verboseLines(withBasis))))
expect_true(any(grepl(
  "sharing one prior total",
  verboseLines("forest(x1) + forest(x2) + forest(x2, basis = z)")
)))
