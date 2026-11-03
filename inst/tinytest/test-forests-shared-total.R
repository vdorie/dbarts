# Several forests with no basis (dec-B292, dec-B503, dec-B505, dec-B509): two
# or more that state no sd and no student() share one prior total, the fitting
# function's leaf prior, each an equal part of its variance whatever its tree
# count, held at 1, with one k for the group; they split the fitting
# function's tree count as blocks() splits a forest's; a stated sd or count is
# the forest's own, outside the share.

forest <- dbartsForests$forest
normal <- dbartsPriors$normal
fixed <- dbartsPriors$fixed

set.seed(20261030)
n <- 90L
x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
z <- factor(rbinom(n, 1L, 0.5))
y <- sin(3 * x1) + x2 + x3^2 + (z == "1") * 0.5 + rnorm(n, sd = 0.3)
yb <- as.double(y > stats::median(y))
x <- cbind(x1 = x1, x2 = x2, x3 = x3)
d <- data.frame(y, yb, x1, x2, x3, z)

ctl <- function(...) {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.samples = 5L,
    n.burn = 0L,
    updateState = FALSE,
    seed = 3L,
    ...
  )
}
familyArgs <- function(family) {
  if (family == "logistic") list(family = "logistic") else list()
}
response <- function(family) if (family == "gaussian") "y" else "yb"
sds <- function(sampler) {
  unname(vapply(sampler$getLeafPrior(), function(prior) unname(prior$sd), 0))
}
kScales <- function(sampler) {
  unname(vapply(sampler$getLeafPrior(), function(prior) prior$k.scale, 0))
}
# the one forest's prior sd at the same leaf prior: its k.scale over its k
loneSd <- function(family, ...) {
  lone <- do.call(
    dbarts,
    c(
      list(as.formula(paste(response(family), "~ x1 + x2 + x3")), d),
      familyArgs(family),
      list(control = ctl(), ...)
    )
  )
  prior <- lone$getLeafPrior()
  c(sd = prior$k.scale / as.vector(lone$getK())[1L], k.scale = prior$k.scale)
}

# --- every door and family: created, each sharer its part of the total ---
doors <- list(
  formula = function(family, ...) {
    do.call(
      dbarts,
      c(
        list(
          as.formula(paste(
            response(family),
            "~ forest(x1) + forest(x2) + forest(x3)"
          )),
          d
        ),
        familyArgs(family),
        list(control = ctl(), ...)
      )
    )
  },
  list = function(family, ...) {
    do.call(
      dbarts,
      c(
        list(x, d[[response(family)]]),
        familyArgs(family),
        list(
          forests = list(forest(x1), forest(x2), forest(x3)),
          control = ctl(),
          ...
        )
      )
    )
  },
  data = function(family, ...) {
    do.call(
      dbarts,
      c(
        list(dbartsData(
          x,
          d[[response(family)]],
          bases = list(NULL, NULL, NULL)
        )),
        familyArgs(family),
        list(control = ctl(), ...)
      )
    )
  }
)
for (family in c("gaussian", "probit", "logistic")) {
  lone <- loneSd(family)
  for (door in names(doors)) {
    info <- paste(family, door)
    sampler <- doors[[door]](family)
    priors <- sampler$getLeafPrior()
    expect_identical(length(priors), 3L, info = info)
    expect_true(
      all(vapply(
        priors,
        function(prior) identical(prior$prior.sd.of, "forest, shared"),
        NA
      )),
      info = info
    )
    # the total's k.scale over sqrt(3) on each: an equal part whatever the
    # count, 25 trees each here
    expect_equal(
      kScales(sampler),
      rep(lone[["k.scale"]] / sqrt(3), 3L),
      tolerance = 1e-12,
      check.attributes = FALSE,
      info = info
    )
    expect_identical(sampler$model@forests[[2L]]@n.trees, 25L)
    if (family == "gaussian") {
      expect_equal(
        sds(sampler),
        rep(lone[["sd"]] / sqrt(3), 3L),
        tolerance = 1e-12,
        check.attributes = FALSE,
        info = info
      )
      expect_identical(
        priors[[1L]]$leaf.prior,
        normal(k = 2),
        info = info
      )
    } else {
      # k is drawn, so the sd is no one number
      expect_true(all(is.na(sds(sampler))), info = info)
      expect_true(
        is(priors[[1L]]$leaf.prior@k, "dbartsChiHyperprior"),
        info = info
      )
    }
  }
}

# --- the parts are per forest at unequal counts: 38 and 37 trees, and 10, 38
# and 37 with a count stated; a part per tree would make them differ ---
total <- loneSd("gaussian")[["k.scale"]]
two <- dbarts(y ~ forest(x1) + forest(x2), d, control = ctl())
expect_identical(vapply(two$model@forests, slot, 0L, "n.trees"), c(38L, 37L))
expect_equal(kScales(two), rep(total / sqrt(2), 2L), tolerance = 1e-12)
stated <- dbarts(
  y ~ forest(x1, n.trees = 10) + forest(x2) + forest(x3),
  d,
  control = ctl()
)
expect_identical(
  vapply(stated$model@forests, slot, 0L, "n.trees"),
  c(10L, 38L, 37L)
)
expect_equal(kScales(stated), rep(total / sqrt(3), 3L), tolerance = 1e-12)
# the control reads the fitting function's count: the parts of the forests with
# no basis, summed
expect_identical(stated$control@n.trees, 85L)

# --- beside a forest with a basis and one stating sd: those keep their own,
# and a stated k or leaf.prior = normal(sd = ) is the group's ---
beside <- dbarts(
  y ~ forest(x1) + forest(x2) + forest(x3, sd = 0.5) + forest(x1, basis = z),
  d,
  leaf.prior = normal(k = 3),
  control = ctl()
)
besidePriors <- beside$getLeafPrior()
expect_identical(
  vapply(besidePriors, `[[`, "", "prior.sd.of"),
  c(
    forest1 = "forest, shared",
    forest2 = "forest, shared",
    forest3 = "forest",
    z = "level difference (at a coefficient of 1)"
  )
)
expect_equal(
  unname(sds(beside)[1:2]),
  rep(total / 3 / sqrt(2), 2L),
  tolerance = 1e-12
)
expect_identical(unname(sds(beside)[3L]), 0.5)
expect_identical(beside$getK()[1:2, 1L], c(3, 3))
expect_identical(beside$getK()[3:4, 1L], c(1, 1))
named <- dbarts(
  y ~ forest(x1) + forest(x2),
  d,
  leaf.prior = normal(sd = 1),
  control = ctl()
)
expect_equal(sds(named), rep(1 / sqrt(2), 2L), tolerance = 1e-12)
expect_identical(named$getLeafPrior(1L)$leaf.prior, normal(sd = 1))
expect_true(is(
  bart(
    y ~ forest(x1) + forest(x2) + forest(x3, basis = z),
    d,
    k = 3,
    n.samples = 3L,
    n.burn = 0L,
    n.chains = 1L,
    verbose = FALSE
  ),
  "bart"
))

# --- held at 1, and fixed() on a sharer changes nothing ---
heldControl <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.samples = 4L,
  n.burn = 2L,
  seed = 5L
)
plainFit <- dbarts(
  y ~ forest(x1) + forest(x2) + forest(x3, basis = z),
  d,
  control = heldControl
)
heldFit <- dbarts(
  y ~ forest(x1, amplitude = fixed()) + forest(x2) + forest(x3, basis = z),
  d,
  control = heldControl
)
plainRun <- plainFit$run()
expect_identical(heldFit$run()$train, plainRun$train)
amplitudes <- plainFit$getForestAmplitudes()
expect_identical(unname(amplitudes[1:2, ]), c(1, 1))

# --- one k for the group, drawn under probit ---
probit <- dbarts(
  yb ~ forest(x1) + forest(x2) + forest(x3, basis = z),
  d,
  control = dbartsControl(
    n.chains = 2L,
    n.threads = 1L,
    n.samples = 5L,
    n.burn = 5L,
    seed = 9L
  )
)
probitRun <- probit$run()
expect_identical(dim(probitRun$k), c(5L, 2L))
k <- probit$getK()
expect_identical(k[1L, ], k[2L, ])
expect_identical(k[2L, ], probitRun$k[5L, ])
expect_identical(k[3L, ], c(1, 1))
probitFit <- bart(
  yb ~ forest(x1) + forest(x2),
  d,
  n.samples = 5L,
  n.burn = 0L,
  n.chains = 2L,
  verbose = FALSE
)
expect_identical(length(probitFit$k), 10L)
expect_null(probitFit$fixed$k)

# --- counts: bart()'s n.trees reaches the forests stating none ---
counted <- bart(
  y ~ forest(x1, n.trees = 10) + forest(x2) + forest(x3),
  d,
  n.trees = 50,
  n.samples = 3L,
  n.burn = 0L,
  n.chains = 1L,
  verbose = FALSE
)
expect_identical(attr(counted, "forest.trees"), c(10L, 25L, 25L))
printed <- capture.output(print(counted))
expect_true("  trees: 10, 25, 25; 60 in all" %in% printed)
expect_true(
  "  forests 1, 2, 3 share one prior total, k 2" %in% printed
)
expect_true(any(grepl(
  "^  1 forest1: no basis; sd [0-9.]+, a share of the total$",
  printed
)))
expect_error(
  bart(
    y ~ forest(x1, n.trees = 10) + forest(x2, n.trees = 5),
    d,
    n.trees = 50,
    verbose = FALSE
  ),
  "'n.trees' is given to the fitting function and to every forest with no basis",
  fixed = TRUE
)
expect_error(
  dbarts(
    x,
    y,
    forests = list(forest(), forest(), forest()),
    control = dbartsControl(n.chains = 1L),
    n.trees = 2L
  ),
  "the fitting function's 2 tree(s) split among 3 forests with no basis",
  fixed = TRUE
)

# --- refusals, each by name ---
expect_error(
  dbarts(
    y ~ forest(x1) + forest(x2, amplitude = student(1)),
    d,
    control = ctl()
  ),
  "forest 2: amplitude = student() draws one forest's coefficient; these forests share one total",
  fixed = TRUE
)
ownCoefficient <- dbarts(
  y ~ forest(x1) + forest(x2) + forest(x3, sd = 0.5, amplitude = student(1)),
  d,
  control = ctl()
)
expect_identical(
  ownCoefficient$getLeafPrior(3L)$prior.sd.of,
  "forest (at a coefficient of 1)"
)
expect_error(
  bart(
    y ~ forest(x1 + x2) + forest(x3),
    d,
    interactions = interactions(max.order = 1),
    verbose = FALSE
  ),
  "'interactions' given to the fitting function is a constraint on the forest with no basis, and this model has 2",
  fixed = TRUE
)
expect_error(
  dbarts(
    x,
    y,
    forests = list(forest(), forest()),
    blocks = blocks(list("x1", c("x2", "x3"))),
    control = ctl()
  ),
  "'blocks' given to the fitting function is a constraint on the forest with no basis",
  fixed = TRUE
)
expect_error(
  dbarts(
    round(exp(y)) ~ forest(x1) + forest(x2),
    d,
    family = "nbinom",
    control = ctl()
  ),
  "does not support family \"nbinom\"",
  fixed = TRUE
)
# where no forest shares, the fitting function's leaf prior has no forest: the
# default restated is taken, anything else refused
expect_true(is(
  dbarts(
    y ~ forest(x1, sd = 1) + forest(x2, basis = z),
    d,
    leaf.prior = normal(2),
    control = ctl()
  ),
  "dbartsSampler"
))
expect_error(
  dbarts(
    y ~ forest(x1, sd = 1) + forest(x2, basis = z),
    d,
    leaf.prior = normal(k = 3),
    control = ctl()
  ),
  "no forest of this model takes a leaf prior's 'k'",
  fixed = TRUE
)
# a sharer's part is fixed at creation; the group's leaf prior changes
live <- dbarts(
  y ~ forest(x1) + forest(x2) + forest(x3, basis = z),
  d,
  control = ctl()
)
expect_error(
  live$setForest(1, sd = 1),
  "forest 1 shares one prior total with forests 1, 2",
  fixed = TRUE
)
expect_silent(live$setLeafPrior(normal(k = 3)))
expect_identical(unname(live$getK()[1:2, 1L]), c(3, 3))
live$setForest(3, sd = 0.4)
expect_identical(unname(live$getLeafPrior(3L)$sd), 0.4)

# --- a copy and a reload continue as one, and as the live sampler ---
kept <- dbarts(
  yb ~ forest(x1) + forest(x2) + forest(x3, basis = z),
  d,
  control = dbartsControl(
    n.chains = 2L,
    n.threads = 1L,
    n.samples = 5L,
    n.burn = 5L,
    seed = 13L,
    updateState = TRUE
  )
)
invisible(kept$run())
kept$storeState()
copied <- kept$copy()
path <- tempfile(fileext = ".rds")
saveRDS(kept, path)
reloaded <- readRDS(path)
unlink(path)
liveNext <- kept$run(0L, 5L)
copyNext <- copied$run(0L, 5L)
expect_identical(reloaded$run(0L, 5L), copyNext)
expect_equal(liveNext$train, copyNext$train, tolerance = 1e-10)
expect_identical(copied$getK(), reloaded$getK())

# --- a drawn k wherever the sharers stand ---
chi <- dbartsPriors$chi
invchi <- dbartsPriors$invchi
drawnFirst <- dbarts(
  x,
  yb,
  forests = list(
    forest(basis = z, sd = 0.5, amplitude = fixed()),
    forest(),
    forest()
  ),
  control = ctl()
)
drawnFirstRun <- drawnFirst$run()
expect_false(is.null(drawnFirstRun$k))
expect_identical(
  as.numeric(utils::tail(drawnFirstRun$k, 1L)),
  drawnFirst$getK()[2L]
)
expect_true(length(unique(as.numeric(drawnFirstRun$k))) > 1L)
statedFirst <- function(formula, ...) {
  bart(
    formula,
    d,
    n.samples = 12L,
    n.burn = 4L,
    n.chains = 2L,
    n.threads = 1L,
    seed = 5L,
    verbose = FALSE,
    ...
  )
}
statedFirstFit <- statedFirst(yb ~ forest(x1, sd = 1) + forest(x2) + forest(x3))
expect_false(is.null(statedFirstFit$k))
# extract names a drawn k per forest, as leaf.prior.sd does: the group's draws
# on a part of the share and the 1 a forest outside it is pinned at, and
# 'forest' selects among them as it does there
statedK <- extract(statedFirstFit, type = "k")
expect_identical(names(statedK), c("forest1", "forest2", "forest3"))
expect_identical(statedK$forest1, 1)
expect_identical(length(statedK$forest2), 24L)
expect_identical(as.vector(statedK$forest2), as.vector(statedFirstFit$k))
expect_identical(statedK$forest3, statedK$forest2)
expect_identical(extract(statedFirstFit, "k", forest = "forest1"), 1)
expect_identical(extract(statedFirstFit, "k", forest = 3L), statedK$forest3)
expect_identical(
  extract(statedFirstFit, "k", forest = c(1L, 3L)),
  statedK[c("forest1", "forest3")]
)
expect_identical(
  dim(extract(statedFirstFit, "k", combineChains = FALSE, forest = 2L)),
  c(2L, 12L)
)
# a name or an index that is no forest is refused, as the held sibling and
# leaf.prior.sd refuse it
for (type in c("k", "leaf.prior.sd")) {
  expect_error(
    extract(statedFirstFit, type, forest = "zz"),
    "'forest' names no forest of this model: \"zz\"",
    fixed = TRUE
  )
  expect_error(
    extract(statedFirstFit, type, forest = 9L),
    "'forest' index must be between 1 and 3",
    fixed = TRUE
  )
}

# the gaussian sweep and the grow-from-root sweep draw the group's k
sharedChi <- dbarts(
  y ~ forest(x1) + forest(x2),
  d,
  leaf.prior = normal(k = chi(1.25, Inf)),
  control = ctl()
)
expect_true(length(unique(as.numeric(sharedChi$run()$k))) > 1L)
sharedGrown <- dbarts(
  y ~ forest(x1) + forest(x2),
  d,
  leaf.prior = normal(k = chi(1.25, Inf)),
  control = ctl()
)
kBefore <- as.numeric(sharedGrown$getK())
invisible(sharedGrown$growFromRoot(5L))
kAfter <- as.numeric(sharedGrown$getK())
expect_true(kAfter[1L] != kBefore[1L])
expect_identical(kAfter[1L], kAfter[2L])

# a drawn named total starts at, and has its prior scale at, the statement
namedDrawn <- dbarts(
  y ~ forest(x1) + forest(x2),
  d,
  leaf.prior = normal(sd = invchi(1.25, 0.7)),
  control = ctl()
)
expect_equal(
  as.numeric(namedDrawn$getK()),
  rep(namedDrawn$getLeafPrior(1L)$k.scale * sqrt(2) / 0.7, 2L)
)

# a held k stays through the probit rescaling step
heldProbit <- dbarts(
  yb ~ forest(x1) + forest(x2),
  d,
  leaf.prior = normal(k = 3),
  control = ctl()
)
expect_true(is.null(heldProbit$run()$k))
expect_true(all(heldProbit$getK() == 3))

# a forest with no basis beside one stating sd shares with none, under probit
# as under gaussian
expect_true(is(
  dbarts(yb ~ forest(x1) + forest(x2, sd = 1), d, control = ctl()),
  "dbartsSampler"
))

# 'test' is refused by the formula's own message
expect_error(
  dbarts(y ~ forest(x1) + forest(x2), d, test = d[1:5, ], control = ctl()),
  "a forest() formula term does not support 'test'",
  fixed = TRUE
)

# the verbose creation summary names the group
verboseLines <- capture.output(dbarts(
  y ~ forest(x1) + forest(x2),
  d,
  control = dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.samples = 2L,
    n.burn = 1L,
    seed = 3L,
    verbose = TRUE
  )
))
expect_true(any(grepl(
  "forests sharing one prior total: 1 2, each 1/2 of its variance",
  verboseLines,
  fixed = TRUE
)))

# --- what a fit reads of a drawn k ---
binaryFit <- statedFirst(yb ~ forest(x1) + forest(x2))
gaussianFit <- statedFirst(y ~ forest(x1) + forest(x2))
gaussianChiFit <- statedFirst(
  y ~ forest(x1) + forest(x2),
  leaf.prior = normal(k = chi(1.25, Inf))
)
gaussianInvchiFit <- statedFirst(
  y ~ forest(x1) + forest(x2),
  leaf.prior = normal(sd = invchi(3, 1))
)
# summary shows the k row where it is drawn, wherever the sharers stand, as a
# fit of one forest does, and a held k keeps each forest's size
for (fit in list(
  binaryFit,
  gaussianChiFit,
  statedFirstFit
)) {
  expect_true("k" %in% summary(fit)$stats$variable)
  expect_false("leaf.prior.sd" %in% summary(fit)$vars)
}
# summary reads the fit's k, which extract's per-forest list leaves as it was
expect_equal(
  summary(statedFirstFit, vars = "k")$stats$mean,
  mean(statedFirstFit$k)
)
# a drawn total named by its sd shows leaf.prior.sd, as one forest does, and
# the row asked for by name is the parts' k.scale over each draw of k
expect_true("leaf.prior.sd" %in% summary(gaussianInvchiFit)$stats$variable)
expect_false("k" %in% summary(gaussianInvchiFit)$vars)
for (fit in list(
  binaryFit,
  gaussianChiFit,
  gaussianInvchiFit,
  statedFirstFit
)) {
  row <- summary(fit, vars = "leaf.prior.sd")$stats
  expect_identical(row$variable, "leaf.prior.sd")
  parts <- extract(fit, type = "leaf.prior.sd")
  expect_equal(row$mean, mean(parts[[length(parts)]]))
}
expect_false("k" %in% summary(gaussianFit)$vars)
expect_length(summary(gaussianFit)$fixed$leaf.prior.sd, 2L)
# leaf.prior.sd is each draw's part, k.scale over the group's k, and the
# stated forest's number beside it
partSd <- extract(binaryFit, type = "leaf.prior.sd")
expect_identical(names(partSd), c("forest1", "forest2"))
expect_identical(
  partSd$forest1,
  binaryFit$leaf.prior$forest1$k.scale /
    extract(binaryFit, type = "k", forest = 1L)
)
expect_identical(
  dim(partSd$forest2),
  dim(extract(binaryFit, type = "k", forest = 2L))
)
statedSd <- extract(statedFirstFit, type = "leaf.prior.sd")
expect_equivalent(statedSd$forest1, 1)
expect_identical(
  dim(statedSd$forest2),
  dim(extract(statedFirstFit, type = "k", forest = 2L))
)
# the chains kept apart, k and each part
for (type in c("k", "leaf.prior.sd")) {
  apart <- extract(statedFirstFit, type, combineChains = FALSE)
  expect_identical(dim(apart$forest2), c(2L, 12L), info = type)
  expect_identical(dim(apart$forest3), c(2L, 12L), info = type)
  expect_identical(length(apart$forest1), 1L, info = type)
}
expect_identical(
  extract(statedFirstFit, type = "leaf.prior.sd", forest = "forest3"),
  statedSd$forest3
)
# held, one number each, as before
expect_equal(
  unname(unlist(extract(gaussianFit, type = "leaf.prior.sd"))),
  unname(vapply(gaussianFit$leaf.prior, function(prior) unname(prior$sd), 0))
)

# fixed() on a sharer restates what the reader says: two spellings, one record
statedHold <- statedFirst(
  y ~ forest(x1, amplitude = fixed()) + forest(x2, amplitude = fixed())
)
expect_identical(gaussianFit$fixed$amplitude, c(forest1 = 1, forest2 = 1))
expect_identical(statedHold$fixed$amplitude, gaussianFit$fixed$amplitude)
beside <- statedFirst(y ~ forest(x1) + forest(x2) + forest(x3, sd = 1))
expect_identical(
  beside$fixed$amplitude,
  c(forest1 = 1, forest2 = 1, forest3 = 1)
)
