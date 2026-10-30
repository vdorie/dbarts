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
    expect_identical(
      attr(sampler$control, "bartcore.forests", exact = TRUE)$params[[2L]][1L],
      25
    )
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
expect_identical(
  vapply(
    attr(two$control, "bartcore.forests")$params,
    function(p) p[[1L]],
    0
  ),
  c(38, 37)
)
expect_equal(kScales(two), rep(total / sqrt(2), 2L), tolerance = 1e-12)
stated <- dbarts(
  y ~ forest(x1, n.trees = 10) + forest(x2) + forest(x3),
  d,
  control = ctl()
)
expect_identical(
  vapply(
    attr(stated$control, "bartcore.forests")$params,
    function(p) p[[1L]],
    0
  ),
  c(10, 38, 37)
)
expect_equal(kScales(stated), rep(total / sqrt(3), 3L), tolerance = 1e-12)
# the control the next fit inherits keeps the fitting function's count
expect_identical(stated$control@n.trees, 10L)
expect_identical(
  attr(stated$control, "bartcore.forests")$control.n.trees,
  75L
)

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
    forest3 = "forest (prior median)",
    z = "level difference (prior median)"
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
    control = dbartsControl(n.trees = 2L, n.chains = 1L)
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
  "forest (prior median)"
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
  "a non-default 'k'",
  fixed = TRUE
)
# a sharer's part is fixed at creation
live <- dbarts(
  y ~ forest(x1) + forest(x2) + forest(x3, basis = z),
  d,
  control = ctl()
)
expect_error(
  live$setLeafPrior(forests = list(forest(sd = 1))),
  "states forest 1's 'sd', but it shares one prior total with forests 1, 2",
  fixed = TRUE
)
expect_error(
  live$setLeafPrior(normal(k = 3)),
  "takes no 'leaf.prior': the total is the leaf prior they were created under",
  fixed = TRUE
)
live$setLeafPrior(forests = list(forest(), forest(), forest(sd = 0.4)))
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
