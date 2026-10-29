# A forest's own proposal.probs, n.perturb.cuts and treeShift: each one stated
# governs its forest and each one left out takes the control's, which
# setControl changes on the forests that state none. Twins pin each path
# bitwise: a forest stating X under a control holding Y draws what the
# control holding X draws.

set.seed(17L)
n <- 60L
x <- matrix(runif(n * 3L), n, 3L, dimnames = list(NULL, c("a", "b", "c")))
dat <- data.frame(x, z = factor(sample(c("t", "c"), n, TRUE)))
dat$y <- sin(3 * dat$a) + dat$b - dat$c + 0.3 * rnorm(n)
labels <- cut(dat$a + 0.2 * rnorm(n), 3L, labels = FALSE)
counts <- sapply(1:3, function(k) as.integer(labels == k))

shipped <- c(birth_death = 0.6, change = 0.4)
moving <- c(
  birth_death = 0.4,
  swap = 0.1,
  change = 0.2,
  perturb = 0.2,
  rule_gibbs = 0.1,
  birth = 0.45
)
# rule_gibbs does nothing on the variance forest's leaves, which refuses it
varianceMoving <- c(birth_death = 0.5, swap = 0.1, change = 0.2, perturb = 0.2)
frozen <- c(birth_death = 0, swap = 0, change = 0, perturb = 0, rule_gibbs = 0)

forest <- dbarts::dbartsForests$forest
varianceForest <- dbarts::dbartsForests$varianceForest
ctl <- function(...) {
  dbarts::dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 4L,
    n.burn = 0L,
    n.samples = 3L,
    updateState = FALSE,
    ...
  )
}
# `mean` states the first forest's settings, `second` the basis forest's and
# `variance` the variance forest's; NULL is not stated
build <- function(
  kind,
  control,
  mean = list(),
  second = mean,
  variance = list()
) {
  set.seed(23L)
  first <- forest(
    proposal.probs = mean$proposal.probs,
    n.perturb.cuts = mean$n.perturb.cuts,
    treeShift = mean$treeShift
  )
  switch(
    kind,
    one = dbarts::dbarts(
      y ~ a + b + c,
      dat,
      forests = list(first),
      control = control
    ),
    hetero = dbarts::dbarts(
      y ~ a + b + c,
      dat,
      forests = list(first),
      variance = varianceForest(
        n.trees = 3L,
        proposal.probs = variance$proposal.probs,
        n.perturb.cuts = variance$n.perturb.cuts
      ),
      control = control
    ),
    two = dbarts::dbarts(
      y ~ a + b + c,
      dat,
      forests = list(
        first,
        forest(
          a + b,
          basis = z,
          proposal.probs = second$proposal.probs,
          n.perturb.cuts = second$n.perturb.cuts,
          treeShift = second$treeShift
        )
      ),
      control = control
    ),
    multi = dbarts::dbarts(
      dbarts::dbartsData(x, counts = counts),
      family = "multinomial",
      forests = list(first),
      control = control
    ),
    dart = dbarts::dbarts(
      y ~ a + b + c,
      dat,
      forests = list(first),
      split.probs = dart(),
      control = control
    )
  )
}
# a run, then the live sampler's continuation, which its copy and its reload
# after storeState run as well, the two re-created ones bitwise and the live
# one to its fits' reassociation; `setup` acts between the run and the rest
draws <- function(sampler, setup = NULL) {
  set.seed(29L)
  first <- sampler$run(20L, 3L)
  if (!is.null(setup)) {
    setup(sampler)
  }
  sampler$storeState()
  duplicate <- sampler$copy()
  file <- tempfile(fileext = ".rds")
  saveRDS(sampler, file)
  reloaded <- readRDS(file)
  unlink(file)
  continued <- lapply(list(sampler, duplicate, reloaded), function(s) {
    set.seed(31L)
    s$run(0L, 3L)
  })
  expect_identical(continued[[3L]], continued[[2L]])
  expect_equal(continued[[2L]], continued[[1L]])
  sampler$storeState()
  list(first = first, continued = continued, state = sampler$state)
}
# the control's three settings set from `args`, the rest kept
restate <- function(sampler, args) {
  given <- do.call(ctl, args)
  control <- sampler$control
  control@proposal.probs <- given@proposal.probs
  control@n.perturb.cuts <- given@n.perturb.cuts
  control@levelGibbs <- given@levelGibbs
  sampler$setControl(control)
}
restated <- function(sampler, args) {
  restate(sampler, args)
  sampler
}

# ---- cross-twins: a forest's own is the control's, on every path ----

cases <- list(
  mixture = list(
    x = list(proposal.probs = moving),
    y = list(proposal.probs = shipped)
  ),
  window3 = list(
    x = list(n.perturb.cuts = 3),
    y = list(proposal.probs = moving, n.perturb.cuts = 1)
  ),
  windowInf = list(
    x = list(n.perturb.cuts = Inf),
    y = list(proposal.probs = moving, n.perturb.cuts = 1)
  ),
  always = list(x = list(treeShift = "always"), y = list(treeShift = "never")),
  never = list(x = list(treeShift = "never"), y = list(treeShift = "always")),
  frozenAlways = list(
    x = list(proposal.probs = frozen, treeShift = "always"),
    y = list(proposal.probs = shipped, treeShift = "never")
  ),
  frozenNever = list(
    x = list(proposal.probs = frozen, treeShift = "never"),
    y = list(proposal.probs = shipped, treeShift = "always")
  )
)
for (kind in c("one", "hetero", "two", "multi")) {
  for (case in names(cases)) {
    X <- cases[[case]]$x
    Y <- cases[[case]]$y
    stated <- draws(build(kind, do.call(ctl, Y), mean = X))
    info <- paste(kind, case)
    # the twin's variance forest states what the control held
    twin <- build(
      kind,
      do.call(ctl, modifyList(Y, X)),
      variance = Y[setdiff(names(X), "treeShift")]
    )
    expect_identical(stated, draws(twin), info = info)
    expect_false(
      identical(stated, draws(build(kind, do.call(ctl, Y)))),
      info = info
    )
  }
}
# the variance forest's, against a twin whose mean forest states the control's
cases$mixture$x <- list(proposal.probs = varianceMoving)
for (case in c("mixture", "window3", "windowInf")) {
  X <- cases[[case]]$x
  Y <- cases[[case]]$y
  stated <- draws(build("hetero", do.call(ctl, Y), variance = X))
  twin <- build("hetero", do.call(ctl, modifyList(Y, X)), mean = Y[names(X)])
  expect_identical(stated, draws(twin), info = case)
  expect_false(
    identical(stated, draws(build("hetero", do.call(ctl, Y)))),
    info = case
  )
}

# a stated "auto" is not an unstated one: on a forest frozen by its own
# mixture it takes the step under a control's "never"
autoFrozen <- list(proposal.probs = frozen, treeShift = "auto")
expect_identical(
  draws(build("one", ctl(treeShift = "never"), mean = autoFrozen)),
  draws(build("one", ctl(proposal.probs = frozen)))
)
expect_false(identical(
  draws(build("one", ctl(treeShift = "never"), mean = autoFrozen)),
  draws(build(
    "one",
    ctl(treeShift = "never"),
    mean = autoFrozen["proposal.probs"]
  ))
))

# ---- forests apart ----

meanShape <- function(sampler) {
  trees <- sampler$getTrees()
  trees$value[trees$var == -1L] <- NA
  trees
}
varianceShape <- function(sampler) {
  sampler$storeState()
  chainState <- sampler$state[[1L]]
  splits <- chainState$variance.vars >= 0L
  list(
    chainState$variance.vars,
    chainState$variance.sizes,
    matrix(chainState$variance.values, 8L)[, splits, drop = FALSE]
  )
}
sweepsMoved <- function(sampler, shapeOf) {
  at <- shapeOf(sampler)
  moved <- 0L
  for (sweep in seq_len(15L)) {
    invisible(sampler$run(0L, 1L))
    moved <- moved + !identical(shapeOf(sampler), at)
  }
  moved
}
# a control frozen after some sweeps holds the forest that states no mixture
# and leaves the one that states its own moving
apart <- build(
  "hetero",
  ctl(),
  variance = list(proposal.probs = varianceMoving)
)
invisible(apart$run(20L, 1L))
restate(apart, list(proposal.probs = frozen))
expect_identical(sweepsMoved(apart, meanShape), 0L)
expect_true(sweepsMoved(apart, varianceShape) > 0L)
apart <- build("hetero", ctl(), mean = list(proposal.probs = moving))
invisible(apart$run(20L, 1L))
restate(apart, list(proposal.probs = frozen))
expect_true(sweepsMoved(apart, meanShape) > 0L)
expect_identical(sweepsMoved(apart, varianceShape), 0L)
# a second forest frozen by its own stays stumps while the first grows
apart <- build("two", ctl(), second = list(proposal.probs = frozen))
invisible(apart$run(20L, 1L))
trees <- apart$getTrees()
expect_true(all(trees$var[trees$forest == 2L] == -1L))
expect_true(any(trees$var[trees$forest == 1L] != -1L))
rm(apart, trees)

# ---- setControl reaches exactly the forests that state none ----

changes <- list(
  mixture = list(
    x = list(proposal.probs = moving),
    y = list(),
    z = list(proposal.probs = frozen)
  ),
  window = list(
    x = list(n.perturb.cuts = 3),
    y = list(proposal.probs = moving),
    z = list(proposal.probs = moving, n.perturb.cuts = Inf)
  ),
  treeShift = list(
    x = list(treeShift = "never"),
    y = list(),
    z = list(treeShift = "always")
  )
)
for (change in names(changes)) {
  X <- changes[[change]]$x
  Y <- changes[[change]]$y
  Z <- changes[[change]]$z
  # before any run: the first forest states X, the other takes Z
  expect_identical(
    draws(restated(
      build("two", do.call(ctl, Y), mean = X, second = list()),
      Z
    )),
    draws(build("two", do.call(ctl, Y), mean = X, second = Z[names(X)])),
    info = change
  )
  # created under Y, then given Z, as created under Z
  for (kind in c("one", "hetero", "two", "multi", "dart")) {
    expect_identical(
      draws(restated(build(kind, do.call(ctl, Y)), Z)),
      draws(build(kind, do.call(ctl, Z))),
      info = paste(kind, change)
    )
  }
  if (change != "treeShift") {
    expect_identical(
      draws(restated(build("hetero", do.call(ctl, Y), mean = X), Z)),
      draws(build("hetero", do.call(ctl, Y), mean = X, variance = Z[names(X)])),
      info = change
    )
  }
  # one entry states every category forest's, which setControl then leaves
  expect_identical(
    draws(restated(build("multi", do.call(ctl, Y), mean = X), Z)),
    draws(build("multi", do.call(ctl, Y), mean = X)),
    info = change
  )
}
# a moving mixture and treeShift changed together take the prior install,
# which installs the mode with them
for (kind in c("one", "hetero")) {
  expect_identical(
    draws(restated(
      build(kind, ctl(treeShift = "never")),
      list(proposal.probs = moving, treeShift = "always")
    )),
    draws(build(kind, ctl(proposal.probs = moving, treeShift = "always"))),
    info = kind
  )
}
# a DART sampler takes a change after some sweeps, its copy and its reload
# running what it runs
dartSampler <- build("dart", ctl())
invisible(draws(dartSampler, function(s) {
  expect_null(restate(s, list(proposal.probs = moving, treeShift = "always")))
}))
expect_identical(dartSampler$control@proposal.probs, moving)
expect_true(dartSampler$control@levelGibbs)
rm(dartSampler)

# treeShift alone takes the move install, which leaves a sigma a host wrote
# to the engine alone; with the mixture the prior install puts the model's
# fixed value back
fixedSigma <- function(args) {
  sampler <- dbarts::dbarts(
    y ~ a + b + c,
    dat,
    family = gaussian(sigma = fixed(1)),
    control = ctl()
  )
  .Call(dbarts:::C_dbarts_bartcore_setSigma, sampler$getPointer(), 0.3)
  restate(sampler, args)
  sampler$run(0L, 1L)$sigma
}
expect_equal(fixedSigma(list(treeShift = "always")), 0.3)
expect_equal(fixedSigma(list(treeShift = "always", proposal.probs = moving)), 1)

# ---- no partial writes ----

# three entries, for a sampler of one forest
wrongCount <- function(sampler) {
  control <- sampler$control
  entry <- list(proposal.probs = moving, n.perturb.cuts = 3, levelGibbs = TRUE)
  attr(control, "bartcore.moves") <- list(
    forests = rep(list(entry), 3L),
    variance = NULL
  )
  control
}
handed <- build("hetero", ctl())
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_setMoveSettings,
    handed$getPointer(),
    wrongCount(handed)
  ),
  "do not match this sampler's forests"
)
model <- handed$model
model@tree.prior@power <- 3
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_setModel,
    handed$getPointer(),
    model,
    handed$data,
    wrongCount(handed)
  ),
  "do not match this sampler's forests"
)
expect_identical(draws(handed), draws(build("hetero", ctl())))
rm(handed, model)
# a refused install through setControl rolls the stored control back
handed <- build("one", ctl())
recorded <- handed$control
handed$control <- wrongCount(handed)
expect_error(
  restate(handed, list(proposal.probs = frozen)),
  "do not match this sampler's forests"
)
expect_identical(
  unname(handed$control@proposal.probs),
  c(0.6, 0, 0.4, 0, 0, 0.5)
)
handed$control <- recorded
expect_identical(draws(handed), draws(build("one", ctl())))
rm(handed, recorded)

# ---- the prior installs keep a forest's own ----

all3 <- list(proposal.probs = moving, n.perturb.cuts = 3, treeShift = "always")
expect_identical(
  draws(build("one", ctl(), mean = all3), function(s) s$setModel(s$model)),
  draws(build("one", ctl(), mean = all3))
)
newHyperprior <- function(s) s$setLeafPrior(normal(k = chi(df = 2, scale = 3)))
expect_identical(
  draws(build("one", ctl(), mean = all3), newHyperprior),
  draws(build("one", do.call(ctl, all3)), newHyperprior)
)

# ---- refusals ----

refusals <- list(
  list(
    list(proposal.probs = c(birth_death = 0.9, change = 0.4)),
    "'proposal.probs'"
  ),
  list(
    list(proposal.probs = c(birth_death = 0.9, change = 0.1, birth = 2)),
    "'proposal.probs'"
  ),
  list(list(proposal.probs = c(walk = 1)), "'proposal.probs'"),
  list(list(n.perturb.cuts = 0), "'n.perturb.cuts'"),
  list(list(n.perturb.cuts = 1.5), "'n.perturb.cuts'"),
  list(list(n.perturb.cuts = 70000), "'n.perturb.cuts'"),
  list(list(n.perturb.cuts = "3"), "'n.perturb.cuts'"),
  list(list(treeShift = "sometimes"), "'treeShift'"),
  list(list(treeShift = TRUE), "'treeShift'")
)
for (refusal in refusals) {
  bad <- refusal[[1L]]
  expect_error(
    build("two", ctl(), second = bad),
    paste0("forest 2's ", refusal[[2L]]),
    fixed = TRUE
  )
  expect_error(
    dbarts::dbarts(
      y ~ forest(
        a + b,
        n.perturb.cuts = bad$n.perturb.cuts,
        treeShift = bad$treeShift,
        proposal.probs = bad$proposal.probs
      ) +
        forest(a, basis = z),
      dat,
      control = ctl()
    ),
    paste0("forest 1's ", refusal[[2L]]),
    fixed = TRUE
  )
}
expect_error(
  dbarts::dbarts(
    y ~ a + b + c,
    dat,
    forests = list(
      main = forest(n.perturb.cuts = 0),
      modifier = forest(a, basis = z)
    ),
    control = ctl()
  ),
  "forest \"main\"'s 'n.perturb.cuts'",
  fixed = TRUE
)
expect_error(
  build("hetero", ctl(), variance = list(proposal.probs = moving)),
  "varianceForest()'s 'proposal.probs' gives 'rule_gibbs' a positive share",
  fixed = TRUE
)
expect_error(
  build("hetero", ctl(), variance = list(n.perturb.cuts = -1)),
  "varianceForest()'s 'n.perturb.cuts'",
  fixed = TRUE
)
expect_error(varianceForest(treeShift = "always"), "unused argument")
leafPriorSampler <- build("two", ctl())
for (knob in names(all3)) {
  expect_error(
    leafPriorSampler$setLeafPrior(
      forests = list(do.call(forest, c(list(sd = 1), all3[knob])))
    ),
    paste0("'", knob, "' is fixed at creation"),
    fixed = TRUE
  )
}
rm(leafPriorSampler)
# a monotone forest's stated mixture is never rewritten, the package
# default's included
for (stated in list(moving, shipped)) {
  expect_error(
    dbarts::dbarts(
      y ~ a + b + c,
      dat,
      monotone = monotone(c(a = 1)),
      forests = list(forest(proposal.probs = stated)),
      control = ctl()
    ),
    "forest 1's 'proposal.probs' proposes other moves",
    fixed = TRUE
  )
}

# ---- a fit's control handed on, and the verbose summary ----

# the record rides the control; dbarts() drops it with the other records
stated <- build("two", ctl(), mean = list(proposal.probs = moving))
expect_false(is.null(attr(stated$control, "bartcore.moves")))
expect_null(attr(build("two", ctl())$control, "bartcore.moves"))
handedOn <- build("two", stated$control)
expect_null(attr(handedOn$control, "bartcore.moves"))
expect_identical(draws(handedOn), draws(build("two", ctl())))
rm(stated, handedOn)

printed <- capture.output(invisible(
  build(
    "two",
    ctl(verbose = TRUE),
    second = list(
      proposal.probs = moving,
      n.perturb.cuts = 4,
      treeShift = "never"
    )
  )$run(0L, 1L)
))
expect_true(any(grepl("^\\s*level fiber gibbs step: auto", printed)))
expect_true(any(grepl(
  "^\\s*proposal probabilities: .*; birth 0.50; perturb window 1$",
  printed
)))
expect_true(any(grepl(
  "forest 2: level fiber gibbs step: false",
  printed,
  fixed = TRUE
)))
expect_true(any(grepl(
  "forest 2: proposal probabilities: birth/death 0.40, swap 0.10, change 0.20, perturb 0.20, rule_gibbs 0.10; birth 0.45; perturb window 4",
  printed,
  fixed = TRUE
)))
# forest 1's own mode and window, and a forest differing from it only in its
# window
printed <- capture.output(invisible(
  build(
    "two",
    ctl(verbose = TRUE),
    mean = list(n.perturb.cuts = 5, treeShift = "always"),
    second = list(n.perturb.cuts = 4, treeShift = "always")
  )$run(0L, 1L)
))
expect_true(any(grepl("^\\s*level fiber gibbs step: true", printed)))
expect_true(any(grepl(
  "^\\s*proposal probabilities: .*; perturb window 5$",
  printed
)))
expect_true(any(grepl(
  "^\\s*forest 2: proposal probabilities: .*; perturb window 4$",
  printed
)))
