# The door a live sampler changes one forest through: $setForest and $getForest,
# $setVarianceForest and $getVarianceForest. The main oracle is the twin: a call
# before any run leaves the records and, on two chains, the next draws of a
# sampler created with the statement it leaves. Every case is named, an error in
# one is a failed expectation carrying that name, and the file runs to its end.

forest <- dbartsForests$forest
varianceForest <- dbartsForests$varianceForest
fixed <- dbartsPriors$fixed
student <- dbartsPriors$student
normal <- dbartsPriors$normal
invchi <- dbartsPriors$invchi
linear <- dbartsPriors$linear
gp <- dbartsPriors$gp

set.seed(11L)
n <- 60L
z <- rbinom(n, 1L, 0.5)
dose <- runif(n, 1, 5)
x <- data.frame(x1 = runif(n), x2 = runif(n))
y <- x$x1 + z * x$x2 + 0.2 * dose + rnorm(n, 0, 0.3)
yBinary <- as.double(y > stats::median(y))
category <- factor(sample(c("a", "b", "c"), n, TRUE))
frame <- cbind(x, z = factor(z), dose = dose, y = y)

# a mixture that proposes perturb moves, so that n.perturb.cuts moves a draw
shipped <- c(
  birth_death = 0.5,
  swap = 0,
  change = 0.3,
  perturb = 0.2,
  rule_gibbs = 0,
  birth = 0.5
)
moved <- c(
  birth_death = 0.4,
  swap = 0,
  change = 0.2,
  perturb = 0.4,
  rule_gibbs = 0,
  birth = 0.5
)
ctl <- function(...) {
  dbartsControl(
    n.chains = 2L,
    n.threads = 1L,
    n.samples = 4L,
    n.burn = 0L,
    seed = 5L,
    updateState = FALSE,
    verbose = FALSE,
    ...
  )
}
factorBasis <- list(basis = ~ factor(z), n.trees = 6L)
numericBasis <- list(basis = ~dose, n.trees = 6L)

# Each kind of forest the door's table separates: the forests of a sampler and
# the one the door is aimed at. `sharer` marks a forest whose sd is refused until
# a sharer can leave its group; `exempt` the names whose draws keep creation's
# start, which reads what creation was given; `drawn` a forest whose k is drawn.
kinds <- list(
  lone = list(forests = list(list()), target = 1L),
  lone.probit = list(
    forests = list(list()),
    target = 1L,
    y = yBinary,
    drawn = TRUE
  ),
  sharer = list(
    forests = list(mu1 = list(), mu2 = list(), tau = factorBasis),
    target = 2L,
    sharer = TRUE
  ),
  sharer.probit = list(
    forests = list(mu1 = list(), mu2 = list(), tau = factorBasis),
    target = 1L,
    sharer = TRUE,
    y = yBinary,
    drawn = TRUE
  ),
  no.basis.pair = list(
    forests = list(list(), list()),
    target = 1L,
    sharer = TRUE
  ),
  lone.sharer = list(
    forests = list(mu = list(), tau = factorBasis),
    target = 1L
  ),
  lone.sharer.probit = list(
    forests = list(mu = list(), tau = factorBasis),
    target = 1L,
    y = yBinary,
    drawn = TRUE
  ),
  stated = list(
    forests = list(mu = list(sd = 1.2), tau = factorBasis),
    target = 1L
  ),
  student.stated = list(
    forests = list(
      mu = list(sd = 1.2, amplitude = student(3)),
      tau = factorBasis
    ),
    target = 1L,
    exempt = "sd"
  ),
  normal.stated = list(
    forests = list(
      mu = list(sd = 1.2, amplitude = normal()),
      tau = factorBasis
    ),
    target = 1L,
    exempt = "sd"
  ),
  student.no.sd = list(
    forests = list(mu = list(amplitude = student(3)), tau = factorBasis),
    target = 1L,
    exempt = "sd"
  ),
  factor = list(forests = list(mu = list(), tau = factorBasis), target = 2L),
  factor.first = list(
    forests = list(tau = factorBasis, mu = list()),
    target = 1L
  ),
  factor.student = list(
    forests = list(
      mu = list(),
      tau = c(factorBasis, list(amplitude = student(2)))
    ),
    target = 2L
  ),
  factor.held = list(
    forests = list(
      mu = list(),
      tau = c(factorBasis, list(amplitude = fixed()))
    ),
    target = 2L
  ),
  numeric = list(forests = list(mu = list(), tau = numericBasis), target = 2L),
  numeric.student = list(
    forests = list(
      mu = list(),
      tau = c(numericBasis, list(amplitude = student(2)))
    ),
    target = 2L
  ),
  numeric.held = list(
    forests = list(
      mu = list(),
      tau = c(numericBasis, list(amplitude = fixed()))
    ),
    target = 2L
  ),
  three = list(
    forests = list(mu = list(), tau = factorBasis, gamma = numericBasis),
    target = 3L
  )
)

# A sampler of a kind with `stated` written into its target's statement, at
# position `target`, or NULL for the kind's.
build <- function(kind, stated = list(), target = NULL) {
  spec <- kinds[[kind]]
  entries <- spec$forests
  position <- if (is.null(target)) spec$target else target
  entries[[position]][names(stated)] <- stated
  forests <- lapply(entries, function(arguments) do.call(forest, arguments))
  response <- if (is.null(spec$y)) y else spec$y
  dbarts(
    x,
    response,
    forests = forests,
    n.trees = 8L,
    control = ctl(proposal.probs = shipped)
  )
}
buildMultinomial <- function(stated = list(), response = category) {
  dbarts(
    x,
    response,
    family = "multinomial",
    forests = list(do.call(forest, stated)),
    n.trees = 8L,
    control = ctl(proposal.probs = shipped)
  )
}
buildVariance <- function(mean = list(), variance = list()) {
  dbarts(
    x,
    y,
    forests = list(do.call(forest, mean)),
    variance = do.call(varianceForest, variance),
    n.trees = 8L,
    control = ctl(proposal.probs = shipped)
  )
}
draws <- function(sampler) {
  sampler$run(5L, 5L)[c("sigma", "train")]
}
# a model's records with what a formula holds of its frame left out
modelOf <- function(sampler) {
  model <- sampler$model
  for (index in seq_along(model@forests)) {
    model@forests[[index]]@basis.term <- NULL
  }
  model
}
door <- function(sampler, target, arguments) {
  do.call(sampler$setForest, c(list(target), arguments))
}
roundTrip <- function(sampler, target) {
  door(sampler, target, sampler$getForest(target))
}
reloaded <- function(sampler) {
  sampler$storeState()
  file <- tempfile(fileext = ".rds")
  saveRDS(sampler, file)
  loaded <- readRDS(file)
  unlink(file)
  loaded
}
refusal <- function(expr) {
  tryCatch(
    {
      expr
      "taken"
    },
    error = conditionMessage
  )
}
# a case's body runs here, and what stops it fails by the case's name
case <- function(name, body) {
  tryCatch(
    body,
    error = function(e) {
      expect_true(
        FALSE,
        info = paste0(name, " stopped: ", conditionMessage(e))
      )
    }
  )
  invisible(NULL)
}
expectFound <- function(text, pattern, info) {
  expect_true(
    grepl(pattern, text, fixed = TRUE),
    info = paste0(info, ": ", text)
  )
}
# the sampler refused, or left inert, reads, records and draws as `untouched`
expectAsUntouched <- function(sampler, untouched, target, info) {
  # a multinomial reader, left out, is the same with `target` NULL
  if (is.null(target) && length(sampler$model@forests) > 1L) {
    target <- 1L
  }
  expect_identical(
    sampler$getForest(target),
    untouched$getForest(target),
    info = paste(info, "reader")
  )
  expect_identical(
    modelOf(sampler),
    modelOf(untouched),
    info = paste(info, "records")
  )
  expect_identical(
    draws(sampler),
    draws(untouched),
    info = paste(info, "draws")
  )
}

# what each name is written to, each differing from every default it meets
written <- list(
  sd = 0.7,
  power = 0.5,
  base = 0.6,
  proposal.probs = moved,
  n.perturb.cuts = 5,
  treeShift = "always"
)
# whether a setting moves a draw on a forest here: a treeShift has a level step
# to act on only on a forest with a basis or a leaf to shift
expectMoves <- function(kind, name) {
  name != "treeShift"
}

## --- each name on each kind ----------------------------------------------------
for (kind in names(kinds)) {
  spec <- kinds[[kind]]
  for (name in names(written)) {
    if (name == "sd" && isTRUE(spec$sharer)) {
      next
    }
    info <- paste(kind, name)
    stated <- written[name]
    case(info, {
      # before any run: the twin
      live <- build(kind)
      expect_identical(live$control@n.chains, 2L, info = info)
      expect_null(door(live, spec$target, stated), info = info)
      twin <- build(kind, stated)
      expect_identical(
        modelOf(live),
        modelOf(twin),
        info = paste(info, "records")
      )
      liveDraws <- draws(live)
      twinDraws <- draws(twin)
      expect_identical(
        ncol(liveDraws$sigma),
        2L,
        info = paste(info, "two chains")
      )
      if (name %in% spec$exempt) {
        expect_false(
          identical(liveDraws, twinDraws),
          info = paste(info, "exempt: the start creation read")
        )
      } else {
        expect_identical(liveDraws, twinDraws, info = paste(info, "draws"))
      }
      if (expectMoves(kind, name)) {
        expect_false(
          identical(twinDraws, draws(build(kind))),
          info = paste(info, "moves a draw")
        )
      }

      # after a run: the reader and the engine, and the next draws of a twin
      # given the state, of a copy and of a reload
      live <- build(kind)
      invisible(live$run(0L, 2L))
      door(live, spec$target, stated)
      twin <- build(kind, stated)
      expect_identical(
        live$getForest(spec$target),
        twin$getForest(spec$target),
        info = paste(info, "reader")
      )
      if (name != "proposal.probs") {
        expect_identical(
          live$getForest(spec$target)[[name]],
          stated[[name]],
          info = paste(info, "reads what was written")
        )
      }
      if (name == "sd") {
        expect_identical(
          live$getLeafPrior(spec$target)$sd,
          twin$getLeafPrior(spec$target)$sd,
          info = paste(info, "getLeafPrior")
        )
      }
      live$storeState()
      twin$setState(live$state)
      expect_identical(
        live$getK(spec$target),
        twin$getK(spec$target),
        info = paste(info, "getK")
      )
      copied <- live$copy()
      loaded <- reloaded(live)
      expect_identical(
        modelOf(copied),
        modelOf(live),
        info = paste(info, "copy records")
      )
      expect_identical(
        modelOf(loaded),
        modelOf(live),
        info = paste(info, "reload records")
      )
      # an install of the state keeps it to its last digits, so what the state
      # restarts agrees bitwise and the sampler that goes on, to those digits
      after <- draws(twin)
      expect_identical(draws(copied), after, info = paste(info, "copy draws"))
      expect_identical(draws(loaded), after, info = paste(info, "reload draws"))
      expect_equal(
        draws(live),
        after,
        tolerance = 1e-10,
        info = paste(info, "goes on as the state twin")
      )
    })
  }
}

# A forest's own move settings and depth on every other forest of three: the
# others' records and reader do not move.
for (name in c("power", "base", "n.perturb.cuts")) {
  for (target in 1:3) {
    info <- paste("three, forest", target, name)
    case(info, {
      live <- build("three")
      door(live, target, written[name])
      twin <- build("three", written[name], target)
      expect_identical(modelOf(live), modelOf(twin), info = info)
      expect_identical(draws(live), draws(twin), info = info)
      untouched <- build("three")
      for (other in setdiff(1:3, target)) {
        expect_identical(
          live$getForest(other),
          untouched$getForest(other),
          info = paste(info, "leaves forest", other)
        )
      }
    })
  }
}

# Several names at once, in one call, as creation takes them.
case("several names", {
  stated <- list(sd = 0.7, power = 3.5, base = 0.6, n.perturb.cuts = 5)
  live <- build("factor")
  door(live, 2L, stated)
  twin <- build("factor", stated)
  expect_identical(modelOf(live), modelOf(twin), info = "several names")
  expect_identical(draws(live), draws(twin), info = "several names")
})

# The other leaf models and a stated k on a lone forest.
loneLeaves <- list(
  fixed.k = list(leaf.prior = normal(k = 3)),
  infinite.k = list(leaf.prior = normal(k = Inf)),
  linear = list(leaf.prior = linear(sd = 1, columns = "x2")),
  gp = list(leaf.prior = gp(sd = 1, columns = "x2"))
)
for (leaf in names(loneLeaves)) {
  info <- paste("lone leaf", leaf, "sd")
  case(info, {
    made <- function() {
      do.call(
        dbarts,
        c(
          list(x, y, n.trees = 8L, control = ctl(proposal.probs = shipped)),
          loneLeaves[[leaf]]
        )
      )
    }
    live <- made()
    live$setForest(1L, sd = 0.7)
    expect_identical(live$getLeafPrior()$sd, 0.7, info = info)
    expect_identical(live$getForest()$sd, 0.7, info = info)
    expect_identical(
      live$model@leaf.prior@prior.sd,
      0.7,
      info = paste(info, "records")
    )
    copied <- live$copy()
    liveDraws <- draws(live)
    expect_identical(draws(copied), liveDraws, info = paste(info, "copy"))
    expect_false(
      identical(draws(made()), liveDraws),
      info = paste(info, "moves a draw")
    )
  })
}

## --- the multinomial model's one statement ------------------------------------
multinomialWritten <- written[names(written) != "sd"]
for (name in names(multinomialWritten)) {
  info <- paste("multinomial", name)
  case(info, {
    stated <- multinomialWritten[name]
    live <- buildMultinomial()
    expect_null(do.call(live$setForest, stated), info = info)
    twin <- buildMultinomial(stated)
    expect_identical(
      modelOf(live),
      modelOf(twin),
      info = paste(info, "records")
    )
    expect_identical(
      live$getForest(),
      twin$getForest(),
      info = paste(info, "reader")
    )
    copied <- live$copy()
    liveDraws <- draws(live)
    expect_identical(liveDraws, draws(twin), info = paste(info, "draws"))
    expect_identical(draws(copied), liveDraws, info = paste(info, "copy"))
  })
}

## --- the reader ---------------------------------------------------------------
everyName <- c(
  "sd",
  "power",
  "base",
  "amplitude",
  "proposal.probs",
  "n.perturb.cuts",
  "treeShift"
)
expectPlainList <- function(read, expected, info) {
  expect_identical(class(read), "list", info = paste(info, "class"))
  expect_identical(
    names(attributes(read)),
    "names",
    info = paste(info, "attributes")
  )
  expect_identical(names(read), expected, info = paste(info, "names"))
}
for (kind in names(kinds)) {
  case(paste("reader", kind), {
    sampler <- build(kind)
    for (position in seq_along(kinds[[kind]]$forests)) {
      expectPlainList(
        sampler$getForest(position),
        everyName,
        paste("reader", kind, position)
      )
    }
  })
}
amplitudeRead <- list(
  lone = fixed(),
  sharer = fixed(),
  lone.sharer = fixed(),
  stated = fixed(),
  student.stated = student(3),
  normal.stated = normal(),
  student.no.sd = student(3),
  factor = normal(),
  factor.student = student(2),
  factor.held = fixed(),
  numeric = normal(),
  numeric.student = student(2),
  numeric.held = fixed()
)
for (kind in names(amplitudeRead)) {
  info <- paste("reader amplitude", kind)
  case(info, {
    expect_identical(
      build(kind)$getForest(kinds[[kind]]$target)$amplitude,
      amplitudeRead[[kind]],
      info = info
    )
  })
}
# no one number is in force where k is drawn: the name holds NULL, not NA or 0
for (kind in names(kinds)) {
  if (isTRUE(kinds[[kind]]$drawn)) {
    info <- paste("reader drawn k", kind)
    case(info, {
      read <- build(kind)$getForest(kinds[[kind]]$target)
      expect_true("sd" %in% names(read), info = paste(info, "name"))
      expect_true(is.null(read$sd), info = paste(info, "NULL"))
    })
  }
}
case("reader sd stated", {
  expect_identical(build("stated")$getForest(1L)$sd, 1.2, info = "no basis")
  expect_identical(
    build("factor", list(sd = 0.7))$getForest(2L)$sd,
    0.7,
    info = "factor"
  )
  expect_identical(
    build("numeric", list(sd = 0.7))$getForest(2L)$sd,
    0.7,
    info = "numeric reads no column name"
  )
  expect_identical(
    unname(build("numeric")$getForest(2L)$sd),
    unname(build("numeric")$getLeafPrior(2L)$sd),
    info = "numeric default is getLeafPrior's"
  )
  expect_null(
    names(build("numeric")$getForest(2L)$sd),
    info = "numeric default reads no column name"
  )
})
case("reader depth from the fitting function", {
  bartForm <- dbarts(
    x,
    y,
    power = 3,
    base = 0.7,
    control = ctl(proposal.probs = shipped)
  )
  sampler <- list(
    live = bartForm,
    copy = bartForm$copy(),
    new = new("dbartsSampler", bartForm$control, bartForm$model, bartForm$data),
    reload = reloaded(bartForm)
  )
  for (how in names(sampler)) {
    read <- sampler[[how]]$getForest()
    expect_identical(read$power, 3, info = paste(how, "power"))
    expect_identical(read$base, 0.7, info = paste(how, "base"))
    expect_identical(
      sampler[[how]]$model@forests[[1L]]@stated,
      bartForm$model@forests[[1L]]@stated,
      info = paste(how, "stated")
    )
  }
  # every forest of several takes the fitting function's depth at creation
  several <- dbarts(
    x,
    y,
    forests = list(forest(), forest(basis = ~ factor(z))),
    power = 3,
    control = ctl(proposal.probs = shipped)
  )
  expect_identical(several$getForest(1L)$power, 3, info = "forest 1 of two")
  expect_identical(several$getForest(2L)$power, 3, info = "forest 2 of two")
})
case("reader takes a control's value and a forest's own", {
  sampler <- dbarts(
    x,
    y,
    forests = list(forest(n.perturb.cuts = 3)),
    control = ctl(proposal.probs = shipped)
  )
  expect_identical(sampler$getForest()$n.perturb.cuts, 3, info = "its own")
  expect_identical(
    sampler$getForest()$proposal.probs,
    shipped,
    info = "the control's"
  )
  expect_identical(
    sampler$getForest()$treeShift,
    "auto",
    info = "the control's"
  )
})
case("reader forest argument", {
  sampler <- build("three")
  expectFound(
    refusal(sampler$getForest()),
    "'forest' must be given: this sampler has 3 forests",
    "left out"
  )
  expectFound(
    refusal(sampler$setForest(power = 3)),
    "'forest' must be given: this sampler has 3 forests",
    "left out, setter"
  )
  expectFound(refusal(sampler$getForest(4L)), "forest", "out of range")
  expectFound(
    refusal(sampler$getForest("nope")),
    "names no forest of this model",
    "unknown label"
  )
  expect_identical(
    sampler$getForest("gamma"),
    sampler$getForest(3L),
    info = "label is its position"
  )
  expect_identical(
    sampler$getForest("forest2"),
    sampler$getForest(2L),
    info = "forest<i> is position i"
  )
  lone <- build("lone")
  expect_identical(
    lone$getForest("forest1"),
    lone$getForest(),
    info = "one forest, forest1"
  )
  expect_identical(lone$getForest(1L), lone$getForest(), info = "one forest, 1")
})

## --- inert: a call that changes nothing has the draws and records of no call ----
inertKinds <- c(
  names(kinds),
  "k.infinite",
  "k.stated",
  "multinomial",
  "variance"
)
buildAny <- function(kind) {
  switch(
    kind,
    k.infinite = dbarts(
      x,
      yBinary,
      leaf.prior = normal(k = Inf),
      n.trees = 8L,
      control = ctl(proposal.probs = shipped)
    ),
    k.stated = dbarts(
      x,
      y,
      leaf.prior = normal(k = 3),
      n.trees = 8L,
      control = ctl(proposal.probs = shipped)
    ),
    multinomial = buildMultinomial(),
    variance = buildVariance(),
    build(kind)
  )
}
forestPositions <- function(kind) {
  if (kind %in% names(kinds)) seq_along(kinds[[kind]]$forests) else 1L
}
# a call on a kind's forest, the multinomial's left out and the others by position
doorAny <- function(sampler, kind, target, arguments) {
  if (kind == "multinomial") {
    do.call(sampler$setForest, arguments)
  } else {
    door(sampler, target, arguments)
  }
}
readerAny <- function(sampler, kind, target) {
  if (kind == "multinomial") sampler$getForest() else sampler$getForest(target)
}
for (kind in inertKinds) {
  case(paste("inert", kind), {
    sampler <- buildAny(kind)
    untouched <- buildAny(kind)
    before <- sampler$model
    for (target in forestPositions(kind)) {
      read <- readerAny(sampler, kind, target)
      # the list whole, handed back
      doorAny(sampler, kind, target, read)
      expect_identical(
        sampler$model,
        before,
        info = paste("round trip", kind, target)
      )
      for (name in names(read)) {
        info <- paste("inert", kind, target, name)
        # a subset restated, and a NULL on the name
        doorAny(sampler, kind, target, read[name])
        expect_identical(sampler$model, before, info = paste(info, "restated"))
        cleared <- list(NULL)
        names(cleared) <- name
        doorAny(sampler, kind, target, cleared)
        expect_identical(sampler$model, before, info = paste(info, "NULL"))
        expect_identical(
          readerAny(sampler, kind, target),
          read,
          info = paste(info, "reader")
        )
      }
    }
    expect_identical(
      draws(sampler),
      draws(untouched),
      info = paste("inert", kind, "draws")
    )
  })
}
# a NULL on a sd in force moves nothing on a forest that states one
case("inert NULL on a stated sd", {
  sampler <- build("stated")
  sampler$setForest(1L, sd = NULL)
  expect_identical(sampler$getForest(1L)$sd, 1.2, info = "reader")
  expect_identical(draws(sampler), draws(build("stated")), info = "draws")
})
case("inert on a run sampler with updateState", {
  sampler <- build("factor")
  invisible(sampler$run(0L, 2L))
  stored <- sampler$state
  sampler$setForest(
    2L,
    power = sampler$getForest(2L)$power,
    updateState = TRUE
  )
  expect_identical(sampler$state, stored, info = "no change stores nothing")
  sampler$setForest(2L, power = 0.5, updateState = TRUE)
  expect_false(
    identical(sampler$state, stored),
    info = "a change stores when told"
  )
  expectFound(
    refusal(sampler$setForest(2L, power = 3, updateState = "a")),
    "updateState",
    "updateState checked"
  )
})

# A value taken from above goes on following what it follows after a round trip.
setPerturbCuts <- function(sampler, value) {
  control <- sampler$control
  control@n.perturb.cuts <- value
  sampler$setControl(control)
}
for (kind in c("lone", "factor", "sharer")) {
  info <- paste("follows setControl after a round trip,", kind)
  case(info, {
    sampler <- build(kind)
    plain <- build(kind)
    target <- kinds[[kind]]$target
    roundTrip(sampler, target)
    setPerturbCuts(sampler, 7L)
    setPerturbCuts(plain, 7L)
    expect_identical(
      sampler$getForest(target)$n.perturb.cuts,
      7L,
      info = paste(info, "reader")
    )
    expect_identical(draws(sampler), draws(plain), info = info)
    # a forest that stated its own does not follow, the later write in force
    stating <- build(kind)
    door(stating, target, list(n.perturb.cuts = 5))
    setPerturbCuts(stating, 7L)
    expect_identical(
      stating$getForest(target)$n.perturb.cuts,
      5,
      info = paste(info, "a statement holds")
    )
  })
}
case("follows the fitting function's door after a round trip", {
  sampler <- build("lone")
  plain <- build("lone")
  roundTrip(sampler, 1L)
  sampler$setLeafPrior(normal(k = 3))
  plain$setLeafPrior(normal(k = 3))
  expect_identical(sampler$model, plain$model, info = "leaf prior reaches")
  expect_identical(draws(sampler), draws(plain), info = "leaf prior draws")
  sampler <- build("lone")
  plain <- build("lone")
  roundTrip(sampler, 1L)
  deeper <- function(s) {
    model <- s$model
    model@tree.prior@power <- 3
    s$setModel(model)
  }
  deeper(sampler)
  deeper(plain)
  expect_identical(sampler$getForest()$power, 3, info = "setModel reaches")
  expect_identical(draws(sampler), draws(plain), info = "setModel draws")
})
# what a door wrote goes back as it read, the mixture's canonical form included
for (kind in c("lone", "factor", "numeric", "normal.stated")) {
  info <- paste("inert after a write,", kind)
  case(info, {
    target <- kinds[[kind]]$target
    sampler <- build(kind)
    door(sampler, target, written[names(written) != "sd"])
    before <- sampler$model
    roundTrip(sampler, target)
    expect_identical(sampler$model, before, info = info)
  })
}
case("a statement left out of a round trip stays out", {
  for (kind in c("three", "lone.sharer", "factor")) {
    sampler <- build(kind)
    for (target in forestPositions(kind)) {
      stated <- sampler$model@forests[[target]]@stated
      roundTrip(sampler, target)
      expect_identical(
        sampler$model@forests[[target]]@stated,
        stated,
        info = paste(kind, target, "stated")
      )
    }
  }
})
case("a stated depth is stated, and a write states it once", {
  sampler <- build("three")
  sampler$setForest(2L, base = 0.6)
  sampler$setForest(2L, power = 0.5)
  expect_identical(
    sampler$model@forests[[2L]]@stated,
    build("three", list(power = 0.5, base = 0.6), 2L)$model@forests[[
      2L
    ]]@stated,
    info = "stated in creation's order"
  )
  sampler$setForest(2L, power = 3)
  expect_identical(
    sum(sampler$model@forests[[2L]]@stated == "power"),
    1L,
    info = "stated once"
  )
})
case("a numeric forest's default re-derives after a round trip", {
  sampler <- build("numeric")
  plain <- build("numeric")
  roundTrip(sampler, 2L)
  sampler$setForestBasis(2L, 3 * dose, updateBasisScale = TRUE)
  plain$setForestBasis(2L, 3 * dose, updateBasisScale = TRUE)
  expect_identical(
    sampler$getForest(2L)$sd,
    plain$getForest(2L)$sd,
    info = "the default"
  )
  expect_false(
    identical(sampler$getForest(2L)$sd, build("numeric")$getForest(2L)$sd),
    info = "the default moved"
  )
  expect_identical(draws(sampler), draws(plain), info = "draws")
  stating <- build("numeric")
  stating$setForest(2L, sd = 0.7)
  expectFound(
    refusal(stating$setForestBasis(2L, 3 * dose, updateBasisScale = TRUE)),
    "$setForest",
    "a stated sd does not re-derive"
  )
})

## --- the later write wins, on a lone forest, whichever door wrote first --------
agree <- function(sampler, info) {
  expect_identical(
    unname(sampler$getForest()$sd),
    {
      sd <- sampler$getLeafPrior()$sd
      if (is.na(sd)) NULL else unname(sd)
    },
    info = paste(info, "reader and engine")
  )
  expect_identical(
    sampler$copy()$getLeafPrior(),
    sampler$getLeafPrior(),
    info = paste(info, "a copy reads it")
  )
}
case("later write, sd after a leaf prior's k", {
  live <- dbarts(
    x,
    y,
    leaf.prior = normal(k = 3),
    n.trees = 8L,
    control = ctl(proposal.probs = shipped)
  )
  live$setForest(1L, sd = 0.7)
  twin <- build("lone", list(sd = 0.7))
  expect_identical(modelOf(live), modelOf(twin), info = "records")
  agree(live, "sd after k")
  expect_identical(draws(live), draws(twin), info = "draws")
})
case("later write, leaf prior after forest(sd = )", {
  live <- build("lone", list(sd = 2))
  live$setLeafPrior(normal(k = 3))
  twin <- dbarts(
    x,
    y,
    leaf.prior = normal(k = 3),
    n.trees = 8L,
    control = ctl(proposal.probs = shipped)
  )
  expect_identical(modelOf(live), modelOf(twin), info = "records")
  agree(live, "k after sd")
  expect_identical(draws(live), draws(twin), info = "draws")
})
case("later write, a probit forest(sd = ) returned to a drawn k", {
  live <- build("lone.probit", list(sd = 2))
  live$setLeafPrior(normal())
  expect_null(live$getForest()$sd, info = "the drawn k has no sd")
  expect_identical(
    live$model@leaf.prior,
    build("lone.probit")$model@leaf.prior,
    info = "the leaf prior of a drawn k"
  )
  # the drawn spread starts where the stated one left it, the one handoff
  expect_true(
    all(live$getK() == live$getLeafPrior()$k.scale / 2),
    info = "from the spread in force"
  )
  roundTrip(live, 1L)
  agree(live, "drawn after sd")
  # it is drawn: a run moves it, and a door's sd holds it again
  start <- live$getK()
  invisible(live$run(0L, 3L))
  expect_false(identical(live$getK(), start), info = "a run draws k")
  live$setForest(1L, sd = 0.7)
  expect_identical(live$getForest()$sd, 0.7, info = "held again")
})
case("later write, power and setModel", {
  live <- build("lone", list(power = 0.5))
  model <- live$model
  model@tree.prior@power <- 3
  live$setModel(model)
  expect_identical(
    live$getForest()$power,
    3,
    info = "setModel after a statement"
  )
  twin <- build("lone", list(power = 3))
  expect_identical(draws(live), draws(twin), info = "setModel draws")
  live <- dbarts(
    x,
    y,
    power = 3,
    n.trees = 8L,
    control = ctl(proposal.probs = shipped)
  )
  live$setForest(1L, power = 0.5)
  twin <- build("lone", list(power = 0.5))
  expect_identical(live$getForest()$power, 0.5, info = "a door after the fit's")
  expect_identical(
    live$model@tree.prior@power,
    0.5,
    info = "the record the fit's wrote"
  )
  expect_identical(draws(live), draws(twin), info = "a door's draws")
})
case("setLeafPrior takes one argument less", {
  live <- build("lone")
  expectFound(refusal(live$setLeafPrior(NULL)), "leaf", "a NULL leaf prior")
  expectFound(
    refusal(live$setLeafPrior(forests = list(1))),
    "unused argument",
    "forests = gone"
  )
  several <- build("factor")
  expectFound(
    refusal(several$setLeafPrior(normal(k = 3), forests = list(1))),
    "unused argument",
    "forests = gone on several"
  )
})

## --- what creation refuses of a statement, the door refuses for the same reason --
# Each row: the statement, the kind it is made on and the words both doors use.
refusals <- list(
  list("power negative", "factor", list(power = -1), "power"),
  list("power long", "factor", list(power = c(1, 2)), "power"),
  list("power NA", "factor", list(power = NA_real_), "power"),
  list("power text", "lone", list(power = "a"), "power"),
  list("base zero", "factor", list(base = 0), "base"),
  list("base above 1", "factor", list(base = 1.5), "base"),
  list("base long", "lone", list(base = c(0.5, 0.6)), "base"),
  list("base NA", "lone", list(base = NA_real_), "base"),
  list("sd negative", "factor", list(sd = -1), "sd"),
  list("sd zero", "numeric", list(sd = 0), "sd"),
  list("sd infinite", "factor", list(sd = Inf), "sd"),
  list("sd NA", "stated", list(sd = NA_real_), "sd"),
  list("sd long", "factor", list(sd = c(1, 2)), "sd"),
  list("sd text", "lone", list(sd = "a"), "sd"),
  list("sd invchi", "factor", list(sd = invchi(1, 1)), "sd"),
  list("sd student", "factor", list(sd = student(1)), "sd"),
  list("cuts zero", "factor", list(n.perturb.cuts = 0L), "n.perturb.cuts"),
  list("cuts part", "lone", list(n.perturb.cuts = 1.5), "n.perturb.cuts"),
  list("cuts text", "lone", list(n.perturb.cuts = "a"), "n.perturb.cuts"),
  list("shift word", "factor", list(treeShift = "sometimes"), "treeShift"),
  list("shift number", "lone", list(treeShift = 2), "treeShift"),
  list(
    "probs short",
    "factor",
    list(proposal.probs = c(1, 2)),
    "proposal.probs"
  ),
  list(
    "probs negative",
    "lone",
    list(proposal.probs = c(birth_death = -1, change = 2)),
    "proposal.probs"
  ),
  list(
    "probs unknown",
    "lone",
    list(proposal.probs = c(birth_death = 0.5, grow = 0.5)),
    "proposal.probs"
  ),
  list("amplitude word", "factor", list(amplitude = "a"), "amplitude"),
  list("amplitude number", "factor", list(amplitude = 2), "amplitude"),
  list(
    "student on a lone forest",
    "lone",
    list(amplitude = student(1)),
    "no coefficient to give a prior"
  ),
  list(
    "normal on a lone forest",
    "lone",
    list(amplitude = normal()),
    "no coefficient to give a prior"
  ),
  list(
    "student on a sharer",
    "sharer",
    list(amplitude = student(1)),
    "share one total"
  ),
  list(
    "normal on a sharer",
    "sharer",
    list(amplitude = normal()),
    "share one total"
  )
)
for (row in refusals) {
  info <- paste("two doors,", row[[1L]])
  case(info, {
    kind <- row[[2L]]
    target <- kinds[[kind]]$target
    atCreation <- refusal(suppressWarnings(build(kind, row[[3L]])))
    live <- build(kind)
    atDoor <- refusal(door(live, target, row[[3L]]))
    expect_false(identical(atCreation, "taken"), info = paste(info, "creation"))
    expect_false(identical(atDoor, "taken"), info = paste(info, "door"))
    expectFound(atDoor, row[[4L]], paste(info, "the reason"))
    expect_identical(atDoor, atCreation, info = paste(info, "the same words"))
    expectAsUntouched(live, build(kind), target, info)
  })
}
# refused at one door only, each for its stated reason: a live sampler keeps its
# coefficient prior and a sharer its part, until later plans
doorOnly <- list(
  list(
    "a changed coefficient on a basis forest",
    "factor",
    list(amplitude = student(2)),
    "a forest's coefficient prior cannot be changed live yet: forest \"tau\" keeps normal()"
  ),
  list(
    "normal() on a student forest",
    "factor.student",
    list(amplitude = normal()),
    "keeps student(2)"
  ),
  list(
    "another df",
    "numeric.student",
    list(amplitude = student(3)),
    "keeps student(2)"
  ),
  list(
    "drawn on a held forest",
    "factor.held",
    list(amplitude = normal()),
    "keeps fixed()"
  ),
  list(
    "held on a drawn forest",
    "factor",
    list(amplitude = fixed()),
    "keeps normal()"
  ),
  list(
    "drawn on a stated forest",
    "stated",
    list(amplitude = normal()),
    "keeps fixed()"
  ),
  list(
    "held on a drawn forest with no basis",
    "normal.stated",
    list(amplitude = fixed()),
    "keeps normal()"
  ),
  list(
    "a sharer's size",
    "sharer",
    list(sd = 0.7),
    "this forest shares one prior total with forests 1, 2; a forest's own size cannot be changed live yet"
  ),
  list(
    "a sharer's size where k is drawn",
    "sharer.probit",
    list(sd = 0.7),
    "shares one prior total"
  )
)
for (row in doorOnly) {
  info <- paste("door only,", row[[1L]])
  case(info, {
    kind <- row[[2L]]
    target <- kinds[[kind]]$target
    live <- build(kind)
    expectFound(refusal(door(live, target, row[[3L]])), row[[4L]], info)
    expectAsUntouched(live, build(kind), target, info)
  })
}
# a statement needs what it names beside it, merged and not name by name: student
# with a stated sd is not refused for the sd it lacks, which it has
case("merged statement", {
  stating <- build("stated")
  expect_identical(stating$getForest(1L)$sd, 1.2, info = "the sd in force")
  expectFound(
    refusal(stating$setForest(1L, amplitude = student(1))),
    "cannot be changed live yet",
    "student beside a stated sd reads the coefficient refusal"
  )
  sharing <- build("sharer")
  expectFound(
    refusal(sharing$setForest(1L, sd = 1, amplitude = student(1))),
    "cannot be changed live yet",
    "student stated with its sd reads the coefficient refusal"
  )
  expectFound(
    refusal(sharing$setForest(1L, amplitude = student(1))),
    "share one total",
    "student with no sd reads the sharers'"
  )
})
case("monotone mixture, at both doors", {
  made <- function(...) {
    dbarts(
      x,
      y,
      monotone = monotone(c(x1 = 1)),
      forests = list(forest(...)),
      n.trees = 8L,
      control = ctl()
    )
  }
  wording <- "forest 1's 'proposal.probs' proposes other moves"
  expectFound(refusal(made(proposal.probs = moved)), wording, "creation")
  live <- made()
  expectFound(
    refusal(live$setForest(1L, proposal.probs = moved)),
    wording,
    "door"
  )
  roundTrip(live, 1L)
  expect_identical(live$model, made()$model, info = "monotone round trip")
  live$setForest(1L, power = 0.5)
  expect_identical(
    live$getForest(1L)$power,
    0.5,
    info = "monotone takes a depth"
  )
})

# Structure by name, whatever the value, the value in force included.
structureRows <- list(
  list("n.trees", 8L, "a forest's tree count is part of the model"),
  list("n.trees", 3L, "a forest's tree count is part of the model"),
  list("vars", "x1", "the columns a forest may split on are part of the model"),
  list(
    "vars",
    c("x1", "x2"),
    "the columns a forest may split on are part of the model"
  ),
  list(
    "interactions",
    list(c("x1", "x2")),
    "a forest's interaction limits are part of the model"
  ),
  list("blocks", list("x1"), "a forest's blocks are part of the model")
)
for (kind in c("lone", "sharer", "factor", "numeric", "multinomial")) {
  for (row in structureRows) {
    info <- paste("structure,", kind, row[[1L]], deparse(row[[2L]])[1L])
    case(info, {
      live <- buildAny(kind)
      arguments <- list(row[[2L]])
      names(arguments) <- row[[1L]]
      arguments$power <- 0.5
      message <- refusal(doorAny(live, kind, 1L, arguments))
      expectFound(
        message,
        paste0("'", row[[1L]], "' is fixed at creation"),
        info
      )
      expectFound(message, row[[3L]], info)
      expectAsUntouched(live, buildAny(kind), NULL, info)
    })
  }
}
case("basis and an unknown name", {
  live <- build("factor")
  expectFound(
    refusal(live$setForest(2L, basis = ~ factor(z), power = 0.5)),
    "'basis' is data, not a setting of the forest; change it with $setForestBasis",
    "a basis names its door"
  )
  expectFound(
    refusal(live$setForest(2L, nonsense = 1, power = 0.5)),
    "unused argument 'nonsense': $setForest takes forest, sd, power",
    "an unknown name"
  )
  expectFound(
    refusal(live$setForest(2L, 1, 2, 3, 4, 5, 6, 7, NULL, 9)),
    "unused argument",
    "one too many, unnamed"
  )
  expectAsUntouched(live, build("factor"), 2L, "basis and an unknown name")
})

# A refused call leaves everything as it was: the bad name comes after and before
# a good one, and the good one is not written.
atomicRows <- list(
  list("factor", list(power = 0.5, base = -1)),
  list("factor", list(base = 0.6, power = -1)),
  list("factor", list(sd = 0.7, power = -1)),
  list("factor", list(power = 0.5, sd = -1)),
  list("factor", list(n.perturb.cuts = 5, sd = -1)),
  list("factor", list(proposal.probs = moved, treeShift = "sometimes")),
  list("factor", list(treeShift = "always", proposal.probs = c(1, 2))),
  list("factor", list(power = 0.5, amplitude = student(2))),
  list("factor", list(sd = 0.7, amplitude = fixed())),
  list("factor", list(n.perturb.cuts = 5, amplitude = student(2))),
  list("sharer", list(power = 0.5, sd = 0.7)),
  list("sharer", list(n.perturb.cuts = 5, sd = 0.7)),
  list("lone", list(power = 0.5, sd = -1)),
  list("lone", list(power = 0.5, n.perturb.cuts = 5, amplitude = normal())),
  list("lone", list(power = 0.5, sd = 1e308))
)
for (row in atomicRows) {
  kind <- row[[1L]]
  info <- paste("atomic,", kind, paste(names(row[[2L]]), collapse = " "))
  case(info, {
    live <- build(kind)
    target <- kinds[[kind]]$target
    expect_false(
      identical(refusal(door(live, target, row[[2L]])), "taken"),
      info = paste(info, "refused")
    )
    expectAsUntouched(live, build(kind), target, info)
  })
}
# The checks come first, and what a write then refuses undoes what went before
# it: here the engine's own refusal, which a check ahead of it answers on every
# statement a sampler can now reach.
refusingWrite <- function(sampler, forest, ...) {
  writer <- dbarts:::writeForestDoor
  environment(writer) <- list2env(
    list(
      .Call = function(symbol, ...) {
        if (identical(symbol, dbarts:::C_dbarts_bartcore_setForestSd)) {
          stop("the engine refused")
        }
        base::.Call(symbol, ...)
      }
    ),
    parent = environment(dbarts:::writeForestDoor)
  )
  writer(sampler, forest, list(...), character(), NULL)
}
case("a write refused after writes", {
  live <- build("factor")
  message <- refusal(refusingWrite(
    live,
    2L,
    sd = 0.7,
    power = 0.5,
    base = 0.6,
    n.perturb.cuts = 5,
    proposal.probs = moved
  ))
  expectFound(message, "the engine refused", "it refused")
  expectAsUntouched(live, build("factor"), 2L, "the write undone")
})

## --- the variance forest ---------------------------------------------------------
varianceWritten <- list(
  power = 0.5,
  base = 0.6,
  proposal.probs = c(
    birth_death = 0.4,
    swap = 0.1,
    change = 0.2,
    perturb = 0.3,
    birth = 0.5
  ),
  n.perturb.cuts = 5
)
for (name in names(varianceWritten)) {
  info <- paste("variance door,", name)
  case(info, {
    stated <- varianceWritten[name]
    live <- buildVariance()
    expect_null(do.call(live$setVarianceForest, stated), info = info)
    twin <- buildVariance(variance = stated)
    expect_identical(
      modelOf(live),
      modelOf(twin),
      info = paste(info, "records")
    )
    expect_identical(
      live$getVarianceForest(),
      twin$getVarianceForest(),
      info = paste(info, "reader")
    )
    # the mean forest does not move
    expect_identical(
      live$getForest(),
      buildVariance()$getForest(),
      info = paste(info, "leaves the mean forest")
    )
    liveDraws <- draws(live)
    expect_identical(liveDraws, draws(twin), info = paste(info, "draws"))
    expect_false(
      identical(liveDraws, draws(buildVariance())),
      info = paste(info, "moves a draw")
    )
    live$storeState()
    twin$setState(live$state)
    copied <- live$copy()
    loaded <- reloaded(live)
    expect_identical(
      modelOf(copied),
      modelOf(live),
      info = paste(info, "copy records")
    )
    expect_identical(
      modelOf(loaded),
      modelOf(live),
      info = paste(info, "reload records")
    )
    after <- draws(twin)
    expect_identical(draws(copied), after, info = paste(info, "copy"))
    expect_identical(draws(loaded), after, info = paste(info, "reload"))
  })
}
case("the variance reader", {
  sampler <- buildVariance()
  read <- sampler$getVarianceForest()
  expectPlainList(
    read,
    c("power", "base", "proposal.probs", "n.perturb.cuts"),
    "variance reader"
  )
  before <- sampler$model
  do.call(sampler$setVarianceForest, read)
  expect_identical(sampler$model, before, info = "variance round trip")
  sampler$setVarianceForest(power = NULL)
  expect_identical(sampler$model, before, info = "variance NULL")
  expect_identical(
    draws(sampler),
    draws(buildVariance()),
    info = "variance inert draws"
  )
  expect_null(
    build("lone")$getVarianceForest(),
    info = "no variance forest reads"
  )
})
case("the mean door leaves the variance forest", {
  live <- buildVariance(variance = list(power = 2))
  live$setForest(1L, power = 0.5)
  twin <- buildVariance(mean = list(power = 0.5), variance = list(power = 2))
  expect_identical(modelOf(live), modelOf(twin), info = "records")
  expect_identical(draws(live), draws(twin), info = "draws")
  expect_identical(
    live$getVarianceForest()$power,
    2,
    info = "reader"
  )
  # it takes the mean forest's depth at creation, and is on its own after
  live <- buildVariance()
  live$setForest(1L, base = 0.6)
  twin <- buildVariance(mean = list(base = 0.6), variance = list(base = 0.95))
  expect_identical(draws(live), draws(twin), info = "the default depth stays")
})
case("setModel leaves the variance forest", {
  live <- buildVariance(variance = list(power = 2))
  model <- live$model
  model@tree.prior@power <- 0.5
  live$setModel(model)
  expect_identical(live$getVarianceForest()$power, 2, info = "reader")
  twin <- buildVariance(mean = list(power = 0.5), variance = list(power = 2))
  expect_identical(draws(live), draws(twin), info = "draws")
})
case("a sampler with no variance forest", {
  live <- build("lone")
  expectFound(
    refusal(live$setVarianceForest(power = 3)),
    "this sampler has no variance forest: whether there is one is part of the model; make a new sampler",
    "refuses"
  )
  expectFound(
    refusal(live$setVarianceForest()),
    "this sampler has no variance forest",
    "refuses with no name"
  )
})
varianceRefusals <- list(
  list("n.trees", list(n.trees = 8L), "'n.trees' is fixed at creation"),
  list("vars", list(vars = "x1"), "'vars' is fixed at creation"),
  list("sd", list(sd = 1), "unused argument 'sd'"),
  list("amplitude", list(amplitude = fixed()), "unused argument 'amplitude'"),
  list("treeShift", list(treeShift = "auto"), "unused argument 'treeShift'"),
  list("basis", list(basis = 1), "unused argument 'basis'"),
  list(
    "interactions",
    list(interactions = list()),
    "unused argument 'interactions'"
  ),
  list("power", list(power = -1), "forest 'power' must be"),
  list("base with power", list(power = 0.5, base = 2), "forest 'base' must be"),
  list("rule_gibbs", list(proposal.probs = shipped + c(0, 0, 0, 0, 0.1, 0)), "")
)
for (row in varianceRefusals) {
  info <- paste("variance refuses,", row[[1L]])
  case(info, {
    live <- buildVariance()
    message <- refusal(do.call(live$setVarianceForest, row[[2L]]))
    expect_false(identical(message, "taken"), info = info)
    expectFound(message, row[[3L]], info)
    expectAsUntouched(live, buildVariance(), NULL, info)
    expect_identical(
      live$getVarianceForest(),
      buildVariance()$getVarianceForest(),
      info = paste(info, "variance reader")
    )
  })
}

## --- labels and the multinomial model --------------------------------------------
case("a label, where a formula and a list put it elsewhere", {
  byList <- function(...) {
    dbarts(
      x,
      y,
      forests = list(...),
      n.trees = 8L,
      control = ctl(proposal.probs = shipped)
    )
  }
  mu <- forest()
  tau <- forest(basis = ~ factor(z), n.trees = 6L)
  first <- byList(mu = mu, tau = tau)
  second <- byList(tau = tau, mu = mu)
  first$setForest("tau", power = 0.5)
  second$setForest("tau", power = 0.5)
  expect_identical(first$getForest(2L)$power, 0.5, info = "position 2")
  expect_identical(second$getForest(1L)$power, 0.5, info = "position 1")
  expect_false(
    identical(second$getForest(2L)$power, 0.5),
    info = "the other forest stays"
  )
  twin <- dbarts(
    x,
    y,
    forests = list(
      tau = forest(basis = ~ factor(z), n.trees = 6L, power = 0.5),
      mu = mu
    ),
    n.trees = 8L,
    control = ctl(proposal.probs = shipped)
  )
  expect_identical(modelOf(second), modelOf(twin), info = "forest 1's record")
  expect_identical(draws(second), draws(twin), info = "forest 1's draws")
  byFormula <- dbarts(
    y ~ forest(x1) + forest(x2, basis = z),
    frame,
    n.trees = 8L,
    control = ctl(proposal.probs = shipped)
  )
  byFormula$setForest("z", power = 0.5)
  expect_identical(
    byFormula$getForest(2L)$power,
    0.5,
    info = "a formula's label"
  )
  expect_identical(
    byFormula$getForest("z"),
    byFormula$getForest(2L),
    info = "read"
  )
  byFormula$setForest("forest1", base = 0.6)
  expect_identical(byFormula$getForest(1L)$base, 0.6, info = "forest<i>")
})
case("the multinomial door takes the forest left out", {
  live <- buildMultinomial()
  read <- live$getForest()
  expectPlainList(
    read,
    c("power", "base", "proposal.probs", "n.perturb.cuts", "treeShift"),
    "multinomial reader"
  )
  for (which in list(1L, 2L, 3L, "a", "b", "c", "forest2")) {
    expect_identical(
      live$getForest(which),
      read,
      info = paste("multinomial, one list", which)
    )
  }
  for (which in list(1L, 2L, "b", "forest2")) {
    message <- refusal(live$setForest(which, power = 4))
    expectFound(
      message,
      "a multinomial model's forests take one statement, so this would change every category; write s$setForest(power = 4)",
      paste("multinomial refuses", which)
    )
  }
  expectFound(
    refusal(live$getForest(4L)),
    "forest",
    "multinomial, out of range"
  )
  expectAsUntouched(live, buildMultinomial(), NULL, "refused setters")
  # the depth reaches every category, from any reader
  live$setForest(power = 0.5)
  for (which in list(1L, 2L, "c")) {
    expect_identical(
      live$getForest(which)$power,
      0.5,
      info = paste("multinomial, written", which)
    )
  }
})
case("the multinomial statement refuses what it has none of", {
  for (stated in list(list(sd = 1), list(amplitude = fixed()))) {
    name <- names(stated)
    live <- buildMultinomial()
    message <- refusal(do.call(live$setForest, stated))
    expect_false(identical(message, "taken"), info = paste("door", name))
    expectFound(message, "multinomial", paste("door reason", name))
    expectAsUntouched(live, buildMultinomial(), NULL, paste("refused", name))
  }
  expectFound(
    refusal(buildMultinomial(list(sd = 1))),
    "multinomial",
    "creation refuses sd"
  )
})
case("a name on the multinomial entry", {
  expectFound(
    refusal(dbarts(
      x,
      category,
      family = "multinomial",
      forests = list(mu = forest()),
      n.trees = 8L,
      control = ctl()
    )),
    "a multinomial model's one forest entry states every category forest and takes no name",
    "refused at creation"
  )
})
case("an entry not in a list states every category too", {
  sampler <- dbarts(
    x,
    category,
    family = "multinomial",
    forests = forest(power = 0.5),
    n.trees = 8L,
    control = ctl(proposal.probs = shipped)
  )
  expect_identical(sampler$getForest("b")$power, 0.5, info = "one statement")
})
case("a level read by every forest-indexed reader", {
  sampler <- buildMultinomial()
  invisible(sampler$run(0L, 2L))
  readers <- list(
    getLeafPrior = sampler$getLeafPrior,
    getK = sampler$getK,
    getForestFits = sampler$getForestFits,
    getForest = sampler$getForest
  )
  for (reader in names(readers)) {
    for (position in 1:3) {
      level <- c("a", "b", "c")[position]
      expect_identical(
        readers[[reader]](level),
        readers[[reader]](position),
        info = paste(reader, level)
      )
      expect_identical(
        readers[[reader]](paste0("forest", position)),
        readers[[reader]](position),
        info = paste(reader, "forest<i>", position)
      )
    }
  }
  expect_identical(
    sampler$getTrees(forest = "b"),
    sampler$getTrees(forest = 2L),
    info = "getTrees"
  )
  expect_identical(
    names(sampler$getLeafPrior()),
    c("a", "b", "c"),
    info = "a NULL read names its list by levels"
  )
  colliding <- factor(sample(c("forest2", "b", "c"), n, TRUE))
  # a level spelled forest<i> names position i, or is refused where it is another
  spelled <- buildMultinomial(response = colliding)
  levelsThere <- levels(colliding)
  expect_identical(
    levelsThere,
    c("b", "c", "forest2"),
    info = "the fixture puts forest2 last"
  )
  expectFound(
    refusal(spelled$getForest("forest2")),
    "forest2",
    "a level at another position's name is refused by the conflict"
  )
})

## --- amplitude = normal() --------------------------------------------------------
for (kind in c("factor", "numeric")) {
  info <- paste("amplitude normal,", kind)
  case(info, {
    statedNormal <- build(kind, list(amplitude = normal()))
    left <- build(kind)
    expect_identical(
      modelOf(statedNormal),
      modelOf(left),
      info = paste(info, "records")
    )
    expect_identical(
      draws(statedNormal),
      draws(left),
      info = paste(info, "draws")
    )
    live <- build(kind)
    before <- live$model
    live$setForest(2L, amplitude = normal())
    expect_identical(live$model, before, info = paste(info, "door, in force"))
    expect_identical(
      draws(live),
      draws(build(kind)),
      info = paste(info, "door draws")
    )
  })
}
case("amplitude in force by its spelling", {
  live <- build("factor.student")
  before <- live$model
  live$setForest(2L, amplitude = student(2))
  live$setForest(2L, amplitude = student(2L))
  expect_identical(live$model, before, info = "student, the same df")
  live <- build("factor.held")
  before <- live$model
  live$setForest(2L, amplitude = fixed())
  expect_identical(live$model, before, info = "held")
  expectFound(
    refusal(live$setForest(2L, amplitude = normal())),
    "keeps fixed()",
    "held, normal refused"
  )
  live <- build("normal.stated")
  before <- live$model
  live$setForest(1L, amplitude = normal())
  expect_identical(live$model, before, info = "a no-basis forest's normal")
  expect_identical(
    live$getForest(1L)$amplitude,
    normal(),
    info = "a no-basis forest's normal reads normal(), not a t"
  )
})
