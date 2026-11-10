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
# the one the door is aimed at. `sharer` marks a forest whose sd takes it out of
# its group, which no creation reaches (test-forest-leave.R); `exempt` the names whose draws keep creation's
# start, which reads what creation was given; `drawn` a forest whose k is drawn;
# `sigest` a starting sigma stated where the estimate, which reads the bases,
# leaves a setting's ten sweeps here drawing as the default's.
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
  three = list(
    forests = list(mu = list(), tau = factorBasis, gamma = numericBasis),
    target = 3L,
    sigest = 0.4
  ),
  bases.only = list(
    forests = list(tau = factorBasis, gamma = numericBasis),
    target = 1L
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
    n.trees = if (all(vapply(entries, function(e) !is.null(e$basis), NA))) {
      NULL
    } else {
      8L
    },
    sigest = spec$sigest,
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
# a case's body runs here, and what stops it fails by the case's name; `same`
# labels an expectation by the case and what it pins
caseName <- ""
same <- function(current, target, ...) {
  expect_identical(current, target, info = paste(caseName, ...))
}
case <- function(name, body) {
  caseName <<- name
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
      same(live$control@n.chains, 2L)
      expect_null(door(live, spec$target, stated), info = info)
      twin <- build(kind, stated)
      same(modelOf(live), modelOf(twin), "records")
      liveDraws <- draws(live)
      twinDraws <- draws(twin)
      same(ncol(liveDraws$sigma), 2L, "two chains")
      if (name %in% spec$exempt) {
        expect_false(
          identical(liveDraws, twinDraws),
          info = paste(info, "exempt: the start creation read")
        )
      } else {
        same(liveDraws, twinDraws, "draws")
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
      same(live$getForest(spec$target), twin$getForest(spec$target), "reader")
      if (name != "proposal.probs") {
        same(
          live$getForest(spec$target)[[name]],
          stated[[name]],
          "reads what was written"
        )
      }
      if (name == "sd") {
        same(
          live$getLeafPrior(spec$target)$sd,
          twin$getLeafPrior(spec$target)$sd,
          "getLeafPrior"
        )
      }
      live$storeState()
      twin$setState(live$state)
      copied <- live$copy()
      loaded <- reloaded(live)
      # an install of the state keeps it to its last digits, so what the state
      # restarts agrees bitwise and the sampler that goes on, to those digits
      after <- draws(twin)
      same(draws(copied), after, "copy draws")
      same(draws(loaded), after, "reload draws")
      expect_equal(
        draws(live),
        after,
        tolerance = 1e-10,
        info = paste(info, "goes on as the state twin")
      )
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
  fixed.k = list(
    list(leaf.prior = normal(k = 3)),
    list(forests = list(forest(sd = 0.7)))
  ),
  infinite.k = list(
    list(leaf.prior = normal(k = Inf)),
    list(forests = list(forest(sd = 0.7)))
  ),
  linear = list(
    list(leaf.prior = linear(sd = 1, columns = "x2")),
    list(leaf.prior = linear(sd = 0.7, columns = "x2"))
  ),
  gp = list(
    list(leaf.prior = gp(sd = 1, columns = "x2")),
    list(leaf.prior = gp(sd = 0.7, columns = "x2"))
  )
)
for (leaf in names(loneLeaves)) {
  info <- paste("lone leaf", leaf, "sd")
  case(info, {
    made <- function(stating) {
      do.call(
        dbarts,
        c(
          list(x, y, n.trees = 8L, control = ctl(proposal.probs = shipped)),
          stating
        )
      )
    }
    live <- made(loneLeaves[[leaf]][[1L]])
    live$setForest(1L, sd = 0.7)
    twin <- made(loneLeaves[[leaf]][[2L]])
    same(live$model, twin$model, "twin records")
    same(draws(live), draws(twin), "twin draws")
    live <- made(loneLeaves[[leaf]][[1L]])
    live$setForest(1L, sd = 0.7)
    same(live$getLeafPrior()$sd, 0.7)
    same(live$getForest()$sd, 0.7)
    same(live$model@leaf.prior@prior.sd, 0.7, "records")
    copied <- live$copy()
    liveDraws <- draws(live)
    same(draws(copied), liveDraws, "copy")
    expect_false(
      identical(draws(made(loneLeaves[[leaf]][[1L]])), liveDraws),
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
    same(modelOf(live), modelOf(twin), "records")
    same(live$getForest(), twin$getForest(), "reader")
    copied <- live$copy()
    liveDraws <- draws(live)
    same(liveDraws, draws(twin), "draws")
    same(draws(copied), liveDraws, "copy")
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
  numeric.student = student(2)
)
for (kind in names(amplitudeRead)) {
  info <- paste("reader amplitude", kind)
  case(info, {
    same(
      build(kind)$getForest(kinds[[kind]]$target)$amplitude,
      amplitudeRead[[kind]],
      ""
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
  same(build("stated")$getForest(1L)$sd, 1.2, "no basis")
  same(build("factor", list(sd = 0.7))$getForest(2L)$sd, 0.7, "factor")
  same(
    build("numeric", list(sd = 0.7))$getForest(2L)$sd,
    0.7,
    "numeric reads no column name"
  )
  same(
    unname(build("numeric")$getForest(2L)$sd),
    unname(build("numeric")$getLeafPrior(2L)$sd),
    "numeric default is getLeafPrior's"
  )
  expect_null(
    names(build("numeric")$getForest(2L)$sd),
    info = "reader sd stated: numeric default reads no column name"
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
  }
  # beside a forest with a basis, it is forest 1's and the other keeps its kind's
  several <- dbarts(
    x,
    y,
    forests = list(forest(), forest(basis = ~ factor(z))),
    power = 2.5,
    control = ctl(proposal.probs = shipped)
  )
  same(several$getForest(1L)$power, 2.5, "forest 1 of two")
  same(several$getForest(2L)$power, 3, "forest 2 of two")
})
case("reader takes a control's value and a forest's own", {
  sampler <- dbarts(
    x,
    y,
    forests = list(forest(n.perturb.cuts = 3)),
    control = ctl(proposal.probs = shipped)
  )
  same(sampler$getForest()$n.perturb.cuts, 3, "its own")
  same(sampler$getForest()$proposal.probs, shipped, "the control's")
  same(sampler$getForest()$treeShift, "auto", "the control's")
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
  expectFound(
    refusal(sampler$getForest(4L)),
    "forest index out of range",
    "out of range"
  )
  expectFound(
    refusal(sampler$getForest("nope")),
    "names no forest of this model",
    "unknown label"
  )
  same(
    sampler$getForest("gamma"),
    sampler$getForest(3L),
    "label is its position"
  )
  same(
    sampler$getForest("forest2"),
    sampler$getForest(2L),
    "forest<i> is position i"
  )
  lone <- build("lone")
  same(lone$getForest("forest1"), lone$getForest(), "one forest, forest1")
  same(lone$getForest(1L), lone$getForest(), "one forest, 1")
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
        # a subset restated, and a NULL on the name
        doorAny(sampler, kind, target, read[name])
        same(sampler$model, before, target, name, "restated")
        cleared <- list(NULL)
        names(cleared) <- name
        doorAny(sampler, kind, target, cleared)
        same(sampler$model, before, target, name, "NULL")
        same(readerAny(sampler, kind, target), read, target, name, "reader")
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
  same(sampler$getForest(1L)$sd, 1.2, "reader")
  same(draws(sampler), draws(build("stated")), "draws")
})
case("inert on a run sampler with updateState", {
  sampler <- build("factor")
  invisible(sampler$run(0L, 2L))
  stored <- sampler$state
  sampler$setForest(
    2L,
    power = sampler$getForest(2L)$power,
    updateState = FALSE
  )
  same(sampler$state, stored, "no change, told not to store, stores nothing")
  sampler$setForest(
    2L,
    power = sampler$getForest(2L)$power,
    updateState = TRUE
  )
  expect_false(
    identical(sampler$state, stored),
    info = "inert on a run sampler with updateState: no change stores when told"
  )
  stored <- sampler$state
  invisible(sampler$run(0L, 2L))
  sampler$setForest(2L, power = 0.5, updateState = TRUE)
  expect_false(
    identical(sampler$state, stored),
    info = "inert on a run sampler with updateState: a change stores when told"
  )
  expectFound(
    refusal(sampler$setForest(2L, power = 3, updateState = "a")),
    "updateState",
    "updateState checked"
  )
})
# an inert call stores state when told to at every door that takes the argument:
# $setLeafPrior's, $setForest's and $setVarianceForest's alike
case("inert call stores state at each door", {
  stores <- function(live, call) {
    invisible(live$run(0L, 2L))
    stored <- live$state
    call(live)
    !identical(live$state, stored)
  }
  expect_true(
    stores(build("lone"), function(live) {
      live$setLeafPrior(live$getLeafPrior()$leaf.prior, updateState = TRUE)
    }),
    info = "setLeafPrior"
  )
  expect_true(
    stores(build("lone"), function(live) {
      live$setForest(1L, power = live$getForest(1L)$power, updateState = TRUE)
    }),
    info = "setForest, lone"
  )
  expect_true(
    stores(build("lone"), function(live) live$setForest(updateState = TRUE)),
    info = "setForest, naming nothing"
  )
  expect_true(
    stores(buildVariance(), function(live) {
      live$setVarianceForest(
        power = live$getVarianceForest()$power,
        updateState = TRUE
      )
    }),
    info = "setVarianceForest"
  )
  expect_false(
    stores(build("lone"), function(live) {
      live$setForest(1L, power = live$getForest(1L)$power, updateState = FALSE)
    }),
    info = "setForest, told not to"
  )
})

# a value only near the one in force is another value, and is written
for (name in c("sd", "power", "base")) {
  info <- paste("a value near the one in force is a change,", name)
  case(info, {
    live <- build("stated")
    inForce <- live$getForest(1L)[[name]]
    near <- inForce * (1 + 1e-10)
    stated <- list(near)
    names(stated) <- name
    door(live, 1L, stated)
    same(live$getForest(1L)[[name]], near, "reader")
    twin <- build("stated", stated)
    same(modelOf(live), modelOf(twin), "records")
    same(draws(live), draws(twin), "draws")
  })
}

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
    same(sampler$getForest(target)$n.perturb.cuts, 7L, "reader")
    same(draws(sampler), draws(plain))
    # a forest that stated its own does not follow, the later write in force
    stating <- build(kind)
    door(stating, target, list(n.perturb.cuts = 5))
    if (length(stating$model@forests) == 1L) {
      # no forest would take the window: refused, naming the door that does
      expect_error(
        setPerturbCuts(stating, 7L),
        "$setForest(forest, n.perturb.cuts = )",
        fixed = TRUE
      )
    } else {
      setPerturbCuts(stating, 7L)
    }
    same(stating$getForest(target)$n.perturb.cuts, 5, "a statement holds")
  })
}
case("follows the fitting function's door after a round trip", {
  sampler <- build("lone")
  plain <- build("lone")
  roundTrip(sampler, 1L)
  sampler$setLeafPrior(normal(k = 3))
  plain$setLeafPrior(normal(k = 3))
  same(sampler$model, plain$model, "leaf prior reaches")
  same(draws(sampler), draws(plain), "leaf prior draws")
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
  same(sampler$getForest()$power, 3, "setModel reaches")
  same(draws(sampler), draws(plain), "setModel draws")
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
    same(sampler$model, before)
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
    info = "a stated depth is stated, and a write states it once: stated in creation's order"
  )
  sampler$setForest(2L, power = 3)
  same(sum(sampler$model@forests[[2L]]@stated == "power"), 1L, "stated once")
})
case("a numeric forest's default re-derives after a round trip", {
  sampler <- build("numeric")
  plain <- build("numeric")
  roundTrip(sampler, 2L)
  sampler$setForestBasis(2L, 3 * dose, updateBasisScale = TRUE)
  plain$setForestBasis(2L, 3 * dose, updateBasisScale = TRUE)
  same(sampler$getForest(2L)$sd, plain$getForest(2L)$sd, "the default")
  expect_false(
    identical(sampler$getForest(2L)$sd, build("numeric")$getForest(2L)$sd),
    info = "a numeric forest's default re-derives after a round trip: the default moved"
  )
  same(draws(sampler), draws(plain), "draws")
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
  same(modelOf(live), modelOf(twin), "records")
  agree(live, "k after sd")
  same(draws(live), draws(twin), "draws")
})
case("later write, a probit forest(sd = ) returned to a drawn k", {
  live <- build("lone.probit", list(sd = 2))
  live$setLeafPrior(normal())
  expect_null(
    live$getForest()$sd,
    info = "later write, a probit forest(sd = ) returned to a drawn k: the drawn k has no sd"
  )
  same(
    live$model@leaf.prior,
    build("lone.probit")$model@leaf.prior,
    "the leaf prior of a drawn k"
  )
  # the drawn spread starts where the stated one left it, the one handoff
  expect_true(
    all(live$getK() == live$getLeafPrior()$k.scale / 2),
    info = "later write, a probit forest(sd = ) returned to a drawn k: from the spread in force"
  )
  roundTrip(live, 1L)
  agree(live, "drawn after sd")
  # it is drawn: a run moves it, and a door's sd holds it again
  start <- live$getK()
  invisible(live$run(0L, 3L))
  expect_false(
    identical(live$getK(), start),
    info = "later write, a probit forest(sd = ) returned to a drawn k: a run draws k"
  )
  live$setForest(1L, sd = 0.7)
  same(live$getForest()$sd, 0.7, "held again")
})
case("later write, power and setModel", {
  live <- build("lone", list(power = 0.5))
  model <- live$model
  model@tree.prior@power <- 3
  live$setModel(model)
  same(live$getForest()$power, 3, "setModel after a statement")
  twin <- build("lone", list(power = 3))
  same(draws(live), draws(twin), "setModel draws")
  live <- dbarts(
    x,
    y,
    power = 3,
    n.trees = 8L,
    control = ctl(proposal.probs = shipped)
  )
  live$setForest(1L, power = 0.5)
  twin <- build("lone", list(power = 0.5))
  same(live$getForest()$power, 0.5, "a door after the fit's")
  same(live$model@tree.prior@power, 0.5, "the record the fit's wrote")
  same(draws(live), draws(twin), "a door's draws")
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
  list(
    "sd past the leaf scale",
    "lone",
    list(sd = 1e308),
    "the leaf prior's sd overflows the leaf scale"
  ),
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
    same(atDoor, atCreation, "the same words")
  })
}
# A forest with no name of its own is named by its position at both doors, as
# creation names it, whatever its basis: a basis written as code is the label it
# is found by and not what a refusal calls it.
case("unlabeled forests are named by their label at both doors", {
  unlabeled <- list(
    list(list(), factorBasis),
    list(list(), list()),
    list(numericBasis, factorBasis),
    list(list(), numericBasis, numericBasis),
    list(numericBasis, numericBasis, factorBasis, list())
  )
  make <- function(statements, position = 0L, extra = list()) {
    if (position > 0L) {
      statements[[position]][names(extra)] <- extra
    }
    dbarts(
      x,
      y,
      forests = lapply(statements, function(a) do.call(forest, a)),
      n.trees = if (
        all(vapply(statements, function(a) !is.null(a$basis), NA))
      ) {
        NULL
      } else {
        8L
      },
      control = ctl(proposal.probs = shipped)
    )
  }
  for (index in seq_along(unlabeled)) {
    statements <- unlabeled[[index]]
    live <- make(statements)
    for (position in unique(c(2L, length(statements)))) {
      for (bad in list(
        list(proposal.probs = c(1, 2)),
        list(n.perturb.cuts = 0L),
        list(treeShift = "sometimes")
      )) {
        info <- paste("unlabeled", index, position, names(bad))
        atCreation <- refusal(make(statements, position, bad))
        atDoor <- refusal(door(live, position, bad))
        # the label the sampler gives it: the text of its basis, or its
        # position where it has none
        recorded <- live$model@forests[[position]]@label
        label <- if (identical(recorded, paste0("forest", position))) {
          position
        } else {
          paste0("\"", recorded, "\"")
        }
        expectFound(atCreation, paste0("forest ", label, "'s "), info)
        same(atDoor, atCreation, info, "the same words")
      }
      if (!is.null(statements[[position]]$basis)) {
        held <- make(statements)
        same(
          refusal(door(held, position, list(amplitude = fixed()))),
          "taken",
          paste("unlabeled", index, position, "a coefficient")
        )
        same(
          held$getForest(position)$amplitude,
          fixed(),
          paste("unlabeled", index, position, "a coefficient in force")
        )
      }
    }
  }
})
# a coefficient's prior, taken before a run as creation takes the statement; a
# release starts the coefficient at its held value and creation at its start
# scale, so its draws are the state twin's (test-coefficient-prior-live.R)
coefficientRows <- list(
  list("a changed coefficient on a basis forest", "factor", student(2), FALSE),
  list("normal() on a student forest", "factor.student", normal(), FALSE),
  list("another df", "numeric.student", student(3), FALSE),
  list("drawn on a held forest", "factor.held", normal(), TRUE),
  list("held on a drawn forest", "factor", fixed(), FALSE),
  list("drawn on a stated forest", "stated", normal(), TRUE),
  list("held on a drawn forest with no basis", "normal.stated", fixed(), FALSE)
)
for (row in coefficientRows) {
  info <- paste("a coefficient,", row[[1L]])
  case(info, {
    kind <- row[[2L]]
    target <- kinds[[kind]]$target
    live <- build(kind)
    same(
      refusal(door(live, target, list(amplitude = row[[3L]]))),
      "taken",
      info
    )
    twin <- build(kind, list(amplitude = row[[3L]]))
    same(live$getForest(target), twin$getForest(target), info, "reader")
    same(modelOf(live), modelOf(twin), info, "records")
    if (!row[[4L]]) {
      same(draws(live), draws(twin), info, "draws")
    }
  })
}
# a statement needs what it names beside it, merged and not name by name: student
# with a stated sd is not refused for the sd it lacks, which it has
case("merged statement", {
  stating <- build("stated")
  same(stating$getForest(1L)$sd, 1.2, "the sd in force")
  same(
    refusal(stating$setForest(1L, amplitude = student(1))),
    "taken",
    "student beside a stated sd is taken"
  )
  same(
    modelOf(stating),
    modelOf(build("stated", list(amplitude = student(1)))),
    "student beside a stated sd, the records"
  )
  sharing <- build("sharer")
  same(
    refusal(sharing$setForest(1L, sd = 1, amplitude = student(1))),
    "taken",
    "student stated with its sd is taken"
  )
  same(
    sharing$getForest(1L)[c("sd", "amplitude")],
    list(sd = 1, amplitude = student(1)),
    "student stated with its sd, in force"
  )
  sharing <- build("sharer")
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
  same(live$model, made()$model, "monotone round trip")
  live$setForest(1L, power = 0.5)
  same(live$getForest(1L)$power, 0.5, "monotone takes a depth")
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
  # the ninth value after the forest lands on n.trees, the first of the names
  # fixed at creation, and a fourteenth, unused, comes after those names
  expectFound(
    refusal(live$setForest(2L, 1, 2, 3, 4, 5, 6, 7, NULL, 9)),
    "'n.trees' is fixed at creation",
    "a tenth, unnamed"
  )
  expectFound(
    refusal(live$setForest(
      2L,
      1,
      2,
      3,
      4,
      5,
      6,
      7,
      NULL,
      9,
      10,
      11,
      12,
      13,
      14
    )),
    "'n.trees' is fixed at creation",
    "one too many, unnamed"
  )
  expectAsUntouched(live, build("factor"), 2L, "basis and an unknown name")
})

# A short name matches as it does at creation: the door's names include the ones
# it refuses, so a name that matches several there matches several here, and a
# short name of a refused one is refused as the full name is.
case("short names at the door and at creation", {
  whichArgument <- function(text) sub("argument [0-9]+ ", "argument N ", text)
  live <- build("lone")
  for (short in list(
    list(n = 50),
    list(b = 0.5),
    list(ba = 0.5),
    list(p = 0.5),
    list(i = 2L),
    list(s = 0.7)
  )) {
    info <- paste("short name", names(short))
    atCreation <- refusal(do.call(forest, short))
    atDoor <- refusal(door(build("lone"), 1L, short))
    if (identical(atCreation, "taken")) {
      # a name that is one name at creation is that name at the door
      expect_true(
        atDoor == "taken" || grepl("fixed at creation", atDoor),
        info = info
      )
    } else {
      same(whichArgument(atDoor), whichArgument(atCreation), info)
    }
  }
  expectFound(
    refusal(live$setForest(1L, n = 50)),
    "matches multiple formal arguments",
    "n is ambiguous"
  )
  expectFound(
    refusal(live$setForest(1L, b = 0.5)),
    "matches multiple formal arguments",
    "b is ambiguous"
  )
  same(
    refusal(live$setForest(1L, n.t = 200)),
    refusal(live$setForest(1L, n.trees = 200)),
    "n.t is n.trees"
  )
  same(
    refusal(live$setForest(1L, v = "x1")),
    refusal(live$setForest(1L, vars = "x1")),
    "v is vars"
  )
  same(
    refusal(live$setForest(1L, int = 2L)),
    refusal(live$setForest(1L, interactions = 2L)),
    "int is interactions"
  )
  # a short name that is one name stays that name
  live$setForest(1L, pow = 0.5)
  same(live$getForest(1L)$power, 0.5, "pow is power")
  expectAsUntouched(live, build("lone", list(power = 0.5)), 1L, "short names")
  untouched <- build("lone")
  invisible(refusal(untouched$setForest(1L, n = 50)))
  invisible(refusal(untouched$setForest(1L, b = 0.5)))
  expectAsUntouched(untouched, build("lone"), 1L, "ambiguous short names")
  # the variance door, against the variance forest's creation names
  variance <- buildVariance()
  for (short in list(
    list(n = 50),
    list(n.t = 5L),
    list(v = "x1"),
    list(b = 0.6)
  )) {
    info <- paste("variance short name", names(short))
    atCreation <- refusal(do.call(varianceForest, short))
    atDoor <- refusal(do.call(variance$setVarianceForest, short))
    if (identical(atCreation, "taken")) {
      expect_true(
        atDoor == "taken" || grepl("fixed at creation", atDoor),
        info = info
      )
    } else {
      same(whichArgument(atDoor), whichArgument(atCreation), info)
    }
  }
  same(
    refusal(variance$setVarianceForest(n.t = 5L)),
    refusal(variance$setVarianceForest(n.trees = 5L)),
    "variance, n.t is n.trees"
  )
  expectFound(
    refusal(variance$setVarianceForest(n = 50)),
    "matches multiple formal arguments",
    "variance, n is ambiguous"
  )
  same(variance$getVarianceForest()$base, 0.6, "variance, b is base")
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
  list("factor", list(amplitude = student(2), power = -1)),
  list("factor", list(sd = -1, amplitude = fixed())),
  list("factor", list(amplitude = student(2), n.perturb.cuts = 0L)),
  list("stated", list(amplitude = normal(), base = -1)),
  list("sharer", list(power = 0.5, sd = -1)),
  list("sharer", list(n.perturb.cuts = 5, sd = -1)),
  list("lone", list(power = 0.5, sd = -1)),
  list("lone", list(power = 0.5, n.perturb.cuts = 5, amplitude = normal())),
  list("lone", list(power = 0.5, sd = 1e308)),
  list("lone", list(n.perturb.cuts = 5, sd = 1e308))
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
# it. A lone forest's sd of 1e308 above is refused before any write, as its
# leaf scale overflows; here is the engine's refusal after the writes, on a
# forest whose sd has no such refusal to reach.
refusingWrite <- function(sampler, forest, ..., lone = FALSE) {
  writer <- dbarts:::writeForestDoor
  stubs <- list(
    .Call = function(symbol, ...) {
      if (identical(symbol, dbarts:::C_dbarts_bartcore_setForestSd)) {
        stop("the engine refused")
      }
      base::.Call(symbol, ...)
    }
  )
  # a lone forest's sd goes through the leaf-prior install
  if (lone) {
    stubs$writeLeafPrior <- function(...) stop("the engine refused")
  }
  environment(writer) <- list2env(
    stubs,
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

# a lone forest's sd goes through the leaf-prior install, the last of its writes
case("a lone write refused after writes", {
  live <- build("lone")
  message <- refusal(refusingWrite(
    live,
    1L,
    sd = 0.7,
    power = 0.5,
    n.perturb.cuts = 5,
    proposal.probs = moved,
    lone = TRUE
  ))
  expectFound(message, "the engine refused", "it refused")
  expectAsUntouched(live, build("lone"), 1L, "the lone write undone")
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
    same(modelOf(live), modelOf(twin), "records")
    same(live$getVarianceForest(), twin$getVarianceForest(), "reader")
    # the mean forest does not move
    same(
      live$getForest(),
      buildVariance()$getForest(),
      "leaves the mean forest"
    )
    liveDraws <- draws(live)
    same(liveDraws, draws(twin), "draws")
    expect_false(
      identical(liveDraws, draws(buildVariance())),
      info = paste(info, "moves a draw")
    )
    live$storeState()
    twin$setState(live$state)
    copied <- live$copy()
    loaded <- reloaded(live)
    after <- draws(twin)
    same(draws(copied), after, "copy")
    same(draws(loaded), after, "reload")
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
  same(sampler$model, before, "variance round trip")
  sampler$setVarianceForest(power = NULL)
  same(sampler$model, before, "variance NULL")
  same(draws(sampler), draws(buildVariance()), "variance inert draws")
  expect_null(
    build("lone")$getVarianceForest(),
    info = "the variance reader: no variance forest reads"
  )
})
case("the variance door stores when told", {
  sampler <- buildVariance()
  invisible(sampler$run(0L, 2L))
  stored <- sampler$state
  sampler$setVarianceForest(
    power = sampler$getVarianceForest()$power,
    updateState = TRUE
  )
  expect_false(
    identical(sampler$state, stored),
    info = "the variance door stores when told: variance, no change stores too"
  )
  stored <- sampler$state
  invisible(sampler$run(0L, 2L))
  sampler$setVarianceForest(power = 0.5, updateState = TRUE)
  expect_false(
    identical(sampler$state, stored),
    info = "the variance door stores when told: variance, a change stores"
  )
})
case("the mean door leaves the variance forest", {
  live <- buildVariance(variance = list(power = 2))
  live$setForest(1L, power = 0.5)
  twin <- buildVariance(mean = list(power = 0.5), variance = list(power = 2))
  same(modelOf(live), modelOf(twin), "records")
  same(draws(live), draws(twin), "draws")
  same(live$getVarianceForest()$power, 2, "reader")
  # it takes the mean forest's depth at creation, and is on its own after
  live <- buildVariance()
  live$setForest(1L, base = 0.6)
  twin <- buildVariance(mean = list(base = 0.6), variance = list(base = 0.95))
  same(draws(live), draws(twin), "the default depth stays")
})
case("setModel leaves the variance forest", {
  live <- buildVariance(variance = list(power = 2))
  model <- live$model
  model@tree.prior@power <- 0.5
  live$setModel(model)
  same(live$getVarianceForest()$power, 2, "reader")
  twin <- buildVariance(mean = list(power = 0.5), variance = list(power = 2))
  same(draws(live), draws(twin), "draws")
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
  list(
    "rule_gibbs",
    list(
      proposal.probs = c(
        birth_death = 0.4,
        swap = 0.1,
        change = 0.2,
        perturb = 0.2,
        rule_gibbs = 0.1,
        birth = 0.5
      )
    ),
    "gives 'rule_gibbs' a positive share"
  )
)
# a variance call refused on one name writes none of the others
for (stated in list(
  list(power = 0.5, base = 2),
  list(power = 0.5, n.perturb.cuts = 0L),
  list(n.perturb.cuts = 5, proposal.probs = c(1, 2)),
  list(base = 0.6, power = -1)
)) {
  info <- paste("variance refused whole,", paste(names(stated), collapse = " "))
  case(info, {
    live <- buildVariance()
    expect_false(
      identical(refusal(do.call(live$setVarianceForest, stated)), "taken"),
      info = paste(info, "refused")
    )
    same(
      live$getVarianceForest(),
      buildVariance()$getVarianceForest(),
      "reader"
    )
    same(modelOf(live), modelOf(buildVariance()), "records")
    same(draws(live), draws(buildVariance()), "draws")
  })
}
for (row in varianceRefusals) {
  info <- paste("variance refuses,", row[[1L]])
  case(info, {
    live <- buildVariance()
    message <- refusal(do.call(live$setVarianceForest, row[[2L]]))
    expect_false(identical(message, "taken"), info = info)
    expectFound(message, row[[3L]], info)
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
  same(first$getForest(2L)$power, 0.5, "position 2")
  same(second$getForest(1L)$power, 0.5, "position 1")
  expect_false(
    identical(second$getForest(2L)$power, 0.5),
    info = "a label, where a formula and a list put it elsewhere: the other forest stays"
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
  same(modelOf(second), modelOf(twin), "forest 1's record")
  same(draws(second), draws(twin), "forest 1's draws")
  byFormula <- dbarts(
    y ~ forest(x1) + forest(x2, basis = z),
    frame,
    n.trees = 8L,
    control = ctl(proposal.probs = shipped)
  )
  byFormula$setForest("z", power = 0.5)
  same(byFormula$getForest(2L)$power, 0.5, "a formula's label")
  same(byFormula$getForest("z"), byFormula$getForest(2L), "read")
  byFormula$setForest("forest1", base = 0.6)
  same(byFormula$getForest(1L)$base, 0.6, "forest<i>")
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
    "forest index out of range",
    "multinomial, out of range"
  )
  expectAsUntouched(live, buildMultinomial(), NULL, "refused setters")
  # the reader gives the one statement, whichever category was named and whatever
  # a hand write left on one of them
  handWritten <- buildMultinomial()
  .Call(
    dbarts:::C_dbarts_bartcore_setForestTreePrior,
    handWritten$getPointer(),
    1L,
    0.6,
    0.5
  )
  same(handWritten$getForest(2L), read, "one statement")
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
  live <- buildMultinomial()
  message <- refusal(live$setForest(sd = 1))
  expect_false(identical(message, "taken"), info = "door sd")
  expectFound(message, "multinomial", "door reason sd")
  expectAsUntouched(live, buildMultinomial(), NULL, "refused sd")
  same(
    refusal(live$setForest(sd = 1)),
    refusal(buildMultinomial(list(sd = 1))),
    "an sd in the words of creation"
  )
})
# no category has a coefficient: fixed() states what is so, at creation and at
# the door as an inert write, and student() and normal() are refused by
# creation's words at both
case("the multinomial door and amplitude", {
  live <- buildMultinomial()
  before <- live$model
  live$setForest(amplitude = fixed())
  same(live$model, before, "fixed() writes nothing")
  expect_false(
    "amplitude" %in% names(live$getForest()),
    info = "the multinomial reader has no amplitude"
  )
  roundTripped <- buildMultinomial()
  do.call(
    roundTripped$setForest,
    c(roundTripped$getForest(), amplitude = fixed())
  )
  same(
    modelOf(roundTripped),
    modelOf(buildMultinomial()),
    "fixed() with a read"
  )
  same(
    draws(live),
    draws(buildMultinomial(list(amplitude = fixed()))),
    "draws as creation with fixed()"
  )
  expectAsUntouched(
    buildMultinomial(),
    buildMultinomial(list(amplitude = fixed())),
    NULL,
    "as created with fixed()"
  )
  for (prior in list(student(2), normal())) {
    word <- if (is(prior, "dbartsNormalPrior")) "normal()" else "student(2)"
    live <- buildMultinomial()
    atDoor <- refusal(live$setForest(amplitude = prior))
    atCreation <- refusal(buildMultinomial(list(amplitude = prior)))
    expect_false(identical(atDoor, "taken"), info = paste("door", word))
    same(atDoor, atCreation, "the words of creation,", word)
    expectAsUntouched(live, buildMultinomial(), NULL, paste("refused", word))
    # with a good name beside it, nothing is written
    atomic <- buildMultinomial()
    invisible(refusal(atomic$setForest(power = 0.5, amplitude = prior)))
    expectAsUntouched(atomic, buildMultinomial(), NULL, paste("atomic", word))
  }
})
# the reader gives a category's moves from the statement, not the control, and
# the door writes them to every category
case("the multinomial door and the moves a statement gave", {
  stated <- list(n.perturb.cuts = 5, proposal.probs = moved)
  live <- buildMultinomial(stated)
  read <- live$getForest()
  same(read$n.perturb.cuts, 5, "the cuts stated")
  same(read$proposal.probs[["perturb"]], 0.4, "the mixture stated")
  for (which in list(1L, 2L, 3L, "b", "c")) {
    same(live$getForest(which), read, "one list stating,", which)
  }
  live$setForest(n.perturb.cuts = 7)
  twin <- buildMultinomial(list(n.perturb.cuts = 7, proposal.probs = moved))
  same(modelOf(live), modelOf(twin), "records after a write")
  same(draws(live), draws(twin), "draws after a write")
  for (which in list(1L, 3L, "c")) {
    same(live$getForest(which)$n.perturb.cuts, 7, "written,", which)
  }
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
  same(sampler$getForest("b")$power, 0.5, "one statement")
})
case("a level read by every forest-indexed reader", {
  sampler <- buildMultinomial()
  invisible(sampler$run(0L, 2L))
  readers <- list(
    getForestFits = sampler$getForestFits,
    getForestVariableCounts = sampler$getForestVariableCounts
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
  same(
    sampler$getTrees(forest = "b"),
    sampler$getTrees(forest = 2L),
    "getTrees"
  )
  same(
    names(sampler$getLeafPrior()),
    c("a", "b", "c"),
    "a NULL read names its list by levels"
  )
  # a level spelled forest<i> at another position names two forests and is refused
  colliding <- factor(
    sample(c("forest2", "b", "c"), n, TRUE),
    levels = c("forest2", "b", "c")
  )
  spelled <- buildMultinomial(response = colliding)
  expectFound(
    refusal(spelled$getForest("forest2")),
    "'forest' (\"forest2\") is the label of forest 1 and the name of position 2",
    "a level at another position's name, by the conflict"
  )
  same(
    spelled$getForest("b"),
    spelled$getForest(2L),
    "the level beside it reads by its position"
  )
})
