# A live sampler's forest takes another coefficient prior through
# $setForest(amplitude = ): fixed() holds the coefficient at its held value at
# once, normal() or student(df) on a held forest releases it, drawn from that
# value, and between normal() and student(df), or two dfs, the coefficient keeps
# its value, a t's variance starting at the normal's. No leaf or trees' scale
# moves but where the fitting function's leaf prior sizes the forest, which
# then takes the size creation gives the statement.
#
# Oracles, every one on two chains: the creation twin before any run, where the
# start is creation's; the state twin otherwise, the live sampler's state
# installed into it and into a creation of the statement the call leaves; the
# reader-to-creation twin for the last forest left of a group; and later writes
# on both. Every case is named, and an error in one is a failed expectation.

forest <- dbartsForests$forest
varianceForest <- dbartsForests$varianceForest
fixed <- dbartsPriors$fixed
student <- dbartsPriors$student
normal <- dbartsPriors$normal
chi <- dbartsPriors$chi

set.seed(20261103L)
n <- 100L
x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
z <- factor(rbinom(n, 1L, 0.5))
w <- rnorm(n)
y <- sin(3 * x1) + (z == "1") * x2 + 0.5 * w * x3 + rnorm(n, sd = 0.3)
yb <- as.integer(y > stats::median(y))
d <- data.frame(y, yb, x1, x2, x3, z, w)

ctl <- function() {
  dbartsControl(
    n.chains = 2L,
    n.threads = 1L,
    n.samples = 5L,
    n.burn = 0L,
    updateState = FALSE,
    verbose = FALSE
  )
}
make <- function(formula, ...) {
  dbarts(
    stats::as.formula(formula),
    d,
    n.trees = 24L,
    control = ctl(),
    seed = 7L,
    ...
  )
}
# a model's records with what a formula holds of its frame left out
recordsOf <- function(sampler) {
  model <- sampler$model
  for (index in seq_along(model@forests)) {
    model@forests[[index]]@basis.term <- NULL
  }
  model
}
draws <- function(sampler) {
  sampler$run(0L, 5L)[c("sigma", "train", "k")]
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
caseName <- ""
same <- function(current, target, ...) {
  expect_identical(current, target, info = paste(caseName, ...))
}
case <- function(name, body) {
  caseName <<- name
  tryCatch(
    body,
    error = function(e) {
      expect_true(FALSE, info = paste0(name, " stopped: ", conditionMessage(e)))
    }
  )
  invisible(NULL)
}
# each chain's coefficient blocks and variances off its stored state: the
# state's glue holds the forest count, the widths, the blocks and a variance
# for each forest
glueOf <- function(sampler, index) {
  sampler$storeState()
  lapply(sampler$state, function(chain) {
    glue <- chain$glue
    numForests <- glue[[1L]]
    widths <- glue[1L + seq_len(numForests)]
    start <- 1L + numForests + sum(widths[seq_len(index - 1L)])
    list(
      block = glue[start + seq_len(widths[[index]])],
      variance = glue[[1L + numForests + sum(widths) + index]]
    )
  })
}
# the live sampler's state installed into it and into `twin`, then both run
installBoth <- function(live, twin) {
  live$storeState()
  state <- live$state
  twin$setState(state)
  live$setState(state)
  invisible(NULL)
}
# the normal's variance on each kind: 1 / z^2 with no basis, kappa / 0.674^2
# with one, kappa 1/2 on a level
normalVariance <- c(
  none = 1 / stats::qnorm(0.75)^2,
  numeric = 1 / 0.674^2,
  factor = 0.5 / 0.674^2
)
expect_equal(
  unname(normalVariance),
  c(2.1981, 2.2013, 1.1007),
  tolerance = 1e-4
)

# --- each transition, from a statement to another ---
# `from` and `to` are the forests' statements; `exempt` where a release starts
# the block at its held value and creation at its start scale, or a t's
# variance with no basis at the normal's and creation at 1 / m^2
factorModel <- function(statement) {
  paste0("y ~ forest(x1) + forest(x2, basis = ~ z", statement, ")")
}
numericModel <- function(statement) {
  paste0("y ~ forest(x1) + forest(x2, basis = ~ w", statement, ")")
}
statingModel <- function(statement) {
  paste0("y ~ forest(x1, sd = 0.7", statement, ") + forest(x2, basis = ~ z)")
}
transitions <- list(
  list(
    "factor hold",
    factorModel(""),
    2L,
    fixed(),
    factorModel(", amplitude = fixed()"),
    FALSE,
    "factor",
    c(0, 1)
  ),
  list(
    "factor release",
    factorModel(", amplitude = fixed()"),
    2L,
    normal(),
    factorModel(""),
    TRUE,
    "factor",
    c(0, 1)
  ),
  list(
    "factor normal to student",
    factorModel(""),
    2L,
    student(3),
    factorModel(", amplitude = student(3)"),
    FALSE,
    "factor",
    NULL
  ),
  list(
    "factor student to normal",
    factorModel(", amplitude = student(3)"),
    2L,
    normal(),
    factorModel(""),
    FALSE,
    "factor",
    NULL
  ),
  list(
    "factor df to df",
    factorModel(", amplitude = student(3)"),
    2L,
    student(5),
    factorModel(", amplitude = student(5)"),
    FALSE,
    "factor",
    NULL
  ),
  list(
    "numeric hold",
    numericModel(""),
    2L,
    fixed(),
    numericModel(", amplitude = fixed()"),
    FALSE,
    "numeric",
    1
  ),
  list(
    "numeric release to student",
    numericModel(", amplitude = fixed()"),
    2L,
    student(3),
    numericModel(", amplitude = student(3)"),
    TRUE,
    "numeric",
    1
  ),
  list(
    "numeric student to normal",
    numericModel(", amplitude = student(3)"),
    2L,
    normal(),
    numericModel(""),
    FALSE,
    "numeric",
    NULL
  ),
  list(
    "stated hold",
    statingModel(", amplitude = normal()"),
    1L,
    fixed(),
    statingModel(", amplitude = fixed()"),
    FALSE,
    "none",
    1
  ),
  list(
    "stated release",
    statingModel(""),
    1L,
    normal(),
    statingModel(", amplitude = normal()"),
    TRUE,
    "none",
    1
  ),
  list(
    "stated normal to student",
    statingModel(", amplitude = normal()"),
    1L,
    student(3),
    statingModel(", amplitude = student(3)"),
    TRUE,
    "none",
    NULL
  ),
  list(
    "stated student to normal",
    statingModel(", amplitude = student(3)"),
    1L,
    normal(),
    statingModel(", amplitude = normal()"),
    FALSE,
    "none",
    NULL
  ),
  list(
    "stated df to df",
    statingModel(", amplitude = student(1)"),
    1L,
    student(3),
    statingModel(", amplitude = student(3)"),
    FALSE,
    "none",
    NULL
  )
)
for (row in transitions) {
  case(paste("before a run,", row[[1L]]), {
    live <- make(row[[2L]])
    live$setForest(row[[3L]], amplitude = row[[4L]])
    twin <- make(row[[5L]])
    same(recordsOf(live), recordsOf(twin), "records")
    same(live$getForest(row[[3L]]), twin$getForest(row[[3L]]), "the reader")
    if (!is.null(row[[8L]])) {
      blocks <- live$getForestAmplitudes(row[[3L]])
      same(
        unname(blocks),
        matrix(row[[8L]], length(row[[8L]]), 2L),
        "the block at its held value"
      )
    }
    if (!row[[6L]]) {
      same(draws(live), draws(twin), "draws")
    }
  })
  case(paste("after a run, the state twin,", row[[1L]]), {
    live <- make(row[[2L]])
    invisible(live$run(0L, 5L))
    before <- glueOf(live, row[[3L]])
    live$setForest(row[[3L]], amplitude = row[[4L]])
    after <- glueOf(live, row[[3L]])
    toStudent <- is(row[[4L]], "dbartsStudentPrior")
    for (chain in 1:2) {
      if (row[[6L]] && !is.null(row[[8L]])) {
        same(
          after[[chain]]$block,
          row[[8L]],
          "a release starts at the held value"
        )
      } else if (!is.null(row[[8L]])) {
        same(after[[chain]]$block, row[[8L]], "a hold sets the block")
      } else {
        same(after[[chain]]$block, before[[chain]]$block, "the block kept")
      }
      fromStudent <- grepl("student", row[[2L]], fixed = TRUE)
      if (toStudent && !fromStudent) {
        same(
          after[[chain]]$variance,
          normalVariance[[row[[7L]]]],
          "a t's variance starts at the normal's"
        )
      } else if (toStudent) {
        same(after[[chain]]$variance, before[[chain]]$variance, "a df keeps v")
      }
    }
    twin <- make(row[[5L]])
    installBoth(live, twin)
    same(recordsOf(live), recordsOf(twin), "records")
    same(draws(live), draws(twin), "draws")
  })
}
case("a hold under probit, before a run", {
  live <- make("yb ~ forest(x1) + forest(x2, basis = ~ z)")
  live$setForest(2L, amplitude = fixed())
  twin <- make("yb ~ forest(x1) + forest(x2, basis = ~ z, amplitude = fixed())")
  same(recordsOf(live), recordsOf(twin), "records")
  same(draws(live), draws(twin), "draws")
})
case("a df changed keeps the variance drawn; student to normal resets it", {
  live <- make(statingModel(", amplitude = student(3)"))
  invisible(live$run(0L, 5L))
  drawn <- vapply(glueOf(live, 1L), function(g) g$variance, 0)
  live$setForest(1L, amplitude = student(5))
  same(vapply(glueOf(live, 1L), function(g) g$variance, 0), drawn, "kept")
  live$setForest(1L, amplitude = normal())
  same(
    vapply(glueOf(live, 1L), function(g) g$variance, 0),
    rep(normalVariance[["none"]], 2L),
    "reset"
  )
})

# --- a forest the fitting function's leaf prior sizes whole ---
wholeModels <- c(
  alone = "~ forest(x1) + forest(x2, basis = ~ z)",
  beside = "~ forest(x1) + forest(x2, sd = 0.5)",
  both = "~ forest(x1) + forest(x2, sd = 0.5) + forest(x3, basis = ~ z)"
)
withFirst <- function(model, statement) {
  sub("forest(x1)", paste0("forest(x1", statement, ")"), model, fixed = TRUE)
}
for (response in c("y", "yb")) {
  for (shape in names(wholeModels)) {
    formula <- paste(response, wholeModels[[shape]])
    case(paste("whole total released and held again,", response, shape), {
      live <- make(formula)
      live$setForest(1L, amplitude = student(3))
      released <- make(withFirst(formula, ", amplitude = student(3)"))
      same(recordsOf(live), recordsOf(released), "records released")
      same(live$getForest(1L), released$getForest(1L), "the reader released")
      same(unname(live$getK(1L)), c(1, 1), "k pinned at 1")
      same(
        unname(live$getForestAmplitudes(1L)),
        matrix(1, 1L, 2L),
        "the block at its held value"
      )
      ran <- live$run(0L, 3L)
      same(is.null(ran$k), TRUE, "no forest draws k, so the run records none")
      installBoth(live, released)
      same(draws(live), draws(released), "the state twin released")
      kWas <- live$getLeafPrior(1L)$sd
      live$setForest(1L, amplitude = fixed())
      held <- make(withFirst(formula, ", amplitude = fixed()"))
      same(recordsOf(live), recordsOf(held), "records held again")
      same(live$getForest(1L), held$getForest(1L), "the reader held again")
      if (response == "yb") {
        # a drawn k starts at the k that is the spread in force
        same(
          unname(live$getLeafPrior(1L)$k.scale / live$getK(1L)),
          rep(unname(kWas), 2L),
          "a drawn k starts at the spread in force"
        )
      }
      installBoth(live, held)
      same(draws(live), draws(held), "the state twin held again")
    })
  }
}
case("whole total released and held again before a run, held k", {
  for (formula in paste("y", wholeModels)) {
    live <- make(formula)
    live$setForest(1L, amplitude = normal())
    live$setForest(1L, amplitude = fixed())
    same(
      draws(live),
      draws(make(withFirst(formula, ", amplitude = fixed()"))),
      formula
    )
  }
  # placed as creation places it under a held k other than 2 or a named sd,
  # each of which gives a group of one where k 2 gives the held entry
  leafPriors <- list("k 3" = normal(k = 3), "named sd" = normal(sd = 0.6))
  for (label in names(leafPriors)) {
    for (formula in paste("y", wholeModels[c("beside", "both")])) {
      live <- make(formula, leaf.prior = leafPriors[[label]])
      live$setForest(1L, amplitude = student(3))
      live$setForest(1L, amplitude = fixed())
      twin <- make(
        withFirst(formula, ", amplitude = fixed()"),
        leaf.prior = leafPriors[[label]]
      )
      same(recordsOf(live), recordsOf(twin), "records,", label, formula)
      same(draws(live), draws(twin), "draws,", label, formula)
    }
  }
})

# --- the last forest left of a group, released ---
lastLeftModel <- function(response, s3 = "") {
  paste0(
    response,
    " ~ forest(x1) + forest(x2) + forest(x3",
    s3,
    ") + forest(x2, basis = ~ z)"
  )
}
for (response in c("y", "yb")) {
  case(paste("the last left released,", response), {
    live <- make(lastLeftModel(response))
    live$setForest(1L, sd = 0.4)
    live$setForest(2L, sd = 0.3)
    before <- live$getLeafPrior(3L)
    hyperprior <- live$model@leaf.hyperprior
    invisible(live$run(0L, 3L))
    live$setForest(3L, amplitude = student(3))
    read <- live$getForest(3L)
    part <- read$sd
    same(live$model@forests[[3L]]@sd, part, "its part recorded as its sd")
    same(is.na(live$model@forests[[3L]]@share), TRUE, "no share")
    if (response == "y") {
      expect_equal(part, unname(before$sd), tolerance = 1e-12, info = "held k")
    } else {
      # a drawn k gives no one sd: the part is taken at the prior's median k
      medianK <- hyperprior@scale *
        sqrt(stats::qchisq(0.5, hyperprior@degreesOfFreedom))
      expect_equal(
        part,
        unname(before$k.scale) / medianK,
        tolerance = 1e-12,
        info = "drawn k"
      )
    }
    twin <- dbarts(
      stats::as.formula(paste0(
        response,
        " ~ forest(x1, sd = 0.4) + forest(x2, sd = 0.3) + ",
        "forest(x3, sd = part, amplitude = student(3)) + ",
        "forest(x2, basis = ~ z)"
      )),
      d,
      n.trees = 24L,
      control = ctl(),
      seed = 7L
    )
    same(recordsOf(live), recordsOf(twin), "the reader's statement creates it")
    installBoth(live, twin)
    same(draws(live), draws(twin), "the reader-to-creation twin")
  })
}

# --- an sd and a release in one call on a sharer of three ---
case("an sd with a release on a sharer", {
  live <- make(lastLeftModel("y"))
  same(
    refusal(live$setForest(1L, amplitude = student(3))),
    refusal(make(sub(
      "forest(x1)",
      "forest(x1, amplitude = student(3))",
      lastLeftModel("y"),
      fixed = TRUE
    ))),
    "refused alone, as creation refuses it"
  )
  live$setForest(1L, sd = 0.4, amplitude = student(3))
  recorded <- new("dbartsSampler", live$control, live$model, live$data)
  same(recordsOf(recorded), recordsOf(live), "records")
  installBoth(live, recorded)
  same(draws(live), draws(recorded), "the recorded twin")
  same(
    vapply(live$model@forests, methods::slot, 0, "share"),
    c(NA, 1 / 3, 1 / 3, NA),
    "the others keep their parts"
  )
})

# --- after a run: readers, copy and reload ---
case("after a run, the change in force, a copy and a reload", {
  live <- make(factorModel(""))
  invisible(live$run(0L, 5L))
  live$setForest(2L, amplitude = fixed())
  same(live$getForest(2L)$amplitude, fixed(), "the reader")
  same(live$getLeafPrior(2L)$amplitude, fixed(), "the leaf prior's reader")
  same(unname(live$getForestAmplitudes(2L)), matrix(c(0, 1), 2L, 2L), "held")
  live$setForest(2L, amplitude = student(3))
  invisible(live$run(0L, 2L))
  live$storeState()
  copied <- live$copy()
  file <- tempfile(fileext = ".rds")
  saveRDS(live, file)
  loaded <- readRDS(file)
  unlink(file)
  # each restarts from the stored state, as a creation of the statement in
  # force does with that state installed
  twin <- make(factorModel(", amplitude = student(3)"))
  twin$setState(live$state)
  same(recordsOf(copied), recordsOf(twin), "the copy's records")
  same(recordsOf(loaded), recordsOf(twin), "the reload's records")
  ran <- draws(twin)
  same(draws(copied), ran, "the copy continues")
  same(draws(loaded), ran, "the reload continues")
})

# --- later writes after each route, against the creation twin ---
laterWrites <- list(
  list(
    "a df changed",
    statingModel(", amplitude = normal()"),
    1L,
    student(3),
    statingModel(", amplitude = student(3)")
  ),
  list(
    "whole total released",
    "y ~ forest(x1) + forest(x2, basis = ~ z)",
    1L,
    student(3),
    "y ~ forest(x1, amplitude = student(3)) + forest(x2, basis = ~ z)"
  ),
  list(
    "whole total released beside sd",
    "y ~ forest(x1) + forest(x2, sd = 0.5)",
    1L,
    student(3),
    "y ~ forest(x1, amplitude = student(3)) + forest(x2, sd = 0.5)"
  ),
  list(
    "held again beside sd and a basis forest",
    "y ~ forest(x1, sd = 0.5) + forest(x2, amplitude = student(3)) + forest(x2, basis = ~ z)",
    2L,
    fixed(),
    "y ~ forest(x1, sd = 0.5) + forest(x2, amplitude = fixed()) + forest(x2, basis = ~ z)"
  ),
  list(
    "held again under probit",
    "yb ~ forest(x1, sd = 0.5) + forest(x2, amplitude = student(3)) + forest(x2, basis = ~ z)",
    2L,
    fixed(),
    "yb ~ forest(x1, sd = 0.5) + forest(x2, amplitude = fixed()) + forest(x2, basis = ~ z)"
  )
)
for (row in laterWrites) {
  for (write in c("sd", "leaf prior")) {
    case(paste("a later", write, "after", row[[1L]]), {
      live <- make(row[[2L]])
      live$setForest(row[[3L]], amplitude = row[[4L]])
      twin <- make(row[[5L]])
      apply <- function(sampler) {
        refusal(
          if (write == "sd") {
            sampler$setForest(row[[3L]], sd = 0.6)
          } else {
            sampler$setLeafPrior(normal(k = 3))
          }
        )
      }
      same(apply(live), apply(twin), "taken or refused alike")
      same(recordsOf(live), recordsOf(twin), "records")
      installBoth(live, twin)
      same(draws(live), draws(twin), "draws")
    })
  }
}

# --- inert: the round trip, NULL and the value in force ---
inertModels <- c(
  factorModel(""),
  factorModel(", amplitude = fixed()"),
  numericModel(", amplitude = student(3)"),
  statingModel(", amplitude = normal()"),
  "y ~ forest(x1) + forest(x2, basis = ~ z)",
  "yb ~ forest(x1, amplitude = student(3)) + forest(x2, sd = 0.5)",
  lastLeftModel("y")
)
for (formula in inertModels) {
  case(paste("inert,", formula), {
    untouched <- make(formula)
    for (target in seq_along(untouched$model@forests)) {
      live <- make(formula)
      do.call(live$setForest, c(list(target), live$getForest(target)))
      live$setForest(target, amplitude = NULL)
      live$setForest(target, amplitude = live$getForest(target)$amplitude)
      same(recordsOf(live), recordsOf(untouched), "records", target)
      same(draws(live), draws(make(formula)), "draws", target)
    }
  })
}

# --- names in one call are compared with the values in force before it ---
# an sd handed back beside a hold is no change wherever k is drawn, so the
# forest follows the leaf prior, k drawn, as after fixed() alone; stating the
# same sd in a second call, or a different one in the call, pins it
drawnK <- list(
  probit = list("yb ~ forest(x1) + forest(x2, basis = ~ z)"),
  logistic = list(
    "yb ~ forest(x1) + forest(x2, basis = ~ z)",
    family = "logistic"
  ),
  "gaussian, chi()" = list(
    "y ~ forest(x1) + forest(x2, basis = ~ z)",
    leaf.prior = normal(k = chi())
  )
)
released <- function(arguments, amplitude = student(3)) {
  live <- do.call(make, arguments)
  live$setForest(1L, amplitude = amplitude)
  live
}
releases <- list("student(3)" = student(3), "normal()" = normal())
for (label in names(drawnK)) {
  arguments <- drawnK[[label]]
  case(paste("an sd in force beside a hold is no change,", label), {
    routes <- list(
      "sd beside fixed()" = function(live) {
        live$setForest(1L, sd = live$getForest(1L)$sd, amplitude = fixed())
      },
      "the reader's list" = function(live) {
        given <- live$getForest(1L)
        given$amplitude <- fixed()
        do.call(live$setForest, c(list(1L), given))
      }
    )
    for (release in names(releases)) {
      for (route in names(routes)) {
        at <- paste(release, route)
        live <- released(arguments, releases[[release]])
        live$run(0L, 4L)
        same(is.null(live$getForest(1L)$sd), FALSE, at, "an sd in force")
        routes[[route]](live)
        twin <- released(arguments, releases[[release]])
        twin$run(0L, 4L)
        twin$setForest(1L, amplitude = fixed())
        same(live$model@forests[[1L]]@sd, NA_real_, at, "not stated")
        same(is.null(live$getForest(1L)$sd), TRUE, at, "the reader's sd")
        same(recordsOf(live), recordsOf(twin), at, "records")
        ran <- draws(live)
        same(any(ran$k != 1), TRUE, at, "k drawn")
        same(ran, draws(twin), at, "draws")
      }
    }
  })
  case(paste("the same sd in a second call is stated,", label), {
    live <- released(arguments)
    sd <- live$getForest(1L)$sd
    live$setForest(1L, sd = sd, amplitude = fixed())
    live$setForest(1L, sd = sd)
    same(live$model@forests[[1L]]@sd, sd, "stated")
    same(unname(live$getK(1L)), c(1, 1), "k pinned")
    # creation refuses a leaf prior's k that no forest takes, so the twin
    # takes the family's leaf prior, which the statement leaves in force
    twin <- do.call(
      make,
      c(
        list(withFirst(arguments[[1L]], ", sd = sd, amplitude = fixed()")),
        arguments[setdiff(names(arguments), c("", "leaf.prior"))]
      )
    )
    same(recordsOf(live), recordsOf(twin), "records")
    same(draws(live), draws(twin), "draws")
  })
}
case("a different sd beside a hold under probit is stated", {
  formula <- drawnK$probit[[1L]]
  live <- released(drawnK$probit)
  sd <- 2 * live$getForest(1L)$sd
  live$setForest(1L, sd = sd, amplitude = fixed())
  same(live$model@forests[[1L]]@sd, sd, "stated")
  same(unname(live$getK(1L)), c(1, 1), "k pinned")
  twin <- make(withFirst(formula, ", sd = sd, amplitude = fixed()"))
  same(recordsOf(live), recordsOf(twin), "records")
  same(draws(live), draws(twin), "draws")
})
# under a held k too the sd in force is no change, so the forest goes back to
# the leaf prior keeping its size
heldLeafPriors <- list(
  default = list(),
  "k 3" = list(leaf.prior = normal(k = 3)),
  "named sd" = list(leaf.prior = normal(sd = 0.6))
)
for (label in names(heldLeafPriors)) {
  arguments <- heldLeafPriors[[label]]
  case(paste("an sd beside a hold under gaussian, held k,", label), {
    formula <- "y ~ forest(x1) + forest(x2, basis = ~ z)"
    live <- do.call(make, c(list(formula), arguments))
    live$setForest(1L, amplitude = student(3))
    sd <- live$getForest(1L)$sd
    same(is.null(sd), FALSE, "the released forest's sd in force")
    live$setForest(1L, sd = sd, amplitude = fixed())
    same(live$model@forests[[1L]]@sd, NA_real_, "dropped")
    twin <- do.call(
      make,
      c(list(withFirst(formula, ", amplitude = fixed()")), arguments)
    )
    same(recordsOf(live), recordsOf(twin), "records")
    expect_equal(
      live$getForest(1L)$sd,
      sd,
      tolerance = 1e-12,
      info = paste(caseName, "the size kept")
    )
    same(draws(live), draws(twin), "draws")
  })
}
case("an sd in force beside another name is dropped", {
  formula <- "y ~ forest(x1) + forest(x2, sd = 0.5) + forest(x3, basis = ~ z)"
  live <- make(formula)
  live$setForest(2L, sd = live$getForest(2L)$sd, power = 3)
  other <- make(formula)
  other$setForest(2L, power = 3)
  same(recordsOf(live), recordsOf(other), "records")
  same(draws(live), draws(other), "draws")
})

# --- a later fitting-function leaf prior reaches a forest held again ---
case("the leaf prior after a release and a hold", {
  alone <- make("y ~ forest(x1) + forest(x2, basis = ~ z)")
  alone$setForest(1L, amplitude = student(3))
  alone$setLeafPrior(normal(k = 3))
  twin <- make(
    "y ~ forest(x1, amplitude = student(3)) + forest(x2, basis = ~ z)",
    leaf.prior = normal(k = 3)
  )
  same(alone$getForest(1L)$sd, twin$getForest(1L)$sd, "the one forest's size")
  beside <- make("y ~ forest(x1) + forest(x2, sd = 0.5)")
  beside$setForest(1L, amplitude = student(3))
  expect_true(
    grepl(
      "stating its own sd or a coefficient's prior",
      refusal(beside$setLeafPrior(normal(k = 3))),
      fixed = TRUE
    ),
    info = "beside a forest stating sd, none takes it"
  )
  beside$setForest(1L, amplitude = fixed())
  same(
    refusal(beside$setLeafPrior(normal(k = 3))),
    "taken",
    "held again, reached"
  )
  same(unname(beside$getK(1L)), c(3, 3), "at its k")
  # created so, the words are true of it too
  stating <- make(
    "yb ~ forest(x1, sd = 1) + forest(x2, amplitude = student(3))"
  )
  expect_true(
    grepl(
      "every forest with no basis stating its own sd or a coefficient's prior",
      refusal(stating$setLeafPrior(normal(k = 3))),
      fixed = TRUE
    ),
    info = "a sampler created stating a df"
  )
})

# --- refusals, as creation's, writing nothing ---
refusals <- list(
  list(
    "a sharer of two with no sd",
    "y ~ forest(x1) + forest(x2) + forest(x2, basis = ~ z)",
    1L,
    list(amplitude = student(3)),
    ", amplitude = student(3)",
    "share one total"
  ),
  list(
    "a size the leaf prior does not give",
    "y ~ forest(x1) + forest(x2, basis = ~ z)",
    1L,
    list(amplitude = student(3)),
    ", amplitude = student(3)",
    "gives no size"
  ),
  list(
    "a fixed() other than 1",
    factorModel(""),
    2L,
    list(amplitude = fixed(2)),
    NULL,
    NULL
  )
)
for (row in refusals) {
  case(paste("refused,", row[[1L]]), {
    sizeless <- grepl("size", row[[1L]], fixed = TRUE)
    built <- function(formula) {
      if (sizeless) {
        make(formula, leaf.prior = normal(k = chi(1.25, Inf)))
      } else {
        make(formula)
      }
    }
    live <- built(row[[2L]])
    atDoor <- refusal(do.call(live$setForest, c(list(row[[3L]]), row[[4L]])))
    expect_false(identical(atDoor, "taken"), info = caseName)
    if (!is.null(row[[6L]])) {
      expect_true(grepl(row[[6L]], atDoor, fixed = TRUE), info = atDoor)
      atCreation <- refusal(built(withFirst(row[[2L]], row[[5L]])))
      same(atDoor, atCreation, "creation's words")
    }
    same(recordsOf(live), recordsOf(built(row[[2L]])), "records unchanged")
    same(draws(live), draws(built(row[[2L]])), "draws unchanged")
  })
}
case("refused, the last left under a leaf prior with no size", {
  live <- make(lastLeftModel("y"), leaf.prior = normal(k = chi(1.25, Inf)))
  live$setForest(1L, sd = 0.4)
  live$setForest(2L, sd = 0.3)
  expect_true(
    grepl(
      "gives no size",
      refusal(live$setForest(3L, amplitude = normal())),
      fixed = TRUE
    )
  )
  same(is.na(live$model@forests[[3L]]@amplitude.df), TRUE, "nothing written")
})
case("refused, a model of one forest and a multinomial model", {
  lone <- make("y ~ forest(x1)")
  expect_true(grepl(
    "this model has one forest",
    refusal(lone$setForest(1L, amplitude = normal())),
    fixed = TRUE
  ))
  same(
    refusal(lone$setForest(1L, amplitude = fixed())),
    "taken",
    "fixed() changes nothing"
  )
  category <- factor(sample(c("a", "b", "c"), n, TRUE))
  softmax <- dbarts(
    cbind(x1, x2),
    category,
    family = "multinomial",
    forests = list(forest()),
    n.trees = 8L,
    control = ctl()
  )
  expect_false(identical(
    refusal(softmax$setForest(amplitude = student(3))),
    "taken"
  ))
  same(
    refusal(softmax$setForest(amplitude = fixed())),
    "taken",
    "fixed() alone"
  )
})

# --- $setModel refuses what its records would drop ---
case("$setModel, a coefficient and a tree count", {
  live <- make(factorModel(""))
  saved <- live$model
  same(refusal(live$setModel(live$model)), "taken", "its own model")
  live$setForest(2L, amplitude = fixed())
  expect_true(
    grepl(
      "cannot change forest 2's 'amplitude': change it with $setForest(\"z\", amplitude = normal())",
      refusal(live$setModel(saved)),
      fixed = TRUE
    ),
    info = refusal(live$setModel(saved))
  )
  drawn <- live$model
  drawn@forests[[2L]]@held.amplitude <- FALSE
  drawn@forests[[2L]]@amplitude.df <- 3
  expect_true(grepl(
    "amplitude = student(3))",
    refusal(live$setModel(drawn)),
    fixed = TRUE
  ))
  # a tree count the record states is structure, refused on every sampler; one
  # it does not state is a default, and the sampler keeps its own
  stating <- function(record, count) {
    record@n.trees <- count
    record@stated <- union(record@stated, "n.trees")
    record
  }
  counted <- live$model
  counted@forests[[2L]] <- stating(counted@forests[[2L]], 3L)
  expect_true(grepl(
    "forest 2's tree count",
    refusal(live$setModel(counted)),
    fixed = TRUE
  ))
  counted@forests[[2L]]@stated <- setdiff(
    counted@forests[[2L]]@stated,
    "n.trees"
  )
  count <- live$model@forests[[2L]]@n.trees
  same(refusal(live$setModel(counted)), "taken", "an unstated count")
  same(live$model@forests[[2L]]@n.trees, count, "the sampler's count kept")
  lone <- make("y ~ x1 + x2")
  oneCount <- lone$model
  oneCount@forests[[1L]] <- stating(oneCount@forests[[1L]], 7L)
  expect_true(grepl(
    "forest 1's tree count",
    refusal(lone$setModel(oneCount)),
    fixed = TRUE
  ))
  withVariance <- dbarts(
    y ~ x1 + x2,
    d,
    n.trees = 24L,
    variance = ~x1,
    control = ctl()
  )
  varianceCount <- withVariance$model
  varianceCount@variance <- stating(varianceCount@variance, 3L)
  expect_true(grepl(
    "the variance forest's tree count",
    refusal(withVariance$setModel(varianceCount)),
    fixed = TRUE
  ))
  moved <- live$model
  moved@forests[[2L]]@moves <- list(
    proposal.probs = c(
      birth_death = 0.5,
      swap = 0,
      change = 0.5,
      perturb = 0,
      rule_gibbs = 0,
      birth = 0.5
    )
  )
  same(
    refusal(live$setModel(moved)),
    refusal(live$setModel(live$model)),
    "moves still taken"
  )
})
