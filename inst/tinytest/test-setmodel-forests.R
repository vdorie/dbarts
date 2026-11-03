# The fitting function's level on a live sampler of several forests or a
# multinomial one: $setModel, $setTreePrior and $setLeafPrior. The main oracle
# is the twin: a write before any run leaves the records and, on two chains,
# the next draws of a sampler created with what the write leaves. Every case
# is named, an error in one is a failed expectation carrying that name, and
# the file runs to its end.

forest <- dbartsForests$forest
varianceForest <- dbartsForests$varianceForest
student <- dbartsPriors$student
normal <- dbartsPriors$normal
chi <- dbartsPriors$chi

set.seed(23L)
n <- 60L
z <- rbinom(n, 1L, 0.5)
dose <- runif(n, 1, 5)
x <- data.frame(x1 = runif(n), x2 = runif(n))
y <- x$x1 + z * x$x2 + 0.2 * dose + rnorm(n, 0, 0.3)
yBinary <- as.double(y > stats::median(y))
category <- factor(sample(c("a", "b", "c"), n, TRUE))

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

# a sampler of the forests `entries`, each a list of forest()'s arguments,
# with the fitting function's `...`
build <- function(entries, ..., response = y, n.trees = 8L) {
  dbarts(
    x,
    response,
    forests = lapply(entries, function(arguments) do.call(forest, arguments)),
    n.trees = n.trees,
    control = ctl(),
    ...
  )
}
buildMultinomial <- function(stated = list(), ...) {
  dbarts(
    x,
    category,
    family = "multinomial",
    forests = list(do.call(forest, stated)),
    n.trees = 8L,
    control = ctl(),
    ...
  )
}
buildVariance <- function(mean = list(), variance = list(), ...) {
  dbarts(
    x,
    y,
    forests = list(do.call(forest, mean)),
    variance = do.call(varianceForest, variance),
    n.trees = 8L,
    control = ctl(),
    ...
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
expectFound <- function(text, pattern, ...) {
  expect_true(
    grepl(pattern, text, fixed = TRUE),
    info = paste0(paste(caseName, ...), ": ", text)
  )
}
# the sampler, after a write, records and draws as `twin`
expectTwin <- function(live, twin, ...) {
  same(modelOf(live), modelOf(twin), ..., "records")
  same(draws(live), draws(twin), ..., "draws")
}
# a refused or inert write leaves the records and the next draws
expectUntouched <- function(sampler, untouched, ...) {
  same(modelOf(sampler), modelOf(untouched), ..., "records")
  same(draws(sampler), draws(untouched), ..., "draws")
}
# a copy and a reload continue as the live sampler does: bitwise alike, and
# the live sampler to the reassociation of its fits a state install leaves
expectContinues <- function(live, ...) {
  copied <- live$copy()
  loaded <- reloaded(live)
  liveDraws <- draws(live)
  copyDraws <- draws(copied)
  same(draws(loaded), copyDraws, ..., "reload as copy")
  expect_equal(
    liveDraws,
    copyDraws,
    tolerance = 1e-10,
    info = paste(caseName, ..., "copy as live")
  )
}
depthOf <- function(sampler, index, name) {
  sampler$getForest(index)[[name]]
}

threeSharers <- list(mu1 = list(), mu2 = list(power = 3), mu3 = list())
oneCandidate <- function(...) list(mu = list(...), tau = factorBasis)

## --- the level's depth on three sharers, one stating power -------------------
case("three sharers: $setTreePrior's power", {
  live <- build(threeSharers)
  expect_null(live$setTreePrior(power = 2.5))
  twin <- build(threeSharers, power = 2.5)
  same(depthOf(live, 1L, "power"), 2.5, "forest 1 follows")
  same(depthOf(live, 3L, "power"), 2.5, "forest 3 follows")
  same(depthOf(live, 2L, "power"), 3, "forest 2 keeps its own")
  expectTwin(live, twin)
  expectContinues(live)
})
case("three sharers: $setModel's power", {
  live <- build(threeSharers)
  model <- live$model
  model@tree.prior@power <- 2.5
  live$setModel(model)
  twin <- build(threeSharers, power = 2.5)
  expectTwin(live, twin)
  expectContinues(live)
})
case("three sharers: $setTreePrior and $setModel agree", {
  viaDoor <- build(threeSharers)
  viaModel <- build(threeSharers)
  viaDoor$setTreePrior(base = 0.8, power = 2.5)
  model <- viaModel$model
  model@tree.prior@power <- 2.5
  model@tree.prior@base <- 0.8
  viaModel$setModel(model)
  same(modelOf(viaDoor), modelOf(viaModel), "records")
  same(draws(viaDoor), draws(viaModel), "draws")
})
case("three sharers: a level change leaves the forests following", {
  live <- build(threeSharers)
  live$setTreePrior(power = 2.5)
  live$setTreePrior(power = 1.5)
  same(depthOf(live, 1L, "power"), 1.5, "forest 1 still follows")
  same(depthOf(live, 3L, "power"), 1.5, "forest 3 still follows")
  same("power" %in% live$model@forests[[1L]]@stated, FALSE, "not stated")
  twin <- build(threeSharers, power = 1.5)
  expectTwin(live, twin)
})
case("three sharers: mid-run, a copy and a reload continue", {
  live <- build(threeSharers)
  invisible(draws(live))
  live$setTreePrior(power = 2.5)
  expectContinues(live)
})
case("every candidate stating its own: refused by name", {
  entries <- list(
    mu1 = list(power = 3),
    mu2 = list(power = 2.5),
    tau = factorBasis
  )
  live <- build(entries)
  untouched <- build(entries)
  message <- refusal(live$setTreePrior(power = 1.5))
  expectFound(message, "reaches no forest")
  expectFound(message, "forest \"mu1\" and forest \"mu2\"")
  expectFound(message, "$setForest(\"mu1\", power = 1.5)")
  model <- live$model
  model@tree.prior@power <- 1.5
  message <- refusal(live$setModel(model))
  expectFound(message, "$setForest")
  expectUntouched(live, untouched)
})
case("$setTreePrior: NULL and a left-out name change nothing", {
  live <- build(threeSharers, power = 1.5, base = 0.8)
  untouched <- build(threeSharers, power = 1.5, base = 0.8)
  expect_null(live$setTreePrior(power = NULL))
  live$setTreePrior(base = NULL)
  live$setTreePrior()
  live$setTreePrior(power = 1.5)
  expectUntouched(live, untouched)
})
case("three sharers: the residual prior", {
  live <- build(threeSharers)
  model <- live$model
  model@resid.prior <- dbartsPriors$chisq(10, 0.5)
  live$setModel(model)
  twin <- build(
    threeSharers,
    family = gaussian(sigma = dbartsPriors$chisq(10, 0.5))
  )
  # the twin states its family, where the live model states the prior alone
  same(live$model@resid.prior, twin$model@resid.prior, "the prior")
  same(draws(live), draws(twin), "draws")
})
case("$setTreePrior: a value that is not one is refused", {
  live <- build(threeSharers)
  expectFound(refusal(live$setTreePrior(base = 1)), "'base' must be")
  expectFound(refusal(live$setTreePrior(power = -1)), "'power' must be")
})

## --- the later write wins (dec-B523) -----------------------------------------
case("a lone forest stating power", {
  live <- dbarts(
    x,
    y,
    forests = list(forest(power = 3)),
    n.trees = 8L,
    control = ctl()
  )
  live$setTreePrior(power = 2.5)
  same(live$getForest(1L)$power, 2.5, "written")
  twin <- dbarts(x, y, n.trees = 8L, control = ctl(), power = 2.5)
  same(draws(live), draws(twin), "draws")
  live$setForest(1L, power = 3)
  same(live$getForest(1L)$power, 3, "written back")
})
case("the one candidate stating power", {
  live <- build(oneCandidate(power = 3))
  live$setTreePrior(power = 2.5)
  same(depthOf(live, 1L, "power"), 2.5, "written")
  twin <- build(oneCandidate(), power = 2.5)
  same(draws(live), draws(twin), "draws")
  live$setForest(1L, power = 3)
  same(depthOf(live, 1L, "power"), 3, "written back")
  model <- live$model
  model@tree.prior@power <- 2
  live$setModel(model)
  same(depthOf(live, 1L, "power"), 2, "$setModel writes it")
})
case("the one candidate stating sd: $setLeafPrior writes it", {
  live <- build(oneCandidate(sd = 1.2))
  live$setLeafPrior(normal(k = 3))
  same(live$model@forests[[1L]]@sd, NA_real_, "the record's sd goes")
  same(live$getLeafPrior(1L)$sd.stated, FALSE, "the reader")
  twin <- build(oneCandidate(), leaf.prior = normal(k = 3))
  expectTwin(live, twin)
  live$setForest(1L, sd = 1.2)
  same(live$model@forests[[1L]]@sd, 1.2, "written back")
  same(live$getLeafPrior(1L)$sd, 1.2, "the reader after")
})
case("probit forest(sd = 2): normal() returns k to drawn", {
  live <- build(oneCandidate(sd = 2), response = yBinary)
  spread <- live$getLeafPrior(1L)$sd
  live$setLeafPrior(normal())
  same(live$model@forests[[1L]]@sd, NA_real_, "the record's sd goes")
  prior <- live$getLeafPrior(1L)
  same(is(prior$leaf.prior@k, "dbartsChiHyperprior"), TRUE, "drawn")
  # from the spread in force, call 8's exemption
  expect_equal(
    prior$k.scale / live$getK(1L),
    rep(spread, 2L),
    info = paste(caseName, "the start")
  )
  expectContinues(live)
})

## --- the one candidate that draws a coefficient (dec-B540, dec-B550) ---------
dfEntries <- oneCandidate(amplitude = student(3))
for (written in list(
  list(name = "sd", spec = normal(sd = 0.9)),
  list(name = "k", spec = normal(k = 3)),
  list(name = "chi", spec = normal(k = chi(3, 2)))
)) {
  case(paste("the df candidate:", written$name), {
    live <- build(dfEntries)
    created <- build(dfEntries, leaf.prior = written$spec)
    target <- created$getLeafPrior(1L)$sd
    coefficient <- live$getForestAmplitudes()
    k <- live$getK(1L)
    live$setLeafPrior(written$spec)
    prior <- live$getLeafPrior(1L)
    same(prior$sd, target, "the sd creation's reading gives")
    expect_equal(
      unname(prior$k.scale),
      unname(target),
      tolerance = 1e-12,
      info = paste(caseName, "k.scale")
    )
    same(prior$amplitude, student(3), "the df kept")
    same(live$getK(1L), k, "no k")
    same(live$getForestAmplitudes(), coefficient, "the coefficient unmoved")
    same(live$model@forests[[1L]]@sd, NA_real_, "no sd stated")
    stated <- build(dfEntries)
    stated$setForest(1L, sd = target)
    same(draws(live), draws(stated), "draws as $setForest(1, sd = )")
    expectContinues(live)
  })
}
case("the df candidate stating sd: the statement goes", {
  live <- build(oneCandidate(amplitude = student(3), sd = 1.2))
  live$setLeafPrior(normal(k = 3))
  same(live$model@forests[[1L]]@sd, NA_real_, "the record's sd")
  same(live$getLeafPrior(1L)$sd.stated, FALSE, "the reader")
  target <- build(dfEntries, leaf.prior = normal(k = 3))$getLeafPrior(1L)$sd
  same(live$getLeafPrior(1L)$sd, target, "the sd")
})
case("the df candidate: after $setResponse, creation's scale", {
  live <- build(dfEntries)
  target <- live$getLeafPrior(1L)$k.scale
  live$setResponse(10 * y)
  live$setLeafPrior(normal(k = 3))
  expect_equal(
    live$getLeafPrior(1L)$sd,
    target * 2 / 3,
    info = paste(caseName, "the scale creation set over k")
  )
  live$setLeafPrior(normal(k = 4))
  expect_equal(
    live$getLeafPrior(1L)$sd,
    target * 2 / 4,
    info = paste(caseName, "a second write")
  )
})

## --- the group's leaf prior ---------------------------------------------------
sharerEntries <- list(mu1 = list(), mu2 = list(), tau = factorBasis)
for (written in list(
  list(name = "held", spec = normal(k = 3)),
  list(name = "sd", spec = normal(sd = 0.8))
)) {
  case(paste("the group:", written$name), {
    live <- build(sharerEntries)
    live$setLeafPrior(written$spec)
    twin <- build(sharerEntries, leaf.prior = written$spec)
    k <- live$getK()
    same(k[1L, ], k[2L, ], "every sharer's k equal")
    expectTwin(live, twin)
    viaModel <- build(sharerEntries)
    model <- viaModel$model
    model@leaf.prior <- twin$model@leaf.prior
    model@leaf.hyperprior <- twin$model@leaf.hyperprior
    model@prior.scale <- twin$model@prior.scale
    viaModel$setModel(model)
    expectTwin(
      viaModel,
      build(sharerEntries, leaf.prior = written$spec),
      "$setModel"
    )
  })
}
case("the group: drawn from the k in force", {
  live <- build(sharerEntries)
  k <- live$getK()
  live$setLeafPrior(normal(k = chi(3, 2)))
  same(live$getK(), k, "the k in force")
  expectContinues(live)
})
case("the group: unequal counts keep their own leaf scales", {
  entries <- list(
    mu1 = list(n.trees = 38L),
    mu2 = list(n.trees = 37L),
    tau = factorBasis
  )
  live <- build(entries, n.trees = NULL)
  live$setLeafPrior(normal(k = 3))
  twin <- build(entries, leaf.prior = normal(k = 3), n.trees = NULL)
  same(live$getLeafPrior(1L)$sd, twin$getLeafPrior(1L)$sd, "forest 1")
  same(live$getLeafPrior(2L)$sd, twin$getLeafPrior(2L)$sd, "forest 2")
  same(draws(live), draws(twin), "draws")
})
case("the group: a forest that states sd unmoved", {
  entries <- list(mu1 = list(), mu2 = list(), mu3 = list(sd = 0.5))
  live <- build(entries)
  live$setLeafPrior(normal(k = 3))
  same(live$getLeafPrior(3L)$sd, 0.5, "its own")
  expectTwin(live, build(entries, leaf.prior = normal(k = 3)))
})
case("the lone sharer: a held 2 restated is inert", {
  live <- build(oneCandidate())
  untouched <- build(oneCandidate())
  live$setLeafPrior(normal(k = 2))
  expectUntouched(live, untouched)
})
case("the lone sharer: a held 3 draws as creation, and back", {
  live <- build(oneCandidate())
  live$setLeafPrior(normal(k = 3))
  expectTwin(live, build(oneCandidate(), leaf.prior = normal(k = 3)))
  back <- build(oneCandidate())
  back$setLeafPrior(normal(k = 3))
  back$setLeafPrior(normal(k = 2))
  created <- build(oneCandidate())
  # the default restated reads as a number, where creation's reads as written
  same(
    back$model@leaf.hyperprior,
    created$model@leaf.hyperprior,
    "the hyperprior"
  )
  same(back$getK(), created$getK(), "getK")
  same(draws(back), draws(created), "draws")
})
case("several candidates all stating sd: the default only", {
  entries <- list(mu1 = list(sd = 0.5), mu2 = list(sd = 0.6))
  live <- build(entries)
  untouched <- build(entries)
  live$setLeafPrior(normal(k = 2))
  live$setLeafPrior(normal())
  message <- refusal(live$setLeafPrior(normal(k = 3)))
  expectFound(message, "no forest of this sampler takes")
  expectFound(message, "$setForest(\"mu1\", sd = v)")
  expectUntouched(live, untouched)
})
case("every forest with a basis: the default only", {
  entries <- list(tau = factorBasis, gamma = list(basis = ~dose, n.trees = 6L))
  live <- dbarts(
    x,
    y,
    forests = lapply(entries, function(a) do.call(forest, a)),
    control = ctl()
  )
  message <- refusal(live$setLeafPrior(normal(k = 3)))
  expectFound(message, "every forest having a basis")
  expectFound(refusal(live$setTreePrior(power = 3)), "every forest of")
})

## --- a multinomial sampler ----------------------------------------------------
case("multinomial: the level's depth and k", {
  live <- buildMultinomial()
  live$setTreePrior(power = 2.5)
  live$setLeafPrior(normal(k = 3))
  twin <- buildMultinomial(power = 2.5, leaf.prior = normal(k = 3))
  expectTwin(live, twin)
})
case("multinomial: its one entry stating power is written", {
  live <- buildMultinomial(list(power = 3))
  live$setTreePrior(power = 2.5)
  same(live$model@tree.prior@power, 2.5, "written")
  same(draws(live), draws(buildMultinomial(power = 2.5)), "draws")
})

## --- refused by name, the sampler untouched ---------------------------------
case("$setModel: a forest's own difference", {
  live <- build(threeSharers)
  untouched <- build(threeSharers)
  model <- live$model
  model@forests[[2L]]@power <- 1.5
  message <- refusal(live$setModel(model))
  expectFound(message, "\"mu2\"'s 'power'")
  expectFound(message, "$setForest(\"mu2\", power = 1.5)")
  model <- live$model
  model@forests[[1L]]@stated <- c(model@forests[[1L]]@stated, "power")
  expectFound(refusal(live$setModel(model)), "$setForest(\"mu1\"")
  model <- live$model
  model@forests[[1L]]@sd <- 0.4
  expectFound(refusal(live$setModel(model)), "$setForest(\"mu1\", sd = 0.4)")
  expectUntouched(live, untouched)
})
case("split.probs on several forests", {
  live <- build(threeSharers)
  untouched <- build(threeSharers)
  expectFound(
    refusal(live$setTreePrior(split.probs = c(0.7, 0.3))),
    "does not support 'split.probs'"
  )
  model <- live$model
  model@tree.prior@splitProbabilities <- c(0.7, 0.3)
  expectFound(refusal(live$setModel(model)), "does not support 'split.probs'")
  live$setTreePrior(split.probs = c(0.5, 0.5))
  expectUntouched(live, untouched)
})

## --- inert ---------------------------------------------------------------------
inertKinds <- list(
  sharers = function() build(threeSharers),
  candidate = function() build(oneCandidate(sd = 1.2, power = 3)),
  student = function() build(dfEntries),
  bases = function() {
    dbarts(
      x,
      y,
      forests = list(
        forest(basis = ~ factor(z), n.trees = 6L),
        forest(basis = ~dose, n.trees = 6L)
      ),
      control = ctl()
    )
  },
  multinomial = function() buildMultinomial(list(power = 3)),
  variance = function() buildVariance(variance = list(power = 3))
)
for (kind in names(inertKinds)) {
  case(paste("$setModel(s$model) inert:", kind), {
    live <- inertKinds[[kind]]()
    untouched <- inertKinds[[kind]]()
    live$setModel(live$model)
    expectUntouched(live, untouched)
  })
}
case("$setModel(s$model) inert, then $setForest as before", {
  live <- build(threeSharers)
  untouched <- build(threeSharers)
  live$setModel(live$model)
  live$setForest(1L, power = 1.5)
  untouched$setForest(1L, power = 1.5)
  expectUntouched(live, untouched)
})

## --- carried pins --------------------------------------------------------------
case("student(3): a new depth leaves the coefficient", {
  live <- build(dfEntries)
  invisible(draws(live))
  amplitude <- live$getForest(1L)$amplitude
  coefficient <- live$getForestAmplitudes()
  live$setTreePrior(power = 2.5)
  same(live$getForest(1L)$amplitude, amplitude, "the reader's amplitude")
  same(live$getForestAmplitudes(), coefficient, "the coefficient")
  expectContinues(live)
})
case("a basis scale restated, then $setModel(s$model)", {
  entries <- list(mu = list(), gamma = list(basis = ~dose, n.trees = 6L))
  live <- build(entries)
  invisible(draws(live))
  live$setForestBasis(2L, dose * 2, updateBasisScale = TRUE)
  sd <- live$getLeafPrior(2L)$sd
  live$setModel(live$model)
  same(live$getLeafPrior(2L)$sd, sd, "the reader's sd")
  expectContinues(live)
})

## --- the variance forest (dec-B524) -------------------------------------------
for (stated in list(list(), list(power = 3))) {
  case(paste("the variance forest stays:", length(stated)), {
    live <- buildVariance(variance = stated)
    before <- live$getVarianceForest()
    live$setTreePrior(power = 2.5)
    same(live$getVarianceForest(), before, "$setTreePrior")
    model <- live$model
    model@tree.prior@power <- 1.5
    live$setModel(model)
    same(live$getVarianceForest(), before, "$setModel")
    twin <- buildVariance(
      mean = list(),
      variance = list(power = before$power),
      power = 1.5
    )
    same(draws(live), draws(twin), "draws")
    live$setVarianceForest(power = 2)
    same(live$getVarianceForest()$power, 2, "its door")
  })
}
case("the variance forest: a model's record is not a door", {
  # no flag records what a variance forest states, and a spec's record takes
  # the fitting function's depth, so the sampler's record stays
  live <- buildVariance(variance = list(power = 3))
  untouched <- buildVariance(variance = list(power = 3))
  model <- live$model
  model@variance@power <- 1.5
  live$setModel(model)
  expectUntouched(live, untouched)
})

## --- a size the leaf prior cannot give ------------------------------------------
case("the df candidate: no size, refused by name", {
  live <- build(dfEntries)
  untouched <- build(dfEntries)
  message <- refusal(live$setLeafPrior(normal(k = Inf)))
  expectFound(message, "gives no positive finite size")
  expectFound(message, "$setForest(\"mu\", sd = v)")
  expectUntouched(live, untouched)
})
