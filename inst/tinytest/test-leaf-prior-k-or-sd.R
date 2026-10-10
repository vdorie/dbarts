# The leaf prior is named by k (relative to the data's scale) or by sd (on
# the family's scale), never both. The model encodes a named sd beside a
# reference k of 2: sd = x is prior.scale 2x with k fixed at 2, sd =
# invchi(df, c) is prior.scale 2c with k ~ chi(df, 2), and invchi(df, 0) is
# chi(df, Inf) with no prior.scale. The bridge divides the pair back to the
# sd, which the engine states as k against the data's scale. The pins below
# are the model's encoding of each spelling, held bitwise.

hex <- function(value) sprintf("%a", value)
inputsOf <- function(model) {
  hyperprior <- model@leaf.hyperprior
  law <- if (is(hyperprior, "dbartsChiHyperprior")) {
    c(
      "chi",
      hex(hyperprior@degreesOfFreedom),
      if (is.finite(hyperprior@scale)) hex(hyperprior@scale) else "Inf"
    )
  } else {
    c("fixed", hex(hyperprior@k))
  }
  c(if (is.na(model@prior.scale)) "NA" else hex(model@prior.scale), law)
}

set.seed(1L)
n <- 120L
x <- matrix(runif(n * 2L), n, 2L)
colnames(x) <- c("x1", "x2")
yg <- 3 * x[, 1L] + rnorm(n)
yb <- rbinom(n, 1L, pnorm(1.5 * (x[, 1L] - 0.5)))
pinControl <- function(n.chains = 2L) {
  dbartsControl(
    n.chains = n.chains,
    n.threads = 1L,
    n.samples = 15L,
    n.burn = 10L,
    seed = 7L,
    n.trees = 20L
  )
}
modelOf <- function(y, leafPrior, family, ...) {
  eval(bquote(
    dbarts(
      x,
      .(y),
      leaf.prior = .(leafPrior),
      family = .(family),
      control = pinControl(),
      ...
    )
  ))$model
}

# --- translation pins: each new spelling's engine inputs are the recorded
# inputs of the spelling it replaces (normal(k = 2, scale = 1.3),
# normal(scale = 1.3) on probit, normal(k = chi(1.25, 2), scale = 1.3)) ---
expect_identical(
  inputsOf(modelOf(yg, quote(normal(sd = 0.65)), "gaussian")),
  c("0x1.4cccccccccccdp+0", "fixed", "0x1p+1")
)
expect_identical(
  inputsOf(modelOf(yb, quote(normal(sd = invchi(1.5, 0.65))), "probit")),
  c("0x1.4cccccccccccdp+0", "chi", "0x1.8p+0", "0x1p+1")
)
expect_identical(
  inputsOf(modelOf(yb, quote(normal(sd = invchi(1.25, 0.65))), "probit")),
  c("0x1.4cccccccccccdp+0", "chi", "0x1.4p+0", "0x1p+1")
)
# normal(k = 2, sd = 0.7) reached the engine as 2 * 0.7 at a fixed 2
expect_identical(
  inputsOf(modelOf(yg, quote(normal(sd = 0.7)), "gaussian")),
  c("0x1.6666666666666p+0", "fixed", "0x1p+1")
)
# the improper limit has no anchor: it is the k form chi(df, Inf)
expect_identical(
  inputsOf(modelOf(yb, quote(normal(sd = invchi(1.25, 0))), "probit")),
  c("NA", "chi", "0x1.4p+0", "Inf")
)
expect_identical(
  inputsOf(modelOf(yb, quote(normal(k = chi(1.25, Inf))), "probit")),
  c("NA", "chi", "0x1.4p+0", "Inf")
)
# the k forms are unchanged
expect_identical(
  inputsOf(modelOf(yb, quote(normal()), "probit")),
  c("NA", "chi", "0x1.8p+0", "0x1p+1")
)
expect_identical(
  inputsOf(modelOf(yg, quote(normal(k = 4)), "gaussian")),
  c("NA", "fixed", "0x1p+2")
)

# --- live bitwise oracles, two chains, every chain compared: a k spelling and
# the sd spelling of the same prior draw identically ---
drawsOf <- function(y, leafPrior, family = "gaussian", ...) {
  sampler <- eval(bquote(
    dbarts(
      x,
      .(y),
      leaf.prior = .(leafPrior),
      family = .(family),
      control = pinControl(),
      ...
    )
  ))
  run <- sampler$run()
  list(train = run$train, k = run$k)
}
expectSameDraws <- function(a, b, info) {
  expect_identical(a$train, b$train, info = info)
  expect_identical(a$k, b$k, info = info)
}
expectSameDraws(
  drawsOf(yb, quote(normal()), "probit"),
  drawsOf(yb, quote(normal(sd = invchi(1.5, 1.5))), "probit"),
  "probit default"
)
expectSameDraws(
  drawsOf(yb, quote(normal()), "logistic"),
  drawsOf(yb, quote(normal(sd = invchi(1.5, pi * sqrt(3) / 2))), "logistic"),
  "logistic default"
)
expectSameDraws(
  drawsOf(yg, quote(normal())),
  drawsOf(yg, bquote(normal(sd = .(diff(range(yg)) / 4)))),
  "gaussian default"
)
expectSameDraws(
  drawsOf(yg, quote(normal(k = 4))),
  drawsOf(yg, bquote(normal(sd = .(diff(range(yg)) / 8)))),
  "gaussian k = 4"
)
expectSameDraws(
  drawsOf(yb, quote(normal(k = chi(1.25, Inf))), "probit"),
  drawsOf(yb, quote(normal(sd = invchi(1.25, 0))), "probit"),
  "improper limit"
)
expectSameDraws(
  drawsOf(yg, quote(linear("x1", k = 4))),
  drawsOf(yg, bquote(linear("x1", sd = .(diff(range(yg)) / 8)))),
  "linear"
)
expectSameDraws(
  drawsOf(yg, quote(gp("x1", k = 4))),
  drawsOf(yg, bquote(gp("x1", sd = .(diff(range(yg)) / 8)))),
  "gp"
)
expectSameDraws(
  drawsOf(yg, quote(normal(k = 4)), monotone = c(x1 = 1)),
  drawsOf(
    yg,
    bquote(normal(sd = .(diff(range(yg)) / 8))),
    monotone = c(x1 = 1)
  ),
  "monotone"
)
# hazard: one binary model on the person-period expansion, whose sd is in its
# link's latent units
hazardFit <- function(leafPrior) {
  eval(bquote(bart(
    x,
    cbind(sample.int(4L, n, TRUE), rbinom(n, 1L, 0.5)),
    family = "hazard",
    leaf.prior = .(leafPrior),
    n.chains = 2L,
    n.threads = 1L,
    n.samples = 10L,
    n.burn = 5L,
    n.trees = 20L,
    seed = 7L,
    verbose = FALSE
  )))
}
set.seed(2L)
hazardK <- hazardFit(quote(normal(k = 4)))
set.seed(2L)
hazardSd <- hazardFit(quote(normal(sd = 0.75)))
expect_identical(hazardK$yhat.train, hazardSd$yhat.train)

# --- the refusals ---
priors <- dbarts:::dbartsPriors
for (constructor in c("normal", "linear", "gp")) {
  args <- if (constructor == "normal") list() else list(columns = 1L)
  expect_error(
    do.call(priors[[constructor]], c(args, list(k = 2, sd = 1))),
    "give either 'k'.*or 'sd'.*not both",
    info = constructor
  )
  expect_error(
    do.call(priors[[constructor]], c(args, list(scale = 1))),
    "unused argument",
    info = constructor
  )
}
expect_error(priors$normal(sd = priors$chi(1.5, 2)), "k = chi\\(\\)")
expect_error(priors$normal(k = priors$invchi(1.5, 1)), "sd = invchi\\(\\)")
expect_error(priors$normal(sd = "1"), "takes no string form")
expect_error(priors$invchi(1.5), "requires 'scale'")
expect_error(priors$invchi(0, 1), "'df' must be a single positive finite")
expect_error(priors$invchi(-1, 1), "'df' must be a single positive finite")
expect_error(priors$invchi(Inf, 1), "'df' must be a single positive finite")
expect_error(
  priors$invchi(NA_real_, 1),
  "'df' must be a single positive finite"
)
expect_error(priors$invchi(1.5, -1), "'scale' must be a single non-negative")
expect_error(priors$invchi(1.5, Inf), "'scale' must be a single non-negative")
expect_error(priors$invchi(1.5, NA_real_), "'scale' must be a single non-neg")
expect_error(priors$invchi(1.5, c(1, 2)), "single numbers")
# an sd hyperprior under a monotone constraint, as a k hyperprior is
expect_error(
  dbarts(
    x,
    yg,
    leaf.prior = normal(sd = invchi(1.5, 1)),
    monotone = c(x1 = 1),
    control = pinControl()
  ),
  "'sd' hyperprior is not supported under a monotone"
)
# a numeric sd under monotone is accepted (pinned above)
# a causal forest refuses a named sd, number or hyperprior, pointing at
# forest(sd = )
z <- rbinom(n, 1L, 0.5)
for (namedSd in list(0.5, priors$invchi(1.5, 0.5))) {
  expect_error(
    dbarts(
      x,
      yg,
      forests = list(forest(), forest(basis = ~ factor(z))),
      leaf.prior = normal(sd = namedSd),
      control = pinControl()
    ),
    "leaf-prior 'sd'.*forest\\(sd = \\)"
  )
}
# a hurdle fit refuses a named sd before either half is built
yHurdle <- ifelse(runif(n) < 0.3, 0, exp(yg / 3))
for (leafPrior in list(
  quote(normal(sd = 0.5)),
  quote(normal(sd = invchi(1.5, 0.5)))
)) {
  expect_error(
    eval(bquote(bart(
      x,
      yHurdle,
      family = "hurdle.lognormal",
      leaf.prior = .(leafPrior),
      n.samples = 5L,
      n.burn = 5L,
      n.chains = 1L,
      verbose = FALSE
    ))),
    "one 'sd' cannot state both"
  )
}

# --- chi()'s first argument is df; degreesOfFreedom is a tombstone that
# warns once per session and uses the value ---
onceState <- dbarts:::onceWarnState
onceState[["tombstone.degreesOfFreedom.chi"]] <- NULL
warnings <- character()
viaOldName <- withCallingHandlers(
  priors$chi(degreesOfFreedom = 3, scale = 2),
  warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
expect_equal(length(warnings), 1L)
expect_true(grepl("'degreesOfFreedom' is now 'df'", warnings))
expect_identical(viaOldName, priors$chi(df = 3, scale = 2))
expect_silent(priors$chi(degreesOfFreedom = 3))
expect_identical(priors$chi(3, 2), priors$chi(df = 3, scale = 2))
expect_error(priors$chi(df = 1, degreesOfFreedom = 1), "supply one")

# --- xbart: an sd grid beside k, exclusive with it ---
xbartArgs <- list(
  formula = x,
  data = yg,
  n.samples = 10L,
  n.reps = 2L,
  n.burn = c(10L, 5L),
  n.trees = 10L,
  n.threads = 1L,
  seed = 3L,
  verbose = FALSE
)
sdGrid <- do.call(xbart, c(xbartArgs, list(sd = c(0.5, 1, 0.25))))
expect_identical(names(dimnames(sdGrid)), c("rep", "sd"))
expect_identical(dimnames(sdGrid)$sd, c("0.5", "1", "0.25"))
# each cell is the one-call fit of that sd: the grid is swept most shrunk
# first with warm starts, so compare against a grid in that same order
sweptOrder <- do.call(xbart, c(xbartArgs, list(sd = c(0.25, 0.5, 1))))
expect_identical(unname(sdGrid[, c(3L, 1L, 2L)]), unname(sweptOrder))
# on a binary family, whose anchor is a constant (3 on probit), a grid of
# fixed cells is the same sweep in either spelling; a modelled cell is not,
# since the two start their chains at different spreads
probitArgs <- xbartArgs
probitArgs$data <- yb
expect_identical(
  unname(do.call(xbart, c(probitArgs, list(sd = c(1.5, 0.75, 3))))),
  unname(do.call(xbart, c(probitArgs, list(k = c(2, 4, 1)))))
)
hyperGrid <- do.call(
  xbart,
  c(xbartArgs, list(sd = list(0.5, priors$invchi(1.5, 0.5))))
)
expect_identical(dimnames(hyperGrid)$sd, c("0.5", "invchi(1.5, 0.5)"))
expect_error(
  do.call(xbart, c(xbartArgs, list(k = 2, sd = 1))),
  "give either 'k'.*or 'sd'.*as the grid"
)
expect_error(
  do.call(xbart, c(xbartArgs, list(k = 2, leaf.prior = quote(normal(sd = 1))))),
  "the leaf prior's 'sd' and the 'k' grid both state the spread"
)
expect_error(
  do.call(
    xbart,
    c(xbartArgs, list(sd = 2, leaf.prior = quote(normal(sd = 1))))
  ),
  "the leaf prior's 'sd' and the 'sd' grid both state the spread"
)
expect_error(
  do.call(xbart, c(xbartArgs, list(sd = 2, leaf.prior = quote(normal(k = 3))))),
  "the leaf prior's 'k' and the 'sd' grid both state the spread"
)
expect_error(
  do.call(xbart, c(xbartArgs, list(sd = list(priors$chi(1.5, 2))))),
  "k = chi\\(\\)"
)
# a named sd in the leaf prior is a one-cell sd axis, the same cell as the
# grid's
oneCell <- do.call(
  xbart,
  c(xbartArgs, list(leaf.prior = quote(normal(sd = 0.5)), drop = FALSE))
)
expect_identical(names(dimnames(oneCell))[3L], "sd")
expect_identical(
  as.vector(oneCell),
  as.vector(do.call(xbart, c(xbartArgs, list(sd = 0.5))))
)
# the per-sd recipe: two calls with the same seed and different sd use the
# same folds, so one call per sd reproduces the grid's cells; its first cell
# is fresh, which the grid's most shrunk cell is too
expect_identical(
  unname(sweptOrder[, 1L]),
  as.vector(do.call(xbart, c(xbartArgs, list(sd = 0.25))))
)

# --- a fit stores the sampler's k whatever terms the leaf prior was named in,
# and the leaf prior it ran under; the spread is the anchor over k ---
fitOf <- function(leafPrior) {
  eval(bquote(bart(
    x,
    yb,
    leaf.prior = .(leafPrior),
    n.samples = 10L,
    n.burn = 5L,
    n.chains = 2L,
    n.threads = 1L,
    n.trees = 20L,
    seed = 11L,
    keepTrees = TRUE,
    verbose = FALSE
  )))
}
kNamed <- fitOf(quote(normal()))
sdNamed <- fitOf(quote(normal(sd = invchi(1.5, 1.5))))
# no spread channel stands in for k, drawn or burn-in
expect_identical(grep("^(first[.])?sd$", names(sdNamed)), integer())
# the same chain: both carry the raw k draws and the burn-in's, and the
# spread is the anchor over them
expect_identical(sdNamed[["k"]], kNamed[["k"]])
expect_identical(sdNamed[["first.k"]], kNamed[["first.k"]])
expect_identical(extract(sdNamed, "k"), extract(kNamed, "k"))
expect_equal(sdNamed$leaf.prior$k.scale, 3)
expect_equal(extract(sdNamed, "leaf.prior.sd"), 3 / extract(kNamed, "k"))
expect_equal(extract(kNamed, "leaf.prior.sd"), 3 / extract(kNamed, "k"))
expect_identical(sdNamed$fixed, list())
expect_error(extract(sdNamed, "sd"), "type must be in")
# summary follows the naming: leaf.prior.sd on the sd-named fit, k on the
# k-named one, and either on request
expect_true("leaf.prior.sd" %in% summary(sdNamed)$stats$variable)
expect_false("k" %in% summary(sdNamed)$stats$variable)
expect_true("k" %in% summary(kNamed)$stats$variable)
expect_false("leaf.prior.sd" %in% summary(kNamed)$stats$variable)
expect_true("k" %in% summary(sdNamed, vars = "k")$stats$variable)
expect_true(
  "first.k" %in% summary(sdNamed, vars = "first.k")$stats$variable
)
# and the reader of the fit's sampler agrees about the terms
expect_identical(
  sdNamed$fit$getLeafPrior()$leaf.prior@prior.sd,
  dbartsPriors$invchi(1.5, 1.5)
)

# a model whose prior.scale names an sd beside a chi prior with an infinite
# scale is refused; no public constructor builds one, so the model is built by
# hand
handBuilt <- dbarts(
  x,
  yg,
  leaf.prior = dbartsPriors$normal(k = dbartsPriors$chi(1.5, 2)),
  control = pinControl()
)
badModel <- handBuilt$model
badModel@leaf.hyperprior <- methods::new(
  "dbartsChiHyperprior",
  degreesOfFreedom = 1.5,
  scale = Inf
)
badModel@prior.scale <- 1.3
expect_error(
  handBuilt$setModel(badModel),
  "named prior scale requires a finite k or chi scale"
)

# a model whose prior.scale is set while its leaf prior carries no sd is read
# back with the sd the encoding names
bare <- dbarts(
  x,
  yg,
  leaf.prior = dbartsPriors$normal(sd = 0.65),
  control = pinControl()
)
bareModel <- bare$model
bareModel@leaf.prior@prior.sd <- NULL
bare$model <- bareModel
expect_identical(
  bare$getLeafPrior()$leaf.prior@prior.sd,
  bareModel@prior.scale / bareModel@leaf.hyperprior@k
)
expect_equal(bare$getLeafPrior()$leaf.prior@prior.sd, 0.65)

# and under a drawn sd, the fallback restates the chi prior's scale
bareDrawn <- dbarts(
  x,
  yg,
  leaf.prior = dbartsPriors$normal(sd = dbartsPriors$invchi(1.5, 0.65)),
  control = pinControl()
)
drawnModel <- bareDrawn$model
drawnModel@leaf.prior@prior.sd <- NULL
bareDrawn$model <- drawnModel
expect_identical(
  bareDrawn$getLeafPrior()$leaf.prior@prior.sd,
  dbartsPriors$invchi(
    1.5,
    drawnModel@prior.scale / drawnModel@leaf.hyperprior@scale
  )
)
expect_equal(
  bareDrawn$getLeafPrior()$leaf.prior@prior.sd,
  dbartsPriors$invchi(1.5, 0.65)
)
