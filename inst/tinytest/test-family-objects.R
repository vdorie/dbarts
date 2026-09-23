# The family objects (dec-B98, dec-B101): one 'family' argument carrying
# every setting only one family reads. The constructors, the exported
# vocabulary, the token/object equivalence, and that a setting written on
# the object reaches the sampler specification unchanged.

source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)
source(
  system.file("common", "countWarnings.R", package = "dbarts"),
  local = TRUE
)

x <- testData$x
y <- testData$y

# --- the vocabulary --------------------------------------------------------

# exported as one object, mirroring dbartsPriors, so no generic name
# (gaussian, probit, logistic) enters the search path to mask stats::gaussian
expect_true(is.list(dbartsFamilies))
expect_true(all(vapply(dbartsFamilies, is.function, logical(1L))))
expect_equal(
  sort(names(dbartsFamilies)),
  sort(c(
    "gaussian",
    "student",
    "probit",
    "logistic",
    "multinomial",
    "ordinal",
    "nbinom",
    "aft",
    "hazard",
    "hurdle.lognormal"
  ))
)
for (name in names(dbartsFamilies)) {
  expect_false(name %in% getNamespaceExports("dbarts"))
}

for (name in names(dbartsFamilies)) {
  family <- dbartsFamilies[[name]]()
  expect_inherits(family, "dbartsFamily")
  # hazard's link rides its token, so it is the one whose token is not the
  # constructor's own name
  expect_true(startsWith(family@token, sub("\\..*$", "", name)))
}

expect_equal(dbartsFamilies$hazard()@token, "hazard.probit")
expect_equal(dbartsFamilies$hazard(link = "logistic")@token, "hazard.logistic")
expect_equal(dbartsFamilies$hazard()@settings$max.rows, 1e7)
expect_null(dbartsFamilies$hazard()@settings$breaks)
expect_equal(dbartsFamilies$nbinom(dispersion = 3)@settings$dispersion, 3)
expect_true(is.na(dbartsFamilies$nbinom()@settings$dispersion))

# validation, by name
expect_error(dbartsFamilies$nbinom(dispersion = -1), "positive")
expect_error(dbartsFamilies$nbinom(dispersion = c(1, 2)), "single positive")
expect_error(dbartsFamilies$hazard(max.rows = 0), "positive")
expect_error(dbartsFamilies$hazard(breaks = "five"), "NULL or numeric")
expect_error(dbartsFamilies$hazard(link = "cloglog"), "should be one of")

# printing names the call a caller would write, hazard's link included
expect_stdout(show(dbartsFamilies$student(3)), "student\\(df = 3\\)")
expect_stdout(show(dbartsFamilies$hazard()), "link = \"probit\"")

# --- resolution ------------------------------------------------------------

control <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 20L,
  n.samples = 30L,
  updateState = FALSE
)

# a bare constructor call resolves in the family vocabulary, not the
# caller's frame, exactly as the prior vocabulary does inside tree.prior
expect_false(exists("student", envir = globalenv(), inherits = FALSE))
samplerNSE <- dbarts::dbarts(x, y, control = control, family = student(df = 6))
expect_equal(attr(samplerNSE$model, "resid.df"), 6)

# the same family spelled three ways gives the same specification
samplerToken <- dbarts::dbarts(x, y, control = control, family = "gaussian")
samplerObject <- dbarts::dbarts(x, y, control = control, family = gaussian())
samplerBare <- dbarts::dbarts(x, y, control = control, family = gaussian)
expect_equal(samplerToken$model@family, samplerObject$model@family)
expect_equal(samplerToken$model@family, samplerBare$model@family)

# a variable holding either spelling still resolves, in the caller's frame
heldToken <- "probit"
heldObject <- dbartsFamilies$probit()
yBinary <- as.numeric(y > median(y))
expect_equal(
  dbarts::dbarts(
    x,
    yBinary,
    control = control,
    family = heldToken
  )$model@family,
  "probit"
)
expect_equal(
  dbarts::dbarts(
    x,
    yBinary,
    control = control,
    family = heldObject
  )$model@family,
  "probit"
)

# a family this entry point does not fit is refused by name
expect_error(
  dbarts::dbarts(x, y, control = control, family = hurdle.lognormal()),
  "does not fit family"
)
expect_error(
  dbarts::dbartsSpec(
    dbarts::dbartsData(x, y),
    control = control,
    family = hazard()
  ),
  "does not fit family"
)
expect_error(
  dbarts::bart(x, y, family = 1L),
  "family name or a family object"
)

# --- forwarded through a wrapper's dots -------------------------------------

# a call forwarded through one or more wrappers' dots resolves as it would
# written directly, at every entry point and in the prior arguments too
residDf <- function(sampler) attr(sampler$model, "resid.df")
viaDots <- function(...) dbarts::dbarts(x, y, control = control, ...)
viaNested <- function(...) viaDots(...)
viaBinary <- function(...) dbarts::dbarts(x, yBinary, control = control, ...)
expect_equal(residDf(viaDots(family = student(3))), 3)
expect_equal(residDf(viaNested(family = student(4))), 4)
expect_equal(
  viaNested(family = gaussian(sigma = chisq(5, 0.5)))$model@resid.prior@df,
  5
)
expect_equal(viaBinary(family = probit)$model@family, "probit")
expect_equal(viaNested(tree.prior = cgm(power = 3))$model@tree.prior@power, 3)
viaClosure <- function(...) {
  inner <- function() dbarts::dbarts(x, y, control = control, ...)
  inner()
}
expect_equal(residDf(viaClosure(family = student(6))), 6)
expect_error(
  (function(...) {
    dbarts::dbartsSpec(dbarts::dbartsData(x, y), control = control, ...)
  })(family = hazard(breaks = 4)),
  "does not fit family"
)
expect_error(
  (function(...) dbarts::xbart(x, y, ...))(family = student(3)),
  "does not fit family"
)
viaBart <- function(...) {
  dbarts::bart(
    x,
    y,
    n.trees = 5L,
    n.samples = 5L,
    n.burn = 5L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE,
    ...
  )
}
expect_equal(
  countWarnings(
    fitViaBart <- viaBart(family = student(3), tree.prior = cgm(power = 3)),
    "warning"
  ),
  0L
)
expect_equal(unique(fitViaBart$resid.df), 3)
expect_error(
  (function(...) {
    dbarts::bart(x, factor(y > 0), family = "multinomial", ...)
  })(tree.prior = dart()),
  "DART 'tree.prior'"
)
viaXbart <- function(...) {
  dbarts::xbart(
    x,
    y,
    n.samples = 5L,
    n.reps = 1L,
    n.burn = c(3L, 2L),
    n.threads = 1L,
    ...
  )
}
expect_true(
  is.numeric(viaXbart(tree.prior = cgm(power = 3), node.prior = normal()))
)

# ordinary variables still resolve where the call was written
expect_equal(viaBinary(family = heldToken)$model@family, "probit")
expect_equal(residDf(viaNested(family = dbartsFamilies$student(5))), 5)
expect_equal(
  (function() {
    heldToken <- "logistic"
    viaBinary(family = heldToken)$model@family
  })(),
  "logistic"
)
expect_equal(
  residDf(do.call(viaDots, list(family = dbartsFamilies$student(8)))),
  8
)
# a wrapper re-entered by eval() after do.call(envir = ) has no caller the
# stack can name, so its reference resolves as an ordinary one, in 'e'
e <- new.env()
e$treePrior <- dbartsPriors$cgm(power = 3)
e$heldToken <- "logistic"
treePrior <- dbartsPriors$cgm(power = 5)
expect_equal(
  do.call(
    viaBart,
    list(tree.prior = quote(treePrior), keepTrees = TRUE),
    envir = e
  )$fit$model@tree.prior@power,
  3
)
viaEval <- function(...) {
  eval(quote(dbarts::dbarts(x, yBinary, control = control, ...)))
}
expect_equal(
  do.call(viaEval, list(family = quote(heldToken)), envir = e)$model@family,
  "logistic"
)
# a wrapper's own formal is an ordinary variable: it forwards a token or an
# object, not a constructor call
viaFormal <- function(fam) {
  dbarts::dbarts(x, y, control = control, family = fam)
}
expect_equal(residDf(viaFormal(dbartsFamilies$student(7))), 7)
expect_error(viaFormal(student(7)), "could not find function")

# --- the settings reach the specification ----------------------------------

# a fixed Student-t df written on the object is the df the model carries;
# the estimate mode is the bridge's 0
expect_equal(
  attr(
    dbarts::dbarts(x, y, control = control, family = student(3))$model,
    "resid.df"
  ),
  3
)
expect_equal(
  attr(
    dbarts::dbarts(x, y, control = control, family = student())$model,
    "resid.df"
  ),
  0
)
expect_null(
  attr(
    dbarts::dbarts(x, y, control = control, family = gaussian())$model,
    "resid.df"
  )
)
# and the same through the consumer surface, which shares the resolution
specStudent <- dbarts::dbartsSpec(
  dbarts::dbartsData(x, y),
  control = control,
  family = student(df = 9)
)
expect_equal(attr(specStudent$model, "resid.df"), 9)
expect_equal(specStudent$family, "gaussian")

# a count dispersion written on nbinom() is the dispersion the fit reports
set.seed(93L)
yCount <- rpois(length(y), 4)
fitFixedDispersion <- dbarts::bart(
  x,
  yCount,
  family = nbinom(dispersion = 4),
  n.trees = 10L,
  n.samples = 10L,
  n.burn = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_true(all(fitFixedDispersion$dispersion == 4))

# a hazard grid written on hazard() is the grid the expansion uses
set.seed(94L)
nSurv <- 60L
xSurv <- matrix(rnorm(nSurv * 2L), nSurv, 2L)
timeSurv <- pmax(1, ceiling(rexp(nSurv, 0.2)))
statusSurv <- rbinom(nSurv, 1L, 0.7)
fitHazard <- dbarts::bart(
  xSurv,
  cbind(timeSurv, statusSurv),
  family = hazard(breaks = c(0, 2, 5, max(timeSurv))),
  n.trees = 10L,
  n.samples = 10L,
  n.burn = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_equal(fitHazard$periods, c(2, 5, max(timeSurv)))
# a cap smaller than the expansion is named rather than run
expect_error(
  dbarts::bart(
    xSurv,
    cbind(timeSurv, statusSurv),
    family = hazard(max.rows = 10),
    n.trees = 10L,
    n.samples = 10L,
    n.burn = 5L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  ),
  "max.rows"
)

# --- the stored call is the caller's own ------------------------------------

# the resolved object is a forwarding detail: a fit that named no family
# records none, and one that named a family records the expression written
callArgs <- list(
  n.trees = 5L,
  n.samples = 5L,
  n.burn = 2L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
callDefaulted <- do.call(dbarts::bart, c(list(x, y), callArgs))
expect_false("family" %in% names(callDefaulted$call))
callToken <- do.call(
  dbarts::bart,
  c(list(x, y), callArgs, list(family = "gaussian"))
)
expect_identical(callToken$call$family, "gaussian")
callObject <- dbarts::bart(
  x,
  y,
  family = student(3),
  n.trees = 5L,
  n.samples = 5L,
  n.burn = 2L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_identical(callObject$call$family, quote(student(3)))

# --- the residual prior rides the family ------------------------------------

# the setting completes dec-B98's rule: the residual scale's own prior is a
# gaussian-family setting, so it is written inside the family call and
# nowhere else. The prior vocabulary resolves there, as it does inside
# 'tree.prior' and 'node.prior'.
expect_inherits(
  dbartsFamilies$gaussian(sigma = dbartsPriors$chisq(5, 0.75))@settings$sigma,
  "dbartsChiSqPrior"
)
expect_null(dbartsFamilies$gaussian()@settings$sigma)
expect_inherits(
  dbartsFamilies$student(3, sigma = dbartsPriors$fixed(2))@settings$sigma,
  "dbartsFixedPrior"
)
expect_equal(dbartsFamilies$student(3)@settings$df, 3)
expect_null(dbartsFamilies$aft()@settings$sigma)
# a bare constructor name means its defaults, as everywhere else
expect_equal(
  dbartsFamilies$gaussian(sigma = dbartsPriors$chisq)@settings$sigma@df,
  3
)
expect_error(
  dbartsFamilies$gaussian(sigma = 3),
  pattern = "must be a residual prior"
)
expect_error(
  dbartsFamilies$aft(sigma = dbartsPriors$normal()),
  pattern = "must be a residual prior"
)
# it reads back as the call that wrote it, not as an S4 object dump
expect_stdout(
  show(dbartsFamilies$gaussian(sigma = dbartsPriors$chisq(5, 0.75))),
  pattern = "gaussian(sigma = chisq(5, 0.75))",
  fixed = TRUE
)

residArgs <- list(
  n.trees = 5L,
  n.samples = 7L,
  n.burn = 3L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 217L,
  keepSampler = TRUE,
  verbose = FALSE
)
residFit <- function(...) do.call(dbarts::bart, c(list(x, y), residArgs, ...))

# the shipped default, written out: the same draws as writing nothing
fitDefault <- residFit()
fitDefaultNamed <- dbarts::bart(
  x,
  y,
  family = gaussian(sigma = chisq(3, 0.9)),
  n.trees = 5L,
  n.samples = 7L,
  n.burn = 3L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 217L,
  keepSampler = TRUE,
  verbose = FALSE
)
expect_identical(fitDefault$yhat.train, fitDefaultNamed$yhat.train)
expect_identical(fitDefault$sigma, fitDefaultNamed$sigma)

# a fixed residual scale suppresses the draw, reached through the family
fitFixed <- dbarts::bart(
  x,
  y,
  family = gaussian(sigma = fixed(1)),
  n.trees = 5L,
  n.samples = 7L,
  n.burn = 3L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 217L,
  keepSampler = TRUE,
  verbose = FALSE
)
expect_inherits(fitFixed$fit$model@resid.prior, "dbartsFixedPrior")
expect_true(all(abs(fitFixed$sigma - 1) < 1e-8))

# the sampler constructor reaches the prior the same way: it is the raw prior
# triple's third member, read off the family object the caller wrote it on
residControl <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 5L,
  n.samples = 5L,
  seed = 217L,
  updateState = FALSE
)
samplerViaFamily <- dbarts::dbarts(
  x,
  y,
  family = gaussian(sigma = fixed(2)),
  control = residControl
)
expect_inherits(samplerViaFamily$model@resid.prior, "dbartsFixedPrior")
expect_equal(samplerViaFamily$model@resid.prior@value, 2)
expect_equal(
  dbarts::dbartsSpec(
    dbarts::dbartsData(x, y),
    family = gaussian(sigma = fixed(2))
  )$model@resid.prior@value,
  2
)

# 'sigest' is the estimate a chisq prior calibrates against, so it stands
# beside one; a fixed residual scale has nothing to calibrate and refuses it
expect_silent(dbarts::dbarts(
  x,
  y,
  family = gaussian(sigma = chisq(5, 0.9)),
  sigest = 1.5,
  control = residControl
))
expect_error(
  dbarts::dbarts(
    x,
    y,
    family = gaussian(sigma = fixed(2)),
    sigest = 1.5,
    control = residControl
  ),
  pattern = "no effect under a fixed residual scale"
)
expect_silent(dbarts::bart(
  x,
  y,
  family = gaussian(sigma = chisq(5, 0.9)),
  sigest = 1.5,
  n.trees = 5L,
  n.samples = 7L,
  n.burn = 3L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 217L,
  verbose = FALSE
))
expect_error(
  dbarts::bart(
    x,
    y,
    family = gaussian(sigma = fixed(2)),
    sigest = 1.5,
    n.trees = 5L,
    n.samples = 7L,
    n.burn = 3L,
    n.chains = 1L,
    n.threads = 1L,
    seed = 217L,
    verbose = FALSE
  ),
  pattern = "no effect under a fixed residual scale"
)

# and a binary family, which has no residual scale to give a prior to, has no
# 'sigma' argument to write one in
expect_error(
  dbartsFamilies$probit(sigma = dbartsPriors$fixed(2)),
  pattern = "unused argument"
)
