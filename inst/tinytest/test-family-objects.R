# The family objects (dec-B98, dec-B101): one 'family' argument carrying
# every setting only one family reads. The constructors, the exported
# vocabulary, the token/object equivalence, and that a setting written on
# the object reaches the sampler specification unchanged.

source(
  system.file("common", "friedmanData.R", package = "dbarts"),
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
