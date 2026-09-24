# What a fit is: $family (the family as specified), family() (the same with
# its settings), and the stored descriptors a family or model fixes -
# resid.scale on a family with a residual law, n.forests on every bart fit -
# none of them logical and none inferred from which draw channels a run
# kept.

set.seed(71L)
n <- 60L
x <- matrix(rnorm(2L * n), n, 2L)
y <- x[, 1L] + rnorm(n)
yBinary <- as.integer(y > 0)
time <- sample(1:4, n, TRUE)
status <- rbinom(n, 1L, 0.7)

fitOf <- function(...) {
  suppressWarnings(suppressMessages(bart(
    x,
    ...,
    n.trees = 5L,
    n.samples = 4L,
    n.burn = 2L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  )))
}

fits <- list(
  gaussian = fitOf(y),
  probit = fitOf(yBinary),
  logistic = fitOf(yBinary, family = "logistic"),
  student = fitOf(y, family = dbartsFamilies$student(3)),
  heteroscedastic = fitOf(y, variance = varianceForest(), keepFits = FALSE),
  aft = fitOf(survival::Surv(exp(y), status)),
  hazard = fitOf(survival::Surv(time, status), family = "hazard"),
  multinomial = fitOf(factor(sample(3L, n, TRUE))),
  ordinal = fitOf(factor(sample(3L, n, TRUE), ordered = TRUE)),
  nbinom = fitOf(rpois(n, 3), family = "nbinom"),
  hurdle = fitOf(pmax(y, 0), family = "hurdle.lognormal")
)

# $family is the family as specified, "auto" resolved; family() the same
# family with its settings; the engine family the link follows is looked up
specifiedTokens <- c(
  gaussian = "gaussian",
  probit = "probit",
  logistic = "logistic",
  student = "student",
  heteroscedastic = "gaussian",
  aft = "aft",
  hazard = "hazard.probit",
  multinomial = "multinomial",
  ordinal = "ordinal",
  nbinom = "nbinom",
  hurdle = "hurdle.lognormal"
)
for (name in names(fits)) {
  fit <- fits[[name]]
  expect_identical(fit$family, specifiedTokens[[name]], info = name)
  expect_true(is(family(fit), "dbartsFamily"), info = name)
  expect_identical(family(fit)@token, fit$family, info = name)
  # no logical component on any fit
  expect_false(any(vapply(fit, is.logical, logical(1L))), info = name)
}
expect_equal(family(fits$student)@settings$df, 3)
expect_identical(dbarts:::fitEngineFamily(fits$student), "gaussian")
expect_identical(dbarts:::fitEngineFamily(fits$hazard), "probit")
# a family the lookup does not name stops rather than taking a default link
unknownFamily <- fits$probit
unknownFamily$family <- "cloglog"
expect_error(fitted(unknownFamily), "unknown family 'cloglog'")

# the residual law's scale model exists exactly on the families with one;
# its shape is not stored, since the family fixes it
for (name in c("gaussian", "student", "heteroscedastic", "aft")) {
  expect_true(is.character(fits[[name]]$resid.scale), info = name)
}
for (name in c("probit", "logistic", "hazard")) {
  expect_null(fits[[name]]$resid.scale, info = name)
}
for (name in names(fits)) {
  expect_null(fits[[name]]$resid.dist, info = name)
}
expect_identical(fits$gaussian$resid.scale, "constant")
# the scale model survives keepFits = FALSE, which dropped its draws
expect_identical(fits$heteroscedastic$resid.scale, "forest")
expect_null(fits$heteroscedastic$s.train)

# n.forests on every bart fit, 1 for a single forest
for (name in names(fits)[vapply(fits, inherits, logical(1L), "bart")]) {
  expect_identical(fits[[name]]$n.forests, 1L, info = name)
}

# the hurdle's components carry their own specified families
expect_identical(family(fits$hurdle$occupancy)@token, "probit")
expect_identical(family(fits$hurdle$positive)@token, "gaussian")

# a hazard fit is recognized by its specified family: survivalProbabilities
# takes the hazard branch, which then asks for the trees
expect_error(survivalProbabilities(fits$hazard), "keepTrees")
