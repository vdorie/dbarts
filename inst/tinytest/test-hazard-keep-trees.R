# bart()'s keepTrees = NULL default: a hazard family keeps its trees (its
# survival curves replay them), every other family does not. An explicit
# TRUE or FALSE is taken as given, and keepSampler follows the resolved value.

if (!requireNamespace("survival", quietly = TRUE)) {
  exit_file("survival is not installed")
}
Surv <- survival::Surv

fitDefault <- function(...) {
  dbarts::bart(
    ...,
    n.trees = 5L,
    n.samples = 6L,
    n.burn = 4L,
    n.chains = 2L,
    n.threads = 1L,
    seed = 3L,
    verbose = FALSE
  )
}
keeps <- function(fit) !is.null(fit$fit) && fit$fit$control@keepTrees

set.seed(21)
n <- 60L
df <- data.frame(a = runif(n), b = rnorm(n))
df$t <- rexp(n, exp(df$a - 0.5))
df$s <- rbinom(n, 1L, 0.8)
x <- as.matrix(df[, c("a", "b")])

hazardMatrix <- fitDefault(x, cbind(df$t, df$s), family = "hazard")
hazardFormula <- fitDefault(Surv(t, s) ~ a + b, df, family = "hazard")
hazardLogistic <- fitDefault(
  Surv(t, s) ~ a + b,
  df,
  family = hazard(link = "logistic")
)
expect_true(keeps(hazardMatrix))
expect_true(keeps(hazardFormula))
expect_true(keeps(hazardLogistic))
expect_false(is.null(hazardFormula$fit))

# a default hazard fit gives its curves and its partial dependence with no
# refit warning
for (fit in list(hazardMatrix, hazardFormula)) {
  warned <- character()
  withCallingHandlers(
    {
      curves <- dbarts::survivalProbabilities(fit)
      pd <- dbarts::pdbart(fit)
    },
    warning = function(w) {
      warned <<- c(warned, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(warned, character())
  expect_true(is.array(curves))
  expect_true(inherits(pd, "pdbart"))
}

# every other family keeps no trees by default
expect_false(keeps(fitDefault(Surv(t, s) ~ a + b, df, family = "aft")))
expect_false(keeps(fitDefault(Surv(t, s) ~ a + b, df)))
expect_false(keeps(fitDefault(a ~ b, df)))
expect_false(keeps(fitDefault(a ~ b, df, family = "gaussian")))
expect_false(keeps(fitDefault(
  (a > 0.5) ~ b,
  df,
  family = "probit"
)))

# explicit values are taken as given; NULL is the default
expect_false(keeps(fitDefault(
  x,
  cbind(df$t, df$s),
  family = "hazard",
  keepTrees = FALSE
)))
expect_true(keeps(fitDefault(a ~ b, df, keepTrees = TRUE)))
expect_true(keeps(fitDefault(
  x,
  cbind(df$t, df$s),
  family = "hazard",
  keepTrees = NULL
)))
expect_false(keeps(fitDefault(a ~ b, df, keepTrees = NULL)))

# a control that keeps trees stands for a hazard fit and for any other
expect_true(keeps(fitDefault(
  a ~ b,
  df,
  control = dbartsControl(keepTrees = TRUE, n.chains = 2L, n.threads = 1L)
)))

# keepSampler follows the resolved value, and an explicit one is kept
expect_true(is.null(fitDefault(a ~ b, df)$fit))
expect_false(is.null(fitDefault(a ~ b, df, keepSampler = TRUE)$fit))

# an explicit FALSE on a hazard fit refuses the curves, naming the fix
hazardNoTrees <- fitDefault(
  Surv(t, s) ~ a + b,
  df,
  family = "hazard",
  keepTrees = FALSE
)
expect_error(dbarts::survivalProbabilities(hazardNoTrees), "keepTrees = TRUE")
