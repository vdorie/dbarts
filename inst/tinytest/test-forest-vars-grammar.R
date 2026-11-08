# forest()'s predictors, its basis and varianceForest()'s predictors are read
# alike: bare terms, a one-sided formula in place or held in a variable (the
# tilde accepted, never required), and a character vector. Every spelling of
# the same predictors gives the identical model and identical draws.

forest <- dbartsForests$forest
varianceForest <- dbartsForests$varianceForest

set.seed(5)
n <- 80L
frame <- data.frame(
  y = rnorm(n),
  x1 = rnorm(n),
  x2 = rnorm(n),
  x3 = rnorm(n),
  w = rnorm(n)
)
control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.samples = 8L,
  n.burn = 4L,
  verbose = FALSE
)
draws <- function(expr) {
  set.seed(11L)
  sampler <- eval(expr)
  sampler$run()$train
}
refusal <- function(expr) {
  tryCatch(
    {
      force(expr)
      "accepted"
    },
    error = conditionMessage
  )
}
same <- function(cases) {
  reference <- draws(cases[[1L]])
  for (case in cases[-1L]) {
    expect_identical(draws(case), reference, info = deparse(case))
  }
}

twoTerms <- ~ x1 + x2
names12 <- c("x1", "x2")
heldNames <- c("x2", "x3")
heldTerms <- ~ x2 + x3

## --- forest()'s predictors, in a forests list ---------------------------------
inList <- function(vars) {
  bquote(dbarts(
    y ~ x1 + x2 + x3,
    frame,
    forests = list(forest(.(vars)), forest(x3, basis = w)),
    control = control,
    n.trees = 6L
  ))
}
same(list(
  inList(quote(x1 + x2)),
  inList(quote(~ x1 + x2)),
  inList(quote(twoTerms)),
  inList(quote(names12)),
  inList(c("x1", "x2")),
  inList(quote(. - x3)),
  inList(quote(~ . - x3)),
  inList(1:2)
))

## --- forest()'s predictors, in a formula --------------------------------------
inFormula <- function(vars) {
  bquote(dbarts(
    y ~ forest(.(vars)) + forest(x3, basis = w),
    frame,
    control = control,
    n.trees = 6L
  ))
}
same(list(
  inFormula(quote(x1 + x2)),
  inFormula(quote(~ x1 + x2)),
  inFormula(quote(twoTerms)),
  inFormula(quote(names12)),
  inFormula(c("x1", "x2"))
))

## --- the basis ------------------------------------------------------------------
heldBasis <- ~w
same(list(
  quote(dbarts(
    y ~ forest(x1) + forest(x2, basis = w),
    frame,
    control = control,
    n.trees = 6L
  )),
  quote(dbarts(
    y ~ forest(x1) + forest(x2, basis = ~w),
    frame,
    control = control,
    n.trees = 6L
  )),
  quote(dbarts(
    y ~ forest(x1) + forest(x2, basis = heldBasis),
    frame,
    control = control,
    n.trees = 6L
  ))
))

## --- varianceForest()'s predictors ----------------------------------------------
withVariance <- function(variance) {
  bquote(dbarts(
    y ~ x1 + x2 + x3,
    frame,
    variance = .(variance),
    control = control,
    n.trees = 6L
  ))
}
same(list(
  withVariance(quote(varianceForest(x2 + x3))),
  withVariance(quote(varianceForest(~ x2 + x3))),
  withVariance(quote(varianceForest(heldTerms))),
  withVariance(quote(varianceForest(heldNames))),
  withVariance(quote(varianceForest(c("x2", "x3")))),
  withVariance(quote(varianceForest(2:3))),
  withVariance(quote(varianceForest(. - x1))),
  withVariance(~ x2 + x3)
))
same(list(
  withVariance(quote(varianceForest())),
  withVariance(quote(varianceForest(.))),
  withVariance(quote(varianceForest(x1 + x2 + x3))),
  withVariance(TRUE)
))
everyColumn <- TRUE
expect_identical(
  draws(withVariance(quote(varianceForest(everyColumn)))),
  draws(withVariance(TRUE))
)

## --- a plain variance formula reads as varianceForest() reads it ---------------
# the dot and a transformed term are taken, as varianceForest() takes them
same(list(
  withVariance(quote(varianceForest(. - x1))),
  withVariance(~ . - x1)
))
logFrame <- frame
logFrame$x3 <- abs(logFrame$x3) + 1
logVariance <- function(variance) {
  bquote(dbarts(
    y ~ x1 + x2 + log(x3),
    logFrame,
    variance = .(variance),
    control = control,
    n.trees = 6L
  ))
}
same(list(
  logVariance(quote(varianceForest(~ log(x3)))),
  logVariance(~ log(x3)),
  logVariance("log(x3)")
))
# a star, a colon and an I() term are refused with the message naming the
# fit's predictors, where the variables mentioned were once read silently
for (term in c("x2 * x3", "x2:x3", "I(x2^2)")) {
  expect_identical(
    refusal(eval(withVariance(str2lang(paste0("~ ", term))))),
    paste0(
      "'",
      term,
      "' is not a predictor of this fit (x1, x2, x3); here the 'variance' ",
      "formula selects among them"
    ),
    info = term
  )
  expect_true(
    grepl(
      "is not a predictor",
      refusal(eval(withVariance(
        str2lang(paste0("varianceForest(~ ", term, ")"))
      ))),
      fixed = TRUE
    ),
    info = term
  )
}
expect_identical(
  refusal(eval(withVariance(~nosuch))),
  paste0(
    "'nosuch' is not a predictor of this fit (x1, x2, x3); here the ",
    "'variance' formula selects among them"
  )
)

# ~ 1 names no predictor: refused with the remedy; an intercept term beside a
# predictor is refused as varianceForest() refuses it
expect_true(grepl(
  "omit variance for a constant variance, or write ~ . for every column",
  refusal(eval(withVariance(~1))),
  fixed = TRUE
))
for (written in c("~ x2 - 1", "~ 0 + x2")) {
  expect_true(
    grepl(
      "an intercept term",
      refusal(eval(withVariance(str2lang(written)))),
      fixed = TRUE
    ),
    info = written
  )
  expect_true(
    grepl(
      "an intercept term",
      refusal(eval(withVariance(str2lang(paste0(
        "varianceForest(",
        written,
        ")"
      ))))),
      fixed = TRUE
    ),
    info = written
  )
}

## --- what is refused --------------------------------------------------------------
twoSided <- y ~ x1 + x2
expect_identical(
  refusal(forest(y ~ x1)),
  paste0(
    "forest()'s first argument is a formula with a left-hand side; write the ",
    "predictors as the right-hand side alone, as x1 + x2 or ~ x1 + x2"
  )
)
expect_identical(
  refusal(varianceForest(y ~ x1)),
  paste0(
    "varianceForest()'s first argument is a formula with a left-hand side; ",
    "write the predictors as the right-hand side alone, as x1 + x2 or ",
    "~ x1 + x2"
  )
)
expect_identical(
  refusal(eval(inList(quote(twoSided)))),
  paste0(
    "forest()'s first argument, 'twoSided', holds a formula with a left-hand ",
    "side; write the predictors as the right-hand side alone, as ~ x1 + x2"
  )
)
expect_identical(
  refusal(eval(inFormula(quote(twoSided)))),
  paste0(
    "forest()'s first argument, 'twoSided', holds a formula with a left-hand ",
    "side; write the predictors as the right-hand side alone, as ~ x1 + x2"
  )
)
expect_identical(
  refusal(eval(withVariance(quote(varianceForest(nosuch + x2))))),
  paste0(
    "'nosuch' is not a predictor of this fit (x1, x2, x3); here ",
    "varianceForest()'s first argument selects among them"
  )
)
expect_identical(
  refusal(eval(withVariance(quote(varianceForest(twoSided))))),
  paste0(
    "varianceForest()'s first argument, 'twoSided', holds a formula with a ",
    "left-hand side; write the predictors as the right-hand side alone, as ",
    "~ x1 + x2"
  )
)

## --- the printed specification ------------------------------------------------------
expect_identical(
  format(varianceForest(~ x2 + x3))[[1L]],
  "vars           = x2 + x3"
)
