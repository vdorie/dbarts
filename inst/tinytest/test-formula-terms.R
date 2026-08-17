# forest() written inside a formula - bare, or as one operand of a ':'/'*'
# node - declares an additional amplitude-coupled forest (docs/design/bcf.md),
# reached by rewriting the formula and feeding the existing
# forests = list(forest(), forest(basis = )) channel (R/formulaTerms.R).
# Block A: a formula with no forest() call anywhere is untouched. Block B: the
# ':' sugar and the symbolic vars slot desugar to exactly what the general
# forest(basis =, vars =) form would, at the same seed. Block C: the refusal
# matrix, one cell per shape the grammar or the family/data does not support.

n <- 60L
set.seed(0)
x1 <- runif(n)
x2 <- runif(n)
a <- runif(n)
b <- runif(n)
z <- rbinom(n, 1L, 0.5)
zf <- factor(sample(c("u", "v", "w"), n, replace = TRUE))
zc <- sample(c("p", "q"), n, replace = TRUE)
o <- rep(0.1, n)
w <- rep(1, n)
weirdName <- runif(n)
y <- x1 + z * (1 + x2) + rnorm(n, 0, 0.2)
d <- data.frame(
  y = y,
  x1 = x1,
  x2 = x2,
  a = a,
  b = b,
  z = z,
  zf = zf,
  zc = zc,
  o = o,
  w = w
)
d[["weird name"]] <- weirdName

tinyArgs <- list(
  seed = 1L,
  n.samples = 3L,
  n.burn = 3L,
  n.trees = 3L,
  n.chains = 1L,
  n.threads = 1L,
  keepTrees = FALSE,
  keepSampler = TRUE,
  verbose = FALSE
)

fit <- function(formula, ...) {
  do.call(
    dbarts::bart2,
    c(list(formula = formula, data = d), tinyArgs, list(...))
  )
}
refuses <- function(formula, pattern, ...) {
  # nolint next: object_usage_linter. tinytest attaches expect_* at run time.
  expect_error(fit(formula, ...), pattern, info = deparse(formula))
}

## --- Block A: a formula with no forest() call is untouched -----------------
## (this file's own draw-level check; the cross-build byte-identity gate is
## an external A/B script run at both the pre- and post-change revisions,
## since the change touches no baseline the harness already covers)
expect_true(is.null(dbarts:::walkFormulaTerms(y ~ a + b)))
expect_true(is.null(dbarts:::walkFormulaTerms(y ~ .)))
expect_true(is.null(dbarts:::walkFormulaTerms(y ~ . - z)))
# a ':' whose neither operand is a forest() call is not a term at all
expect_true(is.null(dbarts:::walkFormulaTerms(y ~ a + z:b)))
expect_true(is.null(dbarts:::walkFormulaTerms(y ~ a - 1)))
expect_true(is.null(dbarts:::walkFormulaTerms(y ~ a + offset(o))))
expect_true(is.null(dbarts:::walkFormulaTerms(as.formula("y ~ `weird name`"))))

expect_silent(fit(y ~ a + b))
expect_silent(fit(y ~ .))
expect_silent(fit(y ~ . - z))
expect_silent(fit(y ~ a - 1))
expect_silent(fit(y ~ a + offset(o)))
expect_silent(fit(as.formula("y ~ `weird name`")))
expect_silent(fit(y ~ a, subset = 1:40))
expect_silent(fit(y ~ a, weights = w, offset = o))
expect_silent(dbarts::bart2(
  cbind(a, b),
  y,
  seed = 1L,
  n.samples = 3L,
  n.burn = 3L,
  n.trees = 3L,
  n.chains = 1L,
  n.threads = 1L,
  keepTrees = FALSE,
  verbose = FALSE
))

## --- Block B: sugar desugars to exactly the general named form -------------
bcfAttr <- function(result) attr(result$fit$control, "bartcore.bcf")
expectSameForest <- function(sugarFormula, namedFormula, ...) {
  s <- fit(sugarFormula, ...)
  h <- fit(namedFormula, ...)
  # nolint next: object_usage_linter. tinytest attaches expect_* at run time.
  expect_identical(s$fit$data@bases, h$fit$data@bases)
  # nolint next: object_usage_linter. tinytest attaches expect_* at run time.
  expect_identical(bcfAttr(s), bcfAttr(h))
  # nolint next: object_usage_linter. tinytest attaches expect_* at run time.
  expect_identical(s$yhat.train, h$yhat.train)
}

# (11) numeric z
expectSameForest(
  y ~ x1 + x2 + z:forest(x1 + x2),
  y ~ x1 + x2 + forest(x1 + x2, basis = ~z)
)
# (12) factor zf, no reference level dropped: n x nlevels
zfSugar <- fit(y ~ x1 + x2 + zf:forest(x1 + x2))
expect_identical(dim(zfSugar$fit$data@bases[[2L]]), c(n, 3L))
expectSameForest(
  y ~ x1 + x2 + zf:forest(x1 + x2),
  y ~ x1 + x2 + forest(x1 + x2, basis = ~zf)
)
# (13) factor(z), z numeric 0/1: n x 2
factorSugar <- fit(y ~ x1 + x2 + factor(z):forest(x1 + x2))
expect_identical(dim(factorSugar$fit$data@bases[[2L]]), c(n, 2L))
expectSameForest(
  y ~ x1 + x2 + factor(z):forest(x1 + x2),
  y ~ x1 + x2 + forest(x1 + x2, basis = ~ factor(z))
)
# (14) compound (a + b): ONE forest, n x 2, not the elementwise sum
compoundSugar <- fit(y ~ x1 + x2 + (a + b):forest(x1 + x2))
expect_identical(length(compoundSugar$fit$data@bases), 2L)
expect_identical(dim(compoundSugar$fit$data@bases[[2L]]), c(n, 2L))
expect_false(
  isTRUE(all.equal(
    compoundSugar$fit$data@bases[[2L]][, 1L] +
      compoundSugar$fit$data@bases[[2L]][, 2L],
    a + b
  )) &&
    ncol(compoundSugar$fit$data@bases[[2L]]) == 1L
)
expectSameForest(
  y ~ x1 + x2 + (a + b):forest(x1 + x2),
  y ~ x1 + x2 + forest(x1 + x2, basis = ~ cbind(a, b))
)
# (15) operand order is immaterial
expectSameForest(
  y ~ x1 + x2 + forest(x1 + x2):z,
  y ~ x1 + x2 + z:forest(x1 + x2)
)
# (16) a dbarts::-qualified head reaches the same configuration
expectSameForest(
  y ~ x1 + x2 + z:dbarts::forest(x1 + x2),
  y ~ x1 + x2 + z:forest(x1 + x2)
)
# (17) the symbolic slot vs vars = by name
expectSameForest(
  y ~ x1 + x2 + forest(x1 + x2, basis = ~z),
  y ~ x1 + x2 + forest(vars = c("x1", "x2"), basis = ~z)
)
# (18) a symbolic factor name expands to all its indicator columns
indicatorFit <- fit(
  y ~ x1 + x2 + zf + forest(zf, basis = ~x1),
  factors = "indicators"
)
expect_identical(
  colnames(indicatorFit$fit$data@x)[3:5],
  c("zf.u", "zf.v", "zf.w")
)
expect_identical(bcfAttr(indicatorFit)$vars[[2L]], 3:5)
expectSameForest(
  y ~ x1 + x2 + zf + forest(zf, basis = ~x1),
  y ~ x1 + x2 + zf + forest(vars = c("zf.u", "zf.v", "zf.w"), basis = ~x1),
  factors = "indicators"
)

## --- Block C: the refusal matrix --------------------------------------------
# (19) '*' with a forest operand
refuses(y ~ x1 + x2 + z * forest(x1 + x2), "z \\* forest.*write z:forest")
# (20) the ancestor chain, not just the immediate parent
refuses(y ~ x1 + x2 + I(z:forest(x1 + x2)), "top-level additive term")
refuses(y ~ (x1 + z:forest(x2))^2, "top-level additive term")
# (21) a term on the left-hand side
refuses(forest(x1) ~ x2, "left-hand side")
# (21a) a forest() inside a removal
refuses(y ~ x1 - z:forest(x2), "top-level additive term")
# (22)-(24) the symbolic-vars grammar
refuses(y ~ x1 + x2 + forest(log(x1)), "not a supported forest\\(\\) vars")
refuses(y ~ x1 + x2 + forest(x1 - x2), "not a supported forest\\(\\) vars")
refuses(y ~ x1 + x2 + forest(.), "not a supported forest\\(\\) vars")
# (24a) a symbolic name absent from the rewritten right-hand side
refuses(y ~ x2 + z:forest(x1 + x2), "not a term of the rewritten")
# (24b) a compound operand with a non-numeric, non-logical member
refuses(y ~ x1 + x2 + (a + zf):forest(x1), "numeric or logical")
refuses(y ~ x1 + x2 + (a + zc):forest(x1), "numeric or logical")
# (24c) a multi-way colon chain
refuses(y ~ x1 + x2 + z:a:forest(x1), "more than one operand")
# a forest() nested inside a ':' operand, both associativity directions
refuses(y ~ x1 + forest(x1):z:a, "more than one operand")
refuses(y ~ x1 + a:(b:forest(x1)), "more than one operand")
# (24d) forest() on both sides of ':'
refuses(y ~ x1 + x2 + forest(x1):forest(x2), "both sides of ':'")
refuses(y ~ x1 + x2 + z:forest(x1):forest(x2), "both sides of ':'")
# (25) sugar and basis = together
refuses(y ~ x1 + x2 + z:forest(x1, basis = ~z), "both through ':' and")
# (26) the symbolic slot and vars = together
refuses(
  y ~ x1 + x2 + forest(x1 + x2, vars = c("x1", "x2"), basis = ~z),
  "positionally and by name"
)
# (27) a rewrite leaving no predictors, both no-intercept spellings
refuses(y ~ z:forest(x1), "no predictors left")
refuses(y ~ 0 + z:forest(x1), "no predictors left")
# (27a) an all-zero basis column
refuses(
  y ~ x1 + x2 + forest(x1 + x2, basis = ~ rep(0, n)),
  "all zeros"
)
# (28) 'test' with a term
refuses(y ~ x1 + x2 + z:forest(x1 + x2), "no test-basis channel", test = d)
# (29) restates the existing forests = collision with a pre-built dbartsData
dd <- dbarts::dbartsData(y ~ x1 + x2, d)
expect_error(
  dbarts:::dbarts(dd, forests = list(forest(), forest(basis = ~z))),
  "already a dbartsData"
)
# the forests = / formula-term collision: two spellings of the same
# multi-forest declaration on one call
expect_error(
  dbarts:::dbarts(
    y ~ x1 + x2 + z:forest(x1 + x2),
    d,
    forests = list(forest(), forest(basis = ~z)),
    control = dbarts::dbartsControl(
      n.chains = 1L,
      n.trees = 3L,
      n.burn = 3L,
      n.samples = 3L,
      n.threads = 1L,
      verbose = FALSE
    )
  ),
  "only be declared one way"
)
# (30) the family matrix - refused at family resolution, naming the family
termFormula <- y ~ x1 + x2 + z:forest(x1 + x2)
for (family in c("hazard", "hazard.probit", "hazard.logistic")) {
  refuses(termFormula, paste0("family \"", family, "\""), family = family)
}
refuses(termFormula, "family = \"multinomial\"", family = "multinomial")
refuses(
  termFormula,
  "family = \"hurdle.lognormal\"",
  family = "hurdle.lognormal"
)
refuses(termFormula, "family = \"hurdle.lognormal\"", family = "twopart")
for (family in c("aft", "ordinal", "nbinom")) {
  refuses(termFormula, paste0("family \"", family, "\""), family = family)
}
