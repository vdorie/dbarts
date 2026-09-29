# Base R's family objects as 'family' (dec-B134): accepted spellings resolve
# to the dbarts family with the same likelihood and link; the rest are
# refused by name.

source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

x <- testData$x
yb <- as.integer(testData$y > median(testData$y))
y <- testData$y

resolve <- function(
  expr,
  tokens = c("auto", "gaussian", "probit", "logistic")
) {
  dbarts:::resolveFamily(expr, tokens, "test", environment())
}

# accepted spellings resolve to the ordinary dbarts family
expect_identical(resolve(quote(stats::gaussian())), resolve(quote(gaussian())))
expect_identical(
  resolve(quote(stats::gaussian(link = "identity"))),
  resolve(quote("gaussian"))
)
expect_identical(
  resolve(quote(binomial(link = "probit"))),
  resolve(quote("probit"))
)
expect_identical(
  resolve(quote(binomial(link = "logit"))),
  resolve(quote("logistic"))
)
# the word and the function mean binomial(), the logit link
expect_identical(resolve(quote("binomial")), resolve(quote("logistic")))
expect_identical(resolve(quote(binomial)), resolve(quote("logistic")))
expect_identical(resolve(quote(binomial())), resolve(quote("logistic")))
expect_identical(resolve(quote(stats::binomial)), resolve(quote("logistic")))
# a partial word is still match.arg's, and gaussian keeps its meaning
expect_identical(resolve(quote("gauss")), resolve(quote("gaussian")))

# refusals name the family and link
expect_error(
  resolve(quote(binomial(link = "cloglog"))),
  "binomial(link = \"cloglog\") is not supported",
  fixed = TRUE
)
expect_error(
  resolve(quote(binomial(link = "cauchit"))),
  "cauchit",
  fixed = TRUE
)
expect_error(
  resolve(quote(stats::gaussian(link = "log"))),
  "gaussian(link = \"log\") is not supported",
  fixed = TRUE
)
for (expr in list(
  quote(poisson),
  quote(poisson()),
  quote(Gamma()),
  quote(quasibinomial()),
  quote(quasipoisson()),
  quote(inverse.gaussian())
)) {
  fam <- eval(if (is.symbol(expr)) expr else expr)
  fam <- if (is.function(fam)) fam() else fam
  expect_error(resolve(expr), fam$family, fixed = TRUE)
  expect_error(resolve(expr), "dbartsFamilies", fixed = TRUE)
}
foreign <- structure(
  list(family = "Negative Binomial(3)", link = "log"),
  class = "family"
)
expect_error(resolve(quote(foreign)), "Negative Binomial(3)", fixed = TRUE)
expect_error(resolve(quote(1L)), "family name or a family object")

# the entry point's own list still applies after mapping
expect_error(
  resolve(quote(binomial), tokens = c("auto", "gaussian", "probit")),
  "does not fit family \"logistic\""
)

# through the entry points: the mapped family draws exactly as its dbarts
# spelling, and family() reports the ordinary dbarts family
fit <- function(family, ...) {
  bart(
    x,
    yb,
    family = family,
    verbose = FALSE,
    n.trees = 5L,
    n.chains = 1L,
    n.threads = 1L,
    n.samples = 20L,
    n.burn = 10L,
    seed = 3L
  )
}
a <- fit(binomial(link = "probit"))
b <- fit("probit")
expect_identical(a$yhat.train, b$yhat.train)
expect_identical(family(a), family(b))
expect_identical(fit(binomial)$yhat.train, fit("logistic")$yhat.train)
expect_identical(fit("binomial")$yhat.train, fit("logistic")$yhat.train)
expect_identical(
  bart(
    x,
    y,
    family = stats::gaussian(),
    verbose = FALSE,
    n.trees = 5L,
    n.chains = 1L,
    n.threads = 1L,
    n.samples = 20L,
    n.burn = 10L,
    seed = 3L
  )$yhat.train,
  bart(
    x,
    y,
    family = "gaussian",
    verbose = FALSE,
    n.trees = 5L,
    n.chains = 1L,
    n.threads = 1L,
    n.samples = 20L,
    n.burn = 10L,
    seed = 3L
  )$yhat.train
)
expect_error(fit(poisson), "poisson")
expect_error(
  dbarts(yb ~ x, data = data.frame(yb = yb, x = x[, 1L]), family = Gamma()),
  "Gamma"
)
expect_error(xbart(x, yb, family = binomial(link = "cloglog")), "cloglog")
expect_error(
  dbartsSpec(dbartsData(x, yb), family = quasibinomial()),
  "quasibinomial"
)

# bare gaussian takes link first, glm's order, and refuses other links
expect_identical(
  resolve(quote(gaussian(link = "identity"))),
  resolve(quote("gaussian"))
)
expect_identical(
  resolve(quote(gaussian("identity"))),
  resolve(quote("gaussian"))
)
expect_error(
  resolve(quote(gaussian(link = "log"))),
  "gaussian(link = \"log\") is not supported",
  fixed = TRUE
)
expect_error(resolve(quote(gaussian(sigma = "identity"))), "link =")
expect_error(resolve(quote(gaussian(fixed(1)))), "sigma =")
expect_identical(
  resolve(quote(gaussian(sigma = fixed(1))))@settings,
  dbartsFamilies$gaussian(sigma = dbartsPriors$fixed(1))@settings
)

# stats family names as strings map or are refused by name
expect_identical(resolve(quote("binomial")), resolve(quote("logistic")))
for (nm in c(
  "poisson",
  "Gamma",
  "quasibinomial",
  "quasipoisson",
  "inverse.gaussian",
  "quasi"
)) {
  expect_error(
    resolve(bquote(.(nm))),
    paste0("family \"", nm, "\" is not supported"),
    fixed = TRUE
  )
}

# a continuous response under binomial's logit link says where it came from
expect_error(
  bart(
    x,
    y,
    family = binomial,
    verbose = FALSE,
    n.trees = 5L,
    n.chains = 1L,
    n.threads = 1L,
    n.samples = 20L,
    n.burn = 10L
  ),
  "logit link"
)
