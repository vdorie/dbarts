# The binomial response's R surface: cbind(successes, failures) under family =
# binomial and "logistic" on both interfaces and xbart; the two scales a fit
# reports on, as glm's; the log-likelihood with its binomial coefficient; the
# predictive draws' trials at training and new rows; the zero-trial row; the
# sampler's doors; and Refusals 11 to 18 by their words. The bridge's reading
# of the counts slot is test-binomial-bridge.R's.

set.seed(20261107L)
n <- 40L
bnm <- data.frame(x1 = runif(n), x2 = runif(n))
trials <- rep_len(c(0L, 1L, 3L, 5L, 2L), n)
bnm$k <- rbinom(n, trials, plogis(2 * bnm$x1 - 1))
bnm$f <- trials - bnm$k
X <- as.matrix(bnm[c("x1", "x2")])
counts <- cbind(k = bnm$k, f = bnm$f)
fitArgs <- list(
  n.samples = 20L,
  n.burn = 10L,
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 10L,
  verbose = FALSE,
  seed = 3L
)
quietly <- function(expr) suppressWarnings(suppressMessages(expr))
fitBinomial <- function(...) {
  quietly(do.call(dbarts::bart, c(list(...), fitArgs)))
}

# --- the routes: formula and matrix interfaces, family spellings ---

fit <- fitBinomial(
  cbind(k, f) ~ x1 + x2,
  bnm,
  family = binomial,
  keepTrees = TRUE
)
expect_inherits(fit, "bart")
expect_identical(fit$family, "logistic")
# fit$y keeps the two columns as written
expect_identical(unname(fit$y), unname(cbind(bnm$k, bnm$f)))
expect_identical(colnames(fit$y), c("k", "f"))
for (family in list(stats::binomial(), "logistic")) {
  other <- fitBinomial(
    cbind(k, f) ~ x1 + x2,
    bnm,
    family = family,
    keepTrees = TRUE
  )
  expect_identical(other$yhat.train, fit$yhat.train)
}
onMatrix <- fitBinomial(X, counts, family = binomial, keepTrees = TRUE)
expect_identical(unname(onMatrix$yhat.train), unname(fit$yhat.train))
onData <- quietly(dbarts::dbarts(
  dbarts::dbartsData(X, counts = counts),
  family = "logistic"
))
expect_true(dbarts:::samplerHoldsTrials(onData))
rm(family, other, onMatrix, onData)

# the zero-trial warning names the binomial's rows, once a session
onceState <- dbarts:::onceWarnState
onceState[["multinomialZeroTrials"]] <- NULL
expect_warning(
  suppressMessages(do.call(
    dbarts::bart,
    c(list(cbind(k, f) ~ x1 + x2, bnm, family = binomial), fitArgs)
  )),
  paste0(
    "binomial rows with zero trials (8 of 40) contribute nothing to the ",
    "likelihood and still receive fitted probabilities"
  ),
  fixed = TRUE
)
rm(onceState)

# all trials 1 is the binary fit, draw for draw
y01 <- as.double(runif(n) < plogis(2 * bnm$x1 - 1))
binary <- fitBinomial(X, y01, family = binomial)
asPair <- fitBinomial(X, cbind(y01, 1 - y01), family = binomial)
expect_identical(unname(asPair$yhat.train), unname(binary$yhat.train))
rm(y01, binary, asPair)

# --- the two scales: probabilities, and counts for the draws ---

p <- dbarts::extract(fit, type = "ev")
expect_identical(dim(p), c(fitArgs$n.samples, n))
expect_true(all(p > 0 & p < 1))
expect_equal(p, plogis(dbarts::extract(fit, type = "bart")))
observed <- ifelse(trials > 0L, bnm$k / trials, 0)
expect_equal(
  unname(residuals(fit)),
  unname(observed - fitted(fit))
)
set.seed(5L)
residualsPpd <- residuals(fit, type = "ppd")
set.seed(5L)
fittedPpd <- fitted(fit, type = "ppd")
expect_equal(unname(residualsPpd), unname(bnm$k - fittedPpd))
# the predictive draws are Binomial(n, p) at a training row: whole counts no
# larger than its trials, 0 at a row of no trials
set.seed(6L)
ppd <- dbarts::extract(fit, type = "ppd")
expect_true(all(ppd == round(ppd)))
expect_true(all(sweep(ppd, 2L, trials, "<=")))
expect_true(all(ppd[, trials == 0L] == 0))
# and their mean estimates n p
set.seed(7L)
many <- replicate(200L, colMeans(dbarts::extract(fit, type = "ppd")))
expect_true(max(abs(rowMeans(many) - trials * colMeans(p))) < 0.15)
rm(p, observed, residualsPpd, fittedPpd, ppd, many)

# the zero-trial row: fitted, its residual -p, its draws 0
zeroRows <- which(trials == 0L)
expect_equal(
  unname(residuals(fit)[zeroRows]),
  unname(-fitted(fit)[zeroRows])
)

# --- the log-likelihood: dbinom's, binomial coefficient included ---

logLikelihood <- dbarts::extract(fit, type = "loglik")
p <- dbarts::extract(fit, type = "ev")
byDbinom <- dbinom(
  matrix(bnm$k, nrow(p), n, byrow = TRUE),
  matrix(trials, nrow(p), n, byrow = TRUE),
  p,
  log = TRUE
)
expect_equal(unname(logLikelihood), unname(byDbinom))
expect_true(all(logLikelihood[, zeroRows] == 0))
# a frequency weight w counts copies of the row: w times the row's term
weights <- rep_len(c(1, 2, 3), n)
weighted <- fitBinomial(
  cbind(k, f) ~ x1 + x2,
  bnm,
  family = binomial,
  weights = weights
)
pWeighted <- dbarts::extract(weighted, type = "ev")
expect_equal(
  unname(dbarts::extract(weighted, type = "loglik")),
  unname(sweep(
    dbinom(
      matrix(bnm$k, nrow(pWeighted), n, byrow = TRUE),
      matrix(trials, nrow(pWeighted), n, byrow = TRUE),
      pWeighted,
      log = TRUE
    ),
    2L,
    weights,
    "*"
  ))
)
# and its predictive draw is still one row's, Binomial(n, p)
set.seed(8L)
ppdWeighted <- dbarts::extract(weighted, type = "ppd")
expect_true(all(sweep(ppdWeighted, 2L, trials, "<=")))
rm(logLikelihood, p, byDbinom, weights, weighted, pWeighted, ppdWeighted)

# --- new rows: their trials from predict's weights, one where none given ---

newData <- data.frame(x1 = c(0.1, 0.5, 0.9), x2 = 0.5, k = 20, f = 20)
set.seed(9L)
single <- predict(fit, newData, type = "ppd")
expect_true(all(single %in% c(0, 1)))
set.seed(9L)
tenTrials <- predict(fit, newData, type = "ppd", weights = c(10, 10, 10))
expect_true(all(tenTrials == round(tenTrials) & tenTrials <= 10))
expect_true(max(tenTrials) > 1)
expect_true(all(
  predict(fit, newData, type = "ev") > 0 &
    predict(fit, newData, type = "ev") < 1
))
# a row of the test sample draws one trial likewise
withTest <- fitBinomial(
  cbind(k, f) ~ x1 + x2,
  bnm,
  family = binomial,
  test = newData
)
set.seed(10L)
expect_true(all(
  dbarts::extract(withTest, type = "ppd", sample = "test") %in% c(0, 1)
))
rm(newData, single, tenTrials, withTest)

# --- glm: the posterior mean probability near glm's on a large design ---

set.seed(11L)
nLarge <- 600L
large <- data.frame(x1 = runif(nLarge))
largeTrials <- rep_len(c(2L, 6L, 10L), nLarge)
large$k <- rbinom(nLarge, largeTrials, plogis(2 * large$x1 - 1))
large$f <- largeTrials - large$k
glmFit <- glm(cbind(k, f) ~ x1, family = binomial, data = large)
bartFit <- quietly(dbarts::bart(
  cbind(k, f) ~ x1,
  large,
  family = binomial,
  n.samples = 200L,
  n.burn = 100L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  seed = 12L
))
expect_true(mean(abs(fitted(bartFit) - fitted(glmFit))) < 0.03)
rm(nLarge, large, largeTrials, glmFit, bartFit)

# a hazard family's (time, status) pair is not a binomial's, under either link
for (family in c("hazard", "hazard.logistic")) {
  hazardSampler <- quietly(dbarts::dbarts(
    X,
    cbind(rep_len(1:3, n), rep_len(0:1, n)),
    family = family
  ))
  expect_false(dbarts:::samplerHoldsTrials(hazardSampler), info = family)
}
rm(family, hazardSampler)

# --- "auto" on a counts data object stays multinomial ---

autoSampler <- quietly(dbarts::dbarts(dbarts::dbartsData(X, counts = counts)))
expect_identical(autoSampler$model@family, "multinomial")
rm(autoSampler)

# --- the sampler's doors ---

control <- dbarts::dbartsControl(
  n.samples = 5L,
  n.burn = 0L,
  n.chains = 1L,
  n.threads = 1L,
  updateState = FALSE,
  seed = 1L
)
sampler <- quietly(dbarts::dbarts(
  cbind(k, f) ~ x1 + x2,
  bnm,
  family = binomial,
  control = control,
  n.trees = 5L
))
expect_true(dbarts:::samplerHoldsTrials(sampler))
expect_false(dbarts:::samplerCarriesCounts(sampler))
invisible(sampler$run())
# a successes vector against the trials held, or the pair
fewer <- pmin(bnm$k, 1L)
sampler$setResponse(fewer)
expect_identical(unname(sampler$data@counts[, 1L]), as.integer(fewer))
expect_identical(as.double(sampler$data@y), as.double(trials))
sampler$setResponse(cbind(bnm$k, bnm$f))
expect_identical(unname(sampler$data@counts), unname(counts))
expect_error(
  sampler$setResponse(trials + 1L),
  "hold more successes than trials"
)
sampler$setWeights(rep(2, n))
expect_error(
  sampler$setWeights(rep(1.5, n)),
  "logistic weights are observation counts"
)
sampler$setWeights(rep(1, n))
sampler$setOffset(rep(0.1, n))
sampler$setOffset(NULL)
sampler$setPredictor(X * 0.9)
sampler$setTestPredictor(X[1:3, , drop = FALSE])
# at the current trees, one value per row
expect_identical(length(sampler$predict(X[1:3, , drop = FALSE])), 3L)
# setData takes the successes and failures, the number of rows free
sampler$setData(dbarts::dbartsData(X[1:20, ], counts = counts[1:20, ]))
expect_identical(nrow(sampler$data@counts), 20L)
expect_error(
  sampler$setData(dbarts::dbartsData(X, as.double(bnm$k > 0))),
  "a binomial sampler's data carry its response as the n x 2 counts"
)
expect_error(
  sampler$setCategoryOffset(matrix(0, 20L, 2L)),
  "only a multinomial \\(softmax\\) sampler has one"
)
invisible(sampler$run())
# a binary sampler does not take successes and failures by setData
binarySampler <- quietly(dbarts::dbarts(
  X,
  as.double(bnm$k > 0),
  family = binomial,
  control = control,
  n.trees = 5L
))
expect_error(
  binarySampler$setData(dbarts::dbartsData(X, counts = counts)),
  "'newData' carries counts of successes and failures"
)
rm(fewer, sampler, binarySampler, control)

# --- the fallback text: a binomial row read at its trials times its weight ---

fallbackText <- character()
withCallingHandlers(
  dbarts::dbarts(
    cbind(k, f) ~ x1 + x2,
    bnm,
    family = binomial,
    weights = c(0, 2, rep(0, n - 2L)),
    n.trees = 5L,
    leaf.prior = dbarts::dbartsPriors$linear("x1"),
    control = dbarts::dbartsControl(
      n.samples = 5L,
      n.chains = 1L,
      n.threads = 1L
    )
  ),
  warning = function(condition) {
    if (inherits(condition, "dbartsScaleFallbackWarning")) {
      fallbackText <<- c(fallbackText, conditionMessage(condition))
    }
    invokeRestart("muffleWarning")
  }
)
expect_identical(
  fallbackText,
  paste0(
    "one observation is in the likelihood (1 of 40 is kept by the ",
    "active-row mask), so the center and scale of \"x1\" are taken from all ",
    "40 at their trials times their weights, a count of 0 read as 1"
  )
)
rm(fallbackText)

# --- the refusals, by their words ---

refusal <- function(expr) {
  tryCatch(
    {
      quietly(expr)
      NA_character_
    },
    error = conditionMessage
  )
}
# 11: the proportion spelling
for (family in list(binomial, stats::binomial(), "logistic")) {
  expect_identical(
    refusal(fitBinomial(
      k / pmax(trials, 1L) ~ x1,
      cbind(bnm, trials),
      family = family,
      weights = pmax(trials, 1L)
    )),
    paste0(
      "under family = binomial the response is 0/1 or cbind(successes, ",
      "failures): a proportion with the trials as weights is written ",
      "cbind(k, n - k), with k = round(p * n) from a proportion p and its ",
      "trials n"
    )
  )
}
# 12: the probit binomial, on both interfaces, xbart's and a data object
probitText <- paste0(
  "family = binomial(\"probit\") does not take a cbind(successes, failures) ",
  "response: the probit binomial is not fitted; use family = binomial, the ",
  "logit link, or a 0/1 response"
)
expect_identical(
  refusal(fitBinomial(cbind(k, f) ~ x1, bnm, family = binomial("probit"))),
  probitText
)
expect_identical(
  refusal(fitBinomial(X, counts, family = "probit")),
  probitText
)
expect_identical(
  refusal(dbarts::dbarts(
    dbarts::dbartsData(X, counts = counts),
    family = "probit"
  )),
  probitText
)
xbartArgs <- list(
  n.reps = 1L,
  n.samples = 5L,
  n.burn = 2L,
  n.trees = 5L,
  n.threads = 1L,
  verbose = FALSE
)
expect_identical(
  refusal(do.call(
    dbarts::xbart,
    c(
      list(cbind(k, f) ~ x1, bnm, family = binomial("probit")),
      xbartArgs
    )
  )),
  probitText
)
expect_identical(
  refusal(do.call(
    dbarts::xbart,
    c(list(X, counts, family = "probit"), xbartArgs)
  )),
  probitText
)
# 13: a count negative or fractional, the first five rows named
negativeRows <- which(bnm$k < 3L)
expect_true(length(negativeRows) > 5L)
expect_identical(
  refusal(fitBinomial(cbind(k - 3, f) ~ x1, bnm, family = binomial)),
  paste0(
    "a binomial response's successes and failures must be non-negative ",
    "whole numbers; row(s) ",
    paste(negativeRows[1:5], collapse = ", "),
    ", ..."
  )
)
expect_identical(
  refusal(fitBinomial(X, counts + c(0.5, rep(0, n - 1L)), family = binomial)),
  paste0(
    "a binomial response's successes and failures must be non-negative ",
    "whole numbers; row(s) 1"
  )
)
# a missing count goes to na.action
withMissing <- bnm
withMissing$k[3L] <- NA
expect_identical(
  NROW(fitBinomial(cbind(k, f) ~ x1, withMissing, family = binomial)$y),
  n - 1L
)
# 14: "auto" names the binomial first
expect_identical(
  refusal(fitBinomial(cbind(k, f) ~ x1, bnm)),
  paste0(
    "'y' is an n x 2 matrix, which family = \"auto\" leaves ambiguous - ",
    "successes and failures need family = binomial; per-category counts ",
    "need family = \"multinomial\"; a (time, status) pair needs a ",
    "survival::Surv response (which \"auto\" reads as survival) or family = ",
    "\"aft\" / \"hazard\""
  )
)
# 16: dbartsData's two-column y
expect_identical(
  refusal(dbarts::dbartsData(X, counts)),
  paste0(
    "'y' is an n x 2 matrix; dbartsData() takes a single-column response - ",
    "successes and failures and per-category counts go to dbartsData(counts ",
    "= ), a (time, status) pair to dbarts()/bart() with family = ",
    "\"aft\"/\"hazard\""
  )
)
# 17: logistic on three count columns, a data object's and a response's
threeColumns <- cbind(counts, 1L)
expect_identical(
  refusal(dbarts::dbarts(
    dbarts::dbartsData(X, counts = threeColumns),
    family = "logistic"
  )),
  paste0(
    "the data carry an n x 3 count matrix: family \"logistic\" takes two ",
    "columns, successes and failures, and family \"multinomial\" fits 3 ",
    "categories"
  )
)
expect_identical(
  refusal(fitBinomial(cbind(k, f, f) ~ x1, bnm, family = binomial)),
  paste0(
    "the response is an n x 3 count matrix: family \"logistic\" takes two ",
    "columns, successes and failures, and family \"multinomial\" fits 3 ",
    "categories"
  )
)
# 18: xbart with a matrix response other than a binomial pair, both routes
xbartText <- paste0(
  "xbart() takes a single-column response or, under family = binomial, ",
  "cbind(successes, failures); 'y' is an n x 3 matrix, and it ",
  "cross-validates neither multinomial counts nor a survival response; use ",
  "bart() instead - per-category counts need family = \"multinomial\"; a ",
  "(time, status) pair needs a survival::Surv response (which \"auto\" ",
  "reads as survival) or family = \"aft\" / \"hazard\""
)
expect_identical(
  refusal(do.call(
    dbarts::xbart,
    c(
      list(cbind(k, f, f) ~ x1, bnm, family = binomial),
      xbartArgs
    )
  )),
  xbartText
)
expect_identical(
  refusal(do.call(
    dbarts::xbart,
    c(
      list(X, threeColumns, family = binomial),
      xbartArgs
    )
  )),
  xbartText
)
# xbart takes the pair on both routes, alike
xbartFormula <- quietly(do.call(
  dbarts::xbart,
  c(
    list(
      cbind(k, f) ~ x1 + x2,
      bnm,
      family = binomial,
      loss = "log",
      seed = 2L
    ),
    xbartArgs
  )
))
xbartMatrix <- quietly(do.call(
  dbarts::xbart,
  c(
    list(X, counts, family = binomial, loss = "log", seed = 2L),
    xbartArgs
  )
))
expect_identical(as.vector(xbartMatrix), as.vector(xbartFormula))
rm(
  refusal,
  family,
  probitText,
  negativeRows,
  xbartArgs,
  withMissing,
  threeColumns
)
rm(xbartText, xbartFormula, xbartMatrix)

rm(fit, zeroRows, fitBinomial, quietly, fitArgs, counts, X, trials, bnm, n)
