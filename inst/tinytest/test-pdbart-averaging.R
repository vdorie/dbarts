# pdbart's and pd2bart's scale (type), the rows they average over (newdata,
# n.average.rows, average.weights), and variables of the data in a formula
# fit.

source(
  system.file("common", "captureWarnings.R", package = "dbarts"),
  local = TRUE
)

fitSmall <- function(...) {
  dbarts::bart(
    ...,
    n.trees = 5L,
    n.samples = 10L,
    n.burn = 5L,
    n.chains = 2L,
    n.threads = 1L,
    seed = 3L,
    keepTrees = TRUE,
    verbose = FALSE
  )
}
# the per-draw mean of predict over rows with 'values' set
rowMeanAt <- function(fit, rows, values, ...) {
  for (name in names(values)) {
    if (is.matrix(rows)) {
      rows[, name] <- values[[name]]
      next
    }
    rows[[name]] <- if (is.factor(rows[[name]])) {
      factor(values[[name]], levels = levels(rows[[name]]))
    } else {
      values[[name]]
    }
  }
  rowMeans(predict(fit, rows, ...))
}

set.seed(11)
n <- 60L
df <- data.frame(
  a = runif(n, 1, 3),
  b = rnorm(n),
  f = factor(sample(c("u", "v", "w"), n, TRUE))
)
df$y <- log(df$a) + df$b^2 + as.integer(df$f) + rnorm(n, sd = 0.3)
df$yBinary <- as.numeric(df$y > median(df$y))
df$count <- rpois(n, exp(0.5 * df$b))
df$semi <- ifelse(df$b > 0, exp(df$a / 2 + rnorm(n, sd = 0.2)), 0)
x <- as.matrix(df[, c("a", "b")])

# --- type ---
# "auto" is the link scale on every family but hurdle, and the mean response
# there; each value is the row mean of predict at that type
autoTypes <- list(
  gaussian = list(fitSmall(y ~ a + b, df), "bart"),
  student = list(fitSmall(y ~ a + b, df, family = student(3)), "bart"),
  probit = list(fitSmall(yBinary ~ a + b, df), "bart"),
  logistic = list(fitSmall(yBinary ~ a + b, df, family = "logistic"), "bart"),
  nbinom = list(fitSmall(count ~ a + b, df, family = "nbinom"), "bart"),
  hurdle = list(fitSmall(x, df$semi, family = "hurdle.lognormal"), "ev")
)
for (family in names(autoTypes)) {
  fit <- autoTypes[[family]][[1L]]
  pd <- dbarts::pdbart(fit, xind = "a", levs = list(c(1.5, 2.5)), pl = FALSE)
  expect_identical(pd$type, autoTypes[[family]][[2L]], info = family)
  expect_equal(
    pd$fd[[1L]][, 2L],
    rowMeanAt(
      fit,
      if (family == "hurdle") x else df,
      list(a = 2.5),
      type = pd$type
    ),
    info = family
  )
}
# the probability is averaged per row: transform, then average
probitFit <- autoTypes$probit[[1L]]
pdEv <- dbarts::pdbart(
  probitFit,
  xind = "a",
  levs = list(2),
  type = "ev",
  pl = FALSE
)
expect_identical(pdEv$type, "ev")
expect_equal(pdEv$fd[[1L]][, 1L], rowMeanAt(probitFit, df, list(a = 2)))
pdLink <- dbarts::pdbart(
  probitFit,
  xind = "a",
  levs = list(2),
  type = "link",
  pl = FALSE
)
expect_identical(pdLink$type, "bart")
expect_true(max(abs(pdEv$fd[[1L]] - pnorm(pdLink$fd[[1L]]))) > 1e-3)
# the hurdle's parts, and its result: the two samplers, no copied draws
hurdleFit <- autoTypes$hurdle[[1L]]
pdProb <- dbarts::pdbart(
  hurdleFit,
  xind = "a",
  levs = list(2),
  type = "prob",
  pl = FALSE
)
expect_equal(
  pdProb$fd[[1L]][, 1L],
  rowMeanAt(hurdleFit, x, list(a = 2), type = "prob")
)
expect_identical(names(pdProb$fit), c("zero", "positive"))
expect_null(pdProb$yhat.train)
expect_false(
  "fit" %in%
    names(dbarts::pdbart(
      hurdleFit,
      xind = "a",
      levs = list(2),
      keepSampler = FALSE,
      pl = FALSE
    ))
)
# and a negative-binomial fit's mean count
negbinFit <- autoTypes$nbinom[[1L]]
pdCount <- dbarts::pd2bart(
  negbinFit,
  xind = c("a", "b"),
  levs = list(2, 0),
  type = "ev",
  pl = FALSE
)
expect_equal(
  pdCount$fd[, 1L],
  rowMeanAt(negbinFit, df, list(a = 2, b = 0), type = "ev")
)
# the posterior predictive is drawn per row and averaged
pdPpd <- dbarts::pdbart(
  autoTypes$gaussian[[1L]],
  xind = "a",
  levs = list(2),
  type = "ppd",
  pl = FALSE
)
expect_identical(pdPpd$type, "ppd")
expect_equal(dim(pdPpd$fd[[1L]]), c(20L, 1L))

# "sigma" on a fit with a variance forest, and refused elsewhere
heteroFit <- fitSmall(y ~ a + b, df, variance = TRUE)
pdSigma <- dbarts::pdbart(
  heteroFit,
  xind = "a",
  levs = list(2),
  type = "sigma",
  pl = FALSE
)
expect_equal(
  pdSigma$fd[[1L]][, 1L],
  rowMeanAt(heteroFit, df, list(a = 2), type = "sigma")
)
gaussianFit <- autoTypes$gaussian[[1L]]
expect_error(
  dbarts::pdbart(gaussianFit, type = "sigma", pl = FALSE),
  "variance forest"
)
expect_error(
  dbarts::pdbart(gaussianFit, type = "prob", pl = FALSE),
  "does not take type = \"prob\""
)
expect_error(
  dbarts::pdbart(negbinFit, type = "sigma", pl = FALSE),
  "does not take type = \"sigma\""
)
expect_error(
  dbarts::pdbart(gaussianFit, type = "forest", pl = FALSE),
  "type = \"forest\""
)
# refused before anything is fit
expect_error(
  dbarts::pdbart(x, df$y, type = "evv", n.trees = -1L, pl = FALSE),
  "'type' must be one of"
)
# a sampler takes the link scale only
expect_error(
  dbarts::pdbart(gaussianFit$fit, type = "ev", pl = FALSE),
  "takes only type = \"bart\""
)
expect_error(
  dbarts::pdbart(gaussianFit$fit, newdata = df, pl = FALSE),
  "take a fit"
)
rm(autoTypes, family, fit, pd, pdEv, pdLink, pdProb, pdCount, pdPpd)
rm(heteroFit, pdSigma, hurdleFit)

# --- the plot labels, by scale and family ---
label <- function(
  type,
  family,
  call = quote(bart(formula = y ~ a, data = df))
) {
  dbarts:::pdScaleLabel(list(type = type, family = family, bartcall = call))
}
expect_identical(label("bart", "gaussian"), "y")
expect_identical(label("bart", "probit"), "probit scale")
expect_identical(label("ev", "probit"), "probability")
expect_identical(label("ev", "hurdle.lognormal"), "mean response")
expect_identical(
  label("bart", "gaussian", quote(bart(formula = x, data = yy))),
  "yy"
)
expect_identical(label(NULL, NULL), "partial-dependence")
pdLabelled <- dbarts::pd2bart(
  gaussianFit,
  xind = c("a", "b"),
  pl = FALSE
)
psFile <- tempfile(fileext = ".ps")
postscript(psFile)
plot(pdLabelled)
plot(pdLabelled, main = "mine")
dev.off()
psLines <- readLines(psFile, warn = FALSE)
expect_equal(length(grep("(Median, y)", psLines, fixed = TRUE)), 1L)
expect_equal(length(grep("(mine)", psLines, fixed = TRUE)), 1L)
unlink(psFile)
rm(label, pdLabelled, psFile, psLines)

# --- variables of the data in a formula fit ---
# xind names a variable; every term built from it moves with it
termFit <- fitSmall(y ~ log(a) + poly(a, 2) + b, df)
pdTerms <- dbarts::pdbart(termFit, xind = "a", pl = FALSE)
expect_equal(
  pdTerms$levs[[1L]],
  quantile(df$a, c(0.05, seq(0.1, 0.9, 0.1), 0.95)),
  check.attributes = FALSE
)
expect_identical(pdTerms$xlbs, "a")
expect_equal(
  pdTerms$fd[[1L]][, 3L],
  rowMeanAt(termFit, df, list(a = pdTerms$levs[[1L]][3L]))
)
expect_error(
  dbarts::pdbart(termFit, xind = 1L, pl = FALSE),
  "names variables of the data"
)
expect_error(
  dbarts::pdbart(termFit, xind = "log(a)", pl = FALSE),
  "unrecognized variables"
)
# a factor variable, by level
factorFit <- fitSmall(y ~ a + f, df)
pdFactor <- dbarts::pdbart(factorFit, xind = "f", pl = FALSE)
expect_identical(pdFactor$levs[[1L]], c("u", "v", "w"))
expect_equal(
  pdFactor$fd[[1L]][, 2L],
  rowMeanAt(factorFit, df, list(f = "v"))
)
# an offset expression in the varied variable is rebuilt with it
offsetFit <- fitSmall(y ~ a + b, df, offset = log(a))
pdOffset <- dbarts::pdbart(offsetFit, xind = "a", levs = list(2), pl = FALSE)
expect_equal(
  pdOffset$fd[[1L]][, 1L],
  rowMeanAt(offsetFit, df, list(a = 2))
)
expect_equal(
  pdOffset$fd[[1L]][, 1L],
  as.vector(colMeans(offsetFit$fit$predict(cbind(a = 2, b = df$b)))) + log(2)
)
# a fit passed in reads its stored data; without its call it needs newdata
expect_identical(pdTerms$fd, dbarts::pdbart(termFit, xind = "a", pl = FALSE)$fd)
noCallFit <- fitSmall(y ~ a + b, df, keepCall = FALSE)
expect_error(
  dbarts::pdbart(noCallFit, xind = "a", pl = FALSE),
  "give the rows as 'newdata'"
)
expect_equal(
  dim(dbarts::pdbart(noCallFit, xind = "a", newdata = df, pl = FALSE)$fd[[1L]]),
  c(20L, 11L)
)
# rows dropped for a missing response are not averaged
dfMissing <- df[1:50, ]
dfMissing$y[c(2, 9, 17, 30, 44)] <- NA
missingFit <- fitSmall(y ~ a + b, dfMissing)
pdMissing <- dbarts::pdbart(missingFit, xind = "a", levs = list(2), pl = FALSE)
expect_equal(
  pdMissing$fd[[1L]][, 1L],
  rowMeanAt(missingFit, dfMissing[!is.na(dfMissing$y), ], list(a = 2))
)
# and a data call reads the data it was given
pdData <- dbarts::pdbart(
  y ~ log(a) + poly(a, 2) + b,
  df,
  xind = "a",
  n.trees = 5L,
  n.samples = 10L,
  n.burn = 5L,
  n.chains = 2L,
  n.threads = 1L,
  seed = 3L,
  verbose = FALSE,
  pl = FALSE
)
expect_identical(pdData$fd, pdTerms$fd)
rm(termFit, pdTerms, factorFit, pdFactor, offsetFit, pdOffset, noCallFit)
rm(dfMissing, missingFit, pdMissing, pdData)

# --- newdata ---
# a subgroup's rows, with its grid
subgroup <- df[df$f == "u", ]
pdGroup <- dbarts::pdbart(
  gaussianFit,
  xind = "a",
  newdata = subgroup,
  pl = FALSE
)
expect_equal(
  pdGroup$levs[[1L]],
  quantile(subgroup$a, c(0.05, seq(0.1, 0.9, 0.1), 0.95)),
  check.attributes = FALSE
)
expect_equal(
  pdGroup$fd[[1L]][, 4L],
  rowMeanAt(gaussianFit, subgroup, list(a = pdGroup$levs[[1L]][4L]))
)
# a factor shows every level the fit knows, whatever the subgroup holds
factorFit <- fitSmall(y ~ a + f, df)
pdLevels <- dbarts::pdbart(
  factorFit,
  xind = "f",
  newdata = subgroup,
  pl = FALSE
)
expect_identical(pdLevels$levs[[1L]], c("u", "v", "w"))
# an unseen level in the rows is refused, as predict refuses it
unseen <- subgroup
unseen$f <- factor(c("zz", as.character(unseen$f[-1L])))
expect_error(
  dbarts::pdbart(factorFit, xind = "a", newdata = unseen, pl = FALSE),
  "levels not present in the training data"
)
# a missing predictor is kept, on a fit that learned a route for it
dfNA <- df
dfNA$b[c(4L, 20L, 33L)] <- NA
naFit <- fitSmall(y ~ a + b, dfNA)
withMissing <- dfNA[dfNA$f == "u", ]
withMissing$b[1:2] <- NA
pdKept <- dbarts::pdbart(
  naFit,
  xind = "a",
  levs = list(2),
  newdata = withMissing,
  pl = FALSE
)
expect_false(anyNA(pdKept$fd[[1L]]))
expect_equal(
  pdKept$fd[[1L]][, 1L],
  rowMeanAt(naFit, withMissing, list(a = 2))
)
# a plain-vector offset cannot be evaluated on other rows
offsetValues <- rnorm(n)
vectorFit <- fitSmall(x, df$y, offset = offsetValues)
expect_error(
  dbarts::pdbart(vectorFit, newdata = x[1:10, ], pl = FALSE),
  "write the offset as a column of the data"
)
expect_error(
  dbarts::pdbart(gaussianFit, newdata = df, n.average.rows = 5L, pl = FALSE),
  "cannot be given with 'newdata'"
)
rm(subgroup, pdGroup, pdLevels, unseen, dfNA, naFit, withMissing, pdKept)

# --- n.average.rows and average.weights ---
# a subsample is reproducible and is the rows it drew
set.seed(5)
pdSample <- dbarts::pdbart(
  gaussianFit,
  xind = "a",
  levs = list(2),
  n.average.rows = 12L,
  pl = FALSE
)
set.seed(5)
sampled <- sort(sample.int(n, 12L))
expect_equal(
  pdSample$fd,
  dbarts::pdbart(
    gaussianFit,
    xind = "a",
    levs = list(2),
    newdata = df[sampled, ],
    pl = FALSE
  )$fd
)
expect_error(
  dbarts::pdbart(gaussianFit, n.average.rows = 0, pl = FALSE),
  "'n.average.rows' must be a positive whole number"
)
# weights: a weighted mean, unchanged by scale
weights <- runif(n)
pdWeighted <- dbarts::pdbart(
  gaussianFit,
  xind = "a",
  levs = list(2),
  average.weights = weights,
  pl = FALSE
)
atTwo <- df
atTwo$a <- 2
expect_equal(
  pdWeighted$fd[[1L]][, 1L],
  drop(predict(gaussianFit, atTwo) %*% weights) / sum(weights)
)
expect_equal(
  pdWeighted$fd,
  dbarts::pdbart(
    gaussianFit,
    xind = "a",
    levs = list(2),
    average.weights = 10 * weights,
    pl = FALSE
  )$fd
)
for (bad in list(
  weights[-1L],
  replace(weights, 1L, -1),
  replace(weights, 1L, NA)
)) {
  expect_error(
    dbarts::pdbart(gaussianFit, average.weights = bad, pl = FALSE),
    "'average.weights'"
  )
}
expect_error(
  dbarts::pdbart(gaussianFit, average.weights = numeric(n), pl = FALSE),
  "must not all be zero"
)
# an entry for a row the fit gives weight 0 is ignored
fitWeights <- rep_len(c(1, 0, 2), n)
zeroFit <- withCallingHandlers(
  dbarts::bart(
    y ~ a + b,
    df,
    weights = fitWeights,
    n.trees = 5L,
    n.samples = 10L,
    n.burn = 5L,
    n.chains = 2L,
    n.threads = 1L,
    seed = 3L,
    keepTrees = TRUE,
    verbose = FALSE
  ),
  warning = function(w) {
    if (grepl("'weights' of 0", conditionMessage(w), fixed = TRUE)) {
      invokeRestart("muffleWarning")
    }
  }
)
pdZero <- function(w) {
  dbarts::pdbart(
    zeroFit,
    xind = "a",
    levs = list(2),
    average.weights = w,
    pl = FALSE
  )$fd
}
expect_equal(pdZero(weights), pdZero(replace(weights, fitWeights == 0, 50)))
# weights with a subsample: the sampled rows keep theirs, renormalized
set.seed(9)
pdBoth <- dbarts::pdbart(
  gaussianFit,
  xind = "a",
  levs = list(2),
  n.average.rows = 15L,
  average.weights = weights,
  pl = FALSE
)
set.seed(9)
sampled <- sort(sample.int(n, 15L))
expect_equal(
  pdBoth$fd[[1L]][, 1L],
  drop(predict(gaussianFit, atTwo[sampled, ]) %*% weights[sampled]) /
    sum(weights[sampled])
)
# offsets under a subsample: an expression is evaluated on the rows, a plain
# vector gives each row its own value
expressionFit <- fitSmall(y ~ a + b, df, offset = b / 2)
set.seed(4)
pdExpression <- dbarts::pdbart(
  expressionFit,
  xind = "a",
  levs = list(2),
  n.average.rows = 10L,
  pl = FALSE
)
set.seed(4)
sampled <- sort(sample.int(n, 10L))
expect_equal(
  pdExpression$fd[[1L]][, 1L],
  rowMeans(predict(expressionFit, atTwo[sampled, ]))
)
set.seed(4)
pdVector <- dbarts::pdbart(
  vectorFit,
  xind = 1L,
  levs = list(2),
  n.average.rows = 10L,
  pl = FALSE
)
xAtTwo <- x
xAtTwo[, 1L] <- 2
expect_equal(
  pdVector$fd[[1L]][, 1L],
  rowMeans(predict(
    vectorFit,
    xAtTwo[sampled, ],
    offset = offsetValues[sampled]
  ))
)
expect_equal(
  pdVector$fd[[1L]][, 1L],
  as.vector(colMeans(vectorFit$fit$predict(xAtTwo[sampled, ]))) +
    mean(offsetValues[sampled])
)
rm(pdSample, sampled, weights, pdWeighted, atTwo, bad, fitWeights, zeroFit)
rm(pdZero, pdBoth, expressionFit, pdExpression, pdVector, xAtTwo)
rm(offsetValues, vectorFit)

# --- pd2bart's two-predictor shortcut ---
# it ignores the averaging arguments, saying so
warnings.shortcut <- captureWarnings(
  pdShort <- dbarts::pd2bart(
    gaussianFit,
    levs = list(2, 0),
    newdata = df[1:5, ],
    pl = FALSE
  )
)
expect_equal(length(warnings.shortcut), 1L)
expect_true(grepl("'newdata'", conditionMessage(warnings.shortcut[[1L]])))
expect_equal(
  pdShort$fd[, 1L],
  rowMeanAt(gaussianFit, df, list(a = 2, b = 0))
)
# but a posterior predictive draw takes the general route, averaging each
# row's own draw
set.seed(2)
pdDraw <- dbarts::pd2bart(
  gaussianFit,
  levs = list(2, 0),
  type = "ppd",
  pl = FALSE
)
set.seed(2)
expect_equal(
  pdDraw$fd[, 1L],
  rowMeanAt(gaussianFit, df, list(a = 2, b = 0), type = "ppd")
)
rm(warnings.shortcut, pdShort, pdDraw)

rm(fitSmall, rowMeanAt, n, df, x, probitFit, negbinFit, gaussianFit)
rm(factorFit)
