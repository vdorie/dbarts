# weights.type and the frequency reading: its resolution per family, its
# refusals, the starting sigma read off the rows written out, and the
# posterior predictive draws a weight asks for on every family (dec-B565,
# dec-B566, dec-B580, dec-B594, dec-B595, dec-B596)

set.seed(71)
n <- 40L
x <- matrix(rnorm(n * 2L), n, dimnames = list(NULL, c("x1", "x2")))
y <- x[, 1L] + rnorm(n)
w <- as.double(sample(1:3, n, replace = TRUE))
mask <- rep(c(0, 1, 1, 1), length.out = n)
df <- data.frame(y = y, x1 = x[, 1L], x2 = x[, 2L], w = w, m = mask)
df$yb <- as.double(df$y > 0)
df$yc <- rpois(n, exp(0.5 * df$x1))
df$yo <- factor(cut(df$y, 3L), ordered = TRUE)
df$ym <- factor(sample(letters[1:3], n, replace = TRUE))
expanded <- df[rep(seq_len(n), w), ]
survivalResponse <- cbind(ceiling(exp(y)), rep_len(c(1, 0, 1), n))

# the call is rebuilt in the caller's frame, so a weights column is read from
# the data as written there
quick <- function(..., n.samples = 20L, n.chains = 2L, keepTrees = TRUE) {
  call <- match.call()
  call[[1L]] <- quote(bart)
  call$n.samples <- n.samples
  call$n.chains <- n.chains
  call$keepTrees <- keepTrees
  call$n.trees <- 10L
  call$n.burn <- 10L
  call$n.threads <- 1L
  call$seed <- 3L
  call$verbose <- FALSE
  eval(call, parent.frame())
}
specOf <- function(data, ...) {
  dbartsSpec(
    data,
    dbartsControl(n.samples = 5L, n.chains = 1L, n.threads = 1L),
    ...,
    n.trees = 10L
  )
}

# ---- the refusals, by their words -------------------------------------------

expect_error(
  dbarts(y ~ x1, df, weights.type = "count"),
  "'weights.type' must be NULL, \"precision\" or \"frequency\"",
  fixed = TRUE
)
expect_error(
  bart(y ~ x1, df, weights.type = NA_character_, verbose = FALSE),
  "'weights.type' must be NULL, \"precision\" or \"frequency\"",
  fixed = TRUE
)
expect_error(
  xbart(y ~ x1, df, weights.type = c("precision", "frequency")),
  "'weights.type' must be NULL, \"precision\" or \"frequency\"",
  fixed = TRUE
)
# 2: a stated frequency on a membership family takes 0 and 1 only
for (family in c("probit", "ordinal")) {
  response <- if (family == "probit") df$yb else df$yo
  expect_error(
    dbarts(
      x,
      response,
      weights = w,
      family = family,
      weights.type = "frequency"
    ),
    paste0(
      "weights.type = \"frequency\" is not taken under family = \"",
      family,
      "\" for weights other than 0 and 1"
    ),
    fixed = TRUE,
    info = family
  )
}
expect_silent(dbarts(
  x,
  df$yb,
  weights = mask,
  family = "probit",
  weights.type = "frequency"
))
# 3: precision where there is no residual scale, with or without weights
for (family in c("logistic", "nbinom", "probit", "ordinal", "multinomial")) {
  response <- switch(
    family,
    logistic = ,
    probit = df$yb,
    nbinom = df$yc,
    ordinal = df$yo,
    multinomial = df$ym
  )
  expect_error(
    dbarts(x, response, family = family, weights.type = "precision"),
    paste0(
      "weights.type = \"precision\" needs a residual scale, and family = \"",
      family,
      "\" has none"
    ),
    fixed = TRUE,
    info = family
  )
}
expect_error(
  bart(
    x,
    survivalResponse,
    family = "hazard.logistic",
    weights.type = "precision",
    verbose = FALSE
  ),
  "needs a residual scale"
)
# 4: aft takes no weights, so no reading; the hurdle's 0/1 weights are rows in
expect_error(
  dbarts(x, survivalResponse, family = "aft", weights.type = "frequency"),
  "weights.type = \"frequency\" is not taken under family = \"aft\", which takes no weights",
  fixed = TRUE
)
yHurdle <- ifelse(df$yb > 0, exp(y), 0)
expect_error(
  bart(
    x,
    yHurdle,
    family = "hurdle.lognormal",
    weights.type = "precision",
    verbose = FALSE
  ),
  paste0(
    "weights.type = \"precision\" is not taken under family = ",
    "\"hurdle.lognormal\": its weights of 0 and 1 mark rows in and out of ",
    "both parts"
  ),
  fixed = TRUE
)
# 5, 6, 7: the count families take non-negative whole numbers
expect_error(
  dbarts(x, df$yb, weights = w + 0.5, family = "logistic"),
  "logistic weights are observation counts and must be non-negative whole numbers",
  fixed = TRUE
)
expect_error(
  dbarts(x, df$yc, weights = w + 0.5, family = "nbinom"),
  paste0(
    "nbinom weights are observation counts and must be non-negative whole ",
    "numbers; an exposure belongs in the offset as a log-exposure term"
  ),
  fixed = TRUE
)
expect_error(
  dbarts(x, df$ym, weights = w + 0.5, family = "multinomial"),
  "multinomial weights are observation counts and must be non-negative whole numbers",
  fixed = TRUE
)
# 10: the reading is fixed when the sampler is created
sampler <- dbarts(y ~ x1 + x2, df, weights = w, weights.type = "frequency")
other <- sampler$model
other@weights.type <- "precision"
expect_error(
  sampler$setModel(other),
  paste0(
    "$setModel cannot change weights.type, which is fixed when the sampler ",
    "is created; create a new sampler"
  ),
  fixed = TRUE
)
# a model stating none reads the sampler's own and passes
same <- sampler$model
same@weights.type <- NA_character_
expect_silent(sampler$setModel(same))
expect_identical(sampler$model@weights.type, "frequency")
rm(sampler, other, same)

# ---- the resolution, row by row ----------------------------------------------

resolved <- function(response, family, ...) {
  specOf(dbartsData(x, response, ...), family = family)$model@weights.type
}
expect_identical(resolved(y, "gaussian"), "precision")
expect_identical(resolved(y, "gaussian", weights = w), "precision")
expect_identical(resolved(y, "student"), "precision")
expect_identical(resolved(df$yb, "probit"), "frequency")
expect_identical(resolved(df$yb, "logistic"), "frequency")
expect_identical(resolved(df$yc, "nbinom"), "frequency")
expect_identical(resolved(df$yo, "ordinal"), "frequency")
expect_identical(
  specOf(
    dbartsData(x, y),
    family = "gaussian",
    weights.type = "frequency"
  )$model@weights.type,
  "frequency"
)
expect_identical(
  dbarts(
    x,
    y,
    variance = varianceForest(n.trees = 5L),
    weights = w
  )$model@weights.type,
  "precision"
)
expect_identical(
  dbarts(x, survivalResponse, family = "aft")$model@weights.type,
  "precision"
)
expect_identical(
  dbarts(x, df$ym, family = "multinomial")$model@weights.type,
  "frequency"
)
expect_identical(
  dbarts(x, survivalResponse, family = "hazard")$model@weights.type,
  "frequency"
)
# a model made by hand states none and resolves at the sampler's creation
spec <- specOf(dbartsData(x, y, weights = w), family = "gaussian")
spec$model@weights.type <- NA_character_
handMade <- new("dbartsSampler", spec$control, spec$model, spec$data)
expect_identical(handMade$model@weights.type, "precision")
rm(spec, handMade)
# a gaussian fit records a stated frequency; every other fit's content is as
# it was
expect_identical(
  quick(y ~ x1 + x2, df, weights = w, weights.type = "frequency")$weights.type,
  "frequency"
)
expect_null(quick(y ~ x1 + x2, df, weights = w)$weights.type)

# ---- the starting sigma, read off the rows written out -----------------------

for (rule in c("dense", "sparse")) {
  frequencyFit <- dbarts(
    y ~ x1 + x2,
    df,
    weights = w,
    weights.type = "frequency",
    sigest = rule
  )
  expandedFit <- dbarts(y ~ x1 + x2, expanded, sigest = rule)
  expect_equal(
    frequencyFit$data@sigma,
    expandedFit$data@sigma,
    tolerance = 1e-12,
    info = rule
  )
  # the sampler starts where the chi prior is anchored, at the estimate
  expect_equal(
    frequencyFit$getSigmas(),
    expandedFit$getSigmas(),
    tolerance = 1e-12,
    info = rule
  )
}
# precision is the weighted regression, its df the rows less the rank
precisionFit <- dbarts(y ~ x1 + x2, df, weights = w, sigest = "dense")
expect_equal(
  precisionFit$data@sigma,
  summary(lm(y ~ x1 + x2, df, weights = w))$sigma,
  tolerance = 1e-12
)
# a forest's basis enters the regression, weighted the same way
df$z <- rnorm(n)
expanded$z <- df$z[rep(seq_len(n), w)]
for (rule in c("dense", "sparse")) {
  expect_equal(
    dbarts(
      y ~ forest(x1) + forest(x2, basis = z),
      df,
      weights = w,
      weights.type = "frequency",
      sigest = rule
    )$data@sigma,
    dbarts(
      y ~ forest(x1) + forest(x2, basis = z),
      expanded,
      sigest = rule
    )$data@sigma,
    tolerance = 1e-12,
    info = paste("basis", rule)
  )
}
# a weight of 0 is a row out, as in the expansion
w0 <- w
w0[1:5] <- 0
expect_equal(
  dbarts(y ~ x1 + x2, df, weights = w0, weights.type = "frequency")$data@sigma,
  dbarts(y ~ x1 + x2, df[rep(seq_len(n), w0), ])$data@sigma,
  tolerance = 1e-12
)
# the regression is defined by the distinct rows, as under precision: three
# rows and three columns, counted 2, 3 and 4, leave it none, so the sd of the
# rows written out stands in where the written-out twin's exact fit lands on
# the floor, and the fit and its summary say which sd it is
small <- data.frame(y = c(1, 3, 2), a = c(0, 1, 0), b = c(1, 1, 0))
smallCounts <- c(2, 3, 4)
for (rule in c("dense", "sparse")) {
  expect_equal(
    dbarts(
      y ~ a + b,
      small,
      weights = smallCounts,
      weights.type = "frequency",
      sigest = rule
    )$data@sigma,
    sd(rep(small$y, smallCounts)),
    tolerance = 1e-12,
    info = rule
  )
}
expect_true(
  dbarts(y ~ a + b, small[rep(1:3, smallCounts), ])$data@sigma < 1e-7
)
saturated <- quick(
  y ~ a + b,
  small,
  weights = smallCounts,
  weights.type = "frequency",
  n.chains = 1L
)
expect_true(isTRUE(attr(saturated$sigest.fallback, "weighted")))
countLine <- paste0(
  "(Starting sigma: the sd of the response with each row written out as ",
  "often as its weight;"
)
expect_true(any(grepl(
  countLine,
  utils::capture.output(print(summary(saturated))),
  fixed = TRUE
)))
# counts of 1 write out nothing: the unweighted fit's words
expect_true(any(grepl(
  "(Starting sigma: the sd of the response;",
  utils::capture.output(print(summary(quick(
    y ~ a + b,
    small,
    weights = c(1, 1, 1),
    weights.type = "frequency",
    n.chains = 1L
  )))),
  fixed = TRUE
)))
# the verbose line of bart() and of dbarts() names the count stand-in
verboseLine <- paste0(
  "starting sigma is the sd of the response with each row written out as ",
  "often as its weight"
)
verboseLines <- character()
invisible(utils::capture.output(withCallingHandlers(
  bart(
    y ~ a + b,
    small,
    weights = smallCounts,
    weights.type = "frequency",
    n.samples = 5L,
    n.burn = 5L,
    n.trees = 5L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = TRUE
  ),
  message = function(m) {
    verboseLines <<- c(verboseLines, conditionMessage(m))
    invokeRestart("muffleMessage")
  }
)))
expect_true(any(grepl(verboseLine, verboseLines, fixed = TRUE)))
dbartsLines <- character()
conditions <- character()
invisible(utils::capture.output(withCallingHandlers(
  dbarts(
    y ~ a + b,
    small,
    weights = smallCounts,
    weights.type = "frequency",
    control = dbartsControl(verbose = TRUE, n.chains = 1L, n.threads = 1L)
  ),
  message = function(m) {
    dbartsLines <<- c(dbartsLines, conditionMessage(m))
    invokeRestart("muffleMessage")
  },
  dbartsSigmaFallback = function(condition) {
    conditions <<- c(conditions, conditionMessage(condition))
  }
)))
expect_true(any(grepl(verboseLine, dbartsLines, fixed = TRUE)))
expect_identical(
  conditions,
  paste0(
    "the starting sigma is ",
    "the sd of the ",
    "response with each row written out as often as its weight"
  )
)
# fractional counts whose rows in exceed the rank but whose sum falls below
# it leave the dense regression no degrees of freedom: the stand-in, silently
set.seed(4)
fracX <- matrix(rnorm(40), 20)
fracY <- rnorm(20)
fracWarnings <- character()
fracSigma <- withCallingHandlers(
  dbarts(
    fracX,
    fracY,
    weights = rep(0.1, 20),
    weights.type = "frequency",
    sigest = "dense"
  )$data@sigma,
  warning = function(condition) {
    fracWarnings <<- c(fracWarnings, conditionMessage(condition))
    invokeRestart("muffleWarning")
  }
)
expect_identical(fracWarnings, character())
expect_equal(
  fracSigma,
  sqrt(sum(0.1 * (fracY - mean(fracY))^2) / (sum(rep(0.1, 20)) - 1)),
  tolerance = 1e-12
)
# a fold's estimate, the routine xbart runs on each fold's training rows
trainRows <- sort(sample(n, 30L))
foldExpanded <- rep(trainRows, w[trainRows])
expect_equal(
  dbarts:::residualStandardError(
    y[trainRows],
    cbind(1, x[trainRows, ]),
    w[trainRows],
    NULL,
    frequency = TRUE
  ),
  summary(lm(y[foldExpanded] ~ x[foldExpanded, ]))$sigma,
  tolerance = 1e-12
)
design <- dbarts:::startingSigmaDesign(x, "sparse", NULL)
expect_equal(
  dbarts:::sparseResidualStandardError(
    y,
    design,
    w,
    NULL,
    rows = trainRows,
    frequency = TRUE
  ),
  summary(lm(y[foldExpanded] ~ x[foldExpanded, ]))$sigma,
  tolerance = 1e-12
)
# and xbart runs under the reading
expect_silent(xbart(
  y ~ x1 + x2,
  df,
  weights = w,
  weights.type = "frequency",
  n.reps = 1L,
  n.samples = 10L,
  n.burn = c(5L, 2L),
  n.trees = 5L,
  k = 2,
  n.threads = 1L,
  verbose = FALSE
))
# a fold's training rows saturated on their distinct rows, where all rows are
# not: twelve rows, ten coefficients and three rows held out. Each fold takes
# the sd of its training rows written out, by both routines, and one line
# under verbose says so
xFold <- matrix(rnorm(12L * 9L), 12L)
yFold <- rnorm(12L)
wFold <- as.double(rep(2:4, 4L))
foldStandIns <- function(rule) {
  original <- dbarts:::floorSigmaEstimate
  calls <- list()
  recording <- function(sigma, residual, weights = NULL, frequency = FALSE) {
    result <- original(sigma, residual, weights, frequency)
    calls[[length(calls) + 1L]] <<- list(
      residual = residual,
      weights = weights,
      result = result
    )
    result
  }
  assignInNamespace("floorSigmaEstimate", recording, "dbarts")
  on.exit(assignInNamespace("floorSigmaEstimate", original, "dbarts"))
  lines <- character()
  withCallingHandlers(
    invisible(utils::capture.output(xbart(
      xFold,
      yFold,
      weights = wFold,
      weights.type = "frequency",
      sigest = rule,
      n.reps = 1L,
      n.test = 4L,
      n.samples = 5L,
      n.burn = 5L,
      n.trees = 5L,
      n.threads = 1L,
      seed = 3L,
      verbose = TRUE
    ))),
    message = function(m) {
      lines <<- c(lines, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  folds <- Filter(function(call) length(call$residual) < 12L, calls)
  list(lines = lines, folds = folds)
}
for (rule in c("dense", "sparse")) {
  made <- foldStandIns(rule)
  expect_equal(length(made$folds), 4L, info = rule)
  for (fold in made$folds) {
    expect_equal(
      fold$result,
      sd(rep(fold$residual, fold$weights)),
      tolerance = 1e-12,
      info = rule
    )
  }
  expect_equal(
    sum(grepl(
      paste0(
        "starting sigma is the sd of the response with each row written ",
        "out as often as its weight in 4 of 4"
      ),
      made$lines,
      fixed = TRUE
    )),
    1L,
    info = rule
  )
}
rm(
  frequencyFit,
  expandedFit,
  precisionFit,
  small,
  smallCounts,
  design,
  saturated,
  countLine,
  xFold,
  yFold,
  wFold,
  foldStandIns,
  made
)

# ---- a training row's draw is one observation --------------------------------

# frequency: one observation at sd sigma, whatever the row's weight
freqFit <- quick(
  y ~ x1 + x2,
  df,
  weights = w,
  weights.type = "frequency",
  n.samples = 2000L,
  n.chains = 1L
)
residual <- extract(freqFit, "ppd") - extract(freqFit, "ev")
scaled <- residual / as.vector(freqFit$sigma)
heavy <- w == 3
expect_true(abs(sd(scaled[, heavy]) - 1) < 0.05)
expect_true(abs(sd(scaled[, !heavy]) - 1) < 0.05)
# the log-likelihood is w times one observation's
expect_equal(
  extract(freqFit, "loglik"),
  t(
    w *
      t(dnorm(
        matrix(y, nrow(residual), n, byrow = TRUE),
        extract(freqFit, "ev"),
        as.vector(freqFit$sigma),
        log = TRUE
      ))
  ),
  tolerance = 1e-10,
  check.attributes = FALSE
)
# predict takes 1, a new observation, and 0, a row held out
expect_error(
  predict(freqFit, type = "ppd", weights = 2),
  paste0(
    "on a fit with weights.type = \"frequency\" predict takes weights of 1, ",
    "a new observation, and 0, a row held out"
  ),
  fixed = TRUE
)
expect_true(all(is.na(predict(freqFit, df[1:3, ], type = "ppd", weights = 0))))
rm(freqFit, residual, scaled, heavy)

# logistic: one 0/1 draw at a training row of any count
logisticFit <- quick(yb ~ x1 + x2, df, weights = w, family = "logistic")
logisticDraws <- extract(logisticFit, "ppd")
expect_true(all(logisticDraws %in% c(0, 1)))
# predict's weights are new trials
set.seed(4)
trials <- predict(logisticFit, type = "ppd", weights = w)
expect_true(all(trials >= 0 & trials <= rep(w, each = nrow(trials))))
expect_true(any(trials > 1))

# ---- a row of weight 0 has NA draws on every family, the stream unmoved -------

zeroWeights <- w
zeroWeights[c(2L, 7L)] <- 0
zeroMask <- rep(1, n)
zeroMask[c(2L, 7L)] <- 0
naCases <- list(
  gaussian = list(formula = y ~ x1 + x2, weights = zeroWeights),
  frequency = list(
    formula = y ~ x1 + x2,
    weights = zeroWeights,
    weights.type = "frequency"
  ),
  probit = list(formula = yb ~ x1 + x2, weights = zeroMask),
  logistic = list(
    formula = yb ~ x1 + x2,
    weights = zeroWeights,
    family = "logistic"
  ),
  nbinom = list(
    formula = yc ~ x1 + x2,
    weights = zeroWeights,
    family = "nbinom"
  ),
  ordinal = list(
    formula = yo ~ x1 + x2,
    weights = zeroMask,
    family = "ordinal"
  ),
  multinomial = list(
    formula = ym ~ x1 + x2,
    weights = zeroWeights,
    family = "multinomial"
  )
)
df$zeroWeights <- zeroWeights
df$zeroMask <- zeroMask
for (name in names(naCases)) {
  case <- naCases[[name]]
  arguments <- c(
    list(case$formula, df),
    case[setdiff(names(case), c("formula", "weights"))],
    list(
      weights = if (identical(case$weights, zeroMask)) {
        quote(zeroMask)
      } else {
        quote(zeroWeights)
      }
    )
  )
  fit <- do.call(quick, arguments)
  set.seed(5)
  draws <- extract(fit, "ppd")
  rowMargin <- length(dim(draws)) - if (name == "multinomial") 0L else 0L
  drawsOut <- if (length(dim(draws)) == 3L) {
    draws[, c(2L, 7L), ]
  } else {
    draws[, c(2L, 7L)]
  }
  drawsIn <- if (length(dim(draws)) == 3L) {
    draws[, -c(2L, 7L), ]
  } else {
    draws[, -c(2L, 7L)]
  }
  expect_true(all(is.na(drawsOut)), info = name)
  expect_false(anyNA(drawsIn), info = name)
  # the fitted value is still reported
  expect_false(anyNA(fitted(fit)), info = name)
  # and stating the row's weight draws it as a new observation
  set.seed(5)
  given <- predict(fit, type = "ppd", weights = 1)
  expect_false(anyNA(given), info = name)
  rm(fit, draws, drawsOut, drawsIn, given, rowMargin, arguments, case)
}
# the stream: the row is drawn and then overwritten, so every other row's draw
# is the draw a call that keeps it makes
gaussianZero <- quick(y ~ x1 + x2, df, weights = zeroWeights)
set.seed(6)
withZero <- extract(gaussianZero, "ppd")
set.seed(6)
keptRow <- predict(
  gaussianZero,
  type = "ppd",
  weights = ifelse(zeroWeights == 0, 1, zeroWeights)
)
expect_identical(withZero[, -c(2L, 7L)], keptRow[, -c(2L, 7L)])
probitZero <- quick(yb ~ x1 + x2, df, weights = zeroMask)
set.seed(6)
withZero <- extract(probitZero, "ppd")
set.seed(6)
probitAll <- quick(yb ~ x1 + x2, df)
expect_identical(dim(withZero), dim(extract(probitAll, "ppd")))
rm(gaussianZero, withZero, keptRow, probitZero, probitAll)

# ---- nbinom: predict's weights are new observations summed --------------------

nbinomFit <- quick(yc ~ x1 + x2, df, family = "nbinom", n.chains = 1L)
mu <- extract(nbinomFit, "ev")
r <- as.vector(extract(nbinomFit, "shape"))
set.seed(7)
summed <- predict(nbinomFit, type = "ppd", weights = 3)
set.seed(7)
byHand <- matrix(
  rnbinom(length(mu), size = 3 * rep(r, n), mu = 3 * as.vector(mu)),
  nrow(mu)
)
expect_equal(summed, byHand, check.attributes = FALSE)
# a 0 is no observations, a count of 0, at newdata and at the training rows
expect_true(all(predict(nbinomFit, df[1:4, ], type = "ppd", weights = 0) == 0))
expect_true(all(predict(nbinomFit, type = "ppd", weights = 0) == 0))
expect_error(
  predict(nbinomFit, type = "ppd", weights = 0.5),
  paste0(
    "the posterior predictive 'weights' of a nbinom fit count new ",
    "observations, whose draws are summed, and must be non-negative integers"
  ),
  fixed = TRUE
)
# counts at the fit: a row of weight 2 is the row twice
countedFit <- quick(yc ~ x1 + x2, df, weights = w, family = "nbinom")
expect_identical(countedFit$weights, w)
rm(nbinomFit, mu, r, summed, byHand, countedFit)

# ---- ordinal and multinomial: counts per category over w trials (dec-B594) ---

ordinalFit <- quick(yo ~ x1 + x2, df, family = "ordinal")
categoryFit <- quick(ym ~ x1 + x2, df, family = "multinomial")
counts <- t(sapply(seq_len(n), function(i) {
  tabulate(sample(3L, sample(1:4, 1L), replace = TRUE), 3L)
}))
colnames(counts) <- c("a", "b", "c")
countFit <- bart(
  x,
  counts,
  family = "multinomial",
  n.trees = 10L,
  n.burn = 10L,
  n.samples = 20L,
  n.chains = 2L,
  n.threads = 1L,
  keepTrees = TRUE,
  seed = 3L,
  verbose = FALSE
)
trialWeights <- rep(c(1, 2, 7, 0), length.out = n)
# a (chains x) draws x rows x K probability array's counts over 'trials': each
# category's count a binomial of the trials left at its share of the
# probability left, every draw of a row at once; a row of no trials is drawn
# at one and reported as zeros
sequentialBinomials <- function(probabilities, trials) {
  d <- dim(probabilities)
  K <- d[length(d)]
  numDraws <- prod(d) %/% (d[length(d) - 1L] * K)
  flat <- matrix(probabilities, prod(d) %/% K, K)
  remaining <- rep(as.integer(trials), each = numDraws)
  none <- remaining == 0L
  remaining[none] <- 1L
  left <- rowSums(flat)
  result <- matrix(0L, nrow(flat), K)
  for (k in seq_len(K - 1L)) {
    share <- ifelse(left > 0, pmin(1, flat[, k] / left), 0)
    result[, k] <- rbinom(length(remaining), remaining, share)
    remaining <- remaining - result[, k]
    left <- pmax(0, left - flat[, k])
  }
  result[, K] <- remaining
  result[none, ] <- 0L
  array(result, d)
}
newRows <- df[1:8, ]
for (fitName in c("ordinal", "category", "count")) {
  fit <- switch(
    fitName,
    ordinal = ordinalFit,
    category = categoryFit,
    count = countFit
  )
  newdata <- if (fitName == "count") x[1:8, ] else newRows
  K <- 3L
  for (combine in c(TRUE, FALSE)) {
    info <- paste(fitName, combine)
    # without weights: a category where the fit models categories, a count
    # vector at a count fit's training rows, a one-trial one at its new rows
    plain <- predict(fit, type = "ppd", combineChains = combine)
    atNew <- predict(fit, newdata, type = "ppd", combineChains = combine)
    if (fitName == "count") {
      expect_identical(dim(plain)[length(dim(plain))], K, info = info)
      expect_true(
        all(apply(atNew, seq_len(length(dim(atNew)) - 1L), sum) == 1),
        info = info
      )
    } else {
      expect_true(all(plain %in% seq_len(K)), info = info)
      expect_true(all(atNew %in% seq_len(K)), info = info)
    }
    # with weights: counts summing to w, at the training rows and at newdata
    weighted <- predict(
      fit,
      type = "ppd",
      weights = trialWeights,
      combineChains = combine
    )
    margins <- seq_len(length(dim(weighted)) - 1L)
    totals <- apply(weighted, margins, sum)
    rowAxis <- length(dim(totals))
    expect_identical(dim(weighted)[length(dim(weighted))], K, info = info)
    expect_true(
      all(sweep(totals, rowAxis, trialWeights) == 0),
      info = info
    )
    weightedNew <- predict(
      fit,
      newdata,
      type = "ppd",
      weights = trialWeights[1:8],
      combineChains = combine
    )
    totalsNew <- apply(weightedNew, seq_len(length(dim(weightedNew)) - 1L), sum)
    expect_true(
      all(sweep(totalsNew, length(dim(totalsNew)), trialWeights[1:8]) == 0),
      info = info
    )
  }
  # equal under one seed to the sequential binomials written by hand from the
  # replayed probabilities, chain by chain
  probabilities <- predict(fit, newdata, type = "ev", combineChains = FALSE)
  set.seed(8)
  drawn <- predict(
    fit,
    newdata,
    type = "ppd",
    weights = trialWeights[1:8],
    combineChains = FALSE
  )
  set.seed(8)
  byHand <- sequentialBinomials(probabilities, trialWeights[1:8])
  expect_equal(drawn, byHand, check.attributes = FALSE, info = fitName)
  # types that draw nothing read no weight
  expect_identical(
    predict(
      fit,
      newdata,
      type = "ev",
      weights = trialWeights[1:8],
      combineChains = FALSE
    ),
    probabilities,
    info = fitName
  )
  rm(
    fit,
    newdata,
    plain,
    atNew,
    weighted,
    weightedNew,
    totals,
    totalsNew,
    probabilities,
    drawn,
    byHand
  )
}
# weights = 1 on a fit to categories is the unit vector of the category drawn
set.seed(9)
oneTrial <- predict(categoryFit, newRows, type = "ppd", weights = 1)
set.seed(9)
categories <- predict(categoryFit, newRows, type = "ppd")
expect_true(all(apply(oneTrial, 1:2, sum) == 1))
# and its categories come at the frequencies predict()'s do, over 4,000 draws
# of each row
oneTrialShares <- Reduce(
  `+`,
  lapply(seq_len(100L), function(i) {
    apply(predict(categoryFit, newRows, type = "ppd", weights = 1), 2:3, sum)
  })
) /
  4000
categoryShares <- Reduce(
  `+`,
  lapply(seq_len(100L), function(i) {
    drawn <- predict(categoryFit, newRows, type = "ppd")
    t(apply(drawn, 2L, tabulate, nbins = 3L))
  })
) /
  4000
expect_true(all(
  abs(oneTrialShares - categoryShares) < 4 * sqrt(0.25 / 4000) * sqrt(2)
))
# on a count fit, predict(newdata) and predict(newdata, weights = 1) agree
set.seed(10)
countNew <- predict(countFit, x[1:8, ], type = "ppd")
set.seed(10)
countOne <- predict(countFit, x[1:8, ], type = "ppd", weights = 1)
expect_identical(countNew, countOne)
# a 0 is no trials: K zero counts, at newdata and at a training row
expect_true(all(predict(ordinalFit, newRows, type = "ppd", weights = 0) == 0))
expect_true(all(predict(categoryFit, type = "ppd", weights = 0) == 0))
expect_true(all(predict(countFit, x[1:8, ], type = "ppd", weights = 0) == 0))
# a fraction, a negative value or Inf is refused
for (bad in c(0.5, -1, Inf)) {
  expect_error(
    predict(ordinalFit, type = "ppd", weights = bad),
    paste0(
      "the posterior predictive 'weights' of an ordinal fit are trial ",
      "counts, each draw the count per category over that many new trials, ",
      "and must be non-negative whole numbers"
    ),
    fixed = TRUE,
    info = paste("ordinal", bad)
  )
  expect_error(
    predict(categoryFit, type = "ppd", weights = bad),
    "the posterior predictive 'weights' of a multinomial fit are trial counts",
    fixed = TRUE,
    info = paste("multinomial", bad)
  )
}
# class is unchanged by weights, reading none
expect_identical(
  predict(ordinalFit, newRows, type = "class", weights = 3),
  predict(ordinalFit, newRows, type = "class")
)
rm(
  oneTrial,
  categories,
  oneTrialShares,
  categoryShares,
  countNew,
  countOne,
  bad
)

# ---- a test sample's weights mean the fit's (dec-B596) -----------------------

testRows <- df[1:6, ]
testRows$tw <- c(0, 1, 1, 0, 1, 1)
testRows$tc <- c(0, 3, 1, 0, 2, 1)
train <- df
train$tw <- mask
train$tc <- w
# logistic: a test weight of 3 is one 0/1 draw, a 0 NA
logisticTest <- quick(
  yb ~ x1 + x2,
  train,
  test = testRows,
  weights = tc,
  family = "logistic"
)
testDraws <- extract(logisticTest, "ppd", sample = "test")
expect_true(all(is.na(testDraws[, c(1L, 4L)])))
expect_true(all(testDraws[, -c(1L, 4L)] %in% c(0, 1)))
# probit and nbinom test rows of weight 0 are NA
probitTest <- quick(yb ~ x1 + x2, train, test = testRows, weights = tw)
expect_true(all(is.na(extract(probitTest, "ppd", sample = "test")[, c(
  1L,
  4L
)])))
nbinomTest <- quick(
  yc ~ x1 + x2,
  train,
  test = testRows,
  weights = tc,
  family = "nbinom"
)
nbinomTestDraws <- extract(nbinomTest, "ppd", sample = "test")
expect_true(all(is.na(nbinomTestDraws[, c(1L, 4L)])))
expect_false(anyNA(nbinomTestDraws[, -c(1L, 4L)]))
# ordinal: a 1 a category, a 0 NA
ordinalTest <- quick(
  yo ~ x1 + x2,
  train,
  test = testRows,
  weights = tw,
  family = "ordinal"
)
ordinalTestDraws <- extract(ordinalTest, "ppd", sample = "test")
expect_true(all(is.na(ordinalTestDraws[, c(1L, 4L)])))
expect_true(all(ordinalTestDraws[, -c(1L, 4L)] %in% 1:3))
# multinomial: a positive count one observation, a 0 NA
multinomialTest <- quick(
  ym ~ x1 + x2,
  train,
  test = testRows,
  weights = tc,
  family = "multinomial"
)
multinomialTestDraws <- extract(multinomialTest, "ppd", sample = "test")
expect_true(all(is.na(multinomialTestDraws[, c(1L, 4L)])))
expect_true(all(multinomialTestDraws[, -c(1L, 4L)] %in% 1:3))
# frequency gaussian: a test weight of 2.5 is one observation at sd sigma
train$tf <- w + 0.5
testRows$tf <- c(2.5, 2.5, 0, 1, 2.5, 1)
gaussianTest <- quick(
  y ~ x1 + x2,
  train,
  test = testRows,
  weights = tf,
  weights.type = "frequency"
)
gaussianTestDraws <- extract(gaussianTest, "ppd", sample = "test")
expect_true(all(is.na(gaussianTestDraws[, 3L])))
# each family's rule refuses a test value its training weights refuse
testRows$bad <- c(0.5, 1, 1, 1, 1, 1)
train$bad <- w
expect_error(
  quick(
    yb ~ x1 + x2,
    train,
    test = testRows,
    weights = bad,
    family = "logistic"
  ),
  "logistic weights are observation counts"
)
expect_error(
  quick(yc ~ x1 + x2, train, test = testRows, weights = bad, family = "nbinom"),
  "nbinom weights are observation counts"
)
# on a count family the test draws agree with predict at the test rows given
# no weights, the 0 rows aside
set.seed(11)
fromTest <- extract(nbinomTest, "ppd", sample = "test")
set.seed(11)
fromPredict <- predict(nbinomTest, testRows, type = "ppd")
expect_identical(fromTest[, -c(1L, 4L)], fromPredict[, -c(1L, 4L)])
rm(
  logisticTest,
  testDraws,
  probitTest,
  nbinomTest,
  nbinomTestDraws,
  ordinalTest,
  ordinalTestDraws,
  multinomialTest,
  multinomialTestDraws,
  gaussianTest,
  gaussianTestDraws,
  fromTest,
  fromPredict
)

# ---- a missing predict weight follows na.action, as at newdata (dec-B482) ----

missingWeights <- rep(1, n)
missingWeights[3L] <- NA
expect_error(
  predict(ordinalFit, type = "ppd", weights = missingWeights),
  "missing"
)
passed <- predict(
  ordinalFit,
  type = "ppd",
  weights = missingWeights,
  na.action = na.pass
)
expect_true(all(is.na(passed[, 3L, ])))
expect_false(anyNA(passed[, -3L, ]))
omitted <- predict(
  ordinalFit,
  type = "ppd",
  weights = missingWeights,
  na.action = na.omit
)
expect_identical(dim(omitted)[2L], n - 1L)

# ---- a predict weight of 0 is drawn at one trial and reported as none -------

# moving one row's weight from 1 to 0 leaves every other row's draws as they
# are, at new rows and at the training rows, on every count family
zeroStream <- list(
  probit = quick(yb ~ x1 + x2, df, family = "probit"),
  logistic = quick(yb ~ x1 + x2, df, family = "logistic"),
  nbinom = quick(yc ~ x1 + x2, df, family = "nbinom"),
  ordinal = quick(yo ~ x1 + x2, df, family = "ordinal"),
  multinomial = quick(ym ~ x1 + x2, df, family = "multinomial")
)
firstRow <- function(draws) {
  if (length(dim(draws)) == 2L) draws[, 1L] else draws[, 1L, ]
}
otherRows <- function(draws) {
  if (length(dim(draws)) == 2L) draws[, -1L] else draws[, -1L, , drop = FALSE]
}
for (fitName in names(zeroStream)) {
  for (sample in c("newdata", "train")) {
    drawAt <- function(first) {
      set.seed(9L)
      if (sample == "newdata") {
        predict(
          zeroStream[[fitName]],
          df[1:6, ],
          type = "ppd",
          weights = c(first, 2, 2, 2, 2, 2)
        )
      } else {
        predict(
          zeroStream[[fitName]],
          type = "ppd",
          weights = c(first, rep(2, n - 1L))
        )
      }
    }
    atZero <- drawAt(0)
    atOne <- drawAt(1)
    expect_identical(
      otherRows(atZero),
      otherRows(atOne),
      info = paste(fitName, sample, "other rows")
    )
    expect_true(
      all(firstRow(atZero) == 0),
      info = paste(fitName, sample, "a 0 draws none")
    )
  }
}

# a fit records weights.type only where a reader needs it, a gaussian fit
# under frequency
expect_identical(
  quick(y ~ x1 + x2, df, weights = w, weights.type = "frequency")$weights.type,
  "frequency"
)
expect_null(quick(y ~ x1 + x2, df, weights = w)$weights.type)
for (fitName in names(zeroStream)) {
  expect_null(zeroStream[[fitName]]$weights.type, info = fitName)
}
expect_null(
  quick(yc ~ x1 + x2, df, weights = w, family = "nbinom")$weights.type
)

# a count weight multiplies the row's log density on nbinom and multinomial
countedNbinom <- quick(
  yc ~ x1 + x2,
  df,
  weights = w,
  family = "nbinom",
  n.chains = 1L
)
mu <- extract(countedNbinom, "ev")
r <- as.vector(extract(countedNbinom, "shape"))
expect_equal(
  unname(extract(countedNbinom, "loglik")),
  matrix(
    rep(w, each = nrow(mu)) *
      dnbinom(
        rep(df$yc, each = nrow(mu)),
        size = rep(r, n),
        mu = as.vector(mu),
        log = TRUE
      ),
    nrow(mu)
  )
)
countedCategories <- quick(
  ym ~ x1 + x2,
  df,
  weights = w,
  family = "multinomial",
  n.chains = 1L
)
probs <- extract(countedCategories, "ev")
code <- match(df$ym, countedCategories$levels)
expect_equal(
  unname(extract(countedCategories, "loglik")),
  vapply(
    seq_len(n),
    function(i) w[i] * log(probs[, i, code[i]]),
    0 * probs[, 1L, 1L]
  )
)

# nbinom's counts are read as counts where they are not all 0 and 1: a linear
# leaf's covariate center and scale are the rows written out's
countLeaf <- dbarts(
  x,
  df$yc,
  weights = w,
  family = "nbinom",
  leaf.prior = dbartsPriors$linear("x1"),
  n.trees = 5L,
  control = dbartsControl(n.samples = 5L, n.chains = 1L, n.threads = 1L)
)
expect_equal(
  unname(unlist(countLeaf$data@standardization[[1L]])),
  c(mean(expanded$x1), sd(expanded$x1))
)

# the starting sigma's fallback under frequency is the sd of the rows written
# out, here at fractional counts, which leave the regression no degrees of
# freedom
few <- data.frame(y = c(0.3, 1.7, -0.4), x1 = c(0.1, 0.5, 0.9), x2 = c(1, 0, 2))
fewCounts <- c(0.5, 0.7, 0.6)
fewCenter <- sum(fewCounts * few$y) / sum(fewCounts)
expect_equal(
  dbarts(
    y ~ x1 + x2,
    few,
    weights = fewCounts,
    weights.type = "frequency",
    n.trees = 5L,
    control = dbartsControl(n.samples = 5L, n.chains = 1L, n.threads = 1L)
  )$data@sigma,
  sqrt(sum(fewCounts * (few$y - fewCenter)^2) / (sum(fewCounts) - 1))
)

# the fallback's text says what it reads: the response's scale and the
# prior's center from every row, a leaf covariate's and a basis's numbers
# from every row at its count
oneValue <- df
oneValue$y[w > 0] <- 0
oneValue$y[1:3] <- 3
fallbackIn <- c(2, 3, 1, rep(0, n - 3L))
fallbackText <- function(...) {
  texts <- character()
  withCallingHandlers(
    dbarts(
      ...,
      weights = fallbackIn,
      weights.type = "frequency",
      n.trees = 5L,
      control = dbartsControl(n.samples = 5L, n.chains = 1L, n.threads = 1L)
    ),
    dbartsScaleFallbackWarning = function(condition) {
      texts <<- c(texts, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  texts
}
expect_identical(
  fallbackText(
    y ~ x1 + x2,
    oneValue,
    leaf.prior = dbartsPriors$linear("x1")
  ),
  paste0(
    "fewer than two distinct values of the response are in the likelihood ",
    "(3 of 40 observations are in it), so the response's scale and the ",
    "prior's center are taken from all 40, and the center and scale of ",
    "\"x1\" from all 40 at their weights' counts, a weight of 0 read as 1"
  )
)
expect_identical(
  fallbackText(
    y ~ x1 + x2,
    oneValue,
    forests = list(forest(), forest(basis = ~x2))
  ),
  paste0(
    "fewer than two distinct values of the response are in the likelihood ",
    "(3 of 40 observations are in it), so the response's scale and the ",
    "prior's center are taken from all 40, and the scale of forest ",
    "\"x2\"'s basis from all 40 at their weights' counts, a weight of 0 read ",
    "as 1"
  )
)

# an nbinom sampler reads its weights as counts however it was created
# (dec-B632): a linear leaf's covariate center and scale, over counts
# {1, 0, 3}, are the expanded data's whether the sampler was created with
# them, or unweighted or with 0/1 weights and given them by setWeights
countedMoments <- function(v, counts) {
  center <- sum(counts * v) / sum(counts)
  c(center, sqrt(sum(counts * (v - center)^2) / (sum(counts) - 1)))
}
set.seed(77L)
pathX <- matrix(runif(120L), 60L, 2L, dimnames = list(NULL, c("x1", "x2")))
pathX[, 1L] <- pathX[, 1L] * (1 + 3 * (seq_len(60L) > 40L))
pathY <- as.double(rpois(60L, exp(0.3 + pathX[, 2L])))
pathCounts <- rep(c(1, 0, 3), each = 20L)
pathSampler <- function(weights = NULL) {
  suppressWarnings(dbarts::dbarts(
    pathX,
    pathY,
    family = "nbinom",
    weights = weights,
    control = dbarts::dbartsControl(
      n.chains = 1L,
      n.threads = 1L,
      n.samples = 2L,
      updateState = FALSE,
      seed = 3L
    ),
    n.trees = 4L,
    leaf.prior = dbarts::dbartsPriors$linear("x1")
  ))
}
standardizationOf <- function(sampler) {
  unname(unlist(sampler$data@standardization[[1L]]))
}
pathExpanded <- countedMoments(pathX[, 1L], pathCounts)
expect_equal(standardizationOf(pathSampler(pathCounts)), pathExpanded)
for (start in list(NULL, as.double(pathCounts > 0))) {
  given <- pathSampler(start)
  invisible(given$run(1L, 0L))
  suppressWarnings(given$setWeights(pathCounts, updateStandardization = TRUE))
  expect_equal(standardizationOf(given), pathExpanded)
}

# where a read-off falls back to every row, a count of 0 is read as 1, as the
# mask's fallback reads a masked row: an nbinom linear leaf with one row in,
# and a gaussian frequency basis whose rows in hold one response value
oneCount <- replace(rep(0, 60L), 5L, 3)
expect_equal(
  standardizationOf(pathSampler(oneCount)),
  countedMoments(pathX[, 1L], pmax(oneCount, 1))
)
# the fallback names what it reads: every number at the counts, an nbinom
# center among them
oneCountText <- character()
withCallingHandlers(
  dbarts::dbarts(
    pathX,
    pathY,
    family = "nbinom",
    weights = oneCount,
    control = dbarts::dbartsControl(n.chains = 1L, n.threads = 1L, seed = 3L),
    n.trees = 4L,
    leaf.prior = dbarts::dbartsPriors$linear("x1")
  ),
  warning = function(condition) {
    if (inherits(condition, "dbartsScaleFallbackWarning")) {
      oneCountText <<- c(oneCountText, conditionMessage(condition))
    }
    invokeRestart("muffleWarning")
  }
)
expect_identical(
  oneCountText,
  paste0(
    "fewer than two distinct values of the response are in the likelihood ",
    "(1 of 60 observations is in it), so the response's scale, the prior's ",
    "center and the center and scale of \"x1\" are taken from all 60 at ",
    "their weights' counts, a weight of 0 read as 1"
  )
)
# and at a setData whose new data hold one counted row in
newPathX <- matrix(
  runif(80L),
  40L,
  2L,
  dimnames = list(NULL, c("x1", "x2"))
)
newPathCounts <- replace(rep(0, 40L), 3L, 2)
replaced <- pathSampler(pathCounts)
invisible(replaced$run(1L, 0L))
suppressWarnings(replaced$setData(
  dbarts::dbartsData(
    newPathX,
    as.double(rpois(40L, 2)),
    weights = newPathCounts
  ),
  updateStandardization = TRUE
))
expect_equal(
  standardizationOf(replaced),
  countedMoments(newPathX[, 1L], pmax(newPathCounts, 1))
)
basisData <- data.frame(
  y = replace(rep(0, 60L), 5L, 2.5),
  x1 = pathX[, 1L],
  b = 0.5 + runif(60L)
)
basisCounts <- replace(rep(0, 60L), c(6L, 7L), c(3, 2))
basisSampler <- suppressWarnings(dbarts::dbarts(
  y ~ x1,
  basisData,
  weights = basisCounts,
  weights.type = "frequency",
  n.trees = 4L,
  forests = list(forest(), forest(basis = ~b)),
  control = dbarts::dbartsControl(n.chains = 1L, n.threads = 1L)
))
expect_equal(
  unlist(basisSampler$data@basis.scale),
  countedMoments(basisData$b, pmax(basisCounts, 1))[2L]
)
# the forests' scales, the anchor's among them, are those of the fit with
# every row in at those counts
allIn <- suppressWarnings(dbarts::dbarts(
  y ~ x1,
  basisData,
  weights = pmax(basisCounts, 1),
  weights.type = "frequency",
  n.trees = 4L,
  forests = list(forest(), forest(basis = ~b)),
  control = dbarts::dbartsControl(n.chains = 1L, n.threads = 1L)
))
for (f in 1:2) {
  expect_equal(
    basisSampler$getLeafPrior(f)$k.scale,
    allIn$getLeafPrior(f)$k.scale
  )
}

rm(
  countedMoments,
  pathX,
  pathY,
  pathCounts,
  pathSampler,
  standardizationOf,
  pathExpanded,
  start,
  given,
  oneCount,
  oneCountText,
  newPathX,
  newPathCounts,
  replaced,
  basisData,
  basisCounts,
  basisSampler,
  allIn,
  f,
  n,
  x,
  y,
  w,
  w0,
  mask,
  df,
  expanded,
  quick,
  specOf,
  survivalResponse,
  yHurdle,
  resolved,
  trainRows,
  foldExpanded,
  logisticFit,
  logisticDraws,
  trials,
  zeroWeights,
  zeroMask,
  naCases,
  name,
  ordinalFit,
  categoryFit,
  counts,
  countFit,
  trialWeights,
  newRows,
  fitName,
  sequentialBinomials,
  K,
  combine,
  info,
  margins,
  rowAxis,
  testRows,
  train,
  missingWeights,
  passed,
  omitted,
  family,
  response,
  rule,
  zeroStream,
  firstRow,
  otherRows,
  sample,
  drawAt,
  atZero,
  atOne,
  countedNbinom,
  mu,
  r,
  countedCategories,
  probs,
  code,
  countLeaf,
  few,
  fewCounts,
  fewCenter,
  oneValue,
  fallbackIn,
  fallbackText
)
