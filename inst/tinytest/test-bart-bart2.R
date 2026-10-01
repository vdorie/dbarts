source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# the legacy door and the modern door fit the same data to the same place
n.burn <- 200L
n.sims <- 400L
bartFit <- dbarts::bartBT(
  testData$x,
  testData$y,
  ndpost = n.sims,
  nskip = n.burn,
  ntree = 50L,
  nchain = 4L,
  nthread = 1L,
  verbose = FALSE
)

bart2Fit <- dbarts::bart(
  testData$x,
  testData$y,
  n.samples = n.sims,
  n.burn = n.burn,
  n.trees = 50L,
  n.chains = 4L,
  n.threads = 1L,
  keepTrees = FALSE,
  verbose = FALSE
)

expect_inherits(bartFit, "bart")
expect_inherits(bart2Fit, "bart")
expect_true(
  sqrt(mean((bartFit$yhat.train.mean - bart2Fit$yhat.train.mean)^2)) /
    sd(testData$y) <
    0.1
)

rm(bart2Fit, bartFit, n.sims, n.burn)

rm(testData)

# The legacy door is 0.9-34's argument list exactly: the modern settings
# never on CRAN's 'bart' are gone from it, so each is an ordinary unused
# argument rather than a silently honoured extra.

set.seed(202)
nS10 <- 30L
xS10 <- matrix(rnorm(nS10 * 2L), nS10, dimnames = list(NULL, c("a", "b")))
yS10 <- rnorm(nS10)
quickS10 <- list(
  ndpost = 5L,
  nskip = 2L,
  ntree = 3L,
  nchain = 1L,
  nthread = 1L,
  verbose = FALSE
)

for (extra in list(
  list(subset = 1:10),
  list(storage = "single"),
  list(family = "logistic")
)) {
  expect_error(
    do.call(dbarts::bartBT, c(list(xS10, yS10), extra, quickS10)),
    pattern = "unused argument"
  )
}

# every one of the three is still reachable at the modern door, which is
# where the message that names it points
expect_inherits(
  dbarts::bart(
    xS10,
    yS10,
    subset = 1:20,
    storage = "single",
    n.samples = 5L,
    n.burn = 2L,
    n.trees = 3L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  ),
  "bart"
)

# a factor response of three or more levels was fit as its integer level
# codes by 0.9-x; refused here, naming both remedies
y3S10 <- factor(sample(c("a", "b", "c"), nS10, replace = TRUE))
factorMsg <- tryCatch(
  do.call(dbarts::bartBT, c(list(xS10, y3S10), quickS10)),
  error = function(e) conditionMessage(e)
)
expect_true(grepl("family = \"multinomial\"", factorMsg, fixed = TRUE))
expect_true(grepl("as.integer(y) - 1L", factorMsg, fixed = TRUE))
# both 3+-level factor kinds reach the SAME refusal through the formula
# path too: the ordered one after dbarts() resolves it to ordinal, the
# unordered one from inside dbarts()'s own categorical refusal, whose
# message names family tokens this door has no formal for
dfOrderedS10 <- data.frame(a = xS10[, 1L], b = xS10[, 2L], y = ordered(y3S10))
dfUnorderedS10 <- data.frame(a = xS10[, 1L], b = xS10[, 2L], y = y3S10)
for (dfCase in list(dfOrderedS10, dfUnorderedS10)) {
  formulaMsg <- tryCatch(
    do.call(dbarts::bartBT, c(list(y ~ a + b, dfCase), quickS10)),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl("three or more levels", formulaMsg, fixed = TRUE))
  expect_true(grepl("as.integer(y) - 1L", formulaMsg, fixed = TRUE))
  expect_false(grepl("bart2", formulaMsg, fixed = TRUE))
}
rm(y3S10, factorMsg, dfOrderedS10, dfUnorderedS10, formulaMsg, dfCase)

# a keepSampler fit carries $n.chains too, not only when the sampler itself
# is dropped
keptBartFit <- do.call(
  dbarts::bartBT,
  c(list(xS10, yS10, keepsampler = TRUE), quickS10)
)
expect_equal(keptBartFit$n.chains, 1L)
keptBart2Fit <- dbarts::bart(
  xS10,
  yS10,
  n.samples = 5L,
  n.burn = 2L,
  n.trees = 3L,
  n.chains = 2L,
  n.threads = 1L,
  keepSampler = TRUE,
  verbose = FALSE
)
expect_equal(keptBart2Fit$n.chains, 2L)

# partial matching still reaches every legacy formal
abbrevFit1 <- dbarts::bartBT(
  xS10,
  yS10,
  ntre = 3L,
  ndpost = 5L,
  nskip = 2L,
  nchain = 1L,
  nthread = 1L,
  verbose = FALSE,
  seed = 55L
)
expect_inherits(abbrevFit1, "bart")
abbrevFit2 <- dbarts::bartBT(
  xS10,
  yS10,
  ntree = 3L,
  ndpo = 5L,
  nskip = 2L,
  nchain = 1L,
  nthread = 1L,
  verbose = FALSE,
  seed = 55L
)
expect_inherits(abbrevFit2, "bart")

# the refusal names the family that fits the response: multinomial for an
# unordered factor, ordinal for an ordered one
y3Fam <- factor(sample(c("a", "b", "c"), nS10, replace = TRUE))
famMsg <- function(y) {
  tryCatch(
    do.call(dbarts::bartBT, c(list(xS10, y), quickS10)),
    error = function(e) conditionMessage(e)
  )
}
famMsgFormula <- function(y) {
  df <- data.frame(a = xS10[, 1L], b = xS10[, 2L], y = y)
  tryCatch(
    do.call(dbarts::bartBT, c(list(y ~ a + b, df), quickS10)),
    error = function(e) conditionMessage(e)
  )
}
for (msgFn in list(famMsg, famMsgFormula)) {
  expect_true(grepl("family = \"multinomial\"", msgFn(y3Fam), fixed = TRUE))
  expect_false(grepl("ordinal", msgFn(y3Fam), fixed = TRUE))
  expect_true(
    grepl("family = \"ordinal\"", msgFn(ordered(y3Fam)), fixed = TRUE)
  )
  expect_false(grepl("multinomial", msgFn(ordered(y3Fam)), fixed = TRUE))
}
rm(y3Fam, famMsg, famMsgFormula, msgFn)

rm(nS10, xS10, yS10, quickS10, abbrevFit1, abbrevFit2)

# keeptrees = TRUE, keepsampler = FALSE keeps $fit anyway (keepsampler's
# default IS keeptrees; an explicit FALSE override must not lose it) -
# predict and extract("trees") both work
set.seed(303)
nKT <- 40L
xKT <- matrix(rnorm(nKT * 2L), nKT, 2L)
yKT <- xKT[, 1L] + rnorm(nKT)
fitKT <- dbarts::bartBT(
  xKT,
  yKT,
  ndpost = 5L,
  nskip = 4L,
  ntree = 4L,
  nchain = 1L,
  nthread = 1L,
  verbose = FALSE,
  seed = 66L,
  keeptrees = TRUE,
  keepsampler = FALSE
)
expect_false(is.null(fitKT$fit))
expect_silent(predict(fitKT, xKT[1:3, ]))
expect_silent(extract(fitKT, "trees"))

rm(nKT, xKT, yKT, fitKT)

# bartBT()'s x/y route (factors = "indicators") keeps the training level
# table and recodes a test factor against it by label, as lm does: a test
# factor declaring fewer levels, or more unused ones, predicts as the
# training-level one; a value with no training rows is refused by name
set.seed(404)
nF <- 60L
dF <- data.frame(
  x1 = rnorm(nF),
  f = factor(sample(c("a", "b", "c"), nF, TRUE))
)
dF$y <- dF$x1 *
  2 +
  ifelse(dF$f == "a", 3, ifelse(dF$f == "b", -3, 0)) +
  rnorm(nF, 0, 0.2)
trF <- dF[1:40, c("x1", "f")]
ytrF <- dF$y[1:40]
teSame <- dF[41:60, c("x1", "f")]
teSame <- teSame[teSame$f != "c", ]
teFewer <- teSame
teFewer$f <- droplevels(teFewer$f)
teExtra <- teSame
teExtra$f <- factor(teExtra$f, levels = c("d", "b", "a", "c"))
teMany <- teSame
teMany$f <- factor(teMany$f, levels = c("a", "b", "c", paste0("z", 1:5000)))
fitBT <- function(test) {
  dbarts::bartBT(
    trF,
    ytrF,
    test,
    ndpost = 5L,
    nskip = 4L,
    ntree = 5L,
    nchain = 1L,
    verbose = FALSE,
    seed = 4L
  )$yhat.test
}
expected <- fitBT(teSame)
expect_identical(fitBT(teFewer), expected)
expect_identical(fitBT(teExtra), expected)
expect_identical(fitBT(teMany), expected)
# predict() reaches the same funnel with the fit's stored level table
fitF <- dbarts::bartBT(
  trF,
  ytrF,
  ndpost = 5L,
  nskip = 4L,
  ntree = 5L,
  nchain = 1L,
  verbose = FALSE,
  seed = 4L,
  keeptrees = TRUE
)
expected <- predict(fitF, teSame)
expect_identical(predict(fitF, teFewer), expected)
expect_identical(predict(fitF, teExtra), expected)
expect_identical(predict(fitF, teMany), expected)
# a character test column is recoded the same way, and a value training
# never saw is refused by name
teChar <- teSame
teChar$f <- as.character(teChar$f)
expect_identical(predict(fitF, teChar), expected)
teChar$f[1L] <- "d"
expect_error(
  predict(fitF, teChar),
  pattern = "test data factor 'f' has level 'd' with no training rows"
)
rm(nF, dF, trF, ytrF, teSame, teFewer, teExtra, teMany, teChar, fitF)
rm(fitBT, expected)

# bart(factors = "indicators") keeps the table too: a level subset away from
# training, or declared and never observed, has no indicator column and is
# refused by name rather than predicted as another level
set.seed(405)
dI <- data.frame(
  x = runif(90L),
  g = factor(
    sample(c("a", "b", "c"), 90L, TRUE),
    levels = c("a", "b", "c", "z")
  )
)
dI$y <- dI$x + 2 * (dI$g == "b") - 2 * (dI$g == "c") + rnorm(90L, 0, 0.1)
fitIndicators <- function(rows = TRUE) {
  bart(
    y ~ x + g,
    data = dI,
    subset = rows,
    factors = "indicators",
    n.trees = 10L,
    n.samples = 10L,
    n.burn = 10L,
    n.chains = 1L,
    n.threads = 1L,
    keepTrees = TRUE,
    verbose = FALSE,
    seed = 2L
  )
}
fitSubset <- fitIndicators(dI$g != "c")
newI <- data.frame(x = 0.5, g = factor(c("a", "b"), levels = levels(dI$g)))
expect_error(
  predict(fitSubset, rbind(newI, data.frame(x = 0.5, g = "c"))),
  pattern = "test data factor 'g' has level 'c' with no training rows"
)
fitAll <- fitIndicators()
expect_error(
  predict(fitAll, data.frame(x = 0.5, g = factor("z", levels = levels(dI$g)))),
  pattern = "test data factor 'g' has level 'z' with no training rows"
)
# a test factor declaring only some levels predicts as the full-level one
expect_identical(
  predict(fitAll, data.frame(x = 0.5, g = factor("b"))),
  predict(fitAll, data.frame(x = 0.5, g = factor("b", levels = levels(dI$g))))
)
rm(dI, fitIndicators, fitSubset, newI, fitAll)

# bartBT() keeps 0.9-34's row rule: an incomplete row is dropped, not
# modelled and not refused, and the door takes no na.action of its own
set.seed(505)
nMiss <- 40L
xMiss <- matrix(rnorm(nMiss * 2L), nMiss, 2L)
xMiss[1L, 1L] <- NA_real_
yMiss <- rnorm(nMiss)
expect_equal(
  length(
    dbarts::bartBT(
      xMiss,
      yMiss,
      ndpost = 3L,
      nskip = 2L,
      ntree = 3L,
      verbose = FALSE
    )$y
  ),
  nMiss - 1L
)
expect_error(
  dbarts::bartBT(
    xMiss,
    yMiss,
    ndpost = 3L,
    nskip = 2L,
    ntree = 3L,
    verbose = FALSE,
    na.action = na.pass
  ),
  pattern = "unused argument"
)
rm(nMiss, xMiss, yMiss)
