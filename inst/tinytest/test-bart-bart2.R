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

# The legacy door is 0.9-34's argument list exactly: the five settings it
# briefly carried on this branch and never on CRAN are gone, so each is an
# ordinary unused argument rather than a silently honoured extra.

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
  list(family = "logistic"),
  list(resid.dist = quote(gaussian)),
  list(prior.scale = 2.0)
)) {
  expect_error(
    do.call(dbarts::bartBT, c(list(xS10, yS10), extra, quickS10)),
    pattern = "unused argument"
  )
}

# every one of the five is still reachable at the modern door, which is
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
# an ORDERED 3+-level response reaches the same refusal through the
# formula path, where dbarts() resolves it to ordinal first
dfS10 <- data.frame(a = xS10[, 1L], b = xS10[, 2L], y = ordered(y3S10))
expect_error(
  do.call(dbarts::bartBT, c(list(y ~ a + b, dfS10), quickS10)),
  pattern = "three or more levels"
)
rm(y3S10, factorMsg, dfS10)

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

# bartBT()'s x/y route (factors = "indicators", no stored level table) names
# the mismatching factor and its levels rather than blaming the column
# count when a test factor's levels differ from training's
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
teFewer <- dF[41:60, c("x1", "f")]
teFewer <- teFewer[teFewer$f != "c", ]
teFewer$f <- droplevels(teFewer$f)
expect_error(
  dbarts::bartBT(
    trF,
    ytrF,
    teFewer,
    ndpost = 5L,
    nskip = 4L,
    ntree = 5L,
    nchain = 1L,
    verbose = FALSE,
    seed = 4L
  ),
  pattern = "'test' factor 'f' does not match training's indicator columns"
)
# a test factor declaring MORE levels than training's is named by the level
# tables themselves, before the drop-pattern replay indexes past its end;
# the count of undeclared levels is unbounded and does not change the refusal
teExtra <- dF[41:60, c("x1", "f")]
teExtra$f <- factor(teExtra$f, levels = c("a", "b", "c", "d"))
teMany <- dF[41:60, c("x1", "f")]
teMany$f <- factor(teMany$f, levels = c("a", "b", "c", paste0("z", 1:5000)))
expect_error(
  dbarts::bartBT(
    trF,
    ytrF,
    teExtra,
    ndpost = 5L,
    nskip = 4L,
    ntree = 5L,
    nchain = 1L,
    verbose = FALSE,
    seed = 4L
  ),
  pattern = "'test' factor 'f' declares 4 levels but the training design declared 3"
)
expect_error(
  dbarts::bartBT(
    trF,
    ytrF,
    teMany,
    ndpost = 5L,
    nskip = 4L,
    ntree = 5L,
    nchain = 1L,
    verbose = FALSE,
    seed = 4L
  ),
  pattern = "'test' factor 'f' declares 5003 levels"
)
# predict() reaches the same funnel with the fit's stored drop pattern
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
expect_error(
  predict(fitF, teExtra),
  pattern = "'test' factor 'f' declares 4 levels but the training design declared 3"
)
expect_error(
  predict(fitF, teMany),
  pattern = "'test' factor 'f' declares 5003 levels"
)
expect_error(
  predict(fitF, teFewer),
  pattern = "'test' factor 'f' does not match training's indicator columns"
)
# a character test column is expanded as a factor, so its distinct values
# are compared to the training table the same way
teChar <- dF[41:60, c("x1", "f")]
teChar$f <- as.character(teChar$f)
teChar$f[1L] <- "d"
expect_error(
  predict(fitF, teChar),
  pattern = "'test' factor 'f' declares 4 levels but the training design declared 3"
)
rm(nF, dF, trF, ytrF, teFewer, teExtra, teMany, teChar, fitF)

# bartBT()'s missing-predictor refusal cannot advise an argument bartBT()
# itself rejects (missing = "incorporate", bart()/dbarts() only); it
# points at those front doors instead
set.seed(505)
nMiss <- 40L
xMiss <- matrix(rnorm(nMiss * 2L), nMiss, 2L)
xMiss[1L, 1L] <- NA_real_
yMiss <- rnorm(nMiss)
expect_error(
  dbarts::bartBT(
    xMiss,
    yMiss,
    ndpost = 3L,
    nskip = 2L,
    ntree = 3L,
    verbose = FALSE
  ),
  pattern = "use bart\\(\\) or dbarts\\(\\)"
)
expect_error(
  dbarts::bartBT(
    xMiss,
    yMiss,
    ndpost = 3L,
    nskip = 2L,
    ntree = 3L,
    verbose = FALSE,
    missing = "incorporate"
  ),
  pattern = "unused argument"
)
rm(nMiss, xMiss, yMiss)
