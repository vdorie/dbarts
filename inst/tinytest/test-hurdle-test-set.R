# A hurdle fit takes a test set: both parts are evaluated at its rows, extract
# and fitted read them with sample = "test", and a fit kept without its
# training fits still answers for them.

set.seed(5201L)
n <- 80L
d <- data.frame(a = runif(n), b = runif(n))
d$f <- factor(sample(c("u", "v", "w"), n, TRUE))
d$y <- ifelse(runif(n) < 0.5, 0, exp(rnorm(n)))
d$y[1L] <- 0
d$y[2L] <- 1.2
# the first test rows hold a positive response row
testRows <- 1:12
te <- d[testRows, ]
te$a <- runif(length(testRows))
m <- nrow(te)

fitArgs <- list(
  n.samples = 6L,
  n.burn = 6L,
  n.trees = 8L,
  n.chains = 2L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrees = TRUE,
  seed = 5202L
)
fitHurdle <- function(formula, ...) {
  do.call(
    bart,
    c(list(formula, family = "hurdle.lognormal"), list(...), fitArgs)
  )
}
plain <- fitHurdle(y ~ a + b + f, data = d)
withTest <- fitHurdle(y ~ a + b + f, data = d, test = te)

# --- the training channel is unchanged by a test set ---------------------------
for (type in c("ev", "prob", "bart", "loglik")) {
  expect_identical(
    extract(withTest, type = type),
    extract(plain, type = type),
    info = type
  )
}
expect_identical(fitted(withTest), fitted(plain))
expect_identical(ncol(withTest$positive$yhat.test), n + m)
expect_identical(ncol(plain$positive$yhat.test), n)

# --- extract and fitted at the test rows ----------------------------------------
# stored draws against the replayed trees, to rounding; these types draw nothing
for (type in c("prob", "bart", "ev")) {
  stored <- extract(withTest, type = type, sample = "test")
  expect_identical(dim(stored), c(2L * 6L, m), info = type)
  expect_equal(
    unname(stored),
    unname(predict(withTest, te, type = type)),
    tolerance = 1e-8,
    info = type
  )
}
set.seed(1L)
ppd <- extract(withTest, type = "ppd", sample = "test")
expect_identical(dim(ppd), c(12L, m))
expect_true(all(ppd >= 0))
expect_equal(
  fitted(withTest, type = "ev", sample = "test"),
  colMeans(extract(withTest, type = "ev", sample = "test"))
)
band <- fitted(withTest, type = "prob", sample = "test", ci.level = 0.9)
expect_identical(dim(band), c(m, 3L))
# the order of fitted's arguments: type, sample, ci.level
expect_identical(
  fitted(withTest, "prob", "test", 0.9),
  band
)
expect_identical(
  names(formals(dbarts:::fitted.bartHurdle)),
  c("object", "type", "sample", "ci.level", "...")
)
# the training rows are the default and are not padded or shifted by the test
expect_identical(length(fitted(withTest, type = "ev")), n)

# --- the refusals -----------------------------------------------------------------
expect_error(
  extract(plain, type = "ev", sample = "test"),
  "cannot extract test sample predictions"
)
expect_error(
  fitted(plain, type = "ev", sample = "test"),
  "cannot extract test sample predictions"
)
expect_error(
  extract(withTest, type = "loglik", sample = "test"),
  "no test response exists"
)
# a new missing value in a column complete in training is refused at the fit
teGap <- te
teGap$a[3L] <- NA
expect_error(
  fitHurdle(y ~ a + b + f, data = d, test = teGap),
  "test predictors have missing values"
)

# --- the matrix interface, a frame and a data object agree -------------------------
viaMatrix <- fitHurdle(
  as.matrix(d[, c("a", "b")]),
  d$y,
  test = as.matrix(te[, c("a", "b")])
)
viaFormula <- fitHurdle(y ~ a + b, data = d, test = te)
for (type in c("ev", "prob", "bart")) {
  expect_identical(
    unname(extract(viaMatrix, type = type, sample = "test")),
    unname(extract(viaFormula, type = type, sample = "test")),
    info = type
  )
}
viaObject <- fitHurdle(dbartsData(y ~ a + b + f, d, test = te))
expect_identical(
  unname(extract(viaObject, type = "ev", sample = "test")),
  unname(extract(withTest, type = "ev", sample = "test"))
)

# --- the test rows' names -----------------------------------------------------------
teNamed <- te
rownames(teNamed) <- paste0("t", seq_len(m))
named <- fitHurdle(y ~ a + b + f, data = d, test = teNamed)
expect_identical(named$row.names.test, rownames(teNamed))
expect_identical(
  colnames(extract(named, type = "ev", sample = "test")),
  rownames(teNamed)
)
expect_identical(names(fitted(named, sample = "test")), rownames(teNamed))
expect_identical(names(fitted(named)), rownames(d))

# --- keepTrainingFits = FALSE: the test rows alone ----------------------------------
dropped <- fitHurdle(
  y ~ a + b + f,
  data = d,
  test = te,
  keepTrainingFits = FALSE
)
expect_identical(ncol(dropped$positive$yhat.test), m)
expect_null(dropped$zero$yhat.train)
for (type in c("ev", "prob", "bart")) {
  expect_identical(
    extract(dropped, type = type, sample = "test"),
    extract(withTest, type = type, sample = "test"),
    info = type
  )
}
set.seed(2L)
droppedPpd <- extract(dropped, type = "ppd", sample = "test")
set.seed(2L)
expect_identical(droppedPpd, extract(withTest, type = "ppd", sample = "test"))
expect_identical(
  fitted(dropped, sample = "test"),
  fitted(withTest, sample = "test")
)
# a gaussian fit kept so refuses training reads in these words
trainingReads <- "cannot extract train sample predictions; bart must be called with 'keepTrainingFits' == TRUE"
expect_error(fitted(dropped), trainingReads, fixed = TRUE)
expect_error(residuals(dropped), trainingReads, fixed = TRUE)
expect_error(extract(dropped, type = "loglik"), trainingReads, fixed = TRUE)
expect_error(predict(dropped), trainingReads, fixed = TRUE)
expect_error(
  plot(dropped),
  "plot requires bart to be called with 'keepTrainingFits' == TRUE",
  fixed = TRUE
)
# and what does not read the training rows answers
expect_inherits(summary(dropped), "summary.bartHurdle")
expect_equal(
  unname(predict(dropped, te, type = "bart")),
  unname(extract(dropped, type = "bart", sample = "test")),
  tolerance = 1e-8
)
expect_inherits(
  pdbart(dropped, xind = "a", levs = list(0.5), pl = FALSE),
  "pdbart"
)
# without a test set the training fits are still required
expect_error(
  fitHurdle(y ~ a + b + f, data = d, keepTrainingFits = FALSE),
  "unless a test set is given"
)
expect_error(
  fitHurdle(y ~ a + b + f, data = d, test = te, keepFits = FALSE),
  "keepFits = TRUE"
)

# --- plot on a fit whose first test rows hold a positive row --------------------------
expect_true(any(d$y[testRows] > 0))
grDevices::pdf(NULL)
expect_silent(plot(withTest))
grDevices::dev.off()

# --- the test channel of a stubbed fit is absent, not an error of its own -------------
expect_error(
  extract(
    structure(
      list(
        y = numeric(2L),
        zero = structure(list(yhat.train = matrix(0, 3L, 2L)), class = "bart"),
        positive = structure(
          list(yhat.test = matrix(0, 3L, 2L)),
          class = "bart"
        )
      ),
      class = "bartHurdle"
    ),
    type = "prob",
    sample = "test"
  ),
  "cannot extract test sample predictions"
)
