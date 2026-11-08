# Offsets on a hurdle fit. A plain offset, the argument or an offset() term,
# enters the positive part, whose log y it shifts; offset = list(zero = ,
# positive = ) states each part's. predict re-evaluates each part's recorded
# offset on newdata, and replaces it at the training rows from the stored draws.

set.seed(5401L)
n <- 90L
d <- data.frame(a = runif(n), b = runif(n))
d$t <- runif(n, 1, 3)
d$u <- runif(n, -0.5, 0.5)
d$y <- ifelse(runif(n) < 0.5, 0, exp(rnorm(n)))
d$y[1L] <- 0
d$y[2L] <- 1.2

seed <- 5402L
fitArgs <- list(
  n.samples = 6L,
  n.burn = 6L,
  n.trees = 8L,
  n.chains = 2L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrees = TRUE
)
fitHurdle <- function(formula, ...) {
  fitCall <- match.call()
  fitCall[[1L]] <- quote(bart)
  fitCall$family <- "hurdle.lognormal"
  fitCall[names(fitArgs)] <- fitArgs
  fitCall$seed <- seed
  eval(fitCall, parent.frame())
}
channels <- function(fit, types = c("prob", "bart")) {
  lapply(setNames(types, types), function(type) extract(fit, type = type))
}
partSeeds <- dbarts:::withFixedSeed(seed, sample.int(.Machine$integer.max, 2L))

plain <- fitHurdle(y ~ a + b, data = d)
xm0 <- as.matrix(d[, c("a", "b")])
positiveOnly <- fitHurdle(y ~ a + b, data = d, offset = log(t))
zeroOnly <- fitHurdle(y ~ a + b, data = d, offset = list(zero = u))
both <- fitHurdle(
  y ~ a + b,
  data = d,
  offset = list(zero = u, positive = log(t))
)

# --- each type moves with the part whose offset it carries --------------------------
expect_identical(extract(positiveOnly, "prob"), extract(plain, "prob"))
expect_false(identical(extract(positiveOnly, "bart"), extract(plain, "bart")))
expect_identical(extract(zeroOnly, "bart"), extract(plain, "bart"))
expect_false(identical(extract(zeroOnly, "prob"), extract(plain, "prob")))
expect_identical(extract(both, "prob"), extract(zeroOnly, "prob"))
expect_identical(extract(both, "bart"), extract(positiveOnly, "bart"))
# the zero part has no offset record of a plain offset
expect_null(positiveOnly$zero$offset)
expect_identical(positiveOnly$positive$offset, log(d$t)[d$y > 0])
expect_identical(both$zero$offset, d$u)

# --- a plain offset, a named one and an offset() term enter the positive part -----------
expect_identical(
  channels(fitHurdle(y ~ a + b, data = d, offset = list(positive = log(t)))),
  channels(positiveOnly)
)
expect_identical(
  channels(fitHurdle(y ~ a + b + offset(log(t)), data = d)),
  channels(positiveOnly)
)
# the argument and a term are summed
summed <- fitHurdle(y ~ a + b + offset(u), data = d, offset = log(t) - u)
expect_equal(summed$positive$offset, log(d$t)[d$y > 0])
# a vector and an expression are the same offset
expect_identical(
  lapply(channels(fitHurdle(xm0, d$y, offset = log(d$t))), unname),
  lapply(channels(positiveOnly), unname)
)

# --- each part is the standalone fit with its offset ---------------------------------------
zeroAlone <- do.call(
  bart,
  c(
    list(as.double(y > 0) ~ a + b, d, family = "probit", offset = quote(u)),
    fitArgs,
    seed = partSeeds[1L]
  )
)
expect_identical(both$zero$yhat.train, zeroAlone$yhat.train)
positiveRows <- d$y > 0
positiveAlone <- do.call(
  bart,
  c(
    list(
      log(y) ~ a + b,
      d[positiveRows, ],
      family = "gaussian",
      offset = quote(log(t)),
      test = d
    ),
    fitArgs,
    seed = partSeeds[2L]
  )
)
expect_identical(both$positive$yhat.train, positiveAlone$yhat.train)
expect_identical(both$positive$yhat.test, positiveAlone$yhat.test)

# --- predict at newdata ------------------------------------------------------------------------
# the training rows as newdata give the stored draws, offsets re-evaluated
for (fit in list(positiveOnly, both)) {
  for (type in c("ev", "prob", "bart")) {
    expect_equal(
      unname(predict(fit, d, type = type)),
      unname(extract(fit, type = type)),
      tolerance = 1e-8,
      info = type
    )
  }
}
withTerm <- fitHurdle(y ~ a + b + offset(log(t)), data = d)
expect_equal(
  unname(predict(withTerm, d, type = "ev")),
  unname(extract(withTerm, type = "ev")),
  tolerance = 1e-8
)
# an offset given adds to the recorded one, on its part (dec-B184)
extra <- runif(n)
shifted <- predict(both, d, type = "bart", offset = list(positive = extra))
expect_equal(shifted, predict(both, d, type = "bart") + rep(extra, each = 12L))
expect_identical(
  predict(both, d, type = "prob", offset = list(positive = extra)),
  predict(both, d, type = "prob")
)
# a plain offset given reaches the positive part
expect_identical(
  predict(both, d, type = "bart", offset = extra),
  shifted
)
# a formula's offset() term on the zero part is not there: it is the positive's
expect_equal(
  unname(predict(withTerm$zero, d, type = "ev")),
  unname(predict(withTerm, d, type = "prob")),
  tolerance = 1e-8
)
changed <- d
changed$t <- d$t * 3
expect_identical(
  predict(withTerm$zero, changed, type = "ev"),
  predict(withTerm$zero, d, type = "ev")
)
expect_null(attr(attr(withTerm$zero$fit$data@x, "terms"), "offset"))
expect_false(is.null(attr(
  attr(withTerm$positive$fit$data@x, "terms"),
  "offset"
)))
# a plain vector the fit was given cannot be re-evaluated on newdata, naming the part
vectorFit <- fitHurdle(
  as.matrix(d[, c("a", "b")]),
  d$y,
  offset = list(zero = d$u, positive = log(d$t))
)
expect_error(
  predict(vectorFit, as.matrix(d[1:5, c("a", "b")])),
  "zero part's offset"
)
expect_identical(
  dim(predict(
    vectorFit,
    as.matrix(d[1:5, c("a", "b")]),
    offset = list(zero = 0, positive = 0)
  )),
  c(12L, 5L)
)
expect_error(
  predict(both, d[1:5, ], offset = list(other = 1)),
  "named 'zero' and 'positive'"
)

# --- predict with no newdata -------------------------------------------------------------------
expect_equal(
  unname(predict(both, type = "ev")),
  unname(extract(both, type = "ev"))
)
# an offset DIFFERENT from the fit's, read off the stored draws, equals predict at
# the training rows as newdata with that offset, to rounding
o2 <- runif(n)
o3 <- runif(n, -1, 1)
for (type in c("bart", "prob", "ev")) {
  viaStored <- predict(both, type = type, offset = list(positive = o2))
  viaNewdata <- predict(
    both,
    d,
    type = type,
    offset = list(positive = o2 - log(d$t))
  )
  expect_equal(
    unname(viaStored),
    unname(viaNewdata),
    tolerance = 1e-8,
    info = type
  )
  viaStored <- predict(both, type = type, offset = list(zero = o3))
  viaNewdata <- predict(both, d, type = type, offset = list(zero = o3 - d$u))
  expect_equal(
    unname(viaStored),
    unname(viaNewdata),
    tolerance = 1e-8,
    info = type
  )
}
# and the other part's draws are unmoved
expect_identical(
  predict(both, type = "prob", offset = list(positive = o2)),
  predict(both, type = "prob")
)
expect_identical(
  predict(both, type = "bart", offset = list(zero = o3)),
  predict(both, type = "bart")
)
# one value, or one per training row
expect_equal(
  predict(both, type = "bart", offset = list(positive = 0.5)),
  predict(both, type = "bart", offset = list(positive = rep(0.5, n)))
)
expect_error(
  predict(both, type = "bart", offset = list(positive = 1:3)),
  "one number or one per training row"
)
# a missing value follows na.action (dec-B482)
withGap <- replace(o2, 4L, NA)
expect_identical(
  dim(predict(both, type = "ev", offset = withGap, na.action = na.omit)),
  c(12L, n - 1L)
)
expect_identical(
  dim(predict(both, type = "ev", offset = withGap, na.action = na.pass)),
  c(12L, n)
)
expect_error(predict(both, type = "ev", offset = withGap), "missing values")
expect_identical(
  dim(predict(
    both,
    d[1:5, ],
    offset = list(positive = c(1, NA, 1, 1, 1)),
    na.action = na.omit
  )),
  c(12L, 4L)
)

# --- offset.test ----------------------------------------------------------------------------------
te <- d[1:10, ]
te$a <- runif(10)
te$t <- runif(10, 1, 3)
te$u <- runif(10, -0.5, 0.5)
withTest <- fitHurdle(
  y ~ a + b,
  data = d,
  test = te,
  offset = list(zero = u, positive = log(t))
)
# each element is evaluated on the test frame
expect_equal(withTest$zero$offset.test, te$u)
expect_equal(withTest$positive$offset.test, c(log(d$t), log(te$t)))
expect_equal(
  unname(extract(withTest, type = "ev", sample = "test")),
  unname(predict(withTest, te, type = "ev")),
  tolerance = 1e-8
)
# the training channel is unchanged by a test set
expect_identical(
  channels(withTest, c("prob", "bart")),
  channels(both, c("prob", "bart"))
)
# stated apart
stated <- fitHurdle(
  y ~ a + b,
  data = d,
  test = te,
  offset = list(zero = u, positive = log(t)),
  offset.test = list(zero = rep(0, 10), positive = rep(1, 10))
)
expect_equal(stated$zero$offset.test, rep(0, 10))
expect_equal(stated$positive$offset.test[91:100], rep(1, 10))
# a matrix fit's vector offset is reused for a test set of n rows and refused for
# another number
xm <- as.matrix(d[, c("a", "b")])
expect_error(
  fitHurdle(xm, d$y, test = xm[1:5, ], offset = list(zero = d$u)),
  "offset.test"
)
reused <- fitHurdle(
  xm,
  d$y,
  test = xm,
  offset = list(zero = d$u, positive = d$t)
)
expect_equal(reused$zero$offset.test, d$u)

# --- a missing offset follows na.action on both routes (dec-A89) ----------------------------------
dNA <- d
dNA$u[5L] <- NA
dNA$t[9L] <- NA
for (action in list(na.omit, na.exclude)) {
  viaFormula <- fitHurdle(
    y ~ a + b,
    data = dNA,
    offset = list(zero = u, positive = log(t)),
    na.action = action
  )
  viaMatrix <- fitHurdle(
    dNA[, c("a", "b")],
    dNA$y,
    offset = list(zero = dNA$u, positive = log(dNA$t)),
    na.action = action
  )
  expect_identical(viaMatrix$na.action, viaFormula$na.action)
  expect_identical(
    unname(extract(viaMatrix, "ev")),
    unname(extract(viaFormula, "ev"))
  )
  expect_identical(length(viaFormula$y), n - 2L)
  expect_identical(length(viaFormula$zero$offset), n - 2L)
}
expect_error(
  fitHurdle(y ~ a + b, data = dNA, offset = list(zero = u)),
  "missing|NA"
)
expect_error(
  fitHurdle(dNA[, c("a", "b")], dNA$y, offset = list(zero = dNA$u)),
  "'offset' contains missing values"
)

# --- a data object has one offset slot, and what is written beside it is ignored ------------------------
expect_warning(
  fitHurdle(
    dbartsData(y ~ a + b, d, keep.data = FALSE),
    offset = list(zero = u)
  ),
  "remaining arguments are ignored"
)

# --- partial dependence carries each part's offset -------------------------------------------------------------
manual <- function(fit, type, level) {
  rows <- d
  rows$a <- level
  rowMeans(predict(fit, rows, type = type))
}
for (type in c("prob", "bart", "ev")) {
  pd <- pdbart(
    both,
    xind = "a",
    levs = list(c(0.3, 0.6)),
    type = type,
    pl = FALSE
  )
  expect_equal(
    unname(pd$fd[[1L]][, 1L]),
    unname(manual(both, type, 0.3)),
    tolerance = 1e-8,
    info = type
  )
}
# an offset() term on the positive part leaves the "prob" average where it was
pdTerm <- pdbart(
  withTerm,
  xind = "a",
  levs = list(0.3),
  type = "prob",
  pl = FALSE
)
expect_equal(
  unname(pdTerm$fd[[1L]][, 1L]),
  unname(manual(withTerm, "prob", 0.3)),
  tolerance = 1e-8
)
# a vector offset on rows the average leaves out is each remaining row's own
keepRows <- rep(c(1, 0), length.out = n)
pdVector <- pdbart(
  vectorFit,
  xind = 1L,
  levs = list(0.3),
  type = "ev",
  average.weights = keepRows,
  pl = FALSE
)
index <- which(keepRows > 0)
rows <- xm[index, , drop = FALSE]
rows[, 1L] <- 0.3
expect_equal(
  unname(pdVector$fd[[1L]][, 1L]),
  unname(rowMeans(predict(
    vectorFit,
    rows,
    type = "ev",
    offset = list(zero = d$u[index], positive = log(d$t)[index])
  ))),
  tolerance = 1e-8
)
expect_error(
  pdbart(
    vectorFit,
    xind = 1L,
    levs = list(0.3),
    newdata = xm[1:5, ],
    pl = FALSE
  ),
  "zero part's 'offset'"
)

# --- plain-vector offsets under keepTrainingFits = FALSE ---------------------------------------------
# the fit keeps the positive part's offset at every training row beside y, so
# partial dependence reads it with or without the training channels, and with
# rows cut by a subset
xTest <- xm[11:40, ]
offsetTest <- list(zero = d$u[11:40], positive = log(d$t)[11:40])
keepRows <- rep(c(1, 0), length.out = n)
averaged <- function(fit) {
  unname(pdbart(
    fit,
    xind = 1L,
    levs = list(0.3),
    type = "ev",
    average.weights = rep(c(1, 0), length.out = length(fit$y)),
    pl = FALSE
  )$fd[[1L]][, 1L])
}
for (cut in list(NULL, quote(seq_len(n) > 3L))) {
  byTraining <- lapply(c(TRUE, FALSE), function(trainingFits) {
    fitHurdle(
      xm,
      d$y,
      subset = eval(cut),
      test = xTest,
      offset = list(zero = d$u, positive = log(d$t)),
      offset.test = offsetTest,
      keepTrainingFits = trainingFits
    )
  })
  expect_identical(byTraining[[2L]]$offset, byTraining[[1L]]$offset)
  expect_equal(
    averaged(byTraining[[2L]]),
    averaged(byTraining[[1L]]),
    tolerance = 1e-8
  )
  rows <- xm[if (is.null(cut)) seq_len(n) else which(eval(cut)), , drop = FALSE]
  index <- which(rep(c(1, 0), length.out = nrow(rows)) > 0)
  rows <- rows[index, , drop = FALSE]
  rows[, 1L] <- 0.3
  use <- if (is.null(cut)) seq_len(n) else which(eval(cut))
  expect_equal(
    averaged(byTraining[[2L]]),
    unname(rowMeans(predict(
      byTraining[[2L]],
      rows,
      type = "ev",
      offset = list(zero = d$u[use][index], positive = log(d$t)[use][index])
    ))),
    tolerance = 1e-8
  )
}
expect_true(is.list(pd2bart(
  byTraining[[2L]],
  xind = 1:2,
  levs = list(0.3, 0.4),
  pl = FALSE
)))

# --- a zero-part offset on a data frame whose row names repeat ---------------------------------
# the second reading of the offset must not match rows by name
dDup <- d
dDup$u[20L] <- NA
dDup$a[c(4L, 10L)] <- NA
dUnique <- dDup
attr(dDup, "row.names") <- rep(c("r1", "r2", "r3"), length.out = n)
for (action in list(na.omit, na.exclude)) {
  viaUnique <- fitHurdle(
    y ~ a + b,
    data = dUnique,
    offset = list(zero = u),
    na.action = action
  )
  viaDup <- fitHurdle(
    y ~ a + b,
    data = dDup,
    offset = list(zero = u),
    na.action = action
  )
  expect_identical(unname(viaDup$na.action), unname(viaUnique$na.action))
  expect_identical(
    lapply(channels(viaDup), unname),
    lapply(channels(viaUnique), unname)
  )
}
