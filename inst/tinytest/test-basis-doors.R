# Every door that builds a model reads a forest's basis against the data by
# one reader: a name is read from the data first (a frame's column, a matrix's
# column, a column a block-holding object holds, a held basis by the text it
# was written as) and then where forest() was called; a value has one entry per
# row handed over and is cut with the rows; a forests list states the forests
# in full; a single string names its column; a hazard fit's basis has one
# value per subject.

suppressMessages(library(survival))

forest <- dbartsForests$forest

set.seed(20261108L)
n <- 60L
d <- data.frame(
  x1 = rnorm(n),
  x2 = runif(n),
  z = rbinom(n, 1L, 0.5),
  dose = runif(n, 0.5, 2)
)
d$y <- d$x1 + d$z * d$x2 + d$dose * d$x1 + rnorm(n, 0, 0.3)
X <- as.matrix(d[c("x1", "x2", "dose")])

ctl <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.samples = 5L,
  n.burn = 3L,
  updateState = FALSE
)
draws <- function(sampler) sampler$run(3L, 5L)$train
# the call is made as written, so that 'subset' and 'weights' are read where a
# caller's would be
fitDraws <- function(...) {
  call <- match.call()
  call[[1L]] <- quote(dbarts)
  call[c("n.trees", "control", "seed")] <- list(5L, quote(ctl), 3L)
  draws(eval(call, parent.frame()))
}
specDraws <- function(data, ...) {
  spec <- dbartsSpec(data, ctl, ..., n.trees = 5L, seed = 3L)
  draws(new("dbartsSampler", spec$control, spec$model, spec$data))
}
refusal <- function(expr) {
  tryCatch(
    {
      expr
      "created"
    },
    error = conditionMessage
  )
}

# ---- D1: beside a frame-keeping object a name reads the frame first ----
# a caller's variable of the column's name, unlike it
dose <- rev(d$dose)
frameObject <- dbartsData(y ~ x1 + x2, d)
byFormula <- fitDraws(
  y ~ x1 + x2,
  d,
  forests = list(forest(), forest(basis = dose))
)
expect_identical(
  fitDraws(frameObject, forests = list(forest(), forest(basis = dose))),
  byFormula
)
expect_identical(
  specDraws(frameObject, forests = list(forest(), forest(basis = dose))),
  byFormula
)
byCaller <- fitDraws(
  y ~ x1 + x2,
  d,
  forests = list(forest(), do.call(forest, list(basis = dose)))
)
expect_false(identical(byCaller, byFormula))

# ---- D2: beside a block-holding object a held column reads first ----
x2 <- rev(d$x2)
byHeld <- fitDraws(
  y ~ x1 + x2,
  d,
  forests = list(forest(), do.call(forest, list(basis = d$x2)))
)
blockObject <- dbartsData(y ~ x1 + x2, d, keep.data = FALSE)
expect_identical(
  fitDraws(blockObject, forests = list(forest(), forest(basis = x2))),
  byHeld
)
s1 <- dbarts(y ~ x1 + x2, d, n.trees = 5L, control = ctl, seed = 3L)
expect_identical(
  fitDraws(s1$data, forests = list(forest(), forest(basis = x2))),
  byHeld
)
rm(x2)

# ---- D3 and the plain list: a forests list states the forests in full ----
s2 <- dbarts(
  y ~ x1 + x2,
  d,
  forests = list(forest(), forest(basis = factor(z))),
  n.trees = 5L,
  control = ctl,
  seed = 3L
)
expect_identical(length(s2$data@bases), 2L)
oneForest <- dbarts(s2$data, n.trees = 5L, control = ctl, seed = 3L)
expect_identical(length(oneForest$model@forests), 1L)
expect_null(oneForest$data@bases)
expect_identical(draws(oneForest), fitDraws(y ~ x1 + x2, d))
plainTwo <- list(forest(n.trees = 4L), forest(n.trees = 3L))
twoPlain <- function(data) {
  sampler <- dbarts(data, forests = plainTwo, control = ctl, seed = 3L)
  list(
    vapply(sampler$model@forests, methods::slot, 1L, "n.trees"),
    all(vapply(sampler$data@bases, is.null, NA)),
    draws(sampler)
  )
}
byStated <- {
  sampler <- dbarts(
    y ~ x1 + x2,
    d,
    forests = plainTwo,
    control = ctl,
    seed = 3L
  )
  list(
    vapply(sampler$model@forests, methods::slot, 1L, "n.trees"),
    all(vapply(sampler$data@bases, is.null, NA)),
    draws(sampler)
  )
}
expect_identical(byStated[1:2], list(c(4L, 3L), TRUE))
expect_identical(twoPlain(s2$data), byStated)
# a data object built by the caller, from a frame or a matrix, alike
expect_identical(twoPlain(frameObject), byStated)
expect_identical(
  twoPlain(dbartsData(X[, 1:2], d$y))[1:2],
  list(c(4L, 3L), TRUE)
)
expect_identical(
  twoPlain(dbartsData(X[, 1:2], d$y))[[3L]],
  draws(dbarts(X[, 1:2], d$y, forests = plainTwo, control = ctl, seed = 3L))
)

# ---- the held basis kept by the text it was written as (dec-B706) ----
# no z anywhere but in d
expect_false(exists("z"))
byWritten <- fitDraws(
  y ~ x1 + x2,
  d,
  forests = list(forest(), forest(basis = factor(z)))
)
expect_identical(
  fitDraws(s2$data, forests = list(forest(), forest(basis = factor(z)))),
  byWritten
)
expect_identical(
  specDraws(s2$data, forests = list(forest(), forest(basis = factor(z)))),
  byWritten
)
expect_identical(
  refusal(dbarts(s2$data, forests = list(forest(), forest(basis = z)))),
  paste0(
    "'basis' (z) names 'z', which these data hold only in the basis written ",
    "as forest(basis = factor(z)); write forest(basis = factor(z)) to keep it"
  )
)

# ---- D4: a held basis keeps its recorded scale, matched by name ----
s3 <- dbarts(
  y ~ x1 + x2,
  d,
  forests = list(
    forest(),
    a = forest(basis = dose),
    b = forest(basis = factor(z))
  ),
  n.trees = 5L,
  control = ctl,
  seed = 3L
)
recorded <- dbarts:::dataBasisScale(s3$data)[[2L]]
expect_equal(recorded, sd(d$dose))
s3$setForestBasis(d$dose / 10, forest = 2L)
kept <- dbarts(
  s3$data,
  forests = list(forest(), forest(basis = dose), forest(basis = factor(z))),
  n.trees = 5L,
  control = ctl,
  seed = 3L
)
expect_identical(dbarts:::dataBasisScale(kept$data)[[2L]], recorded)
expect_identical(kept$data@bases[[2L]], s3$data@bases[[2L]])
# the swapped column in force under the recorded prior, as a sampler made from
# the swapped one's model and data runs it
expect_identical(
  draws(kept),
  draws(new("dbartsSampler", kept$control, s3$model, s3$data))
)
# moved in the list, the record follows the name
moved <- dbarts(
  s3$data,
  forests = list(forest(), forest(basis = factor(z)), forest(basis = dose)),
  n.trees = 5L,
  control = ctl,
  seed = 3L
)
expect_identical(dbarts:::dataBasisScale(moved$data)[[3L]], recorded)
expect_null(dbarts:::dataBasisScale(moved$data)[[2L]])
expect_identical(moved$data@bases[[3L]], s3$data@bases[[2L]])
expect_identical(moved$data@bases[[2L]], s3$data@bases[[3L]])

# ---- D5: at the matrix interface a name reads the matrix's column ----
x2 <- rev(d$x2)
byColumn <- fitDraws(
  X,
  d$y,
  forests = list(forest(), do.call(forest, list(basis = X[, "x2"])))
)
expect_identical(
  fitDraws(X, d$y, forests = list(forest(), forest(basis = x2))),
  byColumn
)
expect_identical(
  fitDraws(dbartsData(X, d$y), forests = list(forest(), forest(basis = x2))),
  byColumn
)
expect_identical(
  specDraws(dbartsData(X, d$y), forests = list(forest(), forest(basis = x2))),
  byColumn
)
sm <- dbarts(X, d$y, n.trees = 5L, control = ctl, seed = 3L)
expect_identical(
  fitDraws(sm$data, forests = list(forest(), forest(basis = x2))),
  byColumn
)
rm(x2)
# a column with no caller's variable of its name is read too
expect_false(exists("x2"))
expect_identical(
  fitDraws(X, d$y, forests = list(forest(), forest(basis = x2))),
  byColumn
)

# ---- a single string names its column (dec-B700) ----
byName <- fitDraws(
  y ~ x1 + x2,
  d,
  forests = list(forest(), forest(basis = dose))
)
nm <- "dose"
expect_identical(
  fitDraws(y ~ x1 + x2, d, forests = list(forest(), forest(basis = "dose"))),
  byName
)
expect_identical(
  fitDraws(y ~ x1 + x2, d, forests = list(forest(), forest(basis = nm))),
  byName
)
expect_identical(
  fitDraws(frameObject, forests = list(forest(), forest(basis = nm))),
  byName
)
byTerm <- fitDraws(y ~ x1 + forest(x2, basis = dose), d)
expect_identical(fitDraws(y ~ x1 + forest(x2, basis = "dose"), d), byTerm)
expect_identical(fitDraws(y ~ x1 + forest(x2, basis = nm), d), byTerm)
# a string vector of one value per row is a basis of its levels
arm <- ifelse(d$z == 1, "b", "a")
expect_identical(
  fitDraws(y ~ x1 + x2, d, forests = list(forest(), forest(basis = arm))),
  fitDraws(
    y ~ x1 + x2,
    d,
    forests = list(forest(), forest(basis = factor(arm)))
  )
)
# and at predict, a string given in 'bases' names a column of newdata
termFit <- bart(
  y ~ x1 + forest(x2, basis = dose),
  d,
  n.trees = 5L,
  n.samples = 5L,
  n.burn = 3L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 3L,
  verbose = FALSE,
  keepTrees = TRUE
)
newRows <- d[1:7, ]
newRows$dose <- newRows$dose + 1
byValue <- predict(termFit, newRows, bases = newRows$dose)
expect_identical(predict(termFit, newRows, bases = "dose"), byValue)
expect_identical(predict(termFit, newRows, bases = list(NULL, nm)), byValue)
expect_identical(predict(termFit, newRows), byValue)
stringTerm <- bart(
  y ~ x1 + forest(x2, basis = "dose"),
  d,
  n.trees = 5L,
  n.samples = 5L,
  n.burn = 3L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 3L,
  verbose = FALSE,
  keepTrees = TRUE
)
expect_identical(predict(stringTerm, newRows), byValue)
expect_identical(
  refusal(predict(termFit, newRows, bases = "dosage")),
  paste0(
    "'bases' gives forest 2 ('dose') the name 'dosage', which names no ",
    "column of 'newdata'"
  )
)

# ---- dec-B680: a value has one entry per row handed over ----
keep <- 1:40
cutDose <- d$dose[keep]
byFormulaSubset <- fitDraws(
  y ~ x1 + x2,
  d,
  subset = keep,
  forests = list(forest(), do.call(forest, list(basis = d$dose)))
)
byMatrixSubset <- fitDraws(
  X[, 1:2],
  d$y,
  subset = keep,
  forests = list(forest(), do.call(forest, list(basis = d$dose)))
)
expect_identical(byMatrixSubset, byFormulaSubset)
fullValue <- list(forest(), do.call(forest, list(basis = d$dose)))
cutValue <- list(forest(), do.call(forest, list(basis = cutDose)))
cutText <- paste0(
  "forest 2's 'basis' has 40 rows, matching 'subset' (40) but not the full ",
  "data (60 rows); a 'basis' must cover the full data and is subset with it"
)
for (object in list(
  dbartsData(y ~ x1 + x2, d, subset = keep),
  dbartsData(y ~ x1 + x2, d, subset = keep, keep.data = FALSE),
  dbartsData(X[, 1:2], d$y, subset = keep)
)) {
  expect_identical(fitDraws(object, forests = fullValue), byFormulaSubset)
  expect_identical(specDraws(object, forests = fullValue), byFormulaSubset)
  expect_identical(refusal(fitDraws(object, forests = cutValue)), cutText)
  expect_identical(refusal(specDraws(object, forests = cutValue)), cutText)
}

# a formula's forest() term is cut with the rows as a list's basis is
expect_identical(
  fitDraws(y ~ x1 + forest(x2, basis = dose), d, subset = keep),
  fitDraws(y ~ x1 + forest(x2, basis = dose), d[keep, ])
)
withGap <- d
withGap$x1[c(3L, 11L)] <- NA
expect_identical(
  fitDraws(y ~ x1 + forest(x2, basis = dose), withGap, na.action = na.omit),
  fitDraws(y ~ x1 + forest(x2, basis = dose), withGap[-c(3L, 11L), ])
)

# ---- a hazard fit's basis has one value per subject (dec-B707) ----
set.seed(20261109L)
h <- data.frame(
  time = sample(1:4, n, TRUE),
  status = rbinom(n, 1L, 0.7),
  x1 = rnorm(n),
  x2 = runif(n),
  dose = runif(n, 0.5, 2)
)
subject <- rep(seq_len(n), h$time)
expanded <- h[subject, ]
expanded$period <- sequence(h$time)
expanded$y <- as.integer(
  expanded$period == h$time[subject] & h$status[subject] == 1L
)
rownames(expanded) <- NULL
hazardFit <- fitDraws(
  Surv(time, status) ~ x1 + x2,
  h,
  family = "hazard",
  forests = list(forest(), forest(basis = dose))
)
expect_identical(
  hazardFit,
  fitDraws(
    y ~ x1 + x2 + period,
    expanded,
    family = "probit",
    forests = list(forest(), forest(basis = dose))
  )
)
expect_identical(
  fitDraws(
    as.matrix(h[c("x1", "x2")]),
    Surv(h$time, h$status),
    family = "hazard",
    forests = list(forest(), do.call(forest, list(basis = h$dose)))
  ),
  hazardFit
)
# 'subset' and the na.action cut the subjects first
withMissing <- h
withMissing$x1[c(4L, 9L)] <- NA
keptSubjects <- setdiff(seq_len(50L), c(4L, 9L))
keptSubject <- rep(keptSubjects, h$time[keptSubjects])
keptRows <- h[keptSubject, ]
keptRows$period <- sequence(h$time[keptSubjects])
keptRows$y <- as.integer(
  keptRows$period == h$time[keptSubject] & h$status[keptSubject] == 1L
)
rownames(keptRows) <- NULL
expect_identical(
  fitDraws(
    Surv(time, status) ~ x1 + x2,
    withMissing,
    subset = 1:50,
    na.action = na.omit,
    family = "hazard",
    forests = list(forest(), forest(basis = dose))
  ),
  fitDraws(
    y ~ x1 + x2 + period,
    keptRows,
    family = "probit",
    forests = list(forest(), forest(basis = dose))
  )
)
expect_identical(
  refusal(fitDraws(
    as.matrix(h[c("x1", "x2")]),
    Surv(h$time, h$status),
    family = "hazard",
    forests = list(forest(), do.call(forest, list(basis = expanded$dose)))
  )),
  paste0(
    "forest 2's 'basis' has ",
    nrow(expanded),
    " rows where the data were handed 60 subjects: a hazard fit's basis has ",
    "one value for each subject, before 'subset' and the na.action, and is ",
    "repeated over its periods"
  )
)

# ---- dec-B681: a declaring list beside a data object ----
declaring <- list(forest(), forest(basis = factor(z)))
byDeclared <- fitDraws(y ~ x1 + x2, d, forests = declaring)
expect_identical(fitDraws(frameObject, forests = declaring), byDeclared)
expect_identical(specDraws(frameObject, forests = declaring), byDeclared)

# ---- dec-B682 and dec-B683: a name the data were built from ----
expect_identical(
  refusal(dbarts(
    dbartsData(y ~ x1, d, keep.data = FALSE),
    forests = list(forest(), forest(basis = dose))
  )),
  paste0(
    "'basis' (dose) names 'dose', a column of the data frame these data were ",
    "built from, which they do not hold; build the data object with its ",
    "frame, dbartsData(formula, data), to read it there"
  )
)
environmentData <- list2env(as.list(d))
expect_identical(
  refusal(dbarts(
    dbartsData(y ~ x1, environmentData),
    forests = list(forest(), forest(basis = dose))
  )),
  paste0(
    "'basis' (dose) names 'dose', a variable of the environment these data ",
    "were built from, which they do not hold; build the data from a data ",
    "frame or list, which they keep, to read it there"
  )
)
# a name of neither kind is the caller's
elsewhere <- rev(d$dose)
expect_identical(
  fitDraws(
    dbartsData(y ~ x1 + x2, d, keep.data = FALSE),
    forests = list(forest(), forest(basis = elsewhere))
  ),
  fitDraws(
    y ~ x1 + x2,
    d,
    forests = list(forest(), do.call(forest, list(basis = elsewhere)))
  )
)

# ---- 'bases' is no argument of dbartsData() ----
expect_error(
  dbartsData(y ~ x1, d, bases = list(NULL, d$dose)),
  "unused argument"
)

# ---- guard: a formula with no data reads where forest() was called ----
y <- d$y
x1 <- d$x1
other <- d$dose
expect_identical(
  specDraws(
    dbartsData(y ~ x1),
    forests = list(forest(), forest(basis = other))
  ),
  fitDraws(y ~ x1, forests = list(forest(), forest(basis = other)))
)
