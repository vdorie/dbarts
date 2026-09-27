# predict's na.action (dec-B34), on every predict method and
# survivalProbabilities. Column 'a' carries training NAs, so its missing
# values have a learned route; column 'b' carries none, so its missing values
# are unroutable. The default predicts the first kind and refuses the second;
# na.pass returns NA for the second; the other functions act on any missing
# value, and na.exclude pads what it dropped back as NA.

set.seed(31)
n <- 40L
trainNames <- paste0("r", seq_len(n))
x <- matrix(
  rnorm(n * 2L),
  n,
  2L,
  dimnames = list(trainNames, c("a", "b"))
)
x[1:6, "a"] <- NA
y <- ifelse(is.na(x[, "a"]), 0, x[, "a"]) + rnorm(n, 0, 0.3)

newX <- x[7:12, ]
rownames(newX) <- paste0("t", 1:6)
newX[2L, "a"] <- NA # routable
newX[4L, "b"] <- NA # unroutable
newNames <- rownames(newX)
complete <- c(1L, 3L, 5L, 6L)
routable <- c(1L, 2L, 3L, 5L, 6L)

quick <- function(...) {
  suppressWarnings(do.call(
    bart,
    list(
      ...,
      n.samples = 8L,
      n.burn = 4L,
      n.trees = 5L,
      n.chains = 2L,
      keepTrees = TRUE,
      verbose = FALSE
    )
  ))
}
lastNames <- function(x) {
  if (is.null(dim(x))) names(x) else dimnames(x)[[length(dim(x))]]
}
# the rows a padded result keeps, along margin 'margin'
rowsOf <- function(x, keep, margin) {
  if (is.null(dim(x))) {
    return(x[keep])
  }
  indices <- rep(list(quote(expr = )), length(dim(x)))
  indices[[margin]] <- keep
  do.call(`[`, c(list(x), indices, list(drop = FALSE)))
}
seeded <- function(expr) {
  set.seed(5)
  expr
}

fit <- quick(x, y)

# --- the default, NULL and names ---

expect_error(predict(fit, newX), "'b'.*na.pass.*na.omit")
expect_error(predict(fit, newX, na.action = NULL), "'b'")
routed <- predict(fit, newX[routable, ])
expect_identical(lastNames(routed), newNames[routable])
expect_identical(predict(fit, newX[routable, ], na.action = NULL), routed)

# --- na.pass and na.exclude: positions, names and kept-row identity ---

for (type in c("ev", "ppd", "bart")) {
  for (combine in c(TRUE, FALSE)) {
    margin <- if (combine) 2L else 3L
    passed <- seeded(predict(
      fit,
      newX,
      type,
      combineChains = combine,
      na.action = na.pass
    ))
    expect_identical(lastNames(passed), newNames, info = type)
    expect_true(all(is.na(rowsOf(passed, 4L, margin))), info = type)
    expect_identical(
      rowsOf(passed, routable, margin),
      seeded(predict(fit, newX[routable, ], type, combineChains = combine)),
      info = type
    )
    excluded <- seeded(predict(
      fit,
      newX,
      type,
      combineChains = combine,
      na.action = na.exclude
    ))
    expect_identical(lastNames(excluded), newNames, info = type)
    expect_true(all(is.na(rowsOf(excluded, c(2L, 4L), margin))), info = type)
    expect_identical(
      rowsOf(excluded, complete, margin),
      seeded(predict(fit, newX[complete, ], type, combineChains = combine)),
      info = type
    )
    expect_identical(
      seeded(predict(
        fit,
        newX,
        type,
        combineChains = combine,
        na.action = na.omit
      )),
      seeded(predict(fit, newX[complete, ], type, combineChains = combine)),
      info = type
    )
  }
}
# the interval's rows are its first margin
band <- predict(fit, newX, ci.level = 0.9, na.action = na.exclude)
expect_identical(rownames(band), newNames)
expect_true(all(is.na(band[c(2L, 4L), ])))
expect_identical(
  band[complete, ],
  predict(fit, newX[complete, ], ci.level = 0.9)
)
expect_identical(
  predict(fit, newX, ci.level = 0.9, na.action = "na.omit"),
  predict(fit, newX[complete, ], ci.level = 0.9)
)
expect_null(attr(band, "na.action"))

# --- na.fail and a custom function ---

expect_error(
  predict(fit, newX, na.action = na.fail),
  "'a', 'b', which na.action = na.fail"
)
expect_identical(
  predict(fit, newX[complete, ], na.action = na.fail),
  predict(fit, newX[complete, ])
)
# a custom function sees one column, NA on the incomplete rows
seen <- NULL
dropSecond <- function(object, ...) {
  seen <<- object
  object[-2L, , drop = FALSE]
}
expect_error(predict(fit, newX, na.action = dropSecond), "'b'")
expect_identical(names(seen), "predictors")
expect_identical(which(is.na(seen$predictors)), c(2L, 4L))
dropIncomplete <- function(object, ...) {
  object[!is.na(object$predictors), , drop = FALSE]
}
expect_identical(
  predict(fit, newX, na.action = dropIncomplete),
  predict(fit, newX[complete, ])
)

# --- no surviving row, and zero-row newdata ---

seedState <- function() {
  if (exists(".Random.seed", globalenv())) {
    get(".Random.seed", globalenv())
  }
}
for (hadSeed in c(TRUE, FALSE)) {
  if (hadSeed) {
    set.seed(9)
  } else if (exists(".Random.seed", globalenv())) {
    rm(".Random.seed", envir = globalenv())
  }
  before <- seedState()
  omitted <- predict(fit, newX[4L, , drop = FALSE], "ppd", na.action = na.omit)
  expect_identical(dim(omitted), c(16L, 0L))
  padded <- predict(
    fit,
    newX[c(2L, 4L), ],
    "ppd",
    combineChains = FALSE,
    na.action = na.exclude
  )
  expect_identical(dim(padded), c(2L, 8L, 2L))
  expect_identical(lastNames(padded), newNames[c(2L, 4L)])
  expect_true(all(is.na(padded)))
  empty <- predict(fit, newX[0L, ], ci.level = 0.9)
  expect_identical(dim(empty), c(0L, 3L))
  expect_identical(seedState(), before, info = hadSeed)
}
expect_identical(dim(predict(fit, newX[0L, ], "bart")), c(16L, 0L))

# --- per-row channels ---

expect_identical(
  predict(fit, newX, offset = 0.5, na.action = na.omit),
  predict(fit, newX[complete, ], offset = 0.5)
)
offsetNew <- seq_len(6L) / 10
expect_identical(
  predict(fit, newX, offset = offsetNew, na.action = na.omit),
  predict(fit, newX[complete, ], offset = offsetNew[complete])
)
expect_error(
  predict(fit, newX, offset = 1:3, na.action = na.omit),
  "'offset' must have the same number of rows"
)
expect_identical(
  seeded(predict(fit, newX, "ppd", weights = 1:6, na.action = na.omit)),
  seeded(predict(fit, newX[complete, ], "ppd", weights = (1:6)[complete]))
)
expect_error(
  predict(fit, newX[complete, ], offset = c(0.1, NA, 0.2, 0.3)),
  "'offset' has missing values"
)
expect_error(
  predict(fit, newX[complete, ], "ppd", weights = c(1, 1, NA, 1)),
  "'weights' has missing values"
)

# --- a positional newdata warns once ---

countWarnings <- function(expr) {
  count <- 0L
  withCallingHandlers(
    expr,
    dbartsPositionalArgsWarning = function(w) {
      count <<- count + 1L
      invokeRestart("muffleWarning")
    }
  )
  count
}
expect_identical(
  countWarnings(predict(fit, unname(newX), na.action = na.exclude)),
  1L
)

# --- a formula fit with a factor, and the unseen-level fix ---

df <- data.frame(
  y = y,
  a = x[, "a"],
  b = x[, "b"],
  g = factor(rep(c("u", "v", "w"), length.out = n)),
  row.names = trainNames
)
fitF <- quick(y ~ a + b + g, df)
newDf <- data.frame(newX, g = factor(c("u", "v", "w", "u", NA, "v")))
expect_identical(
  predict(fitF, newDf, na.action = na.omit),
  predict(fitF, newDf[c(1L, 3L, 6L), ])
)
passedF <- predict(fitF, newDf, na.action = na.pass)
expect_true(all(is.na(passedF[, c(4L, 5L)])))
# an unseen level is refused, even beside a missing value in its column
unseen <- newDf
unseen$g <- factor(c("u", "z", "w", "u", NA, "v"))
expect_error(
  predict(fitF, unseen, na.action = na.pass),
  "levels not present in the training data"
)
# and in a container coded on its own level table
expect_error(
  predict(
    fitF,
    dbarts:::makeCategoricalModelMatrix(unseen[c("a", "b", "g")]),
    na.action = na.pass
  ),
  "levels not present in the training data"
)

# --- sparse newdata ---

if (requireNamespace("Matrix", quietly = TRUE)) {
  xSparse <- Matrix::Matrix(x, sparse = TRUE)
  fitS <- quick(xSparse, y)
  newSparse <- Matrix::Matrix(newX, sparse = TRUE)
  passedS <- predict(fitS, newSparse, na.action = na.pass)
  expect_identical(lastNames(passedS), newNames)
  expect_true(all(is.na(passedS[, 4L])))
  expect_identical(
    passedS[, routable],
    predict(fitS, newSparse[routable, ])
  )
}

# --- the amplitude arms ---

z <- rep(c(0, 1), length.out = n)
dfZ <- data.frame(y = y + z, a = x[, "a"], b = x[, "b"], z = z)
fitZ <- quick(y ~ a + b + z:forest(a + b), dfZ)
newZ <- data.frame(newX, z = rep(c(1, 0), 3L))
forestPad <- predict(fitZ, newZ, "forest", na.action = na.exclude)
expect_identical(dimnames(forestPad)[[2L]], newNames)
expect_true(all(is.na(forestPad[, c(2L, 4L), ])))
expect_identical(
  forestPad[, complete, , drop = FALSE],
  predict(fitZ, newZ[complete, ], "forest")
)
blendPad <- predict(fitZ, newZ, "bart", na.action = na.exclude)
expect_identical(
  blendPad[, complete],
  predict(fitZ, newZ[complete, ], "bart")
)
# a caller's own bases are checked against newdata, then cut to the kept rows
basesZ <- list(NULL, newZ$z)
expect_identical(
  predict(fitZ, newZ, "bart", bases = basesZ, na.action = na.omit),
  predict(fitZ, newZ[complete, ], "bart", bases = list(NULL, newZ$z[complete]))
)
expect_identical(
  dim(predict(fitZ, newZ[4L, ], "bart", na.action = na.omit)),
  c(16L, 0L)
)

# --- the other classes ---

category <- factor(rep(c("u", "v", "w"), length.out = n))
fitM <- quick(x, category, family = "multinomial")
for (type in c("ev", "ppd")) {
  margin <- if (type == "ev") 2L else 2L
  excludedM <- seeded(predict(fitM, newX, type, na.action = na.exclude))
  expect_identical(dimnames(excludedM)[[2L]], newNames)
  expect_identical(
    rowsOf(excludedM, complete, margin),
    seeded(predict(fitM, newX[complete, ], type)),
    info = type
  )
}
classM <- predict(fitM, newX, "class", na.action = na.exclude)
expect_identical(names(classM), newNames)
expect_true(all(is.na(classM[c(2L, 4L)])))
expect_identical(levels(classM), levels(category))
bandM <- predict(fitM, newX, ci.level = 0.9, na.action = na.exclude)
expect_identical(dimnames(bandM)[[1L]], newNames)
expect_identical(dim(predict(fitM, newX[0L, ])), c(16L, 0L, 3L))

fitO <- quick(x, factor(category, ordered = TRUE), family = "ordinal")
expect_identical(
  predict(fitO, newX, "bart", na.action = na.omit),
  predict(fitO, newX[complete, ], "bart")
)
classO <- predict(fitO, newX, "class", na.action = na.pass)
expect_true(is.na(classO[[4L]]))
expect_identical(classO[-4L], predict(fitO, newX[-4L, ], "class"))
expect_identical(
  seeded(predict(fitO, newX, "ppd", na.action = na.exclude))[, complete],
  seeded(predict(fitO, newX[complete, ], "ppd"))
)

counts <- rpois(n, 3)
fitN <- quick(x, counts, family = "nbinom")
ppdN <- seeded(predict(fitN, newX, "ppd", na.action = na.exclude))
expect_true(is.integer(ppdN) || is.double(ppdN))
expect_true(all(is.na(ppdN[, c(2L, 4L)])))
expect_identical(
  ppdN[, complete],
  seeded(predict(fitN, newX[complete, ], "ppd"))
)
expect_identical(
  predict(fitN, newX, offset = offsetNew, na.action = na.omit),
  predict(fitN, newX[complete, ], offset = offsetNew[complete])
)

positive <- ifelse(seq_len(n) %% 3L == 0L, 0, exp(y))
fitH <- quick(x, positive, family = "hurdle.lognormal")
expect_error(predict(fitH, newX), "'b'")
ppdH <- seeded(predict(fitH, newX, "ppd", na.action = na.exclude))
expect_identical(lastNames(ppdH), newNames)
expect_identical(
  ppdH[, complete],
  seeded(predict(fitH, newX[complete, ], "ppd"))
)
expect_identical(
  countWarnings(predict(fitH, unname(newX), na.action = na.omit)),
  1L
)
expect_identical(dim(predict(fitH, newX[0L, ])), c(16L, 0L))

# --- survivalProbabilities ---

if (requireNamespace("survival", quietly = TRUE)) {
  time <- rexp(n, 1)
  status <- rep(c(1, 1, 0), length.out = n)
  fitA <- quick(x, survival::Surv(time, status), family = "aft")
  expect_error(survivalProbabilities(fitA, 1, newX), "'b'")
  survA <- survivalProbabilities(fitA, c(0.5, 1), newX, na.action = na.pass)
  expect_identical(lastNames(survA), newNames)
  expect_true(all(is.na(survA[,, 4L])))
  expect_identical(
    survA[,, routable],
    survivalProbabilities(fitA, c(0.5, 1), newX[routable, ])
  )
  expect_identical(
    dim(survivalProbabilities(fitA, 1, newX[0L, ])),
    c(16L, 1L, 0L)
  )

  hazardTime <- rep(c(2.5, 0.5, 1.5), length.out = n)
  fitZH <- quick(
    x,
    survival::Surv(hazardTime, status),
    family = quote(hazard(breaks = c(0, 1, 2, 3)))
  )
  expect_error(survivalProbabilities(fitZH, newdata = newX), "'b'")
  survZ <- survivalProbabilities(
    fitZH,
    newdata = newX,
    na.action = na.exclude
  )
  expect_identical(lastNames(survZ), newNames)
  expect_true(all(is.na(survZ[,, c(2L, 4L)])))
  expect_identical(
    survZ[,, complete],
    survivalProbabilities(fitZH, newdata = newX[complete, ])
  )
  expect_identical(
    survivalProbabilities(fitZH, newdata = newX, na.action = na.omit),
    survivalProbabilities(fitZH, newdata = newX[complete, ])
  )
  emptyZ <- survivalProbabilities(
    fitZH,
    newdata = newX[4L, , drop = FALSE],
    na.action = na.exclude
  )
  expect_identical(dim(emptyZ), c(16L, 3L, 1L))
  expect_true(all(is.na(emptyZ)))

  # fit time: a formula-path hazard fit restates its subject-level record
  # over the person-period rows the dropped subject had, as the matrix path,
  # which expands first, records them
  hazardX <- x[7:12, ]
  hazardX[4L, "b"] <- NA
  hazardFrame <- data.frame(
    time = hazardTime[1:6],
    status = status[1:6],
    a = hazardX[, "a"],
    b = hazardX[, "b"],
    row.names = rownames(hazardX)
  )
  for (action in list(na.exclude, na.omit)) {
    fitMatrix <- quick(
      hazardX,
      survival::Surv(hazardTime[1:6], status[1:6]),
      na.action = action,
      family = quote(hazard(breaks = c(0, 1, 2, 3)))
    )
    fitFormula <- quick(
      survival::Surv(time, status) ~ a + b,
      hazardFrame,
      na.action = action,
      family = quote(hazard(breaks = c(0, 1, 2, 3)))
    )
    expect_identical(fitFormula$na.action, fitMatrix$na.action)
    expect_identical(fitFormula$row.names.train, fitMatrix$row.names.train)
    expect_identical(names(fitted(fitFormula)), names(fitted(fitMatrix)))
  }
  expect_identical(names(fitMatrix$na.action), c("r10", "r10.1", "r10.2"))

  # the default grid comes from the kept subjects on both paths: the dropped
  # subject's time, the only 3.5, adds no period
  hazardFrame$time[4L] <- 3.5
  fitMatrix <- quick(
    hazardX,
    survival::Surv(hazardFrame$time, status[1:6]),
    na.action = na.exclude,
    family = "hazard"
  )
  fitFormula <- quick(
    survival::Surv(time, status) ~ a + b,
    hazardFrame,
    na.action = na.exclude,
    family = "hazard"
  )
  expect_identical(fitMatrix$periods, c(0.5, 1.5, 2.5))
  expect_identical(fitFormula$periods, fitMatrix$periods)
  expect_identical(fitFormula$na.action, fitMatrix$na.action)
  expect_identical(names(fitted(fitFormula)), names(fitted(fitMatrix)))
}
