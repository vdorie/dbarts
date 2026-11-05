#!/usr/bin/env Rscript

# The hurdle reduction gate (docs/design/hurdle.md section 5) - hurdle
# composes two ordinary bart2() fits and never combines them within a sweep
# (docs/design/hurdle.md section 0), so a bart2(family = "hurdle.lognormal")
# fit's two component fits must equal, DRAW FOR DRAW, standalone
# bart2(family = "probit") / bart2(family = "gaussian") fits built by hand at
# the same derived seeds and control. This is a sanity FLOOR: the wrapper
# literally calls bart2() for each component (R/bart.R's bart2Hurdle), so the
# equality holds by construction. It is not the full correctness argument -
# the analytic combine/retransform oracle is the load-bearing gate and lands
# separately (docs/plans/hurdle.md step 2).
#
# Usage: Rscript benchmarks/R/hurdle-reduction.R

suppressPackageStartupMessages(library(dbarts))

set.seed(1729L)
n <- 250L
p <- 3L
x <- matrix(runif(n * p), n, p)
pi.true <- pnorm(0.8 * x[, 1L] - 0.5 * x[, 2L])
mu.true <- 1.2 + 0.6 * x[, 1L] - 0.4 * x[, 3L]
occupied <- rbinom(n, 1L, pi.true) == 1L
y <- numeric(n)
y[occupied] <- exp(rnorm(sum(occupied), mu.true[occupied], 0.4))

args <- list(
  n.trees = 30L,
  n.burn = 40L,
  n.samples = 80L,
  n.chains = 1L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE,
  seed = 90210L
)

# the SAME derivation bart2Hurdle uses internally (R/bart.R), reproduced here
# so the standalone comparison fits are built at the identical seeds
derivedSeeds <- local({
  oldSeed <- .GlobalEnv[[".Random.seed"]]
  set.seed(args$seed)
  seeds <- sample.int(.Machine$integer.max, 2L)
  if (!is.null(oldSeed)) {
    .GlobalEnv$.Random.seed <- oldSeed
  }
  seeds
})

fitHurdle <- do.call(bart2, c(list(x, y, family = "hurdle.lognormal"), args))

split <- dbarts:::splitHurdleResponse(y)
xPositive <- x[split$positive, , drop = FALSE]

zeroArgs <- args
zeroArgs$seed <- derivedSeeds[1L]
fitZero <- do.call(
  bart2,
  c(list(x, split$z, family = "probit"), zeroArgs)
)

positiveArgs <- args
positiveArgs$seed <- derivedSeeds[2L]
fitPositive <- do.call(
  bart2,
  c(
    list(xPositive, split$logPositive, test = x, family = "gaussian"),
    positiveArgs
  )
)

compareComponent <- function(wrapped, standalone, label, hasSigma) {
  results <- c(
    # values only: row names on the outputs are not draws
    yhat.train = identical(
      unname(extract(wrapped, type = "bart", sample = "train")),
      unname(extract(standalone, type = "bart", sample = "train"))
    ),
    varcount = identical(wrapped$varcount, standalone$varcount),
    trees = identical(
      extract(wrapped, type = "trees"),
      extract(standalone, type = "trees")
    )
  )
  if (hasSigma) {
    results <- c(results, sigma = identical(wrapped$sigma, standalone$sigma))
  }
  bitwise <- all(results)
  cat(sprintf(
    "%-10s reduction %s (%s)\n",
    label,
    if (bitwise) "BITWISE IDENTICAL" else "DIVERGED",
    paste(sprintf("%s=%s", names(results), results), collapse = ", ")
  ))
  bitwise
}

okZero <- compareComponent(
  fitHurdle$zero,
  fitZero,
  "zero",
  hasSigma = FALSE
)
okPositive <- compareComponent(
  fitHurdle$positive,
  fitPositive,
  "positive",
  hasSigma = TRUE
)

# the packaged objects differ ONLY in $call (the wrapper's own reconstructed
# component call vs. the by-hand standalone call, the hazard-reduction.R
# markerOnly idea); every other top-level field runs through the identical
# bart2() code path
markerOnly <-
  identical(fitHurdle$zero$family, "probit") &&
  identical(fitHurdle$positive$family, "gaussian") &&
  identical(fitHurdle$family, "hurdle.lognormal")

if (!markerOnly) {
  cat("  marker check FAILED: the objects differ beyond $call\n")
}

# ---- arms for the inputs a hurdle fit takes beyond a matrix and a response ----
#
# A formula fit, a test set, a subset, weights of 0 and 1, and each part's
# offset. The standalone inputs below are derived by hand from the data, not by
# splitHurdleResponse or the arm's own helpers: the rows a subset keeps, the
# 1{y > 0} and log y split, the positive rows' weights and each part's offset
# are written out here, so a defect in the arm's own derivation is not copied
# into its reference. The positive part's test is the n training rows followed
# by the test rows, and a hand-built positive-part fit cut from a data frame
# keeps the same factor levels as the hurdle's one design.

set.seed(4242L)
nArm <- 220L
dArm <- data.frame(a = runif(nArm), b = runif(nArm))
dArm$f <- factor(sample(c("u", "v", "w"), nArm, TRUE))
dArm$t <- runif(nArm, 1, 3)
dArm$u <- runif(nArm, -0.5, 0.5)
piArm <- pnorm(0.9 * dArm$a - 0.5 * dArm$b + 0.3 * (dArm$f == "u") - 0.2)
inArm <- rbinom(nArm, 1L, piArm) == 1L
dArm$y <- 0
dArm$y[inArm] <- exp(rnorm(sum(inArm), 1 + 0.6 * dArm$a[inArm], 0.4))
# a level seen on zero rows alone does not reach the positive part's rows, so
# the positive fit must code it as the hurdle's one design does
dArm$f <- factor(
  replace(as.character(dArm$f), which(dArm$y == 0)[1:2], "z"),
  levels = c("u", "v", "w", "z")
)
# weights of 0 and 1, the 0 rows on both kinds of row
wArm <- rep(1, nArm)
wArm[c(which(dArm$y == 0)[3:6], which(dArm$y > 0)[1:5])] <- 0

nTest <- 40L
tArm <- data.frame(a = runif(nTest), b = runif(nTest))
tArm$f <- factor(sample(c("u", "v", "w"), nTest, TRUE), levels = levels(dArm$f))
tArm$t <- runif(nTest, 1, 3)
tArm$u <- runif(nTest, -0.5, 0.5)

armChecks <- function(label, wrapped, standalone, hasSigma, hasTest) {
  results <- c(
    yhat.train = identical(
      unname(wrapped$yhat.train),
      unname(standalone$yhat.train)
    ),
    varcount = identical(wrapped$varcount, standalone$varcount),
    trees = identical(
      extract(wrapped, type = "trees"),
      extract(standalone, type = "trees")
    )
  )
  if (hasTest) {
    results <- c(
      results,
      yhat.test = identical(
        unname(wrapped$yhat.test),
        unname(standalone$yhat.test)
      )
    )
  }
  if (hasSigma) {
    results <- c(results, sigma = identical(wrapped$sigma, standalone$sigma))
  }
  ok <- all(results)
  cat(sprintf(
    "%-24s %s\n",
    label,
    if (ok) {
      "BITWISE IDENTICAL"
    } else {
      paste(
        "DIVERGED:",
        paste(
          names(results)[!results],
          collapse = ", "
        )
      )
    }
  ))
  ok
}

# fit one arm: the hurdle call, then each part's standalone fit from the hand
# built inputs, at the derived seeds
runArm <- function(label, hurdle, zero, positive, hasTest = FALSE) {
  # test weights are not part of what these arms compare
  quietly <- function(expr) {
    withCallingHandlers(expr, warning = function(w) {
      if (grepl("not found in test data", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    })
  }
  wrapped <- quietly(do.call(
    bart2,
    c(hurdle, list(family = "hurdle.lognormal"), args)
  ))
  zeroAlone <- quietly(do.call(
    bart2,
    c(
      zero,
      list(family = "probit"),
      args[names(args) != "seed"],
      list(seed = derivedSeeds[1L])
    )
  ))
  positiveAlone <- quietly(do.call(
    bart2,
    c(
      positive,
      list(family = "gaussian"),
      args[names(args) != "seed"],
      list(seed = derivedSeeds[2L])
    )
  ))
  c(
    armChecks(paste(label, "zero"), wrapped$zero, zeroAlone, FALSE, hasTest),
    armChecks(
      paste(label, "positive"),
      wrapped$positive,
      positiveAlone,
      TRUE,
      TRUE
    )
  )
}

arms <- list()
isPos <- dArm$y > 0

# formula fit
dZero <- transform(dArm, z = as.double(y > 0))
dPos <- transform(dArm[isPos, ], ly = log(y))
arms$formula <- runArm(
  "formula",
  list(y ~ a + b + f, data = dArm),
  list(z ~ a + b + f, data = dZero),
  list(ly ~ a + b + f, data = dPos, test = dArm)
)

# test set, on the matrix route
xArm <- as.matrix(dArm[, c("a", "b")])
xTest <- as.matrix(tArm[, c("a", "b")])
yZero <- as.double(dArm$y > 0)
arms$test <- runArm(
  "test set",
  list(xArm, dArm$y, test = xTest),
  list(xArm, yZero, test = xTest),
  list(
    xArm[isPos, , drop = FALSE],
    log(dArm$y[isPos]),
    test = rbind(xArm, xTest)
  ),
  hasTest = TRUE
)

# subset, with a missing predictor value that the default keeps
keep <- dArm$a > 0.15 & dArm$b < 0.9
dSub <- dArm[keep, ]
dSub$b[c(1L, 5L)] <- NA
dSubIn <- dArm
dSubIn$b[which(keep)[c(1L, 5L)]] <- NA
isPosSub <- dSub$y > 0
dSubZero <- transform(dSub, z = as.double(y > 0))
dSubPos <- transform(dSub[isPosSub, ], ly = log(y))
arms$subset <- runArm(
  "subset",
  list(y ~ a + b + f, data = dSubIn, subset = keep),
  list(z ~ a + b + f, data = dSubZero),
  list(ly ~ a + b + f, data = dSubPos, test = dSub)
)

# weights of 0 and 1: the positive part takes the positive rows' weights
arms$weights <- runArm(
  "weights 0/1",
  list(y ~ a + b + f, data = dArm, weights = wArm),
  list(z ~ a + b + f, data = dZero, weights = wArm),
  list(ly ~ a + b + f, data = dPos, weights = wArm[isPos], test = dArm)
)

# each part's offset: the zero offset on the probit latent, the positive one on
# log y, and its test channel the training rows' offset beside the test rows'
arms$offset <- runArm(
  "offsets",
  list(
    y ~ a + b + f,
    data = dArm,
    offset = list(zero = dArm$u, positive = log(dArm$t))
  ),
  list(z ~ a + b + f, data = dZero, offset = dArm$u),
  list(
    ly ~ a + b + f,
    data = dPos,
    offset = log(dArm$t[isPos]),
    test = dArm,
    offset.test = log(dArm$t)
  )
)
arms$offsetPlain <- runArm(
  "plain offset",
  list(y ~ a + b + f, data = dArm, offset = log(dArm$t)),
  list(z ~ a + b + f, data = dZero),
  list(
    ly ~ a + b + f,
    data = dPos,
    offset = log(dArm$t[isPos]),
    test = dArm,
    offset.test = log(dArm$t)
  )
)

# all of them at once: a subset, weights, a test set and both offsets
kept <- which(keep)
dAll <- dArm[kept, ]
wAll <- wArm[kept]
isPosAll <- dAll$y > 0
dAllZero <- transform(dAll, z = as.double(y > 0))
dAllPos <- transform(dAll[isPosAll, ], ly = log(y))
arms$combined <- runArm(
  "combined",
  list(
    y ~ a + b + f,
    data = dArm,
    subset = keep,
    weights = wArm,
    offset = list(zero = dArm$u, positive = log(dArm$t)),
    test = tArm,
    offset.test = list(zero = tArm$u, positive = log(tArm$t))
  ),
  list(
    z ~ a + b + f,
    data = dAllZero,
    weights = wAll,
    offset = dAll$u,
    test = tArm,
    offset.test = tArm$u
  ),
  list(
    ly ~ a + b + f,
    data = dAllPos,
    weights = wAll[isPosAll],
    offset = log(dAll$t[isPosAll]),
    test = rbind(dAll[, names(tArm)], tArm),
    offset.test = c(log(dAll$t), log(tArm$t))
  ),
  hasTest = TRUE
)

okArms <- all(unlist(arms))

if (okZero && okPositive && markerOnly && okArms) {
  cat("\nOK: hurdle reduces bitwise to its two standalone component fits\n")
} else {
  cat("\nFAIL: a hurdle component diverged from its standalone reduction\n")
  quit(status = 1L)
}
