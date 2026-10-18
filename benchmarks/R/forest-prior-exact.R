#!/usr/bin/env Rscript

# Exact-posterior gate for the rule a forest's coefficient and size follow in
# a model of several forests (docs/plans/forest-multiplier.md, "The rule"),
# its held arms. One binary predictor and sigma held, so a single-tree forest
# is a stump or one split and its leaves integrate out in closed form: given
# how many of each forest's trees split, the forests' cell values are jointly
# normal, the response is their sum through the held multipliers plus noise,
# and the posterior over the configurations is a finite sum. No quadrature.
#
# The oracle is written from the rule and reads nothing off the sampler but
# the unit and the shift: a held coefficient is 1 with no basis and (0, 1)
# for the two levels of a factor, so the forest enters with multiplier 1 or
# z; the forest's own fit then has prior sd its sd, exactly, each of T trees'
# leaves N(0, sd^2 / T); a held forest with no basis that states no sd has a
# quarter of the response's range, taken here from the response drawn, and a
# held factor one unit. Each arm compares the posterior mean of each forest's
# term in each cell, the second level's for a factor, and the probability
# that each forest's tree splits.
#
# Verdict: an arm fails when a gap exceeds 0.0005 plus 4 standard errors over
# seeds. With about forty functionals a run a lone z between 4 and 5 is
# expected now and then, so such an arm is run once more on a second set of
# seeds and fails only if the same functional is out again.
#
# The numeric arms follow (a numeric multiplier, its coefficient drawn or
# held), and the functionals there are the posterior mean of the fit in each
# class of rows that share a cell and the basis's value, the probability that
# each tree splits and that a drawn coefficient is within 0.5, 1 and 2.
#
# Not seen here: a drawn coefficient on a forest with no basis or a factor,
# probit and logistic, subset, a numeric multiplier under logistic, and n for
# n - 1 in a column's scale on a held arm, which the scale's literals in
# tests/cpp and the tinytest record's 1e-12 hold.
#
# Usage: Rscript forest-prior-exact.R [quick]

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args
numSeeds <- if (quick) 32L else 64L
numKept <- if (quick) 25000L else 50000L
numBurn <- 500L
numCores <- as.integer(Sys.getenv("EQUIVALENCE_CORES", "2"))
gapFloor <- 0.0005
zBound <- 4
zRerun <- 5

# ---- the design ----

set.seed(20261017L)
perGroup <- 4L
cell <- rep(0:1, each = 4L * perGroup)
z1 <- rep(rep(c(0, 1, 0, 1), each = perGroup), 2L)
z2 <- rep(rep(c(0, 0, 1, 1), each = perGroup), 2L)
n <- length(cell)
x <- matrix(as.double(cell), ncol = 1L, dimnames = list(NULL, "x"))
sigma <- 0.5 # the residual sd, held; fixed() takes its square
y <- 1 + 0.8 * cell + z1 * (0.5 + 0.4 * cell) - 0.3 * z2 + rnorm(n, sd = sigma)
base <- 0.5
power <- 2

# ---- the oracle ----

# forests: a list of list(multiplier = <n-vector>, trees = T, sd = <response
# units>). The forest's two cell values G_f = (G_f(0), G_f(1)) have, given k
# of its T trees split, covariance (sd^2 / T) [[T, T - k], [T - k, T]], and k
# is binomial(T, base). y - shift = sum_f m_f G_f(cell) + N(0, sigma^2).
exactPosterior <- function(forests, shift) {
  K <- length(forests)
  design <- matrix(0, n, 2L * K)
  for (f in seq_len(K)) {
    design[cbind(seq_len(n), 2L * (f - 1L) + cell + 1L)] <- forests[[
      f
    ]]$multiplier
  }
  counts <- as.matrix(expand.grid(lapply(forests, function(forest) {
    0:forest$trees
  })))
  residual <- y - shift
  logWeights <- numeric(nrow(counts))
  means <- matrix(0, nrow(counts), 2L * K)
  for (r in seq_len(nrow(counts))) {
    prior <- matrix(0, 2L * K, 2L * K)
    logPrior <- 0
    for (f in seq_len(K)) {
      trees <- forests[[f]]$trees
      k <- counts[r, f]
      leaf <- forests[[f]]$sd^2 / trees
      block <- leaf * matrix(c(trees, trees - k, trees - k, trees), 2L)
      at <- 2L * (f - 1L) + 1:2
      prior[at, at] <- block
      logPrior <- logPrior + dbinom(k, trees, base, log = TRUE)
    }
    marginal <- design %*% prior %*% t(design) + diag(sigma^2, n)
    chol <- chol(marginal)
    solved <- backsolve(chol, forwardsolve(t(chol), residual))
    logWeights[r] <- logPrior -
      sum(log(diag(chol))) -
      0.5 * sum(residual * solved)
    means[r, ] <- prior %*% t(design) %*% solved
  }
  weights <- exp(logWeights - max(logWeights))
  weights <- weights / sum(weights)
  term <- colSums(weights * means)
  splits <- vapply(
    seq_len(K),
    function(f) sum(weights * counts[, f]) / forests[[f]]$trees,
    0
  )
  names(term) <- paste0(
    "term",
    rep(seq_len(K), each = 2L),
    "[",
    rep(0:1, K),
    "]"
  )
  names(splits) <- paste0("split", seq_len(K))
  c(term, splits)
}

# ---- the sampler ----

# `arm`: a list of forests, each list(basis = NULL, z1 or z2, sd = NULL or a
# number of units), and n.trees
heldForest <- function(spec, trees) {
  arguments <- list(
    n.trees = trees,
    base = base,
    power = power,
    amplitude = dbarts::dbartsPriors$fixed()
  )
  if (!is.null(spec$basis)) {
    arguments$basis <- factor(spec$basis, levels = c(0, 1))
  }
  if (!is.null(spec$sd)) {
    arguments$sd <- spec$sd
  }
  do.call(dbarts::dbartsForests$forest, arguments)
}

# the forests' terms in each cell from the combined train draws: a held
# forest's multiplier is 1 or its z, so the rows at z = 0 of every factor
# carry the forest with no basis and each factor adds its rows at z = 1
termDraws <- function(train, shift, arm) {
  K <- length(arm$forests)
  plain <- matrix(0, 2L, ncol(train))
  out <- list()
  baseRows <- function(c) which(cell == c & z1 == 0 & z2 == 0)
  for (c in 0:1) {
    plain[c + 1L, ] <- colMeans(train[baseRows(c), , drop = FALSE]) - shift
  }
  for (f in seq_len(K)) {
    spec <- arm$forests[[f]]
    for (c in 0:1) {
      name <- sprintf("term%d[%d]", f, c)
      if (is.null(spec$basis)) {
        out[[name]] <- plain[c + 1L, ]
      } else {
        on <- if (identical(spec$basis, z1)) {
          z1 == 1 & z2 == 0
        } else {
          z2 == 1 & z1 == 0
        }
        rows <- which(cell == c & on)
        out[[name]] <- colMeans(train[rows, , drop = FALSE]) -
          shift -
          plain[c + 1L, ]
      }
    }
  }
  out
}

samplerFit <- function(seed, arm, kept) {
  set.seed(seed)
  control <- dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.cuts = 1L,
    updateState = FALSE,
    keepTrainingFits = TRUE
  )
  sampler <- dbarts(
    x,
    y,
    family = gaussian(sigma = dbartsPriors$fixed(sigma^2)),
    forests = lapply(arm$forests, heldForest, trees = arm$trees),
    control = control
  )
  result <- sampler$run(numBurn, kept)
  shift <- sampler$getLeafPrior(1L)$response.shift
  terms <- termDraws(result$train, shift, arm)
  counts <- result$varcount
  K <- length(arm$forests)
  splits <- lapply(seq_len(K), function(f) counts[1L, f, ] / arm$trees)
  names(splits) <- paste0("split", seq_len(K))
  vapply(c(terms, splits), mean, 0)
}

# the unit and the shift, the two numbers the oracle takes off a sampler
samplerConstants <- function(arm) {
  sampler <- dbarts(
    x,
    y,
    family = gaussian(sigma = dbartsPriors$fixed(sigma^2)),
    forests = lapply(arm$forests, heldForest, trees = arm$trees),
    control = dbartsControl(n.chains = 1L, n.threads = 1L, n.cuts = 1L)
  )
  list(
    unit = sampler$data@forest.calibration[["unit"]],
    shift = sampler$getLeafPrior(1L)$response.shift
  )
}

# the oracle's forests for an arm: a stated sd as written, a default by the
# rule, from the response drawn
oracleForests <- function(arm, unit) {
  lapply(arm$forests, function(spec) {
    sd <- if (!is.null(spec$sd)) {
      spec$sd
    } else if (is.null(spec$basis)) {
      diff(range(y)) / 4
    } else {
      unit
    }
    list(
      multiplier = if (is.null(spec$basis)) rep(1, n) else spec$basis,
      trees = arm$trees,
      sd = sd
    )
  })
}

runArm <- function(name, arm, exact = NULL, fit = samplerFit) {
  if (is.null(exact)) {
    constants <- samplerConstants(arm)
    exact <- exactPosterior(
      oracleForests(arm, constants$unit),
      constants$shift
    )
  }
  verdict <- function(seeds) {
    fits <- parallel::mclapply(
      seeds,
      fit,
      arm = arm,
      kept = numKept,
      mc.cores = numCores
    )
    fits <- do.call(rbind, fits)
    mean <- colMeans(fits)[names(exact)]
    se <- apply(fits, 2L, sd)[names(exact)] / sqrt(nrow(fits))
    gap <- mean - exact
    list(
      mean = mean,
      se = se,
      gap = gap,
      z = gap / se,
      out = abs(gap) > gapFloor + zBound * se
    )
  }
  started <- proc.time()[["elapsed"]]
  first <- verdict(seq_len(numSeeds))
  failed <- first$out
  rerun <- any(failed) && all(abs(first$z[failed]) < zRerun)
  if (rerun) {
    second <- verdict(numSeeds + seq_len(numSeeds))
    failed <- failed & second$out
  }
  cat(sprintf("%s (%.1f s)\n", name, proc.time()[["elapsed"]] - started))
  cat(sprintf(
    "  %-10s %9s %9s %9s %7s\n",
    "quantity",
    "sampler",
    "exact",
    "se",
    "z"
  ))
  for (label in names(exact)) {
    cat(sprintf(
      "  %-10s %+9.5f %+9.5f %9.2e %+7.2f%s\n",
      label,
      first$mean[[label]],
      exact[[label]],
      first$se[[label]],
      first$z[[label]],
      if (failed[[label]]) "  <- FAIL" else ""
    ))
  }
  if (rerun) {
    cat("  a lone z under", zRerun, "was run again on a second set of seeds\n")
  }
  any(failed)
}

# ---- the numeric arms ----

# A numeric forest enters as a w_i F(x_i): a drawn coefficient a is N(0, 1)
# and F has prior sd s / 0.674, s its sd per unit of w, which is one unit of
# the response's scale for one standard deviation of w where none is stated,
# that standard deviation R's sd() over the rows of positive weight; a held
# coefficient is 1 and F has prior sd s. Given a and how many of each forest's
# trees split the cell values are jointly normal, so the only quadrature is
# over a, on a grid under its N(0, 1) prior. A row of weight zero leaves the
# likelihood; the unit is read over the rows in it.
#
# The columns, by group of the design: a 0/1 column, and one of four values
# whose standard deviation is about 0.11. The weighted arms put a fifth of the
# rows at weight zero, their column far out.
group <- 2L * z2 + z1
binaryColumn <- as.double(group %% 2L == 1L)
fourValues <- c(0, 0.1, 0.2, 0.3)[group + 1L]
zeroWeight <- seq_len(n) %% 5L == 3L
weightedColumn <- ifelse(zeroWeight, 2, fourValues)
caseWeights <- ifelse(zeroWeight, 0, 1)
gridA <- seq(-8, 8, by = 0.005)

# forests: a list of list(multiplier = <n-vector>, drawn = TRUE where the
# coefficient is drawn, trees, sd = the prior sd of the forest's own fit).
# rows: the rows in the likelihood. Returns the posterior mean of the fit in
# each class of rows that share a cell and every multiplier, the probability
# that each forest's tree splits, and where a coefficient is drawn the
# probability that it is within 0.5, 1 and 2.
numericPosterior <- function(forests, shift, rows, classes) {
  K <- length(forests)
  counts <- as.matrix(expand.grid(lapply(forests, function(forest) {
    0:forest$trees
  })))
  residual <- (y - shift)[rows]
  drawn <- which(vapply(forests, `[[`, NA, "drawn"))
  values <- if (length(drawn) > 0L) gridA else 1
  logPriorA <- if (length(drawn) > 0L) dnorm(values, log = TRUE) else 0
  # the design at a = 1: row i's column for forest f's cell value
  designAll <- matrix(0, n, 2L * K)
  for (f in seq_len(K)) {
    designAll[cbind(seq_len(n), 2L * (f - 1L) + cell + 1L)] <-
      forests[[f]]$multiplier
  }
  priors <- lapply(seq_len(nrow(counts)), function(r) {
    prior <- matrix(0, 2L * K, 2L * K)
    for (f in seq_len(K)) {
      trees <- forests[[f]]$trees
      k <- counts[r, f]
      leaf <- forests[[f]]$sd^2 / trees
      at <- 2L * (f - 1L) + 1:2
      prior[at, at] <- leaf * matrix(c(trees, trees - k, trees - k, trees), 2L)
    }
    prior
  })
  logPriorK <- vapply(
    seq_len(nrow(counts)),
    function(r) {
      sum(vapply(
        seq_len(K),
        function(f) dbinom(counts[r, f], forests[[f]]$trees, base, log = TRUE),
        0
      ))
    },
    0
  )
  logWeights <- matrix(0, length(values), nrow(counts))
  classFits <- array(0, c(length(values), nrow(counts), length(classes)))
  rr <- sum(residual^2) / sigma^2
  for (v in seq_along(values)) {
    scaling <- rep(1, 2L * K)
    if (length(drawn) > 0L) {
      scaling[2L * (drawn - 1L) + 1:2] <- values[v]
    }
    design <- sweep(designAll, 2L, scaling, `*`)
    inDesign <- design[rows, , drop = FALSE]
    M <- crossprod(inDesign) / sigma^2
    b <- crossprod(inDesign, residual) / sigma^2
    for (r in seq_len(nrow(counts))) {
      P <- priors[[r]]
      inner <- diag(2L * K) + M %*% P
      solved <- solve(inner, b)
      logWeights[v, r] <- logPriorA[v] +
        logPriorK[r] -
        0.5 * as.numeric(determinant(inner)$modulus) -
        0.5 * (rr - sum(b * (P %*% solved)))
      mean <- P %*% solved
      fits <- design %*% mean
      classFits[v, r, ] <- vapply(classes, function(at) mean(fits[at]), 0)
    }
  }
  weights <- exp(logWeights - max(logWeights))
  # the grid's midpoint rule; the ends carry no mass
  weights <- weights / sum(weights)
  out <- vapply(
    seq_along(classes),
    function(j) sum(weights * classFits[,, j]),
    0
  )
  names(out) <- names(classes)
  splits <- vapply(
    seq_len(K),
    function(f) {
      sum(sweep(weights, 2L, counts[, f], `*`)) / forests[[f]]$trees
    },
    0
  )
  names(splits) <- paste0("split", seq_len(K))
  out <- c(out, splits)
  if (length(drawn) > 0L) {
    # each grid point stands for the interval of one step around it, so a
    # point at a bound counts half
    marginal <- rowSums(weights)
    step <- values[2L] - values[1L]
    within <- vapply(
      c(0.5, 1, 2),
      function(t) {
        sum(marginal * pmin(1, pmax(0, (t - abs(values)) / step + 0.5)))
      },
      0
    )
    names(within) <- paste0("within", c(0.5, 1, 2))
    out <- c(out, within)
  }
  out
}

# the classes of rows sharing a cell and every forest's multiplier
rowClasses <- function(arm) {
  key <- paste(
    cell,
    do.call(
      paste,
      lapply(arm$forests, function(spec) {
        if (is.null(spec$basis)) 1 else spec$basis
      })
    )
  )
  classes <- split(seq_len(n), key)
  names(classes) <- paste0("fit[", seq_along(classes), "]")
  classes
}

numericForest <- function(spec, trees) {
  arguments <- list(n.trees = trees, base = base, power = power)
  if (!is.null(spec$basis)) {
    arguments$basis <- if (isTRUE(spec$levels)) {
      factor(spec$basis, levels = c(0, 1))
    } else {
      spec$basis
    }
  }
  if (!isTRUE(spec$drawn)) {
    arguments$amplitude <- dbarts::dbartsPriors$fixed()
  }
  if (!is.null(spec$sd)) {
    arguments$sd <- spec$sd
  }
  do.call(dbarts::dbartsForests$forest, arguments)
}

numericSampler <- function(arm, control) {
  dbarts(
    x,
    y,
    weights = arm$weights,
    family = gaussian(sigma = dbartsPriors$fixed(sigma^2)),
    forests = lapply(arm$forests, numericForest, trees = arm$trees),
    control = control
  )
}

numericFit <- function(seed, arm, kept) {
  set.seed(seed)
  sampler <- numericSampler(
    arm,
    dbartsControl(
      n.chains = 1L,
      n.threads = 1L,
      n.cuts = 1L,
      updateState = FALSE,
      keepTrainingFits = TRUE
    )
  )
  result <- sampler$run(numBurn, kept)
  shift <- sampler$getLeafPrior(1L)$response.shift
  classes <- rowClasses(arm)
  fits <- vapply(
    classes,
    function(at) mean(result$train[at, , drop = FALSE]) - shift,
    0
  )
  K <- length(arm$forests)
  splits <- vapply(
    seq_len(K),
    function(f) mean(result$varcount[1L, f, ]) / arm$trees,
    0
  )
  names(splits) <- paste0("split", seq_len(K))
  out <- c(fits, splits)
  drawn <- which(vapply(arm$forests, function(spec) isTRUE(spec$drawn), NA))
  if (length(drawn) > 0L) {
    widths <- vapply(
      arm$forests,
      function(spec) if (isTRUE(spec$levels)) 2L else 1L,
      0L
    )
    coefficient <- result$glue[sum(widths[seq_len(drawn)]), ]
    within <- vapply(c(0.5, 1, 2), function(t) mean(abs(coefficient) < t), 0)
    names(within) <- paste0("within", c(0.5, 1, 2))
    out <- c(out, within)
  }
  out
}

# the oracle's forests: a stated sd as written, a default by the rule, from
# the response drawn and R's sd() over the rows in the likelihood
numericOracle <- function(arm) {
  sampler <- numericSampler(
    arm,
    dbartsControl(n.chains = 1L, n.threads = 1L, n.cuts = 1L)
  )
  unit <- sampler$data@forest.calibration[["unit"]]
  shift <- sampler$getLeafPrior(1L)$response.shift
  rows <- if (is.null(arm$weights)) seq_len(n) else which(arm$weights > 0)
  forests <- lapply(arm$forests, function(spec) {
    drawn <- isTRUE(spec$drawn)
    numeric <- !is.null(spec$basis) && !isTRUE(spec$levels)
    size <- if (!is.null(spec$sd)) {
      spec$sd
    } else if (is.null(spec$basis)) {
      diff(range(y[rows])) / 4
    } else if (numeric) {
      unit / sd(spec$basis[rows])
    } else {
      unit
    }
    list(
      multiplier = if (is.null(spec$basis)) rep(1, n) else spec$basis,
      drawn = drawn,
      trees = arm$trees,
      sd = if (drawn) size / 0.674 else size
    )
  })
  numericPosterior(forests, shift, rows, rowClasses(arm))
}

# ---- the arms ----

constants <- samplerConstants(list(
  forests = list(list(), list(basis = z1)),
  trees = 1L
))
unit <- constants$unit
cat(sprintf(
  "forest prior exact gate (%s): %d seeds of %d kept sweeps, unit %.4f\n",
  if (quick) "quick" else "full",
  numSeeds,
  numKept,
  unit
))
arms <- list(
  "no basis held at 0.3 of the unit" = list(
    forests = list(list(sd = 0.3 * unit), list(basis = z1)),
    trees = 1L
  ),
  "no basis held at 2 of the unit" = list(
    forests = list(list(sd = 2 * unit), list(basis = z1)),
    trees = 1L
  ),
  "no basis held at its default" = list(
    forests = list(list(), list(basis = z1, sd = 0.7 * unit)),
    trees = 1L
  ),
  "factor held at 0.6 of the unit" = list(
    forests = list(list(sd = 1.1 * unit), list(basis = z1, sd = 0.6 * unit)),
    trees = 1L
  ),
  "factor held at 0.2 of the unit" = list(
    forests = list(list(sd = 1.1 * unit), list(basis = z1, sd = 0.2 * unit)),
    trees = 1L
  ),
  "factor held at its default" = list(
    forests = list(list(sd = 1.1 * unit), list(basis = z1)),
    trees = 1L
  ),
  "both held at their defaults" = list(
    forests = list(list(), list(basis = z1)),
    trees = 1L
  ),
  "three forests of three trees" = list(
    forests = list(list(), list(basis = z1), list(basis = z2, sd = 0.4 * unit)),
    trees = 3L
  )
)
anyFailure <- FALSE
for (name in names(arms)) {
  anyFailure <- runArm(name, arms[[name]]) || anyFailure
}

# a numeric column at the second forest, beside a forest with no basis held
# at 1.1 of the unit: default and stated, drawn and held
numericArm <- function(column, drawn, sd = NULL, weights = NULL) {
  list(
    forests = list(
      list(sd = 1.1 * unit),
      list(basis = column, drawn = drawn, sd = sd)
    ),
    trees = 1L,
    weights = weights
  )
}
numericArms <- list(
  "a 0/1 column drawn at its default" = numericArm(binaryColumn, TRUE),
  "a 0/1 column held at its default" = numericArm(binaryColumn, FALSE),
  "a 0/1 column drawn at 0.6 of the unit" = numericArm(
    binaryColumn,
    TRUE,
    0.6 * unit
  ),
  "four values drawn at their default" = numericArm(fourValues, TRUE),
  "four values held at their default" = numericArm(fourValues, FALSE),
  "four values drawn at 3 of the unit" = numericArm(
    fourValues,
    TRUE,
    3 * unit
  ),
  "four values held at 3 of the unit" = numericArm(
    fourValues,
    FALSE,
    3 * unit
  ),
  "rows of weight zero, drawn at the default" = numericArm(
    weightedColumn,
    TRUE,
    weights = caseWeights
  ),
  "rows of weight zero, held at the default" = numericArm(
    weightedColumn,
    FALSE,
    weights = caseWeights
  ),
  "three forests of three trees, the number drawn" = list(
    forests = list(
      list(),
      list(basis = fourValues, drawn = TRUE),
      list(basis = z2, levels = TRUE, sd = 0.4 * unit)
    ),
    trees = 3L
  ),
  "three forests of three trees, the number held" = list(
    forests = list(
      list(),
      list(basis = fourValues),
      list(basis = z2, levels = TRUE, sd = 0.4 * unit)
    ),
    trees = 3L
  )
)
for (name in names(numericArms)) {
  arm <- numericArms[[name]]
  anyFailure <- runArm(
    name,
    arm,
    exact = numericOracle(arm),
    fit = numericFit
  ) ||
    anyFailure
}
if (anyFailure) {
  quit(status = 1L)
}
cat("OK: every arm matches the exact posterior\n")
