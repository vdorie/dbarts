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
# The Student-t arms (amplitude = student(df)) follow: a prior arm, with
# sigma held so large that the likelihood is flat, where a drawn
# coefficient's size is within the closed form's quantiles of its prior with
# the probability each names, on a numeric column, a factor's difference and a
# forest with no basis; the drawn numeric arm under student(1) and
# student(3), its oracle's grid on the prior's quantile scale; and a drawn
# factor at its default and under student(3), a two-dimensional quadrature
# over its pair. Every oracle's grid is checked against one twice as wide
# (twice as fine on a quantile scale): its functionals agree to 1e-4.
#
# Not seen here: probit and logistic, subset, a numeric multiplier under
# logistic, and n for n - 1 in a column's scale on a held arm, which the
# scale's literals in tests/cpp and the tinytest record's 1e-12 hold.
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
thresholds <- c(0.5, 1, 2)
normalQuartile <- qnorm(0.75)

# A drawn numeric coefficient's grid under its prior: `values`, their log
# prior weights and, for each value, the share of the mass it stands for
# within each threshold. The normal's is uniform, -8 to 8 by 0.005, times
# `wide`, each point standing for one step around it; a t's, prior =
# list(df =), is the midpoints of `points` times `wide` equal shares of
# probability, a = c qt(u, df) with c the median match's scale, z_.75 /
# t_df,.75, each point standing for its share.
coefficientGrid <- function(prior, wide, points = 4000L) {
  if (is.null(prior)) {
    values <- seq(-8 * wide, 8 * wide, by = 0.005)
    step <- values[2L] - values[1L]
    within <- vapply(
      thresholds,
      function(t) pmin(1, pmax(0, (t - abs(values)) / step + 0.5)),
      values
    )
    return(list(
      values = values,
      logPrior = dnorm(values, log = TRUE),
      within = within
    ))
  }
  size <- points * wide
  scale <- normalQuartile / qt(0.75, prior$df)
  values <- scale * qt((seq_len(size) - 0.5) / size, prior$df)
  within <- vapply(
    thresholds,
    function(t) {
      inside <- pt(c(-t, t) / scale, prior$df) * size
      pmin(
        1,
        pmax(
          0,
          pmin(seq_len(size), inside[2L]) -
            pmax(seq_len(size) - 1L, inside[1L])
        )
      )
    },
    values
  )
  list(values = values, logPrior = rep(0, size), within = within)
}

# forests: a list of list(multiplier = <n-vector>, drawn = TRUE where the
# coefficient is drawn, prior = NULL for the normal or list(df =), trees, sd
# = the prior sd of the forest's own fit). rows: the rows in the likelihood.
# Returns the posterior mean of the fit in each class of rows that share a
# cell and every multiplier, the probability that each forest's tree splits,
# and where a coefficient is drawn the probability that it is within 0.5, 1
# and 2, on the coefficient grid `wide` times the base one.
numericPosterior <- function(forests, shift, rows, classes, wide = 1) {
  K <- length(forests)
  counts <- as.matrix(expand.grid(lapply(forests, function(forest) {
    0:forest$trees
  })))
  residual <- (y - shift)[rows]
  drawn <- which(vapply(forests, `[[`, NA, "drawn"))
  grid <- if (length(drawn) > 0L) {
    coefficientGrid(forests[[drawn]]$prior, wide)
  } else {
    list(values = 1, logPrior = 0)
  }
  values <- grid$values
  logPriorA <- grid$logPrior
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
    # each grid point stands for its share of the coefficient, part of which
    # may lie within a bound
    within <- colSums(rowSums(weights) * grid$within)
    names(within) <- paste0("within", thresholds)
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
  } else if (!is.null(spec$prior)) {
    arguments$amplitude <- dbarts::dbartsPriors$student(spec$prior$df)
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
    last <- sum(widths[seq_len(drawn)])
    # a factor's coefficient is the difference of its pair
    coefficient <- if (widths[[drawn]] == 2L) {
      result$glue[last, ] - result$glue[last - 1L, ]
    } else {
      result$glue[last, ]
    }
    within <- vapply(thresholds, function(t) mean(abs(coefficient) < t), 0)
    names(within) <- paste0("within", thresholds)
    out <- c(out, within)
  }
  out
}

# the oracle's forests: a stated sd as written, a default by the rule, from
# the response drawn and R's sd() over the rows in the likelihood; `posterior`
# integrates them on its grid and on one twice as wide, whose functionals
# must agree to 1e-4
numericOracle <- function(arm, posterior = numericPosterior) {
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
      prior = spec$prior,
      trees = arm$trees,
      sd = if (drawn) size / 0.674 else size
    )
  })
  oracle <- posterior(forests, shift, rows, rowClasses(arm))
  wider <- posterior(forests, shift, rows, rowClasses(arm), wide = 2)
  if (max(abs(oracle - wider[names(oracle)])) > 1e-4) {
    stop("an oracle's grid and one twice as wide disagree past 1e-4")
  }
  oracle
}

# A drawn factor beside a held forest with no basis, one tree each: given the
# pair (b0, b1) and how many trees split, each class's fit, the forest with
# no basis plus b_z times the factor's forest in its cell, is jointly normal,
# and the rows reduce to their classes' means. The pair is integrated on a
# grid of the midpoints of `points` times `wide` equal shares of each of its
# marginal priors, N(0, 1/2) or t_df with scale sqrt(1/2) z_.75 / t_df,.75,
# each point weighted by the pair's joint density over the product of the
# marginals (one variance for the pair, so a bivariate t); the share of each
# point's b1 within t of its b0 is taken on b1's scale. The grid is summed a
# block of b1 at a time against a running maximum of the log weights.
factorPosterior <- function(
  forests,
  shift,
  rows,
  classes,
  wide = 1,
  points = 600L
) {
  size <- points * wide
  prior <- forests[[2L]]$prior
  df <- if (is.null(prior)) Inf else prior$df
  scale <- sqrt(0.5) *
    if (is.null(prior)) 1 else normalQuartile / qt(0.75, df)
  margin <- scale * qt((seq_len(size) - 0.5) / size, df)
  first <- vapply(classes, `[`, 0L, 1L)
  classCell <- cell[first]
  classLevel <- forests[[2L]]$multiplier[first]
  noise <- sigma^2 / lengths(classes)
  means <- vapply(classes, function(at) mean((y - shift)[at]), 0)
  m <- length(classes)
  leaf <- vapply(forests, function(forest) forest$sd^2, 0)
  configs <- as.matrix(expand.grid(0:1, 0:1))
  # the running sums: the largest log weight seen and, scaled to it, the
  # total weight, each class's fit, each forest's split and each threshold
  top <- -Inf
  sums <- numeric(1L + m + 2L + length(thresholds))
  for (block in split(seq_len(size), ceiling(seq_len(size) / 100L))) {
    b0 <- rep(margin, times = length(block))
    b1 <- rep(margin[block], each = size)
    index <- rep(block, each = size)
    logRatio <- if (is.null(prior)) {
      0
    } else {
      lgamma((df + 2) / 2) -
        lgamma(df / 2) -
        log(df * pi * scale^2) -
        (df + 2) / 2 * log1p((b0^2 + b1^2) / (df * scale^2)) -
        dt(b0 / scale, df, log = TRUE) -
        dt(b1 / scale, df, log = TRUE) +
        2 * log(scale)
    }
    within <- vapply(
      thresholds,
      function(t) {
        lower <- pt((b0 - t) / scale, df) * size
        upper <- pt((b0 + t) / scale, df) * size
        pmin(1, pmax(0, pmin(index, upper) - pmax(index - 1L, lower)))
      },
      b0
    )
    multiplier <- lapply(classLevel, function(z) if (z == 1) b1 else b0)
    for (r in seq_len(nrow(configs))) {
      k <- configs[r, ]
      # the cell values' covariance within each forest at k of its tree's
      # splits: 1 on the diagonal and 1 - k across the cells
      across <- function(f, g, h) {
        leaf[f] * if (classCell[g] == classCell[h]) 1 else 1 - k[f]
      }
      C <- matrix(list(), m, m)
      for (g in seq_len(m)) {
        for (h in seq_len(g)) {
          C[[g, h]] <- across(1L, g, h) +
            multiplier[[g]] * multiplier[[h]] * across(2L, g, h)
          C[[h, g]] <- C[[g, h]]
        }
      }
      # the class means' covariance, C plus the noise, factored at every
      # point: its log determinant and its inverse times the means
      L <- matrix(list(), m, m)
      logDet <- 0
      for (j in seq_len(m)) {
        d <- C[[j, j]] + noise[j]
        for (l in seq_len(j - 1L)) {
          d <- d - L[[j, l]]^2
        }
        L[[j, j]] <- sqrt(d)
        logDet <- logDet + 2 * log(L[[j, j]])
        for (i in seq_len(m)[seq_len(m) > j]) {
          e <- C[[i, j]]
          for (l in seq_len(j - 1L)) {
            e <- e - L[[i, l]] * L[[j, l]]
          }
          L[[i, j]] <- e / L[[j, j]]
        }
      }
      forward <- vector("list", m)
      for (i in seq_len(m)) {
        e <- means[i]
        for (l in seq_len(i - 1L)) {
          e <- e - L[[i, l]] * forward[[l]]
        }
        forward[[i]] <- e / L[[i, i]]
      }
      alpha <- vector("list", m)
      quadratic <- 0
      for (i in rev(seq_len(m))) {
        e <- forward[[i]]
        for (l in seq_len(m)[seq_len(m) > i]) {
          e <- e - L[[l, i]] * alpha[[l]]
        }
        alpha[[i]] <- e / L[[i, i]]
        quadratic <- quadratic + means[i] * alpha[[i]]
      }
      logWeights <- sum(dbinom(k, 1L, base, log = TRUE)) +
        logRatio -
        0.5 * logDet -
        0.5 * quadratic
      if (max(logWeights) > top) {
        sums <- sums * exp(top - max(logWeights))
        top <- max(logWeights)
      }
      weights <- exp(logWeights - top)
      fits <- vapply(
        seq_len(m),
        function(g) {
          fit <- 0
          for (h in seq_len(m)) {
            fit <- fit + C[[g, h]] * alpha[[h]]
          }
          sum(weights * fit)
        },
        0
      )
      total <- sum(weights)
      sums <- sums +
        c(total, fits, total * k, colSums(weights * within))
    }
  }
  out <- sums[-1L] / sums[1L]
  names(out) <- c(
    names(classes),
    paste0("split", 1:2),
    paste0("within", thresholds)
  )
  out
}

# The prior arm's sampler: sigma held at 1e6, so the likelihood is flat to
# 1e-9 and each drawn coefficient is a draw from its prior. `read` names the
# glue rows of each coefficient whose size is read, the difference of a
# factor's pair, and the prior median of that size; the functionals are the
# probability that the size is within the prior's quantile p, which is p.
priorProbabilities <- c(0.25, 0.5, 0.9)
priorQuantile <- function(p, prior, median) {
  df <- if (is.null(prior)) Inf else prior$df
  median * qt((1 + p) / 2, df) / qt(0.75, df)
}
priorFit <- function(seed, arm, kept) {
  set.seed(seed)
  sampler <- dbarts(
    x,
    y,
    family = gaussian(sigma = dbartsPriors$fixed(1e12)),
    forests = lapply(arm$forests, numericForest, trees = 1L),
    control = dbartsControl(
      n.chains = 1L,
      n.threads = 1L,
      n.cuts = 1L,
      updateState = FALSE
    )
  )
  glue <- sampler$run(numBurn, kept)$glue
  unlist(lapply(names(arm$read), function(name) {
    read <- arm$read[[name]]
    size <- abs(
      if (length(read$rows) == 2L) {
        glue[read$rows[2L], ] - glue[read$rows[1L], ]
      } else {
        glue[read$rows, ]
      }
    )
    out <- vapply(
      priorProbabilities,
      function(p) mean(size <= priorQuantile(p, read$prior, read$median)),
      0
    )
    names(out) <- paste0(name, "@", priorProbabilities)
    out
  }))
}
priorExact <- function(arm) {
  out <- rep(priorProbabilities, length(arm$read))
  names(out) <- paste0(
    rep(names(arm$read), each = length(priorProbabilities)),
    "@",
    priorProbabilities
  )
  out
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
numericArm <- function(
  column,
  drawn,
  sd = NULL,
  weights = NULL,
  prior = NULL
) {
  list(
    forests = list(
      list(sd = 1.1 * unit),
      list(basis = column, drawn = drawn, sd = sd, prior = prior)
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

# ---- the Student-t arms ----

studentT <- function(df) list(df = df)
read <- function(rows, prior, median = normalQuartile) {
  list(rows = rows, prior = prior, median = median)
}
coefficientArm <- function(prior) {
  list(
    forests = list(
      list(sd = 1.1 * unit),
      list(basis = fourValues, drawn = TRUE, prior = prior),
      list(basis = z2, levels = TRUE, drawn = TRUE, prior = prior)
    ),
    read = list(number = read(2L, prior), difference = read(3:4, prior))
  )
}
noBasisArm <- function(prior, sd = NULL) {
  list(
    forests = list(
      list(drawn = TRUE, prior = prior, sd = sd),
      list(basis = z1, levels = TRUE)
    ),
    read = list(
      coefficient = read(1L, prior, if (is.null(sd)) 2 else sd / unit)
    )
  )
}
priorArms <- list(
  "prior: a number and a factor at the default" = coefficientArm(NULL),
  "prior: a number and a factor under student(1)" = coefficientArm(
    studentT(1)
  ),
  "prior: a number and a factor under student(3)" = coefficientArm(
    studentT(3)
  ),
  "prior: no basis under student(1) at its default" = noBasisArm(studentT(1)),
  "prior: no basis under student(3) at 0.7 of the unit" = noBasisArm(
    studentT(3),
    0.7 * unit
  )
)
for (name in names(priorArms)) {
  arm <- priorArms[[name]]
  anyFailure <- runArm(name, arm, exact = priorExact(arm), fit = priorFit) ||
    anyFailure
}

factorArm <- function(prior) {
  list(
    forests = list(
      list(sd = 1.1 * unit),
      list(basis = z1, levels = TRUE, drawn = TRUE, prior = prior)
    ),
    trees = 1L
  )
}
studentArms <- list(
  "four values drawn under student(1)" = numericArm(
    fourValues,
    TRUE,
    prior = studentT(1)
  ),
  "four values drawn under student(3)" = numericArm(
    fourValues,
    TRUE,
    prior = studentT(3)
  ),
  "a factor drawn at its default" = factorArm(NULL),
  "a factor drawn under student(3)" = factorArm(studentT(3))
)
for (name in names(studentArms)) {
  arm <- studentArms[[name]]
  levels <- isTRUE(arm$forests[[2L]]$levels)
  anyFailure <- runArm(
    name,
    arm,
    exact = numericOracle(
      arm,
      if (levels) factorPosterior else numericPosterior
    ),
    fit = numericFit
  ) ||
    anyFailure
}
if (anyFailure) {
  quit(status = 1L)
}
cat("OK: every arm matches the exact posterior\n")
