#!/usr/bin/env Rscript

# The negative binomial default-prior study: is there a better default prior
# on the forest's spread for family = "nbinom" fits? The design, the criterion
# and the decision rule are fixed in docs/plans/nbinom-default-study.md; this
# file implements its stages 0 and 1 (the checks and the prior predictive)
# and decides nothing the plan leaves open except where its comments say so.
#
# BUILDS. One R process sees one build: whatever dbarts R_LIBS finds first,
# named by NBS_BUILD ("today", the shipped sampler; "move", the scratch build
# carrying the collapsed scale move's negative binomial half; "mutant", the
# move build with one power of r too many in the move's Jacobian, for the
# drift test's power). Actions that the plan assigns to a build refuse to run
# under another.
#
# RUN DIRECTORY. NBS_RUN holds every output: data/ (the exported real
# datasets and their truths), checks/, prior/ and the summaries. Nothing is
# read from the working directory.
#
# Usage (each action one process):
#   Rscript nbinom-default-study.R export          (data and additive truths)
#   Rscript nbinom-default-study.R truth bart <dataset>
#   Rscript nbinom-default-study.R drift <seed> <reps>
#   Rscript nbinom-default-study.R sbc <first> <last>
#   Rscript nbinom-default-study.R compare <SBC dataset> <seed>
#   Rscript nbinom-default-study.R reweight <prior1..3 | null1..2>
#   Rscript nbinom-default-study.R prior-checks
#   Rscript nbinom-default-study.R prior <cell>
#   Rscript nbinom-default-study.R summarize <stage>
# where <stage> is drift, sbc, compare, reweight or prior.
#
# BLIND. Nothing here scores a fit against a truth: the Part B scoring is
# not written until the pre-registration commit, and no run length, weight
# or acceptance below reads a truth. The stage 0 checks compare the sampler
# with itself and with the prior, which is not scoring.

suppressPackageStartupMessages(library(dbarts))

## ---------------------------------------------------------------- settings

`%not_in_set%` <- function(x, table) !(x %in% table)

runDir <- Sys.getenv("NBS_RUN")
if (!nzchar(runDir)) {
  stop("NBS_RUN must name the run directory")
}
build <- Sys.getenv("NBS_BUILD")
if (build %not_in_set% c("today", "move", "mutant")) {
  stop("NBS_BUILD must be today, move or mutant")
}

requireBuild <- function(builds) {
  if (build %not_in_set% builds) {
    stop(
      "this action runs under build ",
      paste(builds, collapse = " or "),
      ", not ",
      build
    )
  }
}

outPath <- function(...) {
  path <- file.path(runDir, ...)
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  path
}

# The engine's shape grid and its gamma(2, 0.1) weights (NBShapePrior), and
# C6's log-uniform weights, fixed by the plan.
shapeGrid <- c(1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 30, 50)
shapeWeights <- local({
  w <- shapeGrid * exp(-0.1 * shapeGrid)
  w / sum(w)
})
shapeWeightsC6 <- c(
  0.0886,
  0.1404,
  0.0886,
  0.0653,
  0.0518,
  0.0601,
  0.0653,
  0.0518,
  0.0518,
  0.0653,
  0.0886,
  0.1171,
  0.0653
)
drawShape <- function(weights = shapeWeights) {
  shapeGrid[sample.int(length(shapeGrid), 1L, prob = weights)]
}

# The leaf scale A of s = A / k; checked against the sampler by prior-checks.
leafScale <- 3
countCap <- 1e6
numTrees <- 75L

# Neighboring columns correlated 0.5, as the designed cells' covariates.
correlatedNormals <- function(n, p) {
  sigma <- 0.5^abs(outer(seq_len(p), seq_len(p), "-"))
  x <- matrix(rnorm(n * p), n, p) %*% chol(sigma)
  colnames(x) <- paste0("x", seq_len(p))
  x
}

samplerControl <- function(numChains = 1L, numSamples = 1L) {
  dbartsControl(
    n.chains = numChains,
    n.threads = 1L,
    n.samples = numSamples,
    n.burn = 0L,
    n.trees = numTrees,
    updateState = FALSE,
    verbose = FALSE
  )
}

# Trees, k0 and leaves from the prior, with k0 drawn as the default's
# chi(1.5, 2); the sampler then holds k0 as its state. Returns k0 and the
# log mean without the offset, c + f0.
drawFromPrior <- function(
  sampler,
  x,
  drawK = function() 2 * sqrt(rchisq(1L, 1.5))
) {
  sampler$sampleTreesFromPrior()
  k0 <- drawK()
  sampler$setLeafPrior(dbartsPriors$normal(k = k0))
  sampler$sampleLeafParametersFromPrior()
  list(k0 = k0, eta = as.numeric(sampler$predict(x)))
}

# Counts given the log mean and offset; NULL when one passes the cap or
# their total passes totalCap. Either conditions the joint draw on an event
# of y alone, so the posterior given y, and so a check's exactness, is kept;
# the total's cap, above every real dataset's but Migration's, bounds a
# check's cost, which the Polya-Gamma draws make linear in the total.
totalCap <- 50000
drawCounts <- function(eta, offset, shape, cap = totalCap) {
  y <- suppressWarnings(rnbinom(
    length(eta),
    size = shape,
    mu = exp(eta + offset)
  ))
  if (anyNA(y) || max(y) > countCap || sum(y) > cap) NULL else y
}

## ------------------------------------------------- drift test from the truth

# The start is an exact posterior draw: trees, k0 and leaves from the prior,
# the shape r0 from the grid prior and held fixed (its own draw is shipped
# code), counts from the model at the sampler's center c (set from a
# placeholder response at construction and held by updateScale = FALSE), and
# the latents drawn given the counts by setResponse. log k averaged over lags
# 100 to 1,000 must not drift from log k0.
driftLags <- c(1, 10, 100, 200, 300, 400, 500, 700, 1000)
driftWindow <- driftLags >= 100

runDrift <- function(seed, numReps) {
  requireBuild(c("move", "mutant", "today"))
  set.seed(seed)
  n <- 80L
  results <- vector("list", numReps)
  started <- proc.time()[[3L]]
  for (r in seq_len(numReps)) {
    x <- correlatedNormals(n, 3L)
    offset <- runif(n, 0, log(100))
    shape <- drawShape()
    placeholder <- rpois(n, 0.1 * exp(offset))
    sampler <- dbarts(
      x,
      placeholder,
      offset = offset,
      family = nbinom(shape = shape),
      control = samplerControl()
    )
    defaultPrior <- sampler$getLeafPrior()$leaf.prior
    rejected <- 0L
    repeat {
      truth <- drawFromPrior(sampler, x)
      y <- drawCounts(truth$eta, offset, shape)
      if (!is.null(y)) {
        break
      }
      rejected <- rejected + 1L
    }
    sampler$setLeafPrior(defaultPrior)
    sampler$setResponse(y, updateScale = FALSE)
    k <- numeric(length(driftLags))
    previous <- 0
    for (i in seq_along(driftLags)) {
      fit <- sampler$run(as.integer(driftLags[i] - previous) - 1L, 1L)
      k[i] <- as.numeric(fit$k)
      previous <- driftLags[i]
    }
    results[[r]] <- list(
      k0 = truth$k0,
      k = k,
      shape = shape,
      shapeHeld = as.numeric(fit$shape) == shape,
      rejected = rejected,
      total = sum(y)
    )
  }
  seconds <- proc.time()[[3L]] - started
  cat(sprintf(
    "%s seed %d: %d replications, %.0f s\n",
    build,
    seed,
    numReps,
    seconds
  ))
  saveRDS(
    list(build = build, lags = driftLags, results = results, seconds = seconds),
    outPath("checks", sprintf("drift-%s-%d.rds", build, seed))
  )
}

summarizeDrift <- function() {
  files <- list.files(file.path(runDir, "checks"), "^drift-", full.names = TRUE)
  runs <- lapply(files, readRDS)
  for (b in unique(vapply(runs, `[[`, "", "build"))) {
    results <- do.call(
      c,
      lapply(Filter(function(run) run$build == b, runs), `[[`, "results")
    )
    k <- t(vapply(results, `[[`, numeric(length(driftLags)), "k"))
    k0 <- vapply(results, `[[`, 0, "k0")
    drift <- rowMeans(log(k[, driftWindow, drop = FALSE])) - log(k0)
    se <- sd(drift) / sqrt(length(drift))
    byLag <- colMeans(log(k)) - mean(log(k0))
    cat(sprintf(
      "%s: R = %d, drift of log k over lags 100 to 1,000 %+.4f (se %.4f), z %+.2f; shape held in %d of %d; %d rejected draws\n",
      b,
      length(drift),
      mean(drift),
      se,
      mean(drift) / se,
      sum(vapply(results, `[[`, NA, "shapeHeld")),
      length(results),
      sum(vapply(results, `[[`, 0L, "rejected"))
    ))
    cat(
      "  mean log k - mean log k0 by lag:",
      sprintf("%d %+.3f", driftLags, byLag),
      "\n"
    )
  }
}

## ------------------------------------------------- real datasets and truths

# Each loader returns the data frame, the response's column and the offset's
# (NULL when none). The covariates are every other column, as bart() sees
# them. All are in base R, recommended packages or CRAN packages installed on
# the laptop; export writes them to data/ for the box.
loadData <- function(name, package) {
  e <- new.env()
  data(list = name, package = package, envir = e)
  get(name, envir = e)
}
longForm <- function(species, env, columns) {
  d <- data.frame(
    env[rep(seq_len(nrow(env)), ncol(species)), columns, drop = FALSE],
    species = factor(rep(colnames(species), each = nrow(species))),
    y = as.vector(as.matrix(species))
  )
  rownames(d) <- NULL
  d
}
realLoaders <- list(
  roaches = function() {
    d <- loadData("roaches", "rstanarm")
    d$off <- log(d$exposure2)
    list(
      d = d[, c("y", "roach1", "treatment", "senior", "off")],
      y = "y",
      off = "off"
    )
  },
  solder = function() {
    list(d = loadData("solder", "survival"), y = "skips", off = NULL)
  },
  aids = function() {
    d <- loadData("aids", "boot")
    list(
      d = d[d$dud == 0, c("y", "delay", "time", "quarter")],
      y = "y",
      off = NULL
    )
  },
  mite = function() {
    d <- longForm(loadData("mite", "vegan"), loadData("mite.env", "vegan"), 1:5)
    list(d = d, y = "y", off = NULL)
  },
  Arabidopsis = function() {
    d <- loadData("Arabidopsis", "lme4")
    columns <- c(
      "total.fruits",
      "reg",
      "popu",
      "nutrient",
      "amd",
      "rack",
      "status"
    )
    list(d = d[, columns], y = "total.fruits", off = NULL)
  },
  grouseticks = function() {
    d <- loadData("grouseticks", "lme4")
    list(d = d[, c("TICKS", "YEAR", "HEIGHT")], y = "TICKS", off = NULL)
  },
  epil = function() {
    d <- loadData("epil", "MASS")
    list(d = d[, c("y", "trt", "base", "age", "V4")], y = "y", off = NULL)
  },
  Insurance = function() {
    d <- loadData("Insurance", "MASS")
    d$off <- log(d$Holders)
    for (v in c("District", "Group", "Age")) {
      d[[v]] <- factor(d[[v]], ordered = FALSE)
    }
    list(
      d = d[, c("Claims", "District", "Group", "Age", "off")],
      y = "Claims",
      off = "off"
    )
  },
  crimtab = function() {
    d <- as.data.frame(as.table(loadData("crimtab", "datasets")))
    names(d) <- c("finger", "height", "y")
    d$finger <- as.numeric(as.character(d$finger))
    d$height <- as.numeric(as.character(d$height))
    list(d = d, y = "y", off = NULL)
  },
  BCI = function() {
    columns <- c("UTM.EW", "UTM.NS", "EnvHet", "Habitat")
    d <- longForm(
      loadData("BCI", "vegan"),
      loadData("BCI.env", "vegan"),
      columns
    )
    list(d = d, y = "y", off = NULL)
  },
  Traffic = function() {
    d <- loadData("Traffic", "MASS")
    list(d = d[, c("y", "limit", "year", "day")], y = "y", off = NULL)
  },
  quakes = function() {
    list(d = loadData("quakes", "datasets"), y = "stations", off = NULL)
  },
  Migration = function() {
    d <- loadData("Migration", "carData")
    d <- data.frame(
      migrants = d$migrants,
      lp66 = log(d$pops66),
      ld66 = log(d$popd66),
      ldist = log(d$distance)
    )
    list(d = d, y = "migrants", off = NULL)
  }
)
# Fitted only for the reported percentiles (Part A); no Part B cell.
reportedOnly <- "Migration"

# The four fits of the plan's rule: a quadratic NB GLM; the additive NB GAM
# (REML s() of each numeric covariate with more than 10 distinct values, k =
# min(10, distinct values - 1), the rest as factors, no tensor or by-terms);
# and, with two or more such covariates, s() and te() of the first two plus
# additive smooths of the rest. The additive fit is the additive truth. The
# Poisson GLM on main effects, offset in the fit, gives C7's spread.
exportData <- function() {
  requireBuild(c("today", "move"))
  suppressPackageStartupMessages({
    library(MASS)
    library(mgcv)
  })
  spreads <- NULL
  for (name in names(realLoaders)) {
    started <- proc.time()[[3L]]
    s <- realLoaders[[name]]()
    d <- s$d
    offset <- if (is.null(s$off)) rep(0, nrow(d)) else d[[s$off]]
    covariates <- setdiff(names(d), c(s$y, s$off))
    numeric <- covariates[vapply(
      covariates,
      function(v) is.numeric(d[[v]]) && length(unique(d[[v]])) > 10L,
      NA
    )]
    factors <- setdiff(covariates, numeric)
    for (v in factors) {
      if (!is.factor(d[[v]])) d[[v]] <- factor(d[[v]])
    }
    offsetTerm <- if (is.null(s$off)) {
      character()
    } else {
      sprintf("offset(%s)", s$off)
    }
    rhs <- function(terms) {
      paste(c(terms, factors, offsetTerm), collapse = " + ")
    }
    basisSize <- function(v) min(10L, length(unique(d[[v]])) - 1L)
    smooth <- function(v) sprintf("s(%s, k = %d)", v, vapply(v, basisSize, 1L))
    fits <- list(
      quadratic = list(glm = TRUE, f = rhs(sprintf("poly(%s, 2)", numeric))),
      additive = list(glm = FALSE, f = rhs(smooth(numeric)))
    )
    if (length(numeric) >= 2L) {
      rest <- if (length(numeric) > 2L) smooth(numeric[-(1:2)]) else character()
      pair <- paste(numeric[1L], numeric[2L], sep = ", ")
      fits$s2 <- list(glm = FALSE, f = rhs(c(sprintf("s(%s)", pair), rest)))
      fits$te <- list(glm = FALSE, f = rhs(c(sprintf("te(%s)", pair), rest)))
    }
    truth <- NULL
    for (fitName in names(fits)) {
      formula <- as.formula(paste(s$y, "~", fits[[fitName]]$f))
      fit <- tryCatch(
        suppressWarnings(
          if (fits[[fitName]]$glm) {
            glm.nb(formula, data = d)
          } else {
            gam(formula, data = d, family = nb(), method = "REML")
          }
        ),
        error = function(e) NULL
      )
      if (is.null(fit)) {
        next
      }
      eta <- as.numeric(predict(fit, type = "link")) - offset
      theta <- if (fits[[fitName]]$glm) fit$theta else fit$family$getTheta(TRUE)
      identified <- exp(eta + offset) >= 0.01
      spreads <- rbind(
        spreads,
        data.frame(
          data = name,
          fit = fitName,
          S = sd(eta),
          S.id = sd(eta[identified]),
          share.below = mean(!identified),
          theta = theta
        )
      )
      if (fitName == "additive") truth <- list(eta = eta, theta = theta)
    }
    poissonFormula <- as.formula(paste(s$y, "~", rhs(covariates)))
    poisson <- suppressWarnings(glm(
      poissonFormula,
      data = d,
      family = poisson()
    ))
    y <- d[[s$y]]
    saveRDS(
      list(
        name = name,
        x = d[, covariates, drop = FALSE],
        y = y,
        offset = if (is.null(s$off)) NULL else offset,
        meanCount = sum(y) / sum(exp(offset)),
        glmSpread = sd(as.numeric(predict(poisson, type = "link")) - offset),
        truthAdditive = truth,
        reportedOnly = name %in% reportedOnly
      ),
      outPath("data", paste0(name, ".rds"))
    )
    cat(sprintf("%s: %.0f s\n", name, proc.time()[[3L]] - started))
  }
  saveRDS(spreads, outPath("data", "spreads.rds"))
  print(spreads, digits = 3L, row.names = FALSE)
}

readData <- function(name) {
  readRDS(file.path(runDir, "data", paste0(name, ".rds")))
}

## ------------------------------------------------------- the spread priors

# Each candidate prior on s = A / k: draw(m) returns m spreads; above(t) is
# P(s > t), which reads A1 without Monte Carlo error in s.
chiSpread <- function(df, scale) {
  list(
    draw = function(m) leafScale / (scale * sqrt(rchisq(m, df))),
    above = function(t) pchisq((leafScale / (t * scale))^2, df)
  )
}
fixedSpread <- function(k) {
  list(
    draw = function(m) rep(leafScale / k, m),
    above = function(t) as.numeric(leafScale / k > t)
  )
}
# chi(1.5, 2) with k at least kFloor, by inversion as the plan's C4 draws it.
flooredSpread <- function(kFloor, df = 1.5, scale = 2) {
  atFloor <- pchisq((kFloor / scale)^2, df)
  list(
    draw = function(m) {
      leafScale / (scale * sqrt(qchisq(atFloor + (1 - atFloor) * runif(m), df)))
    },
    above = function(t) {
      pmax(pchisq((leafScale / (t * scale))^2, df) - atFloor, 0) / (1 - atFloor)
    }
  )
}
halfNormalSpread <- function(median) {
  sigma <- median / qnorm(0.75)
  list(
    draw = function(m) abs(rnorm(m, 0, sigma)),
    above = function(t) 2 * pnorm(-t / sigma)
  )
}
# C7: chi(4, c) with c putting s's median at 1.7 times the Poisson GLM's
# spread, floored at 0.25.
dataScaledSpread <- function(glmSpread) {
  target <- 1.7 * max(glmSpread, 0.25)
  chiSpread(4, leafScale / (target * sqrt(qchisq(0.5, 4))))
}
spreadPriors <- list(
  C0 = chiSpread(1.5, 2),
  `C1 k=1.5` = fixedSpread(1.5),
  `C1 k=2` = fixedSpread(2),
  `C2 chi(3,1.25)` = chiSpread(3, 1.25),
  `C2 chi(4,1.03)` = chiSpread(4, 1.03),
  `C2 chi(5,0.9)` = chiSpread(5, 0.9),
  `C2 chi(7,0.75)` = chiSpread(7, 0.75),
  `C2 chi(10,0.62)` = chiSpread(10, 0.62),
  `C3 chi(1.5,3)` = chiSpread(1.5, 3),
  `C4 k>=0.6` = flooredSpread(0.6),
  `C4 k>=0.75` = flooredSpread(0.75),
  `C4 k>=1` = flooredSpread(1),
  `C5 hn 1.3` = halfNormalSpread(1.3),
  `C5 hn 1.57` = halfNormalSpread(1.57)
)

## ----------------------------------------------- stage 1: prior predictive

# Cells: the simulated grid's designs at each n (the mean counts are read
# within a cell, since the forest does not depend on them) and the real
# designs at their own mean count.
priorCellNames <- c(
  as.vector(outer(
    c("norm5", "unif20", "mix10off"),
    c(50L, 500L, 2000L, 10000L),
    paste,
    sep = "-"
  )),
  names(realLoaders)
)
simulatedMeans <- c(0.1, 1, 10, 100)
# The bound of A1, A2 and A3, and the threshold of A1.
partABound <- 0.01
spreadThreshold <- 4
identifiedMean <- 0.01

priorCell <- function(cell) {
  index <- match(cell, priorCellNames)
  if (is.na(index)) {
    stop("unknown cell ", cell)
  }
  if (cell %in% names(realLoaders)) {
    data <- readData(cell)
    return(list(
      x = data$x,
      offset = if (is.null(data$offset)) {
        rep(0, length(data$y))
      } else {
        data$offset
      },
      means = data$meanCount,
      build = data$y,
      dataScaled = setNames(data$glmSpread, "own"),
      seed = index
    ))
  }
  parts <- strsplit(cell, "-", fixed = TRUE)[[1L]]
  n <- as.integer(parts[2L])
  set.seed(1000L + index)
  x <- switch(
    parts[1L],
    norm5 = correlatedNormals(n, 5L),
    unif20 = matrix(runif(n * 20L), n, 20L),
    mix10off = cbind(
      correlatedNormals(n, 5L),
      matrix(rbinom(n * 5L, 1L, 0.3), n, 5L)
    )
  )
  # Exposure log-uniform over three orders of magnitude, scaled to mean 1 so
  # that the cell's label is the mean count a row, as in the other designs.
  offset <- if (parts[1L] == "mix10off") {
    exposure <- exp(runif(n, 0, log(1000)))
    log(exposure / mean(exposure))
  } else {
    rep(0, n)
  }
  list(
    x = x,
    offset = offset,
    means = simulatedMeans,
    build = rep(1:2, length.out = n),
    dataScaled = c(`0.25` = 0.25, `1` = 1, `2` = 2),
    seed = index
  )
}

# P(some count above the cap) given rows' log means and the shape, from each
# row's upper tail (the log lower tail warns of underflow near 1).
capProbability <- function(logMean, shape) {
  above <- pnbinom(countCap, shape, mu = exp(logMean), lower.tail = FALSE)
  -expm1(sum(log1p(-above)))
}

logSumExp <- function(v) {
  m <- max(v)
  m + log(sum(exp(v - m)))
}

# One forest at unit spread, u, scored under one spread s at each mean count:
# S, S_id, the probability of a count above the cap (averaged over the shape
# grid under the default's and C6's weights, the counts integrated out), and
# whether the expected share of zeros passes 0.9 at one shape draw.
scoreSpread <- function(u, offset, logExposure, s, means, shape) {
  v <- offset + s * u
  shift <- logExposure - logSumExp(v)
  out <- matrix(NA_real_, 5L, length(means))
  for (j in seq_along(means)) {
    logMean <- v + shift + log(means[j])
    identified <- logMean >= log(identifiedMean)
    top <- logMean[logMean > log(1e3)]
    capped <- if (length(top) == 0L) {
      c(0, 0)
    } else {
      p <- vapply(
        shapeGrid,
        function(r) {
          capProbability(top, r)
        },
        0
      )
      c(sum(shapeWeights * p), sum(shapeWeightsC6 * p))
    }
    zeroShare <- mean(exp(shape * (log(shape) - log(shape + exp(logMean)))))
    out[, j] <- c(
      if (sum(identified) > 1L) s * sd(u[identified]) else NA_real_,
      capped,
      zeroShare > 0.9,
      sum(!identified)
    )
  }
  out
}

priorBlockSize <- 2500L
priorMinimum <- 5000L
priorMaximum <- 20000L

# Draws forests in blocks until every candidate's A1, A2 and A3 verdicts are
# resolved: its estimate at least two standard errors from the bound, or its
# standard error under 5 percent of the bound (a candidate that close sits at
# the bound, and more forests would not move a ruling). At least 5,000 forests
# (50 above a 99th percentile), at most 20,000 (the plan's count).
runPrior <- function(cell) {
  requireBuild(c("today", "move"))
  config <- priorCell(cell)
  n <- length(config$offset)
  sampler <- suppressWarnings(dbarts(
    config$x,
    config$build,
    offset = if (any(config$offset != 0)) config$offset,
    family = "nbinom",
    control = samplerControl()
  ))
  sampler$setLeafPrior(dbartsPriors$normal(k = leafScale))
  shift <- sampler$getLeafPrior()$response.shift
  priors <- spreadPriors
  for (label in names(config$dataScaled)) {
    priors[[paste0("C7 ", label)]] <- dataScaledSpread(config$dataScaled[[
      label
    ]])
  }
  logExposure <- log(sum(exp(config$offset)))
  realTop <- if (cell %in% names(realLoaders)) {
    spreads <- readRDS(file.path(runDir, "data", "spreads.rds"))
    max(spreads$S.id[spreads$data == cell])
  }
  blocks <- list()
  started <- proc.time()[[3L]]
  repeat {
    block <- length(blocks) + 1L
    set.seed(config$seed * 1000L + block)
    m <- priorBlockSize
    unitS <- numeric(m)
    above <- matrix(
      NA_real_,
      m,
      length(priors),
      dimnames = list(NULL, names(priors))
    )
    scored <- array(
      NA_real_,
      c(m, length(priors), 6L, length(config$means)),
      dimnames = list(
        NULL,
        names(priors),
        c("S", "S.id", "capped", "capped.C6", "zeros90", "unidentified"),
        format(config$means)
      )
    )
    for (i in seq_len(m)) {
      sampler$sampleTreesFromPrior()
      sampler$sampleLeafParametersFromPrior()
      u <- as.numeric(sampler$predict(config$x)) - shift
      unitS[i] <- sd(u)
      shape <- drawShape()
      for (p in seq_along(priors)) {
        s <- priors[[p]]$draw(1L)
        above[i, p] <- priors[[p]]$above(spreadThreshold / unitS[i])
        scored[i, p, 1L, ] <- s * unitS[i]
        scored[i, p, -1L, ] <- scoreSpread(
          u,
          config$offset,
          logExposure,
          s,
          config$means,
          shape
        )
      }
    }
    blocks[[block]] <- list(unitS = unitS, above = above, scored = scored)
    total <- block * m
    verdicts <- priorVerdicts(blocks, realTop)
    if (
      total >= priorMaximum || (total >= priorMinimum && all(verdicts$resolved))
    ) {
      break
    }
  }
  seconds <- proc.time()[[3L]] - started
  cat(sprintf("%s: %d forests, %.0f s\n", cell, total, seconds))
  saveRDS(
    list(
      cell = cell,
      n = n,
      means = config$means,
      realTop = realTop,
      dataScaled = config$dataScaled,
      blocks = blocks,
      seconds = seconds
    ),
    outPath("prior", paste0(cell, ".rds"))
  )
}

# Estimates and standard errors of A1 (Rao-Blackwellized in s), A2 (in the
# counts and the shape) and, at a real design, A3 as the share of datasets
# whose S_id passes the real data's widest (A3 passes when that share is at
# least the bound). Resolved as runPrior states.
priorVerdicts <- function(blocks, realTop) {
  above <- do.call(rbind, lapply(blocks, `[[`, "above"))
  scored <- do.call(abind3, lapply(blocks, `[[`, "scored"))
  m <- nrow(above)
  estimate <- function(values) {
    c(mean(values, na.rm = TRUE), sd(values, na.rm = TRUE) / sqrt(m))
  }
  rows <- list()
  for (p in colnames(above)) {
    for (j in seq_len(dim(scored)[4L])) {
      rows[[length(rows) + 1L]] <- data.frame(
        prior = p,
        mean = dimnames(scored)[[4L]][j],
        part = c("A1", "A2", if (!is.null(realTop)) "A3"),
        rbind(
          estimate(above[, p]),
          estimate(scored[, p, "capped", j]),
          if (!is.null(realTop)) {
            estimate(as.numeric(scored[, p, "S.id", j] > realTop))
          }
        )
      )
    }
  }
  result <- do.call(rbind, rows)
  names(result)[4:5] <- c("estimate", "se")
  result$resolved <- abs(result$estimate - partABound) >= 2 * result$se |
    result$se <= 0.05 * partABound
  result
}

# rbind along the first dimension of four-dimensional arrays.
abind3 <- function(...) {
  arrays <- list(...)
  d <- dim(arrays[[1L]])
  out <- array(
    NA_real_,
    c(sum(vapply(arrays, function(a) dim(a)[1L], 0L)), d[-1L]),
    dimnames = c(list(NULL), dimnames(arrays[[1L]])[-1L])
  )
  start <- 0L
  for (a in arrays) {
    out[start + seq_len(dim(a)[1L]), , , ] <- a
    start <- start + dim(a)[1L]
  }
  out
}

## ------------------------------------------------------- stage 0: SBC of k

# Geyer's initial monotone sequence estimate of the integrated autocorrelation
# time, in units of the series' own spacing (as probit-k-mixing.R).
iatGeyer <- function(x) {
  n <- length(x)
  x <- x - mean(x)
  if (!(var(x) > 0)) {
    return(NA_real_)
  }
  m <- 2^ceiling(log2(2 * n))
  fx <- fft(c(x, numeric(m - n)))
  ac <- Re(fft(Mod(fx)^2, inverse = TRUE))[seq_len(n)]
  ac <- ac / ac[1L]
  numPairs <- floor((n - 1) / 2)
  pairSums <- ac[2 * (0:(numPairs - 1)) + 1] + ac[2 * (0:(numPairs - 1)) + 2]
  firstNonPositive <- which(pairSums <= 0)[1L]
  if (!is.na(firstNonPositive)) {
    pairSums <- pairSums[seq_len(firstNonPositive - 1L)]
  }
  if (length(pairSums) == 0L) {
    return(1)
  }
  max(-1 + 2 * sum(cummin(pairSums)), 1 / n)
}

# The SBC design: 50 rows of the drift test's covariates and exposure, fixed
# across datasets; the center set by a placeholder and held, so that the
# counts are drawn at the fitted model's own c.
sbcDesign <- function() {
  set.seed(4242L)
  n <- 50L
  x <- correlatedNormals(n, 3L)
  offset <- runif(n, 0, log(100))
  list(x = x, offset = offset, placeholder = rpois(n, 0.1 * exp(offset)))
}
sbcSampler <- function(design, numChains = 1L) {
  dbarts(
    design$x,
    design$placeholder,
    offset = design$offset,
    family = "nbinom",
    control = samplerControl(numChains)
  )
}
sbcNumRanks <- 99L
sbcBurn <- 2000L
sbcChunk <- 10000L
sbcCap <- 100000L
# The SBC's datasets are held to a total count of 10,000, which keeps the
# check exact (a condition on y) and its slowest datasets at minutes, not
# the quarter hours a total of 50,000 with slow trees costs at the cap.
sbcTotalCap <- 10000

# Per dataset: the truth (trees, k0, leaves, a shape r0 from the grid prior)
# and counts; a start independent of it (fresh prior trees and k); 2,000
# sweeps of burn-in; then chunks of 10,000 recorded sweeps until 99 draws
# spaced at twice k's autocorrelation time fit, capped at 100,000. k0's rank
# among the 99, and r0's with ties broken at random.
runSbc <- function(first, last) {
  requireBuild("move")
  design <- sbcDesign()
  sampler <- sbcSampler(design)
  defaultPrior <- sampler$getLeafPrior()$leaf.prior
  results <- list()
  started <- proc.time()[[3L]]
  for (d in first:last) {
    set.seed(500000L + d)
    repeat {
      truth <- drawFromPrior(sampler, design$x)
      shape <- drawShape()
      y <- drawCounts(truth$eta, design$offset, shape, sbcTotalCap)
      if (!is.null(y)) break
    }
    drawFromPrior(sampler, design$x)
    sampler$setLeafPrior(defaultPrior)
    sampler$setResponse(y, updateScale = FALSE)
    sampler$run(sbcBurn - 1L, 1L)
    k <- shapes <- numeric()
    repeat {
      fit <- sampler$run(0L, sbcChunk)
      k <- c(k, as.numeric(fit$k))
      shapes <- c(shapes, as.numeric(fit$shape))
      tau <- iatGeyer(log(k))
      thin <- max(1L, ceiling(2 * tau))
      if (sbcNumRanks * thin <= length(k) || length(k) >= sbcCap) break
    }
    thin <- min(thin, length(k) %/% sbcNumRanks)
    kept <- length(k) - (seq_len(sbcNumRanks) - 1L) * thin
    ties <- sum(shapes[kept] == shape)
    results[[length(results) + 1L]] <- list(
      d = d,
      k0 = truth$k0,
      shape = shape,
      y = y,
      tau = tau,
      thin = thin,
      sweeps = length(k),
      rankK = sum(k[kept] < truth$k0),
      rankShape = sum(shapes[kept] < shape) + sample.int(ties + 1L, 1L) - 1L
    )
  }
  seconds <- proc.time()[[3L]] - started
  cat(sprintf("sbc %d to %d: %.0f s\n", first, last, seconds))
  saveRDS(
    list(results = results, seconds = seconds),
    outPath("checks", sprintf("sbc-%d-%d.rds", first, last))
  )
}

summarizeSbc <- function() {
  files <- list.files(file.path(runDir, "checks"), "^sbc-", full.names = TRUE)
  results <- do.call(c, lapply(lapply(files, readRDS), `[[`, "results"))
  chiSquared <- function(ranks) {
    counts <- tabulate(ranks %/% 10L + 1L, 10L)
    expected <- length(ranks) / 10
    statistic <- sum((counts - expected)^2 / expected)
    c(statistic = statistic, p = pchisq(statistic, 9L, lower.tail = FALSE))
  }
  rankK <- vapply(results, `[[`, 0L, "rankK")
  rankShape <- vapply(results, `[[`, 0, "rankShape")
  tau <- vapply(results, `[[`, 0, "tau")
  sweeps <- vapply(results, `[[`, 0, "sweeps")
  testK <- chiSquared(rankK)
  testShape <- chiSquared(rankShape)
  meanRankZ <- (mean(rankK / sbcNumRanks) - 0.5) / sqrt(1 / 12 / length(rankK))
  cat(sprintf(
    "SBC: %d datasets; k's ranks chi-squared %.1f on 9 df, p %.3f; mean normalized rank z %+.2f; shape's p %.3f\n",
    length(rankK),
    testK[["statistic"]],
    testK[["p"]],
    meanRankZ,
    testShape[["p"]]
  ))
  cat("  k's rank deciles:", tabulate(rankK %/% 10L + 1L, 10L), "\n")
  cat(sprintf(
    "  k's autocorrelation time median %.0f, 90th %.0f, max %.0f sweeps; %d at the cap\n",
    median(tau),
    quantile(tau, 0.9),
    max(tau),
    sum(sweeps >= sbcCap)
  ))
}

## ------------------------------------- stage 0: the builds on one dataset

# Two of the SBC datasets, fitted on the build this process sees, at the same
# center: 2 chains of 100,000 sweeps after 2,000 of burn-in, k recorded every
# 10th. summarize compares the builds' 5th, 50th and 95th percentiles of k by
# batch means.
compareSweeps <- 100000L
compareThin <- 10L

runCompare <- function(index, seed) {
  requireBuild(c("today", "move"))
  file <- list.files(file.path(runDir, "checks"), "^sbc-", full.names = TRUE)
  results <- do.call(c, lapply(lapply(file, readRDS), `[[`, "results"))
  dataset <- Filter(function(r) r$d == index, results)[[1L]]
  design <- sbcDesign()
  sampler <- sbcSampler(design, 2L)
  sampler$setControl(dbartsControl(
    n.chains = 2L,
    n.threads = 1L,
    n.samples = 1L,
    n.burn = 0L,
    n.thin = compareThin,
    n.trees = numTrees,
    updateState = FALSE,
    verbose = FALSE
  ))
  set.seed(seed)
  sampler$setResponse(dataset$y, updateScale = FALSE)
  started <- proc.time()[[3L]]
  sampler$run(2000L, 1L)
  fit <- sampler$run(0L, compareSweeps %/% compareThin)
  saveRDS(
    list(
      build = build,
      index = index,
      k = fit$k,
      seconds = proc.time()[[3L]] - started
    ),
    outPath("checks", sprintf("compare-%s-%d.rds", build, index))
  )
}

# Batch means: 20 batches a chain, the statistic's standard error from the
# spread of its batch values.
batchStatistic <- function(draws, statistic, numBatches = 20L) {
  draws <- as.matrix(draws)
  values <- unlist(
    lapply(seq_len(ncol(draws)), function(j) {
      batches <- split(
        draws[, j],
        rep(seq_len(numBatches), each = nrow(draws) / numBatches)
      )
      lapply(batches, statistic)
    }),
    recursive = FALSE
  )
  values <- do.call(rbind, values)
  list(
    estimate = statistic(as.vector(draws)),
    se = apply(values, 2L, sd) / sqrt(nrow(values))
  )
}

summarizeCompare <- function() {
  files <- list.files(
    file.path(runDir, "checks"),
    "^compare-",
    full.names = TRUE
  )
  runs <- lapply(files, readRDS)
  probabilities <- c(0.05, 0.5, 0.95)
  for (index in unique(vapply(runs, `[[`, 0, "index"))) {
    pair <- Filter(function(r) r$index == index, runs)
    if (length(pair) != 2L) {
      next
    }
    # Quantiles' Monte Carlo errors as the posterior package states them
    # (from the effective draws of the indicator at each quantile).
    stats <- lapply(pair, function(r) {
      draws <- posterior::as_draws_array(array(r$k, c(dim(as.matrix(r$k)), 1L)))
      list(
        estimate = quantile(r$k, probabilities, names = FALSE),
        se = as.numeric(posterior::mcse_quantile(draws, probabilities))
      )
    })
    z <- (stats[[1L]]$estimate - stats[[2L]]$estimate) /
      sqrt(stats[[1L]]$se^2 + stats[[2L]]$se^2)
    ess <- vapply(
      pair,
      function(r) {
        sum(apply(as.matrix(r$k), 2L, function(k) length(k) / iatGeyer(log(k))))
      },
      0
    )
    cat(sprintf(
      "dataset %d: k's 5th, 50th, 95th %s: %s; %s: %s; z %s; effective draws of k %.0f, %.0f\n",
      index,
      pair[[1L]]$build,
      paste(sprintf("%.3f", stats[[1L]]$estimate), collapse = " "),
      pair[[2L]]$build,
      paste(sprintf("%.3f", stats[[2L]]$estimate), collapse = " "),
      paste(sprintf("%+.2f", z), collapse = " "),
      ess[1L],
      ess[2L]
    ))
  }
}

## ---------------------------------------------------- stage 0: reweighting

# Log densities of k under the spread priors that are reweighted from C0's
# runs: chi(nu, scale) in k, normalized.
logChiDensity <- function(k, df, scale) {
  (df - 1) *
    log(k) -
    k^2 / (2 * scale^2) -
    (df / 2 - 1) * log(2) -
    lgamma(df / 2) -
    df * log(scale)
}

# The long run's rule: 4 chains, 1,000 sweeps of burn-in a chain, then
# chunks of 2,500 sweeps a chain (every 5th recorded) until k, the mean count
# over rows and 90 percent of rows' log mean reach 100 effective draws
# (batchEss, pooled over chains), capped at 100,000 kept sweeps a chain.
longChains <- 4L
longBurn <- 1000L
longChunk <- 2500L
longThin <- 5L
longTarget <- 100
longCap <- 100000L

# Effective draws a chain, the smaller of batch means (10 batches) and Geyer's
# initial monotone sequence, summed over chains. Batch means alone read too
# high while batches are short against the autocorrelation time, which early
# in a run they are; the smaller of the two keeps a run from stopping early.
batchEss <- function(draws, numBatches = 10L) {
  draws <- as.matrix(draws)
  size <- nrow(draws) %/% numBatches
  sum(apply(draws, 2L, function(x) {
    x <- x[seq_len(size * numBatches)]
    means <- colMeans(matrix(x, size))
    v <- var(x)
    if (!(v > 0)) {
      return(length(x))
    }
    batch <- length(x) * v / (size * var(means))
    min(length(x), batch, length(x) / iatGeyer(x))
  }))
}

longRun <- function(x, y, offset, kPrior, seed) {
  set.seed(seed)
  control <- dbartsControl(
    n.chains = longChains,
    n.threads = 1L,
    n.samples = 1L,
    n.burn = 0L,
    n.thin = longThin,
    n.trees = numTrees,
    updateState = FALSE,
    verbose = FALSE
  )
  # C0 is the family's own default, so its run states no leaf prior.
  arguments <- list(x, y, offset = offset, family = "nbinom", control = control)
  if (!is.null(kPrior)) {
    arguments$leaf.prior <- dbartsPriors$normal(k = kPrior)
  }
  sampler <- do.call(dbarts, arguments)
  kPriorUsed <- sampler$getLeafPrior()$leaf.prior@k
  started <- proc.time()[[3L]]
  sampler$run(longBurn, 1L)
  logMean <- NULL
  k <- shapes <- NULL
  repeat {
    fit <- sampler$run(0L, longChunk %/% longThin)
    eta <- fit$train - if (is.null(offset)) 0 else offset
    logMean <- if (is.null(logMean)) eta else abind2(logMean, eta)
    k <- rbind(k, fit$k)
    shapes <- rbind(shapes, fit$shape)
    swept <- nrow(k) * longThin
    essK <- batchEss(log(k))
    essMean <- batchEss(apply(exp(logMean), c(2L, 3L), mean))
    essRows <- apply(logMean, 1L, batchEss)
    settled <- essK >= longTarget &&
      essMean >= longTarget &&
      quantile(essRows, 0.1) >= longTarget
    if (settled || swept >= longCap) break
  }
  list(
    k = k,
    shape = shapes,
    logMean = logMean,
    essK = essK,
    essMean = essMean,
    essRows = quantile(essRows, 0.1, names = FALSE),
    settled = settled,
    sweeps = swept,
    seconds = proc.time()[[3L]] - started,
    kPrior = c(kPriorUsed@degreesOfFreedom, kPriorUsed@scale)
  )
}

# Bind along the draws (second) dimension of rows x draws x chains arrays.
abind2 <- function(a, b) {
  out <- array(NA_real_, c(dim(a)[1L], dim(a)[2L] + dim(b)[2L], dim(a)[3L]))
  out[, seq_len(dim(a)[2L]), ] <- a
  out[, dim(a)[2L] + seq_len(dim(b)[2L]), ] <- b
  out
}

# The cases: three datasets drawn from the default prior (75 trees, k0 from
# chi(1.5, 2), a shape from the grid prior) on 200 rows of the designed
# cells' covariates at a mean count near 2, read for chi(5, 0.9); and two
# replications of D5's null cell (500 rows, mean count 5, shape 5), read for
# chi(1.5, 3) and C7.
reweightCase <- function(case) {
  index <- as.integer(sub("^[a-z]+", "", case))
  set.seed(700000L + if (startsWith(case, "prior")) index else 100L + index)
  if (startsWith(case, "prior")) {
    n <- 200L
    x <- data.frame(correlatedNormals(n, 5L), b = rbinom(n, 1L, 0.5))
    build <- dbarts(
      x,
      rpois(n, 2),
      family = "nbinom",
      control = samplerControl()
    )
    repeat {
      truth <- drawFromPrior(build, x)
      y <- drawCounts(truth$eta, 0, drawShape())
      if (!is.null(y)) break
    }
    candidates <- list(`C2 chi(5,0.9)` = c(5, 0.9))
  } else {
    n <- 500L
    x <- data.frame(correlatedNormals(n, 5L), b = rbinom(n, 1L, 0.5))
    y <- rnbinom(n, size = 5, mu = 5)
    poisson <- glm(y ~ ., data = cbind(x, y = y), family = poisson())
    spread <- sd(predict(poisson, type = "link"))
    target <- 1.7 * max(spread, 0.25)
    candidates <- list(
      `C3 chi(1.5,3)` = c(1.5, 3),
      C7 = c(4, leafScale / (target * sqrt(qchisq(0.5, 4))))
    )
  }
  list(x = x, y = y, candidates = candidates)
}

weightedQuantile <- function(x, w, p) {
  o <- order(x)
  cw <- cumsum(w[o]) / sum(w)
  x[o][pmin(findInterval(p, cw) + 1L, length(x))]
}

# The reweighted reading of each candidate against its own run: k's 5th,
# 50th and 95th percentiles and each row's 90 percent interval endpoints,
# by batch means; and the weights' Kish share and Pareto k-hat.
runReweight <- function(case) {
  requireBuild("move")
  data <- reweightCase(case)
  base <- longRun(data$x, data$y, NULL, NULL, 1L)
  logK <- log(base$k)
  out <- list(
    case = case,
    base = base[c(
      "kPrior",
      "essK",
      "essMean",
      "essRows",
      "settled",
      "sweeps",
      "seconds"
    )]
  )
  probabilities <- c(0.05, 0.5, 0.95)
  readRun <- function(k, logMean, w) {
    # 5 batches a chain: at the run's target of 100 effective draws over 4
    # chains each batch then holds about 5, enough for a batch's quantile.
    numBatches <- 5L
    size <- nrow(k) %/% numBatches
    batchIndex <- function(j, b) (b - 1L) * size + seq_len(size)
    one <- function(rows, chains) {
      kk <- as.vector(k[rows, chains])
      ww <- as.vector(w[rows, chains])
      lm <- matrix(
        logMean[, rows, chains, drop = FALSE],
        nrow = dim(logMean)[1L]
      )
      c(
        weightedQuantile(kk, ww, probabilities),
        apply(lm, 1L, weightedQuantile, w = ww, p = 0.05),
        apply(lm, 1L, weightedQuantile, w = ww, p = 0.95)
      )
    }
    estimate <- one(seq_len(nrow(k)), seq_len(ncol(k)))
    batches <- do.call(
      rbind,
      lapply(seq_len(ncol(k)), function(j) {
        t(vapply(
          seq_len(numBatches),
          function(b) one(batchIndex(j, b), j),
          estimate
        ))
      })
    )
    list(estimate = estimate, se = apply(batches, 2L, sd) / sqrt(nrow(batches)))
  }
  for (label in names(data$candidates)) {
    prior <- data$candidates[[label]]
    logW <- logChiDensity(base$k, prior[1L], prior[2L]) -
      logChiDensity(base$k, 1.5, 2)
    w <- exp(logW - max(logW))
    kish <- sum(w)^2 / sum(w^2) / length(w)
    khat <- posterior::pareto_khat(as.vector(w), tail = "right", r_eff = 1)
    own <- longRun(
      data$x,
      data$y,
      NULL,
      dbartsPriors$chi(prior[1L], prior[2L]),
      2L
    )
    reweighted <- readRun(base$k, base$logMean, w)
    direct <- readRun(own$k, own$logMean, array(1, dim(own$k)))
    z <- (reweighted$estimate - direct$estimate) /
      sqrt(reweighted$se^2 + direct$se^2)
    out[[label]] <- list(
      prior = prior,
      kish = kish,
      khat = as.numeric(khat),
      own = own[c(
        "kPrior",
        "essK",
        "essMean",
        "essRows",
        "settled",
        "sweeps",
        "seconds"
      )],
      kReweighted = reweighted$estimate[1:3],
      kOwn = direct$estimate[1:3],
      zK = z[1:3],
      zRows = z[-(1:3)]
    )
  }
  saveRDS(out, outPath("checks", sprintf("reweight-%s.rds", case)))
}

summarizeReweight <- function() {
  files <- list.files(
    file.path(runDir, "checks"),
    "^reweight-",
    full.names = TRUE
  )
  for (run in lapply(files, readRDS)) {
    cat(sprintf(
      "%s: C0's run (chi(%g, %g)) %d sweeps a chain, effective draws k %.0f, mean %.0f, rows' 10th %.0f%s\n",
      run$case,
      run$base$kPrior[1L],
      run$base$kPrior[2L],
      run$base$sweeps,
      run$base$essK,
      run$base$essMean,
      run$base$essRows,
      if (run$base$settled) "" else " (unsettled)"
    ))
    for (label in setdiff(names(run), c("case", "base"))) {
      r <- run[[label]]
      cat(sprintf(
        "  %s: Kish share %.2f, k-hat %.2f; k's 5/50/95 reweighted %s, own %s, z %s; rows' endpoints within 2 se %.3f (max |z| %.1f); own run (chi(%.3g, %.3g)) %d sweeps\n",
        label,
        r$kish,
        r$khat,
        paste(sprintf("%.2f", r$kReweighted), collapse = " "),
        paste(sprintf("%.2f", r$kOwn), collapse = " "),
        paste(sprintf("%+.1f", r$zK), collapse = " "),
        mean(abs(r$zRows) <= 2),
        max(abs(r$zRows)),
        r$own$kPrior[1L],
        r$own$kPrior[2L],
        r$own$sweeps
      ))
    }
  }
}

## ------------------------------------------- stage 0: the prior draw's checks

# The unit spread (A = 3, so k = 3 gives each row's f a prior sd of 1); cut
# points as bart() sets them; the sampler's c as the formula states it, and a
# fit's c on Part A's datasets at the cell's mean count, as the convention
# means it to be (the realized total varies about its expectation, so this
# reads the spread of the realized mean count); and the top-rows shortcut for a count
# above the cap against full count draws, 1,000 forests a cell.
runPriorChecks <- function() {
  requireBuild(c("today", "move"))
  set.seed(31L)
  lines <- character()
  say <- function(...) {
    lines[length(lines) + 1L] <<- sprintf(...)
    cat(lines[length(lines)], "\n")
  }
  for (cell in c("norm5-500", "mix10off-2000", "norm5-10000", "roaches")) {
    config <- priorCell(cell)
    offset <- if (any(config$offset != 0)) config$offset
    sampler <- suppressWarnings(dbarts(
      config$x,
      config$build,
      offset = offset,
      family = "nbinom",
      control = samplerControl()
    ))
    say("%s: k.scale %g", cell, sampler$getLeafPrior()$k.scale)
    fit <- suppressWarnings(bart(
      config$x,
      config$build,
      offset = offset,
      family = "nbinom",
      n.chains = 1L,
      n.burn = 0L,
      n.samples = 1L,
      keepTrees = TRUE,
      verbose = FALSE
    ))
    say(
      "  cut points as bart() sets them: %s",
      identical(
        sampler$data@cutPoints,
        fit$fit$data@cutPoints
      )
    )
    shift <- sampler$getLeafPrior()$response.shift
    formula <- log(max(sum(config$build), 0.5) / sum(exp(config$offset)))
    say("  c as the formula states it: %s", isTRUE(all.equal(shift, formula)))
    sampler$setLeafPrior(dbartsPriors$normal(k = leafScale))
    logExposure <- log(sum(exp(config$offset)))
    numForests <- 1000L
    u <- matrix(NA_real_, length(config$offset), numForests)
    for (i in seq_len(numForests)) {
      sampler$sampleTreesFromPrior()
      sampler$sampleLeafParametersFromPrior()
      u[, i] <- as.numeric(sampler$predict(config$x)) - shift
    }
    rowVariance <- mean(apply(u, 1L, var))
    say(
      "  unit spread: rows' prior variance at k = 3 %.3f (1 expected)",
      rowVariance
    )
    prior <- spreadPriors$C0
    for (m in config$means) {
      gap <- shortcut <- full <- numeric(numForests)
      for (i in seq_len(numForests)) {
        s <- prior$draw(1L)
        shape <- drawShape()
        v <- config$offset + s * u[, i]
        logMean <- v + logExposure - logSumExp(v) + log(m)
        y <- suppressWarnings(rnbinom(
          length(logMean),
          size = shape,
          mu = exp(logMean)
        ))
        fitted <- log(max(sum(y), 0.5) / sum(exp(config$offset)))
        gap[i] <- fitted - log(m)
        full[i] <- anyNA(y) || max(y) > countCap
        top <- logMean[logMean > log(1e3)]
        shortcut[i] <- if (length(top) == 0L) {
          0
        } else {
          capProbability(top, shape)
        }
      }
      se <- sqrt(var(full - shortcut) / numForests)
      say(
        "  mean %g: a fit's c less log of the cell's mean count, median %+.3f (central 90%% %+.3f to %+.3f); count above the cap, full draws %.4f, shortcut %.4f, z %+.2f",
        m,
        median(gap),
        quantile(gap, 0.05),
        quantile(gap, 0.95),
        mean(full),
        mean(shortcut),
        if (se > 0) (mean(full) - mean(shortcut)) / se else 0
      )
    }
  }
  writeLines(lines, outPath("checks", "prior-checks.txt"))
}

## -------------------------------------------- stage 1: the Part A table

# One row per cell, candidate and mean count: A1, A2 (and A2 under C6's
# shape weights), A3 at the real designs, with standard errors; S's and
# S_id's quantiles; the share of datasets with more than 90 percent zeros
# expected; the real data's percentile and Traffic's at the low end.
summarizePrior <- function() {
  files <- list.files(file.path(runDir, "prior"), "\\.rds$", full.names = TRUE)
  files <- files[basename(files) != "partA.rds"]
  rows <- list()
  for (run in lapply(files, readRDS)) {
    above <- do.call(rbind, lapply(run$blocks, `[[`, "above"))
    scored <- do.call(abind3, lapply(run$blocks, `[[`, "scored"))
    m <- nrow(above)
    # A3 is also read at the middle of the data's range across the four fits,
    # an option offered with the Part A table.
    realMiddle <- if (!is.null(run$realTop)) {
      spreads <- readRDS(file.path(runDir, "data", "spreads.rds"))
      mean(range(spreads$S.id[spreads$data == run$cell]))
    }
    for (p in colnames(above)) {
      for (j in seq_along(run$means)) {
        sId <- scored[, p, "S.id", j]
        capped <- scored[, p, "capped", j]
        a3 <- if (!is.null(run$realTop)) as.numeric(sId > run$realTop)
        a3Middle <- if (!is.null(run$realTop)) as.numeric(sId > realMiddle)
        rows[[length(rows) + 1L]] <- data.frame(
          cell = run$cell,
          n = run$n,
          forests = m,
          prior = p,
          mean = run$means[j],
          A1 = mean(above[, p]),
          A1.se = sd(above[, p]) / sqrt(m),
          A2 = mean(capped),
          A2.se = sd(capped) / sqrt(m),
          A2.C6 = mean(scored[, p, "capped.C6", j]),
          A3.share = if (is.null(a3)) NA else mean(a3, na.rm = TRUE),
          A3.se = if (is.null(a3)) NA else sd(a3, na.rm = TRUE) / sqrt(m),
          real.top = if (is.null(run$realTop)) NA else run$realTop,
          A3.middle.share = if (is.null(a3)) {
            NA
          } else {
            mean(a3Middle, na.rm = TRUE)
          },
          real.middle = if (is.null(run$realTop)) NA else realMiddle,
          S.50 = quantile(scored[, p, "S", j], 0.5, names = FALSE),
          S.90 = quantile(scored[, p, "S", j], 0.9, names = FALSE),
          S.99 = quantile(scored[, p, "S", j], 0.99, names = FALSE),
          S.id.50 = quantile(sId, 0.5, na.rm = TRUE, names = FALSE),
          S.id.99 = quantile(sId, 0.99, na.rm = TRUE, names = FALSE),
          zeros90 = mean(scored[, p, "zeros90", j]),
          unidentified = mean(scored[, p, "unidentified", j]) / run$n,
          traffic.low = if (run$cell == "Traffic") {
            mean(sId < run$realTop, na.rm = TRUE)
          } else {
            NA
          },
          seconds = run$seconds
        )
      }
    }
  }
  table <- do.call(rbind, rows)
  saveRDS(table, outPath("prior", "partA.rds"))
  write.csv(table, outPath("prior", "partA.csv"), row.names = FALSE)
  simulated <- table[table$cell %not_in_set% names(realLoaders), ]
  real <- table[table$cell %in% setdiff(names(realLoaders), reportedOnly), ]
  priors <- unique(table$prior)
  worst <- function(d, column) {
    vapply(
      priors,
      function(p) {
        v <- d[d$prior == p, column]
        if (length(v) == 0L || all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
      },
      0
    )
  }
  summary <- data.frame(
    prior = priors,
    A1.worst = worst(table, "A1"),
    A2.worst = worst(table, "A2"),
    A2.C6.worst = worst(table, "A2.C6"),
    A3.lowest = vapply(
      priors,
      function(p) {
        v <- real$A3.share[real$prior == p]
        if (length(v) == 0L) NA_real_ else min(v)
      },
      0
    ),
    A3.failing = vapply(
      priors,
      function(p) {
        d <- real[real$prior == p, ]
        paste(d$cell[d$A3.share < partABound], collapse = " ")
      },
      ""
    )
  )
  print(summary, digits = 3L, row.names = FALSE)
  cat(sprintf(
    "forests a cell: %d to %d; %.1f core-hours\n",
    min(table$forests),
    max(table$forests),
    sum(unique(table[c("cell", "seconds")])$seconds) / 3600
  ))
}

## ------------------------------------------------- truths fitted by BART

# The posterior mean log mean of a dbarts fit at k fixed at 1. It is a fixed
# function once stored, so its Monte Carlo error does not enter the study:
# 2 chains of 1,000 burn-in and 1,000 kept are enough to give it a BART
# fit's interactions and spread.
runTruthBart <- function(name) {
  requireBuild("today")
  data <- readData(name)
  started <- proc.time()[[3L]]
  fit <- bart(
    data$x,
    data$y,
    offset = data$offset,
    family = "nbinom",
    k = 1,
    n.chains = 2L,
    n.threads = 1L,
    n.burn = 1000L,
    n.samples = 1000L,
    verbose = FALSE,
    seed = 99L
  )
  eta <- colMeans(extract(fit, "bart")) -
    if (is.null(data$offset)) 0 else data$offset
  identified <- exp(eta + if (is.null(data$offset)) 0 else data$offset) >=
    identifiedMean
  saveRDS(
    list(
      eta = eta,
      shape = median(fit$shape),
      S = sd(eta),
      S.id = sd(eta[identified])
    ),
    outPath("data", paste0(name, "-truth-bart.rds"))
  )
  cat(sprintf(
    "%s: S %.2f, S_id %.2f, shape %g, %.0f s\n",
    name,
    sd(eta),
    sd(eta[identified]),
    median(fit$shape),
    proc.time()[[3L]] - started
  ))
}

## ------------------------------------------------------------ smoke runs

# NBS_QUICK=1 shortens every run and count for a smoke run of the pipeline;
# its output says nothing and belongs in a run directory of its own.
if (nzchar(Sys.getenv("NBS_QUICK"))) {
  longChunk <- 250L
  longTarget <- 10
  longCap <- 2000L
  sbcBurn <- 200L
  sbcChunk <- 1000L
  sbcCap <- 5000L
  compareSweeps <- 2000L
  priorBlockSize <- 100L
  priorMinimum <- 200L
  priorMaximum <- 400L
}

## ---------------------------------------------------------------- dispatch

args <- commandArgs(trailingOnly = TRUE)
action <- if (length(args) > 0L) args[1L] else ""
switch(
  action,
  export = exportData(),
  drift = runDrift(as.integer(args[2L]), as.integer(args[3L])),
  prior = runPrior(args[2L]),
  sbc = runSbc(as.integer(args[2L]), as.integer(args[3L])),
  compare = runCompare(as.integer(args[2L]), as.integer(args[3L])),
  reweight = runReweight(args[2L]),
  `prior-checks` = runPriorChecks(),
  truth = switch(
    args[2L],
    bart = runTruthBart(args[3L]),
    stop("truth takes bart; the additive truth is written by export")
  ),
  summarize = switch(
    args[2L],
    drift = summarizeDrift(),
    sbc = summarizeSbc(),
    compare = summarizeCompare(),
    reweight = summarizeReweight(),
    prior = summarizePrior(),
    stop("unknown stage to summarize")
  ),
  stop("unknown action; see the usage at the top of this file")
)
