#!/usr/bin/env Rscript

# The negative binomial default-prior study: is there a better default prior
# on the forest's spread for family = "nbinom" fits? The design, the criterion
# and the decision rule are fixed in docs/plans/nbinom-default-study.md; this
# file implements its stages 0 to 2 (the checks, the prior predictive and the
# long reading's check) and Part B's main run and scoring, and decides
# nothing the plan leaves open except where its comments say so.
#
# BUILDS. One R process sees one build: whatever dbarts R_LIBS finds first,
# named by NBS_BUILD ("today", the shipped sampler; "move", the scratch build
# carrying the collapsed scale move's negative binomial half; "mutant", the
# move build with one power of r too many in the move's Jacobian, for the
# drift test's power; "kprior", the move build whose k step and move take the
# prior on k that NBS_KPRIOR names, read at each use: unset for the shipped
# chi(1.5, 2), "floor:f" for C4 and "halfnormal:m" for C5, with NBS_MOVE=off
# turning the move off, which leaves today's sampler with only the k step).
# Actions that the plan assigns to a build refuse to run under another.
#
# RUN DIRECTORY. NBS_RUN holds every output: data/ (the exported real
# datasets and their truths), checks/, prior/ and the summaries. Nothing is
# read from the working directory.
#
# Usage (each action one process):
#   Rscript nbinom-default-study.R export          (data and additive truths)
#   Rscript nbinom-default-study.R truth bart <dataset>
#   Rscript nbinom-default-study.R truth cf <dataset>
#   Rscript nbinom-default-study.R drift <seed> <reps>
#   Rscript nbinom-default-study.R identity <seed>
#   Rscript nbinom-default-study.R sbc <first> <last>
#   Rscript nbinom-default-study.R compare <SBC dataset> <seed>
#   Rscript nbinom-default-study.R reweight <prior1..3 | null1..2>
#   Rscript nbinom-default-study.R prior-checks
#   Rscript nbinom-default-study.R prior <cell>
#   Rscript nbinom-default-study.R stage2 <cell> <replication>
#   Rscript nbinom-default-study.R correction <cell> <replication>
#   Rscript nbinom-default-study.R plan
#   Rscript nbinom-default-study.R long <cell> <replication>
#   Rscript nbinom-default-study.R default <cell> <replication> [<candidate>]
#   Rscript nbinom-default-study.R guardplan <candidate>
#   Rscript nbinom-default-study.R calibrate <cell>
#   Rscript nbinom-default-study.R uncoded
#   Rscript nbinom-default-study.R summarize <stage> [<look>]
# where <stage> is drift, identity, sbc, compare, reweight, prior, stage2,
# correction, calibration or main, and <candidate> a fitted candidate's file code (its
# name without spaces, parentheses or commas, as C2chi41.03).
#
# BLIND. Nothing here scores a fit against a truth until the pre-registration
# commit sets scoringOn: before it the Part B scoring runs only on a quick
# run's toy fits and prints only the shape of what it computed. No run
# length, weight or acceptance below reads a truth. The stage 0 checks compare the sampler
# with itself and with the prior, which is not scoring. One exception, the
# amendment after stage 1's: stage2Blind computes a candidate's coverage of
# the truth less C0's on the spread cells, under coded labels, and the
# summary reads only its spread within cells, as blinded sample-size
# re-estimation does.

suppressPackageStartupMessages(library(dbarts))

## ---------------------------------------------------------------- settings

`%not_in_set%` <- function(x, table) !(x %in% table)

runDir <- Sys.getenv("NBS_RUN")
if (!nzchar(runDir)) {
  stop("NBS_RUN must name the run directory")
}
build <- Sys.getenv("NBS_BUILD")
if (build %not_in_set% c("today", "move", "mutant", "kprior")) {
  stop("NBS_BUILD must be today, move, mutant or kprior")
}
if (build != "kprior" && nzchar(Sys.getenv("NBS_KPRIOR"))) {
  stop("NBS_KPRIOR is read only by the kprior build")
}
# A kprior process's outputs carry its prior and whether the move is off.
buildLabel <- if (build == "kprior") {
  paste0(
    "kprior",
    gsub("[^A-Za-z0-9.]", "", Sys.getenv("NBS_KPRIOR")),
    if (Sys.getenv("NBS_MOVE") == "off") "-moveoff" else ""
  )
} else {
  build
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

# k0 from the prior that an NBS_KPRIOR setting names, as the kprior build
# samples it: chi(1.5, 2), truncated below at the floor, or 3 over a
# half-normal of the given median.
priorKDraw <- function(setting) {
  if (!nzchar(setting)) {
    return(function() 2 * sqrt(rchisq(1L, 1.5)))
  }
  parts <- strsplit(setting, ":", fixed = TRUE)[[1L]]
  value <- as.numeric(parts[2L])
  switch(
    parts[1L],
    floor = function() {
      upper <- pchisq((value / 2)^2, 1.5, lower.tail = FALSE)
      2 * sqrt(qchisq(runif(1L) * upper, 1.5, lower.tail = FALSE))
    },
    halfnormal = function() {
      leafScale / abs(rnorm(1L, 0, value / qnorm(0.75)))
    },
    stop("NBS_KPRIOR names no prior")
  )
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
  requireBuild(c("move", "mutant", "today", "kprior"))
  drawK <- priorKDraw(Sys.getenv("NBS_KPRIOR"))
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
      truth <- drawFromPrior(sampler, x, drawK)
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
    buildLabel,
    seed,
    numReps,
    seconds
  ))
  saveRDS(
    list(
      build = buildLabel,
      lags = driftLags,
      results = results,
      seconds = seconds
    ),
    outPath("checks", sprintf("drift-%s-%d.rds", buildLabel, seed))
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

## ------------------------------------------- the kprior build's identity

# The kprior build must draw as the builds it extends: with NBS_KPRIOR unset,
# as the move build; with NBS_MOVE=off as well, as today's sampler; and under
# a floor that cannot bind (floor:1e-6), which takes the shipped gamma draw
# and accepts it, as either. One short fit at a fixed seed, its k and fits
# saved for summarize identity to compare draw for draw.
runIdentity <- function(seed) {
  set.seed(seed)
  n <- 200L
  x <- correlatedNormals(n, 3L)
  offset <- runif(n, 0, log(100))
  eta <- 0.8 * sin(x[, 1L]) + 0.5 * x[, 2L] * x[, 3L] - 2
  y <- rnbinom(n, size = 2, mu = exp(eta + offset))
  control <- dbartsControl(
    n.chains = 2L,
    n.threads = 1L,
    n.samples = 300L,
    n.burn = 300L,
    n.trees = numTrees,
    updateState = FALSE,
    verbose = FALSE
  )
  sampler <- dbarts(
    x,
    y,
    offset = offset,
    family = "nbinom",
    control = control
  )
  fit <- sampler$run()
  saveRDS(
    list(build = buildLabel, k = fit$k, train = fit$train, shape = fit$shape),
    outPath("checks", sprintf("identity-%s-%d.rds", buildLabel, seed))
  )
}

summarizeIdentity <- function() {
  files <- list.files(
    file.path(runDir, "checks"),
    "^identity-.*[.]rds$",
    full.names = TRUE
  )
  runs <- lapply(files, readRDS)
  names(runs) <- sub("^identity-(.*)[.]rds$", "\\1", basename(files))
  pairs <- list(
    c("kprior", "move"),
    c("kpriorfloor1e6", "move"),
    c("kprior-moveoff", "today"),
    c("kpriorfloor1e6-moveoff", "today")
  )
  for (seed in unique(sub(".*-", "", names(runs)))) {
    for (pair in pairs) {
      a <- runs[[paste0(pair[1L], "-", seed)]]
      b <- runs[[paste0(pair[2L], "-", seed)]]
      if (is.null(a) || is.null(b)) {
        next
      }
      same <- identical(a$k, b$k) &&
        identical(a$train, b$train) &&
        identical(a$shape, b$shape)
      cat(sprintf(
        "seed %s: %s against %s: %s (largest k difference %.3g)\n",
        seed,
        pair[1L],
        pair[2L],
        if (same) "identical" else "DIFFERENT",
        max(abs(a$k - b$k))
      ))
    }
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

## ------------------------------------------- counterfactual truths for R

# R on a real dataset reads the truth at both levels of one binary column:
# the additive truth refitted by the export's rule (deterministic, and
# checked to reproduce the stored truth) and the BART-fitted truth refitted at
# its seed with the flipped rows as test rows (prediction draws nothing, so
# its training fit is the stored one, also checked). Written to
# data/<dataset>-cf.rds: the column, its two levels, and each truth at each.
binaryColumns <- c(
  roaches = "treatment",
  epil = "trt",
  Traffic = "limit",
  Arabidopsis = "amd",
  solder = "Solder",
  mite = "Topo"
)
runTruthCounterfactual <- function(name) {
  requireBuild("today")
  suppressPackageStartupMessages(library(mgcv))
  column <- binaryColumns[[name]]
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
  basisSize <- function(v) min(10L, length(unique(d[[v]])) - 1L)
  terms <- c(
    sprintf("s(%s, k = %d)", numeric, vapply(numeric, basisSize, 1L)),
    factors,
    if (is.null(s$off)) character() else sprintf("offset(%s)", s$off)
  )
  formula <- as.formula(paste(s$y, "~", paste(terms, collapse = " + ")))
  fit <- suppressWarnings(gam(
    formula,
    data = d,
    family = nb(),
    method = "REML"
  ))
  stored <- readData(name)
  levels <- levels(d[[column]])
  at <- function(level) {
    flipped <- d
    flipped[[column]] <- factor(level, levels = levels)
    as.numeric(predict(fit, newdata = flipped, type = "link")) - offset
  }
  additive <- cbind(at(levels[1L]), at(levels[2L]))
  additiveCheck <- max(abs(
    as.numeric(predict(fit, type = "link")) - offset - stored$truthAdditive$eta
  ))
  x0 <- x1 <- stored$x
  x0[[column]] <- factor(levels[1L], levels = levels(stored$x[[column]]))
  x1[[column]] <- factor(levels[2L], levels = levels(stored$x[[column]]))
  bartFit <- bart(
    stored$x,
    stored$y,
    test = rbind(x0, x1),
    offset = stored$offset,
    offset.test = if (is.null(stored$offset)) NULL else rep(stored$offset, 2L),
    family = "nbinom",
    k = 1,
    n.chains = 2L,
    n.threads = 1L,
    n.burn = 1000L,
    n.samples = 1000L,
    verbose = FALSE,
    seed = 99L
  )
  o <- if (is.null(stored$offset)) 0 else stored$offset
  bartTrain <- colMeans(extract(bartFit, "bart")) - o
  bartTest <- colMeans(extract(bartFit, "bart", sample = "test")) -
    rep(o, length.out = 2L * nrow(stored$x))
  storedBart <- readRDS(file.path(
    runDir,
    "data",
    paste0(name, "-truth-bart.rds")
  ))
  bartCheck <- max(abs(bartTrain - storedBart$eta))
  n <- nrow(stored$x)
  out <- list(
    column = column,
    levels = levels,
    additive = additive,
    bart = cbind(bartTest[seq_len(n)], bartTest[n + seq_len(n)]),
    reproduced = c(additive = additiveCheck, bart = bartCheck)
  )
  saveRDS(out, outPath("data", paste0(name, "-cf.rds")))
  cat(sprintf(
    "%s, %s: stored truths reproduced to %.1e (additive) and %.1e (BART)\n",
    name,
    column,
    additiveCheck,
    bartCheck
  ))
}

## ------------------------------------------- stage 2: the long reading's check

# The five candidates that pass Part A under dec-B502, read from C0's run by
# the weights pi_c(k) / pi_0(k): each a log density in k up to a constant,
# whether its weight grows without bound as k grows (the half-normal, whose
# density in k falls as 1 / k^2 against C0's exp(-k^2 / 8)), and its own run's
# chi prior where the move build can run it. C4's and C5's own runs need a
# scratch build: stage 2 counts where they would, and runs none.
logHalfNormalSpread <- function(k, median) {
  sigma <- median / qnorm(0.75)
  -2 * log(k) - (leafScale / k)^2 / (2 * sigma^2)
}
logFloored <- function(k, kFloor) {
  ifelse(k >= kFloor, logChiDensity(k, 1.5, 2), -Inf)
}
stage2Candidates <- list(
  `C2 chi(4,1.03)` = list(
    logDensity = function(k) logChiDensity(k, 4, 1.03),
    unbounded = FALSE,
    own = c(4, 1.03)
  ),
  `C2 chi(5,0.9)` = list(
    logDensity = function(k) logChiDensity(k, 5, 0.9),
    unbounded = FALSE,
    own = c(5, 0.9)
  ),
  `C4 k>=0.6` = list(
    logDensity = function(k) logFloored(k, 0.6),
    unbounded = FALSE,
    own = NULL
  ),
  `C4 k>=0.75` = list(
    logDensity = function(k) logFloored(k, 0.75),
    unbounded = FALSE,
    own = NULL
  ),
  `C5 hn 1.3` = list(
    logDensity = function(k) logHalfNormalSpread(k, 1.3),
    unbounded = TRUE,
    own = NULL
  )
)
# The amended rules: a candidate's weights hold on a replication when they
# keep 30 percent of the run's effective draws, their Pareto k-hat is at most
# 0.7, and, for a weight without bound, the likelihood is not flat at large k
# (the narrow band, the null cell, or more than 5 percent of C0's draws above
# k = 4). A run is extended until C0 and every candidate whose weights hold
# reach the target after reweighting; a candidate below half the target after
# reweighting gets its own run.
kishFloor <- 0.3
khatCeiling <- 0.7
ownBelow <- 50

candidateWeights <- function(candidate, k, flat) {
  logW <- candidate$logDensity(k) - logChiDensity(k, 1.5, 2)
  if (all(logW == -Inf)) {
    return(list(w = array(0, dim(k)), kish = 0, khat = Inf, holds = FALSE))
  }
  w <- exp(logW - max(logW))
  kish <- sum(w)^2 / sum(w^2) / length(w)
  # Weights taking one positive value (a truncation) have no tail to fit.
  khat <- if (length(unique(w[w > 0])) > 1L) {
    as.numeric(posterior::pareto_khat(as.vector(w), tail = "right", r_eff = 1))
  } else {
    -Inf
  }
  holds <- kish >= kishFloor &&
    khat <= khatCeiling &&
    !(candidate$unbounded && flat)
  list(w = array(w, dim(k)), kish = kish, khat = khat, holds = holds)
}

# A smooth truth with interactions, for the designed cells, scaled to sd S
# over the cell's rows and centered at its mean count.
designedCell <- function(n, spread, meanCount, shape, seed) {
  set.seed(seed)
  x <- data.frame(correlatedNormals(n, 5L), b = rbinom(n, 1L, 0.5))
  g <- with(x, sin(x1) + 0.5 * x2 * x3 + 0.5 * x4^2 + 0.5 * b * x5)
  f <- if (spread > 0) spread * (g - mean(g)) / sd(g) else rep(0, n)
  list(
    x = x,
    offset = NULL,
    truth = f + log(meanCount) - log(mean(exp(f))),
    shape = shape
  )
}

spreadIdentified <- function(eta, offset) {
  identified <- exp(eta + if (is.null(offset)) 0 else offset) >= identifiedMean
  sd(eta[identified])
}

# Widens a truth to S_id of target by scaling its departure from its mean,
# at the smallest scale that reaches the target: S_id is not monotone in the
# scale, since rows falling below a mean of 0.01 leave it. With meanCount,
# the center keeps the mean count; otherwise the center is lowered until the
# expected total is at most 20,000 and the largest expected count at most
# 10^5, the plan's rule.
widenTruth <- function(eta, offset, target, meanCount = NULL) {
  o <- if (is.null(offset)) 0 else offset
  m <- mean(eta)
  at <- function(a, shift) {
    e <- m + shift + a * (eta - m)
    if (!is.null(meanCount)) {
      e <- e - log(mean(exp(e + o))) + log(meanCount)
    }
    e
  }
  scales <- exp(seq(log(0.01), log(50), length.out = 400L))
  for (shift in seq(0, -20, by = -0.1)) {
    reached <- vapply(
      scales,
      function(a) spreadIdentified(at(a, shift), offset) >= target,
      NA
    )
    first <- which(reached %in% TRUE)[1L]
    if (is.na(first) || first == 1L) {
      next
    }
    a <- uniroot(
      function(a) spreadIdentified(at(a, shift), offset) - target,
      scales[first - 1:0]
    )$root
    e <- at(a, shift)
    mu <- exp(e + o)
    if (!is.null(meanCount) || (sum(mu) <= 20000 && max(mu) <= 1e5)) {
      # the affine map taken, for a counterfactual truth: m + shift + a (eta - m)
      attr(e, "map") <- c(
        center = m,
        shift = e[1L] - m - a * (eta[1L] - m),
        scale = a
      )
      return(e)
    }
  }
  stop("no center meets the caps")
}

# Stage 2's cells: the real datasets as fitted but BCI, Insurance and Traffic
# (additive truths); three widened; D4, D5 and D8; and small-sample cells,
# whose rows are drawn afresh each replication from the dataset's truth
# widened by the widened cells' rule, at shape 3. The plan's mean count of 5
# for them is not used: widened to S_id 3 at that mean, only 8 to 30 percent
# of a dataset's rows keep a mean of 0.01, and a subsample of 30 is mostly
# zeros; the pre-registration settles it. The spread cells are run to
# five replications for the paired spread, blind; today's sampler runs on the
# narrow and middle cells but solder, whose trees settle within the cap on
# neither build.
stage2CellNames <- c(
  "quakes",
  "D5",
  "D8",
  "roaches",
  "solder",
  "Arabidopsis",
  "grouseticks",
  "epil",
  "aids",
  "mite",
  "crimtab",
  "roaches-w3",
  "crimtab-w3",
  "quakes-w2",
  "D4-3",
  "roaches-n30-3",
  "quakes-n100-2",
  "aids-n30-2"
)
stage2SpreadCells <- c(
  "quakes",
  "D5",
  "D8",
  "D4-3",
  "roaches-n30-3",
  "quakes-n100-2",
  "aids-n30-2"
)
stage2TodayCells <- c(
  "quakes",
  "D5",
  "D8",
  "roaches",
  "Arabidopsis",
  "grouseticks",
  "epil"
)
stage2SpreadReps <- 5L

stage2Cell <- function(cell, rep) {
  index <- match(cell, stage2CellNames)
  if (is.na(index)) {
    stop("unknown stage 2 cell ", cell)
  }
  parts <- strsplit(cell, "-", fixed = TRUE)[[1L]]
  subsample <- NULL
  if (parts[1L] %in% c("D4", "D5", "D8")) {
    base <- switch(
      parts[1L],
      D4 = designedCell(60L, 3, 5, 3, 810000L),
      D5 = designedCell(500L, 0, 5, 5, 850000L),
      D8 = designedCell(1000L, 0.3, 5, 3, 880000L)
    )
    null <- parts[1L] == "D5"
  } else {
    d <- readData(parts[1L])
    base <- list(
      x = d$x,
      offset = d$offset,
      truth = as.numeric(d$truthAdditive$eta),
      shape = max(d$truthAdditive$theta, 1)
    )
    null <- FALSE
    if (length(parts) == 2L) {
      base$truth <- widenTruth(
        base$truth,
        base$offset,
        as.numeric(sub("w", "", parts[2L]))
      )
    } else if (length(parts) == 3L) {
      base$truth <- widenTruth(
        base$truth,
        base$offset,
        as.numeric(parts[3L])
      )
      base$shape <- 3
      subsample <- as.integer(sub("n", "", parts[2L]))
    }
  }
  # Rows: fixed for a cell, its held-out fifth fixed too; drawn afresh each
  # replication for a small-sample cell.
  set.seed(900000L + 1000L * index + if (is.null(subsample)) 0L else rep)
  rows <- if (is.null(subsample)) {
    seq_len(nrow(base$x))
  } else {
    sample.int(nrow(base$x), subsample)
  }
  heldOut <- sample.int(length(rows), round(length(rows) / 5))
  # Counts: drawn per replication at the truth.
  set.seed(950000L + 1000L * index + rep)
  offset <- if (is.null(base$offset)) NULL else base$offset[rows]
  mu <- exp(base$truth[rows] + if (is.null(offset)) 0 else offset)
  repeat {
    y <- rnbinom(length(rows), size = base$shape, mu = mu)
    if (max(y) <= countCap) break
  }
  x <- base$x[rows, , drop = FALSE]
  # The band and S_id are the cell's, over all its rows; a small-sample
  # replication's rows are a draw from them.
  sId <- if (null) 0 else spreadIdentified(base$truth, base$offset)
  identifiedShare <- mean(
    exp(base$truth + if (is.null(base$offset)) 0 else base$offset) >=
      identifiedMean
  )
  band <- if (null || sId < 0.5) {
    "narrow"
  } else if (length(rows) <= 130L && sId >= 2 - 1e-6) {
    "small"
  } else if (sId < 1.5) {
    "middle"
  } else {
    "wide"
  }
  train <- setdiff(seq_along(rows), heldOut)
  list(
    cell = cell,
    rep = rep,
    band = band,
    sId = sId,
    identifiedShare = identifiedShare,
    shape = base$shape,
    xTrain = x[train, , drop = FALSE],
    yTrain = y[train],
    offsetTrain = if (is.null(offset)) NULL else offset[train],
    xTest = x[heldOut, , drop = FALSE],
    offsetTest = if (is.null(offset)) NULL else offset[heldOut],
    # Read only by the blind paired spread (stage2Blind).
    truth = base$truth[rows][c(train, heldOut)],
    offsetAll = if (is.null(offset)) NULL else offset[c(train, heldOut)]
  )
}

bindDraws <- function(a, b) if (is.null(a)) b else abind2(a, b)

# One long run on the build this process sees: C0's prior when kPrior is
# NULL, else k's prior as given; extended in chunks until C0's effective
# draws of k, the mean and the 10th percentile of rows' log mean (on up to
# 200 rows) reach the target scaled by 1 / Kish of every candidate whose
# weights hold, or the cap. The trajectory records each chunk's numbers, from
# which the length any other rule would have taken is read.
stage2Run <- function(data, kPrior, seed, candidates, flat, extend = TRUE) {
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
  arguments <- list(
    data$xTrain,
    data$yTrain,
    test = data$xTest,
    family = "nbinom",
    control = control
  )
  if (!is.null(data$offsetTrain)) {
    arguments$offset <- data$offsetTrain
    arguments$offset.test <- data$offsetTest
  }
  if (!is.null(kPrior)) {
    arguments$leaf.prior <- dbartsPriors$normal(k = kPrior)
  }
  sampler <- do.call(dbarts, arguments)
  offsetTrain <- if (is.null(data$offsetTrain)) 0 else data$offsetTrain
  offsetTest <- if (is.null(data$offsetTest)) 0 else data$offsetTest
  checkRows <- sort(sample.int(
    length(data$yTrain),
    min(200L, length(data$yTrain))
  ))
  started <- proc.time()[[3L]]
  sampler$run(longBurn, 1L)
  train <- test <- k <- shapes <- trajectory <- NULL
  repeat {
    fit <- sampler$run(0L, longChunk %/% longThin)
    train <- bindDraws(train, fit$train - offsetTrain)
    test <- bindDraws(test, fit$test - offsetTest)
    k <- rbind(k, fit$k)
    shapes <- rbind(shapes, fit$shape)
    swept <- nrow(k) * longThin
    essK <- batchEss(log(k))
    essMean <- batchEss(apply(exp(train), c(2L, 3L), mean))
    essRows <- quantile(
      apply(train[checkRows, , , drop = FALSE], 1L, batchEss),
      0.1,
      names = FALSE
    )
    flatHere <- flat || mean(k > 4) > 0.05
    weights <- lapply(candidates, candidateWeights, k = k, flat = flatHere)
    kish <- vapply(weights, `[[`, 0, "kish")
    khat <- vapply(weights, `[[`, 0, "khat")
    holds <- vapply(weights, `[[`, NA, "holds")
    prefixed <- function(v, prefix) {
      if (length(v)) as.list(setNames(v, paste0(prefix, names(v)))) else list()
    }
    row <- c(
      list(
        sweeps = swept,
        seconds = proc.time()[[3L]] - started,
        essK = essK,
        essMean = essMean,
        essRows = essRows
      ),
      prefixed(kish, "kish:"),
      prefixed(khat, "khat:"),
      prefixed(holds, "holds:")
    )
    trajectory <- rbind(trajectory, as.data.frame(row, check.names = FALSE))
    need <- if (extend) max(c(1, 1 / kish[holds])) else 1
    if (min(essK, essMean, essRows) >= longTarget * need || swept >= longCap) {
      break
    }
  }
  base <- min(essK, essMean, essRows)
  list(
    k = k,
    shape = shapes,
    train = train,
    test = test,
    weights = weights,
    flat = flatHere,
    trajectory = trajectory,
    essK = essK,
    essMean = essMean,
    essRows = essRows,
    settled = base >= longTarget * need,
    settledC0 = base >= longTarget,
    sweeps = swept,
    seconds = proc.time()[[3L]] - started,
    checkRows = checkRows
  )
}

# k's 5th, 50th and 95th percentiles and each row's 90 percent interval, by
# weights w, with standard errors by batch means (5 batches a chain).
readQuantiles <- function(k, logMean, w, numBatches = 5L) {
  size <- nrow(k) %/% numBatches
  one <- function(rows, chains) {
    kk <- as.vector(k[rows, chains])
    ww <- as.vector(w[rows, chains])
    lm <- matrix(logMean[, rows, chains, drop = FALSE], nrow = dim(logMean)[1L])
    c(
      weightedQuantile(kk, ww, c(0.05, 0.5, 0.95)),
      as.vector(t(apply(lm, 1L, weightedQuantile, w = ww, p = c(0.05, 0.95))))
    )
  }
  estimate <- one(seq_len(nrow(k)), seq_len(ncol(k)))
  batches <- do.call(
    rbind,
    lapply(seq_len(ncol(k)), function(j) {
      t(vapply(
        seq_len(numBatches),
        function(b) one((b - 1L) * size + seq_len(size), j),
        estimate
      ))
    })
  )
  list(estimate = estimate, se = apply(batches, 2L, sd) / sqrt(nrow(batches)))
}

# The expected-coverage correction's check: on the check rows, an interval
# from chains 1 and 2 (all their draws, and their first quarter) covers the
# other two chains' draws, weighted, against 0.9 N / (N + 1), N the first
# chains' effective draws of that row after reweighting.
correctionCheck <- function(run, w) {
  first <- 1:2
  second <- 3:4
  out <- NULL
  for (share in c(1, 0.25)) {
    draws <- seq_len(floor(nrow(run$k) * share))
    w1 <- w[draws, first, drop = FALSE]
    if (!(sum(w1) > 0) || !(sum(w[, second]) > 0)) {
      next
    }
    kish <- sum(w1)^2 / sum(w1^2) / length(w1)
    w2 <- as.vector(w[, second])
    for (i in run$checkRows) {
      x1 <- run$train[i, draws, first, drop = FALSE]
      ends <- weightedQuantile(as.vector(x1), as.vector(w1), c(0.05, 0.95))
      x2 <- as.vector(run$train[i, , second])
      out <- rbind(
        out,
        c(
          share = share,
          n = batchEss(matrix(x1, ncol = 2L)) * kish,
          cover = sum(w2 * (x2 >= ends[1L] & x2 <= ends[2L])) / sum(w2)
        )
      )
    }
  }
  as.data.frame(out)
}

# Blind: each candidate's coverage of the truth's log mean by its 90 percent
# interval, over the rows whose true mean is at least 0.01, less C0's, under
# coded labels. Only the summary's spread within cells is read from it.
stage2Blind <- function(data, run, own) {
  logMean <- abind1(run$train, run$test)
  identified <- exp(
    data$truth + if (is.null(data$offsetAll)) 0 else data$offsetAll
  ) >=
    identifiedMean
  coverage <- function(draws, w) {
    ends <- apply(
      matrix(draws, nrow = dim(draws)[1L]),
      1L,
      weightedQuantile,
      w = as.vector(w),
      p = c(0.05, 0.95)
    )
    inside <- data$truth >= ends[1L, ] & data$truth <= ends[2L, ]
    mean(inside[identified])
  }
  base <- coverage(logMean, array(1, dim(run$k)))
  codes <- setNames(paste0("c", blindCodes), names(stage2Candidates))
  difference <- vapply(
    names(stage2Candidates),
    function(name) {
      if (run$weights[[name]]$holds) {
        coverage(logMean, run$weights[[name]]$w) - base
      } else if (!is.null(own[[name]])) {
        o <- own[[name]]
        coverage(abind1(o$train, o$test), array(1, dim(o$k))) - base
      } else {
        NA_real_
      }
    },
    0
  )
  setNames(difference, codes[names(difference)])
}
blindCodes <- c(3L, 5L, 1L, 4L, 2L)

# Bind along the rows (first) dimension.
abind1 <- function(a, b) {
  if (is.null(b) || dim(b)[1L] == 0L) {
    return(a)
  }
  out <- array(NA_real_, c(dim(a)[1L] + dim(b)[1L], dim(a)[2L], dim(a)[3L]))
  out[seq_len(dim(a)[1L]), , ] <- a
  out[dim(a)[1L] + seq_len(dim(b)[1L]), , ] <- b
  out
}

runStage2 <- function(cell, rep) {
  requireBuild(c("move", "today"))
  data <- stage2Cell(cell, rep)
  flat <- data$band == "narrow"
  candidates <- if (build == "move") stage2Candidates else list()
  run <- stage2Run(data, NULL, 960000L + rep, candidates, flat)
  logMean <- abind1(run$train, run$test)
  reading <- readQuantiles(run$k, logMean, array(1, dim(run$k)))
  essAll <- apply(logMean, 1L, batchEss)
  summary <- list(
    cell = cell,
    rep = rep,
    build = build,
    band = data$band,
    sId = data$sId,
    identifiedShare = data$identifiedShare,
    shape = data$shape,
    n = length(data$yTrain) + nrow(data$xTest),
    total = sum(data$yTrain),
    trajectory = run$trajectory,
    essK = run$essK,
    essMean = run$essMean,
    essRows = quantile(essAll, 0.1, names = FALSE),
    settled = run$settled,
    settledC0 = run$settledC0,
    sweeps = run$sweeps,
    seconds = run$seconds,
    kQuantiles = reading$estimate[1:3],
    kSe = reading$se[1:3],
    ends = reading$estimate[-(1:3)],
    endsSe = reading$se[-(1:3)],
    flat = run$flat
  )
  if (build == "move") {
    base <- min(run$essK, run$essMean, summary$essRows)
    own <- list()
    summary$candidates <- lapply(names(candidates), function(name) {
      weights <- run$weights[[name]]
      ess <- base * weights$kish
      needsOwn <- !weights$holds || ess < ownBelow
      out <- list(
        kish = weights$kish,
        khat = weights$khat,
        holds = weights$holds,
        ess = ess,
        needsOwn = needsOwn,
        ownScratch = needsOwn && is.null(candidates[[name]]$own)
      )
      if (needsOwn && !is.null(candidates[[name]]$own)) {
        prior <- candidates[[name]]$own
        ownRun <- stage2Run(
          data,
          dbartsPriors$chi(prior[1L], prior[2L]),
          970000L + rep,
          list(),
          flat
        )
        if (cell %in% stage2SpreadCells) {
          own[[name]] <<- ownRun
        }
        out[["own"]] <- ownRun[c("sweeps", "seconds", "settled")]
      }
      out
    })
    names(summary$candidates) <- names(candidates)
    summary$correction <- c(
      list(C0 = correctionCheck(run, array(1, dim(run$k)))),
      lapply(run$weights, function(w) {
        if (w$holds) correctionCheck(run, w$w) else NULL
      })
    )
    if (cell %in% stage2SpreadCells) {
      summary$blind <- stage2Blind(data, run, own)
    }
  }
  saveRDS(
    summary,
    outPath("stage2", sprintf("%s-%s-%d.rds", build, cell, rep))
  )
  cat(sprintf(
    "%s %s %d (%s): %d sweeps a chain, %.0f s, settled %s\n",
    build,
    cell,
    rep,
    data$band,
    run$sweeps,
    run$seconds,
    run$settled
  ))
}

# The correction's check on a tail count, for the stage 2 ruling: runs as
# the main run would after the ruling (no extension for the candidates) and,
# on the check rows, reads the coverage of the other chains' draws by an
# interval from chains 1 and 2 (all, half and a quarter of their draws)
# against 0.9 N / (N + 1) with N counted two ways: the bulk count of the run
# rule (batchEss of the row's log mean) and a tail count, the smaller of the
# two endpoints' counts, each the effective draws of the indicator of a draw
# below the endpoint (Geyer). Both are times the weights' Kish share.
tailEss <- function(draws, ends) {
  draws <- as.matrix(draws)
  min(vapply(
    ends,
    function(e) {
      sum(apply(draws <= e, 2L, function(z) {
        z <- as.numeric(z)
        tau <- iatGeyer(z)
        if (is.na(tau)) length(z) else length(z) / max(tau, 1)
      }))
    },
    0
  ))
}

correctionCheckTail <- function(run, w) {
  first <- 1:2
  second <- 3:4
  out <- NULL
  w2 <- as.vector(w[, second])
  if (!(sum(w2) > 0)) {
    return(NULL)
  }
  for (share in c(1, 0.5, 0.25)) {
    draws <- seq_len(floor(nrow(run$k) * share))
    w1 <- w[draws, first, drop = FALSE]
    if (!(sum(w1) > 0)) {
      next
    }
    kish <- sum(w1)^2 / sum(w1^2) / length(w1)
    for (i in run$checkRows) {
      x1 <- matrix(run$train[i, draws, first], ncol = 2L)
      ends <- weightedQuantile(as.vector(x1), as.vector(w1), c(0.05, 0.95))
      x2 <- as.vector(run$train[i, , second])
      out <- rbind(
        out,
        c(
          share = share,
          bulk = batchEss(x1) * kish,
          tail = tailEss(x1, ends) * kish,
          cover = sum(w2 * (x2 >= ends[1L] & x2 <= ends[2L])) / sum(w2)
        )
      )
    }
  }
  as.data.frame(out)
}

runCorrection <- function(cell, rep) {
  requireBuild("move")
  data <- stage2Cell(cell, rep)
  flat <- data$band == "narrow"
  run <- stage2Run(
    data,
    NULL,
    980000L + rep,
    stage2Candidates,
    flat,
    extend = FALSE
  )
  weights <- c(
    list(C0 = list(w = array(1, dim(run$k)), kish = 1, holds = TRUE)),
    run$weights
  )
  # Each candidate's counts over the whole run, for the own-run rules: the
  # 10th percentile over the check rows of the bulk and the tail count.
  counts <- lapply(weights, function(w) {
    if (!(w$kish > 0)) {
      return(c(bulk = 0, tail = 0))
    }
    wv <- as.vector(w$w)
    rows <- t(vapply(
      run$checkRows,
      function(i) {
        x <- matrix(run$train[i, , ], ncol = ncol(run$k))
        ends <- weightedQuantile(as.vector(x), wv, c(0.05, 0.95))
        c(bulk = batchEss(x), tail = tailEss(x, ends))
      },
      c(bulk = 0, tail = 0)
    ))
    c(
      bulk = quantile(rows[, "bulk"], 0.1, names = FALSE) * w$kish,
      tail = quantile(rows[, "tail"], 0.1, names = FALSE) * w$kish,
      essK = run$essK * w$kish,
      essMean = run$essMean * w$kish
    )
  })
  checks <- lapply(names(weights), function(name) {
    w <- weights[[name]]
    if (!w$holds) {
      return(NULL)
    }
    d <- correctionCheckTail(run, w$w)
    if (is.null(d)) NULL else cbind(d, prior = name)
  })
  out <- list(
    cell = cell,
    rep = rep,
    band = data$band,
    sweeps = run$sweeps,
    seconds = run$seconds,
    holds = vapply(weights, `[[`, NA, "holds"),
    counts = counts,
    check = do.call(rbind, checks)
  )
  saveRDS(out, outPath("correction", sprintf("%s-%d.rds", cell, rep)))
  cat(sprintf(
    "correction %s %d (%s): %d sweeps a chain, %.0f s\n",
    cell,
    rep,
    data$band,
    run$sweeps,
    run$seconds
  ))
}

summarizeCorrection <- function() {
  files <- list.files(
    file.path(runDir, "correction"),
    "[.]rds$",
    full.names = TRUE
  )
  runs <- lapply(files, readRDS)
  d <- do.call(
    rbind,
    lapply(seq_along(runs), function(j) {
      x <- runs[[j]]$check
      if (is.null(x)) NULL else cbind(x, run = j)
    })
  )
  cuts <- c(0, 25, 50, 75, 100, 150, 200, Inf)
  for (count in c("bulk", "tail")) {
    cat(sprintf("== error of 0.9 N / (N + 1), N the %s count ==\n", count))
    n <- d[[count]]
    d$error <- d$cover - 0.9 * n / (n + 1)
    d$bin <- cut(n, cuts)
    for (group in list(
      C0 = "C0",
      candidates = setdiff(unique(d$prior), "C0")
    )) {
      x <- d[d$prior %in% group, ]
      rows <- lapply(split(x, x$bin), function(b) {
        if (nrow(b) == 0L) {
          return(NULL)
        }
        # Standard error over runs: rows within a run share its draws.
        byRun <- tapply(b$error, b$run, mean)
        data.frame(
          count = levels(x$bin)[b$bin[1L]],
          rows = nrow(b),
          runs = length(byRun),
          error = round(mean(b$error), 4),
          se = round(
            if (length(byRun) > 1L) sd(byRun) / sqrt(length(byRun)) else NA,
            4
          )
        )
      })
      cat(
        " ",
        if (identical(group, "C0")) "C0" else "candidates, reweighted",
        "\n"
      )
      print(do.call(rbind, rows), row.names = FALSE)
    }
  }
  cat("== own runs a replication by rule (candidates whose weights hold) ==\n")
  rules <- do.call(
    rbind,
    lapply(runs, function(r) {
      counts <- r$counts[setdiff(names(r$counts), "C0")]
      holds <- r$holds[names(counts)]
      floor <- vapply(counts, function(cc) min(cc), 0)
      data.frame(
        band = r$band,
        failing = sum(!holds),
        tail50 = sum(
          !holds |
            vapply(
              counts,
              function(cc) {
                min(cc[["tail"]], cc[["essK"]], cc[["essMean"]]) < 50
              },
              NA
            )
        ),
        bulk50 = sum(
          !holds |
            vapply(
              counts,
              function(cc) {
                min(cc[["bulk"]], cc[["essK"]], cc[["essMean"]]) < 50
              },
              NA
            )
        ),
        bulk100 = sum(
          !holds |
            vapply(
              counts,
              function(cc) {
                min(cc[["bulk"]], cc[["essK"]], cc[["essMean"]]) < 100
              },
              NA
            )
        ),
        lowest = round(min(floor), 1)
      )
    })
  )
  print(
    aggregate(cbind(failing, tail50, bulk50, bulk100) ~ band, rules, mean),
    digits = 3
  )
}

# The first chunk at which C0's own target was met: the length the plan's
# rule before the amendment would have run.
lengthC0 <- function(trajectory) {
  met <- which(
    pmin(trajectory$essK, trajectory$essMean, trajectory$essRows) >= longTarget
  )
  if (length(met)) trajectory$sweeps[met[1L]] else NA_real_
}

summarizeStage2 <- function() {
  files <- list.files(
    file.path(runDir, "stage2"),
    "^(move|today)-.*[.]rds$",
    full.names = TRUE
  )
  runs <- lapply(files, readRDS)
  move <- Filter(function(r) r$build == "move", runs)
  today <- Filter(function(r) r$build == "today", runs)

  cat("== The builds' agreement on C0's posterior ==\n")
  zK <- rowsWithin <- NULL
  for (t in today) {
    m <- Filter(function(r) r$cell == t$cell && r$rep == t$rep, move)
    if (length(m) == 0L) {
      next
    }
    m <- m[[1L]]
    z <- (m$kQuantiles - t$kQuantiles) / sqrt(m$kSe^2 + t$kSe^2)
    zRows <- (m$ends - t$ends) / sqrt(m$endsSe^2 + t$endsSe^2)
    both <- m$settledC0 && t$settledC0
    cat(sprintf(
      "  %-12s %-6s k's 5/50/95 move %s, today %s; z %s; rows' ends within 2 se %.3f (max |z| %.1f); sweeps a chain %d, %d; %s\n",
      t$cell,
      t$band,
      paste(sprintf("%.3f", m$kQuantiles), collapse = " "),
      paste(sprintf("%.3f", t$kQuantiles), collapse = " "),
      paste(sprintf("%+.2f", z), collapse = " "),
      mean(abs(zRows) <= 2),
      max(abs(zRows)),
      m$sweeps,
      t$sweeps,
      if (both) "both settled" else "not both settled, not read"
    ))
    if (both) {
      zK <- c(zK, z)
      rowsWithin <- c(rowsWithin, mean(abs(zRows) <= 2))
    }
  }
  if (length(zK)) {
    cat(sprintf(
      "  read: %d z's of k, chi-squared %.1f on %d, p %.3f, max |z| %.2f; rows within 2 se, least %.3f\n",
      length(zK),
      sum(zK^2),
      length(zK),
      pchisq(sum(zK^2), length(zK), lower.tail = FALSE),
      max(abs(zK)),
      min(rowsWithin)
    ))
  }

  cat("\n== Runs: length, time and weights (move build) ==\n")
  table <- do.call(
    rbind,
    lapply(move, function(r) {
      own <- vapply(r$candidates, function(c) c$needsOwn, NA)
      scratch <- vapply(r$candidates, function(c) c$ownScratch, NA)
      ownSeconds <- sum(vapply(
        r$candidates,
        function(c) if (is.null(c[["own"]])) 0 else c[["own"]]$seconds,
        0
      ))
      data.frame(
        cell = r$cell,
        rep = r$rep,
        band = r$band,
        sId = round(r$sId, 2),
        n = r$n,
        sweeps = r$sweeps,
        sweepsC0 = lengthC0(r$trajectory),
        settled = r$settled,
        seconds = round(r$seconds),
        own = paste(which(own), collapse = ""),
        scratch = paste(which(scratch), collapse = ""),
        ownSeconds = round(ownSeconds),
        minKish = round(min(vapply(r$candidates, `[[`, 0, "kish")), 2),
        maxKhat = round(max(vapply(r$candidates, `[[`, 0, "khat")), 2)
      )
    })
  )
  table <- table[
    order(
      match(table$band, c("narrow", "middle", "wide", "small")),
      table$cell,
      table$rep
    ),
  ]
  print(table, row.names = FALSE)
  cat(
    "  own: candidates (1 to 5, in the plan's order) needing their own run; scratch: of those, needing a scratch build\n"
  )
  saveRDS(table, outPath("stage2", "runs-table.rds"))

  cat("\n== The expected-coverage correction ==\n")
  correction <- do.call(
    rbind,
    lapply(move, function(r) {
      do.call(
        rbind,
        lapply(names(r$correction), function(name) {
          d <- r$correction[[name]]
          if (is.null(d) || nrow(d) == 0L) NULL else cbind(d, prior = name)
        })
      )
    })
  )
  correction$bin <- cut(correction$n, c(0, 15, 25, 50, 100, 200, Inf))
  for (group in list(
    C0 = "C0",
    candidates = setdiff(unique(correction$prior), "C0")
  )) {
    d <- correction[correction$prior %in% group, ]
    rows <- lapply(split(d, d$bin), function(b) {
      if (nrow(b) == 0L) {
        return(NULL)
      }
      expected <- 0.9 * b$n / (b$n + 1)
      data.frame(
        effective = levels(d$bin)[b$bin[1L]],
        rows = nrow(b),
        covered = round(mean(b$cover), 4),
        expected = round(mean(expected), 4),
        se = round(sd(b$cover - expected) / sqrt(nrow(b)), 4)
      )
    })
    cat(
      " ",
      if (identical(group, "C0")) "C0" else "candidates, reweighted",
      "\n"
    )
    print(do.call(rbind, rows), row.names = FALSE)
  }

  cat(
    "\n== Paired spread, blind (pooled within-cell sd of a candidate's F coverage less C0's) ==\n"
  )
  spread <- Filter(function(r) !is.null(r$blind), move)
  if (length(spread)) {
    for (bands in list("narrow", "small")) {
      cells <- unique(vapply(
        Filter(function(r) r$band %in% bands, spread),
        `[[`,
        "",
        "cell"
      ))
      sumSquares <- df <- 0
      perCode <- NULL
      for (cell in cells) {
        d <- do.call(
          rbind,
          lapply(Filter(function(r) r$cell == cell, spread), `[[`, "blind")
        )
        if (is.null(d) || nrow(d) < 2L) {
          next
        }
        centered <- sweep(d, 2L, colMeans(d, na.rm = TRUE))
        ok <- colSums(!is.na(d)) >= 2L
        sumSquares <- sumSquares + sum(centered[, ok]^2, na.rm = TRUE)
        df <- df + sum(colSums(!is.na(d[, ok, drop = FALSE])) - 1L)
        perCode <- rbind(perCode, apply(d, 2L, sd, na.rm = TRUE))
      }
      cat(sprintf(
        "  %s: %d cells, pooled sd %.3f on %d df; largest candidate's within-cell sd %.3f\n",
        bands,
        length(cells),
        sqrt(sumSquares / df),
        df,
        max(perCode, na.rm = TRUE)
      ))
    }
  }
}

## ------------------------------------------------- Part B: the main run

# The cells, fixed for the pre-registration: each real dataset but BCI
# (dropped by the amendment after stage 1) with its additive truth as
# fitted, its BART-fitted truth, and its additive truth widened to S_id 2
# and 3 where narrower; the five datasets whose shape is below 1 at that
# shape, the shape-below-1 guard; the small-sample cells; and the designed
# cells. Each real outcome is scored on 5 folds.
mainReal <- c(
  "roaches",
  "solder",
  "aids",
  "mite",
  "Arabidopsis",
  "grouseticks",
  "epil",
  "Insurance",
  "crimtab",
  "Traffic",
  "quakes"
)
mainShapeGuard <- c("roaches", "mite", "Arabidopsis", "grouseticks", "crimtab")
mainCellNames <- c(
  mainReal,
  paste0(mainReal, "-bart"),
  paste0(setdiff(mainReal, "crimtab"), "-w2"),
  paste0(mainReal, "-w3"),
  paste0(mainShapeGuard, "-r"),
  paste0(
    rep(c("roaches", "aids", "epil", "quakes"), each = 4L),
    c("-n30-2", "-n30-3", "-n100-2", "-n100-3")
  ),
  "D1",
  "D2",
  "D3",
  "D4-2",
  "D4-3",
  "D5",
  "D6",
  "D7",
  "D8"
)
mainFoldNames <- paste0(rep(mainReal, each = 5L), "-fold", 1:5)

# Replications by band, from the amendments after stages 1 and 2 and the
# fix round after the Part B review (the small band 10); crimtab's cells 10;
# the folds one repetition. A confirmation block adds 5.
bandReps <- c(narrow = 5L, middle = 5L, wide = 5L, small = 10L)
cellReps <- function(cell, band) {
  if (grepl("-fold", cell, fixed = TRUE)) {
    1L
  } else if (startsWith(cell, "crimtab")) {
    10L
  } else {
    bandReps[[band]]
  }
}
blockReps <- 5L

# The widest cells' guard: crimtab's cells and the cells widened to S_id 3.
widestCells <- c("crimtab", "crimtab-bart", paste0(mainReal, "-w3"))

designedSpecs <- list(
  D1 = list(n = 500L, spread = 1, mean = 0.05, shape = 3, seed = 801000L),
  D2 = list(n = 1000L, spread = 1, mean = 0.2, shape = 2, seed = 802000L),
  D3 = list(n = 300L, spread = 1, mean = 100, shape = 10, seed = 803000L),
  D4 = list(n = 60L, spread = NA, mean = 5, shape = 3, seed = 810000L),
  D5 = list(n = 500L, spread = 0, mean = 5, shape = 5, seed = 850000L),
  D6 = list(n = 5000L, spread = 1.5, mean = 2, shape = 3, seed = 806000L),
  D7 = list(n = 1000L, spread = 1, mean = 5, shape = 3, seed = 807000L),
  D8 = list(n = 1000L, spread = 0.3, mean = 5, shape = 3, seed = 880000L)
)

# A designed cell: covariates, offset and truth, and the truth at either
# value of the binary column b for R (truthAt(1) at b = 0, truthAt(2) at
# b = 1). D4, D5 and D8 draw as designedCell
# did in stage 2. D1 adds an exposure log-uniform over 1 to 1,000, its mean
# the rate times the exposure. D7's treatment depends on the covariates and
# its log rate ratio runs from 0.3 to 1, the rest scaled so that S is 1.
designedBase <- function(name, spread = NULL) {
  spec <- designedSpecs[[name]]
  if (!is.null(spread)) {
    spec$spread <- spread
  }
  set.seed(spec$seed)
  x <- data.frame(correlatedNormals(spec$n, 5L), b = rbinom(spec$n, 1L, 0.5))
  offset <- if (name == "D1") runif(spec$n, 0, log(1000)) else NULL
  if (name == "D7") {
    x$b <- rbinom(spec$n, 1L, plogis(0.5 * (x$x1 - x$x2)))
    prognosis <- with(x, sin(x1) + 0.5 * x2 * x3 + 0.5 * x4^2)
    prognosis <- prognosis - mean(prognosis)
    effect <- 0.3 + 0.7 * plogis(1.5 * x$x5)
    a <- uniroot(
      function(a) sd(a * prognosis + x$b * effect) - spec$spread,
      c(0, 10)
    )$root
    at <- function(b) a * prognosis + b * effect
  } else {
    g <- function(b) {
      with(x, sin(x1) + 0.5 * x2 * x3 + 0.5 * x4^2 + 0.5 * b * x5)
    }
    observed <- g(x$b)
    at <- function(b) {
      if (spec$spread > 0) {
        spec$spread * (g(b) - mean(observed)) / sd(observed)
      } else {
        rep(0, spec$n)
      }
    }
  }
  center <- log(spec$mean) - log(mean(exp(at(x$b))))
  list(
    x = x,
    offset = offset,
    truth = at(x$b) + center,
    shape = spec$shape,
    null = spec$spread == 0,
    rColumn = "b",
    rValues = list(0L, 1L),
    truthAt = function(j) at(j - 1L) + center
  )
}

# A cell's fixed parts: a designed cell, or a real dataset's covariates and
# a truth (as fitted, BART-fitted, widened, at its own shape below 1, or
# widened on a small sample, whose rows are drawn per replication), or its
# real outcome on one fold. On a real dataset with a binary column
# (binaryColumns), the truth at either level of it for R, from
# data/<dataset>-cf.rds: a widened truth's map, affine in the log mean, is
# applied to the counterfactual truth as to the truth.
mainBase <- function(cell) {
  parts <- strsplit(cell, "-", fixed = TRUE)[[1L]]
  if (grepl("^D[0-9]$", parts[1L])) {
    spread <- if (length(parts) == 2L) as.numeric(parts[2L]) else NULL
    return(designedBase(parts[1L], spread))
  }
  d <- readData(parts[1L])
  base <- list(
    x = d$x,
    offset = d$offset,
    truth = as.numeric(d$truthAdditive$eta),
    shape = max(d$truthAdditive$theta, 1),
    null = FALSE
  )
  kind <- if (length(parts) == 1L) "fitted" else sub("[0-9]+$", "", parts[2L])
  if (kind == "bart") {
    path <- file.path(runDir, "data", paste0(parts[1L], "-truth-bart.rds"))
    base$truth <- as.numeric(readRDS(path)$eta)
  } else if (kind == "r") {
    base$shape <- d$truthAdditive$theta
  } else if (kind == "w") {
    target <- as.numeric(sub("w", "", parts[2L]))
    base$truth <- widenTruth(base$truth, base$offset, target)
  } else if (kind == "n") {
    base$truth <- widenTruth(base$truth, base$offset, as.numeric(parts[3L]))
    base$shape <- 3
    base$subsample <- as.integer(sub("n", "", parts[2L]))
  } else if (kind == "fold") {
    base$fold <- as.integer(sub("fold", "", parts[2L]))
    base$y <- d$y
  } else if (kind != "fitted") {
    stop("unknown cell ", cell)
  }
  if (kind != "fold" && parts[1L] %in% names(binaryColumns)) {
    path <- file.path(runDir, "data", paste0(parts[1L], "-cf.rds"))
    cf <- readRDS(path)
    at <- if (kind == "bart") cf$bart else cf$additive
    map <- attr(base$truth, "map")
    if (!is.null(map)) {
      at <- map[["center"]] +
        map[["shift"]] +
        map[["scale"]] *
          (at - map[["center"]])
    }
    base$rColumn <- cf$column
    base$rValues <- lapply(cf$levels, factor, levels = levels(d$x[[cf$column]]))
    base$truthAt <- function(j) at[, j]
  }
  base$truth <- as.numeric(base$truth)
  base
}

# G's groups over a replication's rows: the levels of the first column with
# 2 to 10 distinct values, and the quartiles of the first numeric column with
# more.
mainGroups <- function(x) {
  distinct <- vapply(x, function(v) length(unique(v)), 0L)
  out <- list()
  categorical <- which(distinct >= 2L & distinct <= 10L)
  if (length(categorical)) {
    j <- categorical[1L]
    out[[names(x)[j]]] <- factor(x[[j]])
  }
  continuous <- which(vapply(x, is.numeric, NA) & distinct > 10L)
  if (length(continuous)) {
    j <- continuous[1L]
    q <- quantile(x[[j]], c(0.25, 0.5, 0.75), names = FALSE)
    quartile <- findInterval(x[[j]], q, left.open = TRUE) + 1L
    out[[paste(names(x)[j], "quartile")]] <- factor(quartile)
  }
  out
}

# One replication of a cell. Rows are ordered training then held-out
# ("all rows"); the test rows are the held-out rows and, for a cell with a
# binary column, up to 200 of the rows at its first and at its second value,
# for R. The truth and the groups are carried for the scoring, which is off
# until the pre-registration; no run reads them.
mainCell <- function(cell, rep) {
  index <- match(cell, c(mainCellNames, mainFoldNames))
  if (is.na(index)) {
    stop("unknown cell ", cell)
  }
  base <- mainBase(cell)
  isFold <- !is.null(base$fold)
  n <- nrow(base$x)
  set.seed(5000000L + 1000L * index + if (is.null(base$subsample)) 0L else rep)
  rows <- if (is.null(base$subsample)) {
    seq_len(n)
  } else {
    sample.int(n, base$subsample)
  }
  heldOut <- if (isFold) {
    set.seed(5900000L + match(sub("-fold.*", "", cell), mainReal))
    which(sample(rep_len(1:5, n)) == base$fold)
  } else {
    sample.int(length(rows), round(length(rows) / 5))
  }
  train <- setdiff(seq_along(rows), heldOut)
  all <- rows[c(train, heldOut)]
  rRows <- if (is.null(base$truthAt)) {
    integer()
  } else {
    sort(sample.int(length(all), min(200L, length(all))))
  }
  offsetAll <- if (is.null(base$offset)) NULL else base$offset[all]
  o <- if (is.null(offsetAll)) 0 else offsetAll
  if (isFold) {
    y <- base$y[all]
  } else {
    set.seed(6000000L + 1000L * index + rep)
    mu <- exp(base$truth[all] + o)
    repeat {
      y <- rnbinom(length(all), size = base$shape, mu = mu)
      if (max(y) <= countCap) break
    }
  }
  sId <- if (base$null) 0 else spreadIdentified(base$truth, base$offset)
  identifiedShare <- mean(
    exp(base$truth + if (is.null(base$offset)) 0 else base$offset) >=
      identifiedMean
  )
  # D4 is a small-sample cell at both spreads, as the plan's design table
  # names it, although at S 2 its S_id is 1.80.
  band <- if (base$null || sId < 0.5) {
    "narrow"
  } else if (startsWith(cell, "D4-")) {
    "small"
  } else if (length(rows) <= 130L && sId >= 2 - 1e-6) {
    "small"
  } else if (sId < 1.5) {
    "middle"
  } else {
    "wide"
  }
  numTrain <- length(train)
  numHeld <- length(heldOut)
  xAll <- base$x[all, , drop = FALSE]
  held <- numTrain + seq_len(numHeld)
  xTest <- xAll[held, , drop = FALSE]
  offsetTest <- if (is.null(offsetAll)) NULL else offsetAll[held]
  truthR <- NULL
  if (length(rRows)) {
    x0 <- x1 <- xAll[rRows, , drop = FALSE]
    x0[[base$rColumn]] <- base$rValues[[1L]]
    x1[[base$rColumn]] <- base$rValues[[2L]]
    xTest <- rbind(xTest, x0, x1)
    if (!is.null(offsetAll)) {
      offsetTest <- c(offsetTest, offsetAll[rRows], offsetAll[rRows])
    }
    truthR <- cbind(base$truthAt(1L)[all[rRows]], base$truthAt(2L)[all[rRows]])
  }
  list(
    cell = cell,
    rep = rep,
    band = band,
    sId = sId,
    identifiedShare = identifiedShare,
    shape = base$shape,
    isFold = isFold,
    n = length(all),
    S = if (isFold) NA_real_ else sd(base$truth),
    rColumn = if (length(rRows)) base$rColumn else NA_character_,
    numTrain = numTrain,
    numHeld = numHeld,
    rRows = rRows,
    xTrain = xAll[seq_len(numTrain), , drop = FALSE],
    yTrain = y[seq_len(numTrain)],
    offsetTrain = if (is.null(offsetAll)) {
      NULL
    } else {
      offsetAll[seq_len(numTrain)]
    },
    xTest = xTest,
    offsetTest = offsetTest,
    yTest = y[held],
    offsetAll = offsetAll,
    groups = mainGroups(xAll),
    truth = if (isFold) NULL else base$truth[all],
    truthR = truthR
  )
}

# The Poisson GLM's spread for C7: main effects, the offset in the fit and
# out of the spread; columns constant over the training rows left out.
glmSpreadOf <- function(data) {
  x <- data$xTrain
  x <- x[, vapply(x, function(v) length(unique(v)) > 1L, NA), drop = FALSE]
  o <- if (is.null(data$offsetTrain)) rep(0, nrow(x)) else data$offsetTrain
  d <- data.frame(x, y.count = data$yTrain)
  fit <- suppressWarnings(glm(
    y.count ~ .,
    data = d,
    offset = o,
    family = poisson()
  ))
  sd(as.numeric(predict(fit, type = "link")) - o)
}

# The candidates read on Part B: the five that pass Part A under dec-B502,
# fitted, with own runs where their weights fail or keep too few draws; the
# rest of the list, reported where their weights hold; C6 is the shape's
# weights on C0's run (mainShapeWeights). code orders the code changes for
# step 3's tie: C2 1, C4 2, C5 3.
mainCandidates <- function(glmSpread) {
  chi <- function(df, scale, fitted) {
    list(
      logDensity = function(k) logChiDensity(k, df, scale),
      unbounded = scale > 2,
      fitted = fitted,
      own = if (fitted) list(chi = c(df, scale)) else NULL,
      code = 1L
    )
  }
  floored <- function(kFloor, fitted) {
    list(
      logDensity = function(k) logFloored(k, kFloor),
      unbounded = FALSE,
      fitted = fitted,
      own = if (fitted) list(kprior = paste0("floor:", kFloor)) else NULL,
      code = 2L
    )
  }
  halfNormal <- function(median, fitted) {
    list(
      logDensity = function(k) logHalfNormalSpread(k, median),
      unbounded = TRUE,
      fitted = fitted,
      own = if (fitted) list(kprior = paste0("halfnormal:", median)) else NULL,
      code = 3L
    )
  }
  c7 <- leafScale / (1.7 * max(glmSpread, 0.25) * sqrt(qchisq(0.5, 4)))
  list(
    `C2 chi(4,1.03)` = chi(4, 1.03, TRUE),
    `C2 chi(5,0.9)` = chi(5, 0.9, TRUE),
    `C4 k>=0.6` = floored(0.6, TRUE),
    `C4 k>=0.75` = floored(0.75, TRUE),
    `C5 hn 1.3` = halfNormal(1.3, TRUE),
    `C2 chi(3,1.25)` = chi(3, 1.25, FALSE),
    `C2 chi(7,0.75)` = chi(7, 0.75, FALSE),
    `C2 chi(10,0.62)` = chi(10, 0.62, FALSE),
    `C3 chi(1.5,3)` = chi(1.5, 3, FALSE),
    `C4 k>=1` = floored(1, FALSE),
    `C5 hn 1.57` = halfNormal(1.57, FALSE),
    C7 = chi(4, c7, FALSE)
  )
}
candidateCode <- vapply(mainCandidates(1), `[[`, 0L, "code")
fittedNames <- c(
  "C2 chi(4,1.03)",
  "C2 chi(5,0.9)",
  "C4 k>=0.6",
  "C4 k>=0.75",
  "C5 hn 1.3"
)
fileCode <- function(name) gsub("[^A-Za-z0-9.]", "", name)

# C6's per-draw weights on the drawn shape, against the default's grid
# weights.
mainShapeWeights <- function(shape) {
  j <- match(shape, shapeGrid)
  array((shapeWeightsC6 / shapeWeights)[j], dim(shape))
}

# The own-run rule, from the stage 2 ruling: no extension for the
# candidates; at C0's stop, a fitted candidate whose weights fail, or whose
# effective draws after reweighting (the bulk count, the smaller of batch
# means and Geyer's estimate, times the Kish share) are under 100, gets its
# own run. At 100 or more, 0.9 N / (N + 1) predicted the coverage within
# 0.002.
ownBelowMain <- 100

## the reading of a run: what the scoring needs, computed from all its draws

readProbs <- c(0.025, 0.05, 0.1, 0.25, 0.75, 0.9, 0.95, 0.975)

# Each item's (row's, group's, rate ratio's) weighted quantiles at readProbs
# and weighted mean under every reading's weights, and its effective draws
# (the bulk count, unweighted, times each reading's Kish share). items is
# items x draws, the draws ordered chain by chain.
readItems <- function(items, numChains, weights) {
  m <- nrow(items)
  numDraws <- ncol(items)
  ess <- apply(items, 1L, function(x) batchEss(matrix(x, ncol = numChains)))
  numReadings <- length(weights)
  quantiles <- array(NA_real_, c(m, length(readProbs), numReadings))
  means <- matrix(NA_real_, m, numReadings)
  totals <- vapply(weights, sum, 0)
  for (i in seq_len(m)) {
    x <- items[i, ]
    o <- order(x)
    xo <- x[o]
    for (j in seq_len(numReadings)) {
      w <- weights[[j]]
      cw <- cumsum(w[o]) / totals[j]
      quantiles[i, , j] <- xo[pmin(findInterval(readProbs, cw) + 1L, numDraws)]
      means[i, j] <- sum(w * x) / totals[j]
    }
  }
  out <- lapply(seq_len(numReadings), function(j) {
    w <- weights[[j]]
    kish <- sum(w)^2 / sum(w^2) / length(w)
    list(
      quantiles = matrix(quantiles[,, j], m),
      mean = means[, j],
      ess = ess * kish
    )
  })
  names(out) <- names(weights)
  out
}

# A main-run long run: stage2Run's run without extension, its draws kept as
# a list of chunks (all rows, training then test, offset excluded, rows x
# draws x chains) rather than appended by copy, so that a run near the cap
# holds its draws once. The per-draw mean count over the training rows and
# the check rows' draws, which the stopping rule reads every chunk, are
# kept beside. Draws the same as stage2Run's at the same seed.
mainRun <- function(data, kPrior, seed, candidates, flat) {
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
  arguments <- list(
    data$xTrain,
    data$yTrain,
    test = data$xTest,
    family = "nbinom",
    control = control
  )
  if (!is.null(data$offsetTrain)) {
    arguments$offset <- data$offsetTrain
    arguments$offset.test <- data$offsetTest
  }
  if (!is.null(kPrior)) {
    arguments$leaf.prior <- dbartsPriors$normal(k = kPrior)
  }
  sampler <- do.call(dbarts, arguments)
  offsetTrain <- if (is.null(data$offsetTrain)) 0 else data$offsetTrain
  offsetTest <- if (is.null(data$offsetTest)) 0 else data$offsetTest
  numTrain <- length(data$yTrain)
  checkRows <- sort(sample.int(numTrain, min(200L, numTrain)))
  started <- proc.time()[[3L]]
  sampler$run(longBurn, 1L)
  chunks <- checkChunks <- list()
  k <- shapes <- meanCount <- trajectory <- NULL
  repeat {
    fit <- sampler$run(0L, longChunk %/% longThin)
    train <- fit$train - offsetTrain
    chunks[[length(chunks) + 1L]] <- abind1(train, fit$test - offsetTest)
    checkChunks[[length(checkChunks) + 1L]] <- train[
      checkRows,
      ,
      ,
      drop = FALSE
    ]
    meanCount <- rbind(meanCount, apply(exp(train), c(2L, 3L), mean))
    k <- rbind(k, fit$k)
    shapes <- rbind(shapes, fit$shape)
    rm(fit, train)
    swept <- nrow(k) * longThin
    essK <- batchEss(log(k))
    essMean <- batchEss(meanCount)
    checkDraws <- bindChunks(checkChunks)
    essRows <- quantile(apply(checkDraws, 1L, batchEss), 0.1, names = FALSE)
    rm(checkDraws)
    # Collected each chunk, so that a long run's garbage does not grow with
    # its draws.
    invisible(gc())
    flatHere <- flat || mean(k > 4) > 0.05
    weights <- lapply(candidates, candidateWeights, k = k, flat = flatHere)
    kish <- vapply(weights, `[[`, 0, "kish")
    khat <- vapply(weights, `[[`, 0, "khat")
    holds <- vapply(weights, `[[`, NA, "holds")
    prefixed <- function(v, prefix) {
      if (length(v)) as.list(setNames(v, paste0(prefix, names(v)))) else list()
    }
    row <- c(
      list(
        sweeps = swept,
        seconds = proc.time()[[3L]] - started,
        essK = essK,
        essMean = essMean,
        essRows = essRows
      ),
      prefixed(kish, "kish:"),
      prefixed(khat, "khat:"),
      prefixed(holds, "holds:")
    )
    trajectory <- rbind(trajectory, as.data.frame(row, check.names = FALSE))
    if (min(essK, essMean, essRows) >= longTarget || swept >= longCap) {
      break
    }
  }
  list(
    k = k,
    shape = shapes,
    chunks = chunks,
    meanCount = meanCount,
    weights = weights,
    flat = flatHere,
    trajectory = trajectory,
    essK = essK,
    essMean = essMean,
    essRows = essRows,
    settled = min(essK, essMean, essRows) >= longTarget,
    sweeps = swept,
    seconds = proc.time()[[3L]] - started
  )
}

# Binds chunks (rows x draws x chains) along the draws, into one array
# allocated once; rows and draws, when given, select from the result.
bindChunks <- function(chunks, rows = NULL, draws = NULL) {
  first <- chunks[[1L]]
  if (is.null(rows)) {
    rows <- seq_len(dim(first)[1L])
  }
  sizes <- vapply(chunks, function(x) dim(x)[2L], 0L)
  ends <- cumsum(sizes)
  if (is.null(draws)) {
    draws <- seq_len(ends[length(ends)])
  }
  out <- array(NA_real_, c(length(rows), length(draws), dim(first)[3L]))
  for (i in seq_along(chunks)) {
    inChunk <- which(draws > ends[i] - sizes[i] & draws <= ends[i])
    if (length(inChunk)) {
      out[, inChunk, ] <- chunks[[i]][
        rows,
        draws[inChunk] - (ends[i] - sizes[i]),
        ,
        drop = FALSE
      ]
    }
  }
  out
}

# A default run's draws, one chunk.
chunkedRun <- function(run) {
  run$chunks <- list(abind1(run$train, run$test))
  run$meanCount <- apply(exp(run$train), c(2L, 3L), mean)
  run$train <- run$test <- NULL
  run
}

# The rows read at a time: about 4 million doubles a block, whose working
# copies stay a few hundred MB.
blockRows <- function(run) {
  numDraws <- sum(vapply(run$chunks, function(x) dim(x)[2L], 0L)) *
    dim(run$chunks[[1L]])[3L]
  max(1L, floor(4e6 / numDraws))
}

# Blocks of a run's readings, items in order, bound into one reading.
bindReadings <- function(blocks) {
  blocks <- Filter(Negate(is.null), blocks)
  if (length(blocks) == 0L) {
    return(NULL)
  }
  out <- lapply(names(blocks[[1L]]), function(name) {
    parts <- lapply(blocks, `[[`, name)
    list(
      quantiles = do.call(rbind, lapply(parts, `[[`, "quantiles")),
      mean = unlist(lapply(parts, `[[`, "mean")),
      ess = unlist(lapply(parts, `[[`, "ess"))
    )
  })
  names(out) <- names(blocks[[1L]])
  out
}

# P on held-out rows: per draw the predictive cdf at y and at y - 1 and the
# log density at y, read under each reading as the posterior mean of the
# cdf and the log of the mean density. logMean is rows x draws, the offset
# included.
readPredictive <- function(logMean, shape, y, weights) {
  numRows <- nrow(logMean)
  numDraws <- ncol(logMean)
  y <- rep(y, numDraws)
  size <- rep(as.vector(shape), each = numRows)
  mu <- as.vector(exp(logMean))
  cdf <- matrix(pnbinom(y, size = size, mu = mu), numRows)
  below <- matrix(pnbinom(y - 1, size = size, mu = mu), numRows)
  density <- matrix(dnbinom(y, size = size, mu = mu, log = TRUE), numRows)
  top <- apply(density, 1L, max)
  lapply(weights, function(w) {
    w <- as.vector(w) / sum(w)
    list(
      cdf = as.vector(cdf %*% w),
      below = as.vector(below %*% w),
      logDensity = top + log(as.vector(exp(density - top) %*% w))
    )
  })
}

# One run's readings under each named weight array (draws x chains): F's
# rows, G's groups and R's rate ratios (each row's, then the average over
# R's rows) as items, and P. The draws are read a block of rows at a time
# from the run's chunks; G's and R's average's per-draw sums accumulate
# over the blocks.
readRun <- function(run, data, weights) {
  numChains <- ncol(run$k)
  numDraws <- nrow(run$k) * numChains
  flat <- lapply(weights, as.vector)
  o <- if (is.null(data$offsetAll)) rep(0, data$n) else data$offsetAll
  size <- blockRows(run)
  blocksOf <- function(m) split(seq_len(m), (seq_len(m) - 1L) %/% size)
  indicators <- lapply(data$groups, function(levels) {
    t(model.matrix(~ levels - 1))
  })
  groupSums <- lapply(indicators, function(x) matrix(0, nrow(x), numDraws))
  held <- data$numTrain + seq_len(data$numHeld)
  fBlocks <- pBlocks <- list()
  for (rows in blocksOf(data$n)) {
    f <- matrix(bindChunks(run$chunks, rows), length(rows))
    fBlocks[[length(fBlocks) + 1L]] <- readItems(f, numChains, flat)
    meanCount <- exp(f + o[rows])
    for (name in names(indicators)) {
      groupSums[[name]] <- groupSums[[name]] +
        indicators[[name]][, rows, drop = FALSE] %*% meanCount
    }
    rm(meanCount)
    # P's working copies are several, so its rows go a tenth of a block at
    # a time.
    heldHere <- which(rows %in% held)
    for (part in split(
      heldHere,
      (seq_along(heldHere) - 1L) %/% max(1L, size %/% 10L)
    )) {
      pBlocks[[length(pBlocks) + 1L]] <- readPredictive(
        f[part, , drop = FALSE] + o[rows[part]],
        run$shape,
        data$yTest[rows[part] - data$numTrain],
        weights
      )
    }
    rm(f)
    invisible(gc())
  }
  g <- if (length(indicators)) {
    do.call(
      rbind,
      lapply(names(indicators), function(name) {
        out <- log(groupSums[[name]] / rowSums(indicators[[name]]))
        rownames(out) <- paste0(name, "=", levels(data$groups[[name]]))
        out
      })
    )
  } else {
    NULL
  }
  numR <- length(data$rRows)
  r <- NULL
  if (numR) {
    oR <- if (is.null(data$offsetAll)) rep(0, numR) else o[data$rRows]
    sums <- matrix(0, 2L, numDraws)
    rBlocks <- list()
    for (rows in blocksOf(numR)) {
      r0 <- matrix(bindChunks(run$chunks, data$n + rows), length(rows))
      r1 <- matrix(bindChunks(run$chunks, data$n + numR + rows), length(rows))
      sums[1L, ] <- sums[1L, ] + colSums(exp(r0 + oR[rows]))
      sums[2L, ] <- sums[2L, ] + colSums(exp(r1 + oR[rows]))
      rBlocks[[length(rBlocks) + 1L]] <- readItems(r1 - r0, numChains, flat)
      rm(r0, r1)
      invisible(gc())
    }
    average <- matrix(log(sums[2L, ]) - log(sums[1L, ]), 1L)
    r <- bindReadings(c(rBlocks, list(readItems(average, numChains, flat))))
  }
  classes <- list(
    F = bindReadings(fBlocks),
    G = if (is.null(g)) NULL else readItems(g, numChains, flat),
    R = r
  )
  out <- lapply(names(weights), function(name) {
    list(
      F = classes$F[[name]],
      G = if (is.null(classes$G)) NULL else classes$G[[name]],
      R = if (is.null(classes$R)) NULL else classes$R[[name]],
      P = list(
        cdf = unlist(lapply(pBlocks, function(b) b[[name]]$cdf)),
        below = unlist(lapply(pBlocks, function(b) b[[name]]$below)),
        logDensity = unlist(lapply(pBlocks, function(b) b[[name]]$logDensity))
      )
    )
  })
  names(out) <- names(weights)
  out
}

splitRhat <- function(draws) {
  h <- nrow(draws) %/% 2L
  halves <- cbind(
    draws[seq_len(h), , drop = FALSE],
    draws[h + seq_len(h), , drop = FALSE]
  )
  within <- mean(apply(halves, 2L, var))
  between <- h * var(colMeans(halves))
  sqrt(((h - 1) / h * within + between / h) / within)
}

# Single precision, half the storage of doubles: enough for stored draws.
toSingle <- function(x) {
  list(dim = dim(x), bytes = writeBin(as.numeric(x), raw(), size = 4L))
}
fromSingle <- function(s) {
  array(readBin(s$bytes, "double", n = prod(s$dim), size = 4L), s$dim)
}

# A run's stored draws, 125 a chain evenly spaced (500 in all), so that any
# prior on the spread or the shape can be read later by reweighting; and its
# diagnostics.
storeRun <- function(run, logWeights = NULL) {
  numDraws <- nrow(run$k)
  keep <- unique(round(seq(1, numDraws, length.out = min(125L, numDraws))))
  list(
    k = run$k[keep, , drop = FALSE],
    shape = run$shape[keep, , drop = FALSE],
    eta = toSingle(bindChunks(run$chunks, draws = keep)),
    logWeights = lapply(logWeights, function(w) w[keep, , drop = FALSE]),
    rhatK = splitRhat(log(run$k)),
    rhatMean = splitRhat(run$meanCount),
    essK = run$essK,
    essMean = run$essMean,
    essRows = run$essRows,
    sweeps = run$sweeps,
    seconds = run$seconds,
    settled = run$settled
  )
}

# What the scoring reads beside the readings, carried so that it need not
# rebuild the cell.
cellRecord <- function(data) {
  data[c(
    "cell",
    "rep",
    "band",
    "sId",
    "identifiedShare",
    "shape",
    "isFold",
    "n",
    "S",
    "rColumn",
    "numTrain",
    "numHeld",
    "rRows",
    "yTest",
    "offsetAll",
    "groups",
    "truth",
    "truthR"
  )]
}

cellIndex <- function(cell) match(cell, c(mainCellNames, mainFoldNames))

# The long reading of one replication, on the kprior build with its prior
# unset and the move on (the move build): C0's run to its own target; each
# candidate's reading by weights where they hold and keep 100 effective
# draws; else, for a fitted candidate, its own run (a chi prior's on this
# build, C4's and C5's with NBS_KPRIOR set for the run).
runLong <- function(cell, rep) {
  requireBuild("kprior")
  if (buildLabel != "kprior") {
    stop("the long reading runs with NBS_KPRIOR unset and the move on")
  }
  started <- proc.time()[[3L]]
  data <- mainCell(cell, rep)
  index <- cellIndex(cell)
  flat <- data$band == "narrow"
  candidates <- mainCandidates(glmSpreadOf(data))
  run <- mainRun(data, NULL, 3000000L + 1000L * index + rep, candidates, flat)
  base <- min(run$essK, run$essMean, run$essRows)
  shapeWeights <- mainShapeWeights(run$shape)
  c6 <- local({
    w <- as.vector(shapeWeights)
    kish <- sum(w)^2 / sum(w^2) / length(w)
    khat <- if (length(unique(w)) > 1L) {
      as.numeric(posterior::pareto_khat(w, tail = "right", r_eff = 1))
    } else {
      -Inf
    }
    list(
      w = shapeWeights,
      kish = kish,
      khat = khat,
      holds = kish >= kishFloor && isTRUE(khat <= khatCeiling)
    )
  })
  allWeights <- c(run$weights, list(C6 = c6))
  status <- lapply(allWeights, function(w) {
    ess <- base * w$kish
    list(
      kish = w$kish,
      khat = w$khat,
      holds = w$holds,
      ess = ess,
      reweighted = w$holds && ess >= ownBelowMain
    )
  })
  reweighted <- names(status)[vapply(status, `[[`, NA, "reweighted")]
  readings <- readRun(
    run,
    data,
    c(
      list(C0 = array(1, dim(run$k))),
      lapply(allWeights[reweighted], `[[`, "w")
    )
  )
  source <- c(
    C0 = "C0",
    setNames(rep("weights", length(reweighted)), reweighted)
  )
  stored <- list(
    C0 = storeRun(
      run,
      c(
        lapply(run$weights, function(w) log(w$w)),
        list(C6 = log(shapeWeights))
      )
    )
  )
  # C0's draws are freed before the own runs, so that a worker holds one
  # run's draws at a time.
  run$chunks <- NULL
  invisible(gc())
  for (j in seq_along(fittedNames)) {
    name <- fittedNames[j]
    if (name %in% reweighted) {
      next
    }
    own <- candidates[[name]][["own"]]
    seed <- 20000000L + 1000000L * j + 1000L * index + rep
    if (!is.null(own$kprior)) {
      Sys.setenv(NBS_KPRIOR = own$kprior)
      ownRun <- tryCatch(
        mainRun(data, NULL, seed, list(), flat),
        finally = Sys.unsetenv("NBS_KPRIOR")
      )
    } else {
      prior <- dbartsPriors$chi(own$chi[1L], own$chi[2L])
      ownRun <- mainRun(data, prior, seed, list(), flat)
    }
    readings[[name]] <- readRun(
      ownRun,
      data,
      list(own = array(1, dim(ownRun$k)))
    )$own
    source[[name]] <- "own"
    stored[[name]] <- storeRun(ownRun)
    rm(ownRun)
    invisible(gc())
  }
  out <- list(
    cell = cellRecord(data),
    reading = "long",
    status = status,
    source = source,
    readings = readings,
    stored = stored,
    trajectory = run$trajectory,
    flat = run$flat,
    seconds = proc.time()[[3L]] - started
  )
  saveRDS(out, outPath("main", sprintf("long-%s-%d.rds", cell, rep)))
  cat(sprintf(
    "long %s %d (%s): %d sweeps a chain, %d own runs, %.0f s sampling, %.0f s in all\n",
    cell,
    rep,
    data$band,
    run$sweeps,
    sum(source == "own"),
    sum(vapply(stored, `[[`, 0, "seconds")),
    out$seconds
  ))
}

# The default reading: 4 chains of 500 burn-in and 500 kept at 75 trees.
# C0's on today's sampler; for the guard, a chi candidate's on today's
# sampler, C4's or C5's on the kprior build with the move off and its prior
# set in NBS_KPRIOR.
defaultBurn <- 500L
defaultSamples <- 500L
runDefault <- function(cell, rep, candidate = "C0") {
  own <- if (candidate == "C0") {
    NULL
  } else {
    mainCandidates(1)[[candidate]][["own"]]
  }
  if (is.null(own$kprior)) {
    requireBuild("today")
  } else {
    requireBuild("kprior")
    wanted <- paste0(
      "kprior",
      gsub("[^A-Za-z0-9.]", "", own$kprior),
      "-moveoff"
    )
    if (buildLabel != wanted) {
      stop("this default reading runs as ", wanted, ", not ", buildLabel)
    }
  }
  data <- mainCell(cell, rep)
  set.seed(8000000L + 1000L * cellIndex(cell) + rep)
  control <- dbartsControl(
    n.chains = longChains,
    n.threads = 1L,
    n.burn = defaultBurn,
    n.samples = defaultSamples,
    n.trees = numTrees,
    updateState = FALSE,
    verbose = FALSE
  )
  arguments <- list(
    data$xTrain,
    data$yTrain,
    test = data$xTest,
    family = "nbinom",
    control = control
  )
  if (!is.null(data$offsetTrain)) {
    arguments$offset <- data$offsetTrain
    arguments$offset.test <- data$offsetTest
  }
  if (!is.null(own$chi)) {
    arguments$leaf.prior <- dbartsPriors$normal(
      k = dbartsPriors$chi(own$chi[1L], own$chi[2L])
    )
  }
  started <- proc.time()[[3L]]
  sampler <- do.call(dbarts, arguments)
  fit <- sampler$run()
  offsetTrain <- if (is.null(data$offsetTrain)) 0 else data$offsetTrain
  offsetTest <- if (is.null(data$offsetTest)) 0 else data$offsetTest
  run <- list(
    k = fit$k,
    shape = fit$shape,
    train = fit$train - offsetTrain,
    test = fit$test - offsetTest,
    seconds = proc.time()[[3L]] - started,
    sweeps = defaultBurn + defaultSamples,
    settled = NA
  )
  run <- chunkedRun(run)
  run$essK <- batchEss(log(run$k))
  run$essMean <- batchEss(run$meanCount)
  run$essRows <- NA_real_
  out <- list(
    cell = cellRecord(data),
    reading = "default",
    candidate = candidate,
    seconds = run$seconds,
    readings = readRun(
      run,
      data,
      setNames(list(array(1, dim(run$k))), candidate)
    ),
    stored = setNames(list(storeRun(run)), candidate)
  )
  saveRDS(
    out,
    outPath(
      "main",
      sprintf("default-%s-%s-%d.rds", fileCode(candidate), cell, rep)
    )
  )
  cat(sprintf(
    "default %s %s %d (%s): %.0f s\n",
    candidate,
    cell,
    rep,
    data$band,
    run$seconds
  ))
}

# The cells' table for the pre-registration (band, rows, S and S_id,
# identified share, mean count, shape, G's columns, R's column and rows,
# replications), R's composition by band, and look 1's job queue: the main
# run's replications, numbered from repOffset + 1. Later looks' queues are
# written by the decision (summarize main). Builds no run.
repOffset <- 0L
runPlan <- function() {
  cells <- do.call(
    rbind,
    lapply(c(mainCellNames, mainFoldNames), function(cell) {
      data <- mainCell(cell, repOffset + 1L)
      data.frame(
        cell = cell,
        band = data$band,
        n = data$n,
        S = round(data$S, 2),
        sId = round(data$sId, 2),
        identified = round(data$identifiedShare, 2),
        meanCount = round(mean(c(data$yTrain, data$yTest)), 2),
        total = sum(data$yTrain),
        shape = signif(data$shape, 3),
        G = paste(
          sub(" quartile$", " (quartiles)", names(data$groups)),
          collapse = "; "
        ),
        groups = length(unlist(lapply(data$groups, levels))),
        R = if (is.na(data$rColumn)) "" else data$rColumn,
        rRows = length(data$rRows),
        reps = cellReps(cell, data$band)
      )
    })
  )
  write.table(
    cells,
    outPath("main-cells.txt"),
    quote = FALSE,
    row.names = FALSE,
    sep = "\t"
  )
  print(cells, row.names = FALSE)
  scored <- cells[
    !grepl("-fold", cells$cell) &
      !(cells$cell %in% paste0(mainShapeGuard, "-r")),
  ]
  cat("R's cells by band (the shape-below-1 cells apart):\n")
  for (band in names(bandReps)) {
    x <- scored[scored$band == band, ]
    cat(sprintf(
      "  %s: %d of %d cells: %s\n",
      band,
      sum(x$rRows > 0L),
      nrow(x),
      paste(x$cell[x$rRows > 0L], collapse = ", ")
    ))
  }
  queued <- cells[!grepl("-fold", cells$cell), ]
  folds <- cells[grepl("-fold", cells$cell), ]
  writeQueue(
    1L,
    c(queued$cell, folds$cell),
    c(
      lapply(queued$reps, function(m) repOffset + seq_len(m)),
      lapply(folds$reps, function(m) repOffset + seq_len(m))
    )
  )
}

# A look's job queue: each cell's replications, a long and a C0 default
# reading each.
writeQueue <- function(look, cells, reps) {
  jobs <- unlist(lapply(seq_along(cells), function(i) {
    c(
      sprintf("kprior long %s %d", cells[i], reps[[i]]),
      sprintf("today default %s %d", cells[i], reps[[i]])
    )
  }))
  writeLines(jobs, outPath(sprintf("main-queue-%d.txt", look)))
  cat(sprintf(
    "look %d: %d replications, %d jobs\n",
    look,
    length(jobs) / 2L,
    length(jobs)
  ))
}

# The replications of looks 1 to look, from their queues: cell, rep, look.
readQueues <- function(look) {
  do.call(
    rbind,
    lapply(seq_len(look), function(j) {
      path <- file.path(runDir, sprintf("main-queue-%d.txt", j))
      if (!file.exists(path)) {
        return(NULL)
      }
      jobs <- strsplit(grep(" long ", readLines(path), value = TRUE), " ")
      data.frame(
        cell = vapply(jobs, `[`, "", 3L),
        rep = as.integer(vapply(jobs, `[`, "", 4L)),
        look = j
      )
    })
  )
}

# The guard's queue (step 7): a candidate's default readings on every
# replication the main run's long readings carry, named by its file code.
runGuardPlan <- function(code) {
  candidates <- mainCandidates(1)
  name <- names(candidates)[match(code, fileCode(names(candidates)))]
  if (is.na(name) || !(name %in% fittedNames)) {
    stop("no fitted candidate has the code ", code)
  }
  own <- candidates[[name]][["own"]]
  runBuild <- if (is.null(own$kprior)) {
    "today"
  } else {
    paste0("kprior-off@", own$kprior)
  }
  files <- list.files(file.path(runDir, "main"), "^long-.*[.]rds$")
  cells <- sub("^long-(.*)-([0-9]+)[.]rds$", "\\1 \\2", files)
  jobs <- sprintf("%s default %s %s", runBuild, cells, code)
  writeLines(jobs, outPath(sprintf("guard-queue-%s.txt", code)))
  cat(sprintf("guard %s: %d jobs\n", name, length(jobs)))
}

## ------------------------------------------------- Part B: the scoring

# Turned on by the pre-registration commit. Before it, a quick run
# (NBS_QUICK) exercised the scoring on its toy runs, printing only the shape
# of what it computed, never a value.
scoringOn <- TRUE

# The folded normal's mean at d >= 0, in units of its standard error, and
# its inverse: a cell's distance from nominal corrected for noise, as
# dec-B455 defines it, zero below the mean under no distance.
foldedGrid <- seq(0, 40, by = 0.001)
foldedMean <- sqrt(2 / pi) *
  exp(-foldedGrid^2 / 2) +
  foldedGrid * (1 - 2 * pnorm(-foldedGrid))
correctedDistance <- function(observed, se) {
  t <- observed / se
  out <- ifelse(
    t >= max(foldedMean),
    observed,
    se * approx(foldedMean, foldedGrid, pmax(t, foldedMean[1L]))$y
  )
  out[!is.finite(se) | se == 0] <- observed[!is.finite(se) | se == 0]
  out
}

# One reading's scores on one replication, by class: the 90 percent
# coverage and its expected value (0.9 N / (N + 1) on the long reading, 0.9
# on the default), the 50, 80 and 95 percent coverages, the 90 percent
# interval score and width; RMSE for F; P by randomized PIT, its uniforms v
# the cell's, replication's and row's, shared by every candidate and both
# readings; and the held-out log density. F's and P's rows whose true mean
# is below 0.01 are apart.
scoreReading <- function(reading, cell, expected, v) {
  o <- if (is.null(cell$offsetAll)) 0 else cell$offsetAll
  scoreItems <- function(read, truth, rows = seq_along(truth)) {
    if (is.null(read) || length(rows) == 0L) {
      return(NULL)
    }
    q <- read$quantiles[rows, , drop = FALSE]
    truth <- truth[rows]
    inside <- function(lower, upper) truth >= q[, lower] & truth <= q[, upper]
    n <- read$ess[rows]
    width <- q[, 7L] - q[, 2L]
    score <- width +
      20 * pmax(q[, 2L] - truth, 0) +
      20 * pmax(truth - q[, 7L], 0)
    c(
      cover90 = mean(inside(2L, 7L)),
      expected = mean(if (expected) 0.9 * n / (n + 1) else rep(0.9, length(n))),
      cover50 = mean(inside(4L, 5L)),
      cover80 = mean(inside(3L, 6L)),
      cover95 = mean(inside(1L, 8L)),
      intervalScore = mean(score),
      width = mean(width),
      rmse = sqrt(mean((read$mean[rows] - truth)^2)),
      items = length(rows)
    )
  }
  out <- list()
  held <- cell$numTrain + seq_len(cell$numHeld)
  if (!cell$isFold) {
    identified <- exp(cell$truth + o) >= identifiedMean
    out$F <- scoreItems(reading$F, cell$truth, which(identified))
    out$Flow <- scoreItems(reading$F, cell$truth, which(!identified))
    groupTruth <- unlist(lapply(cell$groups, function(levels) {
      indicator <- t(model.matrix(~ levels - 1))
      as.vector(log(indicator %*% exp(cell$truth + o) / rowSums(indicator)))
    }))
    out$G <- scoreItems(reading$G, groupTruth)
    if (!is.null(cell$truthR)) {
      oR <- if (is.null(cell$offsetAll)) 0 else cell$offsetAll[cell$rRows]
      truthR <- c(
        cell$truthR[, 2L] - cell$truthR[, 1L],
        log(mean(exp(cell$truthR[, 2L] + oR))) -
          log(mean(exp(cell$truthR[, 1L] + oR)))
      )
      numR <- length(cell$rRows)
      perRow <- scoreItems(reading$R, truthR, seq_len(numR))
      average <- scoreItems(reading$R, truthR, numR + 1L)
      out$R <- (perRow + average) / 2
    }
  }
  u <- reading$P$below + v * (reading$P$cdf - reading$P$below)
  covered <- u >= 0.05 & u <= 0.95
  lowP <- if (cell$isFold) {
    rep(FALSE, cell$numHeld)
  } else {
    exp(cell$truth[held] + o[if (length(o) > 1L) held else 1L]) < identifiedMean
  }
  pScore <- function(rows) {
    if (!any(rows)) {
      return(NULL)
    }
    c(
      cover90 = mean(covered[rows]),
      expected = 0.9,
      logDensity = mean(reading$P$logDensity[rows]),
      items = sum(rows)
    )
  }
  out$P <- pScore(!lowP)
  out$Plow <- pScore(lowP)
  out
}

# The randomized PIT's uniforms for a replication's held-out rows.
pitUniform <- function(cell, rep, m) {
  set.seed(9500000L + 1000L * cellIndex(cell) + rep)
  runif(m)
}

# Every reading of the replications queued at looks 1 to look, scored: a
# long data frame of cell, band, replication, reading kind (long or
# default), candidate, class and the measures; and each file's seconds. A
# dataset's folds are one cell, "<dataset>-folds", each fold a replication.
scoreAll <- function(look) {
  queued <- readQueues(look)
  files <- list.files(file.path(runDir, "main"), "[.]rds$", full.names = TRUE)
  costs <- NULL
  scores <- do.call(
    rbind,
    lapply(files, function(file) {
      x <- readRDS(file)
      cell <- x$cell$cell
      rep <- x$cell$rep
      if (!any(queued$cell == cell & queued$rep == rep)) {
        return(NULL)
      }
      costs <<- rbind(
        costs,
        data.frame(
          cell = cell,
          band = x$cell$band,
          rep = rep,
          kind = x$reading,
          candidate = if (x$reading == "long") "C0" else x$candidate,
          seconds = x$seconds
        )
      )
      v <- pitUniform(cell, rep, x$cell$numHeld)
      if (x$cell$isFold) {
        rep <- as.integer(sub(".*-fold", "", cell))
        cell <- sub("-fold[0-9]+$", "-folds", cell)
      }
      do.call(
        rbind,
        lapply(names(x$readings), function(candidate) {
          scores <- scoreReading(
            x$readings[[candidate]],
            x$cell,
            expected = x$reading == "long",
            v = v
          )
          do.call(
            rbind,
            lapply(names(scores), function(class) {
              s <- scores[[class]]
              if (is.null(s)) {
                return(NULL)
              }
              measures <- c(
                cover90 = NA,
                expected = NA,
                cover50 = NA,
                cover80 = NA,
                cover95 = NA,
                intervalScore = NA,
                width = NA,
                rmse = NA,
                logDensity = NA,
                items = NA
              )
              measures[names(s)] <- s
              data.frame(
                cell = cell,
                band = x$cell$band,
                rep = rep,
                kind = x$reading,
                candidate = candidate,
                source = if (x$reading == "long") {
                  x$source[[candidate]]
                } else {
                  "default"
                },
                class = class,
                t(measures),
                check.names = FALSE
              )
            })
          )
        })
      )
    })
  )
  list(scores = scores, costs = costs)
}

bootstrapReps <- 2000L
bootstrapSeed <- 9200000L

# The cells' bootstrap draws: whole replications resampled within a cell,
# by the cell's and replication count's own seed, so that every candidate,
# class and statistic on that cell resamples the same replications.
scoreCellNames <- c(mainCellNames, paste0(mainReal, "-folds"))
bootIndex <- local({
  cache <- new.env()
  function(cell, m) {
    key <- paste(cell, m)
    if (is.null(cache[[key]])) {
      set.seed(bootstrapSeed + 1000L * match(cell, scoreCellNames) + m)
      cache[[key]] <- matrix(
        sample.int(m, m * bootstrapReps, replace = TRUE),
        m
      )
    }
    cache[[key]]
  }
})

# The rows of one candidate's scores in one class and reading kind over the
# given cells, paired with C0's on the replications both carry, split by
# cell, each in replication order.
pairScores <- function(scores, candidate, class, kind, cells) {
  pick <- function(who) {
    s <- scores[
      scores$candidate == who &
        scores$class == class &
        scores$kind == kind &
        scores$cell %in% cells,
    ]
    s[order(s$cell, s$rep), ]
  }
  a <- pick(candidate)
  b <- pick("C0")
  common <- intersect(paste(a$cell, a$rep), paste(b$cell, b$rep))
  a <- a[paste(a$cell, a$rep) %in% common, ]
  b <- b[paste(b$cell, b$rep) %in% common, ]
  if (nrow(a) == 0L) {
    return(NULL)
  }
  lapply(split(seq_len(nrow(a)), a$cell), function(i) {
    list(a = a[i, ], b = b[i, ], band = a$band[i[1L]], m = length(i))
  })
}

# A cell's corrected gap, dec-B455: the folded normal's correction of
# |mean miss| at standard error sd / sqrt(m); on the bootstrap's
# replications when idx is given (one column a draw).
cellGap <- function(miss, idx = NULL) {
  m <- length(miss)
  if (is.null(idx)) {
    se <- if (m > 1L) sd(miss) / sqrt(m) else NA_real_
    return(correctedDistance(abs(mean(miss)), se))
  }
  x <- matrix(miss[idx], m)
  mu <- colMeans(x)
  se <- if (m > 1L) {
    sqrt(pmax(colSums(x^2) - m * mu^2, 0) / (m - 1L) / m)
  } else {
    rep(NA_real_, length(mu))
  }
  correctedDistance(abs(mu), se)
}

# The cells' weights in an average: equal over the bands present, cells
# equal within a band ("band"); or equal over cells ("cells").
cellWeights <- function(bands, how = "band") {
  if (how == "cells") {
    return(rep(1 / length(bands), length(bands)))
  }
  size <- table(bands)
  1 / (length(size) * as.vector(size[bands]))
}

# A statistic: its estimate, bootstrap draws and standard error, and its
# cells' weighted contributions to the draws (cells x draws, linearized for
# a ratio), from which the block allocation reads each band's share of its
# variance.
makeStatistic <- function(estimate, draws, contributions, bands, reps) {
  list(
    estimate = estimate,
    draws = draws,
    se = sd(draws),
    contributions = contributions,
    bands = bands,
    reps = reps
  )
}

# The candidate's corrected gap less C0's, averaged over the cells by
# weights; per cell, beside it, the point difference and its bound's
# standard error: the per-replication difference of the signed misses,
# sign(mean) times miss, which linearizes the gap, sd / sqrt(m).
gapStatistic <- function(pair, how = "band") {
  cells <- names(pair)
  w <- cellWeights(vapply(pair, `[[`, "", "band"), how)
  per <- lapply(pair, function(x) {
    missA <- x$a$cover90 - x$a$expected
    missB <- x$b$cover90 - x$b$expected
    idx <- bootIndex(x$a$cell[1L], x$m)
    linear <- sign(mean(missA)) * missA - sign(mean(missB)) * missB
    list(
      point = cellGap(missA) - cellGap(missB),
      draws = cellGap(missA, idx) - cellGap(missB, idx),
      boundSe = if (x$m > 1L) sd(linear) / sqrt(x$m) else NA_real_
    )
  })
  point <- vapply(per, `[[`, 0, "point")
  draws <- t(vapply(per, `[[`, numeric(bootstrapReps), "draws"))
  contributions <- w * draws
  out <- makeStatistic(
    sum(w * point),
    colSums(contributions),
    contributions,
    vapply(pair, `[[`, "", "band"),
    vapply(pair, `[[`, 0L, "m")
  )
  out$cells <- data.frame(
    cell = cells,
    band = out$bands,
    m = out$reps,
    value = point,
    boundSe = vapply(per, `[[`, 0, "boundSe"),
    row.names = NULL
  )
  out
}

# The ratio of the band-weighted averages of a measure, candidate over C0.
ratioStatistic <- function(pair, measure) {
  w <- cellWeights(vapply(pair, `[[`, "", "band"))
  means <- lapply(pair, function(x) {
    idx <- bootIndex(x$a$cell[1L], x$m)
    a <- x$a[[measure]]
    b <- x$b[[measure]]
    list(
      a = mean(a),
      b = mean(b),
      aDraws = colMeans(matrix(a[idx], x$m)),
      bDraws = colMeans(matrix(b[idx], x$m))
    )
  })
  a <- vapply(means, `[[`, 0, "a")
  b <- vapply(means, `[[`, 0, "b")
  aDraws <- t(vapply(means, `[[`, numeric(bootstrapReps), "aDraws"))
  bDraws <- t(vapply(means, `[[`, numeric(bootstrapReps), "bDraws"))
  ratio <- sum(w * a) / sum(w * b)
  makeStatistic(
    ratio,
    colSums(w * aDraws) / colSums(w * bDraws),
    w * (aDraws - ratio * bDraws) / sum(w * b),
    vapply(pair, `[[`, "", "band"),
    vapply(pair, `[[`, 0L, "m")
  )
}

# The log density per row on the real outcomes, C0's less the candidate's,
# averaged over the datasets, each dataset's folds its replications.
logScoreStatistic <- function(pair) {
  w <- cellWeights(vapply(pair, `[[`, "", "band"), "cells")
  per <- lapply(pair, function(x) {
    idx <- bootIndex(x$a$cell[1L], x$m)
    d <- x$b$logDensity - x$a$logDensity
    list(point = mean(d), draws = colMeans(matrix(d[idx], x$m)))
  })
  draws <- t(vapply(per, `[[`, numeric(bootstrapReps), "draws"))
  makeStatistic(
    sum(w * vapply(per, `[[`, 0, "point")),
    colSums(w * draws),
    w * draws,
    vapply(pair, `[[`, "", "band"),
    vapply(pair, `[[`, 0L, "m")
  )
}

# The width guard's cells: per narrow cell, F's interval score ratio,
# candidate over C0, its standard error by the delta method on the
# replications, sd(IS_c - ratio IS_0) / (sqrt(m) mean(IS_0)).
widthCells <- function(pair) {
  do.call(
    rbind,
    lapply(pair, function(x) {
      a <- x$a$intervalScore
      b <- x$b$intervalScore
      ratio <- mean(a) / mean(b)
      data.frame(
        cell = x$a$cell[1L],
        band = x$band,
        m = x$m,
        value = ratio,
        boundSe = if (x$m > 1L) {
          sd(a - ratio * b) / (sqrt(x$m) * mean(b))
        } else {
          NA_real_
        }
      )
    })
  )
}

# The group-sequential boundaries: two-sided at 0.05, symmetric, by
# Lan-DeMets spending of O'Brien-Fleming type on each side, alpha_1(t) =
# 2 - 2 Phi(z_0.9875 / sqrt(t)), at the information fractions t of the looks
# so far (the final look's 1). Each look's crossing probability is
# integrated exactly in the new increment over the previous look's
# continuation density, carried by Simpson's rule on its region. At equal
# thirds: 3.71, 2.51, 1.99.
ldBoundaries <- function(t, alpha = 0.05) {
  spend <- 2 * (2 - 2 * pnorm(qnorm(1 - alpha / 4) / sqrt(t)))
  out <- numeric(length(t))
  out[1L] <- qnorm(1 - spend[1L] / 2)
  simpson <- function(b, m = 2000L) {
    u <- seq(-b, b, length.out = m + 1L)
    w <- c(1, rep(c(4, 2), length.out = m - 1L), 1) * (2 * b / m) / 3
    list(u = u, w = w)
  }
  g <- simpson(out[1L])
  density <- dnorm(g$u)
  for (j in seq_along(t)[-1L]) {
    delta <- t[j] - t[j - 1L]
    target <- spend[j] - spend[j - 1L]
    mean <- sqrt(t[j - 1L]) * g$u
    crossing <- function(b) {
      above <- pnorm((b * sqrt(t[j]) - mean) / sqrt(delta), lower.tail = FALSE)
      below <- pnorm((-b * sqrt(t[j]) - mean) / sqrt(delta))
      sum(g$w * density * (above + below))
    }
    out[j] <- if (target <= 0) {
      Inf
    } else {
      uniroot(
        function(b) crossing(b) - target,
        c(0.1, 20),
        tol = 1e-10
      )$root
    }
    region <- simpson(min(out[j], 20))
    kernel <- dnorm(outer(sqrt(t[j]) * region$u, mean, `-`) / sqrt(delta)) *
      sqrt(t[j] / delta)
    density <- as.vector(kernel %*% (g$w * density))
    g <- region
  }
  out
}

# The decision rule's constants: the bars, the looks, the stop. The stop is
# in the sizing table's core-hours, the laptop's; the box's measured seconds
# are converted by the calibration (calibrate) written to the run directory.
scoreClasses <- c("F", "G", "R", "P")
gapBar <- 0.02
tieMargin <- 0.005
logScoreBar <- 0.01
intervalScoreBar <- 1.05
widthBar <- 1.1
lastLook <- 3L
mainStop <- 110
speedFactor <- function() {
  path <- file.path(runDir, "calibration.rds")
  if (file.exists(path)) readRDS(path)$factor else 1
}

# A one-sided statistic against its bar at a look's boundary b: "pass" when
# (bar - estimate) / se is at least b, "fail" at -b or below, else
# "pending".
sequentialStatus <- function(statistic, bar, b) {
  z <- (bar - statistic$estimate) / statistic$se
  if (!is.finite(z)) {
    z <- if (statistic$estimate < bar) Inf else -Inf
  }
  list(
    z = z,
    status = if (z >= b) {
      "pass"
    } else if (z <= -b) {
      "fail"
    } else {
      "pending"
    }
  )
}

# Steps 2 and 3 on one look's scores. Step 2, per fitted candidate, on the
# long reading: each class's gap less C0's, the widest cells' gap and each
# class's interval score ratio (but P's) against the look's boundary, both
# sides; the log score on the folds, which no block adds to, at its 95
# percent bound; no cell worse by the bar and none of the shape-below-1
# cells (a cell is worse when its gap difference less the bar exceeds
# t_0.975 on m - 1 df times its bound's standard error); and each narrow
# cell's F interval score ratio at most 1.10 at the same t bound, failing
# when its lower bound exceeds it. Step 3: per pair of candidates not
# failing step 2, the lead in the rank (the mean of F's and R's gaps less
# C0's) less the tie margin against the boundary: untied at b, tied at -b.
lookStatistics <- function(scores, b) {
  truthCells <- setdiff(mainCellNames, paste0(mainShapeGuard, "-r"))
  foldCells <- paste0(mainReal, "-folds")
  narrowCells <- unique(scores$cell[scores$band == "narrow"])
  classCells <- function(class) {
    if (class == "P") c(truthCells, foldCells) else truthCells
  }
  tBound <- function(cells, bar, side = "upper") {
    q <- qt(0.975, pmax(cells$m - 1L, 1L))
    if (side == "upper") {
      cells$value + q * cells$boundSe
    } else {
      cells$value - q * cells$boundSe
    }
  }
  one <- function(candidate) {
    pairs <- lapply(scoreClasses, function(class) {
      pairScores(scores, candidate, class, "long", classCells(class))
    })
    names(pairs) <- scoreClasses
    pairs <- Filter(Negate(is.null), pairs)
    gaps <- lapply(pairs, gapStatistic)
    widest <- lapply(names(pairs), function(class) {
      p <- pairScores(scores, candidate, class, "long", widestCells)
      if (is.null(p)) NULL else gapStatistic(p, "cells")
    })
    names(widest) <- names(pairs)
    widest <- Filter(Negate(is.null), widest)
    shape <- lapply(names(pairs), function(class) {
      p <- pairScores(
        scores,
        candidate,
        class,
        "long",
        paste0(mainShapeGuard, "-r")
      )
      if (is.null(p)) NULL else gapStatistic(p, "cells")$cells
    })
    intervalScore <- lapply(
      pairs[intersect(names(pairs), c("F", "G", "R"))],
      ratioStatistic,
      measure = "intervalScore"
    )
    widthPair <- pairScores(scores, candidate, "F", "long", narrowCells)
    width <- if (is.null(widthPair)) NULL else widthCells(widthPair)
    logPair <- pairScores(scores, candidate, "P", "long", foldCells)
    logScore <- if (is.null(logPair)) NULL else logScoreStatistic(logPair)
    sequential <- c(
      lapply(gaps, sequentialStatus, bar = gapBar, b = b),
      setNames(
        lapply(widest, sequentialStatus, bar = gapBar, b = b),
        paste0("widest ", names(widest))
      ),
      setNames(
        lapply(intervalScore, sequentialStatus, bar = intervalScoreBar, b = b),
        paste0("interval score ", names(intervalScore))
      )
    )
    statistics <- c(
      gaps,
      setNames(widest, paste0("widest ", names(widest))),
      setNames(intervalScore, paste0("interval score ", names(intervalScore)))
    )
    worse <- function(cells) {
      if (is.null(cells)) {
        return(character())
      }
      flagged <- cells$value - gapBar >
        qt(0.975, pmax(cells$m - 1L, 1L)) *
          cells$boundSe
      cells$cell[!is.na(flagged) & flagged]
    }
    cellsWorse <- unique(unlist(lapply(gaps, function(g) worse(g$cells))))
    shapeWorse <- unique(unlist(lapply(shape, worse)))
    widthStatus <- if (is.null(width)) {
      character()
    } else {
      upper <- tBound(width, widthBar)
      lower <- tBound(width, widthBar, "lower")
      status <- rep("pending", nrow(width))
      status[!is.na(upper) & upper <= widthBar] <- "pass"
      status[!is.na(lower) & lower > widthBar] <- "fail"
      setNames(status, width$cell)
    }
    logStatus <- if (is.null(logScore)) {
      "pass"
    } else if (logScore$estimate + qnorm(0.975) * logScore$se < logScoreBar) {
      "pass"
    } else {
      "fail"
    }
    states <- c(
      vapply(sequential, `[[`, "", "status"),
      widthStatus,
      log = logStatus,
      cells = if (length(cellsWorse)) "fail" else "pass",
      shape = if (length(shapeWorse)) "fail" else "pass"
    )
    stepTwo <- if (any(states == "fail")) {
      "fail"
    } else if (all(states == "pass")) {
      "pass"
    } else {
      "pending"
    }
    rankContributions <- function(class) {
      if (is.null(gaps[[class]])) NULL else gaps[[class]]$contributions / 2
    }
    list(
      statistics = statistics,
      sequential = sequential,
      width = width,
      widthStatus = widthStatus,
      logScore = logScore,
      logStatus = logStatus,
      cellsWorse = cellsWorse,
      shapeWorse = shapeWorse,
      stepTwo = stepTwo,
      rank = mean(c(gaps$F$estimate, gaps$R$estimate)),
      rankContributions = list(
        F = rankContributions("F"),
        R = rankContributions("R")
      ),
      narrowScore = if (is.null(width)) Inf else mean(width$value)
    )
  }
  verdicts <- lapply(fittedNames, one)
  names(verdicts) <- fittedNames
  # Step 3's pairs, the better ranked first.
  open <- fittedNames[vapply(verdicts, `[[`, "", "stepTwo") != "fail"]
  pairs <- list()
  for (c1 in open) {
    for (c2 in open) {
      if (
        c1 == c2 ||
          verdicts[[c1]]$rank > verdicts[[c2]]$rank ||
          (verdicts[[c1]]$rank == verdicts[[c2]]$rank && c1 > c2)
      ) {
        next
      }
      contributions <- NULL
      for (class in c("F", "R")) {
        x <- verdicts[[c2]]$rankContributions[[class]]
        y <- verdicts[[c1]]$rankContributions[[class]]
        if (is.null(x) || is.null(y)) {
          next
        }
        rownames(x) <- paste(class, rownames(x))
        rownames(y) <- paste(class, rownames(y))
        common <- intersect(rownames(x), rownames(y))
        contributions <- rbind(
          contributions,
          x[common, , drop = FALSE] -
            y[common, , drop = FALSE]
        )
      }
      bands <- c(
        verdicts[[c1]]$statistics$F$bands,
        verdicts[[c1]]$statistics$R$bands
      )
      names(bands) <- c(
        paste("F", names(verdicts[[c1]]$statistics$F$bands)),
        paste("R", names(verdicts[[c1]]$statistics$R$bands))
      )
      reps <- c(
        verdicts[[c1]]$statistics$F$reps,
        verdicts[[c1]]$statistics$R$reps
      )
      names(reps) <- names(bands)
      lead <- makeStatistic(
        verdicts[[c2]]$rank - verdicts[[c1]]$rank,
        colSums(contributions),
        contributions,
        bands[rownames(contributions)],
        reps[rownames(contributions)]
      )
      z <- (lead$estimate - tieMargin) / lead$se
      pairs[[paste(c1, c2, sep = " | ")]] <- list(
        better = c1,
        worse = c2,
        lead = lead,
        z = z,
        status = if (z >= b) {
          "untied"
        } else if (z <= -b) {
          "tied"
        } else {
          "pending"
        }
      )
    }
  }
  list(verdicts = verdicts, pairs = pairs)
}

# Stage 1's results that the decision's code assumes: the default fails A,
# so step 4 never applies, and no candidate passes A with margin, so step
# 3's first tie-break never applies. Neither path is coded, and the decision
# stops if either constant says otherwise.
defaultPassesA <- FALSE
passesAWithMargin <- character()
checkUncodedPaths <- function(
  defaultPasses = defaultPassesA,
  withMargin = passesAWithMargin
) {
  if (!identical(defaultPasses, FALSE)) {
    stop(
      "the default passes A, so step 4 applies; ",
      "the decision does not code step 4"
    )
  }
  if (length(withMargin)) {
    stop(
      "a candidate passes A with margin (",
      paste(withMargin, collapse = ", "),
      "), so step 3's first tie-break applies; the decision does not code it"
    )
  }
  invisible(TRUE)
}

# The guard's constructed cases: stage 1's values pass; the default passing
# A, a candidate passing it with margin, or both, each stop.
runUncodedCheck <- function() {
  cases <- list(
    `stage 1's` = list(FALSE, character(), TRUE),
    `default passes A` = list(TRUE, character(), FALSE),
    `one with margin` = list(FALSE, "C4 k>=0.6", FALSE),
    both = list(TRUE, "C2 chi(4,1.03)", FALSE)
  )
  for (name in names(cases)) {
    x <- cases[[name]]
    passes <- tryCatch(
      checkUncodedPaths(x[[1L]], x[[2L]]),
      error = function(e) {
        FALSE
      }
    )
    if (!identical(passes, x[[3L]])) {
      stop("the uncoded paths' guard fails its case: ", name)
    }
    cat(sprintf(
      "%s: %s, as it should\n",
      name,
      if (isTRUE(passes)) "passes" else "stops"
    ))
  }
}

# Step 3's choice among the candidates passing step 2: the smallest rank;
# within the tie margin of it (by the pairs' relations), the smaller code
# change, then the smaller interval score on the narrow band. The tie-break
# for passing A with margin is not coded (checkUncodedPaths).
chooseCandidate <- function(passing, tied, verdicts) {
  if (length(passing) == 0L) {
    return("C0")
  }
  rank <- vapply(verdicts[passing], `[[`, 0, "rank")
  best <- passing[which.min(rank)]
  ties <- passing[vapply(
    passing,
    function(d) {
      d == best ||
        isTRUE(tied[paste(best, d, sep = " | ")]) ||
        isTRUE(tied[paste(d, best, sep = " | ")])
    },
    NA
  )]
  code <- candidateCode[ties]
  ties <- ties[code == min(code)]
  narrow <- vapply(verdicts[ties], `[[`, 0, "narrowScore")
  ties[which.min(narrow)]
}

# Step 3's order: the choice, then the choice among the rest, and so on.
candidateOrder <- function(passing, tied, verdicts) {
  out <- character()
  while (length(passing)) {
    chosen <- chooseCandidate(passing, tied, verdicts)
    out <- c(out, chosen)
    passing <- setdiff(passing, chosen)
  }
  out
}

# The recommendation under every way the pending results could resolve:
# a pending step 2 as a pass or a fail, a pending pair as tied or untied. At
# the final look, pending reads as a fail and as tied.
enumerateDecisions <- function(look, final) {
  verdicts <- look$verdicts
  stepTwo <- vapply(verdicts, `[[`, "", "stepTwo")
  relation <- vapply(look$pairs, `[[`, "", "status")
  if (final) {
    stepTwo[stepTwo == "pending"] <- "fail"
    relation[relation == "pending"] <- "tied"
  }
  pendingTwo <- names(stepTwo)[stepTwo == "pending"]
  pendingPairs <- names(relation)[relation == "pending"]
  grid <- as.matrix(expand.grid(
    rep(list(c(TRUE, FALSE)), length(pendingTwo) + length(pendingPairs))
  ))
  outcomes <- vapply(
    seq_len(max(1L, nrow(grid))),
    function(i) {
      passes <- stepTwo == "pass"
      tied <- relation == "tied"
      if (length(pendingTwo)) {
        passes[pendingTwo] <- grid[i, seq_along(pendingTwo)]
      }
      if (length(pendingPairs)) {
        tied[pendingPairs] <- grid[
          i,
          length(pendingTwo) + seq_along(pendingPairs)
        ]
      }
      chooseCandidate(names(stepTwo)[passes], tied, verdicts)
    },
    ""
  )
  passes <- stepTwo == "pass"
  tied <- relation == "tied"
  list(
    outcomes = sort(unique(outcomes)),
    resolved = length(unique(outcomes)) == 1L,
    order = if (all(stepTwo != "pending") && all(relation != "pending")) {
      candidateOrder(names(stepTwo)[passes], tied, verdicts)
    } else {
      NULL
    }
  )
}

# The information fraction of each look so far: the variance of a
# band-weighted average over the truth cells at the planned maximum (each
# cell's look-1 replications and two blocks), over its variance at that
# look. A cell's variance a replication is m times its gap difference's
# bootstrap variance, averaged over the classes and fitted candidates.
informationFractions <- function(look, queued, scores) {
  truthCells <- setdiff(mainCellNames, paste0(mainShapeGuard, "-r"))
  perRep <- NULL
  for (candidate in fittedNames) {
    for (class in scoreClasses) {
      p <- pairScores(scores, candidate, class, "long", truthCells)
      if (is.null(p)) {
        next
      }
      g <- gapStatistic(p, "cells")
      perRep <- rbind(
        perRep,
        data.frame(
          cell = g$cells$cell,
          v = g$reps * apply(g$contributions / (1 / nrow(g$cells)), 1L, var)
        )
      )
    }
  }
  v <- tapply(perRep$v, perRep$cell, mean)
  # A cell without a spread (one replication, or none differing) takes the
  # others' average; with none, every cell weighs the same.
  known <- is.finite(v) & v > 0
  v[!known] <- if (any(known)) mean(v[known]) else 1
  cells <- names(v)
  bands <- vapply(
    cells,
    function(cell) scores$band[match(cell, scores$cell)],
    ""
  )
  w <- cellWeights(bands)
  count <- function(j) {
    x <- queued[queued$look <= j & queued$cell %in% cells, ]
    as.vector(table(factor(x$cell, levels = cells)))
  }
  first <- count(1L)
  variance <- function(n) sum(w^2 * v / n)
  maximum <- variance(first + 2L * blockReps)
  vapply(seq_len(look), function(j) maximum / variance(count(j)), 0)
}

# Which bands get the next look's block: per band, the variance a block of
# blockReps replications would remove from each pending statistic, as a
# share of that statistic's variance, summed, per core-hour of the block
# (the band's cells' measured seconds a replication); every band within a
# quarter of the best band's rate, in that order, while the stop allows,
# the last band cut to its cheapest cells that fit. everything counts every
# statistic as pending (a quick run's, to exercise it).
allocateBlock <- function(look, costs, queued, everything = FALSE) {
  pending <- list()
  for (candidate in names(look$verdicts)) {
    v <- look$verdicts[[candidate]]
    if (v$stepTwo == "fail" && !everything) {
      next
    }
    for (name in names(v$sequential)) {
      if (v$sequential[[name]]$status == "pending" || everything) {
        pending[[length(pending) + 1L]] <- v$statistics[[name]]
      }
    }
    for (cell in names(v$widthStatus)[
      v$widthStatus == "pending" | everything
    ]) {
      x <- v$width[v$width$cell == cell, ]
      contributions <- matrix(x$boundSe, 1L, 2L, dimnames = list(cell, NULL))
      pending[[length(pending) + 1L]] <- list(
        contributions = contributions,
        bands = setNames(x$band, cell),
        reps = setNames(x$m, cell),
        variance = x$boundSe^2
      )
    }
  }
  for (pair in look$pairs) {
    if (pair$status == "pending" || everything) {
      pending[[length(pending) + 1L]] <- pair$lead
    }
  }
  scoreBands <- setNames(numeric(length(bandReps)), names(bandReps))
  for (x in pending) {
    perCell <- if (is.null(x$variance)) {
      apply(x$contributions, 1L, var)
    } else {
      x$variance
    }
    perCell[is.na(perCell)] <- 0
    total <- sum(perCell)
    if (!(total > 0)) {
      next
    }
    cells <- sub("^[FR] ", "", names(x$bands))
    grows <- !grepl("-folds$", cells)
    gain <- perCell * blockReps / (x$reps + blockReps) * grows / total
    for (band in names(scoreBands)) {
      scoreBands[band] <- scoreBands[band] + sum(gain[x$bands == band])
    }
  }
  laptop <- costs$seconds / speedFactor() / 3600
  spent <- sum(laptop)
  main <- costs$candidate == "C0"
  perRep <- tapply(laptop[main], costs$cell[main], sum) /
    tapply(costs$kind[main] == "long", costs$cell[main], sum)
  cellBand <- tapply(costs$band, costs$cell, `[`, 1L)[names(perRep)]
  growing <- !grepl("-fold", names(perRep))
  blockCost <- vapply(
    names(bandReps),
    function(band) blockReps * sum(perRep[growing & cellBand == band]),
    0
  )
  rate <- scoreBands / blockCost
  rate[!is.finite(rate)] <- 0
  chosen <- names(sort(
    rate[rate > 0 & rate >= max(rate) / 4],
    decreasing = TRUE
  ))
  funded <- list()
  for (band in chosen) {
    cells <- names(perRep)[growing & cellBand == band]
    cost <- blockReps * perRep[cells]
    if (spent + sum(cost) <= mainStop) {
      funded[[band]] <- cells
      spent <- spent + sum(cost)
      next
    }
    cells <- cells[order(cost[cells])]
    fits <- cells[spent + cumsum(cost[cells]) <= mainStop]
    if (length(fits)) {
      funded[[band]] <- fits
      spent <- spent + sum(cost[fits])
    }
    break
  }
  list(
    scores = scoreBands,
    costs = blockCost,
    rates = rate,
    funded = funded,
    spent = sum(laptop),
    projected = spent
  )
}

# Step 5: the order's first candidate whose default reading is not worse
# than C0's by the bar in any class (each class's gap less C0's below 0.02
# at its 95 percent bound, read once, on every replication run); a
# candidate whose default readings are missing waits for its guard queue.
stepFive <- function(scores, order, queued) {
  long <- unique(queued[c("cell", "rep")])
  for (name in order) {
    code <- fileCode(name)
    have <- scores[scores$kind == "default" & scores$candidate == name, ]
    if (nrow(unique(have[c("cell", "rep")])) < nrow(long)) {
      return(list(chosen = NA_character_, waiting = name, checked = list()))
    }
    gaps <- lapply(scoreClasses, function(class) {
      truthCells <- setdiff(mainCellNames, paste0(mainShapeGuard, "-r"))
      cells <- if (class == "P") {
        c(truthCells, paste0(mainReal, "-folds"))
      } else {
        truthCells
      }
      p <- pairScores(scores, name, class, "default", cells)
      if (is.null(p)) NULL else gapStatistic(p)
    })
    names(gaps) <- scoreClasses
    gaps <- Filter(Negate(is.null), gaps)
    upper <- vapply(gaps, function(g) g$estimate + qnorm(0.975) * g$se, 0)
    if (all(upper < gapBar)) {
      return(list(chosen = name, waiting = NA_character_, upper = upper))
    }
  }
  list(chosen = "C0", waiting = NA_character_)
}

# The decision at a look: the boundary from the information fractions, steps
# 2 and 3, whether every way the pending results could resolve gives one
# recommendation; if not and a block can run, the next look's queue; once
# decided, step 5, writing the guard queue it waits for. Saved to
# main-decision-<look>.rds. Scoring is off until the pre-registration: a
# quick run computes all of it and prints only its shape, and writes the
# next look's queue whatever the decision, to exercise it.
summarizeMain <- function(lookNumber = 1L) {
  checkUncodedPaths()
  quick <- nzchar(Sys.getenv("NBS_QUICK"))
  if (!scoringOn && !quick) {
    stop("the scoring is off until the pre-registration commit")
  }
  queued <- readQueues(lookNumber)
  scored <- scoreAll(lookNumber)
  scores <- scored$scores
  fractions <- cummax(pmin(
    informationFractions(lookNumber, queued, scores),
    1
  ))
  decide <- function(final) {
    t <- fractions
    if (final) {
      t[lookNumber] <- 1
    }
    boundaries <- ldBoundaries(t)
    look <- lookStatistics(scores, boundaries[lookNumber])
    list(
      t = t,
      boundaries = boundaries,
      look = look,
      decision = enumerateDecisions(look, final)
    )
  }
  # A look is final at the last look, when the stop is spent, or when it
  # can fund no block.
  final <- lookNumber >= lastLook ||
    sum(scored$costs$seconds) / speedFactor() / 3600 >= mainStop
  state <- decide(final)
  allocation <- NULL
  if (
    (!state$decision$resolved && !final) || (quick && lookNumber < lastLook)
  ) {
    allocation <- allocateBlock(state$look, scored$costs, queued, quick)
    if (length(allocation$funded)) {
      runs <- queued[!grepl("-fold", queued$cell), ]
      cells <- unlist(allocation$funded, use.names = FALSE)
      last <- tapply(runs$rep, runs$cell, max)[cells]
      writeQueue(
        lookNumber + 1L,
        cells,
        lapply(last, function(r) r + seq_len(blockReps))
      )
    } else if (!final) {
      final <- TRUE
      state <- decide(TRUE)
    }
  }
  t <- state$t
  boundaries <- state$boundaries
  b <- boundaries[lookNumber]
  look <- state$look
  decision <- state$decision
  guard <- NULL
  if (decision$resolved || final || quick) {
    order <- decision$order
    if (is.null(order)) {
      order <- enumerateDecisions(look, TRUE)$order
    }
    if (quick && length(order) == 0L) {
      order <- fittedNames
    }
    guard <- stepFive(scores, order, queued)
    if (!is.na(guard$waiting)) {
      runGuardPlan(fileCode(guard$waiting))
    }
  }
  saveRDS(
    list(
      look = lookNumber,
      fractions = t,
      boundaries = boundaries,
      final = final,
      statistics = look,
      decision = decision,
      allocation = allocation,
      guard = guard,
      costs = scored$costs
    ),
    outPath(sprintf("main-decision-%d.rds", lookNumber))
  )
  if (!scoringOn) {
    v <- look$verdicts
    cat(sprintf(
      "quick: %d scored rows, %d cells, classes %s; %d NA coverages; %d boundaries; %d candidates with %d sequential statistics, %d width cells, %d pairs; %d outcomes enumerated; next queue %s; guard %s\n",
      nrow(scores),
      length(unique(scores$cell)),
      paste(sort(unique(scores$class)), collapse = " "),
      sum(is.na(scores$cover90)),
      length(boundaries),
      length(v),
      sum(vapply(v, function(x) length(x$sequential), 0L)),
      sum(vapply(v, function(x) length(x$widthStatus), 0L)),
      length(look$pairs),
      length(decision$outcomes),
      if (is.null(allocation)) {
        "none"
      } else {
        sprintf("%d bands", length(allocation$funded))
      },
      if (is.null(guard)) {
        "not reached"
      } else if (!is.na(guard$waiting)) {
        "waiting"
      } else {
        "read"
      }
    ))
    return(invisible(NULL))
  }
  printDecision(
    lookNumber,
    t,
    b,
    final,
    look,
    decision,
    allocation,
    guard,
    scores
  )
}

# The decision's printout, after the pre-registration.
printDecision <- function(
  lookNumber,
  t,
  b,
  final,
  look,
  decision,
  allocation,
  guard,
  scores
) {
  cat(sprintf(
    "== look %d%s: information fractions %s, boundary %.2f ==\n",
    lookNumber,
    if (final) " (final)" else "",
    paste(sprintf("%.3f", t), collapse = ", "),
    b
  ))
  for (candidate in names(look$verdicts)) {
    v <- look$verdicts[[candidate]]
    cat(sprintf("%s: step 2 %s, rank %+.4f\n", candidate, v$stepTwo, v$rank))
    for (name in names(v$sequential)) {
      x <- v$statistics[[name]]
      bar <- if (startsWith(name, "interval")) intervalScoreBar else gapBar
      cat(sprintf(
        "  %s: %+.4f (se %.4f), z for %.2f %.2f: %s\n",
        name,
        x$estimate,
        x$se,
        bar,
        v$sequential[[name]]$z,
        v$sequential[[name]]$status
      ))
    }
    if (!is.null(v$logScore)) {
      cat(sprintf(
        "  log score, C0's less the candidate's %+.4f (95 percent upper %+.4f): %s\n",
        v$logScore$estimate,
        v$logScore$estimate + qnorm(0.975) * v$logScore$se,
        v$logStatus
      ))
    }
    cat(sprintf(
      "  cells worse by the bar: %s; shape below 1: %s\n",
      if (length(v$cellsWorse)) {
        paste(v$cellsWorse, collapse = ", ")
      } else {
        "none"
      },
      if (length(v$shapeWorse)) paste(v$shapeWorse, collapse = ", ") else "none"
    ))
    if (!is.null(v$width)) {
      q <- qt(0.975, pmax(v$width$m - 1L, 1L))
      cat(sprintf(
        "  width guard, F interval score ratio (t upper): %s\n",
        paste(
          sprintf(
            "%s %.3f (%.3f) %s",
            v$width$cell,
            v$width$value,
            v$width$value + q * v$width$boundSe,
            v$widthStatus
          ),
          collapse = "; "
        )
      ))
    }
  }
  for (pair in look$pairs) {
    cat(sprintf(
      "pair %s ahead of %s: lead %+.4f (se %.4f), z for %.3f %.2f: %s\n",
      pair$better,
      pair$worse,
      pair$lead$estimate,
      pair$lead$se,
      tieMargin,
      pair$z,
      pair$status
    ))
  }
  cat(sprintf(
    "recommendations under every resolution: %s; %s\n",
    paste(decision$outcomes, collapse = ", "),
    if (decision$resolved) "resolved" else "unresolved"
  ))
  if (!is.null(allocation)) {
    cat(sprintf(
      "block allocation: rate by band %s; funded %s; core-hours spent %.1f, after the block %.1f\n",
      paste(
        sprintf("%s %.3g", names(allocation$rates), allocation$rates),
        collapse = ", "
      ),
      if (length(allocation$funded)) {
        paste(
          sprintf(
            "%s (%d cells)",
            names(allocation$funded),
            lengths(allocation$funded)
          ),
          collapse = ", "
        )
      } else {
        "none"
      },
      allocation$spent,
      allocation$projected
    ))
  }
  if (!is.null(guard)) {
    if (!is.na(guard$waiting)) {
      cat(sprintf("step 5 waits for %s's default readings\n", guard$waiting))
    } else {
      cat(sprintf("step 5: recommended %s\n", guard$chosen))
    }
  }
  # Reported beside, not read: rows below 0.01, 50, 80 and 95 percent
  # coverage, RMSE.
  beside <- scores[scores$kind == "long", ]
  print(
    aggregate(
      cbind(cover50, cover80, cover90, cover95, rmse) ~ candidate +
        class +
        band,
      data = beside,
      FUN = mean
    ),
    digits = 3
  )
}

## ------------------------------------------------------- box calibration

# The box's speed against the laptop's, before the main run's queue starts:
# C0's long run on three cells at a fixed length (500 sweeps of burn-in and
# one chunk of 2,500, 4 chains) on the kprior build with its prior unset,
# replication 999, which no main or smoke replication uses. The laptop's
# seconds for the same runs are recorded below; their median ratio rescales
# the sizing table and converts the box's measured seconds for the stop.
calibrationCells <- c("roaches-w3", "crimtab-w3", "mite")
calibrationLaptop <- c(`roaches-w3` = 95.2, `crimtab-w3` = 91.5, mite = 77.3)
runCalibrate <- function(cell) {
  requireBuild("kprior")
  if (buildLabel != "kprior") {
    stop("the calibration runs with NBS_KPRIOR unset and the move on")
  }
  data <- mainCell(cell, 999L)
  longBurn <<- 500L
  longTarget <<- Inf
  longCap <<- longChunk
  run <- mainRun(data, NULL, 999L, list(), FALSE)
  saveRDS(
    list(cell = cell, seconds = run$seconds, sweeps = run$sweeps),
    outPath("checks", sprintf("calibrate-%s.rds", cell))
  )
  cat(sprintf("calibrate %s: %.1f s\n", cell, run$seconds))
}

# The sizing table in the laptop's core-hours (the plan's Box readiness),
# printed rescaled by the box's speed, with the wall time on 4 workers.
sizingLaptop <- c(
  `C0's long runs` = 18.0,
  `own runs` = 23.3,
  readings = 2.6,
  `C0's default readings` = 2.2,
  `the guard's default readings` = 2.2
)
summarizeCalibration <- function() {
  seconds <- vapply(
    calibrationCells,
    function(cell) {
      path <- file.path(runDir, "checks", sprintf("calibrate-%s.rds", cell))
      if (file.exists(path)) readRDS(path)$seconds else NA_real_
    },
    0
  )
  ratio <- seconds / calibrationLaptop[calibrationCells]
  print(data.frame(
    cell = calibrationCells,
    here = round(seconds, 1),
    laptop = calibrationLaptop[calibrationCells],
    ratio = round(ratio, 3),
    row.names = NULL
  ))
  factor <- median(ratio, na.rm = TRUE)
  saveRDS(list(factor = factor, ratio = ratio), outPath("calibration.rds"))
  cat(sprintf("speed factor (seconds here a laptop second): %.3f\n", factor))
  sized <- c(sizingLaptop, total = sum(sizingLaptop))
  print(data.frame(
    component = names(sized),
    laptop = round(sized, 1),
    here = round(sized * factor, 1),
    wall4 = round(sized * factor / 4, 1),
    row.names = NULL
  ))
  cat(sprintf(
    "the stop: %d core-hours of the sizing table's, %.1f here\n",
    mainStop,
    mainStop * factor
  ))
}

## ------------------------------------------------------------ smoke runs

# NBS_QUICK=1 shortens every run and count for a smoke run of the pipeline;
# its output says nothing and belongs in a run directory of its own. Its
# replications are numbered from 900, so that no seed it draws with is one
# the main run's replications (at most 30) draw with; each cell runs two at
# look 1 (a fold one), so that its statistics have a spread, and one a
# block.
if (nzchar(Sys.getenv("NBS_QUICK"))) {
  repOffset <- 899L
  cellReps <- function(cell, band) if (grepl("-fold", cell)) 1L else 2L
  blockReps <- 1L
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
  longBurn <- 100L
  ownBelowMain <- 10
  defaultBurn <- 50L
  defaultSamples <- 50L
  bootstrapReps <- 50L
}

## ---------------------------------------------------------------- dispatch

args <- commandArgs(trailingOnly = TRUE)
action <- if (length(args) > 0L) args[1L] else ""
switch(
  action,
  export = exportData(),
  drift = runDrift(as.integer(args[2L]), as.integer(args[3L])),
  identity = runIdentity(as.integer(args[2L])),
  prior = runPrior(args[2L]),
  sbc = runSbc(as.integer(args[2L]), as.integer(args[3L])),
  compare = runCompare(as.integer(args[2L]), as.integer(args[3L])),
  reweight = runReweight(args[2L]),
  stage2 = runStage2(args[2L], as.integer(args[3L])),
  correction = runCorrection(args[2L], as.integer(args[3L])),
  plan = runPlan(),
  guardplan = runGuardPlan(args[2L]),
  calibrate = runCalibrate(args[2L]),
  uncoded = runUncodedCheck(),
  long = runLong(args[2L], as.integer(args[3L])),
  default = runDefault(
    args[2L],
    as.integer(args[3L]),
    if (length(args) > 3L) {
      candidates <- mainCandidates(1)
      names(candidates)[match(args[4L], fileCode(names(candidates)))]
    } else {
      "C0"
    }
  ),
  `prior-checks` = runPriorChecks(),
  truth = switch(
    args[2L],
    bart = runTruthBart(args[3L]),
    cf = runTruthCounterfactual(args[3L]),
    stop("truth takes bart or cf; the additive truth is written by export")
  ),
  summarize = switch(
    args[2L],
    drift = summarizeDrift(),
    identity = summarizeIdentity(),
    sbc = summarizeSbc(),
    compare = summarizeCompare(),
    reweight = summarizeReweight(),
    prior = summarizePrior(),
    stage2 = summarizeStage2(),
    correction = summarizeCorrection(),
    calibration = summarizeCalibration(),
    main = summarizeMain(if (length(args) > 2L) as.integer(args[3L]) else 1L),
    stop("unknown stage to summarize")
  ),
  stop("unknown action; see the usage at the top of this file")
)
