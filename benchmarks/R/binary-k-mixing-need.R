#!/usr/bin/env Rscript

# The k-mixing study: does slow mixing of the leaf-prior scale k harm the
# estimands binary-fit users report, and where? Truths drawn from the model's
# own prior on real covariates (simulation-based calibration, SBC) are read
# against the shipped sampler, an oracle with k fixed at its true value and,
# on reference replications, the collapsed move and k fixed at its posterior's
# quantiles; long references that both kernels reach score each kernel's
# readings against the other kernel's half. The design, the measures, the
# bars and the reading of a cell are fixed in docs/plans/binary-k-mixing-need.md;
# this file implements them and decides nothing the plan leaves open except
# where its comments say so.
#
# BUILDS. One R process sees one build: whatever dbarts R_LIBS finds first,
# named by KMN_BUILD ("shipped", "move" or "main"). An arm refuses to run under
# a build that is not its own. A fixed-k fit runs neither the probit step nor
# the move, so the oracle and the fixed-k arms run under the shipped build.
#
# DATA. DBARTS_BIS_DATA is the binary interval study's data root (twins/,
# ihdp/, causal/, uci/); ACIC 2016 comes from the aciccomp2016 package.
# KMN_BIS_RUN, the binary interval study's run directory, is read by `screen`
# (its pilot references) and `timing summarize` (its time ratios), and by
# `check pilot-stored`.
#
# RUN DIRECTORY. KMN_RUN holds every output: design/, prop/, fits/, scores/,
# ref/, fixedk/, refscore/, timing/, checks/, bands.rds and the readings.
#
# Usage (each action one process):
#   Rscript binary-k-mixing-need.R cells
#   Rscript binary-k-mixing-need.R jobs <phase> <file>     phases below
#   Rscript binary-k-mixing-need.R design <regime> <link> <rep>
#   Rscript binary-k-mixing-need.R fit <regime> <link> <rep> <arm>[,<arm>]
#   Rscript binary-k-mixing-need.R score <regime> <link> <rep>
#   Rscript binary-k-mixing-need.R reference run <regime> <link> <rep> <part>
#   Rscript binary-k-mixing-need.R reference accept <regime> <link> <rep> [final]
#   Rscript binary-k-mixing-need.R reference escalate <regime> <link> <rep> <part>
#   Rscript binary-k-mixing-need.R fixedk <regime> <link> <rep>
#   Rscript binary-k-mixing-need.R refscore <regime> <link> <rep>
#   Rscript binary-k-mixing-need.R timing <regime> <link> <rep> <t>
#   Rscript binary-k-mixing-need.R timing summarize
#   Rscript binary-k-mixing-need.R screen <outfile>
#   Rscript binary-k-mixing-need.R check <what> [args]   (see runChecks)
#   Rscript binary-k-mixing-need.R read [tranche]
# KMN_QUICK=1 runs every length at a tenth (the smoke); quick output is marked
# and readers refuse to mix it with full output.

suppressPackageStartupMessages(library(dbarts))

## ---------------------------------------------------------------- settings

quick <- nzchar(Sys.getenv("KMN_QUICK"))
tenth <- function(n) {
  if (quick) max(2L, as.integer(round(n / 10))) else as.integer(n)
}

# MCMC lengths, in sweeps. A fit runs past the default and is read at several
# lengths: a longer run's first draws are a shorter run's.
mcmc <- list(
  burn = tenth(500), # every study fit
  defaultDraws = tenth(500), # the default reading
  draws4 = tenth(2000), # a 4-chain fit's draws a chain (the long reading)
  draws10 = tenth(1000), # a 10-chain fit's
  eqUnit = tenth(1000), # burn + default draws: equal time is eqUnit * r
  longFactor = 10L, # the long oracle: 10 times every length
  propBurn = tenth(500),
  propDraws = tenth(500),
  refBurn = tenth(5000),
  refRecord = tenth(20000), # recorded sweeps, kept every refThin
  refThin = 20L,
  refTail = tenth(500), # unthinned sweeps after the record: the tails
  escalateRecord = tenth(100000),
  s4bWarmup = tenth(1000), # stan4bart's own default warmup
  s4bDefault = tenth(1000), # its default kept iterations
  s4bLong = tenth(2000),
  s4bRefWarmup = tenth(5000),
  s4bRefIter = tenth(20500),
  timingWarm = tenth(500),
  timingRun = tenth(500)
)
nTrees <- 75L
kDf <- 1.5
kScale <- 2
refChains <- function(chains) if (chains == 10L) 10L else 8L # a kernel's half
refPartChains <- 2L # chains a reference part job runs
pooledRuns <- 200L # pooled stuck runs a reference replication
stuckRowCap <- 500L # quantities a row class keeps in the pooled stuck runs
bootstrapResamples <- 4000L
bar <- 0.02
zProblem <- qnorm(1 - 0.000156) # 3.60
zNoProblem <- qnorm(0.975)
heldOutRows <- 100L
deterministicScore <- 0.995
fixedQuantiles <- c(
  q05 = 0.05,
  q10 = 0.10,
  q30 = 0.30,
  q50 = 0.50,
  q70 = 0.70,
  q90 = 0.90,
  q95 = 0.95
)
ruledQuantiles <- c("q10", "q50", "q90")
# the share of k's posterior each quantile stands for in a pooled stuck run
quantileBands <- c(
  q05 = 0.075,
  q10 = 0.125,
  q30 = 0.2,
  q50 = 0.2,
  q70 = 0.2,
  q90 = 0.125,
  q95 = 0.075
)

## ------------------------------------------------------------------ paths

harnessDir <- function() {
  # a script that sources this file names its directory here
  env <- Sys.getenv("KMN_HARNESS_DIR")
  if (nzchar(env)) {
    return(env)
  }
  flag <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(flag) == 0L) {
    "benchmarks/R"
  } else {
    dirname(normalizePath(sub("^--file=", "", flag[1L])))
  }
}

source(file.path(harnessDir(), "binary-study-common.R"))
source(file.path(harnessDir(), "binary-datasets.R"))

runDir <- function(...) {
  root <- Sys.getenv("KMN_RUN")
  if (!nzchar(root)) {
    stop("KMN_RUN must name the run directory")
  }
  file.path(root, ...)
}

build <- function() {
  b <- Sys.getenv("KMN_BUILD")
  if (!b %in% c("shipped", "move", "main")) {
    stop("KMN_BUILD must be shipped, move or main")
  }
  b
}

# The build each arm runs under; a fixed k runs neither move, so its arms run
# under the shipped build.
armBuild <- c(
  shipped = "shipped",
  move = "move",
  oracle = "shipped",
  longOracle = "shipped",
  fixedk = "shipped",
  mixed = "shipped" # D2: the move's propensity, the shipped outcome stage
)

requireBuild <- function(which) {
  if (build() != which) {
    stop("this belongs to the ", which, " build, not ", build())
  }
}

repTag <- function(rep) sprintf("r%03d", as.integer(rep))

loadChecked <- function(path) {
  if (!file.exists(path)) {
    stop("missing input ", path)
  }
  x <- readRDS(path)
  if (!identical(x$quick, quick)) {
    stop("quick and full output mixed: ", path)
  }
  x
}

## ---------------------------------------------------------------- regimes

# A regime: its kind ("sbc" truth drawn from the prior; "designed" truth for
# D1 and D2; "real" outcomes for B), task (how fits are made and scored),
# links, chains, replications and references by tranche, tier (1 bulk of use,
# 2 common in health data, 3 niche) and selection.
regime <- function(
  id,
  kind,
  task,
  links,
  chains,
  reps,
  refs,
  tier,
  select = "none",
  n = NA,
  datasets = NULL,
  big = FALSE,
  ...
) {
  list(
    id = id,
    kind = kind,
    task = task,
    links = links,
    chains = chains,
    reps = reps,
    refs = refs,
    tier = tier,
    select = select,
    n = n,
    datasets = datasets,
    big = big,
    par = list(...)
  )
}

both <- c("probit", "logistic")
regimes <- list(
  S1 = regime(
    "S1",
    "sbc",
    "prob",
    both,
    4L,
    c(100L, 400L),
    c(12L, 24L),
    3L,
    "nearsep",
    500L,
    c("tictactoe", "banknote", "wdbc", "ionosphere")
  ),
  S2 = regime(
    "S2",
    "sbc",
    "prob",
    both,
    4L,
    c(100L, 400L),
    c(12L, 24L),
    2L,
    "rare",
    500L,
    c("pima", "german", "creditapproval"),
    ratio = TRUE
  ),
  S3 = regime(
    "S3",
    "sbc",
    "prob",
    both,
    4L,
    c(40L, 80L),
    c(6L, 8L),
    2L,
    "rare",
    10000L,
    "adult",
    big = TRUE,
    ratio = TRUE
  ),
  S4 = regime(
    "S4",
    "sbc",
    "prob",
    both,
    4L,
    c(100L, 400L),
    c(12L, 24L),
    2L,
    "strong",
    200L,
    c("kyphosis", "birthwt", "haberman", "cleveland", "infert")
  ),
  S5 = regime(
    "S5",
    "sbc",
    "prob",
    both,
    4L,
    c(100L, 200L),
    c(12L, 20L),
    2L,
    "none",
    1000L,
    "acic2016hd"
  ),
  S6 = regime(
    "S6",
    "sbc",
    "prob",
    both,
    4L,
    c(100L, 400L),
    c(12L, 24L),
    1L,
    "moderate",
    500L,
    c("pima", "german")
  ),
  S7 = regime(
    "S7",
    "sbc",
    "effect",
    both,
    10L,
    c(100L, 300L),
    c(12L, 24L),
    1L,
    "none",
    747L,
    "ihdp"
  ),
  S8 = regime(
    "S8",
    "sbc",
    "effect",
    both,
    10L,
    c(100L, 200L),
    c(12L, 20L),
    2L,
    "rare",
    1000L,
    "acic2016"
  ),
  S9 = regime(
    "S9",
    "sbc",
    "s9",
    "probit",
    10L,
    c(100L, 300L),
    c(12L, 24L),
    1L,
    "overlap",
    747L,
    "ihdp"
  ),
  E1 = regime(
    "E1",
    "sbc",
    "s4b",
    "probit",
    4L,
    c(100L, 200L),
    c(12L, 20L),
    1L,
    "e1",
    1000L,
    "twins",
    group = "brstate"
  ),
  E2 = regime(
    "E2",
    "sbc",
    "s4b",
    "probit",
    4L,
    c(100L, 200L),
    c(12L, 20L),
    1L,
    "e2",
    1000L,
    "twins",
    group = "brstate_reg"
  ),
  E3 = regime(
    "E3",
    "sbc",
    "loop",
    both,
    4L,
    c(100L, 200L),
    c(12L, 20L),
    1L,
    "moderate",
    500L,
    c("pima", "german"),
    masked = 0.3
  ),
  D1 = regime(
    "D1",
    "designed",
    "multi",
    both,
    10L,
    c(4L, 4L),
    c(4L, 4L),
    1L,
    cells = c("arr1-mod", "arr1-sep", "arr2-mod", "null-mod"),
    big = TRUE
  ),
  D2 = regime(
    "D2",
    "designed",
    "d2",
    "probit",
    10L,
    c(6L, 6L),
    c(6L, 6L),
    1L,
    cells = c("r30", "r05"),
    n = 1000L
  ),
  B = regime("B", "real", "real", both, 4L, c(6L, 6L), c(1L, 1L), NA_integer_)
)

# D1's cells are the binary interval study's A6 cells; D2's are ACIC 2016
# setting 1 at base rates 0.3 and 0.05; B's problems are that study's B1 and B2.
regimeCells <- function(reg) {
  if (reg$id == "D1") {
    return(paste0("D1:", reg$par$cells))
  }
  if (reg$id == "D2") {
    return(paste0("D2:", reg$par$cells))
  }
  if (reg$id == "B") {
    return(c(
      paste0("B1:", names(realDatasets)),
      paste0("B2:", names(causalDatasets))
    ))
  }
  reg$id
}

cellRegime <- function(cell) {
  id <- sub(":.*$", "", cell)
  id <- if (id %in% c("B1", "B2")) "B" else id
  reg <- regimes[[id]]
  if (is.null(reg)) {
    stop("no regime for cell ", cell)
  }
  reg
}

cellFile <- function(cell) gsub("[:]", "_", cell)

cellChains <- function(cell) {
  reg <- cellRegime(cell)
  if (grepl("^B2:", cell)) 10L else reg$chains
}

cellTask <- function(cell) {
  reg <- cellRegime(cell)
  if (grepl("^B1:", cell)) {
    "prob"
  } else if (grepl("^B2:", cell)) {
    "effect"
  } else {
    reg$task
  }
}

drawsPerChain <- function(cell) {
  if (cellTask(cell) == "s4b") {
    return(mcmc$s4bLong)
  }
  if (cellChains(cell) == 10L) mcmc$draws10 else mcmc$draws4
}

repsFor <- function(cell, tranche = 1L) cellRegime(cell)$reps[tranche]
refsFor <- function(cell, tranche = 1L) cellRegime(cell)$refs[tranche]
isRefRep <- function(cell, rep, tranche = 1L) {
  as.integer(rep) <= refsFor(cell, tranche)
}

## ----------------------------------------------------------- selection bands

# The separation score: the 5-fold cross-validated AUC of a classification tree
# (rpart) of the outcome on the covariates, a function of the data alone; the
# folds are seeded from the data's own key, never from a truth.
separationScore <- function(x, y, key) {
  n <- length(y)
  folds <- withSeed(seedOf("folds", key), sample(rep_len(1:5, n)))
  pred <- numeric(n)
  d <- data.frame(x, .y = factor(y, levels = 0:1))
  for (f in 1:5) {
    test <- folds == f
    train <- d[!test, , drop = FALSE]
    if (length(unique(train$.y)) < 2L) {
      pred[test] <- mean(as.integer(as.character(train$.y)))
      next
    }
    fit <- rpart::rpart(
      .y ~ .,
      data = train,
      method = "class",
      control = rpart::rpart.control(xval = 0L)
    )
    pred[test] <- predict(fit, d[test, , drop = FALSE], type = "prob")[, 2L]
  }
  auc(pred, y)
}

auc <- function(score, y) {
  pos <- y == 1L
  n1 <- sum(pos)
  n0 <- sum(!pos)
  if (n1 == 0L || n0 == 0L) {
    return(NA_real_)
  }
  r <- rank(score)
  a <- (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  max(a, 1 - a)
}

# The bands, set in stage 1 (`check bands`) and fixed here by the
# pre-registration commit; until then they are read from the run directory.
fixedBands <- NULL

bands <- function() {
  if (!is.null(fixedBands)) {
    return(fixedBands)
  }
  path <- runDir("bands.rds")
  if (!file.exists(path)) {
    stop("no selection bands; run `check bands`")
  }
  readRDS(path)$bands
}

# Whether (x, y) passes a regime's selection; reads the data alone. The rate
# test runs first, so the tree score is computed only on rate passes.
selects <- function(reg, x, y, key) {
  rate <- mean(y)
  out <- list(ok = FALSE, rate = rate, score = NA_real_)
  if (length(unique(y)) < 2L) {
    return(out)
  }
  rateOk <- switch(
    reg$select,
    none = ,
    overlap = ,
    strong = ,
    nearsep = TRUE,
    e1 = rate >= 0.1 && rate <= 0.9,
    rare = (rate >= 0.01 && rate <= 0.05) || (rate >= 0.95 && rate <= 0.99),
    e2 = (rate >= 0.02 && rate <= 0.08) || (rate >= 0.92 && rate <= 0.98),
    moderate = rate >= 0.2 && rate <= 0.8,
    stop("no selection ", reg$select)
  )
  if (!rateOk) {
    return(out)
  }
  if (reg$select %in% c("none", "overlap", "e1")) {
    out$ok <- TRUE
    return(out)
  }
  b <- bands()
  score <- separationScore(x, y, key)
  out$score <- score
  out$ok <- switch(
    reg$select,
    rare = ,
    e2 = score <= b$moderate[2L],
    moderate = score >= b$moderate[1L] && score <= b$moderate[2L],
    strong = score > b$moderate[2L] && score < b$nearsep[1L],
    nearsep = score >= b$nearsep[1L] &&
      score <= b$nearsep[2L] &&
      score < deterministicScore
  )
  out
}

# The k study's simulated cells (benchmarks/R/binary-hyperprior.R, its
# generator restated here): six processes, three sizes, three predictor
# counts, probit truth. The cells whose truth-based sign agreement (y equal to
# 1{f > 0}) is within 0.03 of 0.76 and of 0.90 set the moderate and
# near-separated bands, each the central 80 percent of those cells' scores.
# Only the cells at base rate 0.5 enter: the two signals were defined on
# symmetric data, and at a base rate of 0.05 or 0.2 the agreement measures the
# rate (predicting no event agrees on 95 percent of rows), which puts weak
# rare cells among the near-separated and leaves no strong band.
kStudySurfaces <- list(
  friedman = function(x) {
    10 *
      sin(pi * x[, 1L] * x[, 2L]) +
      20 * (x[, 3L] - 0.5)^2 +
      10 * x[, 4L] +
      5 * x[, 5L]
  },
  linear = function(x) 2 * x[, 1L] + x[, 2L] - x[, 3L],
  strong = function(x) {
    sin(2 * pi * x[, 1L]) + 2 * (x[, 2L] - 0.5)^2 + 0.5 * x[, 3L]
  },
  weak = function(x) {
    sin(2 * pi * x[, 1L]) + 2 * (x[, 2L] - 0.5)^2 + 0.5 * x[, 3L]
  },
  separable = function(x) x[, 1L] + x[, 2L],
  interaction = function(x) {
    3 *
      (x[, 1L] > 0.5) *
      (x[, 2L] > 0.5) -
      3 * (x[, 3L] > 0.5) * (x[, 4L] > 0.5) +
      2 * (2 * x[, 1L] - 1) * (2 * x[, 5L] - 1)
  }
)
kStudySignal <- c(
  friedman = 1.5,
  linear = 1.5,
  strong = 3,
  weak = 0.35,
  separable = 6,
  interaction = 1.5
)

kStudyCell <- function(dgp, n, p, rate, seed) {
  cal <- withSeed(90000L + p, {
    pool <- matrix(runif(20000L * p), 20000L, p)
    raw <- kStudySurfaces[[dgp]](pool)
    center <- mean(raw)
    spread <- sd(raw)
    mu <- kStudySignal[[dgp]] * (raw - center) / spread
    list(
      center = center,
      spread = spread,
      intercept = uniroot(
        function(a) mean(pnorm(mu + a)) - rate,
        c(-30, 30),
        tol = 1e-8
      )$root
    )
  })
  withSeed(seed, {
    x <- matrix(runif(n * p), n, p)
    colnames(x) <- paste0("x", seq_len(p))
    f <- kStudySignal[[dgp]] *
      (kStudySurfaces[[dgp]](x) - cal$center) /
      cal$spread +
      cal$intercept
    y <- rbinom(n, 1L, pnorm(f))
    list(x = as.data.frame(x), y = y, agreement = mean(y == as.integer(f > 0)))
  })
}

computeBands <- function(reps = 5L) {
  rows <- list()
  for (dgp in names(kStudySurfaces)) {
    for (n in c(100L, 500L, 2000L)) {
      for (p in c(5L, 20L, 50L)) {
        for (rate in 0.5) {
          for (r in seq_len(reps)) {
            key <- paste("band", dgp, n, p, rate, r)
            d <- kStudyCell(dgp, n, p, rate, seedOf(key))
            if (length(unique(d$y)) < 2L) {
              next
            }
            rows[[length(rows) + 1L]] <- data.frame(
              dgp = dgp,
              n = n,
              p = p,
              rate = rate,
              rep = r,
              agreement = d$agreement,
              score = separationScore(d$x, d$y, key)
            )
          }
        }
      }
    }
  }
  tab <- do.call(rbind, rows)
  cellMean <- aggregate(agreement ~ dgp + n + p + rate, tab, mean)
  pick <- function(target) {
    keys <- cellMean[
      abs(cellMean$agreement - target) <= 0.03,
      c("dgp", "n", "p", "rate")
    ]
    m <- merge(tab, keys)
    list(
      cells = nrow(keys),
      scores = m$score,
      band = unname(quantile(m$score, c(0.1, 0.9)))
    )
  }
  moderate <- pick(0.76)
  nearsep <- pick(0.90)
  out <- list(
    bands = list(
      moderate = moderate$band,
      nearsep = c(nearsep$band[1L], min(nearsep$band[2L], deterministicScore))
    ),
    cells = c(moderate = moderate$cells, nearsep = nearsep$cells),
    table = tab,
    quick = quick
  )
  if (out$bands$moderate[2L] >= out$bands$nearsep[1L]) {
    stop("the moderate and near-separated bands overlap: no strong band")
  }
  out
}

## ---------------------------------------------------------------- designs

# The covariates of a replication: training rows x and up to heldOutRows
# held-out rows xh (the dataset's rows outside the training sample; a dataset
# all of whose rows train has none), rotating over the regime's datasets with a
# fresh subsample each replication.
realCovariates <- function(dataset, n, rep, tag) {
  raw <- realDatasets[[dataset]]()
  N <- nrow(raw$x)
  nTrain <- min(n, N)
  nHeld <- min(heldOutRows, N - nTrain)
  rows <- withSeed(seedOf(tag, dataset, rep), sample.int(N, nTrain + nHeld))
  x <- raw$x[rows[seq_len(nTrain)], , drop = FALSE]
  xh <- if (nHeld > 0L) {
    raw$x[rows[nTrain + seq_len(nHeld)], , drop = FALSE]
  } else {
    NULL
  }
  rownames(x) <- NULL
  if (!is.null(xh)) {
    rownames(xh) <- NULL
  }
  list(x = x, xh = xh, dataset = dataset)
}

acic2016Hd <- function(rep, n) {
  x <- acic2016Input()
  x <- cbind(x, permutedCopies(x, 192L, seedOf("S5-hd", rep)))
  rows <- withSeed(seedOf("S5-sub", rep), sample.int(nrow(x), n + heldOutRows))
  list(
    x = x[rows[seq_len(n)], , drop = FALSE],
    xh = x[rows[n + seq_len(heldOutRows)], , drop = FALSE],
    dataset = "acic2016hd"
  )
}

twinsCovariates <- function(reg, rep) {
  d <- twinsData()
  x <- d$x
  # the pairs stan4bart keeps: complete on every covariate used
  keep <- which(stats::complete.cases(x))
  rows <- keep[withSeed(
    seedOf("E-sub", reg$id, rep),
    sample.int(length(keep), reg$n)
  )]
  x <- x[rows, , drop = FALSE]
  group <- factor(x[[reg$par$group]])
  x <- x[,
    setdiff(
      names(x),
      c("brstate", "brstate_reg", "stoccfipb", "stoccfipb_reg")
    ),
    drop = FALSE
  ]
  rownames(x) <- NULL
  list(x = x, xh = NULL, group = group, dataset = "twins")
}

# A dbarts sampler on the design, built by bart() exactly as a fit would be,
# one chain: its prior draw is the model the arms fit.
priorSampler <- function(x, link, seed, nTreesPrior = nTrees) {
  d <- data.frame(x, .y = rep_len(c(0, 1), nrow(x)))
  bart(
    .y ~ .,
    data = d,
    family = link,
    n.trees = nTreesPrior,
    n.chains = 1L,
    n.threads = 1L,
    seed = seed,
    verbose = FALSE,
    samplerOnly = TRUE
  )
}

kDraw <- function() kScale * sqrt(rchisq(1L, kDf))

# A draw of (k0, trees, leaves) from the prior, with f0 at every row set given.
priorDraw <- function(sampler, rowSets) {
  sampler$sampleTreesFromPrior()
  k0 <- kDraw()
  sampler$setLeafPrior(dbartsPriors$normal(k = k0))
  sampler$sampleLeafParametersFromPrior()
  f <- lapply(rowSets, function(rows) as.numeric(sampler$predict(rows)))
  list(k0 = k0, f = f)
}

# The SBC truth for one replication of an SBC regime: covariates, prior draw,
# outcome, selection, retried from the prior draw on rejection. `discriminate`
# names a deliberate mismatch for the discrimination runs ("link": the outcome
# drawn through the other link; "trees50": the truth from a 50-tree prior).
sbcDesign <- function(reg, link, rep, discriminate = "") {
  maxAttempts <- if (reg$select %in% c("rare", "e2")) 200000L else 5000L
  key <- paste(reg$id, link, rep, discriminate)
  base <- switch(
    reg$task,
    prob = ,
    loop = if (reg$datasets[1L] == "acic2016hd") {
      acic2016Hd(rep, reg$n)
    } else {
      realCovariates(
        reg$datasets[(as.integer(rep) - 1L) %% length(reg$datasets) + 1L],
        reg$n,
        rep,
        paste0(reg$id, "-sub")
      )
    },
    effect = effectCovariates(reg, rep),
    s9 = list(x = ihdpData()$x, xh = NULL, dataset = "ihdp"),
    s4b = twinsCovariates(reg, rep)
  )
  truthLink <- if (discriminate == "link") setdiff(both, link) else link
  invLink <- linkFunction(truthLink)
  nTreesPrior <- if (discriminate == "trees50") 50L else nTrees
  stage <- list()
  if (reg$task == "s9") {
    # the treatment from a probit prior propensity, selected for weak overlap by
    # an auxiliary logistic regression on the covariates (data only)
    ps <- priorSampler(base$x, "probit", seedOf("S9-ps-sampler", rep))
    for (attempt in seq_len(maxAttempts)) {
      e <- withSeed(seedOf("S9-ps", rep, attempt), priorDraw(ps, list(base$x)))
      z <- withSeed(
        seedOf("S9-z", rep, attempt),
        rbinom(nrow(base$x), 1L, pnorm(e$f[[1L]]))
      )
      aux <- suppressWarnings(stats::glm(
        z ~ .,
        data = data.frame(base$x, z = z),
        family = binomial()
      ))
      outside <- mean(aux$fitted.values < 0.05 | aux$fitted.values > 0.95)
      if (outside >= 0.10 && length(unique(z)) == 2L) break
    }
    if (outside < 0.10) {
      stop("no S9 treatment passed after ", maxAttempts)
    }
    stage <- list(
      e0 = e$f[[1L]],
      ke0 = e$k0,
      z = z,
      overlapOutside = outside,
      psAttempts = attempt
    )
    base$z <- z
  }
  if (reg$task %in% c("effect", "s9")) {
    trainX <- outcomeFrame(base, base$z)
    rowSets <- list(t1 = outcomeFrame(base, 1L), t0 = outcomeFrame(base, 0L))
  } else {
    trainX <- base$x
    rowSets <- list(train = base$x)
    if (!is.null(base$xh)) {
      rowSets$held <- base$xh
    }
  }
  sampler <- priorSampler(
    trainX,
    link,
    seedOf("prior-sampler", key),
    nTreesPrior
  )
  groups <- if (reg$task == "s4b") base$group else NULL
  # E3's masked rows are fixed before the truth, and selection reads only the
  # observed rows: the fit sees those alone, so the truth stays a posterior
  # draw given what the fit sees
  mask <- if (reg$task == "loop") {
    withSeed(
      seedOf("E3-mask", rep),
      sort(sample.int(nrow(trainX), round(reg$par$masked * nrow(trainX))))
    )
  } else {
    integer()
  }
  seen <- setdiff(seq_len(nrow(trainX)), mask)
  for (attempt in seq_len(maxAttempts)) {
    draw <- withSeed(seedOf("prior", key, attempt), {
      d <- priorDraw(sampler, rowSets)
      if (reg$task == "s4b") {
        # stan4bart's prior: decov with shape and scale 1 makes the one group
        # standard deviation unit exponential; group intercepts normal
        tau0 <- rexp(1L)
        b0 <- rnorm(nlevels(groups), 0, tau0)
        d$tau0 <- tau0
        d$b0 <- b0
        d$f$train <- d$f$train + b0[as.integer(groups)]
      }
      d
    })
    if (reg$task %in% c("effect", "s9")) {
      latent <- ifelse(base$z == 1L, draw$f$t1, draw$f$t0)
    } else {
      latent <- draw$f$train
    }
    y <- withSeed(
      seedOf("y", key, attempt),
      rbinom(length(latent), 1L, invLink(latent))
    )
    sel <- selects(
      reg,
      trainX[seen, , drop = FALSE],
      y[seen],
      paste(key, attempt)
    )
    if (sel$ok) break
  }
  if (!sel$ok) {
    stop("no replication passed selection after ", maxAttempts, " attempts")
  }
  out <- list(
    regime = reg$id,
    cell = reg$id,
    link = link,
    rep = as.integer(rep),
    discriminate = discriminate,
    task = reg$task,
    dataset = base$dataset,
    x = base$x,
    xh = base$xh,
    z = base$z,
    ps = base$ps,
    group = groups,
    y = y,
    truth = c(
      list(k0 = draw$k0, f = draw$f, tau0 = draw$tau0, b0 = draw$b0),
      stage[c("e0", "ke0")]
    ),
    attempts = attempt,
    rate = sel$rate,
    score = sel$score,
    stage = stage[setdiff(names(stage), c("e0", "ke0"))],
    quick = quick
  )
  if (reg$task == "loop") {
    out$mask <- mask
  }
  out
}

# Effect designs: S7 is IHDP with its real treatment, S8 a fresh ACIC 2016
# subsample with setting 31's treatment. The propensity is the shipped build's
# 10-chain probit fit, made once per design before the truth is drawn.
effectCovariates <- function(reg, rep) {
  if (reg$id == "S7") {
    d <- ihdpData()
    base <- list(x = d$x, z = d$z, dataset = "ihdp", propKey = "S7-fixed")
  } else {
    rows <- subsampleRows(4802L, reg$n, seedOf("S8-sub", rep))
    x <- acic2016Input()[rows, , drop = FALSE]
    rownames(x) <- NULL
    base <- list(
      x = x,
      z = acic2016Surfaces(31L, rep)$z[rows],
      dataset = "acic2016",
      propKey = sprintf("S8-r%03d", as.integer(rep))
    )
  }
  base$ps <- propensity(base$propKey, base$x, base$z)
  base
}

outcomeFrame <- function(base, z) {
  z <- rep_len(as.integer(z), nrow(base$x))
  if (is.null(base$ps)) {
    data.frame(base$x, trt = z)
  } else {
    data.frame(base$x, trt = z, pscore = base$ps)
  }
}

# bartCause's default treatment model, the posterior mean of a 10-chain probit
# BART fit, by the shipped build, stored once per key; a missing propensity
# made under another build is an error, never a refit.
propensity <- function(key, x, z, which = "shipped") {
  path <- runDir(
    "prop",
    paste0(key, if (which == "shipped") "" else paste0("-", which), ".rds")
  )
  if (file.exists(path)) {
    return(loadChecked(path)$ps)
  }
  requireBuild(which)
  fit <- bart(
    trt ~ .,
    data = data.frame(x, trt = z),
    family = "probit",
    n.chains = 10L,
    n.threads = 1L,
    n.burn = mcmc$propBurn,
    n.samples = mcmc$propDraws,
    n.trees = nTrees,
    seed = seedOf("propensity", key),
    verbose = FALSE,
    keepTrainingFits = TRUE
  )
  ps <- colMeans(pnorm(fit$yhat.train))
  saveAtomic(
    list(
      ps = ps,
      key = key,
      which = which,
      lib = find.package("dbarts"),
      quick = quick
    ),
    path
  )
  ps
}

# D1, the binary interval study's A6: four arms, the generalized propensity the
# shipped build's 4-chain multinomial posterior mean, one per replication.
d1Design <- function(cell, link, rep) {
  arrangement <- switch(
    sub("^D1:", "", cell),
    "arr1-mod" = c(1L, 1),
    "arr1-sep" = c(1L, 2.5),
    "arr2-mod" = c(2L, 1),
    "null-mod" = c(0L, 1)
  )
  b <- a6Base(rep)
  P <- a6Probabilities(
    b,
    list(par = list(arrangement = arrangement[1L], a = arrangement[2L])),
    link
  )
  y <- withSeed(
    seedOf("y", cell, link, rep),
    as.integer(runif(nrow(P)) < P[cbind(seq_len(nrow(P)), b$arm + 1L)])
  )
  path <- runDir("prop", sprintf("D1-r%03d.rds", as.integer(rep)))
  gps <- if (file.exists(path)) {
    loadChecked(path)$ps
  } else {
    requireBuild("shipped")
    fit <- bart(
      arm ~ .,
      data = data.frame(b$x, arm = factor(b$arm, levels = 0:3)),
      family = "multinomial",
      n.chains = 4L,
      n.threads = 1L,
      n.samples = mcmc$propDraws,
      n.burn = mcmc$propBurn,
      n.trees = nTrees,
      seed = seedOf("D1-gps", rep),
      verbose = FALSE
    )
    ps <- apply(fit$yhat.train, c(2L, 3L), mean)
    colnames(ps) <- paste0("gps", 0:3)
    saveAtomic(list(ps = ps, quick = quick), path)
    ps
  }
  list(
    regime = "D1",
    cell = cell,
    link = link,
    rep = as.integer(rep),
    task = "multi",
    x = data.frame(b$x, gps),
    arm = b$arm,
    y = y,
    designed = list(P = P),
    quick = quick
  )
}

# D2: ACIC 2016 setting 1 (linear, weak overlap) at 1,000 units, the surfaces
# mapped to probabilities as the binary interval study's A3 cells were.
d2Design <- function(cell, link, rep) {
  rate <- if (cell == "D2:r30") 0.3 else 0.05
  s <- acic2016Surfaces(1L, rep)
  rows <- subsampleRows(4802L, 1000L, seedOf("D2-sub", rep))
  x <- acic2016Input()[rows, , drop = FALSE]
  rownames(x) <- NULL
  z <- s$z[rows]
  p <- mapSurfaces(list(s$mu0[rows], s$mu1[rows]), z == 0L, 1, rate, link)
  y <- withSeed(
    seedOf("y", cell, link, rep),
    as.integer(runif(length(z)) < ifelse(z == 1L, p[[2L]], p[[1L]]))
  )
  # bartCause's default: the propensity a covariate of the outcome model, each
  # pipeline its own; this is the shipped pipeline's, which the reference also
  # uses (the move's and the mixed pipeline's are made by their fit jobs)
  propKey <- sprintf("D2-r%03d", as.integer(rep))
  ps <- propensity(propKey, x, z)
  list(
    regime = "D2",
    cell = cell,
    link = link,
    rep = as.integer(rep),
    task = "d2",
    x = x,
    z = z,
    ps = ps,
    propKey = propKey,
    y = y,
    designed = list(p0 = p[[1L]], p1 = p[[2L]]),
    quick = quick
  )
}

# B: the binary interval study's real problems and splits; B2's propensity is
# the shipped build's, one per problem.
bDesign <- function(cell, link) {
  name <- sub("^B[12]:", "", cell)
  if (grepl("^B1:", cell)) {
    raw <- realDatasets[[name]]()
    n <- nrow(raw$x)
    seed <- 10L * (8000L + match(name, names(realDatasets))) + 1L
    split <- withSeed(1L, kStudySplit(n, raw$y, seed, n > subsampleAbove))
    held <- split$test[seq_len(min(heldOutRows, length(split$test)))]
    x <- raw$x[split$train, , drop = FALSE]
    xh <- raw$x[held, , drop = FALSE]
    rownames(x) <- NULL
    rownames(xh) <- NULL
    return(list(
      regime = "B",
      cell = cell,
      link = link,
      rep = 1L,
      task = "prob",
      x = x,
      xh = xh,
      y = raw$y[split$train],
      dataset = name,
      quick = quick
    ))
  }
  d <- causalDatasets[[name]]()
  base <- list(x = d$x, z = d$z)
  base$ps <- propensity(paste0("B2-", name), d$x, d$z)
  list(
    regime = "B",
    cell = cell,
    link = link,
    rep = 1L,
    task = "effect",
    x = d$x,
    z = d$z,
    ps = base$ps,
    y = d$y,
    dataset = name,
    quick = quick
  )
}

designFile <- function(cell, link, rep) {
  runDir("design", link, cellFile(cell), paste0(repTag(rep), ".rds"))
}

# Every design is made once, by the shipped build, before any arm runs.
makeDesign <- function(cell, link, rep, discriminate = "") {
  path <- designFile(cell, link, rep)
  if (nzchar(discriminate)) {
    path <- runDir(
      "design-disc",
      discriminate,
      link,
      cellFile(cell),
      paste0(repTag(rep), ".rds")
    )
  }
  if (file.exists(path)) {
    return(invisible(loadChecked(path)))
  }
  requireBuild("shipped")
  reg <- cellRegime(cell)
  t0 <- proc.time()[["elapsed"]]
  d <- switch(
    reg$kind,
    sbc = sbcDesign(reg, link, rep, discriminate),
    designed = if (reg$id == "D1") {
      d1Design(cell, link, rep)
    } else {
      d2Design(cell, link, rep)
    },
    real = bDesign(cell, link)
  )
  d$seconds <- proc.time()[["elapsed"]] - t0
  saveAtomic(d, path)
  cat(sprintf(
    "design %s %s %d: %d attempts, rate %.3f, score %.3f, %.1f s\n",
    cell,
    link,
    as.integer(rep),
    if (is.null(d$attempts)) 1L else d$attempts,
    if (is.null(d$rate)) mean(d$y) else d$rate,
    if (is.null(d$score)) NA_real_ else d$score,
    d$seconds
  ))
  invisible(d)
}

loadDesign <- function(cell, link, rep, discriminate = "") {
  if (nzchar(discriminate)) {
    return(loadChecked(runDir(
      "design-disc",
      discriminate,
      link,
      cellFile(cell),
      paste0(repTag(rep), ".rds")
    )))
  }
  loadChecked(designFile(cell, link, rep))
}

# The design as the arms see it: no truth. Every label, subgroup and score
# below reads this, never design$truth or design$designed.
blind <- function(design) {
  design$truth <- NULL
  design$designed <- NULL
  design
}

## ---------------------------------------------------------------- subgroups

# Subgroups defined from covariates in advance: for IHDP the birth weight
# quartiles and the mother's education levels; elsewhere the quartiles of the
# numeric covariate with the most distinct values and the levels (at most six)
# of the first factor covariate, if any.
mostDistinct <- function(x) {
  num <- names(x)[vapply(x, is.numeric, logical(1L))]
  num[which.max(vapply(x[num], function(c) length(unique(c)), integer(1L)))]
}

quartiles <- function(v) {
  b <- unique(quantile(
    v,
    c(0, 0.25, 0.5, 0.75, 1),
    na.rm = TRUE,
    names = FALSE
  ))
  as.integer(cut(v, b, include.lowest = TRUE))
}

subgroupsOf <- function(design) {
  x <- design$x
  if (identical(design$dataset, "ihdp")) {
    edu <- ifelse(
      x$mom.lths == 1,
      1L,
      ifelse(x$mom.hs == 1, 2L, ifelse(x$mom.scoll == 1, 3L, 4L))
    )
    g <- c(
      lapply(1:4, function(q) quartiles(x$bw) == q),
      lapply(1:4, function(e) edu == e)
    )
    names(g) <- c(paste0("bw", 1:4), paste0("edu", 1:4))
  } else {
    q <- quartiles(x[[mostDistinct(x)]])
    g <- lapply(sort(unique(q)), function(v) q == v)
    names(g) <- paste0("q", sort(unique(q)))
    fac <- names(x)[vapply(
      x,
      function(c) is.factor(c) || is.character(c),
      logical(1L)
    )]
    if (length(fac)) {
      lv <- table(x[[fac[1L]]])
      lv <- names(sort(lv, decreasing = TRUE))[seq_len(min(6L, length(lv)))]
      gf <- lapply(lv, function(l) x[[fac[1L]]] == l)
      names(gf) <- paste0(fac[1L], "=", lv)
      g <- c(g, gf)
    }
  }
  g <- Filter(function(s) sum(s, na.rm = TRUE) >= 5L, g)
  lapply(g, function(s) s & !is.na(s))
}

## -------------------------------------------------------- scored quantities

# Classes of scored quantities and their tiers (a cell's tier is the rarer of
# its regime's and its estimand's); "latent" classes are intervals on the
# latent scale (a probability's interval maps through the link), "row" classes
# are many quantities, the rest scalar estimands.
classInfo <- list(
  fInt = list(tier = 1L, row = TRUE, latent = TRUE),
  fTail = list(tier = 3L, row = TRUE, latent = TRUE),
  bandF = list(tier = 3L, row = FALSE, band = TRUE),
  ite = list(tier = 2L, row = TRUE),
  bandIte = list(tier = 3L, row = FALSE, band = TRUE),
  ate = list(tier = 1L, row = FALSE),
  att = list(tier = 1L, row = FALSE),
  subgroup = list(tier = 1L, row = FALSE),
  logRR = list(tier = 2L, row = FALSE),
  logOR = list(tier = 2L, row = FALSE),
  wate = list(tier = 1L, row = FALSE),
  ipw = list(tier = 1L, row = FALSE),
  ps = list(tier = 1L, row = TRUE, latent = TRUE),
  iw = list(tier = 1L, row = TRUE),
  essW = list(tier = 1L, row = FALSE),
  tau = list(tier = 1L, row = FALSE),
  ranef = list(tier = 1L, row = TRUE),
  contrast = list(tier = 1L, row = TRUE),
  contrastAvg = list(tier = 1L, row = FALSE),
  bandUnit = list(tier = 1L, row = FALSE, band = TRUE),
  bandAll = list(tier = 1L, row = FALSE, band = TRUE),
  maxIw = list(tier = 1L, row = FALSE)
)

# The raw latent columns a task's fits produce, as one draws-by-columns matrix
# (chains pooled chain-major), and the classes computed from it in double
# precision. `parts` carries what a class needs beside the latent draws: the
# propensity draws (S9), the group intercepts and scale (E1, E2).
classesFrom <- function(L, design, link, parts = NULL, labels = NULL) {
  invLink <- linkFunction(link)
  task <- design$task
  out <- list()
  if (task %in% c("prob", "loop", "s4b")) {
    held <- if (is.null(design$xh)) 0L else nrow(design$xh)
    n <- ncol(L) - held
    if (is.null(labels)) {
      labels <- rep(TRUE, ncol(L))
    }
    out$fInt <- L[, labels, drop = FALSE]
    out$fTail <- L[, !labels, drop = FALSE]
    attr(out$fInt, "cols") <- which(labels)
    attr(out$fTail, "cols") <- which(!labels)
    if (held > 0L) {
      out$bandF <- L[, n + seq_len(held), drop = FALSE]
    }
    if (isTRUE(cellRegime(design$cell)$par$ratio)) {
      q <- quartiles(design$x[[mostDistinct(design$x)]])
      P <- invLink(L[, seq_len(n), drop = FALSE])
      out$logRR <- matrix(
        log(
          rowMeans(P[, q == max(q), drop = FALSE]) /
            rowMeans(P[, q == min(q), drop = FALSE])
        ),
        ncol = 1L
      )
    }
    if (task == "s4b") {
      out$ranef <- parts$ranef
      out$tau <- matrix(parts$tau, ncol = 1L)
    }
    return(out)
  }
  if (task %in% c("effect", "s9", "d2")) {
    n <- ncol(L) %/% 2L
    P1 <- invLink(L[, seq_len(n), drop = FALSE])
    P0 <- invLink(L[, n + seq_len(n), drop = FALSE])
    E <- P1 - P0
    z <- design$z
    out$ite <- E
    out$bandIte <- E
    out$ate <- matrix(rowMeans(E), ncol = 1L)
    out$att <- matrix(rowMeans(E[, z == 1L, drop = FALSE]), ncol = 1L)
    sg <- subgroupsOf(design)
    out$subgroup <- vapply(
      sg,
      function(s) rowMeans(E[, s, drop = FALSE]),
      numeric(nrow(E))
    )
    if (task %in% c("effect", "d2") && design$regime != "B") {
      m1 <- rowMeans(P1)
      m0 <- rowMeans(P0)
      out$logRR <- matrix(log(m1 / m0), ncol = 1L)
      out$logOR <- matrix(qlogis(m1) - qlogis(m0), ncol = 1L)
    }
    if (task == "d2") {
      e <- design$ps
      out$wate <- matrix(drop(E %*% e) / sum(e), ncol = 1L)
      out$maxIw <- matrix(
        max(1 / pmin(pmax(ifelse(z == 1L, e, 1 - e), 0.025), 0.975)),
        nrow(E),
        1L
      )
    }
    if (task == "s9") {
      # the propensity draws pair one to one with the outcome draws
      Le <- parts$psLatent
      e <- pnorm(Le)
      eb <- pmin(pmax(e, 0.025), 0.975)
      y <- design$y
      Z <- matrix(z, nrow(e), ncol(e), byrow = TRUE)
      Y <- matrix(y, nrow(e), ncol(e), byrow = TRUE)
      w1 <- Z / eb
      w0 <- (1 - Z) / (1 - eb)
      ipwAte <- rowSums(w1 * Y) / rowSums(w1) - rowSums(w0 * Y) / rowSums(w0)
      wt <- eb / (1 - eb) * (1 - Z)
      ipwAtt <- rowMeans(Y[, z == 1L, drop = FALSE]) -
        rowSums(wt * Y) / rowSums(wt)
      out$ipw <- cbind(ate = ipwAte, att = ipwAtt)
      out$wate <- matrix(rowSums(e * E) / rowSums(e), ncol = 1L)
      out$ps <- Le
      W <- w1 + w0
      out$iw <- W
      out$essW <- matrix(rowSums(W)^2 / rowSums(W^2), ncol = 1L)
    }
    return(out)
  }
  if (task == "multi") {
    n <- ncol(L) %/% 4L
    probs <- lapply(0:3, function(a) {
      invLink(L[, a * n + seq_len(n), drop = FALSE])
    })
    C <- matrix(0, nrow(L), 6L * n)
    for (j in seq_len(nrow(contrastPairs))) {
      C[, (seq_len(n) - 1L) * 6L + j] <- probs[[contrastPairs[j, 1L] + 1L]] -
        probs[[contrastPairs[j, 2L] + 1L]]
    }
    out$contrast <- C
    out$contrastAvg <- vapply(
      seq_len(6L),
      function(j) rowMeans(C[, (seq_len(n) - 1L) * 6L + j, drop = FALSE]),
      numeric(nrow(C))
    )
    out$bandUnit <- C
    out$bandAll <- C
    return(out)
  }
  stop("no classes for task ", task)
}

# The truth's value of every class, from the same functionals at the truth.
truthClasses <- function(design, labels) {
  t <- design$truth
  if (is.null(t)) {
    return(NULL)
  }
  task <- design$task
  L <- if (task %in% c("effect", "s9")) {
    matrix(c(t$f$t1, t$f$t0), nrow = 1L)
  } else {
    matrix(c(t$f$train, t$f$held), nrow = 1L)
  }
  parts <- NULL
  if (task == "s9") {
    parts <- list(psLatent = matrix(t$e0, nrow = 1L))
  }
  if (task == "s4b") {
    parts <- list(ranef = matrix(t$b0, nrow = 1L), tau = t$tau0)
  }
  out <- classesFrom(L, design, design$link, parts, labels)
  lapply(out, function(m) as.numeric(if (is.matrix(m)) m[1L, ] else m))
}

## --------------------------------------------------------------- summaries

probsKept <- sort(c(probsPointwise, 0.5))

# Per class: quantiles at the four levels and the median, mean and standard
# deviation per quantity, each column with its own handling; a band class
# also keeps the rank-based band at 0.90 and 0.95 (and, for D1's unit band,
# one over each unit's six contrasts).
summarizeClass <- function(Q, name) {
  info <- classInfo[[name]]
  if (ncol(Q) == 0L) {
    return(list(
      q = matrix(0, length(probsKept), 0L),
      mean = numeric(),
      sd = numeric(),
      S = nrow(Q)
    ))
  }
  q <- apply(
    Q,
    2L,
    stats::quantile,
    probs = probsKept,
    names = FALSE,
    type = 7L
  )
  if (!is.matrix(q)) {
    q <- matrix(q, ncol = ncol(Q))
  }
  out <- list(
    q = q,
    mean = colMeans(Q),
    sd = apply(Q, 2L, stats::sd),
    S = nrow(Q)
  )
  if (isTRUE(info$band)) {
    bs <- summarizeDraws(Q, groupSize = if (name == "bandUnit") 6L else NULL)
    out$bands <- if (name == "bandUnit") {
      lapply(bs$unitBands, function(b) list(lo = b$lo, hi = b$hi))
    } else {
      lapply(bs$bands, function(b) list(lo = b$lo, hi = b$hi))
    }
    names(out$bands) <- c("0.9", "0.95")
  }
  out
}

intervalOf <- function(s, level) {
  lo <- match(round((1 - level) / 2, 6), round(probsKept, 6))
  hi <- match(round((1 + level) / 2, 6), round(probsKept, 6))
  list(lo = s$q[lo, ], hi = s$q[hi, ])
}

summarizeAll <- function(classes) {
  out <- lapply(names(classes), function(nm) summarizeClass(classes[[nm]], nm))
  names(out) <- names(classes)
  for (nm in c("fInt", "fTail")) {
    if (!is.null(out[[nm]])) {
      out[[nm]]$cols <- attr(classes[[nm]], "cols")
    }
  }
  out
}

## -------------------------------------------------------------------- arms

# The training frame and test rows of a dbarts fit of the design, in the order
# every reading and the reference use.
fitData <- function(design) {
  task <- design$task
  if (task %in% c("prob", "loop")) {
    return(list(train = data.frame(design$x, y = design$y), test = design$xh))
  }
  if (task %in% c("effect", "d2", "s9")) {
    base <- list(x = design$x, ps = if (task == "s9") NULL else design$ps)
    return(list(
      train = data.frame(outcomeFrame(base, design$z), y = design$y),
      test = rbind(outcomeFrame(base, 1L), outcomeFrame(base, 0L))
    ))
  }
  if (task == "multi") {
    arm <- factor(design$arm, levels = 0:3)
    return(list(
      train = data.frame(design$x, arm = arm, y = design$y),
      test = do.call(
        rbind,
        lapply(0:3, function(a) {
          data.frame(
            design$x,
            arm = factor(rep(a, nrow(design$x)), levels = 0:3)
          )
        })
      )
    ))
  }
  stop("no dbarts fit for task ", task)
}

# A fit's draws as chains x draws x columns, plus k (draws x chains or NULL)
# and seconds. `kFixed` fixes k; `lengths` = c(burn, draws a chain).
fitBart <- function(design, link, chains, lengths, seed, kFixed = NULL) {
  task <- design$task
  fd <- fitData(design)
  train <- fd$train
  test <- fd$test
  args <- list(
    y ~ .,
    data = train,
    family = link,
    n.trees = nTrees,
    n.chains = chains,
    n.threads = 1L,
    n.burn = lengths[1L],
    n.samples = lengths[2L],
    seed = seed,
    verbose = FALSE,
    combineChains = FALSE,
    keepTrainingFits = task %in% c("prob", "loop")
  )
  if (!is.null(test)) {
    args$test <- test
  }
  if (!is.null(kFixed)) {
    args$k <- kFixed
  }
  seconds <- system.time(fit <- do.call(bart, args))[["elapsed"]]
  arr <- if (task %in% c("prob", "loop")) {
    a <- as3(fit$yhat.train, chains)
    if (!is.null(test)) {
      a <- bindColumns(a, as3(fit$yhat.test, chains))
    }
    a
  } else {
    as3(fit$yhat.test, chains)
  }
  if (any(!is.finite(arr))) {
    stop("non-finite draws")
  }
  list(arr = arr, k = kMatrix(fit[["k"]], chains), seconds = seconds)
}

# bart()'s per-chain draws as chains x draws x columns, one chain included.
as3 <- function(a, chains) {
  d <- dim(a)
  if (length(d) == 3L) {
    return(array(as.numeric(a), d))
  }
  array(
    as.numeric(a),
    c(chains, length(a) %/% (chains * d[length(d)]), d[length(d)])
  )
}

bindColumns <- function(a, b) {
  out <- array(0, c(dim(a)[1L], dim(a)[2L], dim(a)[3L] + dim(b)[3L]))
  out[,, seq_len(dim(a)[3L])] <- a
  out[,, dim(a)[3L] + seq_len(dim(b)[3L])] <- b
  out
}

# k as draws x chains.
kMatrix <- function(k, chains) {
  if (is.null(k)) {
    return(NULL)
  }
  if (is.matrix(k)) {
    return(matrix(t(k), ncol = chains))
  }
  matrix(as.numeric(k), ncol = chains)
}

# S9's propensity stage under the same build and lengths, its draws paired one
# to one with the outcome's.
fitPropensityS9 <- function(design, chains, lengths, seed, kFixed = NULL) {
  args <- list(
    trt ~ .,
    data = data.frame(design$x, trt = design$z),
    family = "probit",
    n.trees = nTrees,
    n.chains = chains,
    n.threads = 1L,
    n.burn = lengths[1L],
    n.samples = lengths[2L],
    seed = seed,
    verbose = FALSE,
    combineChains = FALSE
  )
  if (!is.null(kFixed)) {
    args$k <- kFixed
  }
  seconds <- system.time(fit <- do.call(bart, args))[["elapsed"]]
  list(
    arr = as3(fit$yhat.train, chains),
    k = kMatrix(fit[["k"]], chains),
    seconds = seconds
  )
}

# stan4bart with its own defaults (a varying intercept on the design's group),
# warmup fixed so that a longer run's draws extend a shorter one's.
fitStan4bart <- function(
  design,
  chains,
  lengths,
  seed,
  kFixed = NULL,
  callback = NULL
) {
  suppressPackageStartupMessages(library(stan4bart))
  d <- data.frame(design$x, grp = design$group, y = design$y)
  args <- list(
    y ~ bart(. - grp) + (1 | grp),
    data = d,
    verbose = -1L,
    chains = chains,
    cores = 1L,
    iter = lengths[1L] + lengths[2L],
    warmup = lengths[1L],
    seed = seed,
    keep_fits = is.null(callback)
  )
  if (!is.null(kFixed)) {
    args$bart_args <- list(k = kFixed)
  }
  if (!is.null(callback)) {
    args$callback <- callback
  }
  seconds <- system.time(fit <- do.call(stan4bart::stan4bart, args))[[
    "elapsed"
  ]]
  list(fit = fit, seconds = seconds)
}

s4bArrays <- function(fit, chains) {
  bartL <- dbarts::extract(fit, "indiv.bart", combine_chains = FALSE)
  ranL <- dbarts::extract(fit, "indiv.ranef", combine_chains = FALSE)
  L <- bartL + ranL # rows x draws x chains
  arr <- aperm(L, c(3L, 2L, 1L))
  re <- dbarts::extract(fit, "ranef", combine_chains = FALSE)[[1L]] # 1 x groups x draws x chains
  ranef <- aperm(array(re, dim(re)[-1L]), c(3L, 2L, 1L))
  Sig <- dbarts::extract(fit, "Sigma", combine_chains = FALSE)[[1L]]
  tau <- t(matrix(sqrt(Sig), ncol = chains))
  # a fixed k is not a modeled parameter, and stan4bart then has none to extract
  k <- tryCatch(
    dbarts::extract(fit, "k", combine_chains = FALSE),
    error = function(e) NULL
  )
  list(
    arr = arr,
    ranef = ranef,
    tau = tau,
    k = if (is.null(k)) NULL else matrix(k, ncol = chains)
  )
}

# E3: the dbartsSampler loop, one single-chain sampler a chain, each sweep's
# missing outcomes redrawn from the current fit and set as the response.
fitLoop <- function(design, link, chains, lengths, seed, kFixed = NULL) {
  invLink <- linkFunction(link)
  mask <- design$mask
  y <- design$y
  t0 <- proc.time()[["elapsed"]]
  per <- lapply(seq_len(chains), function(ch) {
    yy <- y
    yy[mask] <- withSeed(
      seedOf("E3-start", seed, ch),
      rbinom(length(mask), 1L, mean(y[-mask]))
    )
    args <- list(
      y ~ .,
      data = data.frame(design$x, y = yy),
      test = design$xh,
      family = link,
      n.trees = nTrees,
      n.chains = 1L,
      n.threads = 1L,
      seed = seed + ch,
      verbose = FALSE,
      samplerOnly = TRUE
    )
    if (!is.null(kFixed)) {
      args$k <- kFixed
    }
    s <- do.call(bart, args)
    total <- lengths[1L] + lengths[2L]
    nCols <- nrow(design$x) + if (is.null(design$xh)) 0L else nrow(design$xh)
    draws <- matrix(0, lengths[2L], nCols)
    k <- numeric(lengths[2L])
    set.seed(seed + 1000L * ch)
    for (it in seq_len(total)) {
      r <- s$run(0L, 1L)
      f <- as.numeric(r$train)
      yy[mask] <- rbinom(length(mask), 1L, invLink(f[mask]))
      s$setResponse(yy)
      if (it > lengths[1L]) {
        draws[it - lengths[1L], ] <- c(
          f,
          if (is.null(r$test)) NULL else as.numeric(r$test)
        )
        k[it - lengths[1L]] <- if (is.null(r$k)) {
          NA_real_
        } else {
          as.numeric(r$k)[1L]
        }
      }
    }
    list(draws = draws, k = k)
  })
  arr <- aperm(simplify2array(lapply(per, `[[`, "draws")), c(3L, 1L, 2L))
  k <- if (is.null(kFixed)) {
    vapply(per, `[[`, numeric(lengths[2L]), "k")
  } else {
    NULL
  }
  list(arr = arr, k = k, seconds = proc.time()[["elapsed"]] - t0)
}

# Pooled draws (draws x columns, chain-major) at some draws of every chain.
poolDraws <- function(arr, idx, chains = seq_len(dim(arr)[1L])) {
  a <- arr[chains, idx, , drop = FALSE]
  matrix(aperm(a, c(2L, 1L, 3L)), dim(a)[1L] * dim(a)[2L], dim(a)[3L])
}

poolMatrix <- function(m, idx, chains = seq_len(dim(m)[1L])) {
  if (is.null(m)) {
    return(NULL)
  }
  if (length(dim(m)) == 2L) {
    return(as.vector(t(m[chains, idx, drop = FALSE])))
  }
  poolDraws(m, idx, chains)
}

mixingOf <- function(series) {
  if (is.null(series) || nrow(series) < 4L || any(!is.finite(series))) {
    return(c(rhat = NA_real_, ess = NA_real_))
  }
  c(rhat = posterior::rhat(series), ess = posterior::ess_bulk(series))
}

timingFile <- function() runDir("timing", "ratios.rds")

timeRatio <- function(cell, link) {
  fixed <- Sys.getenv("KMN_FIXED_R")
  if (nzchar(fixed)) {
    return(as.numeric(fixed))
  }
  tab <- readRDS(timingFile())
  r <- tab$r[tab$link == link & tab$cell == cell]
  if (length(r) != 1L || !is.finite(r)) {
    stop("no time ratio for ", cell, " ", link)
  }
  r
}

readingIndices <- function(cell, link, arm) {
  n <- drawsPerChain(cell)
  if (cellTask(cell) == "s4b") {
    default <- seq_len(mcmc$s4bDefault)
    unit <- mcmc$s4bWarmup + mcmc$s4bDefault
    burn <- mcmc$s4bWarmup
  } else {
    default <- seq_len(mcmc$defaultDraws)
    unit <- mcmc$eqUnit
    burn <- mcmc$burn
  }
  out <- list(default = default, long = seq_len(n))
  if (arm %in% c("shipped")) {
    L <- round(unit * timeRatio(cell, link)) - burn
    if (L > n) {
      stop("equal time past the run for ", cell)
    }
    out$eqDraws <- seq_len(L)
  }
  out
}

armSeed <- function(cell, link, rep, arm) {
  # the shipped and move arms (and the oracle) start from one seed, so the
  # move's default is the k-fast control from the same start
  seedOf(
    "fit",
    cell,
    link,
    rep,
    if (arm %in% c("shipped", "move", "oracle", "mixed")) "arm" else arm
  )
}

fitFile <- function(cell, link, rep, arm, discriminate = "") {
  if (nzchar(discriminate)) {
    return(runDir(
      "fits-disc",
      discriminate,
      link,
      cellFile(cell),
      sprintf("%s-%s.rds", repTag(rep), arm)
    ))
  }
  runDir("fits", link, cellFile(cell), sprintf("%s-%s.rds", repTag(rep), arm))
}

# One arm of one replication, under the arm's build, read at every length.
# Stored: the class summaries of each reading and k's mixing; no draws.
runArm <- function(cell, link, rep, arm, discriminate = "") {
  requireBuild(armBuild[[arm]])
  design <- blind(loadDesign(cell, link, rep, discriminate))
  if (cellRegime(cell)$kind == "real" || cellRegime(cell)$kind == "designed") {
    if (arm %in% c("oracle", "longOracle")) {
      stop("no oracle without a drawn truth")
    }
  }
  truthK <- if (arm %in% c("oracle", "longOracle")) {
    t <- loadDesign(cell, link, rep, discriminate)$truth
    list(k0 = t$k0, ke0 = t$ke0)
  } else {
    list()
  }
  chains <- cellChains(cell)
  n <- drawsPerChain(cell)
  task <- cellTask(cell)
  lengths <- if (arm == "longOracle") {
    if (task == "s4b") {
      c(mcmc$s4bWarmup, mcmc$s4bDefault) * mcmc$longFactor
    } else {
      c(mcmc$burn, mcmc$defaultDraws) * mcmc$longFactor
    }
  } else if (task == "s4b") {
    c(mcmc$s4bWarmup, n)
  } else {
    c(mcmc$burn, n)
  }
  seed <- armSeed(cell, link, rep, arm)
  if (task == "d2" && arm %in% c("move", "mixed")) {
    full <- loadDesign(cell, link, rep)
    design$ps <- propensity(full$propKey, full$x, full$z, which = "move")
  }
  parts <- NULL
  if (task == "s4b") {
    f <- fitStan4bart(design, chains, lengths, seed, truthK$k0)
    a <- s4bArrays(f$fit, chains)
    res <- list(arr = a$arr, k = a$k, seconds = f$seconds)
    parts <- list(ranef = a$ranef, tau = a$tau)
  } else if (task == "loop") {
    res <- fitLoop(design, link, chains, lengths, seed, truthK$k0)
  } else {
    res <- fitBart(design, link, chains, lengths, seed, truthK$k0)
    if (task == "s9") {
      pr <- fitPropensityS9(
        design,
        chains,
        lengths,
        seedOf("S9-prop", cell, rep, arm),
        truthK$ke0
      )
      parts <- list(psArr = pr$arr)
      res$seconds <- res$seconds + pr$seconds
      res$kPs <- pr$k
    }
  }
  indices <- if (arm == "longOracle") {
    list(default = seq_len(lengths[2L]))
  } else {
    readingIndices(cell, link, arm)
  }
  scoring <- system.time({
    readings <- lapply(indices, function(idx) {
      L <- poolDraws(res$arr, idx)
      p <- NULL
      if (task == "s9") {
        p <- list(psLatent = poolDraws(parts$psArr, idx))
      } else if (task == "s4b") {
        p <- list(
          ranef = poolDraws(parts$ranef, idx),
          tau = poolMatrix(parts$tau, idx)
        )
      }
      cl <- classesFrom(L, design, link, p)
      list(
        summary = summarizeAll(cl),
        mixing = c(
          k = mixingOf(if (is.null(res$k)) NULL else res$k[idx, , drop = FALSE])
        ),
        draws = length(idx)
      )
    })
  })[["elapsed"]]
  out <- list(
    cell = cell,
    link = link,
    rep = as.integer(rep),
    arm = arm,
    build = build(),
    lib = find.package("dbarts"),
    seconds = res$seconds,
    scoringSeconds = scoring,
    readings = readings,
    quick = quick
  )
  saveAtomic(out, fitFile(cell, link, rep, arm, discriminate))
  cat(sprintf(
    "fit %s %s %d %s: %.1f s, scored %.1f s\n",
    cell,
    link,
    as.integer(rep),
    arm,
    res$seconds,
    out$scoringSeconds
  ))
  invisible(out)
}

## -------------------------------------------------------------- SBC scores

# The interior label of every row: inside [0.01, 0.99] on the probability
# scale by the shipped arm's long reading (its 5th and 95th percentiles), pooled
# with the move's where the move runs (the two arms' percentiles averaged, as
# only summaries are kept), a function of the data and the chains, never of the
# truth or the oracle.
interiorLabels <- function(shippedFit, moveFit = NULL) {
  ends <- function(fit) {
    s <- fit$readings$long$summary
    q <- cbind(s$fInt$q, s$fTail$q)
    cols <- c(s$fInt$cols, s$fTail$cols)
    o <- order(cols)
    list(
      lo = q[match(0.05, round(probsKept, 6)), o],
      hi = q[match(0.95, round(probsKept, 6)), o],
      cols = cols[o]
    )
  }
  e <- ends(shippedFit)
  lo <- e$lo
  hi <- e$hi
  cols <- e$cols
  if (!is.null(moveFit)) {
    m <- ends(moveFit)
    lo <- (lo + m$lo) / 2
    hi <- (hi + m$hi) / 2
  }
  inv <- linkFunction(shippedFit$link)
  lab <- logical(length(cols))
  lab[cols] <- inv(lo) >= interiorRange[1L] & inv(hi) <= interiorRange[2L]
  lab
}

# Re-split a reading's f summary by labels (fits summarize f with every row in
# fInt, as their labels are not yet known).
relabel <- function(summary, labels) {
  if (is.null(summary$fInt)) {
    return(summary)
  }
  all <- summary$fInt
  pick <- function(keep) {
    list(
      q = all$q[, keep, drop = FALSE],
      mean = all$mean[keep],
      sd = all$sd[keep],
      S = all$S,
      cols = which(keep)
    )
  }
  summary$fInt <- pick(labels)
  summary$fTail <- pick(!labels)
  summary
}

# Coverage of the truth per class: counts at 90 and 95 percent, the
# dispersion sum and the interval score (probability scale for latent classes),
# and for a band class whether the whole truth lies inside.
coverageOf <- function(summary, truth, link) {
  inv <- linkFunction(link)
  out <- list()
  for (nm in intersect(names(summary), names(truth))) {
    s <- summary[[nm]]
    t <- truth[[nm]]
    if (isTRUE(classInfo[[nm]]$band)) {
      b <- s$bands
      out[[nm]] <- c(
        n = 1,
        cov90 = as.numeric(all(b[["0.9"]]$lo <= t & t <= b[["0.9"]]$hi)),
        cov95 = as.numeric(all(b[["0.95"]]$lo <= t & t <= b[["0.95"]]$hi)),
        disp = NA_real_,
        is90 = NA_real_
      )
      next
    }
    if (length(t) == 0L) {
      out[[nm]] <- c(n = 0, cov90 = 0, cov95 = 0, disp = 0, is90 = 0)
      next
    }
    i90 <- intervalOf(s, 0.9)
    i95 <- intervalOf(s, 0.95)
    sc <- if (isTRUE(classInfo[[nm]]$latent)) {
      intervalScore(inv(i90$lo), inv(i90$hi), inv(t), 0.9)
    } else {
      intervalScore(i90$lo, i90$hi, t, 0.9)
    }
    v <- s$sd^2
    out[[nm]] <- c(
      n = length(t),
      cov90 = sum(i90$lo <= t & t <= i90$hi),
      cov95 = sum(i95$lo <= t & t <= i95$hi),
      disp = sum(ifelse(v > 0, (t - s$mean)^2 / v, NA), na.rm = TRUE),
      is90 = sum(sc)
    )
  }
  out
}

scoreFile <- function(cell, link, rep, discriminate = "") {
  if (nzchar(discriminate)) {
    return(runDir(
      "scores-disc",
      discriminate,
      link,
      cellFile(cell),
      paste0(repTag(rep), ".rds")
    ))
  }
  runDir("scores", link, cellFile(cell), paste0(repTag(rep), ".rds"))
}

# The SBC readings of one replication: every arm's coverage of the truth,
# split by the shipped arm's labels. A missing arm stops the job.
runScore <- function(cell, link, rep, discriminate = "") {
  design <- loadDesign(cell, link, rep, discriminate)
  # the shipped arm and the oracle always; the move and the long oracle where
  # their jobs ran (the score job waits for every fit job it was given)
  arms <- c("shipped", "oracle")
  if (file.exists(fitFile(cell, link, rep, "move", discriminate))) {
    arms <- c(arms, "move")
  }
  if (file.exists(fitFile(cell, link, rep, "longOracle", discriminate))) {
    arms <- c(arms, "longOracle")
  }
  fits <- lapply(arms, function(a) {
    loadChecked(fitFile(cell, link, rep, a, discriminate))
  })
  names(fits) <- arms
  labels <- if (cellTask(cell) %in% c("prob", "loop", "s4b")) {
    interiorLabels(fits$shipped, fits$move)
  } else {
    NULL
  }
  truth <- truthClasses(design, labels)
  res <- list()
  for (a in arms) {
    for (r in names(fits[[a]]$readings)) {
      summ <- if (is.null(labels)) {
        fits[[a]]$readings[[r]]$summary
      } else {
        relabel(fits[[a]]$readings[[r]]$summary, labels)
      }
      res[[paste(a, r, sep = ".")]] <- coverageOf(summ, truth, link)
    }
  }
  out <- list(
    cell = cell,
    link = link,
    rep = as.integer(rep),
    coverage = res,
    labels = labels,
    rate = design$rate,
    score = design$score,
    k0 = design$truth$k0,
    kMixing = lapply(fits, function(f) f$readings$default$mixing),
    seconds = vapply(fits, function(f) f$seconds, 0),
    quick = quick
  )
  saveAtomic(out, scoreFile(cell, link, rep, discriminate))
  cat(sprintf(
    "score %s %s %d: %s\n",
    cell,
    link,
    as.integer(rep),
    paste(names(res), collapse = " ")
  ))
  invisible(out)
}

## ------------------------------------------------------------- references

referenceDir <- function(cell, link, rep) {
  runDir("ref", link, sprintf("%s-%s", cellFile(cell), repTag(rep)))
}

refPartsPerKernel <- function(cell) {
  refChains(cellChains(cell)) %/% refPartChains
}

# One part of a reference: refPartChains chains of this build's kernel through
# the sampler API: refBurn sweeps of burn-in, refRecord sweeps recorded every
# refThin, then refTail sweeps unthinned (the tails). Latent single-precision
# draws of every raw column; k; the shipped parts keep their sampler so an
# escalation can continue the chains.
runReferencePart <- function(cell, link, rep, part) {
  which <- build()
  if (!which %in% c("shipped", "move")) {
    stop("a reference is made by the two kernels' builds")
  }
  design <- blind(loadDesign(cell, link, rep))
  task <- cellTask(cell)
  dir <- ensureDir(referenceDir(cell, link, rep))
  seed <- seedOf("reference", cell, link, rep, which, part)
  t0 <- proc.time()[["elapsed"]]
  if (task == "s4b") {
    res <- referenceStan4bart(design, seed)
  } else {
    res <- referenceSampler(design, link, task, seed, part, which, dir)
  }
  res$seconds <- proc.time()[["elapsed"]] - t0
  writeSingle(
    res$record,
    file.path(dir, sprintf("%s-part%d-record.f32", which, part))
  )
  writeSingle(
    res$tail,
    file.path(dir, sprintf("%s-part%d-tail.f32", which, part))
  )
  meta <- res[setdiff(names(res), c("record", "tail"))]
  meta$dimsRecord <- dim(res$record)
  meta$dimsTail <- dim(res$tail)
  meta$which <- which
  meta$part <- as.integer(part)
  meta$lib <- find.package("dbarts")
  meta$quick <- quick
  saveAtomic(meta, file.path(dir, sprintf("%s-part%d.rds", which, part)))
  cat(sprintf(
    "reference %s %s %d %s part %d: %.1f s\n",
    cell,
    link,
    as.integer(rep),
    which,
    as.integer(part),
    res$seconds
  ))
}

# The raw columns of a reference as chains x draws x columns: the sampler's
# training fits and test fits in a fit's own order. S9 adds the propensity
# sampler, paired draw for draw; E3's reference is plain fits on the observed
# rows, its masked and held-out rows predicted.
referenceSampler <- function(design, link, task, seed, part, which, dir) {
  build1 <- function(train, test, family, s) {
    args <- list(
      y ~ .,
      data = train,
      family = family,
      n.trees = nTrees,
      n.chains = refPartChains,
      n.threads = 1L,
      n.burn = mcmc$refBurn,
      n.samples = mcmc$refRecord,
      n.thin = mcmc$refThin,
      seed = s,
      verbose = FALSE,
      samplerOnly = TRUE,
      keepTrainingFits = TRUE
    )
    if (!is.null(test)) {
      args$test <- test
    }
    do.call(bart, args)
  }
  colsOf <- function(r, rowsTrain = TRUE) {
    a <- if (rowsTrain) r$train else NULL
    if (!is.null(r$test)) {
      a <- if (is.null(a)) r$test else abind3(a, r$test)
    }
    a
  }
  order <- NULL
  if (task %in% c("prob", "loop")) {
    if (task == "loop") {
      obs <- setdiff(seq_along(design$y), design$mask)
      train <- data.frame(design$x[obs, , drop = FALSE], y = design$y[obs])
      test <- rbind(design$x[design$mask, , drop = FALSE], design$xh)
      order <- order(c(
        obs,
        design$mask,
        length(design$y) +
          seq_len(if (is.null(design$xh)) 0L else nrow(design$xh))
      ))
    } else {
      train <- data.frame(design$x, y = design$y)
      test <- design$xh
    }
    s <- build1(train, test, link, seed)
    pick <- function(r) colsOf(r, TRUE)
  } else {
    base <- list(x = design$x, ps = if (task == "s9") NULL else design$ps)
    if (task == "multi") {
      arm <- factor(design$arm, levels = 0:3)
      train <- data.frame(design$x, arm = arm, y = design$y)
      test <- do.call(
        rbind,
        lapply(0:3, function(a) {
          data.frame(
            design$x,
            arm = factor(rep(a, nrow(design$x)), levels = 0:3)
          )
        })
      )
    } else {
      train <- data.frame(outcomeFrame(base, design$z), y = design$y)
      test <- rbind(outcomeFrame(base, 1L), outcomeFrame(base, 0L))
    }
    s <- build1(train, test, link, seed)
    pick <- function(r) colsOf(r, FALSE)
  }
  sp <- NULL
  if (task == "s9") {
    sp <- build1(data.frame(design$x, y = design$z), NULL, "probit", seed + 7L)
  }
  r1 <- s$run()
  ctrl <- s$control
  ctrl@n.thin <- 1L
  s$setControl(ctrl)
  r2 <- s$run(0L, mcmc$refTail)
  record <- pick(r1)
  tail <- pick(r2)
  extra <- list()
  if (!is.null(sp)) {
    p1 <- sp$run()
    ctrl <- sp$control
    ctrl@n.thin <- 1L
    sp$setControl(ctrl)
    p2 <- sp$run(0L, mcmc$refTail)
    record <- abind3(record, p1$train)
    tail <- abind3(tail, p2$train)
    extra$kPs <- p1$k
  }
  if (!is.null(order)) {
    record <- record[order, , , drop = FALSE]
    tail <- tail[order, , , drop = FALSE]
  }
  if (which == "shipped") {
    s$storeState()
    saveRDS(s, file.path(dir, sprintf("shipped-part%d-sampler.rds", part)))
    if (!is.null(sp)) {
      sp$storeState()
      saveRDS(sp, file.path(dir, sprintf("shipped-part%d-psampler.rds", part)))
    }
  }
  c(
    list(
      record = aperm(record, c(3L, 2L, 1L)),
      tail = aperm(tail, c(3L, 2L, 1L)),
      k = r1$k,
      kTail = r2$k,
      chains = refPartChains
    ),
    extra
  )
}

abind3 <- function(a, b) {
  # rows x draws x chains, bound along rows
  out <- array(0, c(dim(a)[1L] + dim(b)[1L], dim(a)[2L], dim(a)[3L]))
  out[seq_len(dim(a)[1L]), , ] <- a
  out[dim(a)[1L] + seq_len(dim(b)[1L]), , ] <- b
  out
}

# stan4bart has no continuation: a reference part is one long stan4bart run of
# refPartChains chains, every iteration's latent fits (BART plus group
# intercepts), group intercepts and scale kept, thinned afterwards; the last
# refTail iterations are the tail. The plan records them through stan4bart's
# callback; a fit with a callback keeps nothing else, k included, so the run
# keeps its fits instead (the same draws).
referenceStan4bart <- function(design, seed) {
  f <- fitStan4bart(
    design,
    refPartChains,
    c(mcmc$s4bRefWarmup, mcmc$s4bRefIter),
    seed
  )
  a <- s4bArrays(f$fit, refPartChains)
  if (is.null(a$k)) {
    stop("a stan4bart reference draws k")
  }
  iters <- dim(a$arr)[2L]
  recIdx <- seq(mcmc$refThin, iters - mcmc$refTail, by = mcmc$refThin)
  tailIdx <- seq.int(iters - mcmc$refTail + 1L, iters)
  G <- dim(a$ranef)[3L]
  pick <- function(idx) {
    out <- array(0, c(refPartChains, length(idx), dim(a$arr)[3L] + G + 1L))
    out[,, seq_len(dim(a$arr)[3L])] <- a$arr[, idx, , drop = FALSE]
    out[,, dim(a$arr)[3L] + seq_len(G)] <- a$ranef[, idx, , drop = FALSE]
    out[,, dim(a$arr)[3L] + G + 1L] <- a$tau[, idx, drop = FALSE]
    out
  }
  list(
    record = pick(recIdx),
    tail = pick(tailIdx),
    k = a$k[recIdx, , drop = FALSE],
    kTail = a$k[tailIdx, , drop = FALSE],
    chains = refPartChains,
    groups = G
  )
}

# Escalation: the shipped chains continue from their saved samplers to
# escalateRecord recorded sweeps, then a new tail; only when the acceptance
# projected a pass. Otherwise a no-op.
escalateReferencePart <- function(cell, link, rep, part) {
  requireBuild("shipped")
  dir <- referenceDir(cell, link, rep)
  acc <- loadChecked(file.path(dir, "acceptance.rds"))
  if (!identical(acc$decision, "escalate")) {
    cat("no escalation:", acc$decision, "\n")
    return(invisible(NULL))
  }
  if (cellTask(cell) == "s4b") {
    stop(
      "a stan4bart reference escalates by a longer rerun of its shipped half; not built"
    )
  }
  metaPath <- file.path(dir, sprintf("shipped-part%d.rds", part))
  meta <- loadChecked(metaPath)
  if (isTRUE(meta$escalated)) {
    return(invisible(NULL))
  }
  s <- readRDS(file.path(dir, sprintf("shipped-part%d-sampler.rds", part)))
  ctrl <- s$control
  ctrl@n.thin <- mcmc$refThin
  s$setControl(ctrl)
  more <- (mcmc$escalateRecord - mcmc$refRecord) %/% mcmc$refThin
  r1 <- s$run(0L, more)
  ctrl@n.thin <- 1L
  s$setControl(ctrl)
  r2 <- s$run(0L, mcmc$refTail)
  task <- cellTask(cell)
  pick <- function(r) {
    if (task %in% c("prob", "loop")) {
      if (is.null(r$test)) r$train else abind3(r$train, r$test)
    } else {
      r$test
    }
  }
  if (task == "loop") {
    stop("E3 escalation not built")
  }
  rec <- aperm(pick(r1), c(3L, 2L, 1L))
  old <- readSingle(file.path(dir, partFile(meta, "record")), meta$dimsRecord)
  both <- array(0, c(dim(old)[1L], dim(old)[2L] + dim(rec)[2L], dim(old)[3L]))
  both[, seq_len(dim(old)[2L]), ] <- old
  both[, dim(old)[2L] + seq_len(dim(rec)[2L]), ] <- rec
  # new files under new names, the metadata switched to them last, so a killed
  # job leaves the first pass readable
  meta$recordFile <- sprintf("shipped-part%d-record-esc.f32", part)
  meta$tailFile <- sprintf("shipped-part%d-tail-esc.f32", part)
  writeSingle(both, file.path(dir, meta$recordFile))
  writeSingle(aperm(pick(r2), c(3L, 2L, 1L)), file.path(dir, meta$tailFile))
  meta$dimsRecord <- dim(both)
  meta$dimsTail <- dim(aperm(pick(r2), c(3L, 2L, 1L)))
  meta$k <- rbind(meta$k, r1$k)
  meta$kTail <- r2$k
  meta$escalated <- TRUE
  saveAtomic(meta, metaPath)
}

partFile <- function(meta, what) {
  f <- meta[[paste0(what, "File")]]
  if (is.null(f)) {
    sprintf("%s-part%d-%s.f32", meta$which, meta$part, what)
  } else {
    f
  }
}

# A reference's draws: per kernel, chains x draws x raw columns (record and
# tail), k, and the parts' metadata.
loadReference <- function(cell, link, rep) {
  dir <- referenceDir(cell, link, rep)
  out <- list()
  for (which in c("shipped", "move")) {
    rec <- list()
    tl <- list()
    ks <- list()
    kt <- list()
    for (part in seq_len(refPartsPerKernel(cell))) {
      meta <- loadChecked(file.path(dir, sprintf("%s-part%d.rds", which, part)))
      rec[[part]] <- readSingle(
        file.path(dir, partFile(meta, "record")),
        meta$dimsRecord
      )
      tl[[part]] <- readSingle(
        file.path(dir, partFile(meta, "tail")),
        meta$dimsTail
      )
      ks[[part]] <- as.matrix(meta$k)
      kt[[part]] <- as.matrix(meta$kTail)
    }
    bindChains <- function(l) {
      d <- dim(l[[1L]])
      a <- array(0, c(sum(vapply(l, function(x) dim(x)[1L], 0L)), d[2L], d[3L]))
      at <- 0L
      for (x in l) {
        a[at + seq_len(dim(x)[1L]), , ] <- x
        at <- at + dim(x)[1L]
      }
      a
    }
    out[[which]] <- list(
      record = bindChains(rec),
      tail = bindChains(tl),
      k = do.call(cbind, ks),
      kTail = do.call(cbind, kt)
    )
  }
  out
}

# The reference's classes for a set of chains and draws, from the raw columns.
referenceClasses <- function(
  ref,
  which,
  design,
  link,
  labels = NULL,
  tail = FALSE,
  chains = NULL
) {
  r <- ref[[which]]
  A <- if (tail) r$tail else r$record
  if (is.null(chains)) {
    chains <- seq_len(dim(A)[1L])
  }
  L <- poolDraws(A, seq_len(dim(A)[2L]), chains)
  task <- design$task
  parts <- NULL
  if (task == "s9") {
    n <- nrow(design$x)
    parts <- list(psLatent = L[, 2L * n + seq_len(n), drop = FALSE])
    L <- L[, seq_len(2L * n), drop = FALSE]
  }
  if (task == "s4b") {
    n <- nrow(design$x)
    G <- nlevels(design$group)
    parts <- list(
      ranef = L[, n + seq_len(G), drop = FALSE],
      tau = L[, n + G + 1L]
    )
    L <- L[, seq_len(n), drop = FALSE]
  }
  classesFrom(L, design, link, parts, labels)
}

# Labels for a reference replication: the pooled reference's own intervals.
referenceLabels <- function(ref, design, link) {
  if (!design$task %in% c("prob", "loop", "s4b")) {
    return(NULL)
  }
  cl <- lapply(c("shipped", "move"), function(w) {
    referenceClasses(ref, w, design, link)$fInt
  })
  P <- linkFunction(link)(rbind(cl[[1L]], cl[[2L]]))
  interiorColumns(P)
}

acceptanceFile <- function(cell, link, rep) {
  file.path(referenceDir(cell, link, rep), "acceptance.rds")
}

# Acceptance (the plan's reading of dec-B462's "agree and converge"): the
# halves agree on every scored quantity's mean, 5th and 95th percentiles and
# on k's 5th to 95th percentiles within 4.5 combined Monte Carlo standard
# errors; split-Rhat over all chains at most 1.01 for every scalar estimand and
# 1.02 at the 99th percentile of row quantities; tail effective draws at least
# 400 for scalars and at the 1st percentile of row quantities; k's tail
# effective draws at least 100. A failing first pass escalates only when its
# effective draws per sweep project a pass at escalateRecord sweeps.
acceptReference <- function(cell, link, rep, final = FALSE) {
  design <- blind(loadDesign(cell, link, rep))
  ref <- loadReference(cell, link, rep)
  labels <- referenceLabels(ref, design, link)
  cs <- referenceClasses(ref, "shipped", design, link, labels)
  cm <- referenceClasses(ref, "move", design, link, labels)
  nS <- dim(ref$shipped$record)[1L]
  perChainS <- dim(ref$shipped$record)[2L]
  perChainM <- dim(ref$move$record)[2L]
  arr3 <- function(m, chains, per) {
    a <- array(m, c(per, chains, ncol(m)))
    a
  }
  rows <- list()
  for (nm in names(cs)) {
    if (isTRUE(classInfo[[nm]]$band) || ncol(cs[[nm]]) == 0L) {
      next
    }
    a <- arr3(cs[[nm]], nS, perChainS)
    b <- arr3(cm[[nm]], dim(ref$move$record)[1L], perChainM)
    for (j in seq_len(dim(a)[3L])) {
      x <- a[,, j]
      y <- b[,, j]
      st <- function(v) {
        c(
          mean(v),
          quantile(v, c(0.05, 0.95), names = FALSE),
          posterior::mcse_mean(v),
          posterior::mcse_quantile(v, 0.05),
          posterior::mcse_quantile(v, 0.95)
        )
      }
      sx <- st(x)
      sy <- st(y)
      z <- abs(sx[1:3] - sy[1:3]) / sqrt(sx[4:6]^2 + sy[4:6]^2)
      both <- cbind(x, y)
      rows[[length(rows) + 1L]] <- data.frame(
        class = nm,
        j = j,
        row = isTRUE(classInfo[[nm]]$row),
        zMax = max(z, na.rm = TRUE),
        rhat = posterior::rhat(both),
        essTail = posterior::ess_tail(both),
        essS = posterior::ess_tail(x)
      )
    }
  }
  tab <- do.call(rbind, rows)
  kS <- ref$shipped$k
  kM <- ref$move$k
  kq <- c(0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95)
  if (is.null(kS)) {
    kZ <- 0
    kEss <- Inf
    kRhat <- NA_real_
    kEssS <- Inf
  } else {
    kDiff <- quantile(as.vector(kS), kq, names = FALSE) -
      quantile(as.vector(kM), kq, names = FALSE)
    kSe <- sqrt(
      vapply(kq, function(p) posterior::mcse_quantile(kS, p), 0)^2 +
        vapply(kq, function(p) posterior::mcse_quantile(kM, p), 0)^2
    )
    kZ <- max(abs(kDiff) / kSe)
    kEss <- posterior::ess_tail(cbind(kS, kM))
    kEssS <- posterior::ess_tail(kS)
    kRhat <- posterior::rhat(cbind(kS, kM))
  }
  scal <- tab[!tab$row, ]
  rowq <- tab[tab$row, ]
  checks <- c(
    halvesAgree = all(tab$zMax <= 4.5, na.rm = TRUE),
    kAgree = is.finite(kZ) && kZ <= 4.5,
    rhatScalar = nrow(scal) == 0L || all(scal$rhat <= 1.01, na.rm = TRUE),
    rhatRows = nrow(rowq) == 0L ||
      quantile(rowq$rhat, 0.99, na.rm = TRUE) <= 1.02,
    essScalar = nrow(scal) == 0L || all(scal$essTail >= 400, na.rm = TRUE),
    essRows = nrow(rowq) == 0L ||
      quantile(rowq$essTail, 0.01, na.rm = TRUE) >= 400,
    kEss = kEss >= 100
  )
  accepted <- all(checks)
  decision <- if (accepted) {
    "accepted"
  } else if (final) {
    "failed"
  } else {
    # projection at escalateRecord recorded sweeps: the shipped half's
    # effective draws scale with its sweeps, its standard errors by their root
    f <- mcmc$escalateRecord / mcmc$refRecord
    projEss <- function(total, shipped) total + (f - 1) * shipped
    pass <- c(
      kAgree = checks[["kAgree"]] ||
        (is.finite(kZ) && kZ * sqrt(2 / (1 + 1 / f)) <= 4.5),
      kEss = checks[["kEss"]] || projEss(kEss, kEssS) >= 100,
      essScalar = checks[["essScalar"]] ||
        all(projEss(scal$essTail, scal$essS) >= 400, na.rm = TRUE),
      essRows = checks[["essRows"]] ||
        quantile(projEss(rowq$essTail, rowq$essS), 0.01, na.rm = TRUE) >= 400,
      halvesAgree = checks[["halvesAgree"]],
      rhatScalar = checks[["rhatScalar"]],
      rhatRows = checks[["rhatRows"]]
    )
    if (
      all(pass) && cellTask(cell) %in% c("prob", "effect", "multi", "d2", "s9")
    ) {
      "escalate"
    } else {
      "failed"
    }
  }
  out <- list(
    cell = cell,
    link = link,
    rep = as.integer(rep),
    final = final,
    accepted = accepted,
    decision = decision,
    checks = checks,
    kZ = kZ,
    kEss = kEss,
    kRhat = kRhat,
    labels = labels,
    quantities = nrow(tab),
    maxZ = max(tab$zMax, na.rm = TRUE),
    rhat99 = if (nrow(rowq)) quantile(rowq$rhat, 0.99, na.rm = TRUE) else NA,
    quick = quick
  )
  if (final || decision != "escalate") {
    out$final <- TRUE
  }
  # the quantiles the fixed-k arms use: the pooled reference's k, or the move
  # half alone when the reference fails (flagged)
  if (isTRUE(out$final)) {
    kPool <- if (accepted) c(as.vector(kS), as.vector(kM)) else as.vector(kM)
    out$kQuantiles <- if (is.null(kPool)) {
      NULL
    } else {
      quantile(kPool, fixedQuantiles, names = FALSE)
    }
    names(out$kQuantiles) <- names(fixedQuantiles)
    out$kFromMoveHalf <- !accepted
  }
  saveAtomic(out, acceptanceFile(cell, link, rep))
  cat(sprintf(
    "accept %s %s %d: %s (%s)\n",
    cell,
    link,
    as.integer(rep),
    decision,
    paste(names(checks)[!checks], collapse = ",")
  ))
  invisible(out)
}

finalAcceptance <- function(cell, link, rep) {
  acc <- loadChecked(acceptanceFile(cell, link, rep))
  if (!isTRUE(acc$final)) {
    stop("reference ", cell, " ", link, " ", rep, " has no final acceptance")
  }
  acc
}

## --------------------------------------------------------- reference mass

# Each reference class column sorted once; the mass of an interval is then two
# binary searches.
sortedColumns <- function(M) {
  if (ncol(M) == 0L) {
    return(M)
  }
  matrix(apply(M, 2L, sort.int, method = "quick"), nrow(M), ncol(M))
}

massOf <- function(sorted, lo, hi) {
  S <- nrow(sorted)
  vapply(
    seq_len(ncol(sorted)),
    function(j) {
      (findInterval(hi[j], sorted[, j]) -
        findInterval(lo[j], sorted[, j], left.open = TRUE)) /
        S
    },
    0
  )
}

# A reading's readings against a reference made of one or both kernels'
# halves (each half its classes and their sorted columns; a pooled reference
# is never copied, its counts are the halves' counts summed): per class, the
# mean reference mass at 90 and 95 percent, the location error and the width
# ratio (for a pooled reference, against the halves' averaged percentiles), and
# for a band class the simultaneous mass.
massCount <- function(sorted, lo, hi) {
  vapply(
    seq_len(ncol(sorted)),
    function(j) {
      findInterval(hi[j], sorted[, j]) -
        findInterval(lo[j], sorted[, j], left.open = TRUE)
    },
    0
  )
}

massReading <- function(summary, halves) {
  out <- list()
  for (nm in intersect(names(summary), names(halves[[1L]]$classes))) {
    s <- summary[[nm]]
    Rs <- lapply(halves, function(h) h$classes[[nm]])
    m <- ncol(Rs[[1L]])
    S <- sum(vapply(Rs, nrow, 0L))
    if (m == 0L) {
      out[[nm]] <- c(n = 0, mass90 = NA, mass95 = NA, loc = NA, width = NA)
      next
    }
    if (isTRUE(classInfo[[nm]]$band)) {
      inside <- function(b) {
        sum(vapply(
          Rs,
          function(R) {
            ok <- sweep(R, 2L, b$lo, ">=") & sweep(R, 2L, b$hi, "<=")
            if (nm == "bandUnit") {
              a <- array(ok, c(nrow(R), 6L, ncol(R) %/% 6L))
              u <- a[, 1L, , drop = FALSE]
              for (g in 2:6) {
                u <- u & a[, g, , drop = FALSE]
              }
              sum(u)
            } else {
              sum(rowSums(!ok) == 0L)
            }
          },
          0
        )) /
          (if (nm == "bandUnit") S * (m %/% 6L) else S)
      }
      out[[nm]] <- c(
        n = 1,
        mass90 = inside(s$bands[["0.9"]]),
        mass95 = inside(s$bands[["0.95"]]),
        loc = NA,
        width = NA
      )
      next
    }
    i90 <- intervalOf(s, 0.9)
    i95 <- intervalOf(s, 0.95)
    count <- function(iv) {
      Reduce(
        `+`,
        lapply(halves, function(h) massCount(h$sorted[[nm]], iv$lo, iv$hi))
      ) /
        S
    }
    sums <- Reduce(`+`, lapply(Rs, colSums))
    squares <- Reduce(`+`, lapply(Rs, function(R) colSums(R^2)))
    rm <- sums / S
    rsd <- sqrt(pmax(squares / S - rm^2, 0) * S / (S - 1))
    rq <- Reduce(
      `+`,
      lapply(halves, function(h) {
        srt <- h$sorted[[nm]]
        srt[pmax(1L, round(c(0.05, 0.95) * nrow(srt))), , drop = FALSE]
      })
    ) /
      length(halves)
    out[[nm]] <- c(
      n = m,
      mass90 = mean(count(i90)),
      mass95 = mean(count(i95)),
      loc = mean(ifelse(rsd > 0, abs(s$mean - rm) / rsd, NA), na.rm = TRUE),
      width = mean(
        ifelse(
          rq[2L, ] > rq[1L, ],
          (i90$hi - i90$lo) / (rq[2L, ] - rq[1L, ]),
          NA
        ),
        na.rm = TRUE
      )
    )
  }
  out
}

# Both kernels' halves of a reference, each with its classes' sorted columns
# (band classes are read from the draws and not sorted).
halfOf <- function(cl) {
  srt <- lapply(names(cl), function(nm) {
    if (isTRUE(classInfo[[nm]]$band)) NULL else sortedColumns(cl[[nm]])
  })
  names(srt) <- names(cl)
  list(classes = cl, sorted = srt)
}

halfClasses <- function(ref, design, link, labels) {
  list(
    shipped = halfOf(referenceClasses(ref, "shipped", design, link, labels)),
    move = halfOf(referenceClasses(ref, "move", design, link, labels))
  )
}

# The halves a reference named "shipped", "move" or "pooled" is made of.
halvesOf <- function(H, which) {
  if (which == "pooled") H[c("shipped", "move")] else H[which]
}

# The same reference restricted to some columns of each class.
subsetHalves <- function(halves, keepCols) {
  lapply(halves, function(h) {
    cl <- lapply(names(keepCols), function(nm) {
      h$classes[[nm]][, keepCols[[nm]], drop = FALSE]
    })
    srt <- lapply(names(keepCols), function(nm) {
      if (is.null(h$sorted[[nm]])) {
        NULL
      } else {
        h$sorted[[nm]][, keepCols[[nm]], drop = FALSE]
      }
    })
    names(cl) <- names(keepCols)
    names(srt) <- names(keepCols)
    list(classes = cl, sorted = srt)
  })
}

refscoreFile <- function(cell, link, rep) {
  runDir("refscore", link, cellFile(cell), paste0(repTag(rep), ".rds"))
}

# The reference readings of one reference replication, after its final
# acceptance: each kernel's default, equal-time and long readings and its
# stationary tails against the other kernel's half; the oracle against the
# pooled reference (the move half when the reference failed); every reading
# also against the move half alone (the 80 percent rule's sensitivity reading)
# and the pooled reference, and each tail against the pooled reference less its
# own chains (cells of at most pooledLessCap quantities). B's six seeds are all
# read against its one reference. Every entry records the half it was scored
# against, for `check halves`.
pooledLessCap <- 3000L

runRefscore <- function(cell, link, rep, seedCount = NA_integer_) {
  acc <- finalAcceptance(cell, link, rep)
  design <- blind(loadDesign(cell, link, rep))
  ref <- loadReference(cell, link, rep)
  labels <- acc$labels
  H <- halfClasses(ref, design, link, labels)
  other <- c(shipped = "move", move = "shipped", mixed = "move")
  scoreAgainst <- function(summ, vs) {
    list(
      vsOther = massReading(summ, halvesOf(H, vs)),
      vsMove = massReading(summ, halvesOf(H, "move")),
      vsPooled = massReading(summ, halvesOf(H, "pooled")),
      against = vs
    )
  }
  res <- list()
  seeds <- if (cellRegime(cell)$id == "B") {
    seq_len(if (is.na(seedCount)) repsFor(cell) else seedCount)
  } else {
    as.integer(rep)
  }
  for (seedRep in seeds) {
    for (arm in c("shipped", "move", "oracle", "mixed")) {
      path <- fitFile(cell, link, seedRep, arm)
      if (arm %in% c("oracle", "mixed") && !file.exists(path)) {
        next
      }
      f <- loadChecked(path)
      for (r in names(f$readings)) {
        summ <- if (is.null(labels)) {
          f$readings[[r]]$summary
        } else {
          relabel(f$readings[[r]]$summary, labels)
        }
        vs <- if (arm == "oracle") {
          (if (acc$accepted) "pooled" else "move")
        } else {
          other[[arm]]
        }
        key <- paste(arm, r, sep = ".")
        if (cellRegime(cell)$id == "B") {
          key <- paste(key, repTag(seedRep), sep = ".")
        }
        res[[key]] <- scoreAgainst(summ, vs)
      }
    }
  }
  # the tails: default-length runs of each kernel begun at stationarity, as
  # many runs of the user's chain count as the half holds
  C <- cellChains(cell)
  totalCols <- sum(vapply(H$shipped$classes, ncol, 0L))
  for (which in c("shipped", "move")) {
    nC <- dim(ref[[which]]$tail)[1L]
    runs <- split(seq_len(nC), ceiling(seq_len(nC) / C))
    runs <- Filter(function(r) length(r) == C, runs)
    if (length(runs) == 0L) {
      # a half of fewer chains than the user's fit (none at full size) is one run
      runs <- list(seq_len(nC))
    }
    vals <- lapply(runs, function(ch) {
      cl <- referenceClasses(
        ref,
        which,
        design,
        link,
        labels,
        tail = TRUE,
        chains = ch
      )
      summ <- summarizeAll(cl)
      out <- scoreAgainst(summ, other[[which]])
      if (totalCols <= pooledLessCap) {
        keep <- setdiff(seq_len(dim(ref[[which]]$record)[1L]), ch)
        # a 10-chain cell's tail run is its whole half: the rest is the other half
        rest <- halvesOf(H, other[[which]])
        if (length(keep)) {
          rest <- c(
            list(
              less = halfOf(referenceClasses(
                ref,
                which,
                design,
                link,
                labels,
                chains = keep
              ))
            ),
            rest
          )
        }
        out$vsPooledLess <- massReading(summ, rest)
      }
      out
    })
    res[[paste0(which, ".tails")]] <- averageReadings(vals)
  }
  out <- list(
    cell = cell,
    link = link,
    rep = as.integer(rep),
    accepted = acc$accepted,
    readings = res,
    labels = labels,
    quick = quick
  )
  saveAtomic(out, refscoreFile(cell, link, rep))
  cat(sprintf(
    "refscore %s %s %d: %d readings\n",
    cell,
    link,
    as.integer(rep),
    length(res)
  ))
  invisible(out)
}

averageReadings <- function(vals) {
  if (length(vals) == 1L) {
    return(vals[[1L]])
  }
  out <- vals[[1L]]
  for (h in setdiff(names(out), "against")) {
    for (nm in names(out[[h]])) {
      out[[h]][[nm]] <- Reduce(`+`, lapply(vals, function(v) v[[h]][[nm]])) /
        length(vals)
    }
  }
  out
}

## ---------------------------------------------------------------- fixed k

fixedkFile <- function(cell, link, rep) {
  runDir("fixedk", link, cellFile(cell), paste0(repTag(rep), ".rds"))
}

# The fixed-k arms of a reference replication, after its final acceptance:
# k at q10, q50 and q90 with every chain, at q05, q30, q70 and q95 with one,
# default length; the ruled three read alone and 200 pooled stuck runs (as
# many chains as the user's fit, drawn without replacement, each chain's chance
# its quantile's share of k's posterior over its chain count), all against the
# pooled reference. Row classes keep at most stuckRowCap quantities in the
# pooled runs.
runFixedk <- function(cell, link, rep) {
  requireBuild("shipped")
  acc <- finalAcceptance(cell, link, rep)
  if (is.null(acc$kQuantiles)) {
    stop("no k quantiles for ", cell)
  }
  design <- blind(loadDesign(cell, link, rep))
  ref <- loadReference(cell, link, rep)
  labels <- acc$labels
  H <- halfClasses(ref, design, link, labels)
  # a failed reference's fixed-k readings are scored against its move half, flagged
  against <- if (acc$accepted) "pooled" else "move"
  C <- cellChains(cell)
  task <- cellTask(cell)
  lengths <- if (task == "s4b") {
    c(mcmc$s4bWarmup, mcmc$s4bDefault)
  } else {
    c(mcmc$burn, mcmc$defaultDraws)
  }
  chainSets <- list()
  seconds <- 0
  for (q in names(fixedQuantiles)) {
    nC <- if (q %in% ruledQuantiles) C else 1L
    k <- acc$kQuantiles[[q]]
    seed <- seedOf("fixedk", cell, link, rep, q)
    if (task == "s4b") {
      f <- fitStan4bart(design, nC, lengths, seed, k)
      a <- s4bArrays(f$fit, nC)
      res <- list(
        arr = a$arr,
        parts = list(ranef = a$ranef, tau = a$tau),
        seconds = f$seconds
      )
    } else if (task == "loop") {
      res <- fitLoop(design, link, nC, lengths, seed, k)
    } else {
      res <- fitBart(design, link, nC, lengths, seed, k)
      if (task == "s9") {
        pr <- fitPropensityS9(design, nC, lengths, seed + 1L, NULL)
        res$parts <- list(psArr = pr$arr)
        res$seconds <- res$seconds + pr$seconds
      }
    }
    seconds <- seconds + res$seconds
    chainSets[[q]] <- res
  }
  classesOfChains <- function(q, chains) {
    res <- chainSets[[q]]
    idx <- seq_len(dim(res$arr)[2L])
    L <- poolDraws(res$arr, idx, chains)
    p <- NULL
    if (task == "s9") {
      p <- list(psLatent = poolDraws(res$parts$psArr, idx, chains))
    } else if (task == "s4b") {
      p <- list(
        ranef = poolDraws(res$parts$ranef, idx, chains),
        tau = poolMatrix(res$parts$tau, idx, chains)
      )
    }
    classesFrom(L, design, link, p, labels)
  }
  refHalves <- halvesOf(H, against)
  ruled <- lapply(ruledQuantiles, function(q) {
    massReading(summarizeAll(classesOfChains(q, seq_len(C))), refHalves)
  })
  names(ruled) <- ruledQuantiles
  # the pooled stuck runs
  pool <- do.call(
    rbind,
    lapply(names(fixedQuantiles), function(q) {
      data.frame(q = q, chain = seq_len(dim(chainSets[[q]]$arr)[1L]))
    })
  )
  pool$w <- quantileBands[pool$q] / ifelse(pool$q %in% ruledQuantiles, C, 1L)
  # each chain's classes kept only at the stuck runs' columns: band classes
  # out, row classes at most stuckRowCap quantities
  first <- classesOfChains(pool$q[1L], pool$chain[1L])
  bandless <- setdiff(
    names(first),
    names(Filter(function(i) isTRUE(i$band), classInfo))
  )
  keepCols <- lapply(bandless, function(nm) {
    m <- ncol(first[[nm]])
    if (isTRUE(classInfo[[nm]]$row) && m > stuckRowCap) {
      sort(withSeed(
        seedOf("stuck-cols", cell, rep, nm),
        sample.int(m, stuckRowCap)
      ))
    } else {
      seq_len(m)
    }
  })
  names(keepCols) <- bandless
  perChain <- lapply(seq_len(nrow(pool)), function(i) {
    cl <- classesOfChains(pool$q[i], pool$chain[i])
    out <- lapply(bandless, function(nm) {
      cl[[nm]][, keepCols[[nm]], drop = FALSE]
    })
    names(out) <- bandless
    out
  })
  subHalves <- subsetHalves(refHalves, keepCols)
  sets <- withSeed(
    seedOf("stuck-sets", cell, link, rep),
    lapply(seq_len(pooledRuns), function(i) {
      sample.int(nrow(pool), C, prob = pool$w)
    })
  )
  stuck <- lapply(sets, function(set) {
    cl <- lapply(bandless, function(nm) {
      do.call(rbind, lapply(set, function(i) perChain[[i]][[nm]]))
    })
    names(cl) <- bandless
    massReading(summarizeAll(cl), subHalves)
  })
  stuckMean <- stuck[[1L]]
  for (nm in names(stuckMean)) {
    stuckMean[[nm]] <- Reduce(`+`, lapply(stuck, `[[`, nm)) / length(stuck)
  }
  out <- list(
    cell = cell,
    link = link,
    rep = as.integer(rep),
    against = against,
    kQuantiles = acc$kQuantiles,
    kFromMoveHalf = acc$kFromMoveHalf,
    ruled = ruled,
    stuck = stuckMean,
    stuckSpread = lapply(names(stuckMean), function(nm) {
      sd(vapply(stuck, function(s) s[[nm]][["mass90"]], 0))
    }),
    seconds = seconds,
    quick = quick
  )
  saveAtomic(out, fixedkFile(cell, link, rep))
  cat(sprintf(
    "fixedk %s %s %d: %.1f s\n",
    cell,
    link,
    as.integer(rep),
    seconds
  ))
  invisible(out)
}

## ------------------------------------------------------------------ timing

timingPieceFile <- function(cell, link, rep, which, t) {
  runDir(
    "timing",
    link,
    cellFile(cell),
    sprintf("%s-%s-t%d.rds", repTag(rep), which, as.integer(t))
  )
}

# One repetition of the timing pass under this build: a sampler-only run of
# timingRun sweeps after timingWarm, one chain, on the design (S9 adds its
# propensity sampler); stan4bart and the loop are timed by total fit time.
runTiming <- function(cell, link, rep, t) {
  which <- build()
  design <- blind(loadDesign(cell, link, rep))
  task <- cellTask(cell)
  seed <- seedOf("timing", cell, link, rep)
  if (task == "s4b") {
    seconds <- fitStan4bart(
      design,
      1L,
      c(mcmc$timingWarm, mcmc$timingRun),
      seed
    )$seconds
  } else if (task == "loop") {
    seconds <- fitLoop(
      design,
      link,
      1L,
      c(mcmc$timingWarm, mcmc$timingRun),
      seed
    )$seconds
  } else {
    timeSampler <- function(train, test, family) {
      args <- list(
        y ~ .,
        data = train,
        family = family,
        n.trees = nTrees,
        n.chains = 1L,
        n.threads = 1L,
        seed = seed,
        verbose = FALSE,
        keepTrainingFits = FALSE,
        samplerOnly = TRUE
      )
      if (!is.null(test)) {
        args$test <- test
      }
      s <- do.call(bart, args)
      invisible(s$run(mcmc$timingWarm, 0L))
      system.time(invisible(s$run(0L, mcmc$timingRun)))[["elapsed"]]
    }
    fd <- fitData(design)
    seconds <- timeSampler(fd$train, fd$test, link)
    if (task == "s9") {
      seconds <- seconds +
        timeSampler(data.frame(design$x, y = design$z), NULL, "probit")
    }
  }
  load <- if (file.exists("/proc/loadavg")) {
    readLines("/proc/loadavg", warn = FALSE)
  } else {
    NA_character_
  }
  saveAtomic(
    list(
      cell = cell,
      link = link,
      rep = as.integer(rep),
      build = which,
      t = as.integer(t),
      seconds = seconds,
      load = load,
      quick = quick
    ),
    timingPieceFile(cell, link, rep, which, t)
  )
  cat(sprintf(
    "timing %s %s %d %s %d: %.3f s\n",
    cell,
    link,
    as.integer(rep),
    which,
    as.integer(t),
    seconds
  ))
}

# r per cell: the median over replications and repetitions of the move's time
# over the shipped build's, each pair adjacent in time; D1, D2 and B take the
# binary interval study's ratios (D1 its A6 cells, D2 its setting 1 cells at
# 1,000 units, moderate and rare, B its B1 and B2 problems).
bisRatioCell <- c(
  "D1:arr1-mod" = "A6:arr1-mod",
  "D1:arr1-sep" = "A6:arr1-sep",
  "D1:arr2-mod" = "A6:arr2-mod",
  "D1:null-mod" = "A6:null-mod",
  "D2:r30" = "A3:r1-mod",
  "D2:r05" = "A3:r1-rare05"
)

summarizeTiming <- function() {
  files <- list.files(
    runDir("timing"),
    pattern = "-t[0-9]+[.]rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  pieces <- Filter(
    function(p) identical(p$quick, quick),
    lapply(files, readRDS)
  )
  d <- data.frame(
    link = vapply(pieces, `[[`, "", "link"),
    cell = vapply(pieces, `[[`, "", "cell"),
    rep = vapply(pieces, `[[`, 0L, "rep"),
    build = vapply(pieces, `[[`, "", "build"),
    t = vapply(pieces, `[[`, 0L, "t"),
    seconds = vapply(pieces, `[[`, 0, "seconds")
  )
  w <- reshape(
    d,
    idvar = c("link", "cell", "rep", "t"),
    timevar = "build",
    direction = "wide"
  )
  w$ratio <- w$seconds.move / w$seconds.shipped
  rows <- lapply(split(w, paste(w$link, w$cell)), function(g) {
    data.frame(
      link = g$link[1L],
      cell = g$cell[1L],
      n = sum(is.finite(g$ratio)),
      r = median(g$ratio, na.rm = TRUE),
      rMin = min(g$ratio, na.rm = TRUE),
      rMax = max(g$ratio, na.rm = TRUE),
      sShipped = median(g$seconds.shipped),
      sMove = median(g$seconds.move),
      source = "measured"
    )
  })
  tab <- do.call(rbind, rows)
  bisRun <- Sys.getenv("KMN_BIS_RUN")
  if (nzchar(bisRun)) {
    bis <- readRDS(file.path(bisRun, "timing", "ratios.rds"))
    extra <- list()
    for (link in both) {
      cells <- c(names(bisRatioCell), regimeCells(regimes$B))
      for (cell in cells) {
        src <- if (cell %in% names(bisRatioCell)) bisRatioCell[[cell]] else cell
        row <- bis[bis$link == link & bis$cell == src, , drop = FALSE]
        if (nrow(row) == 1L) {
          extra[[length(extra) + 1L]] <- data.frame(
            link = link,
            cell = cell,
            n = row$n,
            r = row$r,
            rMin = row$rMin,
            rMax = row$rMax,
            sShipped = NA_real_,
            sMove = NA_real_,
            source = paste("binary interval study", src)
          )
        }
      }
    }
    tab <- rbind(tab, do.call(rbind, extra))
  }
  saveAtomic(tab, timingFile())
  print(tab, row.names = FALSE)
  invisible(tab)
}

## ------------------------------------------------------------------ screen

# Stage 0: on the binary interval study's pilot references (16 chains of the
# scored quantities, k with every draw), for every interior quantity: the share
# of its posterior variance explained by log k in decile bins; the conditional
# mean and 90 percent interval in the bins of k's draws around its 10th and
# 90th percentiles (from its 5th to 15th, and 85th to 95th), as a shift in
# posterior standard deviations and as reference mass; and a proxy of the
# pooled stuck-run reading: 200 runs of four chains, each chain 500 draws from
# the bin around one of the seven quantiles chosen as the stuck runs choose
# them, their 90 percent interval's mass under the whole reference.
screenReference <- function(dir, link) {
  parts <- list.files(
    dir,
    pattern = "^(shipped|move)-part[0-9]+[.]rds$",
    full.names = TRUE
  )
  mats <- list()
  ks <- list()
  meta1 <- NULL
  for (p in parts) {
    meta <- readRDS(p)
    meta1 <- meta
    mats[[length(mats) + 1L]] <- readSingle(
      sub("[.]rds$", ".f32", p),
      meta$dims
    )
    ks[[length(ks) + 1L]] <- as.vector(t(meta$k)) # chain-major, as the draws
  }
  D <- do.call(rbind, mats)
  k <- unlist(ks)
  rm(mats)
  scale <- if (is.null(meta1$scale)) "scored" else meta1$scale
  if (scale == "latent") {
    stop("a latent pilot reference: not expected among the stored B references")
  }
  isProb <- grepl("/B1_", dir)
  interior <- if (isProb) {
    interiorColumns(D)
  } else {
    q <- apply(D, 2L, stats::quantile, probs = c(0.05, 0.95), names = FALSE)
    q[2L, ] > q[1L, ]
  }
  cls <- list()
  if (isProb) {
    cls$fInt <- D[, interior, drop = FALSE]
  } else {
    cls$ite <- D[, interior, drop = FALSE]
    cls$ate <- matrix(rowMeans(D), ncol = 1L)
  }
  lk <- log(k)
  dec <- cut(lk, unique(quantile(lk, seq(0, 1, 0.1))), include.lowest = TRUE)
  qk <- quantile(k, c(0.05, 0.15, 0.85, 0.95))
  lowBin <- k >= qk[1L] & k <= qk[2L]
  highBin <- k >= qk[3L] & k <= qk[4L]
  binAt <- lapply(fixedQuantiles, function(p) {
    b <- quantile(k, pmin(pmax(c(p - 0.05, p + 0.05), 0), 1))
    which(k >= b[1L] & k <= b[2L])
  })
  sets <- withSeed(
    seedOf("screen", dir),
    lapply(seq_len(pooledRuns), function(i) {
      sample(names(fixedQuantiles), 4L, replace = TRUE, prob = quantileBands)
    })
  )
  rows <- lapply(names(cls), function(nm) {
    M <- cls[[nm]]
    if (ncol(M) == 0L) {
      return(NULL)
    }
    srt <- sortedColumns(M)
    v <- apply(M, 2L, stats::var)
    binMeans <- apply(M, 2L, function(col) tapply(col, dec, mean))
    binN <- as.numeric(table(dec))
    between <- colSums(
      binN * sweep(matrix(binMeans, ncol = ncol(M)), 2L, colMeans(M))^2
    ) /
      length(lk)
    r2 <- ifelse(v > 0, between / v, NA)
    sdAll <- sqrt(v)
    shift <- function(bin) {
      ifelse(
        sdAll > 0,
        abs(colMeans(M[bin, , drop = FALSE]) - colMeans(M)) / sdAll,
        NA
      )
    }
    mass <- function(idx) {
      q <- apply(
        M[idx, , drop = FALSE],
        2L,
        stats::quantile,
        probs = c(0.05, 0.95),
        names = FALSE
      )
      massOf(srt, q[1L, ], q[2L, ])
    }
    stuck <- withSeed(
      seedOf("screen-stuck", dir, nm),
      vapply(
        sets,
        function(set) {
          idx <- unlist(lapply(set, function(q) {
            sample(binAt[[q]], 500L, replace = TRUE)
          }))
          mean(mass(idx))
        },
        0
      )
    )
    data.frame(
      class = nm,
      quantities = ncol(M),
      r2 = mean(r2, na.rm = TRUE),
      r2Max = max(r2, na.rm = TRUE),
      shiftLow = mean(shift(lowBin), na.rm = TRUE),
      shiftHigh = mean(shift(highBin), na.rm = TRUE),
      massLow = mean(mass(which(lowBin))),
      massHigh = mean(mass(which(highBin))),
      stuckMass = mean(stuck)
    )
  })
  do.call(rbind, rows)
}

# The pilot datasets' regimes: by name for the prediction regimes, and the
# three real causal datasets for the effect regimes.
screenRegimes <- list(
  S1 = c("tictactoe", "banknote", "wdbc", "ionosphere"),
  S2 = c("pima", "german", "creditapproval"),
  S3 = "adult",
  S4 = c("kyphosis", "birthwt", "haberman", "cleveland", "infert"),
  S6 = c("pima", "german"),
  "S7, S8, S9" = c("lindner", "actg175", "lalonde")
)

runScreen <- function(outfile) {
  bisRun <- Sys.getenv("KMN_BIS_RUN")
  if (!nzchar(bisRun)) {
    stop("KMN_BIS_RUN must name the binary interval study's run directory")
  }
  rows <- list()
  for (link in both) {
    dirs <- list.dirs(file.path(bisRun, "ref", link), recursive = FALSE)
    for (dir in dirs) {
      acc <- tryCatch(
        readRDS(file.path(dir, "acceptance.rds")),
        error = function(e) NULL
      )
      r <- screenReference(dir, link)
      if (is.null(r)) {
        next
      }
      r$link <- link
      r$problem <- sub("-r001$", "", basename(dir))
      r$accepted <- if (is.null(acc)) NA else isTRUE(acc$accepted)
      rows[[length(rows) + 1L]] <- r
      cat(sprintf(
        "screen %s %s: %s\n",
        link,
        basename(dir),
        paste(sprintf("%s %.4f", r$class, 0.9 - r$stuckMass), collapse = ", ")
      ))
    }
  }
  tab <- do.call(rbind, rows)
  tab$dataset <- sub("^B[12]_", "", tab$problem)
  tab$stuckShortfall <- 0.9 - tab$stuckMass
  tab$lowShortfall <- 0.9 - tab$massLow
  tab$highShortfall <- 0.9 - tab$massHigh
  pred <- list()
  for (reg in names(screenRegimes)) {
    for (link in both) {
      sub <- tab[
        tab$link == link & tab$dataset %in% screenRegimes[[reg]],
        ,
        drop = FALSE
      ]
      for (cl in unique(sub$class)) {
        s <- sub[sub$class == cl, , drop = FALSE]
        pred[[length(pred) + 1L]] <- data.frame(
          regime = reg,
          link = link,
          class = cl,
          problems = nrow(s),
          worstStuck = max(s$stuckShortfall),
          meanStuck = mean(s$stuckShortfall),
          worstEnd = max(pmax(s$lowShortfall, s$highShortfall)),
          prediction = if (max(s$stuckShortfall) < 0.01) {
            "bounded"
          } else {
            "not bounded"
          }
        )
      }
    }
  }
  out <- list(table = tab, prediction = do.call(rbind, pred), quick = quick)
  saveAtomic(out, outfile)
  print(out$prediction, row.names = FALSE)
  invisible(out)
}

## ------------------------------------------------------------------ checks

checkFile <- function(...) runDir("checks", ...)

report <- function(name, pass, detail = NULL) {
  cat(sprintf("CHECK %s: %s\n", name, if (isTRUE(pass)) "PASS" else "FAIL"))
  if (!is.null(detail)) {
    print(detail)
  }
  invisible(isTRUE(pass))
}

# SBC self-consistency on a design: the prior-draw sampler's predictors and
# cut points equal to those of the sampler bart() builds for the fit; recorded
# training and test fits equal to predict() at the same state (logistic and
# probit, every design's factors, missing values and test rows); the fixed k0
# handed back as the drawn hyperprior keeps the chain's k and then moves it.
checkSelf <- function(cell, link, rep) {
  design <- blind(loadDesign(cell, link, rep))
  task <- cellTask(cell)
  if (task == "s4b") {
    trainX <- design$x
    test <- NULL
  } else if (task == "prob" || task == "loop") {
    trainX <- design$x
    test <- design$xh
  } else {
    fd <- fitData(design)
    trainX <- fd$train[, setdiff(names(fd$train), "y"), drop = FALSE]
    test <- fd$test
  }
  seed <- seedOf("self", cell, link, rep)
  prior <- priorSampler(trainX, link, seed)
  fitArgs <- list(
    y ~ .,
    data = data.frame(trainX, y = design$y),
    family = link,
    n.trees = nTrees,
    n.chains = 1L,
    n.threads = 1L,
    seed = seed,
    verbose = FALSE,
    samplerOnly = TRUE,
    keepTrainingFits = TRUE
  )
  if (!is.null(test)) {
    fitArgs$test <- test
  }
  fitSampler <- do.call(bart, fitArgs)
  # the predictor matrices differ only in the response's name in their terms
  px <- prior$data@x
  fx <- fitSampler$data@x
  attr(px, "terms") <- NULL
  attr(fx, "terms") <- NULL
  cutsEqual <- identical(prior$data@cutPoints, fitSampler$data@cutPoints) &&
    identical(px, fx)
  set.seed(seed)
  fitSampler$sampleTreesFromPrior()
  fitSampler$setLeafPrior(dbartsPriors$normal(k = 1.1))
  fitSampler$sampleLeafParametersFromPrior()
  fitSampler$setResponse(as.double(design$y))
  r <- fitSampler$run(0L, 1L)
  dTrain <- max(abs(
    as.numeric(r$train) - as.numeric(fitSampler$predict(trainX))
  ))
  dTest <- if (is.null(test)) {
    0
  } else {
    max(abs(as.numeric(r$test) - as.numeric(fitSampler$predict(test))))
  }
  # the prior sampler's predict at the fit's test rows equals the fit sampler's
  # at the same state: the truth is evaluated where the fit is read
  fitSampler$setLeafPrior(dbartsPriors$normal(k = 0.7))
  fitSampler$setLeafPrior(dbartsPriors$normal(dbartsPriors$chi(kDf, kScale)))
  kBack <- as.numeric(fitSampler$getK())
  r2 <- fitSampler$run(50L, 20L)
  kMoves <- stats::sd(as.numeric(r2$k)) > 0
  res <- c(
    cutsEqual = cutsEqual,
    fitsTrain = dTrain < 1e-8,
    fitsTest = dTest < 1e-8,
    kHandBack = isTRUE(all.equal(kBack, 0.7)),
    kMoves = kMoves
  )
  saveAtomic(
    list(
      cell = cell,
      link = link,
      rep = rep,
      res = res,
      dTrain = dTrain,
      dTest = dTest,
      quick = quick
    ),
    checkFile("self", link, paste0(cellFile(cell), ".rds"))
  )
  report(sprintf("self %s %s", cell, link), all(res), res)
}

# The fits the build checks compare: on pima (numeric) and creditapproval
# (factors, missing values), both links, a drawn-k fit at two lengths from one
# seed (the prefix property), a fixed-k fit; stan4bart on Twins, probit, at two
# lengths from one warmup and with k fixed. Saved whole for `check compare`.
checkBuildFits <- function(out) {
  res <- list(build = build(), lib = find.package("dbarts"))
  for (name in c("pima", "creditapproval")) {
    raw <- realDatasets[[name]]()
    rows <- withSeed(11L, sample.int(nrow(raw$x), 240L))
    train <- data.frame(
      raw$x[rows[1:200], , drop = FALSE],
      y = raw$y[rows[1:200]]
    )
    test <- raw$x[rows[201:240], , drop = FALSE]
    for (link in both) {
      run <- function(n, k = NULL) {
        args <- list(
          y ~ .,
          data = train,
          test = test,
          family = link,
          n.chains = 2L,
          n.threads = 1L,
          n.burn = 100L,
          n.samples = n,
          seed = 7L,
          verbose = FALSE,
          combineChains = FALSE
        )
        if (!is.null(k)) {
          args$k <- k
        }
        f <- do.call(bart, args)
        list(train = f$yhat.train, test = f$yhat.test, k = f$k)
      }
      res[[paste(name, link, "drawn100", sep = ".")]] <- run(100L)
      res[[paste(name, link, "drawn200", sep = ".")]] <- run(200L)
      res[[paste(name, link, "fixed", sep = ".")]] <- run(100L, 1.3)
    }
  }
  suppressPackageStartupMessages(library(stan4bart))
  d <- twinsData()
  keep <- which(stats::complete.cases(d$x))[1:300]
  x <- d$x[keep, , drop = FALSE]
  grp <- factor(x$brstate_reg)
  x <- x[,
    setdiff(
      names(x),
      c("brstate", "brstate_reg", "stoccfipb", "stoccfipb_reg")
    ),
    drop = FALSE
  ]
  y <- d$y0[keep]
  s4 <- function(iter, k = NULL) {
    args <- list(
      y ~ bart(. - grp) + (1 | grp),
      data = data.frame(x, grp = grp, y = y),
      verbose = -1L,
      chains = 1L,
      cores = 1L,
      iter = iter,
      warmup = 100L,
      seed = 5L
    )
    if (!is.null(k)) {
      args$bart_args <- list(k = k)
    }
    f <- do.call(stan4bart::stan4bart, args)
    list(
      fit = dbarts::extract(f, "indiv.bart"),
      ranef = dbarts::extract(f, "ranef")[[1L]]
    )
  }
  res$s4b200 <- s4(200L)
  res$s4b300 <- s4(300L)
  res$s4bFixed <- s4(200L, 1.3)
  saveRDS(res, out)
  cat("build fits saved:", out, "\n")
}

# The build gates on the three builds' fits: the move differs from the shipped
# build in drawn-k fits (both links); fixed-k fits equal on all three builds,
# bit for bit; the shipped build equal to main in drawn-k fits; the prefix
# property on every build (a longer run's first draws are the shorter run's).
checkCompare <- function(shippedFile, moveFile, mainFile) {
  f <- list(
    shipped = readRDS(shippedFile),
    move = readRDS(moveFile),
    main = readRDS(mainFile)
  )
  keys <- grep("[.]drawn100$", names(f$shipped), value = TRUE)
  rows <- list()
  for (key in keys) {
    base <- sub("[.]drawn100$", "", key)
    fixedKey <- paste0(base, ".fixed")
    rows[[base]] <- c(
      moveDiffers = !identical(f$shipped[[key]], f$move[[key]]),
      shippedIsMain = identical(f$shipped[[key]], f$main[[key]]),
      fixedEqual = identical(f$shipped[[fixedKey]], f$move[[fixedKey]]) &&
        identical(f$shipped[[fixedKey]], f$main[[fixedKey]]),
      prefixShipped = prefixHolds(
        f$shipped[[key]],
        f$shipped[[paste0(base, ".drawn200")]]
      ),
      prefixMove = prefixHolds(
        f$move[[key]],
        f$move[[paste0(base, ".drawn200")]]
      )
    )
  }
  s4 <- c(
    moveDiffers = !identical(f$shipped$s4b200, f$move$s4b200),
    shippedIsMain = identical(f$shipped$s4b200, f$main$s4b200),
    fixedEqual = identical(f$shipped$s4bFixed, f$move$s4bFixed) &&
      identical(f$shipped$s4bFixed, f$main$s4bFixed),
    prefixShipped = isTRUE(all(
      f$shipped$s4b300$fit[, 1:100] == f$shipped$s4b200$fit[, 1:100]
    )),
    prefixMove = isTRUE(all(
      f$move$s4b300$fit[, 1:100] == f$move$s4b200$fit[, 1:100]
    ))
  )
  rows$stan4bart <- s4
  tab <- do.call(rbind, rows)
  report("builds compare", all(tab), tab)
}

prefixHolds <- function(a, b) {
  n <- dim(a$train)[2L]
  isTRUE(all(a$train == b$train[, seq_len(n), , drop = FALSE])) &&
    isTRUE(all(a$test == b$test[, seq_len(n), , drop = FALSE])) &&
    isTRUE(all(a$k == b$k[, seq_len(n), drop = FALSE]))
}

# stan4bart's prior against the harness's transcription. stan4bart's probit
# takes only weights of 0 and 1 and refuses all zeros, so one row stays in the
# likelihood: stan4bart then samples the prior tilted by that row's
# likelihood, and the harness's prior draws are weighted by the same
# likelihood (a weight in [0, 1], so the reweighting is stable). Quantiles 5,
# 25, 50, 75 and 95 of the group standard deviation, the active row's group
# intercept, another group's intercept and the BART fit at another row must
# agree within 4.5 combined standard errors (stan4bart's Monte Carlo errors;
# the harness's by bootstrap).
weightedQuantile <- function(v, w, p) {
  o <- order(v)
  cw <- cumsum(w[o]) / sum(w)
  v[o][pmin(length(v), findInterval(p, cw) + 1L)]
}

checkS4bPrior <- function(nDraws = 2000L, nPrior = 20000L) {
  suppressPackageStartupMessages(library(stan4bart))
  d <- twinsData()
  keep <- which(stats::complete.cases(d$x))[1:300]
  x <- d$x[keep, , drop = FALSE]
  grp <- factor(x$brstate_reg)
  x <- x[,
    setdiff(
      names(x),
      c("brstate", "brstate_reg", "stoccfipb", "stoccfipb_reg")
    ),
    drop = FALSE
  ]
  y <- d$y0[keep]
  w <- c(1, rep(0, length(y) - 1L))
  g1 <- as.integer(grp[1L])
  g2 <- if (g1 == 1L) 2L else 1L
  chains <- 4L
  fit <- stan4bart::stan4bart(
    y ~ bart(. - grp) + (1 | grp),
    data = data.frame(x, grp = grp, y = y),
    weights = w,
    verbose = -1L,
    chains = chains,
    cores = 1L,
    iter = 500L + nDraws / chains,
    warmup = 500L,
    seed = 3L
  )
  sFit <- dbarts::extract(fit, "indiv.bart", combine_chains = FALSE)[2L, , ]
  sTau <- matrix(
    sqrt(dbarts::extract(fit, "Sigma", combine_chains = FALSE)[[1L]]),
    ncol = chains
  )
  re <- dbarts::extract(fit, "ranef", combine_chains = FALSE)[[1L]]
  sB1 <- matrix(re[1L, g1, , ], ncol = chains)
  sB2 <- matrix(re[1L, g2, , ], ncol = chains)
  sampler <- priorSampler(x, "probit", 3L)
  h <- withSeed(
    4L,
    t(vapply(
      seq_len(nPrior),
      function(i) {
        p <- priorDraw(sampler, list(x[1:2, , drop = FALSE]))
        tau <- rexp(1L)
        b <- rnorm(2L, 0, tau)
        c(
          f1 = p$f[[1L]][1L] + b[1L],
          fit = p$f[[1L]][2L],
          tau = tau,
          b1 = b[1L],
          b2 = b[2L]
        )
      },
      numeric(5L)
    ))
  )
  lik <- if (y[1L] == 1L) pnorm(h[, "f1"]) else pnorm(-h[, "f1"])
  probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  cmp <- function(sv, hv) {
    qs <- quantile(as.vector(sv), probs, names = FALSE)
    qh <- vapply(probs, function(p) weightedQuantile(hv, lik, p), 0)
    boot <- withSeed(
      5L,
      replicate(200L, {
        i <- sample.int(length(hv), replace = TRUE)
        vapply(probs, function(p) weightedQuantile(hv[i], lik[i], p), 0)
      })
    )
    se <- sqrt(
      vapply(probs, function(p) posterior::mcse_quantile(sv, p), 0)^2 +
        apply(boot, 1L, stats::sd)^2
    )
    data.frame(p = probs, stan4bart = qs, harness = qh, z = (qs - qh) / se)
  }
  tab <- rbind(
    cbind(what = "BART fit, another row", cmp(sFit, h[, "fit"])),
    cbind(what = "group sd", cmp(sTau, h[, "tau"])),
    cbind(what = "active row's group", cmp(sB1, h[, "b1"])),
    cbind(what = "another group", cmp(sB2, h[, "b2"]))
  )
  saveAtomic(
    list(table = tab, ess = sum(lik)^2 / sum(lik^2), quick = quick),
    checkFile("s4b-prior.rds")
  )
  report("stan4bart prior", all(abs(tab$z) <= 4.5), tab)
}

# The selection bands from the k study's simulated cells.
checkBands <- function() {
  b <- computeBands()
  saveAtomic(b, runDir("bands.rds"))
  report("bands", TRUE, list(bands = b$bands, cells = b$cells))
}

# Each selecting regime's cost per accepted replication against its fit cost:
# the design files' seconds against the shipped arm's (from fits when there,
# else from the timing pass: its shipped seconds a sweep times the arm's sweeps
# over all chains).
fitSecondsFromTiming <- function(cell, link) {
  files <- list.files(
    runDir("timing", link, cellFile(cell)),
    pattern = "-shipped-t[0-9]+[.]rds$",
    full.names = TRUE
  )
  if (length(files) == 0L) {
    return(NA_real_)
  }
  perRun <- median(vapply(files, function(f) readRDS(f)$seconds, 0))
  sweeps <- if (cellTask(cell) == "s4b") {
    mcmc$s4bWarmup + mcmc$s4bLong
  } else {
    mcmc$burn + drawsPerChain(cell)
  }
  timed <- if (cellTask(cell) %in% c("s4b", "loop")) {
    mcmc$timingWarm + mcmc$timingRun
  } else {
    mcmc$timingRun
  }
  perRun / timed * sweeps * cellChains(cell)
}

checkSelectionCost <- function() {
  rows <- list()
  for (reg in regimes) {
    if (reg$kind != "sbc") {
      next
    }
    for (link in reg$links) {
      dir <- runDir("design", link, cellFile(reg$id))
      files <- list.files(dir, pattern = "[.]rds$", full.names = TRUE)
      if (length(files) == 0L) {
        next
      }
      ds <- lapply(files, readRDS)
      fits <- list.files(
        runDir("fits", link, cellFile(reg$id)),
        pattern = "shipped[.]rds$",
        full.names = TRUE
      )
      fitSec <- if (length(fits)) {
        median(vapply(fits, function(f) readRDS(f)$seconds, 0))
      } else {
        fitSecondsFromTiming(reg$id, link)
      }
      designSec <- vapply(ds, function(d) d$seconds, 0)
      rows[[length(rows) + 1L]] <- data.frame(
        regime = reg$id,
        link = link,
        designs = length(ds),
        attempts = mean(vapply(ds, function(d) d$attempts, 0)),
        designSeconds = mean(designSec),
        designMax = max(designSec),
        fitSeconds = fitSec,
        within = is.na(fitSec) || mean(designSec) <= fitSec
      )
    }
  }
  tab <- do.call(rbind, rows)
  saveAtomic(list(table = tab, quick = quick), checkFile("selection-cost.rds"))
  report("selection cost", all(tab$within), tab)
}

# The pilot stored D1's shipped default reading for replications 1 to 10.
checkPilotStored <- function() {
  bisRun <- Sys.getenv("KMN_BIS_RUN")
  want <- expand.grid(
    link = both,
    cell = c("arr1-mod", "arr1-sep", "arr2-mod", "null-mod"),
    rep = 1:10,
    stringsAsFactors = FALSE
  )
  want$path <- file.path(
    bisRun,
    "draws",
    want$link,
    paste0("A6_", want$cell),
    sprintf("r%03d-shipped-default.f32", want$rep)
  )
  want$present <- file.exists(want$path)
  report(
    "pilot stored",
    all(want$present),
    table(present = want$present, link = want$link)
  )
}

# No acceptance rule, row label or subgroup reads the truth: the functions
# that make them never name it, and selection is called with (x, y) alone.
checkBlind <- function() {
  fns <- c(
    "selects",
    "separationScore",
    "interiorLabels",
    "subgroupsOf",
    "referenceLabels",
    "relabel",
    "classesFrom",
    "acceptReference",
    "mostDistinct",
    "quartiles",
    "summarizeClass",
    "massReading",
    "halfClasses",
    "referenceClasses",
    "runRefscore",
    "runFixedk",
    "fitBart",
    "fitLoop",
    "fitStan4bart",
    "fitPropensityS9",
    "fitData"
  )
  bad <- vapply(
    fns,
    function(f) {
      txt <- paste(deparse(get(f)), collapse = "\n")
      grepl("truth|\\bk0\\b|\\be0\\b|\\bf0\\b|designed", txt)
    },
    logical(1L)
  )
  call <- paste(deparse(sbcDesign), collapse = "\n")
  selectionCall <- grepl(
    "selects(reg, trainX[seen, , drop = FALSE], y[seen]",
    call,
    fixed = TRUE
  )
  runArmText <- paste(deparse(runArm), collapse = "\n")
  # runArm reads the truth only for the oracle's k
  armUse <- !grepl("truth\\$f|truth\\$b0|truth\\$e0|truth\\$tau0", runArmText)
  report(
    "truth-blind labels",
    !any(bad) && selectionCall && armUse,
    list(
      reading = fns[bad],
      selectionCall = selectionCall,
      oracleOnlyK = armUse
    )
  )
}

# Every kernel's readings scored against the other kernel's half, a tail
# never against its own chains, the oracle against the pooled reference.
checkHalves <- function() {
  files <- list.files(
    runDir("refscore"),
    pattern = "[.]rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  bad <- character()
  for (f in files) {
    rs <- readRDS(f)
    for (key in names(rs$readings)) {
      arm <- sub("[.].*$", "", key)
      want <- switch(
        arm,
        shipped = "move",
        move = "shipped",
        mixed = "move",
        oracle = if (rs$accepted) "pooled" else "move"
      )
      if (!identical(rs$readings[[key]]$against, want)) {
        bad <- c(bad, paste(basename(f), key))
      }
    }
  }
  report(
    "scoring halves",
    length(files) > 0L && length(bad) == 0L,
    list(files = length(files), wrong = bad)
  )
}

# The binary interval study's defects, not repeated, on a job file: every
# replication of every phase has its jobs, propensities in the designs; every
# dependency names an earlier job; no two jobs write one output; fixed-k and
# reference scoring wait for the final acceptance.
checkDefects <- function(jobfile, phase = "tranche1", tranche = 1L) {
  jobs <- readJobs(jobfile)
  ids <- jobs$id
  pos <- setNames(seq_along(ids), ids)
  depBad <- ids[
    !vapply(
      seq_along(ids),
      function(i) {
        d <- jobDeps(jobs$deps[i])
        all(d %in% ids) && all(pos[d] < i)
      },
      logical(1L)
    )
  ]
  depOk <- length(depBad) == 0L && !anyDuplicated(ids)
  outs <- unlist(lapply(jobs$cmds, function(c) {
    vapply(strsplit(c, " ;; ", fixed = TRUE)[[1L]], commandOutput, "")
  }))
  outs <- outs[nzchar(outs)]
  unique <- !anyDuplicated(outs)
  waits <- all(vapply(
    which(grepl("^(fixedk|refscore)-", ids)),
    function(i) {
      any(grepl("^acceptfinal-", jobDeps(jobs$deps[i])))
    },
    logical(1L)
  ))
  missing <- character()
  for (reg in regimes) {
    for (cell in regimeCells(reg)) {
      for (link in reg$links) {
        for (rep in seq_len(
          if (phase == "smoke") 1L else repsFor(cell, tranche)
        )) {
          need <- c(
            jobId("design", cell, link, rep),
            jobId("fit", cell, link, rep, "shipped")
          )
          if (reg$kind == "sbc") {
            need <- c(
              need,
              jobId("fit", cell, link, rep, "oracle"),
              jobId("score", cell, link, rep)
            )
          }
          if (isRefRep(cell, rep, tranche)) {
            need <- c(
              need,
              jobId("fit", cell, link, rep, "move"),
              jobId("acceptfinal", cell, link, rep),
              jobId("fixedk", cell, link, rep),
              jobId("refscore", cell, link, rep)
            )
          }
          if (cellRegime(cell)$id == "B" && rep > 1L) {
            need <- setdiff(
              need,
              c(
                jobId("acceptfinal", cell, link, rep),
                jobId("fixedk", cell, link, rep),
                jobId("refscore", cell, link, rep)
              )
            )
          }
          missing <- c(missing, setdiff(need, ids))
        }
      }
    }
  }
  report(
    "defects",
    depOk && unique && waits && length(missing) == 0L,
    list(
      jobs = length(ids),
      dependenciesEarlier = depOk,
      badDependencies = head(depBad, 10L),
      outputsUnique = unique,
      waitForFinalAcceptance = waits,
      missing = head(missing, 20L)
    )
  )
}

# A command's output path, for the uniqueness check.
commandOutput <- function(cmd) {
  b <- sub("^@([a-z]+) .*$", "\\1", cmd)
  a <- strsplit(sub("^@[a-z]+ ", "", cmd), " ", fixed = TRUE)[[1L]]
  key <- function(v) paste(v[!is.na(v)], collapse = " ")
  switch(
    a[1L],
    design = key(a[1:5]),
    fit = key(a[1:6]),
    score = key(a[1:5]),
    reference = switch(
      a[2L],
      run = key(c(a[1:6], b)),
      accept = key(a[1:6]),
      ""
    ),
    fixedk = key(a[1:4]),
    refscore = key(a[1:4]),
    timing = key(c(a[1:5], b)),
    ""
  )
}

# The SBC null run: on S6, the oracle's and the long oracle's pooled coverage
# within noise of nominal in every class, both links: |coverage - 0.90| at most
# 2.81 bootstrap standard errors (two-sided 0.005, about 0.05 over the run's
# ten or so readings).
checkNull <- function(reps = 100L) {
  tab <- sbcTable("S6", reps, arms = c("oracle.default", "longOracle.default"))
  tab$z <- (tab$cov90 - 0.9) / tab$se90
  saveAtomic(list(table = tab, quick = quick), checkFile("null.rds"))
  report(
    "SBC null run",
    nrow(tab) > 0L && all(abs(tab$z) <= 2.81, na.rm = TRUE),
    tab
  )
}

# The discrimination runs must flag: the outcome drawn through the other link,
# read on the oracle; the truth from a 50-tree prior fit at 75 trees, read on
# the shipped arm. Flagged: the 95 percent bound of interior f's pooled
# coverage excludes 0.90.
checkDiscrimination <- function(reps = 100L) {
  out <- list()
  for (disc in c("link", "trees50")) {
    arm <- if (disc == "link") "oracle.default" else "shipped.default"
    tab <- sbcTable("S6", reps, arms = arm, discriminate = disc)
    tab <- tab[tab$class == "fInt", , drop = FALSE]
    tab$flag <- abs(tab$cov90 - 0.9) - zNoProblem * tab$se90 > 0
    tab$disc <- disc
    out[[disc]] <- tab
  }
  tab <- do.call(rbind, out)
  saveAtomic(list(table = tab, quick = quick), checkFile("discrimination.rds"))
  report("discrimination", nrow(tab) > 0L && all(tab$flag), tab)
}

## ------------------------------------------------------------------- read

# Pooled coverage by class over a regime's scored replications, both links,
# with bootstrap standard errors over replications (the same resample for
# every link and class of the regime).
sbcScores <- function(regimeId, reps, discriminate = "") {
  out <- list()
  for (link in regimes[[regimeId]]$links) {
    for (rep in seq_len(reps)) {
      f <- scoreFile(regimeId, link, rep, discriminate)
      if (file.exists(f)) {
        out[[paste(link, rep)]] <- loadChecked(f)
      }
    }
  }
  out
}

sbcTable <- function(regimeId, reps, arms = NULL, discriminate = "") {
  sc <- sbcScores(regimeId, reps, discriminate)
  if (length(sc) == 0L) {
    return(data.frame())
  }
  repIds <- sort(unique(vapply(sc, function(s) s$rep, 0L)))
  idx <- withSeed(
    seedOf("bootstrap", regimeId, discriminate),
    lapply(seq_len(bootstrapResamples), function(b) {
      sample.int(length(repIds), length(repIds), replace = TRUE)
    })
  )
  rows <- list()
  for (link in unique(vapply(sc, `[[`, "", "link"))) {
    mine <- Filter(function(s) s$link == link, sc)
    byRep <- mine[match(repIds, vapply(mine, function(s) s$rep, 0L))]
    armNames <- if (is.null(arms)) {
      unique(unlist(lapply(mine, function(s) names(s$coverage))))
    } else {
      arms
    }
    for (a in armNames) {
      classes <- unique(unlist(lapply(mine, function(s) {
        names(s$coverage[[a]])
      })))
      for (cl in classes) {
        get1 <- function(field) {
          vapply(
            byRep,
            function(s) {
              if (is.null(s) || is.null(s$coverage[[a]][[cl]])) {
                NA_real_
              } else {
                s$coverage[[a]][[cl]][[field]]
              }
            },
            0
          )
        }
        n <- get1("n")
        c90 <- get1("cov90")
        c95 <- get1("cov95")
        ok <- is.finite(n) & n > 0
        if (!any(ok)) {
          next
        }
        pooled <- function(i, v) sum(v[i][ok[i]]) / sum(n[i][ok[i]])
        b90 <- vapply(idx, function(i) pooled(i, c90), 0)
        b95 <- vapply(idx, function(i) pooled(i, c95), 0)
        floorSe <- sqrt(0.09 / sum(n[ok]))
        rows[[length(rows) + 1L]] <- data.frame(
          regime = regimeId,
          link = link,
          arm = a,
          class = cl,
          reps = sum(ok),
          n = sum(n[ok]),
          cov90 = pooled(seq_along(n), c90),
          se90 = max(sd(b90, na.rm = TRUE), floorSe),
          cov95 = pooled(seq_along(n), c95),
          se95 = max(sd(b95, na.rm = TRUE), sqrt(0.0475 / sum(n[ok])))
        )
      }
    }
  }
  do.call(rbind, rows)
}

# G and D of a regime's SBC cells, with the reading of each.
sbcReading <- function(regimeId, reps) {
  sc <- sbcScores(regimeId, reps)
  if (length(sc) == 0L) {
    return(data.frame())
  }
  repIds <- sort(unique(vapply(sc, function(s) s$rep, 0L)))
  idx <- withSeed(
    seedOf("bootstrap", regimeId),
    lapply(seq_len(bootstrapResamples), function(b) {
      sample.int(length(repIds), length(repIds), replace = TRUE)
    })
  )
  rows <- list()
  for (link in unique(vapply(sc, `[[`, "", "link"))) {
    mine <- Filter(function(s) s$link == link, sc)
    byRep <- mine[match(repIds, vapply(mine, function(s) s$rep, 0L))]
    classes <- unique(unlist(lapply(mine, function(s) {
      names(s$coverage$shipped.default)
    })))
    for (cl in classes) {
      get1 <- function(a, field) {
        vapply(
          byRep,
          function(s) {
            if (is.null(s) || is.null(s$coverage[[a]][[cl]])) {
              NA_real_
            } else {
              s$coverage[[a]][[cl]][[field]]
            }
          },
          0
        )
      }
      n <- get1("shipped.default", "n")
      cs <- get1("shipped.default", "cov90")
      co <- get1("oracle.default", "cov90")
      ok <- is.finite(n) & n > 0 & is.finite(cs) & is.finite(co)
      if (!any(ok)) {
        next
      }
      pooled <- function(i, v) sum(v[i][ok[i]]) / sum(n[i][ok[i]])
      gEst <- pooled(seq_along(n), cs) - 0.9
      dEst <- pooled(seq_along(n), cs) - pooled(seq_along(n), co)
      bG <- vapply(idx, function(i) pooled(i, cs) - 0.9, 0)
      bD <- vapply(idx, function(i) pooled(i, cs) - pooled(i, co), 0)
      floorSe <- sqrt(0.09 / sum(n[ok]))
      seG <- max(sd(bG, na.rm = TRUE), floorSe)
      seD <- max(sd(bD, na.rm = TRUE), floorSe)
      G <- correctedDistanceFloor(abs(gEst), seG)
      band <- isTRUE(classInfo[[cl]]$band)
      problem <- if (band) {
        abs(dEst) >= bar && abs(dEst) - zProblem * seD > 0
      } else {
        G >= bar &&
          abs(gEst) - zProblem * seG > 0 &&
          abs(dEst) >= bar &&
          sign(dEst) == sign(gEst) &&
          abs(dEst) - zProblem * seD > 0
      }
      measured <- (!band && abs(gEst) + zNoProblem * seG < bar) ||
        (dEst - zNoProblem * seD > -bar && dEst + zNoProblem * seD < bar)
      rows[[length(rows) + 1L]] <- data.frame(
        regime = regimeId,
        cell = regimeId,
        link = link,
        class = cl,
        reps = sum(ok),
        coverage = gEst + 0.9,
        gSigned = gEst,
        G = G,
        seG = seG,
        D = dEst,
        seD = seD,
        sbcReading = if (problem) {
          "problem"
        } else if (measured) {
          "measured"
        } else {
          "unresolved"
        }
      )
    }
  }
  do.call(rbind, rows)
}

# dec-B455's noise-corrected distance with a zero spread never read as no
# noise: the standard error passed in already carries the binomial floor.
correctedDistanceFloor <- function(observed, s) {
  d <- observed
  floorAt <- s * sqrt(2 / pi)
  if (s <= 0 || observed <= floorAt) {
    return(0)
  }
  for (iter in seq_len(60L)) {
    g <- s * sqrt(2 / pi) * exp(-d^2 / (2 * s^2)) + d * (1 - 2 * pnorm(-d / s))
    slope <- 1 - 2 * pnorm(-d / s)
    step <- (g - observed) / max(slope, 1e-12)
    d <- max(d - step, 0)
    if (abs(step) < 1e-12) break
  }
  d
}

# The reference readings of a cell: per class, S (0.90 minus the shipped
# default's mass against the move half), the stuck-run shortfall, the
# transients and k's shares, over the cell's passing references (B: over its
# six seeds against its one reference), with bootstrap standard errors and t
# bounds; the 80 percent rule and its sensitivity reading.
refRows <- function(cell, link, tranche) {
  reps <- seq_len(refsFor(cell, tranche))
  rs <- list()
  fk <- list()
  for (rep in reps) {
    f <- refscoreFile(cell, link, rep)
    g <- fixedkFile(cell, link, rep)
    if (file.exists(f)) {
      rs[[as.character(rep)]] <- loadChecked(f)
    }
    if (file.exists(g)) {
      fk[[as.character(rep)]] <- loadChecked(g)
    }
  }
  if (length(rs) == 0L) {
    return(NULL)
  }
  isB <- cellRegime(cell)$id == "B"
  mass <- function(r, key, cl, half = "vsOther") {
    x <- r$readings[[key]][[half]][[cl]]
    if (is.null(x)) NA_real_ else x[["mass90"]]
  }
  classes <- unique(unlist(lapply(rs, function(r) {
    names(r$readings$shipped.tails$vsOther)
  })))
  rows <- list()
  for (cl in classes) {
    per <- do.call(
      rbind,
      lapply(names(rs), function(rep) {
        r <- rs[[rep]]
        # B's six seeds are its replications, all read against its one reference
        keysOf <- function(prefix) {
          if (isB) {
            grep(
              paste0("^", prefix, "[.]default[.]r[0-9]+$"),
              names(r$readings),
              value = TRUE
            )
          } else {
            paste0(prefix, ".default")
          }
        }
        sh <- vapply(keysOf("shipped"), function(k) mass(r, k, cl), 0)
        mv <- vapply(keysOf("move"), function(k) mass(r, k, cl), 0)
        shMove <- vapply(
          keysOf("shipped"),
          function(k) mass(r, k, cl, "vsMove"),
          0
        )
        st <- if (!is.null(fk[[rep]]$stuck[[cl]])) {
          fk[[rep]]$stuck[[cl]][["mass90"]]
        } else {
          NA_real_
        }
        data.frame(
          rep = as.integer(rep),
          seed = seq_along(sh),
          accepted = r$accepted,
          S = 0.9 - sh,
          Smove = 0.9 - shMove,
          moveS = 0.9 - mv,
          stuck = 0.9 - st,
          shippedTransient = mass(r, "shipped.tails", cl) - sh,
          moveTransient = mass(r, "move.tails", cl) - mv,
          stationary = mass(r, "move.tails", cl) - mass(r, "shipped.tails", cl),
          kFromMoveHalf = isTRUE(fk[[rep]]$kFromMoveHalf)
        )
      })
    )
    per$kShare <- per$shippedTransient - per$moveTransient
    rows[[cl]] <- refClassReading(per, cell, link, cl, length(reps), isB)
  }
  do.call(rbind, rows)
}

refClassReading <- function(per, cell, link, cl, attempted, isB) {
  units <- if (isB) seq_len(nrow(per)) else seq_len(nrow(per))
  bound <- function(v, use) {
    v <- v[use]
    v <- v[is.finite(v)]
    if (length(v) < 2L) {
      return(c(
        est = if (length(v)) mean(v) else NA_real_,
        se = NA_real_,
        n = length(v)
      ))
    }
    b <- withSeed(
      seedOf("refboot", cell, link, cl),
      vapply(
        seq_len(bootstrapResamples),
        function(i) mean(sample(v, replace = TRUE)),
        0
      )
    )
    c(est = mean(v), se = max(sd(b), 0.001 / sqrt(length(v))), n = length(v))
  }
  passing <- per$accepted
  passShare <- if (isB) as.numeric(all(passing)) else mean(passing)
  decide <- function(use, sCol) {
    S <- bound(per[[sCol]], use)
    stuck <- bound(per$stuck, use & !per$kFromMoveHalf)
    k <- bound(per$kShare, use)
    df <- max(1, S[["n"]] - 1)
    tP <- qt(1 - 0.000156, df)
    tN <- qt(0.975, df)
    moveTest <- is.finite(S[["est"]]) &&
      mean(per$moveS[use], na.rm = TRUE) < bar &&
      S[["est"]] >= bar
    problem <- is.finite(S[["se"]]) &&
      S[["est"]] >= bar &&
      S[["est"]] - tP * S[["se"]] > 0 &&
      (moveTest ||
        (is.finite(stuck[["se"]]) &&
          stuck[["est"]] >= bar &&
          stuck[["est"]] - tP * stuck[["se"]] > 0))
    bounded <- is.finite(stuck[["se"]]) &&
      is.finite(k[["se"]]) &&
      stuck[["est"]] + tN * stuck[["se"]] < bar &&
      k[["est"]] + tN * k[["se"]] < bar
    measured <- is.finite(S[["se"]]) && S[["est"]] + tN * S[["se"]] < bar
    list(
      reading = if (problem) {
        "problem"
      } else if (bounded) {
        "bounded"
      } else if (measured) {
        "measured"
      } else {
        "unresolved"
      },
      S = S,
      stuck = stuck,
      k = k
    )
  }
  main <- decide(passing, "S")
  sens <- decide(rep(TRUE, nrow(per)), "Smove")
  reading <- main$reading
  reason <- ""
  if (passShare < 0.8) {
    reading <- "unresolved"
    reason <- "failed references"
  } else if (!identical(main$reading, sens$reading)) {
    reading <- "unresolved"
    reason <- "changes when failing references are read against their move half"
  } else if (reading == "unresolved") {
    reason <- "counts"
  }
  data.frame(
    cell = cell,
    link = link,
    class = cl,
    refs = length(unique(per$rep)),
    attempted = attempted,
    passShare = passShare,
    S = main$S[["est"]],
    seS = main$S[["se"]],
    stuck = main$stuck[["est"]],
    seStuck = main$stuck[["se"]],
    kShare = main$k[["est"]],
    seK = main$k[["se"]],
    shippedTransient = mean(per$shippedTransient[passing], na.rm = TRUE),
    moveTransient = mean(per$moveTransient[passing], na.rm = TRUE),
    stationary = mean(per$stationary[passing], na.rm = TRUE),
    moveS = mean(per$moveS[passing], na.rm = TRUE),
    refReading = reading,
    refReason = reason
  )
}

# A cell's tier: the rarer of its regime's and its estimand's; a B problem
# takes the tier of the regime its separation score and outcome rate place it
# in (effects: S7's, or S8's at a rare outcome).
cellTier <- function(cell, link, cl) {
  reg <- cellRegime(cell)
  t <- reg$tier
  if (reg$id == "B") {
    d <- loadDesign(cell, link, 1L)
    rate <- mean(d$y)
    rare <- rate <= 0.05 || rate >= 0.95
    if (d$task == "effect") {
      t <- if (rare) 2L else 1L
    } else {
      score <- separationScore(d$x, d$y, paste("tier", cell))
      b <- bands()
      t <- if (score >= b$nearsep[1L]) {
        3L
      } else if (rare || score > b$moderate[2L]) {
        2L
      } else {
        1L
      }
    }
  }
  max(t, classInfo[[cl]]$tier)
}

combineReading <- function(sbc, ref) {
  if (identical(sbc, "problem") || identical(ref, "problem")) {
    return("problem")
  }
  if (identical(ref, "bounded")) {
    return("no problem, bounded")
  }
  avail <- c(sbc, ref)
  avail <- avail[!is.na(avail)]
  if (length(avail) && all(avail == "measured")) {
    return("no problem, measured")
  }
  "unresolved"
}

runRead <- function(tranche = 1L) {
  all <- list()
  for (reg in regimes) {
    sbc <- if (reg$kind == "sbc") {
      sbcReading(reg$id, reg$reps[tranche])
    } else {
      data.frame()
    }
    for (cell in regimeCells(reg)) {
      for (link in reg$links) {
        ref <- refRows(cell, link, tranche)
        s <- if (nrow(sbc)) {
          sbc[sbc$link == link, , drop = FALSE]
        } else {
          data.frame()
        }
        classes <- union(
          if (nrow(s)) s$class else character(),
          if (!is.null(ref)) ref$class else character()
        )
        for (cl in classes) {
          sr <- if (nrow(s) && cl %in% s$class) {
            s[s$class == cl, , drop = FALSE]
          } else {
            NULL
          }
          rr <- if (!is.null(ref) && cl %in% ref$class) {
            ref[ref$class == cl, , drop = FALSE]
          } else {
            NULL
          }
          all[[length(all) + 1L]] <- data.frame(
            regime = reg$id,
            cell = cell,
            link = link,
            class = cl,
            tier = tryCatch(cellTier(cell, link, cl), error = function(e) {
              NA_integer_
            }),
            G = if (is.null(sr)) NA else sr$G,
            gSigned = if (is.null(sr)) NA else sr$gSigned,
            seG = if (is.null(sr)) NA else sr$seG,
            D = if (is.null(sr)) NA else sr$D,
            seD = if (is.null(sr)) NA else sr$seD,
            sbcReading = if (is.null(sr)) NA else sr$sbcReading,
            S = if (is.null(rr)) NA else rr$S,
            seS = if (is.null(rr)) NA else rr$seS,
            stuck = if (is.null(rr)) NA else rr$stuck,
            seStuck = if (is.null(rr)) NA else rr$seStuck,
            kShare = if (is.null(rr)) NA else rr$kShare,
            seK = if (is.null(rr)) NA else rr$seK,
            refs = if (is.null(rr)) NA else rr$refs,
            passShare = if (is.null(rr)) NA else rr$passShare,
            refReading = if (is.null(rr)) NA else rr$refReading,
            refReason = if (is.null(rr)) NA else rr$refReason,
            reading = combineReading(
              if (is.null(sr)) NA else sr$sbcReading,
              if (is.null(rr)) NA else rr$refReading
            )
          )
        }
      }
    }
  }
  tab <- do.call(rbind, all)
  saveAtomic(
    list(table = tab, tranche = tranche, quick = quick),
    runDir(sprintf("readings-t%d.rds", tranche))
  )
  num <- vapply(tab, is.numeric, logical(1L))
  shown <- tab
  shown[num] <- lapply(tab[num], function(v) round(v, 4L))
  print(shown, row.names = FALSE)
  invisible(tab)
}

## -------------------------------------------------------------------- jobs

# Job files for the box queue: one job a line, tab-separated: id, the ids it
# waits for (comma-separated, or "-"), big (0 or 1), cell key, expected
# seconds or NA, commands joined by " ;; ", each "@<build> <harness args>".
jobId <- function(...) {
  parts <- as.character(c(...))
  paste(gsub("[:]", "_", parts[nzchar(parts)]), collapse = "-")
}

jobDeps <- function(field) {
  if (identical(field, "-")) {
    character()
  } else {
    strsplit(field, ",", fixed = TRUE)[[1L]]
  }
}

readJobs <- function(path) {
  lines <- readLines(path)
  parts <- strsplit(lines, "\t", fixed = TRUE)
  data.frame(
    id = vapply(parts, `[`, "", 1L),
    deps = vapply(parts, `[`, "", 2L),
    big = vapply(parts, `[`, "", 3L),
    cell = vapply(parts, `[`, "", 4L),
    est = vapply(parts, `[`, "", 5L),
    cmds = vapply(parts, `[`, "", 6L)
  )
}

longOracleRegimes <- c("S1", "S5", "S6", "S7", "E1")
longOracleReps <- 25L

# The jobs of one replication of one cell and link, in dependency order.
replicationJobs <- function(
  cell,
  link,
  rep,
  tranche,
  add,
  phases = c("fits", "references"),
  discriminate = "",
  longReps = longOracleReps,
  bSeeds = repsFor(cell, tranche)
) {
  reg <- cellRegime(cell)
  big <- reg$big
  disc <- if (nzchar(discriminate)) paste0(" ", discriminate) else ""
  design <- jobId("design", cell, link, rep, discriminate)
  add(
    design,
    character(),
    big,
    cell,
    sprintf("@shipped design %s %s %d%s", cell, link, rep, disc)
  )
  ref <- isRefRep(cell, rep, tranche) && !nzchar(discriminate)
  if (reg$id == "B" && rep > 1L) {
    # B's seeds share its one reference, made on its first replication
    ref <- FALSE
  }
  fitIds <- character()
  if ("fits" %in% phases) {
    arms <- "shipped"
    if (reg$kind == "sbc") {
      arms <- c(arms, "oracle")
    }
    if (ref || (reg$id == "B")) {
      arms <- c(arms, "move")
    }
    if (
      reg$kind == "sbc" &&
        reg$id %in% longOracleRegimes &&
        rep <= longReps &&
        !nzchar(discriminate)
    ) {
      arms <- c(arms, "longOracle")
    }
    for (arm in arms) {
      id <- jobId("fit", cell, link, rep, arm, discriminate)
      add(
        id,
        design,
        big,
        cell,
        sprintf(
          "@%s fit %s %s %d %s%s",
          armBuild[[arm]],
          cell,
          link,
          rep,
          arm,
          disc
        )
      )
      fitIds <- c(fitIds, id)
    }
    if (reg$id == "D2") {
      id <- jobId("fit", cell, link, rep, "mixed")
      add(
        id,
        jobId("fit", cell, link, rep, "move"),
        big,
        cell,
        sprintf("@shipped fit %s %s %d mixed", cell, link, rep)
      )
      fitIds <- c(fitIds, id)
    }
    if (reg$kind == "sbc") {
      add(
        jobId("score", cell, link, rep, discriminate),
        fitIds,
        FALSE,
        cell,
        sprintf("@shipped score %s %s %d%s", cell, link, rep, disc)
      )
    }
  }
  if (ref && "references" %in% phases) {
    parts <- seq_len(refPartsPerKernel(cell))
    partIds <- character()
    for (which in c("shipped", "move")) {
      for (part in parts) {
        id <- jobId("ref", cell, link, rep, which, part)
        add(
          id,
          design,
          big,
          cell,
          sprintf("@%s reference run %s %s %d %d", which, cell, link, rep, part)
        )
        partIds <- c(partIds, id)
      }
    }
    acc <- jobId("accept", cell, link, rep)
    add(
      acc,
      partIds,
      big,
      cell,
      sprintf("@shipped reference accept %s %s %d", cell, link, rep)
    )
    escIds <- character()
    for (part in parts) {
      id <- jobId("escalate", cell, link, rep, part)
      add(
        id,
        acc,
        big,
        cell,
        sprintf(
          "@shipped reference escalate %s %s %d %d",
          cell,
          link,
          rep,
          part
        )
      )
      escIds <- c(escIds, id)
    }
    fin <- jobId("acceptfinal", cell, link, rep)
    add(
      fin,
      escIds,
      big,
      cell,
      sprintf("@shipped reference accept %s %s %d final", cell, link, rep)
    )
    add(
      jobId("fixedk", cell, link, rep),
      fin,
      big,
      cell,
      sprintf("@shipped fixedk %s %s %d", cell, link, rep)
    )
    seedFits <- if (reg$id == "B") {
      unlist(lapply(seq_len(bSeeds), function(r) {
        c(
          jobId("fit", cell, link, r, "shipped"),
          jobId("fit", cell, link, r, "move")
        )
      }))
    } else {
      fitIds
    }
    add(
      jobId("refscore", cell, link, rep),
      c(fin, jobId("fixedk", cell, link, rep), seedFits),
      big,
      cell,
      sprintf(
        "@shipped refscore %s %s %d%s",
        cell,
        link,
        rep,
        if (reg$id == "B") paste0(" ", bSeeds) else ""
      ),
      last = reg$id == "B"
    )
  }
}

writeJobs <- function(phase, path, tranche = 1L) {
  jobs <- list()
  tail <- list()
  add <- function(id, deps, big, cell, cmd, last = FALSE) {
    line <- paste(
      id,
      if (length(deps)) paste(deps, collapse = ",") else "-",
      as.integer(big),
      cell,
      "NA",
      cmd,
      sep = "\t"
    )
    if (last) {
      tail[[length(tail) + 1L]] <<- line
    } else {
      jobs[[length(jobs) + 1L]] <<- line
    }
  }
  allCells <- unlist(lapply(regimes, regimeCells))
  if (phase == "smoke") {
    for (cell in allCells) {
      for (link in cellRegime(cell)$links) {
        replicationJobs(cell, link, 1L, tranche, add, bSeeds = 1L)
      }
    }
  } else if (phase == "tranche1") {
    order <- c(
      "S6",
      "S1",
      withSeed(seedOf("cell-order"), sample(setdiff(allCells, c("S6", "S1"))))
    )
    for (cell in order) {
      for (link in cellRegime(cell)$links) {
        for (rep in seq_len(repsFor(cell, tranche))) {
          replicationJobs(cell, link, rep, tranche, add)
        }
      }
    }
  } else if (phase == "null") {
    for (link in both) {
      for (rep in seq_len(100L)) {
        replicationJobs(
          "S6",
          link,
          rep,
          tranche,
          add,
          phases = "fits",
          longReps = 100L
        )
      }
    }
  } else if (phase == "disc") {
    for (disc in c("link", "trees50")) {
      for (link in both) {
        for (rep in seq_len(100L)) {
          replicationJobs(
            "S6",
            link,
            rep,
            tranche,
            add,
            phases = "fits",
            discriminate = disc
          )
        }
      }
    }
  } else if (phase == "timing") {
    cells <- unlist(lapply(
      Filter(function(r) r$kind == "sbc", regimes),
      regimeCells
    ))
    for (cell in cells) {
      for (link in cellRegime(cell)$links) {
        for (rep in 1:3) {
          design <- jobId("design", cell, link, rep)
          add(
            design,
            character(),
            cellRegime(cell)$big,
            cell,
            sprintf("@shipped design %s %s %d", cell, link, rep)
          )
        }
      }
    }
    for (t in 1:3) {
      for (cell in cells) {
        for (link in cellRegime(cell)$links) {
          for (rep in 1:3) {
            prev <- jobId("design", cell, link, rep)
            for (which in if (t %% 2L == 1L) {
              c("shipped", "move")
            } else {
              c("move", "shipped")
            }) {
              id <- jobId("timing", cell, link, rep, which, t)
              add(
                id,
                prev,
                FALSE,
                cell,
                sprintf("@%s timing %s %s %d %d", which, cell, link, rep, t)
              )
              prev <- id
            }
          }
        }
      }
    }
  } else if (grepl("^designs:", phase)) {
    # designs of one regime, for the selection cost: designs:<regime>:<reps>
    p <- strsplit(phase, ":", fixed = TRUE)[[1L]]
    for (cell in regimeCells(regimes[[p[2L]]])) {
      for (link in cellRegime(cell)$links) {
        for (rep in seq_len(as.integer(p[3L]))) {
          add(
            jobId("design", cell, link, rep),
            character(),
            cellRegime(cell)$big,
            cell,
            sprintf("@shipped design %s %s %d", cell, link, rep)
          )
        }
      }
    }
  } else {
    stop("no phase ", phase)
  }
  lines <- c(unlist(jobs), unlist(tail))
  writeLines(lines, path)
  cat(sprintf("%d jobs written to %s\n", length(lines), path))
}

## -------------------------------------------------------------------- main

main <- function(args) {
  if (length(args) == 0L) {
    stop("no action")
  }
  action <- args[1L]
  a <- args[-1L]
  opt <- function(i) if (length(a) >= i) a[i] else ""
  switch(
    action,
    cells = {
      for (reg in regimes) {
        cat(reg$id, paste(regimeCells(reg), collapse = " "), "\n")
      }
    },
    jobs = writeJobs(a[1L], a[2L]),
    design = makeDesign(a[1L], a[2L], as.integer(a[3L]), opt(4L)),
    fit = for (arm in strsplit(a[4L], ",", fixed = TRUE)[[1L]]) {
      runArm(a[1L], a[2L], as.integer(a[3L]), arm, opt(5L))
    },
    score = runScore(a[1L], a[2L], as.integer(a[3L]), opt(4L)),
    reference = switch(
      a[1L],
      run = runReferencePart(
        a[2L],
        a[3L],
        as.integer(a[4L]),
        as.integer(a[5L])
      ),
      accept = acceptReference(
        a[2L],
        a[3L],
        as.integer(a[4L]),
        identical(opt(5L), "final")
      ),
      escalate = escalateReferencePart(
        a[2L],
        a[3L],
        as.integer(a[4L]),
        as.integer(a[5L])
      ),
      stop("no reference action ", a[1L])
    ),
    fixedk = runFixedk(a[1L], a[2L], as.integer(a[3L])),
    refscore = runRefscore(
      a[1L],
      a[2L],
      as.integer(a[3L]),
      if (nzchar(opt(4L))) as.integer(opt(4L)) else NA_integer_
    ),
    timing = if (identical(a[1L], "summarize")) {
      summarizeTiming()
    } else {
      runTiming(a[1L], a[2L], as.integer(a[3L]), as.integer(a[4L]))
    },
    screen = runScreen(a[1L]),
    read = runRead(if (length(a)) as.integer(a[1L]) else 1L),
    check = {
      ok <- switch(
        a[1L],
        bands = checkBands(),
        self = checkSelf(a[2L], a[3L], as.integer(a[4L])),
        fits = checkBuildFits(a[2L]),
        compare = checkCompare(a[2L], a[3L], a[4L]),
        s4bprior = checkS4bPrior(),
        selection = checkSelectionCost(),
        pilot = checkPilotStored(),
        blind = checkBlind(),
        halves = checkHalves(),
        defects = checkDefects(
          a[2L],
          if (nzchar(opt(3L))) opt(3L) else "tranche1"
        ),
        null = checkNull(),
        disc = checkDiscrimination(),
        stop("no check ", a[1L])
      )
      if (!isTRUE(ok)) {
        quit(status = 3L)
      }
    },
    stop("no action ", action)
  )
  invisible(NULL)
}

if (sys.nframe() == 0L || identical(Sys.getenv("KMN_MAIN"), "1")) {
  main(commandArgs(TRUE))
}
