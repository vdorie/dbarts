#!/usr/bin/env Rscript

# Pricing the collapsed scale move under logistic: steps 1 to 3, the timing,
# the real-data guard and the convergence measure of
# docs/plans/k-move-pricing.md. The k study's harness
# (binary-k-mixing-need.R) is sourced and its designs, selection, fits,
# classes, summaries, references and acceptance are reused; this file adds
# what pricing reads differently and decides nothing the plan leaves open
# except where its comments say so. The real-data side (the loader, the sealed
# scores and assignment, the guard's references, fits and reading, the
# held-out score and step 3) is in k-move-pricing-guard.R, sourced below.
#
# BUILDS. One R process sees one build, named by KMP_BUILD: "shipped" (the
# rebuilt move's commit with the engine toggle off, k-move-pricing-shipped.patch)
# or "move" (that commit as committed); "mutant1" and "mutant3" are step 1's
# power libraries, run as the move. Designs are made under shipped.
#
# ENVIRONMENT. KMP_RUN, the run directory; KMP_QUICK=1 runs every length at a
# tenth (the smoke; quick output is marked and never mixed with full);
# KMP_SEED_TAG (default "kmp") enters every seed of a truth and an arm, and a
# confirmation run names a fresh one in a run directory of its own;
# DBARTS_BIS_DATA, the binary studies' data root; KMP_KSTUDY_RUN, the k
# study's run directory (Reading 3); KMP_ARCHIVE, the private archive of the
# real datasets' files (outside the tree).
#
# Usage (each action one process):
#   Rscript k-move-pricing.R design <regime> <rep>
#   Rscript k-move-pricing.R timing fit <regime> <rep> <round>
#   Rscript k-move-pricing.R timing segment <cell> <round>
#   Rscript k-move-pricing.R timing state <cell>
#   Rscript k-move-pricing.R timing summarize
#   Rscript k-move-pricing.R fit <regime> <rep> shipped|move
#   Rscript k-move-pricing.R score <regime> <rep>
#   Rscript k-move-pricing.R ess <regime> <rep>
#   Rscript k-move-pricing.R sbc pilot|run <prior> <regime> <rep>
#   Rscript k-move-pricing.R stream <regime> <rep> default|dispersed <L>
#   Rscript k-move-pricing.R r3 fingerprint <tag> | r3 check <ref> <b>... | r3 stream <key> <L>
#   Rscript k-move-pricing.R reference run|extend|accept <cell> logistic 1 ...
#   Rscript k-move-pricing.R guard standins|load <name>|scores|assign
#   Rscript k-move-pricing.R guard design|pilot|more <name> | timing <name> <round>
#   Rscript k-move-pricing.R guard fits <name> <from> <to> [tag] | read [no995] | heldout
#   Rscript k-move-pricing.R step3 features|select [r3 size]|read|fresh
#   Rscript k-move-pricing.R read step1pilot|step1 [build [reps]]|step2 <look>
#   Rscript k-move-pricing.R read convergence|ess|r3
#   Rscript k-move-pricing.R check jobs|outputs <file> | conjugate|blind|eqtime
#   Rscript k-move-pricing.R check scoring|libs|streams
#   Rscript k-move-pricing.R jobs <phase> <file>
# Phases, in order: timing; step1:pilot, step1:power (read before) step1:run;
# step2:1, step2:2; ess; conv, conv:more, conv:long (each written after the
# last is read); r3; guard:prefit, guard:refs, guard:fits; smoke.
# Every library is built by k-move-pricing-install.sh, which stamps its
# commit; every output records it with its shared object's md5 (libId).
# A dataset's name in a command writes each space as %20. KMP_SMOKE=1 puts
# three simulated stand-ins in the sealed list's place, written to
# KMP_LOADER_SOURCE (the smoke reads no passer's file).

## ---------------------------------------------------------------- settings

kmpScriptDir <- function() {
  env <- Sys.getenv("KMP_HARNESS_DIR")
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

kmpBuild <- Sys.getenv("KMP_BUILD", Sys.getenv("KMN_BUILD"))
if (!kmpBuild %in% c("shipped", "move", "mutant1", "mutant3", "kstudy", "")) {
  stop("KMP_BUILD must be shipped, move, mutant1, mutant3 or kstudy")
}
local({
  map <- c(
    KMP_RUN = "KMN_RUN",
    KMP_QUICK = "KMN_QUICK"
  )
  for (nm in names(map)) {
    v <- Sys.getenv(nm)
    if (nzchar(v)) {
      do.call(Sys.setenv, stats::setNames(list(v), map[[nm]]))
    }
  }
  if (nzchar(kmpBuild)) {
    Sys.setenv(
      KMN_BUILD = if (kmpBuild %in% c("shipped", "kstudy")) {
        "shipped"
      } else {
        "move"
      }
    )
  }
  tag <- Sys.getenv("KMP_SEED_TAG", "kmp")
  Sys.setenv(KMN_SEED_TAG = tag, KMN_HARNESS_DIR = kmpScriptDir())
})
source(file.path(kmpScriptDir(), "binary-k-mixing-need.R"))
if (!requireNamespace("posterior", quietly = TRUE)) {
  stop("the posterior package is needed")
}

kmpSmoke <- nzchar(Sys.getenv("KMP_SMOKE"))
# A fixed r replaces the measured one in equalTimeDraws and in the
# equal-time check alike, so the check could not see it: only the smoke may
# set it, and every output records it.
if (nzchar(Sys.getenv("KMN_FIXED_R")) && !kmpSmoke) {
  stop("KMN_FIXED_R is refused outside the smoke (KMP_SMOKE=1)")
}

# The library a process runs: the build's name and commit that
# k-move-pricing-install.sh wrote into it, and the md5 of its shared object.
libId <- local({
  id <- NULL
  function() {
    if (is.null(id)) {
      pkg <- find.package("dbarts")
      dcf <- file.path(pkg, "kmp-build.dcf")
      so <- Sys.glob(file.path(pkg, "libs", "dbarts.*"))
      info <- if (file.exists(dcf)) read.dcf(dcf)[1L, ] else NULL
      id <<- c(
        build = kmpBuild,
        stamped = if (is.null(info)) NA_character_ else unname(info[["Build"]]),
        commit = if (is.null(info)) NA_character_ else unname(info[["Commit"]]),
        patch = if (is.null(info)) NA_character_ else unname(info[["Patch"]]),
        so = if (length(so)) unname(tools::md5sum(so[1L])) else NA_character_,
        fixedR = Sys.getenv("KMN_FIXED_R")
      )
    }
    id
  }
})

# Every output that carries the quick flag carries the library's identity.
kstudySaveAtomic <- saveAtomic
saveAtomic <- function(object, path) {
  if (is.list(object) && !is.data.frame(object) && !is.null(object$quick)) {
    object$libId <- libId()
  }
  kstudySaveAtomic(object, path)
}

kmpLink <- "logistic"
kmpFixedScale <- 2 # chi(1.5, 2), the default; step 1's tilted priors below

# The real-data guard's references: the first job records 16,000 sweeps, the
# checks read at 12,000 and 16,000, a block a round to 20,000, then the
# shipped half alone to 40,000, 80,000 and 100,000; one chain a part job, so
# that every chain starts from its own draw from the prior.
mcmc$refFirst <- tenth(16000)
refPartChains <- 1L

# Regimes: the k study's, logistic only, with step 2's counts and looks (Sizing
# and stops). A is the regimes with strong signal or rare outcomes, B the fits
# with none.
for (id in names(regimes)) {
  regimes[[id]]$links <- intersect(regimes[[id]]$links, kmpLink)
}
step2Looks <- list(
  S1 = c(100L, 200L),
  S2 = c(240L, 480L),
  S3 = c(50L, 100L),
  S4 = c(100L, 200L),
  S5 = c(135L, 270L),
  S6 = c(110L, 220L),
  S7 = c(150L, 300L),
  S8 = c(150L, 300L),
  E3 = c(110L, 220L)
)
setA <- c("S1", "S2", "S3", "S4")
setB <- c("S5", "S6", "S7", "S8", "E3")
step2Regimes <- c(setA, setB)
lookCount <- function(id, look) {
  n <- step2Looks[[id]][look]
  if (quick) max(2L, as.integer(round(n / 50))) else n
}

# The planned standard errors at full size (Sizing and stops), per regime and
# class, on the worse weighting; a first look's is scaled to its count.
plannedSe <- list(
  S1 = c(fInt = 0.014, fTail = 0.014),
  S2 = c(fInt = 0.010),
  S3 = c(fInt = 0.012),
  S4 = c(fInt = 0.005),
  S5 = c(fInt = 0.006, fTail = 0.005),
  S6 = c(fInt = 0.006, fTail = 0.006),
  S7 = c(ite = 0.006, subgroup = 0.006),
  S8 = c(ite = 0.006),
  E3 = c(fInt = 0.005, fTail = 0.005)
)

# The k study's five non-band problem cells, where a gain is read.
gainCells <- list(
  S4 = c("fInt", "fTail"),
  S1 = c("fInt", "fTail"),
  S2 = "fTail"
)
barClasses <- c("fInt", "fTail", "ite", "subgroup") # harmless below the bar
harmOnlyClasses <- c("ate", "att", "logRR", "logOR") # read for shown harm
bandClasses <- c("bandF", "bandIte") # reported, no bar (dec-B625)
step2Bar <- 0.02
zLook <- qnorm(1 - 0.0025) # 2.81: gain and harm, one-sided 0.0025 a look
zHarmless <- qnorm(1 - 0.025) # harmless, one-sided 0.025 a look
ratioBar <- 1.02 # RMSE and interval score, move over shipped
costBar <- 1.25 # any timed logistic cell (dec-B620)

kmpCheckQuick <- function(x, path) {
  if (!identical(x$quick, quick)) {
    stop("quick and full output mixed: ", path)
  }
  x
}
readChecked <- function(path) {
  if (!file.exists(path)) {
    stop("missing input ", path)
  }
  kmpCheckQuick(readRDS(path), path)
}

stampLoad <- function() {
  if (file.exists("/proc/loadavg")) {
    readLines("/proc/loadavg", warn = FALSE)
  } else {
    paste(
      tryCatch(system("uptime", intern = TRUE), error = function(e) NA),
      collapse = ""
    )
  }
}

leafPriorAt <- function(scale) {
  dbartsPriors$normal(dbartsPriors$chi(kDf, scale))
}

# A draw from the model's prior on a sampler: trees, k from chi(1.5, scale),
# leaves given k; with keepK the sampler's k prior is restored after the
# leaves, so its chain starts at the drawn k and draws k from there.
priorStart <- function(sampler, scale = kmpFixedScale, restore = TRUE) {
  sampler$sampleTreesFromPrior()
  k0 <- scale * sqrt(rchisq(1L, kDf))
  sampler$setLeafPrior(dbartsPriors$normal(k = k0))
  sampler$sampleLeafParametersFromPrior()
  if (restore) {
    sampler$setLeafPrior(leafPriorAt(scale))
  }
  k0
}

source(file.path(kmpScriptDir(), "k-move-pricing-guard.R"))

## ------------------------------------------------------------------ timing

# Equal running time (Common to steps 1 to 3, call 20). Per regime, three
# replications and three rounds a replication, the two builds alternating
# within a round (the job list's order). A piece times a full default fit as
# the user runs it and a one-sweep fit on the same data; the one-sweep fit's
# time is the minimum of three, since it is short and its noise is the
# timer's (the study's own choice). r is the median over the nine pairs of the
# move's remainder over the shipped build's, the ratio of seconds per sweep.

timingPiecePath <- function(cell, rep, round) {
  runDir(
    "timing",
    cellFile(cell),
    sprintf("%s-%s-t%d.rds", repTag(rep), kmpBuild, as.integer(round))
  )
}

defaultFitArgs <- function(design, chains, burn, draws, seed) {
  fd <- fitData(design)
  args <- list(
    y ~ .,
    data = fd$train,
    family = kmpLink,
    n.trees = nTrees,
    n.chains = chains,
    n.threads = 1L,
    n.burn = burn,
    n.samples = draws,
    seed = seed,
    verbose = FALSE,
    keepTrainingFits = design$task %in% c("prob", "loop")
  )
  if (!is.null(fd$test)) {
    args$test <- fd$test
  }
  args
}

timeFit <- function(design, chains, lengths, seed) {
  if (design$task == "loop") {
    return(fitLoop(design, kmpLink, chains, lengths, seed)$seconds)
  }
  args <- defaultFitArgs(design, chains, lengths[1L], lengths[2L], seed)
  system.time(do.call(bart, args))[["elapsed"]]
}

runTimingFit <- function(cell, rep, round) {
  design <- blind(loadDesign(cell, kmpLink, rep))
  chains <- cellChains(cell)
  seed <- seedOf("timing", cell, rep)
  one <- min(vapply(
    1:3,
    function(i) timeFit(design, chains, c(0L, 1L), seed),
    0
  ))
  full <- timeFit(design, chains, c(mcmc$burn, mcmc$defaultDraws), seed)
  sweeps <- mcmc$burn + mcmc$defaultDraws
  out <- list(
    cell = cell,
    rep = as.integer(rep),
    round = as.integer(round),
    build = kmpBuild,
    chains = chains,
    sweeps = sweeps,
    full = full,
    oneSweep = one,
    load = stampLoad(),
    quick = quick
  )
  saveAtomic(out, timingPiecePath(cell, rep, round))
  cat(sprintf(
    "timing %s %d round %d %s: full %.2f s, one sweep %.2f s\n",
    cell,
    as.integer(rep),
    as.integer(round),
    kmpBuild,
    full,
    one
  ))
}

# The cost benchmark's time-only cells: the move plan's large-cell generator
# (10 uniform predictors, eta = a (2 sin(pi x1 x2) + 2 (x3 - 0.5)^2 + x4 - 1),
# y = 1{eta + e > 0} with logistic e, a set so the sign of eta matches y on
# about 0.90 or 0.76 of rows), 10,000 and 50,000 rows by 75 and 200 trees;
# person-period rows of hazard.logistic and a count-weighted logistic cell.
# Timed as 100-sweep segments, both builds restored from one stored state
# before every segment, ten interleaved repetitions in alternating order.
timedCells <- c(
  outer(
    c("large:10000", "large:50000"),
    c(":75:sep", ":75:mod", ":200:sep", ":200:mod"),
    paste0
  ),
  "hazard:2000:75",
  "counts:2000:75"
)
segmentSweeps <- 100L

largeSignal <- function(target, seed) {
  withSeed(seed, {
    x <- matrix(runif(20000L * 10L), ncol = 10L)
    base <- 2 *
      sin(pi * x[, 1L] * x[, 2L]) +
      2 * (x[, 3L] - 0.5)^2 +
      x[, 4L] -
      1
    e <- rlogis(nrow(x))
    uniroot(
      function(a) mean(sign(a * base) == sign(a * base + e)) - target,
      c(0.01, 100)
    )$root
  })
}

timedCellData <- function(cell) {
  p <- strsplit(cell, ":", fixed = TRUE)[[1L]]
  n <- as.integer(p[2L])
  if (quick) {
    n <- n %/% 10L
  }
  seed <- seedOf("timed-cell", cell)
  if (p[1L] == "large") {
    a <- largeSignal(
      if (p[4L] == "sep") 0.90 else 0.76,
      seedOf("signal", p[4L])
    )
    return(withSeed(seed, {
      x <- matrix(runif(n * 10L), ncol = 10L)
      eta <- a *
        (2 * sin(pi * x[, 1L] * x[, 2L]) + 2 * (x[, 3L] - 0.5)^2 + x[, 4L] - 1)
      y <- as.integer(eta + rlogis(n) > 0)
      list(
        data = data.frame(x, y = y),
        family = "logistic",
        trees = as.integer(p[3L]),
        agreement = mean(y == as.integer(eta > 0))
      )
    }))
  }
  if (p[1L] == "hazard") {
    # person-period rows: a subject's rows to its event or censoring
    return(withSeed(seed, {
      subjects <- n %/% 5L
      x <- matrix(runif(subjects * 5L), ncol = 5L)
      time <- pmin(rgeom(subjects, plogis(x[, 1L] - 2)) + 1L, 10L)
      status <- as.integer(time < 10L)
      list(
        data = data.frame(x, time = time, status = status),
        family = "hazard.logistic",
        trees = as.integer(p[3L])
      )
    }))
  }
  withSeed(seed, {
    x <- matrix(runif(n * 5L), ncol = 5L)
    size <- sample.int(5L, n, replace = TRUE)
    succ <- rbinom(n, size, plogis(2 * x[, 1L] - 1))
    list(
      data = data.frame(
        x[c(seq_len(n), seq_len(n)), , drop = FALSE],
        y = c(rep(1L, n), rep(0L, n)),
        w = c(succ, size - succ)
      ),
      family = "logistic",
      trees = as.integer(p[3L]),
      weights = TRUE
    )
  })
}

timedSampler <- function(cell) {
  d <- timedCellData(cell)
  args <- list(
    data = d$data,
    family = d$family,
    n.trees = d$trees,
    n.chains = 1L,
    n.threads = 1L,
    seed = seedOf("timed-sampler", cell),
    verbose = FALSE,
    keepTrainingFits = FALSE,
    samplerOnly = TRUE
  )
  if (d$family == "hazard.logistic") {
    args$formula <- survival::Surv(time, status) ~ .
  } else {
    args$formula <- y ~ .
  }
  if (isTRUE(d$weights)) {
    keep <- d$data$w > 0
    args$data <- d$data[keep, , drop = FALSE]
    args$weights <- args$data$w
    args$data$w <- NULL
    args$weights.type <- "frequency"
  }
  do.call(bart, args)
}

timedStatePath <- function(cell) {
  runDir("timing", "state", paste0(cellFile(cell), ".rds"))
}

# The stored state both builds restore before every segment: a shipped
# sampler past its default burn-in.
runTimedState <- function(cell) {
  requireBuild("shipped")
  s <- timedSampler(cell)
  invisible(s$run(mcmc$burn, 0L))
  saveSampler(s, timedStatePath(cell))
  cat("timing state", cell, "\n")
}

runTimedSegment <- function(cell, round) {
  s <- readRDS(timedStatePath(cell))
  seconds <- system.time(invisible(s$run(0L, segmentSweeps)))[["elapsed"]]
  saveAtomic(
    list(
      cell = cell,
      rep = 1L,
      round = as.integer(round),
      build = kmpBuild,
      seconds = seconds,
      sweeps = segmentSweeps,
      load = stampLoad(),
      quick = quick
    ),
    timingPiecePath(cell, 1L, round)
  )
  cat(sprintf(
    "segment %s round %d %s: %.3f s\n",
    cell,
    as.integer(round),
    kmpBuild,
    seconds
  ))
}

# Per cell: r, its range, the share of fixed time in a whole fit, and each
# build's seconds per chain-sweep (the convergence measure's clock).
summarizeTimingPieces <- function() {
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
  if (length(pieces) == 0L) {
    stop("no timing pieces")
  }
  rows <- lapply(pieces, function(p) {
    remainder <- if (is.null(p$full)) p$seconds else p$full - p$oneSweep
    data.frame(
      cell = p$cell,
      rep = p$rep,
      round = p$round,
      build = p$build,
      remainder = remainder,
      full = if (is.null(p$full)) NA_real_ else p$full,
      oneSweep = if (is.null(p$oneSweep)) NA_real_ else p$oneSweep,
      chainSweeps = if (is.null(p$full)) {
        p$sweeps
      } else {
        p$chains * (p$sweeps - 1L)
      }
    )
  })
  d <- do.call(rbind, rows)
  w <- reshape(
    d,
    idvar = c("cell", "rep", "round"),
    timevar = "build",
    direction = "wide"
  )
  w$ratio <- w$remainder.move / w$remainder.shipped
  tab <- do.call(
    rbind,
    lapply(split(w, w$cell), function(g) {
      data.frame(
        link = kmpLink,
        cell = g$cell[1L],
        n = sum(is.finite(g$ratio)),
        r = stats::median(g$ratio, na.rm = TRUE),
        rMin = min(g$ratio, na.rm = TRUE),
        rMax = max(g$ratio, na.rm = TRUE),
        fixedShare = stats::median(
          g$oneSweep.shipped / g$full.shipped,
          na.rm = TRUE
        ),
        secShipped = stats::median(
          g$remainder.shipped / g$chainSweeps.shipped,
          na.rm = TRUE
        ),
        secMove = stats::median(
          g$remainder.move / g$chainSweeps.move,
          na.rm = TRUE
        ),
        costBarExceeded = stats::median(g$ratio, na.rm = TRUE) > costBar
      )
    })
  )
  rownames(tab) <- NULL
  tab$quick <- quick
  saveAtomic(tab, timingFile())
  print(tab, row.names = FALSE)
  invisible(tab)
}

ratioFor <- function(cell) {
  timeRatio(cell, kmpLink)
}

secondsPerSweep <- function(cell, which) {
  tab <- readRDS(timingFile())
  v <- tab[[if (which == "shipped") "secShipped" else "secMove"]][
    tab$cell == cell
  ]
  if (length(v) != 1L) NA_real_ else v
}

# Equal time: the shipped arm's draws a chain after the same burn-in, so that
# its sweeps take the time the move's default sweeps take.
equalTimeDraws <- function(cell) {
  unit <- mcmc$burn + mcmc$defaultDraws
  as.integer(ceiling(unit * ratioFor(cell))) - mcmc$burn
}

## ------------------------------------------------------------- step 2 arms

step2FitPath <- function(cell, rep, arm) {
  runDir("p2fits", cellFile(cell), sprintf("%s-%s.rds", repTag(rep), arm))
}
step2ScorePath <- function(cell, rep) {
  runDir("p2scores", cellFile(cell), paste0(repTag(rep), ".rds"))
}

kMixing <- function(k) {
  if (is.null(k) || nrow(k) < 4L || any(!is.finite(k))) {
    return(c(rhat = NA_real_, essBulk = NA_real_, essTail = NA_real_))
  }
  c(
    rhat = posterior::rhat(k),
    essBulk = posterior::ess_bulk(k),
    essTail = posterior::ess_tail(k)
  )
}

readingOf <- function(res, design, idx) {
  L <- poolDraws(res$arr, idx)
  cl <- classesFrom(L, design, kmpLink)
  list(
    summary = withProbMeans(summarizeAll(cl), cl),
    mixing = kMixing(if (is.null(res$k)) NULL else res$k[idx, , drop = FALSE]),
    draws = length(idx),
    maxAbsLatent = max(abs(L))
  )
}

# One arm of a step 2 replication: the shipped build at equal time (its
# default reading the prefix of the same run) or the move at the default
# length; both after bart()'s 500 burn-in (dec-B621).
runStep2Arm <- function(cell, rep, arm) {
  requireBuild(arm)
  design <- blind(loadDesign(cell, kmpLink, rep))
  chains <- cellChains(cell)
  n <- if (arm == "shipped") equalTimeDraws(cell) else mcmc$defaultDraws
  seed <- armSeed(cell, kmpLink, rep, arm)
  # the shipped arm's default reading needs the default's draws even where
  # equal time gives fewer (r below 1, timing noise at a tenth)
  nRun <- max(n, mcmc$defaultDraws)
  res <- if (design$task == "loop") {
    fitLoop(design, kmpLink, chains, c(mcmc$burn, nRun), seed)
  } else {
    fitBart(design, kmpLink, chains, c(mcmc$burn, nRun), seed)
  }
  idx <- list(default = seq_len(mcmc$defaultDraws))
  if (arm == "shipped") {
    idx$eq <- seq_len(n)
  }
  scoring <- system.time({
    readings <- lapply(idx, function(i) readingOf(res, design, i))
  })[["elapsed"]]
  out <- list(
    cell = cell,
    rep = as.integer(rep),
    arm = arm,
    build = kmpBuild,
    lib = find.package("dbarts"),
    drawsPerChain = n,
    seconds = res$seconds,
    scoringSeconds = scoring,
    readings = readings,
    quick = quick
  )
  saveAtomic(out, step2FitPath(cell, rep, arm))
  cat(sprintf(
    "fit %s %d %s: %d draws a chain, %.1f s, scored %.1f s\n",
    cell,
    as.integer(rep),
    arm,
    n,
    res$seconds,
    scoring
  ))
  invisible(out)
}

# Labels: the shipped arm's equal-time reading pooled with the move's (the two
# arms' 5th and 95th percentiles averaged), a function of the data and the
# chains, never of the truth (the k study's rule).
pricingLabels <- function(shipped, move) {
  ends <- function(s) {
    q <- cbind(s$fInt$q, s$fTail$q)
    cols <- c(s$fInt$cols, s$fTail$cols)
    o <- order(cols)
    list(
      lo = q[match(0.05, round(probsKept, 6)), o],
      hi = q[match(0.95, round(probsKept, 6)), o]
    )
  }
  a <- ends(shipped$readings$eq$summary)
  b <- ends(move$readings$default$summary)
  plogis((a$lo + b$lo) / 2) >= interiorRange[1L] &
    plogis((a$hi + b$hi) / 2) <= interiorRange[2L]
}

# The posterior mean of the probability beside a latent row class's
# summary (pmean), which the probability-scale errors read.
withProbMeans <- function(summary, classes) {
  for (nm in names(summary)) {
    if (isTRUE(classInfo[[nm]]$latent) && ncol(classes[[nm]]) > 0L) {
      summary[[nm]]$pmean <- colMeans(plogis(classes[[nm]]))
    }
  }
  summary
}

# relabel keeping each column's probability mean.
relabelP <- function(summary, labels) {
  pm <- summary$fInt$pmean
  out <- relabel(summary, labels)
  if (!is.null(pm) && !is.null(out$fInt)) {
    out$fInt$pmean <- pm[labels]
    out$fTail$pmean <- pm[!labels]
  }
  out
}

# Per class: coverage counts at 90 percent, the squared error and the signed
# error of the posterior mean, the interval's width and its interval score on
# the scale the class is summarized on (latent for the row classes), and the
# same on the probability scale (sqP, errP, widthP, is90P: the posterior
# mean of the probability, the interval's ends transformed), which rule 2's
# ratio bar reads, as the k study scored them.
pricingCoverage <- function(summary, truth) {
  out <- list()
  for (nm in intersect(names(summary), names(truth))) {
    s <- summary[[nm]]
    t <- truth[[nm]]
    if (isTRUE(classInfo[[nm]]$band)) {
      b <- s$bands[["0.9"]]
      out[[nm]] <- c(n = 1, cov90 = as.numeric(all(b$lo <= t & t <= b$hi)))
      next
    }
    if (length(t) == 0L) {
      out[[nm]] <- c(
        n = 0,
        cov90 = 0,
        sq = 0,
        err = 0,
        width = 0,
        is90 = 0,
        sqP = 0,
        errP = 0,
        widthP = 0,
        is90P = 0
      )
      next
    }
    i90 <- intervalOf(s, 0.9)
    latent <- isTRUE(classInfo[[nm]]$latent)
    pm <- if (latent) s$pmean else s$mean
    tp <- if (latent) plogis(t) else t
    loP <- if (latent) plogis(i90$lo) else i90$lo
    hiP <- if (latent) plogis(i90$hi) else i90$hi
    if (is.null(pm)) {
      pm <- rep(NA_real_, length(t))
    }
    out[[nm]] <- c(
      n = length(t),
      cov90 = sum(i90$lo <= t & t <= i90$hi),
      sq = sum((s$mean - t)^2),
      err = sum(s$mean - t),
      width = sum(i90$hi - i90$lo),
      is90 = sum(intervalScore(i90$lo, i90$hi, t, 0.9)),
      sqP = sum((pm - tp)^2),
      errP = sum(pm - tp),
      widthP = sum(hiP - loP),
      is90P = sum(intervalScore(loP, hiP, tp, 0.9))
    )
  }
  out
}

runStep2Score <- function(cell, rep) {
  design <- loadDesign(cell, kmpLink, rep)
  ship <- readChecked(step2FitPath(cell, rep, "shipped"))
  move <- readChecked(step2FitPath(cell, rep, "move"))
  labels <- if (design$task %in% c("prob", "loop")) {
    pricingLabels(ship, move)
  } else {
    NULL
  }
  truth <- truthClasses(design, labels)
  read1 <- function(fit, r) {
    s <- fit$readings[[r]]$summary
    if (!is.null(labels)) {
      s <- relabelP(s, labels)
    }
    pricingCoverage(s, truth)
  }
  out <- list(
    cell = cell,
    rep = as.integer(rep),
    coverage = list(
      shipped.eq = read1(ship, "eq"),
      shipped.default = read1(ship, "default"),
      move.default = read1(move, "default")
    ),
    labels = labels,
    rate = design$rate,
    score = design$score,
    k0 = design$truth$k0,
    kMixing = list(
      shipped.eq = ship$readings$eq$mixing,
      move.default = move$readings$default$mixing
    ),
    maxAbsLatent = c(
      shipped = ship$readings$eq$maxAbsLatent,
      move = move$readings$default$maxAbsLatent
    ),
    seconds = c(shipped = ship$seconds, move = move$seconds),
    drawsPerChain = c(shipped = ship$drawsPerChain, move = move$drawsPerChain),
    quick = quick
  )
  saveAtomic(out, step2ScorePath(cell, rep))
  cat(sprintf("score %s %d\n", cell, as.integer(rep)))
  invisible(out)
}

## ------------------------------------------- effective draws a second

essRegimes <- c("S1", "S2", "S4", "S6")

# 4 chains of 5,000 draws after 1,000 burn-in, both builds: each row latent's
# bulk and tail effective draws and k's, over the fit's seconds.
runEss <- function(cell, rep) {
  design <- blind(loadDesign(cell, kmpLink, rep))
  res <- fitBart(
    design,
    kmpLink,
    4L,
    c(tenth(1000), tenth(5000)),
    taggedSeedOf("ess", cell, rep)
  )
  S <- dim(res$arr)[2L]
  essCol <- function(j, f) f(matrix(res$arr[,, j], S, 4L, byrow = TRUE))
  nCol <- dim(res$arr)[3L]
  bulk <- vapply(seq_len(nCol), essCol, 0, f = posterior::ess_bulk)
  tailE <- vapply(seq_len(nCol), essCol, 0, f = posterior::ess_tail)
  out <- list(
    cell = cell,
    rep = as.integer(rep),
    build = kmpBuild,
    seconds = res$seconds,
    rowBulk = stats::quantile(bulk, c(0.01, 0.1, 0.5), names = FALSE),
    rowTail = stats::quantile(tailE, c(0.01, 0.1, 0.5), names = FALSE),
    k = kMixing(res$k),
    quick = quick
  )
  saveAtomic(
    out,
    runDir("ess", cellFile(cell), sprintf("%s-%s.rds", repTag(rep), kmpBuild))
  )
  cat(sprintf(
    "ess %s %d %s: %.1f s\n",
    cell,
    as.integer(rep),
    kmpBuild,
    res$seconds
  ))
}

## ------------------------------------------------------------- step 1: SBC

# Step 1 (simulation-based calibration of the move): one chain, 99 draws at
# thin T after burn-in B, on S1, S2 and S4's covariates, under the default
# chi(1.5, 2) with the k study's selection, chi(1.5, 0.1) with the bands and
# the deterministic rejection dropped and S2's rate test kept, and the
# contingent chi(1.5, 0.5) likewise. Four functionals are ranked: k, the mean
# latent fit over training rows and the latent fit at two held-out rows. Two
# rows are always held out (from the dataset's rows outside the training
# sample, or the training sample's size less two where the dataset has no
# more rows), the study's own reading of "fixed by index in advance".
step1Priors <- c(default = 2, tilted = 0.1, third = 0.5)
step1Regimes <- c("S1", "S2", "S4")
step1Draws <- 99L
step1Cutoffs <- list(
  default = c(0.2, 0.5),
  tilted = c(0.02, 0.05),
  third = c(0.2, 0.5)
)
# Replications a regime: the pilot 20, the power run 50 (150 a prior), the
# main run 200 (600 a prior).
step1Reps <- function(phase) {
  n <- c(pilot = 20L, power = 50L, run = 200L)[[phase]]
  if (quick) max(2L, n %/% 20L) else n
}
step1PilotLength <- function() tenth(20000)

step1Design <- function(regimeId, prior, rep) {
  reg <- regimes[[regimeId]]
  scale <- step1Priors[[prior]]
  key <- paste("step1", regimeId, prior, rep, seedTag)
  dataset <- reg$datasets[(as.integer(rep) - 1L) %% length(reg$datasets) + 1L]
  raw <- realDatasets[[dataset]]()
  N <- nrow(raw$x)
  nTrain <- min(reg$n, N - 2L)
  rows <- withSeed(seedOf(key, "rows"), sample.int(N, nTrain + 2L))
  x <- raw$x[rows[seq_len(nTrain)], , drop = FALSE]
  xh <- raw$x[rows[nTrain + 1:2], , drop = FALSE]
  rownames(x) <- NULL
  rownames(xh) <- NULL
  sampler <- priorSampler(x, kmpLink, seedOf(key, "sampler"))
  maxAttempts <- if (reg$select == "rare") 200000L else 5000L
  for (attempt in seq_len(maxAttempts)) {
    draw <- withSeed(seedOf(key, "prior", attempt), {
      k0 <- priorStart(sampler, scale, restore = FALSE)
      list(
        k0 = k0,
        f = as.numeric(sampler$predict(x)),
        fh = as.numeric(sampler$predict(xh))
      )
    })
    y <- withSeed(
      seedOf(key, "y", attempt),
      rbinom(nTrain, 1L, plogis(draw$f))
    )
    if (length(unique(y)) < 2L) {
      next
    }
    ok <- if (prior == "default") {
      selects(reg, x, y, paste(key, attempt))$ok
    } else if (reg$select == "rare") {
      rate <- mean(y)
      (rate >= 0.01 && rate <= 0.05) || (rate >= 0.95 && rate <= 0.99)
    } else {
      TRUE
    }
    if (ok) break
  }
  if (!ok) {
    stop("no step 1 replication passed selection")
  }
  list(
    regime = regimeId,
    prior = prior,
    rep = as.integer(rep),
    dataset = dataset,
    x = x,
    xh = xh,
    y = y,
    truth = draw,
    attempts = attempt,
    score = separationScore(x, y, paste(key, "score")),
    quick = quick
  )
}

step1Path <- function(phase, prior, regimeId, rep) {
  runDir(
    "step1",
    phase,
    kmpBuild,
    prior,
    regimeId,
    paste0(repTag(rep), ".rds")
  )
}

step1Settings <- function(prior) {
  path <- runDir("step1", "settings.rds")
  if (!file.exists(path)) {
    stop("no step 1 thin and burn-in; read the pilot first")
  }
  s <- readChecked(path)$settings
  s[[prior]]
}

# The pilot runs long (no thin, step1PilotLength sweeps) and reports each
# functional's effective draws per sweep and autocorrelation time; the run
# keeps the 99 thinned draws and ranks the truth among them.
runStep1 <- function(phase, prior, regimeId, rep) {
  if (!kmpBuild %in% c("move", "mutant1", "mutant3")) {
    stop("step 1 runs the move or its mutants")
  }
  d <- step1Design(regimeId, prior, rep)
  scale <- step1Priors[[prior]]
  if (phase == "pilot") {
    burn <- 0L
    thin <- 1L
    draws <- step1PilotLength()
  } else {
    st <- step1Settings(prior)
    burn <- st$burn
    thin <- st$thin
    # bart()'s n.samples is the sweeps recorded before thinning
    draws <- step1Draws * thin
  }
  s <- bart(
    y ~ .,
    data = data.frame(d$x, y = d$y),
    test = d$xh,
    family = kmpLink,
    n.trees = nTrees,
    n.chains = 1L,
    n.threads = 1L,
    n.burn = burn,
    n.samples = draws,
    n.thin = thin,
    seed = taggedSeedOf("step1-fit", regimeId, prior, rep),
    verbose = FALSE,
    samplerOnly = TRUE
  )
  s$setLeafPrior(leafPriorAt(scale))
  seconds <- system.time(r <- s$run())[["elapsed"]]
  k <- as.numeric(r$k)
  if (phase != "pilot" && length(k) != step1Draws) {
    stop("step 1 kept ", length(k), " draws, not ", step1Draws)
  }
  fun <- cbind(
    k = k,
    meanFit = colMeans(r$train),
    held1 = r$test[1L, ],
    held2 = r$test[2L, ]
  )
  truth <- c(
    k = d$truth$k0,
    meanFit = mean(d$truth$f),
    held1 = d$truth$fh[1L],
    held2 = d$truth$fh[2L]
  )
  finite <- all(is.finite(r$train)) &&
    all(is.finite(r$test)) &&
    all(is.finite(k))
  out <- list(
    phase = phase,
    prior = prior,
    regime = regimeId,
    rep = as.integer(rep),
    build = kmpBuild,
    score = d$score,
    attempts = d$attempts,
    k0 = d$truth$k0,
    finite = finite,
    maxAbsLatent = max(abs(r$train), abs(r$test)),
    seconds = seconds,
    quick = quick
  )
  if (phase == "pilot") {
    out$essPerSweep <- apply(fun, 2L, function(v) {
      posterior::ess_bulk(matrix(v, ncol = 1L)) / length(v)
    })
    out$tau <- apply(fun, 2L, function(v) {
      length(v) / posterior::ess_basic(matrix(v, ncol = 1L))
    })
  } else {
    out$ranks <- vapply(
      colnames(fun),
      function(nm) sum(fun[, nm] < truth[[nm]]),
      0L
    )
    out$essRatio <- apply(fun, 2L, function(v) {
      posterior::ess_bulk(matrix(v, ncol = 1L)) / length(v)
    })
    cut <- step1Cutoffs[[prior]]
    out$smallK <- vapply(
      cut,
      function(c0) mean(k < c0) - as.numeric(d$truth$k0 < c0),
      0
    )
    names(out$smallK) <- format(cut)
    # rows near 0 or 1 by the fit's own draws, covered by the interval
    # between the 5th and 95th of the 99 draws (order statistics, the
    # plan's interval and its target's); the near-end rows' effective-draw
    # ratio sets the target
    ends <- apply(r$train, 1L, function(v) sort.int(v)[c(5L, 95L)])
    lo <- ends[1L, ]
    hi <- ends[2L, ]
    nearEnd <- !(plogis(lo) >= interiorRange[1L] &
      plogis(hi) <= interiorRange[2L])
    out$tailCover <- c(
      n = sum(nearEnd),
      covered = sum((lo <= d$truth$f & d$truth$f <= hi)[nearEnd])
    )
    out$tailRho <- if (any(nearEnd)) {
      mean(vapply(
        which(nearEnd),
        function(i) {
          min(
            1,
            posterior::ess_bulk(matrix(r$train[i, ], ncol = 1L)) / step1Draws
          )
        },
        0
      ))
    } else {
      NA_real_
    }
    out$draws <- ncol(r$train)
  }
  saveAtomic(out, step1Path(phase, prior, regimeId, rep))
  cat(sprintf(
    "step1 %s %s %s %d %s: %.1f s\n",
    phase,
    prior,
    regimeId,
    as.integer(rep),
    kmpBuild,
    seconds
  ))
}

# The pilot's reading: T the smallest multiple of 50 keeping at least 50
# effective draws of 99 on every functional in every pilot replication, B at
# least nine times the slowest tau, a thin past 500 a mixing finding; the
# separation scores' distribution, the third arm's trigger (over half the
# tilted arm at 0.995 or more).
readStep1Pilot <- function() {
  settings <- list()
  report <- list()
  for (prior in names(step1Priors)) {
    files <- Sys.glob(runDir("step1", "pilot", "move", prior, "*", "r*.rds"))
    if (length(files) == 0L) {
      next
    }
    p <- lapply(files, readChecked)
    eps <- do.call(rbind, lapply(p, `[[`, "essPerSweep"))
    tau <- do.call(rbind, lapply(p, `[[`, "tau"))
    worst <- min(eps)
    thin <- max(50L, as.integer(50L * ceiling((50 / 99) / worst / 50)))
    burn <- as.integer(ceiling(9 * max(tau)))
    if (quick) {
      thin <- max(1L, thin %/% 10L)
      burn <- max(10L, burn %/% 10L)
    }
    scores <- vapply(p, `[[`, 0, "score")
    settings[[prior]] <- list(thin = thin, burn = burn)
    report[[prior]] <- list(
      reps = length(p),
      thin = thin,
      burn = burn,
      mixingFinding = thin > 500L,
      worstEssPerSweep = worst,
      slowestTau = max(tau),
      scores = stats::quantile(scores, c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE),
      shareDeterministic = mean(scores >= deterministicScore, na.rm = TRUE)
    )
  }
  third <- !is.null(report$tilted) && report$tilted$shareDeterministic > 0.5
  if (third) {
    # the contingent arm's jobs wait for this file
    ensureDir(runDir("step1"))
    writeLines("third arm on", runDir("step1", "third-go"))
  }
  out <- list(
    settings = settings,
    report = report,
    runThird = third,
    quick = quick
  )
  saveAtomic(out, runDir("step1", "settings.rds"))
  str(report)
  cat("third arm:", third, "\n")
  invisible(out)
}

# Step 1's rules (fixed in the plan): ranks within the ecdf band and a
# chi-square p over 20 bins at or above 0.05 / 16 (0.05 / 24 with the third
# arm); the small-k share within 3 bootstrap standard errors of 0; coverage
# near 0 or 1 within 3.29 standard errors of the target at the measured
# effective-draw ratio; every draw finite. The early look at 150 a prior goes
# to the confirming rerun on a p below alpha / 10.
step1Target <- function(rho) {
  # coverage, for an exact sampler, of the interval between the 5th and 95th
  # of 99 draws of a stationary AR(1) chain whose effective-draw ratio
  # (1 - phi) / (1 + phi) is rho, the truth an independent draw (the plan's
  # 0.898 independent, about 0.8945 at 0.5 and 0.885 at 0.25)
  rho <- min(max(rho, 0.01), 1)
  phi <- (1 - rho) / (1 + rho)
  B <- if (quick) 20000L else 200000L
  withSeed(seedOf("step1-target", round(rho, 4)), {
    x <- matrix(0, step1Draws, B)
    x[1L, ] <- rnorm(B)
    e <- sqrt(1 - phi^2)
    for (t in 2:step1Draws) {
      x[t, ] <- phi * x[t - 1L, ] + e * rnorm(B)
    }
    truth <- rnorm(B)
    ends <- apply(x, 2L, function(v) sort.int(v)[c(5L, 95L)])
    mean(ends[1L, ] <= truth & truth <= ends[2L, ])
  })
}

sbcEnv <- new.env()
sbcRankUniformity <- function(...) {
  if (is.null(sbcEnv$rankUniformity)) {
    sys.source(file.path(kmpScriptDir(), "sbc.R"), envir = sbcEnv)
  }
  sbcEnv$rankUniformity(...)
}

# A mutant's power run reads rule 1 on the tilted prior only, at the main
# run's alpha; `reps` keeps the first replications of each regime (50 for
# the early look at 150 a prior).
readStep1 <- function(build = "move", reps = NA_integer_) {
  pilot <- readChecked(runDir("step1", "settings.rds"))
  priors <- names(step1Priors)
  if (!pilot$runThird) {
    priors <- setdiff(priors, "third")
  }
  # two rank tests on each of four functionals a prior
  alpha <- 0.05 / (8L * length(priors))
  mutant <- build %in% c("mutant1", "mutant3")
  if (mutant) {
    priors <- "tilted"
  }
  rows <- list()
  verdict <- list()
  for (prior in priors) {
    files <- Sys.glob(runDir("step1", "run", build, prior, "*", "r*.rds"))
    if (length(files) == 0L) {
      next
    }
    p <- lapply(files, readChecked)
    if (is.finite(reps)) {
      p <- Filter(function(x) x$rep <= reps, p)
    }
    ranks <- do.call(rbind, lapply(p, `[[`, "ranks"))
    for (nm in colnames(ranks)) {
      u <- sbcRankUniformity(ranks[, nm], step1Draws, alpha = alpha)
      rows[[length(rows) + 1L]] <- data.frame(
        prior = prior,
        functional = nm,
        reps = nrow(ranks),
        chisqP = u$chisqP,
        ecdfPass = u$pass,
        flag = !u$pass || u$chisqP < alpha,
        early = u$chisqP < alpha / 10
      )
    }
    sk <- do.call(rbind, lapply(p, `[[`, "smallK"))
    boot <- withSeed(seedOf("step1-boot", prior), {
      replicate(
        bootstrapResamples,
        colMeans(sk[sample.int(nrow(sk), replace = TRUE), , drop = FALSE])
      )
    })
    se <- apply(matrix(boot, nrow = ncol(sk)), 1L, stats::sd)
    tc <- do.call(rbind, lapply(p, `[[`, "tailCover"))
    draws <- vapply(p, `[[`, 0L, "draws")
    if (any(draws != step1Draws)) {
      stop("step 1 replications hold other than ", step1Draws, " draws")
    }
    nTail <- sum(tc[, "n"])
    covTail <- if (nTail > 0) sum(tc[, "covered"]) / nTail else NA_real_
    # the near-end rows' effective-draw ratio, weighted by their count
    tr <- vapply(
      p,
      function(x) if (is.null(x$tailRho)) NA_real_ else x$tailRho,
      0
    )
    rho <- if (nTail > 0) {
      stats::weighted.mean(tr, tc[, "n"], na.rm = TRUE)
    } else {
      NA_real_
    }
    target <- if (is.finite(rho)) step1Target(rho) else NA_real_
    # its standard error by a bootstrap over replications, rows of one fit
    # resampled together
    seTail <- if (nTail > 0) {
      bt <- withSeed(seedOf("step1-tail-boot", prior), {
        replicate(bootstrapResamples, {
          i <- sample.int(nrow(tc), replace = TRUE)
          if (sum(tc[i, "n"]) > 0) {
            sum(tc[i, "covered"]) / sum(tc[i, "n"])
          } else {
            NA_real_
          }
        })
      })
      stats::sd(bt, na.rm = TRUE)
    } else {
      NA_real_
    }
    verdict[[prior]] <- list(
      smallK = colMeans(sk),
      smallKSe = se,
      smallKPass = all(abs(colMeans(sk)) <= 3 * se),
      tailCoverage = covTail,
      tailTarget = target,
      tailRho = rho,
      tailSe = seTail,
      tailN = nTail,
      tailPass = is.na(covTail) || abs(covTail - target) <= 3.29 * seTail,
      allFinite = all(vapply(p, `[[`, NA, "finite")),
      maxAbsLatent = max(vapply(p, `[[`, 0, "maxAbsLatent"))
    )
  }
  ranksTab <- do.call(rbind, rows)
  if (mutant) {
    # the power check: the mutant must flag under rule 1
    out <- list(
      build = build,
      alpha = alpha,
      ranks = ranksTab,
      flagged = any(ranksTab$flag),
      reps = reps,
      quick = quick
    )
    saveAtomic(out, runDir("readings", sprintf("step1-%s.rds", build)))
    print(ranksTab, row.names = FALSE)
    cat(
      "step 1 power:",
      build,
      if (out$flagged) "flagged" else "NOT flagged",
      "\n"
    )
    return(invisible(out))
  }
  if (is.finite(reps)) {
    return(readStep1Early(build, reps, ranksTab, verdict, alpha))
  }
  pass <- !any(ranksTab$flag) &&
    all(vapply(
      verdict,
      function(v) v$smallKPass && v$tailPass && v$allFinite,
      NA
    ))
  out <- list(
    build = build,
    alpha = alpha,
    ranks = ranksTab,
    rules = verdict,
    pass = pass,
    quick = quick
  )
  saveAtomic(
    out,
    runDir(
      "readings",
      if (is.finite(reps)) {
        sprintf("step1-%s-first%d.rds", build, as.integer(reps))
      } else {
        sprintf("step1-%s.rds", build)
      }
    )
  )
  print(ranksTab, row.names = FALSE)
  str(verdict)
  cat("step 1", build, if (pass) "PASS" else "FLAG", "\n")
  invisible(out)
}

# Step 1's early look (the plan's stop, at 50 a regime, 150 a prior): a p
# below alpha / 10 on any functional goes to the confirming rerun, a
# non-finite draw ends step 1 as failed; otherwise the run goes on to 200 a
# regime (early-go). Only that rule decides here; the full rules are kept in
# the file and read at the end.
readStep1Early <- function(build, reps, ranksTab, verdict, alpha) {
  finite <- all(vapply(verdict, `[[`, NA, "allFinite"))
  hit <- ranksTab[ranksTab$early, c("prior", "functional"), drop = FALSE]
  decision <- if (!finite) {
    "failed: a non-finite draw"
  } else if (nrow(hit)) {
    "to the confirming rerun"
  } else {
    "go on"
  }
  out <- list(
    build = build,
    reps = reps,
    alpha = alpha,
    ranks = ranksTab,
    rules = verdict,
    replications = tapply(ranksTab$reps, ranksTab$prior, max),
    decision = decision,
    quick = quick
  )
  saveAtomic(
    out,
    runDir("readings", sprintf("step1-%s-first%d.rds", build, as.integer(reps)))
  )
  if (decision == "go on") {
    writeLines("early look passed", runDir("step1", "early-go"))
  }
  print(
    ranksTab[, c("prior", "functional", "reps", "chisqP", "early")],
    row.names = FALSE
  )
  cat(sprintf(
    "step 1 early look (%s at alpha / 10 = %.2g): %s\n",
    paste(
      sprintf("%s %d", names(out$replications), out$replications),
      collapse = ", "
    ),
    alpha / 10,
    decision
  ))
  invisible(out)
}

## ------------------------------------------------------- step 2 reading

step2Scores <- function(id, n) {
  paths <- vapply(seq_len(n), function(r) step2ScorePath(id, r), "")
  present <- file.exists(paths)
  list(scores = lapply(paths[present], readChecked), missing = sum(!present))
}

# One class of one regime at a look, under both weightings: each arm's
# coverage, its noise-corrected distance from 0.90 (dec-B455), the worsening W
# = d(move) - d(shipped at equal time) with its bootstrap standard error, and
# the readings; the ratios of RMSE and interval score; bias and width.
classReading <- function(sc, cl, idx) {
  get1 <- function(arm, field) {
    vapply(
      sc,
      function(s) {
        v <- s$coverage[[arm]][[cl]]
        if (is.null(v) || is.na(v[field])) NA_real_ else v[[field]]
      },
      0
    )
  }
  armM <- "move.default"
  armS <- "shipped.eq"
  n <- get1(armM, "n")
  cm <- get1(armM, "cov90")
  cs <- get1(armS, "cov90")
  ok <- is.finite(n) & n > 0 & is.finite(cm) & is.finite(cs)
  if (!any(ok)) {
    return(NULL)
  }
  # a resample over every replication keeps those holding this class, each
  # mapped to its place among them
  okAt <- cumsum(ok)
  idx <- lapply(idx, function(i) okAt[i[ok[i]]])
  n <- n[ok]
  cm <- cm[ok]
  cs <- cs[ok]
  reps <- length(n)
  band <- cl %in% bandClasses
  weights <- c("pooled", "perfit")
  covOf <- function(i, c, w) {
    if (w == "pooled") sum(c[i]) / sum(n[i]) else mean(c[i] / n[i])
  }
  rows <- list()
  for (w in weights) {
    if (w == "perfit" && (band || cl %in% harmOnlyClasses)) {
      next
    }
    full <- seq_len(reps)
    covM <- covOf(full, cm, w)
    covS <- covOf(full, cs, w)
    bootM <- vapply(idx, function(i) covOf(i, cm, w), 0)
    bootS <- vapply(idx, function(i) covOf(i, cs, w), 0)
    floorSe <- sqrt(0.09 / sum(n))
    # a resample holding no replication with this class is dropped
    sdAtLeast <- function(x) {
      s <- stats::sd(x[is.finite(x)])
      if (is.na(s)) floorSe else max(s, floorSe)
    }
    seM <- sdAtLeast(bootM)
    seS <- sdAtLeast(bootS)
    W <- correctedDistanceFloor(abs(covM - 0.9), seM) -
      correctedDistanceFloor(abs(covS - 0.9), seS)
    seW <- sdAtLeast(abs(bootM - 0.9) - abs(bootS - 0.9))
    harmless <- !band && W + zHarmless * seW < step2Bar
    harm <- !band && W >= step2Bar && W - zLook * seW > 0
    gain <- !band && -W >= step2Bar && -W - zLook * seW > 0
    rows[[w]] <- data.frame(
      class = cl,
      weighting = w,
      reps = reps,
      n = sum(n),
      covMove = covM,
      covShipped = covS,
      ciMove = paste(
        format(
          stats::quantile(bootM, c(0.025, 0.975), na.rm = TRUE),
          digits = 3
        ),
        collapse = " to "
      ),
      ciShipped = paste(
        format(
          stats::quantile(bootS, c(0.025, 0.975), na.rm = TRUE),
          digits = 3
        ),
        collapse = " to "
      ),
      diff = covM - covS,
      diffSe = stats::sd(bootM - bootS, na.rm = TRUE),
      W = W,
      seW = seW,
      harmless = harmless && !(cl %in% harmOnlyClasses),
      harm = harm,
      gain = gain,
      needAtZero = as.integer(ceiling(
        ((zHarmless + qnorm(0.9)) * seW * sqrt(reps) / step2Bar)^2
      ))
    )
  }
  tab <- do.call(rbind, rows)
  if (!band) {
    # the bar on the probability scale; the latent scale reported beside
    ratioOn <- function(fsq, fis) {
      sM <- get1(armM, fsq)[ok]
      sS <- get1(armS, fsq)[ok]
      iM <- get1(armM, fis)[ok]
      iS <- get1(armS, fis)[ok]
      rm <- function(i) {
        sqrt(sum(sM[i]) / sum(n[i])) / sqrt(sum(sS[i]) / sum(n[i]))
      }
      ir <- function(i) sum(iM[i]) / sum(iS[i])
      list(
        rmse = rm(seq_len(reps)),
        rmseUpper = stats::quantile(
          vapply(idx, rm, 0),
          0.975,
          na.rm = TRUE,
          names = FALSE
        ),
        is = ir(seq_len(reps)),
        isUpper = stats::quantile(
          vapply(idx, ir, 0),
          0.975,
          na.rm = TRUE,
          names = FALSE
        )
      )
    }
    lat <- ratioOn("sq", "is90")
    tab$rmseRatioLatent <- lat$rmse
    tab$rmseRatioLatentUpper <- lat$rmseUpper
    tab$isRatioLatent <- lat$is
    tab$isRatioLatentUpper <- lat$isUpper
    sqM <- get1(armM, "sqP")[ok]
    sqS <- get1(armS, "sqP")[ok]
    isM <- get1(armM, "is90P")[ok]
    isS <- get1(armS, "is90P")[ok]
    rmse <- function(i) {
      sqrt(sum(sqM[i]) / sum(n[i])) / sqrt(sum(sqS[i]) / sum(n[i]))
    }
    isr <- function(i) sum(isM[i]) / sum(isS[i])
    bR <- vapply(idx, rmse, 0)
    bI <- vapply(idx, isr, 0)
    tab$rmseRatio <- rmse(seq_len(reps))
    tab$rmseRatioUpper <- stats::quantile(
      bR,
      0.975,
      na.rm = TRUE,
      names = FALSE
    )
    tab$isRatio <- isr(seq_len(reps))
    tab$isRatioUpper <- stats::quantile(bI, 0.975, na.rm = TRUE, names = FALSE)
    tab$ratiosHarmless <- tab$rmseRatioUpper < ratioBar &
      tab$isRatioUpper < ratioBar
    tab$biasMove <- sum(get1(armM, "errP")[ok]) / sum(n)
    tab$biasShipped <- sum(get1(armS, "errP")[ok]) / sum(n)
    tab$widthMove <- sum(get1(armM, "widthP")[ok]) / sum(n)
    tab$widthShipped <- sum(get1(armS, "widthP")[ok]) / sum(n)
    tab$rmseMove <- sqrt(sum(sqM) / sum(n))
    tab$rmseShipped <- sqrt(sum(sqS) / sum(n))
    tab$rmseMoveLatent <- sqrt(sum(get1(armM, "sq")[ok]) / sum(n))
    tab$rmseShippedLatent <- sqrt(sum(get1(armS, "sq")[ok]) / sum(n))
  }
  tab
}

# A class resolves when it is harmless under both weightings with its ratios
# harmless, or shows a harm under either; a one-value estimand resolves only
# by a harm shown, else stays read; a band is reported.
classStatus <- function(tab) {
  cl <- tab$class[1L]
  if (cl %in% bandClasses) {
    return("reported")
  }
  if (any(tab$harm)) {
    return("harm shown")
  }
  if (cl %in% harmOnlyClasses) {
    return("no harm shown")
  }
  ratiosOk <- is.null(tab$ratiosHarmless) || all(tab$ratiosHarmless)
  if (all(tab$harmless) && ratiosOk) {
    if (all(tab$gain)) "improves" else "harmless"
  } else {
    "unresolved"
  }
}

# Rows of data frames with differing columns (a band class carries no
# ratios), a column a frame lacks filled with NA.
rbindFill <- function(l) {
  cols <- unique(unlist(lapply(l, names)))
  do.call(
    rbind,
    lapply(l, function(d) {
      for (cl in setdiff(cols, names(d))) {
        d[[cl]] <- NA
      }
      d[cols]
    })
  )
}

readStep2 <- function(look) {
  out <- list()
  all <- list()
  for (id in step2Regimes) {
    n <- lookCount(id, look)
    got <- step2Scores(id, n)
    if (length(got$scores) == 0L) {
      next
    }
    reps <- length(got$scores)
    idx <- withSeed(seedOf("step2-boot", id, look), {
      lapply(seq_len(bootstrapResamples), function(b) {
        sample.int(reps, reps, replace = TRUE)
      })
    })
    classes <- unique(unlist(lapply(got$scores, function(s) {
      names(s$coverage$move.default)
    })))
    tabs <- Filter(
      Negate(is.null),
      lapply(classes, function(cl) classReading(got$scores, cl, idx))
    )
    for (t in tabs) {
      t$regime <- id
      t$set <- if (id %in% setA) "A" else "B"
      t$status <- classStatus(t)
      t$gainCell <- t$class[1L] %in% gainCells[[id]]
      planned <- plannedSe[[id]][t$class[1L]]
      full <- step2Looks[[id]][2L]
      t$plannedSe <- if (is.na(planned)) NA_real_ else planned * sqrt(full / n)
      t$restatePower <- look == 1L &
        !is.na(t$plannedSe) &
        t$seW > 1.5 * t$plannedSe
      all[[length(all) + 1L]] <- t
    }
    km <- do.call(rbind, lapply(got$scores, function(s) s$kMixing$move.default))
    out[[id]] <- list(
      reps = reps,
      missing = got$missing,
      stops = all(vapply(
        tabs,
        function(t) classStatus(t) != "unresolved",
        NA
      )),
      kMixing = c(
        shareRhat = mean(km[, "rhat"] < 1.01, na.rm = TRUE),
        shareEssBulk = mean(km[, "essBulk"] >= 400, na.rm = TRUE),
        shareEssTail = mean(km[, "essTail"] >= 400, na.rm = TRUE)
      ),
      maxAbsLatent = max(vapply(got$scores, function(s) max(s$maxAbsLatent), 0))
    )
  }
  if (length(all) == 0L) {
    stop("no step 2 scores to read at look ", look)
  }
  tab <- rbindFill(all)
  status <- tapply(
    seq_len(nrow(tab)),
    paste(tab$regime, tab$class),
    function(i) {
      tab$status[i[1L]]
    }
  )
  gainMet <- any(tab$gainCell & tab$weighting == "pooled" & tab$gain) &&
    any(tab$gainCell & tab$weighting == "perfit" & tab$gain) &&
    any(vapply(
      split(tab[tab$gainCell, ], paste(tab$regime, tab$class)[tab$gainCell]),
      function(t) all(t$gain),
      NA
    ))
  harmA <- any(tab$set == "A" & tab$status == "harm shown")
  harmB <- any(tab$set == "B" & tab$status == "harm shown")
  unresolved <- unique(paste(tab$regime, tab$class)[tab$status == "unresolved"])
  # the real-data guard is one more class of rule 2: unresolved until its
  # reading says harmless (a harm it shows stands only once confirmed)
  guardFile <- runDir("readings", "guard.rds")
  guardVerdict <- if (file.exists(guardFile)) {
    readChecked(guardFile)$verdict
  } else {
    "not read"
  }
  if (!identical(guardVerdict, "harmless")) {
    unresolved <- c(unresolved, paste0("R2 guard (", guardVerdict, ")"))
  }
  timing <- if (file.exists(timingFile())) readRDS(timingFile()) else NULL
  costCell <- !is.null(timing) && any(timing$costBarExceeded, na.rm = TRUE)
  costTooHigh <- (harmB && !harmA && gainMet) || costCell
  framing <- if (!gainMet || harmA) {
    "neither a default nor a trigger is supported: parking (step 4)"
  } else if (costTooHigh) {
    "the cost is too high: step 3"
  } else if (length(unresolved)) {
    "no recommendation past what is shown: classes unresolved at these counts"
  } else {
    "the move on for every logistic fit it reaches"
  }
  # a harm frames only once confirmed on fresh seeds (a rerun in a run
  # directory of its own); until then the framing it leads to is preliminary
  if (harmA || harmB) {
    framing <- paste0(
      "preliminary (a harm shown, not yet confirmed on fresh seeds): ",
      framing
    )
  }
  res <- list(
    look = look,
    table = tab,
    regimes = out,
    gainMet = gainMet,
    harmA = harmA,
    harmB = harmB,
    costCell = costCell,
    costTooHigh = costTooHigh,
    unresolved = unresolved,
    guardVerdict = guardVerdict,
    framing = framing,
    status = status,
    quick = quick
  )
  saveAtomic(res, runDir("readings", sprintf("step2-look%d.rds", look)))
  print(
    tab[
      tab$weighting == "pooled" | tab$class %in% barClasses,
      c(
        "regime",
        "class",
        "weighting",
        "reps",
        "covMove",
        "covShipped",
        "W",
        "seW",
        "status"
      )
    ],
    row.names = FALSE,
    digits = 3
  )
  cat("gain:", gainMet, " harm in A:", harmA, " in B:", harmB, "\n")
  cat("framing:", framing, "\n")
  invisible(res)
}

readEss <- function() {
  files <- Sys.glob(runDir("ess", "*", "r*-*.rds"))
  p <- lapply(files, readChecked)
  tab <- do.call(
    rbind,
    lapply(p, function(x) {
      data.frame(
        cell = x$cell,
        rep = x$rep,
        build = x$build,
        seconds = x$seconds,
        rowBulkMedian = x$rowBulk[3L],
        rowTailP01 = x$rowTail[1L],
        kBulk = x$k[["essBulk"]],
        bulkPerSecond = x$rowBulk[3L] / x$seconds,
        tailPerSecond = x$rowTail[1L] / x$seconds,
        kPerSecond = x$k[["essBulk"]] / x$seconds
      )
    })
  )
  agg <- stats::aggregate(
    cbind(bulkPerSecond, tailPerSecond, kPerSecond) ~ cell + build,
    tab,
    stats::median
  )
  saveAtomic(
    list(table = tab, summary = agg, quick = quick),
    runDir("readings", "ess.rds")
  )
  print(agg, row.names = FALSE, digits = 3)
  invisible(agg)
}

## ------------------------------------------------- the convergence measure

# Streams: 20 chains of each build from the start with no burn-in, read as the
# fits a user runs (five fits of 4 chains, two of 10 in the effect regimes),
# every draw from sweep 1 to L kept, at every length of the series. Each
# fit's classes are formed as step 2 forms them, labeled by step 2's score of
# the replication; bands from 250 sweeps. R-hat and effective draws are
# computed for k and a fixed subset of each class's rows (diagRows a class,
# evenly spaced by column), the study's own economy: every row at every
# length would cost several times the streams' own sampling.
convRegimes <- c("S1", "S2", "S4", "S5", "S6", "S7", "S8", "E3")
convLengthsFull <- c(
  25L,
  50L,
  100L,
  200L,
  300L,
  500L,
  700L,
  1000L,
  1500L,
  2000L,
  3000L,
  4000L,
  6000L,
  8000L
)
convLengths <- function() unique(tenth(convLengthsFull))
convBandFrom <- function() tenth(250)
convChains <- 20L
convDispersedReps <- 25L
diagRows <- 20L

streamPath <- function(cell, rep, start, which = kmpBuild) {
  runDir(
    "conv",
    cellFile(cell),
    sprintf("%s-%s-%s.rds", repTag(rep), start, which)
  )
}
diagPath <- function(cell, rep, start, which = kmpBuild) {
  runDir(
    "conv",
    "diag",
    cellFile(cell),
    sprintf("%s-%s-%s.f32", repTag(rep), start, which)
  )
}

# One chain of a stream: the draws (sweeps x columns, the fit's own order) and
# k, from the default start or the chain's own draw from the prior.
streamChain <- function(design, seed, start, L) {
  task <- design$task
  fd <- fitData(design)
  mask <- if (task == "loop") design$mask else integer()
  train <- fd$train
  if (task == "loop") {
    train$y[mask] <- withSeed(
      seedOf("stream-mask-start", seed),
      rbinom(length(mask), 1L, mean(design$y[-mask]))
    )
  }
  args <- list(
    y ~ .,
    data = train,
    family = kmpLink,
    n.trees = nTrees,
    n.chains = 1L,
    n.threads = 1L,
    n.burn = 0L,
    n.samples = L,
    seed = seed,
    verbose = FALSE,
    samplerOnly = TRUE,
    keepTrainingFits = task %in% c("prob", "loop")
  )
  if (!is.null(fd$test)) {
    args$test <- fd$test
  }
  t0 <- proc.time()[["elapsed"]]
  s <- do.call(bart, args)
  if (start == "dispersed") {
    withSeed(seedOf("stream-start", seed), priorStart(s))
  }
  if (task == "loop") {
    nCols <- nrow(design$x) + if (is.null(design$xh)) 0L else nrow(design$xh)
    draws <- matrix(0, L, nCols)
    k <- numeric(L)
    yy <- train$y
    withSeed(seedOf("stream-mask", seed), {
      for (it in seq_len(L)) {
        r <- s$run(0L, 1L)
        f <- as.numeric(r$train)
        yy[mask] <- rbinom(length(mask), 1L, plogis(f[mask]))
        s$setResponse(yy)
        draws[it, ] <- c(f, if (is.null(r$test)) NULL else as.numeric(r$test))
        k[it] <- as.numeric(r$k)[1L]
      }
    })
  } else {
    r <- s$run(0L, L)
    m <- if (task == "prob") {
      if (is.null(r$test)) r$train else rbind(r$train, r$test)
    } else {
      r$test
    }
    draws <- t(matrix(m, nrow = nrow(m)))
    k <- as.numeric(r$k)
  }
  list(draws = draws, k = k, seconds = proc.time()[["elapsed"]] - t0)
}

# The diagnostic quantities of a draws matrix: per class, diagRows columns
# (for the effect regimes, individual effects at diagRows units), each draw's
# value; the column sets are fixed by index from the labels.
diagIndex <- function(design, labels, nCols) {
  if (design$task %in% c("prob", "loop")) {
    pick <- function(cols) {
      if (length(cols) <= diagRows) {
        cols
      } else {
        cols[unique(round(seq(1, length(cols), length.out = diagRows)))]
      }
    }
    list(fInt = pick(which(labels)), fTail = pick(which(!labels)))
  } else {
    n <- nCols %/% 2L
    units <- unique(round(seq(1, n, length.out = min(n, diagRows))))
    list(ite = units)
  }
}

diagValues <- function(draws, design, index) {
  if (design$task %in% c("prob", "loop")) {
    draws[, unlist(index), drop = FALSE]
  } else {
    n <- ncol(draws) %/% 2L
    u <- index$ite
    plogis(draws[, u, drop = FALSE]) - plogis(draws[, n + u, drop = FALSE])
  }
}

# Rank-normalized split R-hat with bulk and tail effective draws of a
# sweeps-by-chains matrix (Vehtari et al. 2021), and the split pooled and
# within-chain variances Brooks and Gelman plot.
mixingAt <- function(m) {
  if (nrow(m) < 4L || any(!is.finite(m))) {
    return(c(rhat = NA, essBulk = NA, essTail = NA, W = NA, varPlus = NA))
  }
  half <- nrow(m) %/% 2L
  sp <- cbind(
    m[seq_len(half), , drop = FALSE],
    m[half + seq_len(half), , drop = FALSE]
  )
  W <- mean(apply(sp, 2L, stats::var))
  B <- half * stats::var(colMeans(sp))
  c(
    rhat = posterior::rhat(m),
    essBulk = posterior::ess_bulk(m),
    essTail = posterior::ess_tail(m),
    W = W,
    varPlus = (half - 1) / half * W + B / half
  )
}

# The scorer of one fit at one length as first written: the pooled draws of
# the prefix, its classes, their summaries (kept for the equality check).
streamSummaryAt <- function(arr, bl, labels, L) {
  Ld <- poolDraws(arr, seq_len(L))
  cl <- classesFrom(Ld, bl, kmpLink)
  if (L < convBandFrom()) {
    cl <- cl[setdiff(names(cl), bandClasses)]
  }
  summ <- withProbMeans(summarizeAll(cl), cl)
  if (!is.null(labels)) {
    summ <- relabelP(summ, labels)
  }
  summ
}

# The same summaries at every length from one pass: every class is a
# per-draw function of the latent draws and a prefix's pooled draws are rows
# of the longest's in the same chain-major order, so the classes are built
# once and sliced; each column is ordered once and a prefix's order is the
# longest's restricted to the prefix's rows (ties keep their row order, as
# order() does). The quantiles are R's type 7, computed as quantile.default
# computes them, and the bands from the ranks as summarizeDraws forms them,
# so every summary is bit-identical to streamSummaryAt's.
streamSummaries <- function(arr, bl, labels, lengths) {
  nChains <- dim(arr)[1L]
  Lmax <- max(lengths)
  clAll <- classesFrom(poolDraws(arr, seq_len(Lmax)), bl, kmpLink)
  sweepOf <- rep.int(seq_len(Lmax), nChains)
  ords <- lapply(clAll, function(Q) {
    if (ncol(Q) == 0L) {
      return(NULL)
    }
    o <- apply(Q, 2L, order)
    if (!is.matrix(o)) matrix(o, nrow(Q), ncol(Q)) else o
  })
  lapply(lengths, function(L) {
    keep <- sweepOf <= L
    rows <- which(keep)
    localOf <- integer(length(keep))
    localOf[rows] <- seq_along(rows)
    names1 <- names(clAll)
    if (L < convBandFrom()) {
      names1 <- setdiff(names1, bandClasses)
    }
    summ <- lapply(names1, function(nm) {
      prefixSummary(clAll[[nm]], nm, rows, keep, localOf, ords[[nm]])
    })
    names(summ) <- names1
    for (nm in c("fInt", "fTail")) {
      if (!is.null(summ[[nm]])) {
        summ[[nm]]$cols <- attr(clAll[[nm]], "cols")
      }
    }
    for (nm in names1) {
      if (isTRUE(classInfo[[nm]]$latent) && ncol(clAll[[nm]]) > 0L) {
        summ[[nm]]$pmean <- colMeans(plogis(clAll[[nm]][rows, , drop = FALSE]))
      }
    }
    if (!is.null(labels)) {
      summ <- relabelP(summ, labels)
    }
    summ
  })
}

# summarizeClass on the prefix rows of Q, given Q's column orders.
prefixSummary <- function(Qall, name, rows, keep, localOf, ordAll) {
  info <- classInfo[[name]]
  Q <- Qall[rows, , drop = FALSE]
  S <- nrow(Q)
  m <- ncol(Q)
  if (m == 0L) {
    return(list(
      q = matrix(0, length(probsKept), 0L),
      mean = numeric(),
      sd = numeric(),
      S = S
    ))
  }
  o <- matrix(localOf[ordAll[keep[ordAll]]], S, m)
  col <- rep(seq_len(m), each = S)
  sorted <- matrix(Q[cbind(as.vector(o), col)], S, m)
  index <- 1 + max(S - 1L, 0) * probsKept
  lo <- floor(index)
  hi <- ceiling(index)
  qs <- sorted[lo, , drop = FALSE]
  xh <- sorted[hi, , drop = FALSE]
  h <- matrix(index - lo, length(probsKept), m)
  i <- which((index > lo) & (xh != qs))
  qs[i] <- (1 - h[i]) * qs[i] + h[i] * xh[i]
  if (!is.null(colnames(Q))) {
    dimnames(qs) <- list(NULL, colnames(Q))
  }
  out <- list(
    q = qs,
    mean = colMeans(Q),
    sd = apply(Q, 2L, stats::sd),
    S = S
  )
  if (isTRUE(info$band)) {
    out$bands <- prefixBands(sorted, o, name)
  }
  out
}

# summarizeDraws' simultaneous bands (and per-unit bands for bandUnit)
# from a prefix's sorted columns and their orders.
prefixBands <- function(sorted, o, name) {
  S <- nrow(sorted)
  m <- ncol(sorted)
  levels <- c(0.9, 0.95)
  tail <- min(S, ceiling(0.06 * S) + 2L)
  ranks <- matrix(0L, S, m)
  ranks[cbind(as.vector(o), rep(seq_len(m), each = S))] <- rep.int(
    seq_len(S),
    m
  )
  ext <- pmax(ranks, S + 1L - ranks)
  if (name == "bandUnit") {
    groupSize <- 6L
    units <- m %/% groupSize
    ext3 <- array(ext, c(S, groupSize, units))
    eg <- ext3[, 1L, ]
    for (k in seq_len(groupSize)[-1L]) {
      eg <- pmax(eg, ext3[, k, ])
    }
    eg <- matrix(eg, S, units)
    bands <- lapply(levels, function(l) {
      pos <- ceiling(l * S)
      kStar <- apply(eg, 2L, function(v) sort.int(v, partial = pos)[pos])
      kCol <- rep(kStar, each = groupSize)
      list(
        lo = sorted[cbind(S + 1L - kCol, seq_len(m))],
        hi = sorted[cbind(kCol, seq_len(m))]
      )
    })
  } else {
    extreme <- ext[cbind(seq_len(S), max.col(ext, "first"))]
    bands <- lapply(levels, function(l) {
      kStar <- sort.int(extreme)[ceiling(l * S)]
      if (S + 1L - kStar > tail) {
        stop("band level outside the stored tails")
      }
      list(lo = sorted[S + 1L - kStar, ], hi = sorted[kStar, ])
    })
  }
  names(bands) <- c("0.9", "0.95")
  bands
}

runStream <- function(cell, rep, start, Lmax) {
  if (!kmpBuild %in% c("shipped", "move")) {
    stop("a stream belongs to the shipped or the move build")
  }
  design <- loadDesign(cell, kmpLink, rep)
  sc <- readChecked(step2ScorePath(cell, rep))
  labels <- sc$labels
  truth <- truthClasses(design, labels)
  bl <- blind(design)
  per <- cellChains(cell)
  nFits <- convChains %/% per
  lengths <- convLengths()
  lengths <- lengths[lengths <= Lmax]
  index <- NULL
  diag <- list()
  kAll <- matrix(0, Lmax, convChains)
  seconds <- numeric(convChains)
  fits <- vector("list", nFits)
  scoring <- logical()
  for (f in seq_len(nFits)) {
    chainIds <- (f - 1L) * per + seq_len(per)
    runs <- lapply(chainIds, function(ch) {
      streamChain(bl, taggedSeedOf("stream", cell, rep, start, ch), start, Lmax)
    })
    arr <- aperm(simplify2array(lapply(runs, `[[`, "draws")), c(3L, 1L, 2L))
    kM <- vapply(runs, `[[`, numeric(Lmax), "k")
    kAll[, chainIds] <- kM
    seconds[chainIds] <- vapply(runs, `[[`, 0, "seconds")
    if (is.null(index)) {
      index <- diagIndex(bl, labels, dim(arr)[3L])
    }
    for (j in seq_along(chainIds)) {
      diag[[chainIds[j]]] <- diagValues(runs[[j]]$draws, bl, index)
    }
    rm(runs)
    summaries <- streamSummaries(arr, bl, labels, lengths)
    if (kmpSmoke) {
      same <- vapply(
        seq_along(lengths),
        function(li) {
          identical(
            summaries[[li]],
            streamSummaryAt(arr, bl, labels, lengths[li])
          )
        },
        NA
      )
      scoring <- c(scoring, same)
    }
    fits[[f]] <- lapply(seq_along(lengths), function(li) {
      L <- lengths[li]
      idx <- seq_len(L)
      summ <- summaries[[li]]
      dv <- simplify2array(lapply(chainIds, function(ch) {
        diag[[ch]][idx, , drop = FALSE]
      }))
      mixQ <- apply(dv, 2L, function(m) mixingAt(matrix(m, L))[1:3])
      list(
        L = L,
        coverage = pricingCoverage(summ, truth),
        k = mixingAt(kM[idx, , drop = FALSE]),
        rows = mixQ
      )
    })
    rm(arr)
  }
  across <- lapply(lengths, function(L) {
    idx <- seq_len(L)
    dv <- simplify2array(lapply(diag, function(d) d[idx, , drop = FALSE]))
    list(
      L = L,
      k = mixingAt(kAll[idx, , drop = FALSE]),
      rows = apply(dv, 2L, function(m) mixingAt(matrix(m, L))[1:3])
    )
  })
  both <- NULL
  if (rep <= convDispersedReps) {
    keep <- simplify2array(diag)
    writeSingle(
      array(c(keep, kAll), c(Lmax, ncol(diag[[1L]]) + 1L, convChains)),
      diagPath(cell, rep, start)
    )
    other <- diagPath(
      cell,
      rep,
      if (start == "default") "dispersed" else "default"
    )
    otherMeta <- streamPath(
      cell,
      rep,
      if (start == "default") "dispersed" else "default"
    )
    if (file.exists(other) && file.exists(otherMeta)) {
      om <- readRDS(otherMeta)
      if (identical(om$quick, quick) && om$Lmax == Lmax) {
        od <- readSingle(other, c(Lmax, ncol(diag[[1L]]) + 1L, convChains))
        mine <- array(c(keep, kAll), dim(od))
        both <- lapply(lengths, function(L) {
          idx <- seq_len(L)
          comb <- array(0, c(L, dim(od)[2L], 2L * convChains))
          comb[,, seq_len(convChains)] <- mine[idx, , , drop = FALSE]
          comb[,, convChains + seq_len(convChains)] <- od[idx, , , drop = FALSE]
          list(
            L = L,
            rhat = apply(comb, 2L, function(m) posterior::rhat(matrix(m, L)))
          )
        })
      }
    }
  }
  out <- list(
    cell = cell,
    rep = as.integer(rep),
    start = start,
    build = kmpBuild,
    Lmax = Lmax,
    lengths = lengths,
    chainsPerFit = per,
    fits = fits,
    across = across,
    both = both,
    seconds = seconds,
    diagIndex = index,
    scoringSame = scoring,
    quick = quick
  )
  saveAtomic(out, streamPath(cell, rep, start))
  cat(sprintf(
    "stream %s %d %s %s to %d: %.1f s\n",
    cell,
    as.integer(rep),
    start,
    kmpBuild,
    as.integer(Lmax),
    sum(seconds)
  ))
  invisible(out)
}

# The regimes whose first look goes on to 100 replications (the plan's
# stop): a row class with a planned standard error whose per-fit coverage
# standard error at 500 or 1,000 sweeps, default start, either library,
# exceeds 1.5 times its planned value at the look's count (step 2's
# planned figure at its full count, scaled to the replications read).
convMoreRegimes <- function() {
  f <- runDir("readings", "convergence.rds")
  if (!file.exists(f)) {
    stop("the convergence first look is read before its second is written")
  }
  tab <- readChecked(f)$coverage
  at <- tab[
    tab$start == "default" & tab$L %in% tenth(c(500L, 1000L)),
    ,
    drop = FALSE
  ]
  more <- vapply(
    convRegimes,
    function(id) {
      t <- at[at$regime == id & at$class %in% names(plannedSe[[id]]), ]
      if (nrow(t) == 0L) {
        return(FALSE)
      }
      planned <- plannedSe[[id]][t$class] * sqrt(step2Looks[[id]][2L] / t$reps)
      any(t$perfitSe > 1.5 * planned, na.rm = TRUE)
    },
    NA
  )
  convRegimes[more]
}

# Whether a regime's streams run on to 8,000 for a build: k's R-hat over a
# replication's 20 chains at 2,000 sweeps at or above 1.01 in more than half
# the replications of the first look (the study's own rule).
streamsRunOn <- function(cell, which, reps) {
  r <- vapply(
    seq_len(reps),
    function(rep) {
      p <- streamPath(cell, rep, "default", which)
      if (!file.exists(p)) {
        return(NA_real_)
      }
      s <- readChecked(p)
      a <- s$across[[match(tenth(2000L), s$lengths)]]
      if (is.null(a)) NA_real_ else a$k[["rhat"]]
    },
    0
  )
  mean(r >= 1.01, na.rm = TRUE) > 0.5
}

# Reading 1 and 2 by regime, start, build and length: coverage pooled and per
# fit with bootstrap intervals over replications, length (width), RMSE and
# bias; R-hat and effective draws of k and of the rows by their distribution
# over replications; the first length from which R-hat stays below 1.01; each
# length in seconds from the timing's seconds per sweep (the move's whole
# cost charged).
readConvergence <- function() {
  rows <- list()
  mix <- list()
  for (cell in convRegimes) {
    for (start in c("default", "dispersed")) {
      for (which in c("shipped", "move")) {
        files <- Sys.glob(runDir(
          "conv",
          cellFile(cell),
          sprintf("r*-%s-%s.rds", start, which)
        ))
        if (length(files) == 0L) {
          next
        }
        st <- lapply(files, readChecked)
        sps <- if (file.exists(timingFile())) {
          secondsPerSweep(cell, which)
        } else {
          NA
        }
        lengths <- sort(unique(unlist(lapply(st, `[[`, "lengths"))))
        # seconds are a user's fit of chainsPer chains on one thread
        chainsPer <- st[[1L]]$chainsPerFit
        for (L in lengths) {
          # a type-7 90 percent interval of n exact independent draws covers
          # 0.9 (n - 1) / (n + 1) in expectation, the curves' floor
          nDraws <- L * chainsPer
          exactCov <- 0.9 * (nDraws - 1) / (nDraws + 1)
          perRep <- lapply(st, function(s) {
            i <- match(L, s$lengths)
            if (is.na(i)) NULL else lapply(s$fits, `[[`, i)
          })
          perRep <- Filter(Negate(is.null), perRep)
          classes <- unique(unlist(lapply(perRep, function(fl) {
            names(fl[[1L]]$coverage)
          })))
          reps <- length(perRep)
          idx <- withSeed(seedOf("conv-boot", cell, start, which, L), {
            lapply(seq_len(bootstrapResamples), function(b) {
              sample.int(reps, reps, replace = TRUE)
            })
          })
          for (cl in classes) {
            agg <- t(vapply(
              perRep,
              function(fl) {
                # a replication may hold no row of a class (no tail row)
                v <- Filter(
                  Negate(is.null),
                  lapply(fl, function(f) f$coverage[[cl]])
                )
                if (length(v) == 0L) {
                  return(c(
                    n = 0,
                    cov = 0,
                    perfit = NA,
                    sq = 0,
                    err = 0,
                    width = 0,
                    sqP = 0,
                    errP = 0,
                    widthP = 0
                  ))
                }
                n <- sum(vapply(v, function(x) x[["n"]], 0))
                cov <- sum(vapply(v, function(x) x[["cov90"]], 0))
                perfit <- mean(
                  vapply(
                    v,
                    function(x) {
                      if (x[["n"]] > 0) x[["cov90"]] / x[["n"]] else NA
                    },
                    0
                  ),
                  na.rm = TRUE
                )
                get <- function(f) {
                  sum(vapply(
                    v,
                    function(x) if (f %in% names(x)) x[[f]] else NA,
                    0
                  ))
                }
                c(
                  n = n,
                  cov = cov,
                  perfit = perfit,
                  sq = get("sq"),
                  err = get("err"),
                  width = get("width"),
                  sqP = get("sqP"),
                  errP = get("errP"),
                  widthP = get("widthP")
                )
              },
              numeric(9L)
            ))
            pooled <- function(i) sum(agg[i, "cov"]) / sum(agg[i, "n"])
            pf <- function(i) mean(agg[i, "perfit"], na.rm = TRUE)
            bP <- vapply(idx, pooled, 0)
            bF <- vapply(idx, pf, 0)
            rows[[length(rows) + 1L]] <- data.frame(
              regime = cell,
              start = start,
              build = which,
              L = L,
              seconds = L * sps * chainsPer,
              class = cl,
              reps = reps,
              exactCoverage = exactCov,
              pooled = pooled(seq_len(reps)),
              pooledLo = stats::quantile(
                bP,
                0.025,
                na.rm = TRUE,
                names = FALSE
              ),
              pooledHi = stats::quantile(
                bP,
                0.975,
                na.rm = TRUE,
                names = FALSE
              ),
              perfit = pf(seq_len(reps)),
              perfitLo = stats::quantile(
                bF,
                0.025,
                na.rm = TRUE,
                names = FALSE
              ),
              perfitHi = stats::quantile(
                bF,
                0.975,
                na.rm = TRUE,
                names = FALSE
              ),
              perfitSe = stats::sd(bF, na.rm = TRUE),
              width = sum(agg[, "widthP"]) / sum(agg[, "n"]),
              rmse = sqrt(sum(agg[, "sqP"]) / sum(agg[, "n"])),
              bias = sum(agg[, "errP"]) / sum(agg[, "n"]),
              widthLatent = sum(agg[, "width"]) / sum(agg[, "n"]),
              rmseLatent = sqrt(sum(agg[, "sq"]) / sum(agg[, "n"])),
              biasLatent = sum(agg[, "err"]) / sum(agg[, "n"])
            )
          }
          kFit <- unlist(lapply(perRep, function(fl) {
            vapply(fl, function(f) f$k[["rhat"]], 0)
          }))
          kEss <- unlist(lapply(perRep, function(fl) {
            vapply(fl, function(f) f$k[["essBulk"]], 0)
          }))
          kAcross <- vapply(
            st,
            function(s) {
              i <- match(L, s$lengths)
              if (is.na(i)) NA_real_ else s$across[[i]]$k[["rhat"]]
            },
            0
          )
          rowAcross <- unlist(lapply(st, function(s) {
            i <- match(L, s$lengths)
            if (is.na(i)) NULL else s$across[[i]]$rows["rhat", ]
          }))
          bothR <- unlist(lapply(st, function(s) {
            if (is.null(s$both)) {
              return(NULL)
            }
            i <- match(L, vapply(s$both, `[[`, 0L, "L"))
            if (is.na(i)) NULL else s$both[[i]]$rhat
          }))
          mix[[length(mix) + 1L]] <- data.frame(
            regime = cell,
            start = start,
            build = which,
            L = L,
            seconds = L * sps * chainsPer,
            kRhatFitMedian = stats::median(kFit, na.rm = TRUE),
            kEssFitMedian = stats::median(kEss, na.rm = TRUE),
            kEssPerSecond = stats::median(kEss, na.rm = TRUE) /
              (L * sps * chainsPer),
            kRhatAcrossMedian = stats::median(kAcross, na.rm = TRUE),
            kRhatAcrossShareBelow = mean(kAcross < 1.01, na.rm = TRUE),
            rowRhatAcrossP90 = stats::quantile(
              rowAcross,
              0.9,
              na.rm = TRUE,
              names = FALSE
            ),
            bothStartsRhatP90 = if (length(bothR)) {
              stats::quantile(bothR, 0.9, na.rm = TRUE, names = FALSE)
            } else {
              NA_real_
            }
          )
        }
      }
    }
  }
  tab <- do.call(rbind, rows)
  if (length(rows) == 0L || length(mix) == 0L) {
    stop("no streams to read")
  }
  mixTab <- do.call(rbind, mix)
  firstBelow <- do.call(
    rbind,
    lapply(
      split(mixTab, paste(mixTab$regime, mixTab$start, mixTab$build)),
      function(g) {
        g <- g[order(g$L), ]
        ok <- g$kRhatAcrossMedian < 1.01
        stay <- rev(cumprod(rev(ok))) == 1
        i <- which(stay)[1L]
        data.frame(
          regime = g$regime[1L],
          start = g$start[1L],
          build = g$build[1L],
          firstL = if (is.na(i)) NA_integer_ else g$L[i],
          firstSeconds = if (is.na(i)) NA_real_ else g$seconds[i]
        )
      }
    )
  )
  out <- list(
    coverage = tab,
    mixing = mixTab,
    firstBelow = firstBelow,
    quick = quick
  )
  saveAtomic(out, runDir("readings", "convergence.rds"))
  print(firstBelow, row.names = FALSE)
  invisible(out)
}

## ------------------------------------------------------------ Reading 3

# The k study's stored datasets with accepted long runs (their compact pooled
# reference: 1,000 quantiles a row, k in full), on the regimes read here. A
# dataset is used only when the model the pricing build makes from its data
# equals the one the k study's build made (r3 fingerprint under each, r3 check).
kstudyRun <- function() {
  root <- Sys.getenv("KMP_KSTUDY_RUN")
  if (!nzchar(root)) {
    stop("KMP_KSTUDY_RUN must name the k study's run directory")
  }
  root
}

r3Keys <- function() {
  keys <- character()
  for (cell in convRegimes) {
    dirs <- Sys.glob(file.path(
      kstudyRun(),
      "ref",
      kmpLink,
      paste0(cell, "-r*")
    ))
    for (d in dirs) {
      acc <- file.path(d, "acceptance.rds")
      if (
        file.exists(acc) &&
          file.exists(file.path(d, "compact.rds")) &&
          isTRUE(readRDS(acc)$accepted)
      ) {
        keys <- c(keys, basename(d))
      }
    }
  }
  keys
}

r3Design <- function(key) {
  cell <- sub("-r[0-9]+$", "", key)
  rep <- as.integer(sub("^.*-r", "", key))
  readRDS(file.path(
    kstudyRun(),
    "design",
    kmpLink,
    cell,
    sprintf("r%03d.rds", rep)
  ))
}

# The model a build makes from a stored dataset: its cut points and the
# priors the sampler states.
r3Fingerprint <- function(tag) {
  out <- list()
  for (key in r3Keys()) {
    d <- blind(r3Design(key))
    fd <- fitData(d)
    s <- bart(
      y ~ .,
      data = fd$train,
      family = kmpLink,
      n.trees = nTrees,
      n.chains = 1L,
      n.threads = 1L,
      seed = 1L,
      verbose = FALSE,
      samplerOnly = TRUE
    )
    slot1 <- function(nm) {
      tryCatch(methods::slot(s$model, nm), error = function(e) NULL)
    }
    out[[key]] <- list(
      cutPoints = s$data@cutPoints,
      power = s$model@tree.prior@power,
      base = s$model@tree.prior@base,
      leafScale = slot1("leaf.scale"),
      priorScale = slot1("prior.scale"),
      leafPrior = utils::capture.output(print(slot1("leaf.prior"))),
      leafHyperprior = utils::capture.output(print(slot1("leaf.hyperprior"))),
      k = as.numeric(s$getK())
    )
  }
  saveAtomic(
    list(tag = tag, lib = libId(), fingerprints = out),
    runDir("r3", sprintf("fingerprint-%s.rds", tag))
  )
  cat("fingerprint", tag, length(out), "datasets\n")
}

# Each pricing library's fingerprints against the reference's (the k
# study's engine, built from its archived tree): a dataset is used when
# every library makes the model the k study's did.
r3Check <- function(reference, others) {
  fp <- function(tag) {
    readRDS(runDir("r3", sprintf("fingerprint-%s.rds", tag)))$fingerprints
  }
  fa <- fp(reference)
  fbs <- lapply(others, fp)
  keys <- unique(c(names(fa), unlist(lapply(fbs, names))))
  same <- vapply(
    keys,
    function(k) {
      !is.null(fa[[k]]) &&
        all(vapply(fbs, function(fb) identical(fa[[k]], fb[[k]]), NA))
    },
    NA
  )
  saveAtomic(
    list(used = keys[same], dropped = keys[!same]),
    runDir("r3", "model-check.rds")
  )
  cat(
    "model check:",
    sum(same),
    "used,",
    sum(!same),
    "dropped:",
    keys[!same],
    "\n"
  )
  TRUE
}

r3Used <- function() {
  p <- runDir("r3", "model-check.rds")
  if (!file.exists(p)) {
    stop("no Reading 3 model check")
  }
  readRDS(p)$used
}

# 20 chains of this build from the default start on a stored dataset, read as
# fits: at each length, each fit's squared error against the reference in the
# posterior mean of k and of each row's latent fit, and in the ends of the 90
# percent interval, each scaled by the reference's posterior variance, by
# class (the reference's labels).
runR3Stream <- function(key, Lmax) {
  if (!key %in% r3Used()) {
    stop(key, " is not a dataset Reading 3 uses")
  }
  d <- r3Design(key)
  d$cell <- sub("-r[0-9]+$", "", key)
  bl <- blind(d)
  dir <- file.path(kstudyRun(), "ref", kmpLink, key)
  comp <- readRDS(file.path(dir, "compact.rds"))
  acc <- readRDS(file.path(dir, "acceptance.rds"))
  refRow <- function(nm) {
    q <- readSingle(file.path(dir, comp$rows[[nm]]$file), comp$rows[[nm]]$dims)
    list(
      cols = comp$rows[[nm]]$cols,
      mean = colMeans(q),
      var = apply(q, 2L, stats::var),
      q05 = q[which.min(abs(comp$probs - 0.05)), ],
      q95 = q[which.min(abs(comp$probs - 0.95)), ]
    )
  }
  # the row classes the reference holds: fInt and fTail (latent, at their
  # columns), or, on the effect regimes, ite (every unit, each draw's
  # difference of probabilities, as classesFrom forms it)
  have <- intersect(c("fInt", "fTail", "ite"), names(comp$rows))
  # a class may hold no row on a dataset (no tail row): it is not read there
  have <- have[vapply(
    have,
    function(nm) {
      is.null(comp$rows[[nm]]$cols) || length(comp$rows[[nm]]$cols) > 0L
    },
    NA
  )]
  ref <- lapply(have, refRow)
  names(ref) <- have
  rowsOf <- function(Ld, nm, cols) {
    if (nm == "ite") {
      n <- ncol(Ld) %/% 2L
      Q <- plogis(Ld[, seq_len(n), drop = FALSE]) -
        plogis(Ld[, n + seq_len(n), drop = FALSE])
    } else {
      Q <- Ld
    }
    if (is.null(cols)) Q else Q[, cols, drop = FALSE]
  }
  kRef <- c(
    as.vector(evenRows(comp$k$shipped$k, nrow(as.matrix(comp$k$move$k)))),
    as.vector(comp$k$move$k)
  )
  per <- cellChains(d$cell)
  lengths <- convLengths()
  lengths <- lengths[lengths <= Lmax]
  res <- list()
  for (f in seq_len(convChains %/% per)) {
    chainIds <- (f - 1L) * per + seq_len(per)
    runs <- lapply(chainIds, function(ch) {
      streamChain(bl, seedOf("r3", key, kmpBuild, ch), "default", Lmax)
    })
    arr <- aperm(simplify2array(lapply(runs, `[[`, "draws")), c(3L, 1L, 2L))
    kM <- vapply(runs, `[[`, numeric(Lmax), "k")
    res[[f]] <- lapply(lengths, function(L) {
      Ld <- poolDraws(arr, seq_len(L))
      kv <- as.vector(kM[seq_len(L), , drop = FALSE])
      byClass <- lapply(have, function(nm) {
        r <- ref[[nm]]
        Q <- rowsOf(Ld, nm, r$cols)
        q <- apply(Q, 2L, stats::quantile, c(0.05, 0.95), names = FALSE)
        c(
          mean = mean((colMeans(Q) - r$mean)^2 / r$var),
          q05 = mean((q[1L, ] - r$q05)^2 / r$var),
          q95 = mean((q[2L, ] - r$q95)^2 / r$var)
        )
      })
      names(byClass) <- have
      list(
        L = L,
        k = (mean(kv) - mean(kRef))^2 / stats::var(kRef),
        rows = byClass
      )
    })
  }
  out <- list(
    key = key,
    build = kmpBuild,
    lengths = lengths,
    fits = res,
    independent100 = 0.01,
    labelsFromReference = !is.null(acc$labels),
    quick = quick
  )
  saveAtomic(out, runDir("r3", "streams", sprintf("%s-%s.rds", key, kmpBuild)))
  cat("r3", key, kmpBuild, "\n")
}

# Reading 3: per build and length, each error (k's mean; per class the
# rows' mean and interval ends) averaged over a dataset's fits, then over
# the datasets used, with the move over shipped ratio of the dataset
# averages; in units of the reference's posterior variance, so E
# independent draws err about 1 / E on a mean.
readR3 <- function() {
  files <- Sys.glob(runDir("r3", "streams", "*.rds"))
  if (length(files) == 0L) {
    stop("no Reading 3 streams to read")
  }
  rows <- list()
  for (f in files) {
    x <- readChecked(f)
    for (li in seq_along(x$lengths)) {
      per <- lapply(x$fits, `[[`, li)
      cls <- unique(unlist(lapply(per, function(p) names(p$rows))))
      vals <- c(k = mean(vapply(per, `[[`, 0, "k")))
      for (cl in cls) {
        m <- rowMeans(vapply(per, function(p) p$rows[[cl]], numeric(3L)))
        vals[paste(cl, names(m), sep = ".")] <- m
      }
      rows[[length(rows) + 1L]] <- data.frame(
        key = x$key,
        regime = sub("-r[0-9]+$", "", x$key),
        build = x$build,
        L = x$lengths[li],
        estimand = names(vals),
        error = unname(vals)
      )
    }
  }
  tab <- do.call(rbind, rows)
  wide <- function(agg, id) {
    w <- stats::reshape(
      agg,
      idvar = id,
      timevar = "build",
      direction = "wide"
    )
    if (all(c("error.move", "error.shipped") %in% names(w))) {
      w$ratio <- w$error.move / w$error.shipped
    }
    w[do.call(order, unname(as.list(w[rev(id)]))), , drop = FALSE]
  }
  # per regime: the average over the regime's datasets; overall: the
  # average of the regime averages, so no regime weighs by its count
  byRegime <- wide(
    stats::aggregate(error ~ build + regime + L + estimand, tab, mean),
    c("regime", "L", "estimand")
  )
  overallAgg <- stats::aggregate(
    error ~ build + L + estimand,
    stats::aggregate(error ~ build + regime + L + estimand, tab, mean),
    mean
  )
  w <- wide(overallAgg, c("L", "estimand"))
  out <- list(
    perDataset = tab,
    byRegime = byRegime,
    table = w,
    used = tryCatch(r3Used(), error = function(e) NULL),
    quick = quick
  )
  saveAtomic(out, runDir("readings", "r3.rds"))
  print(w, row.names = FALSE, digits = 3)
  invisible(out)
}

## ------------------------------------------------------------------ checks

# A job file's form: seven tab-separated fields, unique ids, every
# dependency an earlier job, every command one this file dispatches.
checkJobs <- function(path) {
  lines <- readLines(path, warn = FALSE)
  f <- strsplit(lines, "\t", fixed = TRUE)
  ok <- all(lengths(f) == 7L)
  ids <- vapply(f, `[`, "", 1L)
  ok <- ok && !anyDuplicated(ids) && !any(grepl("[ ,~]", ids))
  seen <- character()
  for (x in f) {
    deps <- if (identical(x[2L], "-")) {
      character()
    } else {
      sub("^~", "", strsplit(x[2L], ",")[[1L]])
    }
    if (!all(deps %in% seen)) {
      cat("job", x[1L], "waits for a later or missing job\n")
      ok <- FALSE
    }
    if (!grepl("^@(shipped|move|mutant1|mutant3|kstudy) ", x[6L])) {
      cat("job", x[1L], "names no build\n")
      ok <- FALSE
    }
    act <- strsplit(sub("^@[a-z0-9]+ ", "", x[6L]), " ", fixed = TRUE)[[1L]][1L]
    if (!act %in% mainActions) {
      cat("job", x[1L], "runs no action of this file:", act, "\n")
      ok <- FALSE
    }
    seen <- c(seen, x[1L])
  }
  cat(
    "jobs",
    path,
    ":",
    length(lines),
    "jobs,",
    if (ok) "ok" else "FAILED",
    "\n"
  )
  ok
}

# The scorers on a posterior known exactly (Verification). The guard: N(0, 1)
# draws at every row stand in for the reference, for each of its halves and
# for both samplers' fits, so a fit of E independent draws errs about 100 / E
# units on a mean, D is about 0 against the pooled reference and against
# either half, and a fit of 100 independent draws errs about 1. The
# convergence measure: truths drawn from a N(0, 1) prior, one observation
# each at unit noise, chains of independent draws from the N(y / 2, 1 / 2)
# posterior, read at several lengths: coverage of the 90 percent interval
# about 0.90, R-hat about 1, bulk effective draws about the draws held, and
# the posterior mean's squared error about that of as many independent
# draws.
# runStream fed exact draws: a stub streamChain returns, for every column,
# independent normal draws about a center that is the truth plus one normal
# error a replication (shared by its chains), so each latent row's interval
# covers its truth with probability 0.9 (n - 1) / (n + 1) at n pooled draws,
# R-hat is about 1 and the posterior mean's error has sd s sqrt(1 + 1 / n).
# Run on one stored replication of a regime (its design and step 2 labels),
# its outputs written to a directory of their own.
checkStreamExact <- function(cell = "S1", reps = 30L, s = 0.5) {
  design <- loadDesign(cell, kmpLink, 1L)
  t <- design$truth
  truth <- if (design$task %in% c("effect", "s9")) {
    c(t$f$t1, t$f$t0)
  } else {
    c(t$f$train, t$f$held)
  }
  center <- truth
  stub <- function(design, seed, start, L) {
    withSeed(seed, {
      z <- matrix(rnorm(L * length(center), sd = s), L, length(center))
      list(
        draws = sweep(z, 2L, center, `+`),
        k = rnorm(L),
        seconds = 0
      )
    })
  }
  dir <- runDir("checks", "stream-exact")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  env <- environment(runStream)
  saved <- mget(c("streamChain", "streamPath", "diagPath"), envir = env)
  on.exit(
    for (nm in names(saved)) {
      assign(nm, saved[[nm]], envir = env)
    },
    add = TRUE
  )
  r <- 0L
  assign("streamChain", stub, envir = env)
  assign(
    "streamPath",
    function(...) file.path(dir, sprintf("r%03d.rds", r)),
    envir = env
  )
  assign(
    "diagPath",
    function(...) file.path(dir, sprintf("r%03d.f32", r)),
    envir = env
  )
  Lmax <- max(convLengths()[convLengths() <= tenth(2000L)])
  out <- vector("list", reps)
  for (r in seq_len(reps)) {
    center <- truth +
      withSeed(seedOf("stream-exact", cell, r), rnorm(length(truth), sd = s))
    out[[r]] <- runStream(cell, 1L, "default", Lmax)
  }
  per <- out[[1L]]$chainsPerFit
  lengths <- out[[1L]]$lengths
  ok <- TRUE
  for (li in seq_along(lengths)) {
    L <- lengths[li]
    n <- L * per
    byRep <- vapply(
      out,
      function(x) {
        v <- lapply(x$fits, function(f) f[[li]]$coverage$fInt)
        c(
          cov = sum(vapply(v, `[[`, 0, "cov90")) / sum(vapply(v, `[[`, 0, "n")),
          sq = sum(vapply(v, `[[`, 0, "sq")) / sum(vapply(v, `[[`, 0, "n")),
          rhat = x$across[[li]]$k[["rhat"]]
        )
      },
      numeric(3L)
    )
    expect <- 0.9 * (n - 1) / (n + 1)
    covSe <- stats::sd(byRep["cov", ]) / sqrt(reps)
    rmse <- sqrt(mean(byRep["sq", ]))
    rmseExpect <- s * sqrt(1 + 1 / n)
    rhat <- stats::median(byRep["rhat", ])
    good <- abs(mean(byRep["cov", ]) - expect) <= 4 * covSe + 1e-12 &&
      abs(rmse / rmseExpect - 1) < 0.1 &&
      (L < 50L || rhat < 1.02)
    cat(sprintf(
      "stream exact L %4d: coverage %.4f (exact %.4f, se %.4f), rmse %.3f (%.3f), k R-hat %.3f %s\n",
      L,
      mean(byRep["cov", ]),
      expect,
      covSe,
      rmse,
      rmseExpect,
      rhat,
      if (good) "ok" else "FAIL"
    ))
    ok <- ok && good
  }
  cat("stream exact:", if (ok) "ok" else "FAILED", "\n")
  if (ok) {
    writeLines(format(Sys.time()), file.path(dir, "passed"))
  }
  ok
}

# That the one-pass stream scorer gave summaries identical to the scorer as
# first written at every length of every fit the smoke ran (runStream
# compares them under KMP_SMOKE).
checkScoring <- function() {
  files <- list.files(
    runDir("conv"),
    "[.]rds$",
    full.names = TRUE,
    recursive = TRUE
  )
  same <- unlist(lapply(files, function(f) readRDS(f)$scoringSame))
  ok <- length(files) > 0L && length(same) > 0L && all(same)
  cat(sprintf(
    "scoring: %d stream files, %d fit lengths compared, %d identical: %s\n",
    length(files),
    length(same),
    sum(same),
    if (ok) "ok" else "FAIL"
  ))
  ok
}

# That every output carrying a library identity names one commit and one
# shared object per build, every identity is stamped, and no output was made
# with a fixed scale (KMN_FIXED_R).
checkLibs <- function() {
  files <- list.files(runDir(), "[.]rds$", full.names = TRUE, recursive = TRUE)
  ids <- list()
  for (f in files) {
    o <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.list(o) && !is.null(o$libId)) {
      ids[[f]] <- o$libId
    }
  }
  if (length(ids) == 0L) {
    cat("libs: no output carries a library identity: FAIL\n")
    return(FALSE)
  }
  tab <- do.call(rbind, lapply(ids, function(x) as.data.frame(as.list(x))))
  ok <- TRUE
  for (b in unique(tab$build)) {
    t <- tab[tab$build == b, , drop = FALSE]
    one <- length(unique(t$commit)) == 1L && length(unique(t$so)) == 1L
    stamped <- all(!is.na(t$stamped) & t$stamped == b)
    fixed <- any(nzchar(t$fixedR))
    cat(sprintf(
      "libs: %s: %d outputs, commit %s, so %s, %s%s\n",
      b,
      nrow(t),
      paste(unique(substr(t$commit, 1L, 9L)), collapse = "/"),
      paste(unique(substr(t$so, 1L, 8L)), collapse = "/"),
      if (stamped) "stamped" else "UNSTAMPED",
      if (fixed) ", FIXED R" else ""
    ))
    ok <- ok && one && stamped && !fixed
  }
  cat(sprintf("libs: %s\n", if (ok) "ok" else "FAIL"))
  ok
}

checkConjugate <- function() {
  withSeed(seedOf("conjugate"), {
    rows <- 30L
    half <- 10000L
    refS <- matrix(rnorm(half * rows), ncol = rows)
    refM <- matrix(rnorm(half * rows), ncol = rows)
    ref <- rbind(refS, refM)
    values <- referenceValues(ref)
    units <- referenceUnits(ref, values, "conjugate")
    labels <- rep(c(TRUE, FALSE), length.out = rows)
    fitOf <- function(E) {
      fitEstimands(matrix(rnorm(E * rows), ncol = rows), rep(1, E), rows)
    }
    E <- 400L
    fits <- list(
      shipped = replicate(20L, fitOf(E), simplify = FALSE),
      move = replicate(20L, fitOf(E), simplify = FALSE)
    )
    errOn <- function(v) {
      lapply(fits, fitErrors, values = v, units = units, labels = labels)
    }
    errs <- errOn(values)
    errS <- errOn(referenceValues(refS))
    errM <- errOn(referenceValues(refM))
    meanErr <- mean(errs$move[, "all.mean"])
    w <- welch(errs$move[, "all.mean"], errs$shipped[, "all.mean"])
    wS <- welch(errS$move[, "all.mean"], errS$shipped[, "all.mean"])
    wM <- welch(errM$move[, "all.mean"], errM$shipped[, "all.mean"])
    err100 <- mean(fitErrors(
      replicate(40L, fitOf(guardIndependent), simplify = FALSE),
      values,
      units,
      labels
    )[, "all.mean"])
    # the convergence measure's scorers
    truths <- 400L
    conv <- t(vapply(
      c(25L, 100L, 1000L),
      function(L) {
        theta <- rnorm(truths)
        y <- theta + rnorm(truths)
        cover <- sq <- rh <- es <- numeric(truths)
        for (i in seq_len(truths)) {
          m <- matrix(rnorm(L * 4L, y[i] / 2, sqrt(0.5)), L, 4L)
          q <- stats::quantile(m, c(0.05, 0.95), names = FALSE)
          cover[i] <- theta[i] >= q[1L] && theta[i] <= q[2L]
          sq[i] <- (mean(m) - y[i] / 2)^2 / 0.5
          if (i <= 40L) {
            mx <- mixingAt(m)
            rh[i] <- mx[["rhat"]]
            es[i] <- mx[["essBulk"]]
          }
        }
        c(
          L = L,
          coverage = mean(cover),
          errRatio = mean(sq) * 4L * L,
          rhat = mean(rh[1:40]),
          essShare = mean(es[1:40]) / (4L * L)
        )
      },
      numeric(5L)
    ))
  })
  checks <- c(
    unitMean = abs(meanErr - 100 / E) < 0.5 * 100 / E,
    welchHarmless = w[["hi"]] < guardMargin,
    dNearZero = abs(w[["D"]]) < 3 * w[["se"]],
    dShippedHalf = abs(wS[["D"]]) < 3 * wS[["se"]],
    dMoveHalf = abs(wM[["D"]]) < 3 * wM[["se"]],
    error100 = abs(err100 - 1) < 0.3,
    coverage = all(abs(conv[, "coverage"] - 0.9) < 3 * sqrt(0.09 / truths)),
    errIndependent = all(abs(conv[, "errRatio"] - 1) < 0.3),
    rhat = all(abs(conv[, "rhat"] - 1) < 0.05),
    ess = all(conv[, "essShare"] > 0.7 & conv[, "essShare"] < 1.4)
  )
  print(c(
    meanErr = meanErr,
    expected = 100 / E,
    D = w[["D"]],
    hi = w[["hi"]],
    DShippedHalf = wS[["D"]],
    DMoveHalf = wM[["D"]],
    err100 = err100
  ))
  print(conv)
  print(checks)
  cat("conjugate:", if (all(checks)) "ok" else "FAILED", "\n")
  all(checks)
}

# Equal time as run: every step 2 shipped fit and guard shipped fit holds
# the draws r from the timing record gives (Verification).
# Whether an output was made with KMN_FIXED_R set (or carries no identity
# to say it was not).
fixedRUsed <- function(x) {
  is.null(x$libId) || nzchar(x$libId[["fixedR"]])
}

checkEqualTime <- function() {
  ok <- TRUE
  files <- Sys.glob(runDir("p2fits", "*", "r*-shipped.rds"))
  for (f in files) {
    x <- readChecked(f)
    if (x$drawsPerChain != equalTimeDraws(x$cell) || fixedRUsed(x)) {
      cat("unequal time or fixed r:", f, "\n")
      ok <- FALSE
    }
  }
  gfiles <- Sys.glob(runDir("guard", "*", "shipped-fit*.rds"))
  for (f in grep("-grid", gfiles, invert = TRUE, value = TRUE)) {
    x <- readChecked(f)
    if (x$drawsPerChain != equalTimeDraws(guardCell(x$name)) || fixedRUsed(x)) {
      cat("unequal time or fixed r:", f, "\n")
      ok <- FALSE
    }
  }
  cat(
    "equal time:",
    length(files) + length(gfiles),
    "fits,",
    if (ok) "ok" else "FAILED",
    "\n"
  )
  ok
}

# A job file's outputs: every fit, score, stream and effective-draws job's
# file present, and no two jobs writing one path (Verification: every
# replication of every arm present, no reading written to another's path).
checkOutputs <- function(path) {
  f <- strsplit(readLines(path, warn = FALSE), "\t", fixed = TRUE)
  outs <- character()
  for (x in f) {
    build <- sub("^@([a-z0-9]+) .*$", "\\1", x[6L])
    a <- strsplit(sub("^@[a-z0-9]+ ", "", x[6L]), " ", fixed = TRUE)[[1L]]
    out <- switch(
      a[1L],
      fit = step2FitPath(a[2L], as.integer(a[3L]), a[4L]),
      score = step2ScorePath(a[2L], as.integer(a[3L])),
      stream = streamPath(a[2L], as.integer(a[3L]), a[4L], build),
      ess = runDir(
        "ess",
        cellFile(a[2L]),
        sprintf("%s-%s.rds", repTag(as.integer(a[3L])), build)
      ),
      NULL
    )
    if (!is.null(out)) {
      outs <- c(outs, out)
    }
  }
  missing <- outs[!file.exists(outs)]
  dup <- unique(outs[duplicated(outs)])
  if (length(missing)) {
    cat("missing:", head(missing, 20L), sep = "\n  ")
  }
  if (length(dup)) {
    cat("written twice:", dup, sep = "\n  ")
  }
  ok <- length(missing) == 0L && length(dup) == 0L
  cat(
    "outputs:",
    length(outs),
    "expected,",
    length(missing),
    "missing,",
    if (ok) "ok" else "FAILED",
    "\n"
  )
  ok
}

# The functions that make arms, labels, fits, streams and guard readings
# never read a truth, and no function but the loader, the design and the
# held-out score reads a sealed outcome; then the k study's own check.
checkPricingBlind <- function() {
  fns <- c(
    "runStep2Arm",
    "pricingLabels",
    "streamChain",
    "diagIndex",
    "runGuardFits",
    "runGuardPilot",
    "guardLabels",
    "referenceValues",
    "referenceUnits",
    "fitEstimands",
    "fitErrors",
    "guardDatasetReading",
    "triggerFeatures",
    "runTimingFit",
    "runEss"
  )
  txt <- function(f) paste(deparse(get(f)), collapse = "\n")
  noTruth <- !vapply(fns, function(f) grepl("truth|designed", txt(f)), NA)
  sealedReaders <- c(
    "passerRows",
    "runGuardDesign",
    "runHeldout",
    "runLoad",
    "runSealedScores"
  )
  others <- setdiff(
    ls(globalenv())[vapply(
      ls(globalenv()),
      function(f) is.function(get(f)),
      NA
    )],
    c(sealedReaders, "sealedPath", "checkPricingBlind")
  )
  noSealed <- !vapply(
    others,
    function(f) grepl("\"outcomes\"|\"heldout\"", txt(f)),
    NA
  )
  if (!all(noTruth)) {
    cat("reads a truth:", fns[!noTruth], "\n")
  }
  if (!all(noSealed)) {
    cat("reads a sealed outcome:", others[!noSealed], "\n")
  }
  ok <- all(noTruth) && all(noSealed) && isTRUE(checkBlind())
  cat("blind:", if (ok) "ok" else "FAILED", "\n")
  ok
}

## -------------------------------------------------------------------- jobs

kmpJobId <- function(...) jobId(...)

# Job files in the k study queue's form: id, ids waited for, peak memory
# (GiB), cell, NA, the command and the gate files (relative to KMP_RUN). A
# command names its build (@shipped, @move, @mutant1, @mutant3). Phases:
# smoke (every action at quick sizes on a few replications and three
# stand-in datasets), timing, step2:<look>, ess, step1:pilot, step1:run,
# conv, r3, guard:prefit, guard:refs, guard:fits, guard:read, step3.
kmpWriteJobs <- function(phase, path) {
  jobs <- list()
  written <- character()
  # a dependency on a job of an earlier phase's file is met by the phases'
  # order (each written after the last is read), so only this file's are kept
  add <- function(id, deps, cell, cmd, gate = character(), mem = 2) {
    deps <- deps[sub("^~", "", deps) %in% written]
    written <<- c(written, id)
    jobs[[length(jobs) + 1L]] <<- paste(
      id,
      if (length(deps)) paste(deps, collapse = ",") else "-",
      mem,
      cell,
      "NA",
      cmd,
      if (length(gate)) paste(gate, collapse = ",") else "-",
      sep = "\t"
    )
    id
  }
  smoke <- phase == "smoke"
  repsOf <- function(id, look) if (smoke) 2L else lookCount(id, look)
  designs <- function(ids, reps) {
    for (id in ids) {
      for (rep in seq_len(reps(id))) {
        add(
          kmpJobId("design", id, rep),
          character(),
          id,
          sprintf("@shipped design %s %d", id, rep)
        )
      }
    }
  }
  timingJobs <- function(ids, rounds, reps) {
    last <- character()
    for (t in seq_len(rounds)) {
      for (id in ids) {
        for (rep in seq_len(reps)) {
          for (which in if (t %% 2L == 1L) {
            c("shipped", "move")
          } else {
            c("move", "shipped")
          }) {
            last <- add(
              kmpJobId("timing", id, rep, which, t),
              c(kmpJobId("design", id, rep), last),
              id,
              sprintf("@%s timing fit %s %d %d", which, id, rep, t)
            )
          }
        }
      }
    }
    last
  }
  timedJobs <- function(cells, rounds) {
    last <- character()
    for (cell in cells) {
      last <- add(
        kmpJobId("state", cell),
        last,
        cell,
        sprintf("@shipped timing state %s", cell)
      )
    }
    for (t in seq_len(rounds)) {
      for (cell in cells) {
        for (which in if (t %% 2L == 1L) {
          c("shipped", "move")
        } else {
          c("move", "shipped")
        }) {
          last <- add(
            kmpJobId("segment", cell, which, t),
            last,
            cell,
            sprintf("@%s timing segment %s %d", which, cell, t)
          )
        }
      }
    }
    last
  }
  step2Jobs <- function(look, after) {
    for (id in step2Regimes) {
      from <- if (look == 1L || smoke) 1L else lookCount(id, 1L) + 1L
      for (rep in seq.int(from, repsOf(id, look))) {
        if (!smoke) {
          add(
            kmpJobId("design", id, rep),
            character(),
            id,
            sprintf("@shipped design %s %d", id, rep)
          )
        }
        fs <- add(
          kmpJobId("fit", id, rep, "shipped"),
          c(kmpJobId("design", id, rep), after),
          id,
          sprintf("@shipped fit %s %d shipped", id, rep)
        )
        fm <- add(
          kmpJobId("fit", id, rep, "move"),
          c(kmpJobId("design", id, rep), after),
          id,
          sprintf("@move fit %s %d move", id, rep)
        )
        add(
          kmpJobId("score", id, rep),
          c(fs, fm),
          id,
          sprintf("@shipped score %s %d", id, rep)
        )
      }
    }
    add(
      kmpJobId("read", "step2", look),
      paste0(
        "~",
        unlist(lapply(step2Regimes, function(id) {
          vapply(
            seq_len(repsOf(id, look)),
            function(r) kmpJobId("score", id, r),
            ""
          )
        }))
      ),
      "read",
      sprintf("@shipped read step2 %d", look)
    )
  }
  essJobs <- function() {
    ids <- character()
    for (id in essRegimes) {
      for (rep in seq_len(if (smoke) 1L else 10L)) {
        for (which in c("shipped", "move")) {
          ids <- c(
            ids,
            add(
              kmpJobId("ess", id, rep, which),
              kmpJobId("design", id, rep),
              id,
              sprintf("@%s ess %s %d", which, id, rep)
            )
          )
        }
      }
    }
    add("read-ess", paste0("~", ids), "read", "@shipped read ess")
  }
  # step 1: the pilot (the move, every prior); the power run (mutants 1
  # and 3, the tilted prior, 50 a regime) read before the main run; the main
  # run (200 a regime, the contingent arm gated on the pilot's trigger), its
  # early look at 150 a prior and its reading
  step1Jobs <- function(phase1) {
    ids <- character()
    sbc <- function(
      b,
      ph,
      prior,
      deps,
      gate = character(),
      reps = seq_len(step1Reps(ph))
    ) {
      for (id in step1Regimes) {
        for (rep in reps) {
          ids <<- c(
            ids,
            add(
              kmpJobId("sbc", ph, b, prior, id, rep),
              deps,
              id,
              sprintf(
                "@%s sbc %s %s %s %d",
                b,
                if (ph == "power") "run" else ph,
                prior,
                id,
                rep
              ),
              gate
            )
          )
        }
      }
    }
    if (phase1 == "pilot") {
      for (prior in names(step1Priors)) {
        sbc("move", "pilot", prior, character())
      }
      add("read-step1pilot", paste0("~", ids), "read", "@move read step1pilot")
    } else if (phase1 == "power") {
      for (b in c("mutant1", "mutant3")) {
        ids <- character()
        sbc(b, "power", "tilted", if (smoke) "read-step1pilot" else character())
        add(
          kmpJobId("read", "step1", b),
          paste0("~", ids),
          "read",
          sprintf("@%s read step1 %s", b, b)
        )
      }
    } else {
      deps <- if (smoke) {
        c("read-step1pilot", "~read-step1-mutant1", "~read-step1-mutant3")
      } else {
        character()
      }
      if (!smoke) {
        for (b in c("mutant1", "mutant3")) {
          f <- runDir("readings", sprintf("step1-%s.rds", b))
          if (!file.exists(f)) {
            stop(
              "step 1's power run (",
              b,
              ") is read before its main run is written"
            )
          }
          if (!isTRUE(readChecked(f)$flagged)) {
            warning(
              b,
              " was not flagged: step 1 is restated before its main run"
            )
          }
        }
      }
      # the early look reads the first 50 a regime (150 a prior) once they
      # are done; the rest run only when it writes early-go (no p below
      # alpha / 10, every draw finite)
      early <- step1Reps("power")
      first <- seq_len(early)
      rest <- seq.int(early + 1L, step1Reps("run"))
      for (prior in c("default", "tilted")) {
        sbc("move", "run", prior, deps, reps = first)
      }
      sbc("move", "run", "third", deps, "step1/third-go", reps = first)
      firstIds <- ids
      el <- add(
        "read-step1-early",
        paste0("~", firstIds),
        "read",
        sprintf("@move read step1 move %d", early)
      )
      for (prior in c("default", "tilted")) {
        sbc("move", "run", prior, el, "step1/early-go", reps = rest)
      }
      sbc(
        "move",
        "run",
        "third",
        el,
        c("step1/third-go", "step1/early-go"),
        reps = rest
      )
      add(
        "read-step1",
        c(paste0("~", ids), "~read-step1-early"),
        "read",
        "@move read step1"
      )
    }
  }
  # the convergence measure: the first look (50 a regime, the dispersed
  # start at 25, to 2,000 sweeps); the second look (to 100 in a regime the
  # first look's standard errors send on, both starts' rule); the long
  # streams (rerun from the start to 8,000, same seeds, where streamsRunOn
  # says so for a library, its dispersed streams lengthened with it)
  convJobs <- function(look = "first") {
    L2 <- max(convLengths()[convLengths() <= tenth(2000L)])
    L8 <- max(convLengths())
    first <- if (smoke) 1L else 50L
    ids <- character()
    stream <- function(id, rep, start, which, Lmax, tag) {
      ids <<- c(
        ids,
        add(
          kmpJobId(tag, id, rep, start, which),
          c(kmpJobId("score", id, rep), if (smoke) "check-streams"),
          id,
          sprintf("@%s stream %s %d %s %d", which, id, rep, start, Lmax)
        )
      )
    }
    dispersed <- min(first, convDispersedReps)
    gate <- runDir("checks", "stream-exact", "passed")
    if (!smoke && look == "first" && !file.exists(gate)) {
      stop("check streams passes before the convergence phase is written")
    }
    if (look == "first") {
      for (id in convRegimes) {
        for (rep in seq_len(first)) {
          for (start in c("default", if (rep <= dispersed) "dispersed")) {
            for (which in c("shipped", "move")) {
              stream(id, rep, start, which, L2, "stream")
            }
          }
        }
      }
    } else if (look == "more") {
      for (id in convMoreRegimes()) {
        for (rep in seq.int(first + 1L, 2L * first)) {
          for (which in c("shipped", "move")) {
            stream(id, rep, "default", which, L2, "stream")
          }
        }
      }
    } else {
      for (id in convRegimes) {
        for (which in c("shipped", "move")) {
          if (!streamsRunOn(id, which, first)) {
            next
          }
          reps <- sum(file.exists(vapply(
            seq_len(2L * first),
            function(r) streamPath(id, r, "default", which),
            ""
          )))
          for (rep in seq_len(reps)) {
            for (start in c("default", if (rep <= dispersed) "dispersed")) {
              stream(id, rep, start, which, L8, "stream8k")
            }
          }
        }
      }
    }
    if (length(ids) == 0L) {
      return(invisible(NULL))
    }
    add(
      paste0("read-conv-", look),
      paste0("~", ids),
      "read",
      "@shipped read convergence"
    )
  }
  r3Jobs <- function(keys) {
    r0 <- add(
      "r3-fp-kstudy",
      character(),
      "r3",
      "@kstudy r3 fingerprint kstudy"
    )
    a <- add("r3-fp-shipped", r0, "r3", "@shipped r3 fingerprint shipped")
    b <- add("r3-fp-move", a, "r3", "@move r3 fingerprint move")
    ck <- add("r3-check", b, "r3", "@shipped r3 check kstudy shipped move")
    ids <- character()
    for (key in keys) {
      for (which in c("shipped", "move")) {
        ids <- c(
          ids,
          add(
            kmpJobId("r3", key, which),
            ck,
            "r3",
            sprintf("@%s r3 stream %s %d", which, key, tenth(2000L))
          )
        )
      }
    }
    add("read-r3", paste0("~", ids), "read", "@shipped read r3")
  }
  guardRefJobs <- function(name, after) {
    cell <- guardCell(name)
    d <- add(
      kmpJobId("gdesign", cell),
      after,
      cell,
      sprintf("@shipped guard design %s", shQuoteName(name))
    )
    rounds <- if (smoke) 2L else 9L
    last <- d
    for (t in seq_len(rounds)) {
      for (which in if (t %% 2L == 1L) {
        c("shipped", "move")
      } else {
        c("move", "shipped")
      }) {
        last <- add(
          kmpJobId("gtiming", cell, which, t),
          last,
          cell,
          sprintf("@%s guard timing %s %d", which, shQuoteName(name), t)
        )
      }
    }
    prevAcc <- NULL
    parts <- seq_len(refPartsPerKernel(cell))
    for (round in refRoundsFor(cell)) {
      L <- refLengths()[round + 1L]
      halves <- if (L <= mcmc$refBlockMax) c("shipped", "move") else "shipped"
      go <- if (round > 0L) {
        refGateRel(cell, kmpLink, 1L, sprintf("go-%d", round))
      } else {
        character()
      }
      partIds <- character()
      for (which in halves) {
        for (part in parts) {
          partIds <- c(
            partIds,
            if (round == 0L) {
              add(
                kmpJobId("ref", cell, which, part),
                d,
                cell,
                sprintf(
                  "@%s reference run %s %s 1 %d",
                  which,
                  cell,
                  kmpLink,
                  part
                )
              )
            } else {
              add(
                kmpJobId("ref", cell, which, part, paste0("r", round)),
                prevAcc,
                cell,
                sprintf(
                  "@%s reference extend %s %s 1 %d %d",
                  which,
                  cell,
                  kmpLink,
                  part,
                  round
                ),
                go
              )
            }
          )
        }
      }
      prevAcc <- add(
        kmpJobId("accept", cell, paste0("r", round)),
        partIds,
        cell,
        sprintf("@shipped reference accept %s %s 1 %d", cell, kmpLink, round),
        go
      )
    }
    c(last = last, accept = prevAcc)
  }
  guardFitJobs <- function(name, after, step3 = FALSE) {
    cell <- guardCell(name)
    n0 <- if (quick) {
      4L
    } else if (step3) {
      step3Fits
    } else {
      guardFitsStart
    }
    ids <- character()
    for (which in c("shipped", "move")) {
      ids <- c(
        ids,
        add(
          kmpJobId("gfits", cell, which),
          after,
          cell,
          sprintf("@%s guard fits %s 1 %d", which, shQuoteName(name), n0)
        )
      )
    }
    # step 3's fits: 20 a library, no pilot raise (20 more only to confirm
    # a miss, written by hand with a tag)
    if (step3) {
      return(ids)
    }
    p <- add(
      kmpJobId("gpilot", cell),
      ids,
      cell,
      sprintf("@shipped guard pilot %s", shQuoteName(name))
    )
    for (which in c("shipped", "move")) {
      ids <- c(
        ids,
        add(
          kmpJobId("gmore", cell, which),
          p,
          cell,
          sprintf("@%s guard more %s", which, shQuoteName(name))
        )
      )
    }
    ids
  }
  guardNames <- function(set) {
    a <- tryCatch(assignment(), error = function(e) NULL)
    if (is.null(a)) {
      stop("the assignment is drawn before the guard's jobs are written")
    }
    a[[set]]
  }
  switch(
    phase,
    smoke = {
      designs(step2Regimes, function(id) 2L)
      tm <- timingJobs(step2Regimes, 2L, 1L)
      td <- timedJobs(
        c("large:10000:75:sep", "hazard:2000:75", "counts:2000:75"),
        2L
      )
      sm <- add(
        "timing-summarize",
        c(tm, td),
        "timing",
        "@shipped timing summarize"
      )
      step2Jobs(1L, sm)
      essJobs()
      step1Jobs("pilot")
      step1Jobs("power")
      step1Jobs("run")
      add(
        "check-streams",
        kmpJobId("score", "S1", 1L),
        "S1",
        "@move check streams"
      )
      convJobs()
      r3Jobs(head(tryCatch(r3Keys(), error = function(e) character()), 1L))
      # the guard on three stand-in datasets, never a passer
      sd <- add("standins", character(), "guard", "@shipped guard standins")
      loads <- vapply(
        passers$name,
        function(nm) {
          add(
            kmpJobId("load", fileTag(nm)),
            sd,
            "guard",
            sprintf("@shipped guard load %s", shQuoteName(nm))
          )
        },
        ""
      )
      sc <- add("gscores", loads, "guard", "@shipped guard scores")
      as <- add("gassign", sc, "guard", "@shipped guard assign")
      f3 <- add(
        "step3-features",
        c(as, kmpJobId("read", "step2", 1L)),
        "step3",
        "@shipped step3 features"
      )
      sel <- add("step3-select", f3, "step3", "@shipped step3 select 1")
      acc <- character()
      for (nm in passers$name) {
        r <- guardRefJobs(nm, sel)
        tsum <- add(
          kmpJobId("gtsum", fileTag(nm)),
          r[["last"]],
          "timing",
          "@shipped timing summarize"
        )
        acc <- c(acc, guardFitJobs(nm, c(tsum, r[["accept"]])))
      }
      gr <- add("guard-read", paste0("~", acc), "read", "@shipped guard read")
      add("guard-read-no995", gr, "read", "@shipped guard read no995")
      ho <- add("held-out-score", gr, "read", "@shipped guard heldout")
      add("step3-read", ho, "read", "@shipped step3 read")
    },
    timing = {
      designs(step2Regimes, function(id) 3L)
      tm <- timingJobs(step2Regimes, 3L, 3L)
      td <- timedJobs(timedCells, 10L)
      add("timing-summarize", c(tm, td), "timing", "@shipped timing summarize")
    },
    "step2:1" = step2Jobs(1L, character()),
    "step2:2" = step2Jobs(2L, character()),
    ess = essJobs(),
    "step1:pilot" = step1Jobs("pilot"),
    "step1:power" = step1Jobs("power"),
    "step1:run" = step1Jobs("run"),
    conv = convJobs("first"),
    "conv:more" = convJobs("more"),
    "conv:long" = convJobs("long"),
    r3 = r3Jobs(r3Keys()),
    "guard:prefit" = {
      loads <- vapply(
        passers$name,
        function(nm) {
          add(
            kmpJobId("load", fileTag(nm)),
            character(),
            "guard",
            sprintf("@shipped guard load %s", shQuoteName(nm))
          )
        },
        ""
      )
      sc <- add("gscores", loads, "guard", "@shipped guard scores")
      add("gassign", sc, "guard", "@shipped guard assign")
    },
    "guard:refs" = {
      for (nm in guardNames("r2")) {
        guardRefJobs(nm, character())
      }
    },
    "guard:fits" = {
      if (!file.exists(runDir("step3", "commit.rds"))) {
        stop("the guard's fits are written after step 3 is committed")
      }
      add(
        "timing-summarize",
        character(),
        "timing",
        "@shipped timing summarize"
      )
      for (nm in guardNames("r2")) {
        guardFitJobs(nm, "timing-summarize")
      }
    },
    "guard:read" = {
      gr <- add("guard-read", character(), "read", "@shipped guard read")
      add("guard-read-no995", gr, "read", "@shipped guard read no995")
      add("held-out-score", gr, "read", "@shipped guard heldout")
    },
    step3 = {
      for (nm in guardNames("r3")) {
        r <- guardRefJobs(nm, character())
        guardFitJobs(nm, r[["accept"]], step3 = TRUE)
      }
    },
    stop("no phase ", phase)
  )
  lines <- as.character(unlist(jobs))
  writeLines(lines, path)
  cat(sprintf("%d jobs written to %s\n", length(lines), path))
}

# A dataset's name inside a job's command, whose fields split on spaces.
shQuoteName <- function(name) gsub(" ", "%20", name, fixed = TRUE)
unquoteName <- function(name) gsub("%20", " ", name, fixed = TRUE)

## -------------------------------------------------------------------- main

mainActions <- c(
  "design",
  "timing",
  "fit",
  "score",
  "ess",
  "sbc",
  "stream",
  "r3",
  "reference",
  "guard",
  "step3",
  "read",
  "check",
  "jobs"
)

kmpMain <- function(args) {
  if (length(args) == 0L) {
    stop("no action")
  }
  action <- args[1L]
  a <- args[-1L]
  opt <- function(i) if (length(a) >= i) a[i] else ""
  int <- function(i) as.integer(a[i])
  if (!action %in% mainActions) {
    stop("no action ", action)
  }
  switch(
    action,
    design = makeDesign(a[1L], kmpLink, int(2L)),
    timing = switch(
      a[1L],
      fit = runTimingFit(a[2L], int(3L), int(4L)),
      segment = runTimedSegment(a[2L], int(3L)),
      state = runTimedState(a[2L]),
      summarize = summarizeTimingPieces(),
      stop("no timing action ", a[1L])
    ),
    fit = runStep2Arm(a[1L], int(2L), a[3L]),
    score = runStep2Score(a[1L], int(2L)),
    ess = runEss(a[1L], int(2L)),
    sbc = runStep1(a[1L], a[2L], a[3L], int(4L)),
    stream = runStream(a[1L], int(2L), a[3L], int(4L)),
    r3 = switch(
      a[1L],
      fingerprint = r3Fingerprint(a[2L]),
      check = r3Check(a[2L], a[-(1:2)]),
      stream = runR3Stream(a[2L], int(3L)),
      stop("no r3 action ", a[1L])
    ),
    reference = switch(
      a[1L],
      run = runReferencePart(a[2L], a[3L], int(4L), int(5L)),
      extend = extendReferencePart(a[2L], a[3L], int(4L), int(5L), int(6L)),
      accept = acceptReference(a[2L], a[3L], int(4L), int(5L)),
      stop("no reference action ", a[1L])
    ),
    guard = switch(
      a[1L],
      standins = writeStandins(),
      load = runLoad(unquoteName(a[2L])),
      scores = runSealedScores(),
      assign = runAssign(),
      design = runGuardDesign(unquoteName(a[2L])),
      timing = runTimingFit(guardCell(unquoteName(a[2L])), 1L, int(3L)),
      fits = runGuardFits(unquoteName(a[2L]), int(3L), int(4L), opt(5L)),
      more = runGuardMore(unquoteName(a[2L])),
      pilot = runGuardPilot(unquoteName(a[2L])),
      read = runGuardRead(identical(opt(2L), "no995")),
      heldout = runHeldout(),
      stop("no guard action ", a[1L])
    ),
    step3 = switch(
      a[1L],
      features = runStep3Features(),
      select = runStep3Select(if (nzchar(opt(2L))) int(2L) else NA_integer_),
      read = runStep3Read(),
      fresh = runStep3Fresh(),
      stop("no step3 action ", a[1L])
    ),
    read = switch(
      a[1L],
      step1pilot = readStep1Pilot(),
      step1 = readStep1(
        if (nzchar(opt(2L))) a[2L] else "move",
        if (nzchar(opt(3L))) int(3L) else NA_integer_
      ),
      step2 = readStep2(int(2L)),
      convergence = readConvergence(),
      ess = readEss(),
      r3 = readR3(),
      stop("no reading ", a[1L])
    ),
    check = {
      ok <- switch(
        a[1L],
        jobs = checkJobs(a[2L]),
        conjugate = checkConjugate(),
        eqtime = checkEqualTime(),
        outputs = checkOutputs(a[2L]),
        blind = checkPricingBlind(),
        scoring = checkScoring(),
        streams = checkStreamExact(),
        libs = checkLibs(),
        stop("no check ", a[1L])
      )
      if (!isTRUE(ok)) {
        quit(status = 3L)
      }
    },
    jobs = kmpWriteJobs(a[1L], a[2L])
  )
  invisible(NULL)
}

if (sys.nframe() == 0L) {
  kmpMain(commandArgs(TRUE))
}
