#!/usr/bin/env Rscript

# The feel study for the monotone default
# (docs/plans/monotone-exact-birth-death.md, Default: feel study): how a user
# would meet monotone(prior = "leaf") and monotone(prior = "joint") as the
# default, beyond run time. Four arms, each fit through bart():
#   leaf        monotone(c(x1 = "increasing"), prior = "leaf"), cgm()
#   joint       monotone(c(x1 = "increasing"), prior = "joint"), cgm()
#   joint.mbart the same under cgm(power = 0.8, base = 0.25), the mBART paper's
#               tree prior
#   free        no constraint, cgm(): the reference a user compares against
# at 200, 50 and 5 trees.
#
# Prior predictive: 500 draws of f per arm and tree count through
# samplePriorPredictive, x1 constrained and x2 free, trained on a 41 x 41 grid
# with a response spanning [-0.5, 0.5], so a rise is in units of the response
# range. Each draw is read along x1 on 201 points at x2 = 0.1, 0.3, ..., 0.9:
# distinct levels, the largest jump as a share of the rise, the share of flat
# intervals, the rise, and the share of curves with no rise at all; and the
# prior sd of f at a point, averaged over the points. The free arm's curves
# are read after sorting (its monotone rearrangement: BART's draw with its
# levels sorted), and the share of its intervals that fall is kept apart. The
# share of splits on x1 comes from a second 500 draws of sampleTreesFromPrior
# on the same sampler (samplePriorPredictive exposes no trees), with the mean
# leaves per tree.
#
# Fits: x1 increasing, x2 free and x3 noise, uniform on [0, 1]^3; noise sd a
# third of the truth's sd under that law. Truths: ramp x1; step 1{x1 > 0.5};
# hinge max(0, x1 - 0.6) / 0.4; interaction x1 (1 + 2 x2) + sin(2 pi x2). n
# 200 and 2000, 4 replicates (the recorded run took 8), one chain of 500
# burn-in and 500 kept sweeps, one thread. A job is one data set (truth, n,
# replicate) at one tree count, fit by all four arms in turn, the order
# rotated by replicate, so sweep times compare within a job. Per fit, on a
# test grid of x1 on 101 points by x2 on 0.05, 0.15, ..., 0.95, with x3
# uniform per x2 row and held along x1, so every curve along x1 holds x2 and
# x3 fixed: RMSE of the posterior mean against the truth, 95% interval
# coverage and mean width; the partial dependence along x1 (the posterior mean
# averaged over the grid's rows) against the truth's, its RMSE, its visible
# steps (increments over 5% of the truth's rise), its largest fall as a share
# of that rise, and the draws' own mean number of steps (distinct increments
# along x1); varcount shares on x1, x2, x3 and splits per tree; held-out log
# predictive score on 1000 fresh points; wall and CPU seconds per sweep.
# Under "leaf" the slow-count threshold is set to -1 s through the unexported
# test hook, so every count the sampler runs is timed into the tally: the
# number of counts, and the largest count's seconds, leaves and down-sets per
# fit, and whether the warning would fire at the default 1 s threshold. The
# timing is the engine's own wall clock and runs at every threshold, so on a
# loaded machine it also counts time the thread was descheduled.
#
# Usage:
#   run out=<dir> [quick] [cores=4] [reps=4] [part=prior,fits]
#     one result file per job in <dir>; a job whose file exists is skipped, so
#     a stopped run resumes by running it again; <dir>/DONE is written last.
#     quick: 20 prior draws, trees 50 and 5, n 200, one replicate, 50 + 50
#     sweeps.
#   report out=<dir> [csv=<file>]
#     summary tables on stdout; csv writes one row per prior job and fit

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1L] else "run"
option <- function(name, default) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^", name, "="), "", hit[1L]) else default
}
quick <- "quick" %in% args
outDir <- option("out", NA_character_)
if (is.na(outDir)) {
  stop("out=<dir> is required")
}

suppressPackageStartupMessages(library(dbarts))

arms <- c("leaf", "joint", "joint.mbart", "free")
treeCounts <- if (quick) c(50L, 5L) else c(200L, 50L, 5L)
truths <- c("ramp", "step", "hinge", "interaction")
sizes <- if (quick) 200L else c(2000L, 200L)
nReps <- if (quick) 1L else as.integer(option("reps", "4"))
nBurn <- if (quick) 50L else 500L
nKept <- if (quick) 50L else 500L
nPriorDraws <- if (quick) 20L else 500L
nHeldOut <- 1000L
parts <- strsplit(option("part", "prior,fits"), ",", fixed = TRUE)[[1L]]

truthFunction <- function(truth, x1, x2) {
  switch(
    truth,
    ramp = x1,
    step = as.numeric(x1 > 0.5),
    hinge = pmax(0, x1 - 0.6) / 0.4,
    interaction = x1 * (1 + 2 * x2) + sin(2 * pi * x2)
  )
}
# the truth's sd under the uniform law, by a 1000 x 1000 midpoint rule
truthSd <- function(truth) {
  mid <- (seq_len(1000L) - 0.5) / 1000
  f <- truthFunction(truth, rep(mid, 1000L), rep(mid, each = 1000L))
  sqrt(mean(f^2) - mean(f)^2)
}

# every monotone arm names its prior; monotone() and cgm() resolve by bare
# name inside bart()
fitArm <- function(arm, data, test, n.trees, seed) {
  common <- list(
    formula = y ~ x1 + x2 + x3,
    data = data,
    test = test,
    n.trees = n.trees,
    n.samples = nKept,
    n.burn = nBurn,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE,
    keepCall = FALSE,
    seed = seed
  )
  call <- switch(
    arm,
    leaf = quote(bart(
      monotone = monotone(c(x1 = "increasing"), prior = "leaf")
    )),
    joint = quote(bart(
      monotone = monotone(c(x1 = "increasing"), prior = "joint")
    )),
    joint.mbart = quote(bart(
      monotone = monotone(c(x1 = "increasing"), prior = "joint"),
      tree.prior = cgm(power = 0.8, base = 0.25)
    )),
    free = quote(bart())
  )
  call <- as.call(c(as.list(call), common))
  eval(call)
}

# ---- prior predictive ------------------------------------------------------

priorSampler <- function(arm, n.trees) {
  side <- seq(0, 1, length.out = 41L)
  data <- data.frame(x1 = rep(side, 41L), x2 = rep(side, each = 41L))
  data$y <- rep_len(c(-0.5, 0.5), nrow(data))
  control <- dbartsControl(n.chains = 1L, n.threads = 1L, n.trees = n.trees)
  switch(
    arm,
    leaf = dbarts(
      y ~ x1 + x2,
      data = data,
      monotone = monotone(c(x1 = "increasing"), prior = "leaf"),
      control = control
    ),
    joint = dbarts(
      y ~ x1 + x2,
      data = data,
      monotone = monotone(c(x1 = "increasing"), prior = "joint"),
      control = control
    ),
    joint.mbart = dbarts(
      y ~ x1 + x2,
      data = data,
      monotone = monotone(c(x1 = "increasing"), prior = "joint"),
      tree.prior = cgm(power = 0.8, base = 0.25),
      control = control
    ),
    free = dbarts(y ~ x1 + x2, data = data, control = control)
  )
}

curveMeasures <- function(curves, sorted) {
  # curves: one row per curve along x1
  tol <- 1e-12
  falls <- if (sorted) mean(t(apply(curves, 1L, diff)) < -tol) else 0
  if (sorted) {
    curves <- t(apply(curves, 1L, sort))
  }
  d <- t(apply(curves, 1L, diff))
  rise <- rowSums(d)
  risen <- rise > tol
  c(
    levels = mean(1 + rowSums(d > tol)),
    levels.median = stats::median(1 + rowSums(d > tol)),
    max.jump.share = mean(
      apply(d[risen, , drop = FALSE], 1L, max) / rise[risen]
    ),
    flat.share = mean(d <= tol),
    rise = mean(rise),
    no.rise = mean(!risen),
    falls = falls
  )
}

runPriorJob <- function(job) {
  set.seed(job$seed)
  sampler <- priorSampler(job$arm, job$trees)
  x1 <- seq(0, 1, length.out = 201L)
  x2 <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  test <- data.frame(x1 = rep(x1, length(x2)), x2 = rep(x2, each = 201L))
  f <- samplePriorPredictive(sampler, x.test = test, n.samples = nPriorDraws)
  curves <- do.call(
    rbind,
    lapply(seq_along(x2), function(j) f[, (j - 1L) * 201L + seq_len(201L)])
  )
  measures <- curveMeasures(curves, sorted = job$arm == "free")
  splitsX1 <- 0
  splits <- 0
  leaves <- 0
  for (i in seq_len(nPriorDraws)) {
    sampler$sampleTreesFromPrior(updateState = FALSE)
    trees <- sampler$getTrees()
    splitsX1 <- splitsX1 + sum(trees$var == 1L)
    splits <- splits + sum(trees$var > 0L)
    leaves <- leaves + sum(trees$var < 0L)
  }
  c(
    measures,
    sd.f = mean(apply(f, 2L, stats::sd)),
    x1.split.share = if (splits > 0) splitsX1 / splits else NA_real_,
    leaves.per.tree = leaves / (nPriorDraws * job$trees)
  )
}

# ---- fits against truth ----------------------------------------------------

makeData <- function(truth, n, seed) {
  set.seed(seed)
  noiseSd <- truthSd(truth) / 3
  draw <- function(n) {
    x <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n))
    x$f <- truthFunction(truth, x$x1, x$x2)
    x$y <- x$f + stats::rnorm(n, 0, noiseSd)
    x
  }
  train <- draw(n)
  x1 <- seq(0, 1, length.out = 101L)
  x2 <- seq(0.05, 0.95, by = 0.1)
  grid <- data.frame(
    x1 = rep(x1, length(x2)),
    x2 = rep(x2, each = length(x1)),
    x3 = rep(runif(length(x2)), each = length(x1))
  )
  grid$f <- truthFunction(truth, grid$x1, grid$x2)
  grid$y <- NA_real_
  list(train = train, grid = grid, heldOut = draw(nHeldOut), nx1 = 101L)
}

logMeanExp <- function(x) {
  m <- max(x)
  m + log(mean(exp(x - m)))
}

fitMeasures <- function(fit, data, n.trees) {
  nGrid <- nrow(data$grid)
  draws <- fit$yhat.test
  gridDraws <- draws[, seq_len(nGrid), drop = FALSE]
  heldDraws <- draws[, nGrid + seq_len(nHeldOut), drop = FALSE]
  truth <- data$grid$f
  postMean <- colMeans(gridDraws)
  lower <- apply(gridDraws, 2L, stats::quantile, 0.025)
  upper <- apply(gridDraws, 2L, stats::quantile, 0.975)
  # partial dependence along x1: the grid's x1 varies fastest
  nx1 <- data$nx1
  pdOf <- function(v) rowMeans(matrix(v, nx1))
  pdHat <- pdOf(postMean)
  pdTrue <- pdOf(truth)
  pdRise <- pdTrue[nx1] - pdTrue[1L]
  pdSteps <- diff(pdHat)
  drawSteps <- apply(gridDraws, 1L, function(v) {
    sum(abs(diff(pdOf(v))) > 1e-12)
  })
  sigma <- fit$sigma
  yHeld <- data$heldOut$y
  score <- vapply(
    seq_along(yHeld),
    function(i) {
      logMeanExp(stats::dnorm(yHeld[i], heldDraws[, i], sigma, log = TRUE))
    },
    numeric(1L)
  )
  varcount <- fit$varcount
  shares <- colSums(varcount) / sum(varcount)
  c(
    rmse = sqrt(mean((postMean - truth)^2)),
    coverage = mean(truth >= lower & truth <= upper),
    width = mean(upper - lower),
    pd.rmse = sqrt(mean((pdHat - pdTrue)^2)),
    pd.visible.steps = sum(pdSteps > 0.05 * pdRise),
    pd.max.fall = max(0, -min(pdSteps)) / pdRise,
    pd.draw.steps = mean(drawSteps),
    share.x1 = shares[["x1"]],
    share.x2 = shares[["x2"]],
    share.x3 = shares[["x3"]],
    splits.per.tree = mean(rowSums(varcount)) / n.trees,
    log.score = mean(score)
  )
}

runFitJob <- function(job) {
  data <- makeData(job$truth, job$n, job$dataSeed)
  test <- rbind(
    data$grid[c("x1", "x2", "x3")],
    data$heldOut[c("x1", "x2", "x3")]
  )
  order <- arms[(seq_along(arms) + job$rep - 2L) %% length(arms) + 1L]
  rows <- lapply(order, function(arm) {
    tally <- NULL
    if (arm == "leaf") {
      previous <- .Call(
        dbarts:::C_dbarts_bartcore_setMonotoneCountHooks,
        -1,
        FALSE,
        NA_integer_
      )
      on.exit(
        .Call(
          dbarts:::C_dbarts_bartcore_setMonotoneCountHooks,
          previous,
          FALSE,
          NA_integer_
        )
      )
    }
    start <- proc.time()
    fit <- withCallingHandlers(
      fitArm(arm, data$train, test, job$trees, job$fitSeed),
      dbartsSlowCountWarning = function(w) {
        tally <<- w$tally
        invokeRestart("muffleWarning")
      }
    )
    used <- proc.time() - start
    if (
      arm != "free" && !identical(fit$monotone.prior, sub("[.].*", "", arm))
    ) {
      stop("arm ", arm, " fit under prior ", format(fit$monotone.prior))
    }
    sweeps <- nBurn + nKept
    counts <- if (is.null(tally)) {
      c(0, 0, 0, 0)
    } else {
      unname(tally[c(
        "counts",
        "slowest.seconds",
        "slowest.leaves",
        "slowest.down.sets"
      )])
    }
    c(
      list(arm = arm, position = match(arm, order)),
      as.list(fitMeasures(fit, data, job$trees)),
      list(
        wall.ms.per.sweep = 1000 * used[["elapsed"]] / sweeps,
        cpu.ms.per.sweep = 1000 *
          (used[["user.self"]] + used[["sys.self"]]) /
          sweeps,
        counts = counts[1L],
        max.count.seconds = counts[2L],
        max.count.leaves = counts[3L],
        max.count.down.sets = counts[4L],
        warns.at.default = arm == "leaf" && counts[2L] > 1
      )
    )
  })
  do.call(rbind, lapply(rows, as.data.frame))
}

# ---- jobs ------------------------------------------------------------------

jobList <- function() {
  jobs <- list()
  if ("prior" %in% parts) {
    for (arm in arms) {
      for (trees in treeCounts) {
        jobs[[length(jobs) + 1L]] <- list(
          kind = "prior",
          id = sprintf("prior-%s-t%03d", arm, trees),
          arm = arm,
          trees = trees,
          seed = 100L * match(arm, arms) + trees
        )
      }
    }
  }
  if ("fits" %in% parts) {
    # the heaviest first, for balance across workers
    for (n in sort(sizes, decreasing = TRUE)) {
      for (trees in treeCounts) {
        for (truth in truths) {
          for (rep in seq_len(nReps)) {
            dataSeed <- 10000L * match(truth, truths) + 100L * n %/% 100L + rep
            jobs[[length(jobs) + 1L]] <- list(
              kind = "fit",
              id = sprintf("fit-%s-n%04d-t%03d-r%d", truth, n, trees, rep),
              truth = truth,
              n = n,
              trees = trees,
              rep = rep,
              dataSeed = dataSeed,
              fitSeed = dataSeed + trees
            )
          }
        }
      }
    }
  }
  jobs
}

runJob <- function(job) {
  file <- file.path(outDir, paste0(job$id, ".rds"))
  if (file.exists(file)) {
    return(invisible(NULL))
  }
  start <- Sys.time()
  result <- tryCatch(
    if (job$kind == "prior") runPriorJob(job) else runFitJob(job),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    writeLines(
      conditionMessage(result),
      file.path(outDir, paste0(job$id, ".err"))
    )
    return(invisible(NULL))
  }
  temp <- paste0(file, ".tmp")
  saveRDS(
    list(job = job, result = result, started = start, finished = Sys.time()),
    temp
  )
  file.rename(temp, file)
  cat(format(Sys.time(), "%H:%M:%S"), job$id, "\n")
  invisible(NULL)
}

if (mode == "run") {
  dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
  unlink(file.path(outDir, "DONE"))
  cores <- as.integer(option("cores", "4"))
  jobs <- jobList()
  cat(length(jobs), "jobs,", cores, "workers\n")
  parallel::mclapply(jobs, runJob, mc.cores = cores, mc.preschedule = FALSE)
  missing <- sum(
    !file.exists(file.path(
      outDir,
      paste0(vapply(jobs, `[[`, "", "id"), ".rds")
    ))
  )
  writeLines(
    sprintf("%d of %d jobs missing", missing, length(jobs)),
    file.path(outDir, "DONE")
  )
  cat(missing, "of", length(jobs), "jobs missing\n")
  quit(status = if (missing) 1L else 0L)
}

# ---- report ----------------------------------------------------------------

if (mode == "report") {
  files <- list.files(outDir, "[.]rds$", full.names = TRUE)
  results <- lapply(files, readRDS)
  kinds <- vapply(results, function(r) r$job$kind, "")
  prior <- do.call(
    rbind,
    lapply(results[kinds == "prior"], function(r) {
      data.frame(arm = r$job$arm, trees = r$job$trees, t(r$result))
    })
  )
  fits <- do.call(
    rbind,
    lapply(results[kinds == "fit"], function(r) {
      data.frame(
        truth = r$job$truth,
        n = r$job$n,
        trees = r$job$trees,
        rep = r$job$rep,
        r$result
      )
    })
  )
  options(width = 200L, digits = 3L)
  if (!is.null(prior)) {
    prior <- prior[order(-prior$trees, match(prior$arm, arms)), ]
    cat("\nprior predictive along x1\n")
    print(prior, row.names = FALSE)
  }
  if (!is.null(fits)) {
    fits$arm <- factor(fits$arm, arms)
    measures <- c(
      "rmse",
      "coverage",
      "width",
      "pd.rmse",
      "pd.visible.steps",
      "pd.max.fall",
      "pd.draw.steps",
      "share.x1",
      "share.x2",
      "share.x3",
      "splits.per.tree",
      "log.score",
      "cpu.ms.per.sweep"
    )
    for (measure in measures) {
      cat("\n", measure, ": mean over replicates\n", sep = "")
      table <- stats::aggregate(
        fits[[measure]],
        fits[c("arm", "trees", "n", "truth")],
        mean
      )
      wide <- stats::reshape(
        table,
        idvar = c("truth", "n", "trees"),
        timevar = "arm",
        direction = "wide"
      )
      names(wide) <- sub("^x[.]", "", names(wide))
      wide <- wide[order(match(wide$truth, truths), -wide$n, -wide$trees), ]
      print(wide, row.names = FALSE)
    }
    # sweep time relative to "joint" within each job, where the arms ran back
    # to back
    job <- interaction(fits$truth, fits$n, fits$trees, fits$rep, drop = TRUE)
    jointCpu <- fits$cpu.ms.per.sweep[fits$arm == "joint"][match(
      job,
      job[fits$arm == "joint"]
    )]
    fits$cpu.vs.joint <- fits$cpu.ms.per.sweep / jointCpu
    cat("\nCPU per sweep over joint's in the same job: median\n")
    print(
      stats::aggregate(cpu.vs.joint ~ arm + trees + n, fits, stats::median),
      row.names = FALSE
    )
    # paired over replicates: each arm minus "leaf" on the same data set,
    # mean and standard error
    for (measure in c("rmse", "coverage", "width", "log.score")) {
      cat(
        "\n",
        measure,
        ": arm minus leaf, mean (se) over replicates\n",
        sep = ""
      )
      leafValue <- fits[[measure]][fits$arm == "leaf"][match(
        job,
        job[fits$arm == "leaf"]
      )]
      fits$diff <- fits[[measure]] - leafValue
      others <- fits[fits$arm != "leaf", ]
      summary <- stats::aggregate(
        diff ~ arm + trees + n + truth,
        others,
        function(d) {
          sprintf("%.4f (%.4f)", mean(d), stats::sd(d) / sqrt(length(d)))
        }
      )
      wide <- stats::reshape(
        summary,
        idvar = c("truth", "n", "trees"),
        timevar = "arm",
        direction = "wide"
      )
      wide <- wide[order(match(wide$truth, truths), -wide$n, -wide$trees), ]
      print(wide, row.names = FALSE)
    }
    fits$diff <- NULL
    leaf <- fits[fits$arm == "leaf", ]
    cat("\nleaf counts per fit: counts run, the largest by down-sets and by")
    cat(" time, and warnings at 1 s\n")
    print(
      stats::aggregate(
        cbind(
          counts,
          max.count.leaves,
          max.count.down.sets,
          max.count.seconds,
          warns.at.default
        ) ~
          trees + n,
        leaf,
        max
      ),
      row.names = FALSE
    )
  }
  csv <- option("csv", NA_character_)
  if (!is.na(csv)) {
    rows <- list()
    if (!is.null(prior)) {
      rows[[1L]] <- data.frame(part = "prior", prior)
    }
    if (!is.null(fits)) {
      fits$cpu.vs.joint <- NULL
      rows[[length(rows) + 1L]] <- data.frame(part = "fit", fits)
    }
    all <- Reduce(
      function(a, b) {
        for (name in setdiff(names(b), names(a))) {
          a[[name]] <- NA
        }
        for (name in setdiff(names(a), names(b))) {
          b[[name]] <- NA
        }
        rbind(a, b[names(a)])
      },
      rows
    )
    numeric <- vapply(all, is.double, TRUE)
    all[numeric] <- lapply(all[numeric], signif, 4L)
    utils::write.csv(all, csv, row.names = FALSE)
  }
}
