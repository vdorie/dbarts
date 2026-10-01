#!/usr/bin/env Rscript

# The tree prior under monotone(prior = "joint")
# (docs/plans/monotone-exact-birth-death.md, Tree prior under "joint"): whether
# mBART's cgm(power = 0.8, base = 0.25) replaces cgm()'s defaults under every
# monotone fit, decided by a rule fixed before the run. Three arms, each fit
# through bart():
#   joint       monotone(<constrained> = "increasing", prior = "joint"), cgm()
#   joint.mbart the same under cgm(power = 0.8, base = 0.25)
#   free        no constraint, cgm(): the reference
# at 200 and 50 trees.
#
# Designs, predictors uniform on [0, 1]^p:
#   friedman    p 10, 10 sin(pi x1 x2 / 2) + 20 (x3 - 0.5)^2 + 10 x4 + 5 x5;
#               x1, x2, x4, x5 constrained, x3 free, x6-x10 noise
#   additive    p 10, 2 x1 + 1{x2 > 0.5} + exp(2 x3) / e^2 - sin(2 pi x4);
#               x1, x2, x3 constrained, x4 free, x5-x10 noise
#   interaction p 5, x1 (1 + 4 x2) + x1 x3 + sin(2 pi x2); x1 constrained, x2
#               and x3 free, x4 and x5 noise
# Each truth is checked nondecreasing in its constrained predictors before a
# run starts. Noise sd is the truth's sd under the design's law (Monte Carlo,
# 10^6 draws) over 3 ("low") or times 1 ("high"). n 200 and 2000, 8
# replicates, one chain of 500 burn-in and 500 kept sweeps, one thread. A job
# is one data set (design, noise, n, replicate) at one tree count, fit by all
# three arms in turn, the order rotated by replicate, so sweep times compare
# within a job. Per fit, on 2000 held-out points from the design's law: RMSE
# of the posterior mean against the truth, 95% interval coverage of the truth
# and mean width, and the log predictive score of fresh noisy responses there
# (mean per point); varcount shares per predictor, summed over the
# constrained, free and noise groups, and splits per tree; wall and CPU
# seconds per sweep. As a check of the constraint, each constrained predictor
# also gets 3 lines of 51 points with the other predictors drawn once per line,
# and the largest fall of any draw along any line is kept as a share of the
# truth's sd (zero for a monotone fit; the free arm's is its own).
#
# Rule: mBART's values become the tree prior under every monotone fit if
#   (a) at 200 trees, in every design, noise and n cell, joint.mbart minus
#       joint is at most 2 paired se in RMSE and at least -2 paired se in
#       score, the se over replicates; and
#   (b) joint.mbart's mean coverage is at least 0.90 in every cell at 200 and
#       50 trees.
# Otherwise cgm()'s defaults stay. report computes the verdict and names the
# cells failing each clause.
#
# Usage:
#   run out=<dir> [quick] [cores=4] [reps=8]
#     one result file per job in <dir>; a job whose file exists is skipped, so
#     a stopped run resumes by running it again; <dir>/DONE is written last.
#     quick: 2 replicates, 50 + 50 sweeps.
#   report out=<dir> [csv=<file>]
#     summary tables and the verdict on stdout; csv writes one row per fit

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

arms <- c("joint", "joint.mbart", "free")
treeCounts <- c(200L, 50L)
noises <- c(low = 1 / 3, high = 1)
sizes <- c(2000L, 200L)
nReps <- as.integer(option("reps", if (quick) "2" else "8"))
nBurn <- if (quick) 50L else 500L
nKept <- if (quick) 50L else 500L
nHeldOut <- 2000L
nLines <- 3L
nLinePoints <- 51L

designs <- list(
  friedman = list(
    p = 10L,
    f = function(x) {
      10 *
        sin(pi * x[, 1L] * x[, 2L] / 2) +
        20 * (x[, 3L] - 0.5)^2 +
        10 * x[, 4L] +
        5 * x[, 5L]
    },
    constrained = c(1L, 2L, 4L, 5L),
    free = 3L
  ),
  additive = list(
    p = 10L,
    f = function(x) {
      2 *
        x[, 1L] +
        (x[, 2L] > 0.5) +
        exp(2 * x[, 3L]) / exp(2) -
        sin(2 * pi * x[, 4L])
    },
    constrained = c(1L, 2L, 3L),
    free = 4L
  ),
  interaction = list(
    p = 5L,
    f = function(x) {
      x[, 1L] * (1 + 4 * x[, 2L]) + x[, 1L] * x[, 3L] + sin(2 * pi * x[, 2L])
    },
    constrained = 1L,
    free = c(2L, 3L)
  )
)
for (name in names(designs)) {
  d <- designs[[name]]
  designs[[name]]$noise <- setdiff(seq_len(d$p), c(d$constrained, d$free))
}

drawX <- function(design, n) {
  x <- matrix(stats::runif(n * design$p), n)
  colnames(x) <- paste0("x", seq_len(design$p))
  x
}

# the truth's sd under the design's law, by Monte Carlo on its own seed
truthSd <- function(design) {
  set.seed(1L)
  stats::sd(design$f(drawX(design, 1000000L)))
}

# nondecreasing in each constrained predictor: f at a point against f with
# that predictor moved up, on 10^5 random pairs
checkTruths <- function() {
  set.seed(2L)
  for (name in names(designs)) {
    design <- designs[[name]]
    x <- drawX(design, 100000L)
    for (j in design$constrained) {
      up <- x
      up[, j] <- x[, j] + stats::runif(nrow(x)) * (1 - x[, j])
      if (any(design$f(up) < design$f(x) - 1e-12)) {
        stop("design ", name, " falls in x", j)
      }
    }
  }
}

# every monotone arm names its prior; monotone() and cgm() resolve by bare
# name inside bart()
fitArm <- function(arm, design, train, test, n.trees, seed) {
  directions <- stats::setNames(
    rep("increasing", length(design$constrained)),
    paste0("x", design$constrained)
  )
  common <- list(
    formula = y ~ .,
    data = train,
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
    joint = bquote(bart(
      monotone = monotone(.(directions), prior = "joint")
    )),
    joint.mbart = bquote(bart(
      monotone = monotone(.(directions), prior = "joint"),
      tree.prior = cgm(power = 0.8, base = 0.25)
    )),
    free = quote(bart())
  )
  call <- as.call(c(as.list(call), common))
  eval(call)
}

makeData <- function(design, noiseShare, n, seed) {
  fSd <- truthSd(design)
  set.seed(seed)
  sd <- noiseShare * fSd
  x <- drawX(design, n)
  train <- data.frame(x, y = design$f(x) + stats::rnorm(n, 0, sd))
  heldX <- drawX(design, nHeldOut)
  heldF <- design$f(heldX)
  heldY <- heldF + stats::rnorm(nHeldOut, 0, sd)
  # lines along each constrained predictor, the rest drawn once per line
  lines <- do.call(
    rbind,
    lapply(design$constrained, function(j) {
      base <- drawX(design, nLines)
      do.call(
        rbind,
        lapply(seq_len(nLines), function(l) {
          line <- base[rep(l, nLinePoints), , drop = FALSE]
          line[, j] <- seq(0, 1, length.out = nLinePoints)
          line
        })
      )
    })
  )
  list(
    train = train,
    test = as.data.frame(rbind(heldX, lines)),
    heldF = heldF,
    heldY = heldY,
    fSd = fSd,
    noiseSd = sd
  )
}

fitMeasures <- function(fit, design, data, n.trees) {
  draws <- fit$yhat.test
  held <- draws[, seq_len(nHeldOut), drop = FALSE]
  postMean <- colMeans(held)
  lower <- apply(held, 2L, stats::quantile, 0.025)
  upper <- apply(held, 2L, stats::quantile, 0.975)
  truth <- data$heldF
  # log predictive score: per point, log of the mean over draws of the
  # normal density at the fresh response; sigma recycles down the rows
  logDensity <- stats::dnorm(
    matrix(data$heldY, nrow(held), nHeldOut, byrow = TRUE),
    held,
    fit$sigma,
    log = TRUE
  )
  peak <- apply(logDensity, 2L, max)
  score <- peak + log(colMeans(exp(logDensity - rep(peak, each = nrow(held)))))
  # largest fall along the lines, each line's points consecutive
  lineDraws <- draws[, -seq_len(nHeldOut), drop = FALSE]
  falls <- vapply(
    seq_len(ncol(lineDraws) / nLinePoints),
    function(l) {
      d <- lineDraws[, (l - 1L) * nLinePoints + seq_len(nLinePoints)]
      max(0, -apply(d, 1L, diff))
    },
    numeric(1L)
  )
  varcount <- fit$varcount
  shares <- colSums(varcount) / sum(varcount)
  perPredictor <- stats::setNames(
    rep(NA_real_, 10L),
    paste0("share.x", seq_len(10L))
  )
  perPredictor[seq_len(design$p)] <- shares
  c(
    rmse = sqrt(mean((postMean - truth)^2)),
    coverage = mean(truth >= lower & truth <= upper),
    width = mean(upper - lower),
    log.score = mean(score),
    max.fall = max(falls) / data$fSd,
    share.constrained = sum(shares[design$constrained]),
    share.free = sum(shares[design$free]),
    share.noise = sum(shares[design$noise]),
    splits.per.tree = mean(rowSums(varcount)) / n.trees,
    perPredictor
  )
}

runFitJob <- function(job) {
  design <- designs[[job$design]]
  data <- makeData(design, noises[[job$noise]], job$n, job$dataSeed)
  order <- arms[(seq_along(arms) + job$rep - 2L) %% length(arms) + 1L]
  rows <- lapply(order, function(arm) {
    start <- proc.time()
    fit <- fitArm(arm, design, data$train, data$test, job$trees, job$fitSeed)
    used <- proc.time() - start
    if (
      arm != "free" && !identical(fit$monotone.prior, sub("[.].*", "", arm))
    ) {
      stop("arm ", arm, " fit under prior ", format(fit$monotone.prior))
    }
    sweeps <- nBurn + nKept
    c(
      list(arm = arm, position = match(arm, order)),
      as.list(fitMeasures(fit, design, data, job$trees)),
      list(
        truth.sd = data$fSd,
        noise.sd = data$noiseSd,
        wall.ms.per.sweep = 1000 * used[["elapsed"]] / sweeps,
        cpu.ms.per.sweep = 1000 *
          (used[["user.self"]] + used[["sys.self"]]) /
          sweeps,
        cpu.s = used[["user.self"]] + used[["sys.self"]]
      )
    )
  })
  do.call(rbind, lapply(rows, as.data.frame))
}

# ---- jobs ------------------------------------------------------------------

jobList <- function() {
  jobs <- list()
  # the heaviest first, for balance across workers
  for (n in sizes) {
    for (trees in treeCounts) {
      for (design in names(designs)) {
        for (noise in names(noises)) {
          for (rep in seq_len(nReps)) {
            dataSeed <- 100000L *
              match(design, names(designs)) +
              10000L * match(noise, names(noises)) +
              100L * n %/% 100L +
              rep
            jobs[[length(jobs) + 1L]] <- list(
              id = sprintf(
                "fit-%s-%s-n%04d-t%03d-r%d",
                design,
                noise,
                n,
                trees,
                rep
              ),
              design = design,
              noise = noise,
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
  result <- tryCatch(runFitJob(job), error = function(e) e)
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
  checkTruths()
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
  if (!length(files)) {
    stop("no results in ", outDir)
  }
  fits <- do.call(
    rbind,
    lapply(files, function(file) {
      r <- readRDS(file)
      data.frame(
        design = r$job$design,
        noise = r$job$noise,
        n = r$job$n,
        trees = r$job$trees,
        rep = r$job$rep,
        r$result
      )
    })
  )
  fits$arm <- factor(fits$arm, arms)
  cells <- c("design", "noise", "n", "trees")
  cellOrder <- function(table) {
    order(
      match(table$design, names(designs)),
      match(table$noise, names(noises)),
      -table$n,
      -table$trees
    )
  }
  options(width = 200L, digits = 3L)
  measures <- c(
    "rmse",
    "coverage",
    "width",
    "log.score",
    "max.fall",
    "share.constrained",
    "share.free",
    "share.noise",
    "splits.per.tree",
    "cpu.ms.per.sweep",
    "wall.ms.per.sweep"
  )
  for (measure in measures) {
    cat("\n", measure, ": mean over replicates\n", sep = "")
    table <- stats::aggregate(fits[[measure]], fits[c("arm", cells)], mean)
    wide <- stats::reshape(
      table,
      idvar = cells,
      timevar = "arm",
      direction = "wide"
    )
    names(wide) <- sub("^x[.]", "", names(wide))
    print(wide[cellOrder(wide), ], row.names = FALSE)
  }

  # paired over replicates: joint.mbart minus joint on the same data set
  job <- interaction(fits[c(cells, "rep")], drop = TRUE)
  paired <- function(measure) {
    joint <- fits[[measure]][fits$arm == "joint"]
    mbart <- fits[[measure]][fits$arm == "joint.mbart"]
    mbart -
      joint[match(job[fits$arm == "joint.mbart"], job[fits$arm == "joint"])]
  }
  mbartRows <- fits[fits$arm == "joint.mbart", c(cells, "rep", "coverage")]
  mbartRows$rmse.diff <- paired("rmse")
  mbartRows$score.diff <- paired("log.score")
  se <- function(d) stats::sd(d) / sqrt(length(d))
  verdict <- do.call(
    rbind,
    lapply(
      split(mbartRows, mbartRows[cells], drop = TRUE),
      function(cell) {
        data.frame(
          cell[1L, cells],
          reps = nrow(cell),
          rmse.diff = mean(cell$rmse.diff),
          rmse.se = se(cell$rmse.diff),
          score.diff = mean(cell$score.diff),
          score.se = se(cell$score.diff),
          mbart.coverage = mean(cell$coverage)
        )
      }
    )
  )
  verdict <- verdict[cellOrder(verdict), ]
  verdict$fails.a.rmse <- verdict$trees == 200L &
    !(verdict$rmse.diff <= 2 * verdict$rmse.se)
  verdict$fails.a.score <- verdict$trees == 200L &
    !(verdict$score.diff >= -2 * verdict$score.se)
  verdict$fails.b <- !(verdict$mbart.coverage >= 0.90)
  cat("\njoint.mbart minus joint, paired over replicates, and mBART coverage\n")
  print(verdict, row.names = FALSE)

  cellName <- function(rows) {
    if (!nrow(rows)) {
      return("none")
    }
    paste(
      sprintf("%s/%s/n%d/t%d", rows$design, rows$noise, rows$n, rows$trees),
      collapse = ", "
    )
  }
  failsA <- verdict$fails.a.rmse | verdict$fails.a.score
  failsB <- verdict$fails.b
  expected <- length(designs) *
    length(noises) *
    length(sizes) *
    length(treeCounts)
  cat("\nRule\n")
  cat(
    "(a) at 200 trees, RMSE diff <= 2 se and score diff >= -2 se, every cell\n"
  )
  cat("    failing RMSE:", cellName(verdict[verdict$fails.a.rmse, ]), "\n")
  cat("    failing score:", cellName(verdict[verdict$fails.a.score, ]), "\n")
  cat("(b) mBART coverage >= 0.90, every cell at 200 and 50 trees\n")
  cat("    failing:", cellName(verdict[failsB, ]), "\n")
  complete <- nrow(verdict) == expected && all(verdict$reps >= 2L)
  cat(
    "verdict:",
    if (!complete) {
      sprintf(
        "incomplete (%d of %d cells, or fewer than 2 reps)",
        nrow(verdict),
        expected
      )
    } else if (any(failsA) || any(failsB)) {
      "cgm() defaults stay"
    } else {
      "mBART's cgm(power = 0.8, base = 0.25) becomes the monotone tree prior"
    },
    "\n"
  )
  cat(
    sprintf(
      "\n%d fits, %.2f CPU-hours in fits\n",
      nrow(fits),
      sum(fits$cpu.s) / 3600
    )
  )

  csv <- option("csv", NA_character_)
  if (!is.na(csv)) {
    numeric <- vapply(fits, is.double, TRUE)
    fits[numeric] <- lapply(fits[numeric], signif, 4L)
    utils::write.csv(fits, csv, row.names = FALSE)
  }
}
