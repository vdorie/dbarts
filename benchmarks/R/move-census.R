#!/usr/bin/env Rscript

# Stage 0 of the tree-mixing falsifier: the move census
# (docs/design/tree-mixing-proposals.md section 6.1). Pure measurement - no
# kill criterion, no engine behaviour change - answering the three questions
# the veto-scaffold counts left open: the distribution of the log-likelihood
# difference among REJECTED structural proposals, change-move proposals and
# acceptances by target node depth, and the cut displacement that puts a
# same-variable cut move at a workable acceptance rate.
#
# The numbers come from the engine's own acceptance terms, appended one line
# per structural proposal by scaffolding in the move kernels that is compiled
# only under -DBARTCORE_MOVE_CENSUS. THIS SCRIPT DOES NOTHING against an
# ordinary build: the library must be built with the flag.
#
# Building the instrumented library. src/Makevars sets PKG_CPPFLAGS with '=',
# so the flag rides CPPFLAGS instead, which R composes into ALL_CPPFLAGS and
# no package makefile touches. R CMD INSTALL reads the user Makevars LAST, so
# an append there survives:
#
#   printf 'CPPFLAGS += -DBARTCORE_MOVE_CENSUS\n' > /tmp/census-makevars
#   mkdir -p /tmp/census-lib
#   R_MAKEVARS_USER=/tmp/census-makevars \
#     R CMD INSTALL --preclean --library=/tmp/census-lib .
#   R_LIBS=/tmp/census-lib Rscript benchmarks/R/move-census.R
#
# A private --library keeps the instrumented build off the ordinary one; the
# flag also reaches misc/external/rc, which ignore it. The run mode stops with
# a diagnostic if the first cell writes no census file, which is what an
# uninstrumented library looks like from R.
#
# Five cells, 200 burn-in plus 500 sampled sweeps, one chain, one thread,
# fixed seeds:
#
#   default   dbarts defaults: Friedman, n = 5000, p = 10, m = 75, sigma = 1
#   lownoise  the same at sigma^2 = 0.1, where structure freezing is
#             established (section 3.2)
#   wide      Friedman p = 50 with 45 noise columns, the regime that killed
#             the warm-start default
#   bcf       the causal-forest strong-scale cell: two forests, prognostic
#             amplitude 8 sigma, which is docs/design/bcf.md's strong-|a|
#             regime (|a|/sigma large, where the mu forest's structure mixes
#             slowly at high SNR)
#   c1        the He-Hahn independent design, n = 10000, p = 30, Trig+poly at
#             kappa = 1, the primary benefit cell of the surface battery
#
# Beyond the per-proposal records the census has always taken, the build logs
# five generator-only probes: the closed rule neighbourhood at a nog node
# (weight entropy, the incumbent's share and rank, both jointly over the
# available variables and restricted to the incumbent's), the informed-death
# weights over the nog nodes against the realized uniform pick, the perturb
# proposal's signed displacement, the rule_gibbs kernel's own cut-scan cost,
# and every tree's settled leaf count. Nothing
# there draws or changes a draw; the record format lives in moves.hpp.
#
# Usage:
#   Rscript move-census.R                          run every cell, summarize
#   Rscript move-census.R run [dir] [cell ...]     run cells, write census
#   Rscript move-census.R summarize [dir]          summarize existing files
# Append 'quick' for a smoke test (fewer sweeps, smaller n; not comparable).
# Append 'perturb' to run the perturb-carrying mixture instead of the shipped
# one, which is the only way the signed-displacement probe records anything, or
# 'gibbs' for the rule_gibbs-carrying one, which is the only way the cut-scan
# cost records anything; both skip the bcf cell, whose treatment forest refuses
# the argument.
#
# One cell runs per R process, spawned by the run mode: the engine opens the
# census file once, on the first record, so a second cell in the same process
# would append to the first cell's file.
#
# Summaries are taken over the SAMPLED sweeps only; burn-in acceptance is
# reported beside them for context.

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args
args <- setdiff(args, "quick")
perturbing <- "perturb" %in% args
args <- setdiff(args, "perturb")
drawingRules <- "gibbs" %in% args
args <- setdiff(args, "gibbs")
modes <- c("run", "summarize", "runcell")
mode <- if (length(args) >= 1L && args[[1L]] %in% modes) args[[1L]] else "both"
args <- setdiff(args, modes)
outputDir <- if (length(args) >= 1L) args[[1L]] else "benchmarks/census"
cellArgs <- if (length(args) >= 2L) args[-1L] else character()

nBurn <- if (quick) 20L else 200L
nSamples <- if (quick) 50L else 500L
nObservations <- if (quick) 500L else 5000L
nObservationsC1 <- if (quick) 500L else 10000L
nTrees <- 75L
dataSeed <- 20260907L
samplerSeed <- 7L

scriptDirectory <- function() {
  file <- sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))
  if (length(file) != 1L) {
    stop("this script needs its own path; invoke it with Rscript")
  }
  dirname(file)
}

# ---------------------------------------------------------------- generators

genFriedman <- function(n, p, sigma) {
  x <- matrix(runif(n * p), n, p)
  f <- 10 *
    sin(pi * x[, 1L] * x[, 2L]) +
    20 * (x[, 3L] - 0.5)^2 +
    10 * x[, 4L] +
    5 * x[, 5L]
  list(x = x, y = f + rnorm(n, sd = sigma))
}

# Hahn, Murray and Carvalho's shape with the prognostic surface standardized
# and rescaled to strength * sigma, so the cell's |a|/sigma is set explicitly;
# the propensity is a function of mu, so the confounding rides the amplitude.
genCausal <- function(n, p, sigma, strength) {
  x <- matrix(runif(n * p), n, p)
  mu <- 2 * sin(pi * x[, 1L] * x[, 2L]) + 2 * (x[, 3L] - 0.5)^2 + x[, 4L]
  mu <- strength * sigma * (mu - mean(mu)) / sd(mu)
  tau <- 1 + 2 * x[, 3L]
  propensity <- 0.05 + 0.8 * pnorm(mu / (strength * sigma) - 0.5 * x[, 1L])
  z <- rbinom(n, 1L, propensity)
  list(x = x, z = z, y = mu + z * tau + rnorm(n, sd = sigma))
}

# The shipped mixture, or one of the two carrying a move that ships at zero:
# neither perturb nor rule_gibbs ever fires at its shipped zero, so neither the
# run-length question nor the cut-scan cost can be asked of the default kernel
# at all.
censusProposalProbs <- function() {
  if (perturbing) {
    c(birth_death = 0.5, swap = 0, change = 0.34, perturb = 0.16, birth = 0.5)
  } else if (drawingRules) {
    c(
      birth_death = 0.5,
      swap = 0,
      change = 0.34,
      rule_gibbs = 0.16,
      birth = 0.5
    )
  } else {
    c(birth_death = 0.6, swap = 0, change = 0.4, perturb = 0, birth = 0.5)
  }
}

censusControl <- function() {
  dbartsControl(
    verbose = FALSE,
    n.trees = nTrees,
    n.chains = 1L,
    n.threads = 1L,
    n.samples = nSamples,
    updateState = FALSE,
    seed = samplerSeed
  )
}

# ---------------------------------------------------------------- the cells

cells <- list(
  default = function() {
    data <- genFriedman(nObservations, 10L, 1)
    dbarts(
      data$x,
      data$y,
      proposal.probs = censusProposalProbs(),
      control = censusControl()
    )
  },
  lownoise = function() {
    data <- genFriedman(nObservations, 10L, sqrt(0.1))
    dbarts(
      data$x,
      data$y,
      proposal.probs = censusProposalProbs(),
      control = censusControl()
    )
  },
  wide = function() {
    data <- genFriedman(nObservations, 50L, 1)
    dbarts(
      data$x,
      data$y,
      proposal.probs = censusProposalProbs(),
      control = censusControl()
    )
  },
  bcf = function() {
    data <- genCausal(nObservations, 10L, 1, strength = 8)
    z <- data$z
    dbarts(
      data$x,
      data$y,
      forests = list(forest(), forest(basis = ~ factor(z))),
      proposal.probs = censusProposalProbs(),
      control = censusControl()
    )
  },
  c1 = function() {
    source(
      file.path(scriptDirectory(), "surfaces", "surfaces-common.R"),
      chdir = FALSE
    )
    data <- surfacesHeHahn(
      nObservationsC1,
      0L,
      30L,
      "trigpoly",
      1,
      design = "independent"
    )
    dbarts(
      data$x,
      data$y,
      proposal.probs = censusProposalProbs(),
      control = censusControl()
    )
  }
)

censusFile <- function(dir, cell) {
  file.path(dir, sprintf("move-census-%s.csv", cell))
}

runCell <- function(dir, cell) {
  file <- censusFile(dir, cell)
  unlink(file)
  Sys.setenv(BARTCORE_MOVE_CENSUS_FILE = file)
  on.exit(Sys.unsetenv("BARTCORE_MOVE_CENSUS_FILE"), add = TRUE)
  set.seed(dataSeed)
  sampler <- cells[[cell]]()
  elapsed <- system.time(invisible(sampler$run(nBurn, nSamples)))[["elapsed"]]
  rm(sampler)
  invisible(gc(FALSE))
  # the engine's stream is still buffered here, so existence is the test: an
  # uninstrumented library never opens the file at all
  if (!file.exists(file)) {
    stop(
      "no census records for cell '",
      cell,
      "': the loaded dbarts was not built with -DBARTCORE_MOVE_CENSUS ",
      "(see this file's header)"
    )
  }
  cat(sprintf(
    "%-9s %6.1f s  %s (%.1f MB)\n",
    cell,
    elapsed,
    file,
    file.size(file) / 1024^2
  ))
}

# ------------------------------------------------------------- the summaries

proposalNames <- c(
  "kind",
  "sweep",
  "forest",
  "tree",
  "move",
  "noop",
  "accepted",
  "nodeDepth",
  "treeDepth",
  "interior",
  "logLik",
  "logPrior",
  "logCorr"
)
cutNames <- c(
  "kind",
  "sweep",
  "forest",
  "tree",
  "nodeDepth",
  "displacement",
  "logRatio"
)
nogNames <- c(
  "kind",
  "sweep",
  "forest",
  "tree",
  "node",
  "nodeDepth",
  "isNog",
  "interior",
  "nog",
  "scanned",
  "jointCandidates",
  "jointEntropy",
  "jointIncumbent",
  "jointMaximum",
  "jointRank",
  "cutCandidates",
  "cutEntropy",
  "cutIncumbent",
  "cutMaximum",
  "cutRank"
)
deathNames <- c(
  "kind",
  "sweep",
  "forest",
  "tree",
  "candidates",
  "entropy",
  "pickWeight",
  "pickRank",
  "maxWeight"
)
runNames <- c(
  "kind",
  "sweep",
  "forest",
  "tree",
  "node",
  "current",
  "target",
  "accepted"
)
gibbsNames <- c(
  "kind",
  "sweep",
  "forest",
  "tree",
  "eligible",
  "scanned",
  "candidates",
  "stratum"
)
shapeNames <- c("kind", "sweep", "forest", "tree", "leaves", "interior", "nog")

# an absent record kind is an empty frame with the right columns, not an
# error: perturb records nothing at the shipped mixture, and the neighbourhood
# probe writes nothing for a leaf model the cut scan cannot score
readRecords <- function(lines, names) {
  if (length(lines) == 0L) {
    frame <- as.data.frame(
      matrix(numeric(0L), 0L, length(names) - 1L),
      stringsAsFactors = FALSE
    )
    names(frame) <- names[-1L]
    return(frame)
  }
  read.csv(text = lines, header = FALSE, col.names = names)[, -1L]
}

readCensus <- function(file) {
  lines <- readLines(file)
  kind <- substr(lines, 1L, 1L)
  list(
    proposals = read.csv(
      text = lines[kind == "p"],
      header = FALSE,
      col.names = proposalNames
    ),
    cuts = read.csv(
      text = lines[kind == "d"],
      header = FALSE,
      col.names = cutNames
    ),
    nog = readRecords(lines[kind == "g"], nogNames),
    deaths = readRecords(lines[kind == "x"], deathNames),
    runs = readRecords(lines[kind == "r"], runNames),
    gibbs = readRecords(lines[kind == "n"], gibbsNames),
    shapes = readRecords(lines[kind == "t"], shapeNames)
  )
}

# proposals, no-op rate, and the two acceptance denominators the veto scaffold
# separated: per proposal made and per proposal that reached a score
moveTable <- function(p) {
  by <- split(p, p$move)
  rows <- lapply(names(by), function(move) {
    m <- by[[move]]
    scored <- sum(m$noop == 0L)
    data.frame(
      move = move,
      proposals = nrow(m),
      noop.pct = 100 * mean(m$noop == 1L),
      accept.pct = 100 * mean(m$accepted == 1L),
      accept.scored.pct = if (scored > 0L) {
        100 * sum(m$accepted) / scored
      } else {
        NA
      },
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  scored <- sum(p$noop == 0L)
  rbind(
    out,
    data.frame(
      move = "all",
      proposals = nrow(p),
      noop.pct = 100 * mean(p$noop == 1L),
      accept.pct = 100 * mean(p$accepted == 1L),
      accept.scored.pct = if (scored > 0L) {
        100 * sum(p$accepted) / scored
      } else {
        NA
      },
      stringsAsFactors = FALSE
    )
  )
}

# The fork the whole program turns on: are rejections close calls (scale, so
# the temperature family is live) or wrong proposals (only better-aimed
# proposals help)?
rejectionTable <- function(p) {
  r <- p[p$noop == 0L & p$accepted == 0L, ]
  rows <- lapply(c(sort(unique(r$move)), "all"), function(move) {
    d <- if (move == "all") r$logLik else r$logLik[r$move == move]
    finite <- d[is.finite(d)]
    q <- quantile(finite, c(0.01, 0.05, 0.25, 0.5, 0.75), names = FALSE)
    data.frame(
      move = move,
      rejected = length(d),
      vetoed.pct = 100 * mean(!is.finite(d)),
      q01 = q[1L],
      q05 = q[2L],
      q25 = q[3L],
      q50 = q[4L],
      q75 = q[5L],
      within1.pct = 100 * mean(abs(finite) <= 1),
      within2.pct = 100 * mean(abs(finite) <= 2),
      within5.pct = 100 * mean(abs(finite) <= 5),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# Section 3.3's never-measured claim: does change acceptance collapse with the
# depth of the node whose rule is redrawn?
changeDepthTable <- function(p) {
  c0 <- p[p$move == "change" & p$noop == 0L, ]
  if (nrow(c0) == 0L) {
    return(NULL)
  }
  by <- split(c0, c0$nodeDepth)
  do.call(
    rbind,
    lapply(names(by), function(depth) {
      d <- by[[depth]]
      data.frame(
        depth = as.integer(depth),
        proposals = nrow(d),
        accepted = sum(d$accepted),
        accept.pct = 100 * mean(d$accepted == 1L),
        stringsAsFactors = FALSE
      )
    })
  )
}

# The window-width input: the acceptance a same-variable cut move would have
# had at each displacement, min(1, exp(logRatio)) averaged over the interior
# nodes a change proposal visited. Signed displacements are pooled by
# magnitude after reporting both signs' counts.
cutTable <- function(d) {
  if (nrow(d) == 0L) {
    return(NULL)
  }
  d$magnitude <- abs(d$displacement)
  by <- split(d, d$magnitude)
  do.call(
    rbind,
    lapply(names(by), function(magnitude) {
      m <- by[[magnitude]]
      alpha <- pmin(1, exp(m$logRatio))
      data.frame(
        displacement = as.integer(magnitude),
        probes = nrow(m),
        left = sum(m$displacement < 0L),
        accept.pct = 100 * mean(alpha),
        vetoed.pct = 100 * mean(!is.finite(m$logRatio)),
        median.logratio = median(m$logRatio[is.finite(m$logRatio)]),
        stringsAsFactors = FALSE
      )
    })
  )
}

# Section 15.3 rows 1 to 3: the closed rule neighbourhood at a nog node, priced at
# every change proposal. A collapsed Gibbs draw there is worth making only if
# the conditional is not already sitting on the incumbent, so the readout is
# the entropy in nats, P(incumbent) and how often the incumbent is the mode.
# "joint" ranges over the available variables and their cuts, "cut" over the
# incumbent variable's cuts alone.
nogTable <- function(g) {
  if (nrow(g) == 0L) {
    return(NULL)
  }
  scanned <- g[g$scanned == 1L, ]
  if (nrow(scanned) == 0L) {
    return(NULL)
  }
  summary <- function(forest, label, f, prefix) {
    data.frame(
      forest = forest,
      neighbourhood = label,
      proposals = nrow(f),
      candidates = median(f[[paste0(prefix, "Candidates")]]),
      entropy = median(f[[paste0(prefix, "Entropy")]]),
      p.incumbent = median(f[[paste0(prefix, "Incumbent")]]),
      incumbent.top.pct = 100 * mean(f[[paste0(prefix, "Rank")]] == 1),
      rank = median(f[[paste0(prefix, "Rank")]]),
      max.weight = median(f[[paste0(prefix, "Maximum")]]),
      stringsAsFactors = FALSE
    )
  }
  by <- split(scanned, scanned$forest)
  do.call(
    rbind,
    lapply(names(by), function(forest) {
      f <- by[[forest]]
      rbind(
        summary(as.integer(forest), "joint", f, "joint"),
        summary(as.integer(forest), "cut-only", f, "cut")
      )
    })
  )
}

# The two shares the neighbourhood table is conditional on: how often a change
# proposal lands on a nog node at all, and how much of the tree is nog.
nogShareTable <- function(g) {
  if (nrow(g) == 0L) {
    return(NULL)
  }
  by <- split(g, g$forest)
  do.call(
    rbind,
    lapply(names(by), function(forest) {
      f <- by[[forest]]
      data.frame(
        forest = as.integer(forest),
        proposals = nrow(f),
        target.nog.pct = 100 * mean(f$isNog == 1L),
        nog.share.pct = 100 * sum(f$nog) / sum(f$interior),
        scanned.pct = 100 * mean(f$scanned == 1L),
        stringsAsFactors = FALSE
      )
    })
  )
}

# Section 15.3 row 3: informed death. The nog nodes weighted by exp of the merged-leaf
# marginal ratio, against the uniform pick the kernel actually made. A uniform
# pick already sitting at the weighted mode leaves the weighting nothing to buy.
deathTable <- function(x) {
  if (nrow(x) == 0L) {
    return(NULL)
  }
  scored <- x[!is.na(x$entropy), ]
  if (nrow(scored) == 0L) {
    return(NULL)
  }
  row <- function(forest, label, f) {
    data.frame(
      forest = forest,
      nog = label,
      proposals = nrow(f),
      candidates = median(f$candidates),
      entropy = median(f$entropy),
      pick.weight = median(f$pickWeight),
      pick.rank = median(f$pickRank),
      pick.top.pct = 100 * mean(f$pickRank == 1),
      max.weight = median(f$maxWeight),
      stringsAsFactors = FALSE
    )
  }
  by <- split(scored, scored$forest)
  do.call(
    rbind,
    lapply(names(by), function(forest) {
      f <- by[[forest]]
      # a tree with one nog node leaves the uniform pick nothing to be wrong
      # about, so the weighting can only pay where there are at least two
      several <- f[f$candidates >= 2, ]
      rows <- row(as.integer(forest), "all", f)
      if (nrow(several) > 0L) {
        rows <- rbind(rows, row(as.integer(forest), ">= 2", several))
      }
      rows
    })
  )
}

# Section 16.3 row 2: whether accepted cut displacements at one node run in the same
# direction. A run continues only while the next accepted displacement starts
# where the last one landed, which is what keeps a reused arena slot or an
# intervening change move from splicing two nodes' histories together. Under a
# reversible walk the sign is a fair coin, so run lengths are geometric with
# P(k) = 2^-k and mean 2.
runPairs <- function(r) {
  accepted <- r[r$accepted == 1L, ]
  if (nrow(accepted) < 2L) {
    return(NULL)
  }
  accepted <- accepted[
    order(accepted$forest, accepted$tree, accepted$node, accepted$sweep),
  ]
  n <- nrow(accepted)
  direction <- sign(accepted$target - accepted$current)
  list(
    # the chain continued: same node, and the next displacement starts where
    # the last one landed
    continues = accepted$forest[-1L] == accepted$forest[-n] &
      accepted$tree[-1L] == accepted$tree[-n] &
      accepted$node[-1L] == accepted$node[-n] &
      accepted$current[-1L] == accepted$target[-n],
    same = direction[-1L] == direction[-n]
  )
}

# The pairwise form of the same question, immune to the censoring below: over
# consecutive accepted displacements that DID continue the chain, a reversible
# walk makes the direction a fair coin.
continuationTable <- function(r) {
  pairs <- runPairs(r)
  if (is.null(pairs) || sum(pairs$continues) == 0L) {
    return(NULL)
  }
  data.frame(
    pairs = sum(pairs$continues),
    same.direction.pct = 100 * mean(pairs$same[pairs$continues]),
    reversible.pct = 50,
    stringsAsFactors = FALSE
  )
}

# Run length as a HAZARD rather than a distribution: a streak the chain broke
# under - an intervening change move, a reused arena slot, the end of the
# record - is censored, and censoring falls on long streaks, so the run-length
# histogram itself is not readable against any null. What is readable is the
# chance a streak already k displacements long extends by one more, taken over
# the pairs where the chain did continue. The geometric null a reversible walk
# implies is 50 percent at every k, and the run-length distribution it carries
# is P(k) = 2^-k.
runTable <- function(r) {
  pairs <- runPairs(r)
  if (is.null(pairs) || sum(pairs$continues) == 0L) {
    return(NULL)
  }
  extends <- pairs$continues & pairs$same
  starts <- which(!c(FALSE, extends))
  position <- sequence(diff(c(starts, length(extends) + 2L)))
  at <- pmin(position[-length(position)], 5L)[pairs$continues]
  extended <- pairs$same[pairs$continues]
  do.call(
    rbind,
    lapply(sort(unique(at)), function(k) {
      data.frame(
        length = if (k < 5L) as.character(k) else "5+",
        streaks = sum(at == k),
        extends.pct = 100 * mean(extended[at == k]),
        reversible.pct = 50,
        stringsAsFactors = FALSE
      )
    })
  )
}

# Section 16.4: the leaf count itself, which section 16.2 could only bound by Jensen.
# One record per tree per sweep, so the quantiles are over trees and sweeps
# together and the mean column is the mean of the per-sweep means.
leafTable <- function(t) {
  if (nrow(t) == 0L) {
    return(NULL)
  }
  by <- split(t, t$forest)
  do.call(
    rbind,
    lapply(names(by), function(forest) {
      f <- by[[forest]]
      q <- quantile(f$leaves, c(0.05, 0.5, 0.95), names = FALSE)
      perSweep <- tapply(f$leaves, f$sweep, mean)
      data.frame(
        forest = as.integer(forest),
        trees = length(unique(f$tree)),
        mean.leaves = mean(perSweep),
        sd.sweep.mean = sd(perSweep),
        q05 = q[1L],
        q50 = q[2L],
        q95 = q[3L],
        max = max(f$leaves),
        stump.pct = 100 * mean(f$leaves == 1L),
        nog.share.pct = 100 * sum(f$nog) / sum(f$interior),
        stringsAsFactors = FALSE
      )
    })
  )
}

# The rule_gibbs kernel's cost, which no other record carries: its own probe
# rides changeMove, the branch a nonzero rule_gibbs share does not take. One
# cut-scan unit is one scanOrdinalCuts pass over the chosen node's members, so
# the move's cost per sweep is the eligible-node reach times the variables
# scanned there, summed over the trees a sweep touches.
gibbsCostTable <- function(g) {
  if (nrow(g) == 0L) {
    return(NULL)
  }
  by <- split(g, g$forest)
  do.call(
    rbind,
    lapply(names(by), function(forest) {
      f <- by[[forest]]
      reached <- f[f$eligible > 0L & f$scanned > 0L, ]
      data.frame(
        forest = as.integer(forest),
        proposals = nrow(f),
        reach.pct = 100 * mean(f$eligible > 0L),
        eligible = median(f$eligible),
        scanned = if (nrow(reached) > 0L) median(reached$scanned) else NA,
        candidates = if (nrow(reached) > 0L) {
          median(reached$candidates)
        } else {
          NA
        },
        vetoed.pct = 100 * mean(f$stratum > 0L),
        scans.per.sweep = sum(f$scanned) / length(unique(f$sweep)),
        stringsAsFactors = FALSE
      )
    })
  )
}

roundFrame <- function(x, digits = 3L) {
  numeric <- vapply(x, is.numeric, logical(1L))
  x[numeric] <- lapply(x[numeric], round, digits)
  x
}

summarizeCell <- function(dir, cell) {
  file <- censusFile(dir, cell)
  if (!file.exists(file)) {
    cat("\n== ", cell, ": no census file at ", file, "\n", sep = "")
    return(invisible(NULL))
  }
  census <- readCensus(file)
  p <- census$proposals
  d <- census$cuts
  g <- census$nog
  x <- census$deaths
  r <- census$runs
  gibbs <- census$gibbs
  shapes <- census$shapes
  sampled <- p$sweep >= nBurn
  cat(
    "\n== ",
    cell,
    " (",
    format(nrow(p), big.mark = ","),
    " proposals, ",
    length(unique(p$forest)),
    " forest(s), ",
    sum(!sampled),
    " burn-in and ",
    sum(sampled),
    " sampled)\n",
    sep = ""
  )
  cat("\nper move, burn-in sweeps:\n")
  print(roundFrame(moveTable(p[!sampled, ]), 2L), row.names = FALSE)
  p <- p[sampled, ]
  d <- d[d$sweep >= nBurn, ]
  g <- g[g$sweep >= nBurn, ]
  x <- x[x$sweep >= nBurn, ]
  r <- r[r$sweep >= nBurn, ]
  gibbs <- gibbs[gibbs$sweep >= nBurn, ]
  shapes <- shapes[shapes$sweep >= nBurn, ]

  cat("\nper move, sampled sweeps:\n")
  print(roundFrame(moveTable(p), 2L), row.names = FALSE)
  cat("\nlog-likelihood difference among rejected proposals:\n")
  print(roundFrame(rejectionTable(p), 2L), row.names = FALSE)
  cat("\nchange proposals by target node depth:\n")
  print(roundFrame(changeDepthTable(p), 2L), row.names = FALSE)
  cat("\nsame-variable cut move, acceptance against displacement:\n")
  print(roundFrame(cutTable(d), 2L), row.names = FALSE)

  cat("\nchange target and the tree's nog share:\n")
  print(roundFrame(nogShareTable(g), 2L), row.names = FALSE)
  cat("\nclosed rule neighbourhood at a nog node:\n")
  print(roundFrame(nogTable(g), 4L), row.names = FALSE)
  cat("\ninformed-death weights against the uniform pick:\n")
  print(roundFrame(deathTable(x), 4L), row.names = FALSE)
  cat("\nperturb: direction of consecutive accepted displacements:\n")
  print(roundFrame(continuationTable(r), 2L), row.names = FALSE)
  cat("\nperturb: streak extension by streak length so far:\n")
  print(roundFrame(runTable(r), 2L), row.names = FALSE)
  cat("\nrule_gibbs: nodes reached and cut scans taken:\n")
  print(roundFrame(gibbsCostTable(gibbs), 2L), row.names = FALSE)
  cat("\nleaves per tree:\n")
  print(roundFrame(leafTable(shapes), 2L), row.names = FALSE)
  invisible(NULL)
}

# ------------------------------------------------------------------- driver

selected <- if (length(cellArgs) > 0L) {
  intersect(names(cells), cellArgs)
} else if (perturbing || drawingRules) {
  # a treatment forest refuses a non-default 'proposal.probs' outright, so a
  # mixture carrying either zero-default move cannot be put to the
  # causal-forest cell at all
  setdiff(names(cells), "bcf")
} else {
  names(cells)
}
if (length(selected) == 0L) {
  stop(
    "no known cell selected; cells are: ",
    paste(names(cells), collapse = ", ")
  )
}

# one cell per process; see the header
spawnCell <- function(dir, cell) {
  file <- file.path(scriptDirectory(), "move-census.R")
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(
      file,
      "runcell",
      dir,
      cell,
      if (quick) "quick",
      if (perturbing) "perturb",
      if (drawingRules) "gibbs"
    )
  )
  if (status != 0L) {
    stop("cell '", cell, "' failed")
  }
}

if (mode %in% c("run", "both")) {
  dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)
  cat(
    "writing census files to ",
    outputDir,
    " (",
    nBurn,
    " burn-in + ",
    nSamples,
    " sampled sweeps, n = ",
    nObservations,
    ")\n",
    sep = ""
  )
  for (cell in selected) {
    spawnCell(outputDir, cell)
  }
}

if (mode == "runcell") {
  if (length(selected) != 1L) {
    stop("runcell takes exactly one cell")
  }
  runCell(outputDir, selected)
}

if (mode %in% c("summarize", "both")) {
  for (cell in selected) {
    summarizeCell(outputDir, cell)
  }
}
