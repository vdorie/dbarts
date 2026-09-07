#!/usr/bin/env Rscript

# P6, the diagonal shelf with targeted selection: Hahn, Murray and Carvalho
# (2020) Example 1. Two uniform covariates, 250 rows, a homogeneous
# treatment effect of -1, and a prognostic surface with a near-step "shelf"
# along x1 = x2 that the propensity tracks, so a single split on the
# treatment indicator can stand in for the many axis-aligned splits the
# shelf would need. The estimand is the outer one: the average treatment
# effect, not the fit.
#
# Published numbers to reproduce (arXiv 1706.09523v4 Table 1, 200
# replications): standard BART bias 0.27, 95% interval coverage 65%, RMSE
# 0.31; the propensity-augmented BCF prior 0.14, 95% and 0.21.
#
# Primary statistic: ATE bias and 95% interval coverage over replications.
#
# The paper does not print the prognostic function. Three reconstructions
# are run. The primary, `figure`, is calibrated to the paper's own Figure 4,
# which plots mu against the propensity for one realization: it reproduces
# that figure's mu range (-1.9 to 3.0), its propensity range (0.10 to 0.88)
# and the value of mu where the propensity crosses one half. The other two
# are the symmetric near-step family the Figure 3 caption describes ("a
# shelf at the line x1 = x2", mu ranging -3 to 3) at two widths, and are
# there to show how far the answer moves with a choice the paper leaves
# open. No arm of this cell can carry a reproduction verdict on its own.
#
# Estimators:
#   bart   the paper's standard BART arm: one forest over (x1, x2, z), ATE
#          read off the counterfactual contrast
#   bcf    the paper's fix, in dbarts' own shape: an estimated propensity
#          added as a covariate and a second forest whose amplitude the
#          treatment indicator modulates. Shipped mixture only, for the
#          reason the move-set arms below give.
#
# Move-set arms, on the same matched seeds and differing only in
# `proposal.probs`, so that every contrast is paired:
#
#   default     proposal.probs unset, which is the shipped mixture
#               (birth_death 0.6, swap 0, change 0.4, perturb 0)
#   birthdeath  birth_death 1, swap 0, change 0
#   swap        birth_death 0.5, swap 0.1, change 0.4, the former default
#
# The mixture applies to EVERY bart2 call an estimator makes, so an arm is
# one sampler throughout. That confines the move-set contrast to the `bart`
# estimator, which is the one the pre-registered primary and the published
# number attach to: a treatment forest refuses a non-default
# `proposal.probs` outright, so the `bcf` estimator has no move-set arm to
# run and stays on the shipped mixture. `default` is re-run beside the other
# two rather than read off an earlier session, so the contrast is paired
# within one run; on matched seeds it reproduces the single-arm run exactly.
#
# Usage: Rscript P6-diagonal-shelf.R [outputDir] [quick] [moveset ...]

source(
  file.path(
    dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))),
    "surfaces-common.R"
  ),
  chdir = FALSE
)

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args

nReplicates <- if (quick) 5L else 200L
n <- 250L
trueTau <- -1

# NULL leaves proposal.probs unset, which is the shipped mixture.
movesets <- list(
  default = NULL,
  birthdeath = c(
    birth_death = 1,
    swap = 0,
    change = 0,
    perturb = 0,
    birth = 0.5
  ),
  swap = c(
    birth_death = 0.5,
    swap = 0.1,
    change = 0.4,
    perturb = 0,
    birth = 0.5
  )
)
movesetNames <- names(movesets)
selectedMovesets <- intersect(movesetNames, args)
if (length(selectedMovesets) > 0L) {
  movesets <- movesets[selectedMovesets]
}
# The control every other arm is read against.
movesetControl <- "default"

outputDir <- surfacesOutputDir(args, flags = c("quick", movesetNames))

# One bart2 argument list with the arm's proposal mixture attached.
movesetCall <- function(call, probs) {
  if (!is.null(probs)) {
    call$proposal.probs <- probs
  }
  call
}

# The Figure-4-calibrated prognostic surface. A gentle ramp in x1 - x2 plus
# a near-step of width 0.08 centred just inside the positive side of the
# diagonal, which is where the paper's text puts the threshold ("positive
# levels of x1 - x2 being a critical threshold").
figureMu <- function(x1, x2) {
  d <- x1 - x2
  -0.87 + 1.2 * d + 2.75 * pnorm((d - 0.15) / 0.08)
}

designs <- list(
  figure = figureMu,
  shelf015 = function(x1, x2) surfacesShelfMu(x1, x2, 0.15),
  shelf040 = function(x1, x2) surfacesShelfMu(x1, x2, 0.40)
)

drawShelf <- function(muFunction) {
  x1 <- runif(n)
  x2 <- runif(n)
  mu <- muFunction(x1, x2)
  propensity <- surfacesShelfPropensity(mu, x1, x2)
  z <- rbinom(n, 1L, propensity)
  data.frame(
    y = mu + trueTau * z + rnorm(n),
    x1 = x1,
    x2 = x2,
    z = z,
    mu = mu,
    propensity = propensity
  )
}

# Section 6.4's separation condition on its own, for the metrics whose
# margin is a ratio of aggregates rather than a mean difference: does the
# one-sided 95% bound of the paired mean, taken in the worse direction,
# exclude zero? `worse` is the sign of a difference that counts as a
# regression.
armSeparated <- function(x, worse) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) {
    return("-")
  }
  se <- sd(x) / sqrt(length(x))
  if (worse * mean(x) - 1.645 * se > 0) "separated" else "not separated"
}

# ATE draws from a counterfactual contrast: the posterior mean surface
# evaluated at every row twice, once with z = 1 and once with z = 0.
ateFromContrast <- function(fitted1, fitted0) {
  rowMeans(fitted1) - rowMeans(fitted0)
}

surfacesUptime("uptime before")

rows <- list()
for (movesetName in names(movesets)) {
  probs <- movesets[[movesetName]]
  for (designName in names(designs)) {
    # bcf runs under the shipped mixture alone: its treatment forest refuses
    # a non-default proposal.probs, so there is no paired arm to give it.
    armNames <- if (designName == "figure" && movesetName == movesetControl) {
      c("bart", "bcf")
    } else {
      "bart"
    }
    for (replicate in seq_len(nReplicates)) {
      set.seed(surfacesDataSeed("P6", designName, replicate))
      data <- drawShelf(designs[[designName]])
      counterfactual <- rbind(
        transform(data, z = 1L),
        transform(data, z = 0L)
      )
      for (armName in armNames) {
        startedAt <- proc.time()
        if (armName == "bart") {
          fit <- do.call(
            bart2,
            movesetCall(
              list(
                y ~ x1 + x2 + z,
                data = data,
                test = counterfactual,
                n.threads = 1L,
                verbose = FALSE,
                seed = surfacesSamplerSeed(replicate)
              ),
              probs
            )
          )
          surface <- extract(fit, type = "ev", sample = "test")
          ate <- ateFromContrast(
            surface[, seq_len(n), drop = FALSE],
            surface[, n + seq_len(n), drop = FALSE]
          )
        } else {
          propensityFit <- do.call(
            bart2,
            movesetCall(
              list(
                z ~ x1 + x2,
                data = data,
                family = "probit",
                n.threads = 1L,
                verbose = FALSE,
                seed = surfacesSamplerSeed(replicate)
              ),
              probs
            )
          )
          augmented <- data
          augmented$pihat <- fitted(propensityFit)
          fit <- do.call(
            bart2,
            movesetCall(
              list(
                y ~ x1 + x2 + pihat + z:forest(x1 + x2 + pihat),
                data = augmented,
                keepTrees = TRUE,
                n.threads = 1L,
                verbose = FALSE,
                seed = surfacesSamplerSeed(replicate)
              ),
              probs
            )
          )
          treated <- augmented
          treated$z <- 1L
          control <- augmented
          control$z <- 0L
          ate <- ateFromContrast(
            predict(fit, newdata = treated),
            predict(fit, newdata = control)
          )
        }
        elapsed <- (proc.time() - startedAt)[["elapsed"]]
        interval <- quantile(ate, c(0.025, 0.975), names = FALSE)
        rows[[length(rows) + 1L]] <- data.frame(
          moveset = movesetName,
          design = designName,
          arm = armName,
          replicate = replicate,
          estimate = mean(ate),
          lower = interval[1L],
          upper = interval[2L],
          covered = trueTau >= interval[1L] & trueTau <= interval[2L],
          width = interval[2L] - interval[1L],
          treatedFraction = mean(data$z),
          wall = elapsed,
          stringsAsFactors = FALSE
        )
      }
      if (replicate %% 25L == 0L) {
        cat(sprintf(
          "%-11s %-9s replicate %d of %d\n",
          movesetName,
          designName,
          replicate,
          nReplicates
        ))
      }
    }
  }
}
results <- do.call(rbind, rows)

surfacesHeader(sprintf(
  "P6 diagonal shelf: %d replications, true ATE %.0f",
  nReplicates,
  trueTau
))
cat(sprintf(
  "%-9s %-5s %-11s %8s %8s %8s %8s %8s\n",
  "design",
  "arm",
  "moveset",
  "bias",
  "|bias|",
  "coverage",
  "rmse",
  "length"
))
summaries <- list()
for (designName in names(designs)) {
  for (armName in unique(results$arm[results$design == designName])) {
    for (movesetName in names(movesets)) {
      keep <- results$design == designName &
        results$arm == armName &
        results$moveset == movesetName
      if (!any(keep)) {
        next
      }
      bias <- mean(results$estimate[keep]) - trueTau
      rmse <- sqrt(mean((results$estimate[keep] - trueTau)^2))
      coverage <- mean(results$covered[keep])
      summaries[[length(summaries) + 1L]] <- data.frame(
        design = designName,
        arm = armName,
        moveset = movesetName,
        bias = bias,
        coverage = coverage,
        rmse = rmse,
        length = mean(results$width[keep]),
        stringsAsFactors = FALSE
      )
      cat(sprintf(
        "%-9s %-5s %-11s %8.3f %8.3f %8.3f %8.3f %8.3f\n",
        designName,
        armName,
        movesetName,
        bias,
        abs(bias),
        coverage,
        rmse,
        mean(results$width[keep])
      ))
    }
  }
}
summaries <- do.call(rbind, summaries)

# The paired quantities. Bias is an aggregate over replications, so its
# paired signal is the per-replication ATE estimate, signed so that a
# positive difference is the one that enlarges the control's own |bias|;
# RMSE's is the per-replication squared error, whose ratio of means is the
# RMSE ratio squared.
contrastMovesets <- setdiff(names(movesets), movesetControl)
if (movesetControl %in% names(movesets) && length(contrastMovesets) > 0L) {
  surfacesHeader(sprintf(
    "P6 move sets: paired difference against %s, replication by replication",
    movesetControl
  ))
  cat(sprintf(
    "%-9s %-5s %-11s %-33s %-33s %s\n",
    "design",
    "arm",
    "moveset",
    "d ATE estimate",
    "d 95% coverage",
    "d squared error"
  ))
  verdicts <- list()
  for (designName in names(designs)) {
    for (armName in unique(results$arm[results$design == designName])) {
      base <- results[
        results$design == designName &
          results$arm == armName &
          results$moveset == movesetControl,
      ]
      if (nrow(base) == 0L) {
        next
      }
      baseBias <- mean(base$estimate) - trueTau
      baseSquared <- (base$estimate - trueTau)^2
      for (movesetName in contrastMovesets) {
        contrast <- results[
          results$design == designName &
            results$arm == armName &
            results$moveset == movesetName,
        ]
        paired <- match(base$replicate, contrast$replicate)
        keepBase <- !is.na(paired)
        contrast <- contrast[paired[keepBase], ]
        pairedBase <- base[keepBase, ]
        if (nrow(contrast) == 0L) {
          next
        }
        dEstimate <- contrast$estimate - pairedBase$estimate
        dCovered <- as.numeric(contrast$covered) -
          as.numeric(pairedBase$covered)
        contrastSquared <- (contrast$estimate - trueTau)^2
        dSquared <- contrastSquared - baseSquared[keepBase]
        cat(sprintf(
          "%-9s %-5s %-11s %-33s %-33s %s\n",
          designName,
          armName,
          movesetName,
          surfacesPairedDifference(dEstimate),
          surfacesPairedDifference(dCovered),
          surfacesPairedDifference(dSquared)
        ))
        biasRatio <- abs(mean(contrast$estimate) - trueTau) / abs(baseBias)
        rmseRatio <- sqrt(mean(contrastSquared) / mean(baseSquared[keepBase]))
        verdicts[[length(verdicts) + 1L]] <- data.frame(
          design = designName,
          arm = armName,
          moveset = movesetName,
          biasRatio = biasRatio,
          rmseRatio = rmseRatio,
          coverageDifference = mean(dCovered),
          coverageVerdict = surfacesMarginVerdict(dCovered, -0.010, -1),
          # A bias or error ratio's margin is 1.02, read on the point
          # estimate, with the paired quantity underneath it carrying
          # section 6.4's separation condition.
          biasSeparation = armSeparated(sign(baseBias) * dEstimate, 1),
          rmseSeparation = armSeparated(dSquared, 1),
          coverageImprovement = surfacesImprovementVerdict(dCovered, 1),
          biasImprovement = surfacesImprovementVerdict(
            sign(baseBias) * dEstimate,
            -1
          ),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  verdicts <- do.call(rbind, verdicts)

  surfacesHeader(
    "P6 move sets: margins (coverage -0.010, bias/RMSE ratio 1.02)"
  )
  cat(sprintf(
    "%-9s %-5s %-11s %-10s %-28s %-12s %-12s %s\n",
    "design",
    "arm",
    "moveset",
    "d coverage",
    "coverage verdict",
    "|bias| ratio",
    "rmse ratio",
    "ratio verdict"
  ))
  for (i in seq_len(nrow(verdicts))) {
    ratioVerdict <- if (
      verdicts$biasRatio[i] > 1.02 || verdicts$rmseRatio[i] > 1.02
    ) {
      sprintf(
        "past margin (%s bias, %s rmse)",
        verdicts$biasSeparation[i],
        verdicts$rmseSeparation[i]
      )
    } else {
      "within margin"
    }
    cat(sprintf(
      "%-9s %-5s %-11s %-10.3f %-28s %-12.3f %-12.3f %s\n",
      verdicts$design[i],
      verdicts$arm[i],
      verdicts$moveset[i],
      verdicts$coverageDifference[i],
      verdicts$coverageVerdict[i],
      verdicts$biasRatio[i],
      verdicts$rmseRatio[i],
      ratioVerdict
    ))
  }

  surfacesHeader("P6 move sets: the primary against the 4x-SE improvement bar")
  cat(sprintf(
    "%-9s %-5s %-11s %-30s %s\n",
    "design",
    "arm",
    "moveset",
    "coverage",
    "|bias|"
  ))
  for (i in seq_len(nrow(verdicts))) {
    cat(sprintf(
      "%-9s %-5s %-11s %-30s %s\n",
      verdicts$design[i],
      verdicts$arm[i],
      verdicts$moveset[i],
      verdicts$coverageImprovement[i],
      verdicts$biasImprovement[i]
    ))
  }
} else {
  verdicts <- NULL
}

surfacesHeader("published reference (Hahn, Murray and Carvalho 2020 Table 1)")
cat("BART  bias 0.27  coverage 0.65  rmse 0.31\n")
cat("BCF   bias 0.14  coverage 0.95  rmse 0.21\n")

surfacesUptime("uptime after")
surfacesSave(
  list(
    results = results,
    summaries = summaries,
    verdicts = verdicts,
    settings = list(
      nReplicates = nReplicates,
      n = n,
      trueTau = trueTau,
      designs = names(designs),
      movesets = movesets,
      movesetControl = movesetControl
    )
  ),
  outputDir,
  "P6-diagonal-shelf"
)
