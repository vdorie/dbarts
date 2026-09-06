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
# Arms:
#   bart   the paper's standard BART arm: one forest over (x1, x2, z), ATE
#          read off the counterfactual contrast
#   bcf    the paper's fix, in dbarts' own shape: an estimated propensity
#          added as a covariate and a second forest whose amplitude the
#          treatment indicator modulates
#
# Usage: Rscript P6-diagonal-shelf.R [outputDir] [quick]

source(
  file.path(
    dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))),
    "surfaces-common.R"
  ),
  chdir = FALSE
)

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args
outputDir <- surfacesOutputDir(args, flags = "quick")

nReplicates <- if (quick) 5L else 200L
n <- 250L
trueTau <- -1

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

# ATE draws from a counterfactual contrast: the posterior mean surface
# evaluated at every row twice, once with z = 1 and once with z = 0.
ateFromContrast <- function(fitted1, fitted0) {
  rowMeans(fitted1) - rowMeans(fitted0)
}

surfacesUptime("uptime before")

rows <- list()
for (designName in names(designs)) {
  armNames <- if (designName == "figure") c("bart", "bcf") else "bart"
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
        fit <- bart2(
          y ~ x1 + x2 + z,
          data = data,
          test = counterfactual,
          n.threads = 1L,
          verbose = FALSE,
          seed = surfacesSamplerSeed(replicate)
        )
        surface <- extract(fit, type = "ev", sample = "test")
        ate <- ateFromContrast(
          surface[, seq_len(n), drop = FALSE],
          surface[, n + seq_len(n), drop = FALSE]
        )
      } else {
        propensityFit <- bart2(
          z ~ x1 + x2,
          data = data,
          family = "probit",
          n.threads = 1L,
          verbose = FALSE,
          seed = surfacesSamplerSeed(replicate)
        )
        augmented <- data
        augmented$pihat <- fitted(propensityFit)
        fit <- bart2(
          y ~ x1 + x2 + pihat + z:forest(x1 + x2 + pihat),
          data = augmented,
          keepTrees = TRUE,
          n.threads = 1L,
          verbose = FALSE,
          seed = surfacesSamplerSeed(replicate)
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
        "%-9s replicate %d of %d\n",
        designName,
        replicate,
        nReplicates
      ))
    }
  }
}
results <- do.call(rbind, rows)

surfacesHeader("P6 diagonal shelf: 200 replications, true ATE -1")
cat(sprintf(
  "%-9s %-5s %8s %8s %8s %8s %8s\n",
  "design",
  "arm",
  "bias",
  "|bias|",
  "coverage",
  "rmse",
  "length"
))
summaries <- list()
for (designName in names(designs)) {
  for (armName in unique(results$arm[results$design == designName])) {
    keep <- results$design == designName & results$arm == armName
    bias <- mean(results$estimate[keep]) - trueTau
    rmse <- sqrt(mean((results$estimate[keep] - trueTau)^2))
    coverage <- mean(results$covered[keep])
    summaries[[length(summaries) + 1L]] <- data.frame(
      design = designName,
      arm = armName,
      bias = bias,
      coverage = coverage,
      rmse = rmse,
      length = mean(results$width[keep]),
      stringsAsFactors = FALSE
    )
    cat(sprintf(
      "%-9s %-5s %8.3f %8.3f %8.3f %8.3f %8.3f\n",
      designName,
      armName,
      bias,
      abs(bias),
      coverage,
      rmse,
      mean(results$width[keep])
    ))
  }
}
summaries <- do.call(rbind, summaries)

surfacesHeader("published reference (Hahn, Murray and Carvalho 2020 Table 1)")
cat("BART  bias 0.27  coverage 0.65  rmse 0.31\n")
cat("BCF   bias 0.14  coverage 0.95  rmse 0.21\n")

surfacesUptime("uptime after")
surfacesSave(
  list(
    results = results,
    summaries = summaries,
    settings = list(
      nReplicates = nReplicates,
      n = n,
      trueTau = trueTau,
      designs = names(designs)
    )
  ),
  outputDir,
  "P6-diagonal-shelf"
)
