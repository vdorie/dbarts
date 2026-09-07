#!/usr/bin/env Rscript

# P1, the low-noise Friedman emulator: the battery's known-positive control
# and its absolute gate. Pratola (2016) built this cell so that the tree
# structure freezes as the noise falls, and the battery's rule is that if it
# does not still freeze here, the harness is mismeasuring and no verdict from
# any other cell is valid.
#
# Two rungs, each an arm list run on the same matched seeds.
#
#   pratola  Pratola's own cell (arXiv 1312.1895 section 2.2): n = 5000,
#            p = 10, m = 200, 5000 burn-in plus 5000 kept, at sigma^2 = 1 and
#            sigma^2 = 0.1
#   house    the cheaper cell the move-set A/B measured, design 2 of
#            docs/design/tree-mixing-proposals.md section 13: n = 2000,
#            p = 10, m = 200, sigma = 0.25, 1000 burn-in plus 2000 kept
#
# Published numbers to reproduce, under birth/death proposals only:
#
#   sigma^2 = 1     acceptance rate around 18%, 90% coverage of eta 81%
#   sigma^2 = 0.1   acceptance rate around 4%,  90% coverage of eta 53.8%
#
# Both are the coverage of the 90% pointwise credible interval for eta(x).
# The section reports them without qualification and its figure plots the
# intervals against the fitted 5000 settings, so they are read here as
# in-sample; the paper's section 6 is where an out-of-sample coverage is
# named as such. The response is eta plus noise, not the deterministic
# simulator output: "we treated the deterministic Friedman function as if it
# were our simulator, i.e. y(x) = eta(x) + eps(x)".
#
# The measured house-rung reference is 90% coverage 0.714, 0.725 and 0.714
# under the three mixtures section 13 ran, on held-out rows. That is the
# absolute gate: the default arm must sit near 0.71.
#
# Primary statistic: 90% pointwise coverage of the true f, at the training
# points and on 1000 held-out rows. Secondaries: RMSE of the posterior mean
# against true f, and minimum ESS over 25 evenly spaced held-out rows. The
# structural acceptance readout is a separate script, P1-friedman-census.R,
# because it needs an instrumented build.
#
# Arms, all reachable at runtime through proposal.probs:
#
#   default     the shipped mixture, left unset
#   birthdeath  birth_death 1.0, swap 0.0, change 0.0, birth 0.5, Pratola's
#               own arm and the one his numbers were taken on
#   swap        birth_death 0.5, swap 0.1, change 0.4, birth 0.5, the
#               historical mixture
#
# Pratola's rung carries a fourth design, run on the birth/death arm alone:
# the same sigma^2 = 0.1 cell on Friedman's own function rather than the
# doubled frequency his paper prints, so the transcription ambiguity is
# measured rather than assumed.
#
# Usage: Rscript P1-friedman.R [outputDir] [quick] [house] [pratola]

source(
  file.path(
    dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))),
    "surfaces-common.R"
  ),
  chdir = FALSE
)

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args

rungNames <- c("house", "pratola")
selectedRungs <- intersect(rungNames, args)
if (length(selectedRungs) == 0L) {
  selectedRungs <- rungNames
}
outputDir <- surfacesOutputDir(args, flags = c("quick", rungNames))

nReplicates <- if (quick) 3L else 20L
nTest <- 1000L
nTrees <- 200L
p <- 10L
level <- 0.9

# NULL leaves proposal.probs unset, which is what the default arm is.
arms <- list(
  default = NULL,
  birthdeath = c(birth_death = 1, swap = 0, change = 0, birth = 0.5),
  swap = c(birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5)
)

rungs <- list(
  house = list(
    n = if (quick) 500L else 2000L,
    nBurn = if (quick) 250L else 1000L,
    nSamples = if (quick) 250L else 2000L,
    designs = list(
      sigma025 = list(sigma = 0.25, frequency = 1, arms = names(arms))
    )
  ),
  pratola = list(
    n = if (quick) 500L else 5000L,
    nBurn = if (quick) 250L else 5000L,
    nSamples = if (quick) 250L else 5000L,
    designs = list(
      variance1 = list(sigma = 1, frequency = 2, arms = names(arms)),
      variance01 = list(sigma = sqrt(0.1), frequency = 2, arms = names(arms)),
      variance01friedman = list(
        sigma = sqrt(0.1),
        frequency = 1,
        arms = "birthdeath"
      )
    )
  )
)

published <- list(
  house = "section 13 design 2, held-out 90% coverage 0.714 / 0.725 / 0.714",
  pratola = paste(
    "Pratola 2016 sec 2.2, birth/death only:",
    "sigma^2 = 1 acceptance ~18%, coverage 0.81;",
    "sigma^2 = 0.1 acceptance ~4%, coverage 0.538"
  )
)

# 25 evenly spaced held-out rows carry the ESS, as the move-set grid does.
essPoints <- as.integer(round(seq(1, nTest, length.out = 25L)))

runRung <- function(rungName) {
  rung <- rungs[[rungName]]
  rows <- list()
  for (designName in names(rung$designs)) {
    design <- rung$designs[[designName]]
    for (armName in design$arms) {
      for (replicate in seq_len(nReplicates)) {
        set.seed(surfacesDataSeed(
          "P1",
          paste0(rungName, designName),
          replicate
        ))
        data <- surfacesFriedman(
          rung$n,
          nTest,
          p,
          sigma = design$sigma,
          frequency = design$frequency
        )
        call <- list(
          data$x,
          data$y,
          test = data$xTest,
          n.trees = nTrees,
          n.chains = 1L,
          n.burn = rung$nBurn,
          n.samples = rung$nSamples,
          n.thin = 1L,
          n.threads = 1L,
          verbose = FALSE,
          seed = surfacesSamplerSeed(replicate)
        )
        if (!is.null(arms[[armName]])) {
          call$proposal.probs <- arms[[armName]]
        }
        startedAt <- proc.time()
        fit <- do.call(bart2, call)
        elapsed <- (proc.time() - startedAt)[["elapsed"]]
        trainDraws <- extract(fit, type = "ev", sample = "train")
        testDraws <- extract(fit, type = "ev", sample = "test")
        trainInterval <- apply(
          trainDraws,
          2L,
          quantile,
          probs = c((1 - level) / 2, 1 - (1 - level) / 2),
          names = FALSE
        )
        rows[[length(rows) + 1L]] <- data.frame(
          rung = rungName,
          design = designName,
          arm = armName,
          replicate = replicate,
          coverageTrain = mean(
            data$f >= trainInterval[1L, ] & data$f <= trainInterval[2L, ]
          ),
          lengthTrain = mean(trainInterval[2L, ] - trainInterval[1L, ]),
          rmseTrain = surfacesRmse(trainDraws, data$f),
          coverageTest = surfacesCoverage(testDraws, data$fTest, level = level),
          rmseTest = surfacesRmse(testDraws, data$fTest),
          minEss = min(surfacesPointEss(testDraws, essPoints)),
          sigmaTruth = design$sigma,
          sigmaPosterior = mean(fit$sigma),
          wall = elapsed,
          stringsAsFactors = FALSE
        )
        cat(sprintf(
          "%-18s %-11s rep %2d  train %.3f  test %.3f  rmse %.3f  %.0fs\n",
          designName,
          armName,
          replicate,
          rows[[length(rows)]]$coverageTrain,
          rows[[length(rows)]]$coverageTest,
          rows[[length(rows)]]$rmseTest,
          elapsed
        ))
        rm(fit, trainDraws, testDraws)
        invisible(gc(FALSE))
      }
    }
  }
  do.call(rbind, rows)
}

for (rungName in selectedRungs) {
  surfacesUptime(sprintf("uptime before %s", rungName))
  results <- runRung(rungName)

  surfacesHeader(sprintf(
    "P1 Friedman, %s rung: mean over seeds (min-max)",
    rungName
  ))
  cat(sprintf(
    "%-18s %-11s %-20s %-20s %-20s %-20s %-12s %s\n",
    "design",
    "arm",
    "90% cover train",
    "90% cover held-out",
    "held-out RMSE",
    "interval length",
    "min ESS",
    "wall s"
  ))
  rung <- rungs[[rungName]]
  for (designName in names(rung$designs)) {
    for (armName in rung$designs[[designName]]$arms) {
      keep <- results$design == designName & results$arm == armName
      cat(sprintf(
        "%-18s %-11s %-20s %-20s %-20s %-20s %-12s %.0f\n",
        designName,
        armName,
        surfacesRange(results$coverageTrain[keep]),
        surfacesRange(results$coverageTest[keep]),
        surfacesRange(results$rmseTest[keep]),
        surfacesRange(results$lengthTrain[keep], digits = 2L),
        surfacesRange(results$minEss[keep], digits = 0L),
        mean(results$wall[keep])
      ))
    }
  }

  surfacesHeader("published reference")
  cat(published[[rungName]], "\n", sep = "")

  surfacesUptime(sprintf("uptime after %s", rungName))
  surfacesSave(
    list(
      results = results,
      published = published[[rungName]],
      settings = list(
        rung = rungName,
        nReplicates = nReplicates,
        nTest = nTest,
        p = p,
        nTrees = nTrees,
        level = level,
        n = rung$n,
        nBurn = rung$nBurn,
        nSamples = rung$nSamples,
        designs = rung$designs,
        arms = arms,
        essPoints = essPoints
      )
    ),
    outputDir,
    paste0("P1-friedman-", rungName)
  )
}
