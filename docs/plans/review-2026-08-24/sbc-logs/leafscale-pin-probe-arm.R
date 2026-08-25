suppressMessages(library(dbarts))
args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]; R <- as.integer(args[2]); nTrees <- as.integer(args[3])
src <- Sys.getenv("SBC_SRC"); out <- Sys.getenv("SBC_OUT")
source(file.path(src, "benchmarks", "R", "sbc.R"))
cfg <- sbcAddGrouping(sbcConfig(family = "gaussian", n = 160L, nTrees = nTrees),
                      nGroups = 8L, relScale = 0.2, zeroWeightFrac = 0.2)
if (mode == "pinned") {
  # OPTION B, the whole harness change: name the anchor calibration so every
  # rebuilt fit carries the build vector's leaf prior instead of range(y0)'s.
  anchor <- dbarts(cfg$x, cfg$yBuild,
                   control = dbartsControl(n.trees = cfg$nTrees, n.chains = 1L,
                                           n.threads = 1L, updateState = FALSE,
                                           verbose = FALSE))
  S <- anchor$getCalibration(1L)[1L, "prior.scale"]
  cat("pinned prior.scale =", S, "\n")
  cfg$nodePrior <- dbartsPriors$normal(cfg$k, scale = S)
}
stopifnot(isTRUE(sbcCheckSigmaPrior(cfg$sigest, cfg$sigDf, cfg$sigQuant)$pass),
          isTRUE(sbcCheckTau(cfg$relScale)$pass))
fit <- runSbcGrouped(cfg, R = R, L = 150L, thin = 30L, report = 25L, swap = identical(mode, "swap"))
saveRDS(fit, file.path(out, paste0("grouped-", mode, "-m", nTrees, "-R", R, ".rds")))
sbcReport(fit, alpha = 0.05)
cat("\nDONE", mode, sprintf("%.1f", fit$elapsed), "s\n")
