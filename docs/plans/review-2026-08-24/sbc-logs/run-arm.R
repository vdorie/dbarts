#!/usr/bin/env Rscript
# Ensemble-scale SBC driver: reuses benchmarks/R/sbc.R's machinery unchanged
# (source() never reaches its CLI block) and only re-parameterises the config
# to the SHIPPED default forest size. Writes the rank matrix as an RDS so the
# bands can be recomputed (per-functional and Bonferroni) after the fact.
args <- commandArgs(trailingOnly = TRUE)
arm <- args[1]
R <- as.integer(args[2])
L <- as.integer(args[3])
thin <- as.integer(args[4])
burnSweeps <- if (length(args) >= 5L && nzchar(args[5])) as.numeric(args[5]) else NULL
nTrees <- if (length(args) >= 6L) as.integer(args[6]) else 200L
src <- Sys.getenv("SBC_SRC")
out <- Sys.getenv("SBC_OUT")
tag <- Sys.getenv("SBC_TAG", arm)
source(file.path(src, "benchmarks", "R", "sbc.R"))

cfg <- switch(arm,
  gaussian = sbcConfig(family = "gaussian", nTrees = nTrees),
  probit = sbcConfig(family = "probit", nTrees = nTrees),
  logistic = sbcConfig(family = "logistic", nTrees = nTrees),
  weighted = {
    c0 <- sbcConfig(family = "gaussian", nTrees = nTrees)
    set.seed(7L)
    c0$weights <- rgamma(c0$n, 2, 2)
    c0
  },
  t = sbcConfig(family = "t", nTrees = nTrees),
  ordinal = sbcConfig(family = "ordinal", numCategories = 4L, nTest = 3L, nTrees = nTrees),
  nbinom = sbcConfig(family = "nbinom", k = 8, nTrees = nTrees),
  multinom = sbcConfig(family = "multinomial", numCategories = 3L, nTrees = nTrees),
  `grouped-gaussian` = sbcAddGrouping(
    sbcConfig(family = "gaussian", n = 160L, nTrees = nTrees),
    nGroups = 8L, relScale = 0.2, zeroWeightFrac = 0.2
  ),
  `grouped-probit` = sbcAddGrouping(
    sbcConfig(family = "probit", n = 160L, nTrees = nTrees),
    nGroups = 8L, relScale = 0.2, zeroWeightFrac = 0
  ),
  bcf = sbcAddBCF(sbcConfig(family = "gaussian", n = 200L, nTrees = nTrees)),
  stop("unknown arm ", arm)
)

cat(sprintf("== arm %s  nTrees=%d n=%d R=%d L=%d thin=%d ==\n", arm, cfg$nTrees, cfg$n, R, L, thin))
ok <- c(sigma = isTRUE(sbcCheckSigmaPrior(cfg$sigest, cfg$sigDf, cfg$sigQuant)$pass))
isTier <- arm %in% c("t", "ordinal", "nbinom", "multinom")
if (isTier) {
  if (cfg$family == "ordinal") ok["cutpoints"] <- isTRUE(sbcCheckOrdinalCutpointPrior()$pass)
  if (cfg$family %in% c("nbinom", "t")) {
    ok["grid"] <- isTRUE(sbcCheckGridPrior(if (cfg$family == "nbinom") sbcNbGrid else sbcTGrid)$pass)
  }
  if (cfg$family == "multinomial") {
    m <- sbcCheckMultinomialProbs(cfg)
    cat(sprintf("  softmax vs reported: %.2e\n", m$maxDiff))
    ok["softmax"] <- isTRUE(m$pass)
  } else {
    lc <- sbcCheckLatentConsistency(cfg)
    cat(sprintf("  predict vs recorded latent: train %.2e test %.2e\n", lc$maxDiff, lc$maxDiffTest))
    ok["latent"] <- isTRUE(lc$pass)
  }
}
if (arm %in% c("grouped-gaussian", "grouped-probit")) {
  ok["tau"] <- isTRUE(sbcCheckTau(cfg$relScale)$pass)
}
if (arm == "bcf") ok["glue"] <- isTRUE(sbcCheckBCFGlue(cfg$sdControl, cfg$bPriorVariance)$pass)
cat("  self-checks: ", paste(names(ok), ifelse(ok, "PASS", "FAIL"), collapse = " "), "\n", sep = "")
if (any(!ok)) stop("harness self-check failed at nTrees=", nTrees)

fit <- if (isTier) {
  if (is.null(burnSweeps)) {
    runSbcFamily(cfg, R = R, L = L, thin = thin)
  } else {
    runSbcFamily(cfg, R = R, L = L, thin = thin, burnSweeps = burnSweeps)
  }
} else if (arm %in% c("grouped-gaussian", "grouped-probit")) {
  if (is.null(burnSweeps)) {
    runSbcGrouped(cfg, R = R, L = L, thin = thin)
  } else {
    runSbcGrouped(cfg, R = R, L = L, thin = thin, burn = as.integer(ceiling(burnSweeps / thin)))
  }
} else if (arm == "bcf") {
  if (is.null(burnSweeps)) {
    runSbcBCF(cfg, R = R, L = L, thin = thin)
  } else {
    runSbcBCF(cfg, R = R, L = L, thin = thin, burn = as.integer(ceiling(burnSweeps / thin)))
  }
} else {
  if (is.null(burnSweeps)) {
    runSbc(cfg, R = R, L = L, thin = thin)
  } else {
    runSbc(cfg, R = R, L = L, thin = thin, burn = as.integer(ceiling(burnSweeps / thin)))
  }
}
saveRDS(fit, file.path(out, paste0("fit-", tag, ".rds")))
sbcReport(fit, alpha = 0.05)
cat("\nDONE ", tag, " elapsed ", sprintf("%.1f", fit$elapsed), "s\n", sep = "")
