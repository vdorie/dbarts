# Packaging the per-forest in-sample channels (forestFits, glue) and
# extract(type = "forest") on top of the sampler-level reporting
# test-bcf-reporting.R already exercises.

set.seed(311)
n <- 50L
x <- matrix(runif(n * 3L), n, 3L)
z <- rbinom(n, 1L, 0.5)
y <- 2 * sin(pi * x[, 1L]) + z * (1 + x[, 2L]) + rnorm(n, sd = 0.2)

fastControl <- function(n.chains = 1L) {
  dbartsControl(
    n.threads = 1L,
    n.trees = 8L,
    n.chains = n.chains,
    n.samples = 6L,
    n.burn = 2L,
    updateState = FALSE,
    seed = 311L
  )
}

packageFrom <- function(forests, n.chains = 1L, combineChains = TRUE) {
  sampler <- dbarts::dbarts(
    x,
    y,
    forests = forests,
    control = fastControl(n.chains)
  )
  burn <- dbarts:::runWithBurnIn(sampler, sampler$control, FALSE)
  dbarts:::packageBartResults(
    sampler,
    burn$samples,
    burn$burnInSigma,
    burn$burnInK,
    combineChains,
    FALSE
  )
}

twoForests <- list(forest(), forest(basis = ~ factor(z)))

# --- shapes and dimnames, combineChains = TRUE ---
fit <- packageFrom(twoForests)
expect_equal(dim(fit$forestFits), c(6L, n, 2L))
expect_identical(dimnames(fit$forestFits)[[3L]], c("forest1", "forest2"))
expect_equal(dim(fit$glue), c(6L, 3L))
expect_identical(attr(fit$glue, "forest"), c("forest1", "forest2", "forest2"))
expect_equal(fit$n.forests, 2L)
expect_null(fit$bases[[1L]])
expect_equal(dim(fit$bases[[2L]]), c(n, 2L))
expect_null(attr(fit, "forest.labels"))

# --- forest.labels rides the declaration's own names ---
namedFit <- packageFrom(list(
  prognostic = forest(),
  treatment = forest(basis = ~ factor(z))
))
expect_identical(attr(namedFit, "forest.labels"), c("prognostic", "treatment"))

# --- combineChains = FALSE keeps the per-chain axis ---
fitUncombined <- packageFrom(twoForests, n.chains = 2L, combineChains = FALSE)
expect_equal(dim(fitUncombined$forestFits), c(2L, 6L, n, 2L))
expect_equal(dim(fitUncombined$glue), c(2L, 6L, 3L))
expect_identical(
  attr(fitUncombined$glue, "forest"),
  c("forest1", "forest2", "forest2")
)

# --- extract(type = "forest"): default slice, by index/name, contribution ---
expect_identical(extract(fit, type = "forest"), fit$forestFits)

byIndex <- extract(fit, type = "forest", forest = 2L)
byName <- extract(fit, type = "forest", forest = "forest2")
expect_identical(byIndex, byName)
expect_equal(dim(byIndex), c(6L, n, 1L))
expect_identical(dimnames(byIndex)[[3L]], "forest2")

# contribution = (basis %*% glue) * forestFits, computed on demand rather
# than stored; re-derived here from the packaged elements alone, one
# draw, both forests
contrib <- extract(fit, type = "forest", contribution = TRUE)
expect_equal(dim(contrib), dim(fit$forestFits))
glueForest <- attr(fit$glue, "forest")
expected1 <- fit$glue[1L, glueForest == "forest1"] * fit$forestFits[1L, , 1L]
expected2 <- as.vector(
  fit$bases[[2L]] %*% fit$glue[1L, glueForest == "forest2"]
) *
  fit$forestFits[1L, , 2L]
expect_equal(contrib[1L, , 1L], expected1)
expect_equal(contrib[1L, , 2L], expected2)

# combineChains = FALSE is a no-op at one chain, as it is for every other
# channel - there is no separate chain axis to split
expect_identical(
  extract(fit, type = "forest", combineChains = FALSE),
  fit$forestFits
)

# extract(combineChains = TRUE) on an uncombined-at-packaging-time (two
# chains) fit reshapes back, bitwise, to the same values a directly-combined
# packaging would carry
recombined <- extract(fitUncombined, type = "forest", combineChains = TRUE)
directFit <- packageFrom(twoForests, n.chains = 2L, combineChains = TRUE)
expect_identical(unname(recombined), unname(directFit$forestFits))

# --- refusals, both by name ---
expect_error(
  extract(fit, type = "forest", sample = "test"),
  pattern = "sample = \"test\""
)
plainFit <- packageFrom(NULL)
expect_error(
  extract(plainFit, type = "forest"),
  pattern = "per-forest reporting"
)
