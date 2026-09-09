# A DONOR warm start on a multi-forest sampler. installForests reassembles
# each chain's trees from a saved slot but takes the amplitudes off the donor's
# LIVE state, so at more than one forest it would pair one draw's forests with
# another's glue and answer with a legal-looking, miscalibrated fit; nothing
# covers that result, so every entrance refuses. One assertion per entrance -
# the modelling argument, the R5 method, the bridge entry - plus the shapes the
# guard must NOT catch. The refusal is keyed on the forest COUNT, so the
# K-forest softmax is pinned beside the two-forest amplitude model.
#
# Grow-from-root is deliberately NOT refused: it composes through the combiner
# inside its own sweep and its two-forest branch is covered from R
# (test-prior-init-composed-law.R) and in tests/cpp, so the arms below pin that
# it still runs at two forests rather than that it errors.

source(
  system.file("common", "captureWarnings.R", package = "dbarts"),
  local = TRUE
)

set.seed(3)
n <- 120L
p <- 3L
x <- matrix(runif(n * p), n, p)
colnames(x) <- paste0("x", seq_len(p))
z <- rbinom(n, 1L, 0.5)
y <- 2 * sin(pi * x[, 1L]) + z * (1 + x[, 2L]) + rnorm(n, sd = 0.2)
zBasis <- cbind(1 - as.double(z), as.double(z))

controlWarmStartRefusal <- function(...) {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 20L,
    n.samples = 5L,
    updateState = FALSE,
    ...
  )
}

# the warm-start donor every arm below offers: an ordinary single-forest fit
# over the same predictors, so nothing but the forest count can be the refusal
donor <- dbarts(x, y, control = controlWarmStartRefusal(keepTrees = TRUE))
donor$sampleTreesFromPrior()
invisible(donor$run(0L, 3L))

bcf <- dbarts(
  x,
  y,
  forests = list(forest(), forest(basis = ~ factor(z))),
  control = controlWarmStartRefusal()
)

# --- the R5 surface: the method, the shape and the count a caller can read
# back for itself ---
expect_error(
  bcf$installTrees(donor),
  "installTrees\\(\\) does not support a multi-forest sampler: this one carries 2 forests"
)

# the donor is not even resolved before the refusal fires, so a caller offering
# a donor this sampler could never take reads the forest count rather than a
# donor complaint
expect_error(
  bcf$installTrees("not a donor"),
  "does not support a multi-forest sampler"
)

# --- the bridge backstop, for the routes that skip the R5 layer: the same
# rule under the C entry point's own name ---
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_installForests,
    bcf$getPointer(),
    donor$state,
    NULL
  ),
  "bartcore_installForests: a multi-forest sampler \\(2 forests\\) has no tested warm start from a donor"
)

# --- the modelling surface: a data object carrying forest bases reaches bart
# as an ordinary fit, so the donor argument refuses there by the name the
# caller wrote ---
basesData <- dbartsData(x, y, bases = list(NULL, zBasis))
expect_error(
  bart(
    basesData,
    warm.start = donor,
    n.samples = 5L,
    n.burn = 1L,
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 20L,
    verbose = FALSE
  ),
  "'warm.start' does not support a multi-forest sampler: this one carries 2 forests"
)

# ... while the OTHER initialization argument is unaffected on the same object:
# grow-from-root starts a two-forest fit and the fit comes back finite
grownFit <- bart(
  basesData,
  n.grow.sweeps = 2L,
  n.samples = 5L,
  n.burn = 1L,
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 20L,
  verbose = FALSE
)
expect_true(all(is.finite(grownFit$yhat.train)))
expect_equal(grownFit$n.forests, 2L)
expect_silent(bcf$growFromRoot(2L))

# --- the count is the key, not the amplitude coupling: a K-forest softmax
# carries no amplitudes and is refused all the same, at its own K ---
labels <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
softmax <- dbarts(
  x,
  labels,
  family = "multinomial",
  control = controlWarmStartRefusal()
)
expect_error(
  softmax$installTrees(donor),
  "installTrees\\(\\) does not support a multi-forest sampler: this one carries 3 forests"
)
expect_silent(softmax$growFromRoot(2L))

# --- and the shapes the guard must not catch. A plain single-forest sampler
# takes a donor, silently, as it always has ---
single <- dbarts(x, y, control = controlWarmStartRefusal())
expect_silent(single$installTrees(donor))

# a heteroscedastic sampler's variance forest is not one of the forests the
# count reports - the mean forest is the only one - so it keeps the donor route
# too, and the refusal above cannot be read as a variance-forest rule
heteroDonor <- dbarts(
  x,
  y,
  variance = TRUE,
  control = controlWarmStartRefusal(keepTrees = TRUE)
)
heteroDonor$sampleTreesFromPrior()
invisible(heteroDonor$run(0L, 3L))
hetero <- dbarts(x, y, variance = TRUE, control = controlWarmStartRefusal())
expect_silent(hetero$installTrees(heteroDonor))

# --- the refusal raises and warns about nothing: an entrance that warned and
# then went on would be the failure this file exists to catch ---
expect_equal(
  length(captureWarnings(try(bcf$installTrees(donor), silent = TRUE))),
  0L
)
