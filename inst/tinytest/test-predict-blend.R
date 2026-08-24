# The combined arm at new rows on an amplitude-coupled fit: predict(type =
# "ev"/"ppd"/"bart") performs the recombination predict(type = "forest")
# leaves to the caller, from the per-forest replay, the fit's own glue, and
# the bases at those rows. The raw per-forest arm itself is
# test-predict-forest.R; what is pinned here is the BLEND - that it reproduces
# the fit in sample on every family and shape that can carry amplitudes, that
# it moves with the basis it is given, that a forest() term's basis is
# re-derived from newdata, and what it refuses.

set.seed(811)
n <- 60L
p <- 3L
x <- matrix(runif(n * p), n, p)
colnames(x) <- paste0("x", seq_len(p))
z <- rbinom(n, 1L, 0.5)
# forest 1's own basis in the both-bases fixture below
w <- runif(n)
g3 <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
y <- 2 * sin(pi * x[, 1L]) + z * (1 + 2 * x[, 2L]) + rnorm(n, sd = 0.2)
yBin <- rbinom(n, 1L, pnorm(-0.4 + 1.5 * x[, 1L] + 0.8 * z))

nNew <- 25L
xNew <- matrix(runif(nNew * p), nNew, p)
colnames(xNew) <- colnames(x)

n.samples <- 6L

# keepTrees is installed for the SAMPLING run only, after burn-in, as bart()
# installs it: the saved-tree store then holds exactly the recorded draws, in
# order, which is what pairs each draw's forests with its own amplitudes below.
keptControl <- function(n.chains = 1L) {
  dbartsControl(
    n.threads = 1L,
    n.trees = 8L,
    n.chains = n.chains,
    n.samples = n.samples,
    n.burn = 2L,
    keepTrees = FALSE,
    updateState = FALSE,
    seed = 811L
  )
}

fitFrom <- function(
  forests,
  response = y,
  family = "auto",
  n.chains = 1L,
  combineChains = TRUE,
  keepTrees = TRUE
) {
  sampler <- dbarts(
    x,
    response,
    forests = forests,
    family = family,
    control = keptControl(n.chains)
  )
  burn <- dbarts:::runWithBurnIn(sampler, sampler$control, keepTrees)
  dbarts:::packageBartResults(
    sampler,
    burn$samples,
    burn$burnInSigma,
    burn$burnInK,
    combineChains,
    TRUE
  )
}

# ragged widths (1, 2): forest 1's implicit all-ones column beside a two-level
# basis, the shape that tells reading the glue by its "forest" attribute apart
# from reading it positionally
twoForests <- list(forest(), forest(basis = ~ factor(z)))
fit <- fitFrom(twoForests)

# --- THE PRIMARY ORACLE. At the TRAINING rows, with the fit's OWN bases, the
# blend is the fit. yhat.train is what the engine's own combiner recorded, so
# this is independent of any R-side recombination. Not bitwise, and the gap is
# accumulation ORDER alone: the engine sums its forests last-first and applies
# the shift as fitScale * f + fitShift, where the blend sums forest 1 first
# into a shift-filled matrix. Measured 1.8e-15 here, three orders under the
# 1e-12 bar the reconstruction identities elsewhere use.
blendAtTraining <- function(fitted, ...) {
  predict(fitted, x, bases = fitted$bases, ...)
}
expect_true(
  max(abs(blendAtTraining(fit, type = "bart") - fit$yhat.train)) < 1e-12
)
# the gaussian fixture is the only one whose shift and scale are both nonzero
# and non-unit, which is what makes it - and not the binary arms - the cell
# that would catch a dropped shift or a doubled response scale
shift <- fit$fit$getCalibration(1L)[1L, "response.shift"]
expect_true(abs(shift) > 1e-8)
expect_true(abs(fit$fit$getCalibration(1L)[1L, "response.scale"] - 1) > 1e-8)

# the same oracle on every other shape that can carry amplitudes. Under a
# latent family the response transform is the identity, so "ev" is the link
# applied to the recorded latent draws.
probitFit <- fitFrom(twoForests, yBin, "probit")
expect_true(
  max(abs(
    blendAtTraining(probitFit, type = "ev") - pnorm(probitFit$yhat.train)
  )) <
    1e-12
)
logisticFit <- fitFrom(twoForests, yBin, "logistic")
expect_true(
  max(abs(
    blendAtTraining(logisticFit, type = "ev") - plogis(logisticFit$yhat.train)
  )) <
    1e-12
)
# two chains stored uncombined: glue and the replay pair draw for draw only in
# the combined layout, so the split result is a reshape of that one
chainFit <- fitFrom(twoForests, n.chains = 2L, combineChains = FALSE)
expect_equal(dim(chainFit$yhat.train), c(2L, n.samples, n))
expect_true(
  max(abs(
    blendAtTraining(chainFit, type = "bart", combineChains = FALSE) -
      chainFit$yhat.train
  )) <
    1e-12
)
expect_identical(
  blendAtTraining(chainFit, type = "bart", combineChains = FALSE),
  dbarts:::uncombineChains(blendAtTraining(chainFit, type = "bart"), 2L)
)
# a three-column basis, so the ragged margin is (1, 3)
wideFit <- fitFrom(list(forest(), forest(basis = ~g3)))
expect_identical(dim(wideFit$bases[[2L]]), c(n, 3L))
expect_true(
  max(abs(blendAtTraining(wideFit, type = "bart") - wideFit$yhat.train)) < 1e-12
)
# forest 1 carrying a basis of its own rather than the implicit column
bothBasesFit <- fitFrom(list(forest(basis = ~w), forest(basis = ~ factor(z))))
expect_identical(dim(bothBasesFit$bases[[1L]]), c(n, 1L))
expect_true(
  max(abs(
    blendAtTraining(bothBasesFit, type = "bart") - bothBasesFit$yhat.train
  )) <
    1e-12
)

# --- SECONDARY: the manual recombination man/bart.Rd documents as the fallback
# (test-predict-forest.R's own closure, verbatim) is what the blend performs,
# bitwise, on BOTH counterfactual arms. The same loop written twice, so what it
# pins is accumulation order; the oracle above is what pins the arithmetic.
glueForest <- attr(fit$glue, "forest")
predicted <- predict(fit, xNew, type = "forest")
recombine <- function(perForest, bases, nRows) {
  out <- matrix(shift, nrow(fit$glue), nRows)
  for (k in seq_len(fit$n.forests)) {
    basis <- if (is.null(bases[[k]])) matrix(1, nRows, 1L) else bases[[k]]
    g <- fit$glue[, glueForest == dimnames(perForest)[[3L]][k], drop = FALSE]
    out <- out + (g %*% t(basis)) * perForest[,, k]
  }
  out
}
zeroBasis <- list(NULL, unname(cbind(rep(1, nNew), rep(0, nNew))))
oneBasis <- list(NULL, unname(cbind(rep(0, nNew), rep(1, nNew))))
atZero <- predict(fit, xNew, type = "bart", bases = zeroBasis)
atOne <- predict(fit, xNew, type = "bart", bases = oneBasis)
expect_equal(dim(atZero), c(n.samples, nNew))
expect_identical(atZero, recombine(predicted, zeroBasis, nNew))
expect_identical(atOne, recombine(predicted, oneBasis, nNew))

# the arms differ by exactly forest 2's own total moved from its first
# amplitude to its second - the whole point of a basis, and what an
# implementation that ignored the basis would collapse to zero
amplitudes <- fit$glue[, glueForest == "forest2", drop = FALSE]
expect_true(max(abs(atOne - atZero)) > 1e-6)
expect_true(
  max(abs(
    (atOne - atZero) - (amplitudes[, 2L] - amplitudes[, 1L]) * predicted[,, 2L]
  )) <
    1e-12
)

# a bare value positions itself when exactly one forest carries a basis, which
# is how a counterfactual arm is written for a Bayesian causal forest
expect_identical(
  atOne,
  predict(fit, xNew, type = "bart", bases = oneBasis[[2L]])
)

# --- the caller's own offset and weights at the predicted rows ---
# an offset enters eta, before any link, length-1 recycled or one per row
expect_identical(
  predict(fit, xNew, type = "bart", bases = zeroBasis, offset = 0.5),
  atZero + 0.5
)
offsetNew <- runif(nNew)
expect_identical(
  predict(fit, xNew, type = "bart", bases = zeroBasis, offset = offsetNew),
  atZero + rep(offsetNew, each = n.samples)
)
# under a gaussian response "ev" is eta itself, and "ppd" is that eta through
# the unchanged sampleFromPPD, so a shared seed makes the equality exact
expect_identical(predict(fit, xNew, type = "ev", bases = zeroBasis), atZero)
set.seed(7)
drawn <- predict(fit, xNew, type = "ppd", bases = zeroBasis)
set.seed(7)
expect_identical(drawn, dbarts:::sampleFromPPD(atZero, fit, NULL, 1L))
# and the band ci.level opts into is taken over the blended draws
band <- predict(fit, xNew, type = "bart", bases = zeroBasis, ci.level = 0.9)
expect_equal(dim(band), c(nNew, 3L))
expect_identical(colnames(band), c("est", "ci.lower", "ci.upper"))
expect_equal(band[, "est"], apply(atZero, 2L, mean))

# --- the bart2 forest() term route: the basis is re-derived from newdata, so
# nothing has to be spelled twice, and an explicit 'bases' overrides that ---
set.seed(812)
nFrame <- 80L
frame <- data.frame(
  x1 = runif(nFrame),
  x2 = runif(nFrame),
  z = rbinom(nFrame, 1L, 0.5)
)
frame$y <- 2 * frame$x1 + frame$z * (1 + frame$x2) + rnorm(nFrame, sd = 0.2)
termFit <- bart2(
  y ~ x1 + x2 + z:forest(x1 + x2),
  data = frame,
  n.trees = 8L,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = 2L,
  n.samples = n.samples,
  keepTrees = TRUE,
  verbose = FALSE
)
# the sugar declares a ONE-column numeric basis, not an indicator pair, so the
# explicit spelling of the same thing is that one column
expect_identical(dim(termFit$bases[[2L]]), c(nFrame, 1L))
newFrame <- frame[1:10, ]
expect_identical(
  predict(termFit, newFrame, type = "bart"),
  predict(termFit, newFrame, type = "bart", bases = matrix(newFrame$z))
)
expect_true(
  max(abs(predict(termFit, frame, type = "bart") - termFit$yhat.train)) < 1e-12
)
# a counterfactual is written into newdata on this route
armed <- newFrame
armed$z <- 1
expect_true(
  max(abs(
    predict(termFit, armed, type = "bart") -
      predict(termFit, newFrame, type = "bart")
  )) >
    1e-6
)
# the variable the basis names is demanded by name, rather than resolved in an
# enclosing scope the way model.frame would resolve it
withoutZ <- newFrame
withoutZ$z <- NULL
expect_error(
  predict(termFit, withoutZ, type = "bart"),
  pattern = "missing variable 'z', required by forest 2"
)

# a factor basis expands to the indicator pair, and the fit's own LEVELS are
# what set the width at the new rows: an all-control arm still gets both
# columns, in the order amplitude j is stated against, rather than the one
# column factor(z) would derive from that newdata alone
factorFit <- bart2(
  y ~ x1 + x2 + forest(x1 + x2, basis = ~ factor(z)),
  data = frame,
  n.trees = 8L,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = 2L,
  n.samples = n.samples,
  keepTrees = TRUE,
  verbose = FALSE
)
expect_identical(dim(factorFit$bases[[2L]]), c(nFrame, 2L))
expect_true(
  max(abs(predict(factorFit, frame, type = "bart") - factorFit$yhat.train)) <
    1e-12
)
allControl <- frame
allControl$z <- 0
expect_identical(
  predict(factorFit, allControl, type = "bart"),
  predict(
    factorFit,
    allControl,
    type = "bart",
    bases = list(NULL, cbind(rep(1, nFrame), rep(0, nFrame)))
  )
)

# --- refusals ---
# 'bases' on a fit with nothing to distribute it over, named by count
plainFit <- fitFrom(NULL)
expect_error(
  predict(plainFit, xNew, bases = list(NULL)),
  pattern = "this fit has 1 forest"
)
# ill-shaped: the wrong list length, the wrong row count, the wrong width
expect_error(
  predict(fit, xNew, type = "bart", bases = list(NULL, NULL, NULL)),
  pattern = "must be a length-2 list"
)
expect_error(
  predict(fit, xNew, type = "bart", bases = list(NULL, cbind(1:3, 3:1))),
  pattern = "same length as 'newdata'"
)
expect_error(
  predict(fit, xNew, type = "bart", bases = list(NULL, rep(1, nNew))),
  pattern = "gives forest 2 1 column; its amplitudes take 2"
)
# a bare value is unambiguous only at one basis-carrying forest
expect_error(
  predict(bothBasesFit, xNew, type = "bart", bases = rep(1, nNew)),
  pattern = "bare value only when exactly one forest carries a basis"
)
# a basis for a forest that declared none, whose amplitude has an implicit
# all-ones column of its own
expect_error(
  predict(
    fit,
    xNew,
    type = "bart",
    bases = list(rep(1, nNew), zeroBasis[[2L]])
  ),
  pattern = "which it declares none of"
)
# and nothing at all, on a fit whose bases arrived as values rather than as a
# term there is any expression to replay
expect_error(
  predict(fit, xNew, type = "bart"),
  pattern = "carries a basis, so the blend needs its 2 columns"
)
# 'bases' does not reach the pre-basis arm
expect_error(
  predict(fit, xNew, type = "forest", bases = zeroBasis),
  pattern = "does not apply to type = \"forest\""
)
# the offset and the weights are checked against the predicted rows
expect_error(
  predict(fit, xNew, type = "bart", bases = zeroBasis, offset = rep(0.5, 3L)),
  pattern = "'offset' must have the same number of rows"
)
expect_error(
  predict(fit, xNew, type = "ppd", bases = zeroBasis, weights = rep(1, 3L)),
  pattern = "'weights' must have the same number of rows"
)
# keepTrees gates BOTH amplitude arms, with one text: without the store the
# replay reports the current trees once for every draw, and the pairing the
# blend is defined by does not exist
noTrees <- fitFrom(twoForests, keepTrees = FALSE)
expect_false(noTrees$fit$control@keepTrees)
expect_error(
  predict(noTrees, xNew, type = "bart", bases = zeroBasis),
  pattern = "requires 'keeptrees'/'keepTrees'"
)
expect_error(
  predict(noTrees, xNew, type = "forest"),
  pattern = "requires 'keeptrees'/'keepTrees'"
)

# --- a single-forest fit's own arms are untouched: no gate, no blend. Without
# the store it still reports the current trees as ONE draw, which is the shape
# the amplitude arms refuse rather than pair against an n.samples glue.
expect_equal(dim(predict(plainFit, xNew, type = "bart")), c(n.samples, nNew))
noTreesPlain <- fitFrom(NULL, keepTrees = FALSE)
plainReplay <- predict(noTreesPlain, xNew, type = "bart")
expect_null(dim(plainReplay))
expect_equal(length(plainReplay), nNew)
