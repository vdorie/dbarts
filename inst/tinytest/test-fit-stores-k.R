# A fit stores the k its sampler recorded under every naming of the leaf prior,
# the leaf prior it ran under (leaf.prior) and the scalars it held fixed
# (fixed); extract answers sigma, shape, k and leaf.prior.sd with draws for a
# parameter the fit sampled and one number for one it held fixed, each checked
# here against the sampler's own readers. The draw channels keep their layout.

set.seed(23L)
n <- 40L
x <- matrix(runif(n * 3L), n, 3L, dimnames = list(NULL, c("x1", "x2", "x3")))
f <- 2 * sin(pi * x[, 1L]) + x[, 2L]
yGaussian <- f + rnorm(n, sd = 0.3)
yBinary <- rbinom(n, 1L, pnorm(f - 1.5))
yCount <- rnbinom(n, size = 3, mu = exp(0.5 * f))
yOrdinal <- cut(
  f + rnorm(n, sd = 0.5),
  3L,
  labels = c("a", "b", "c"),
  ordered_result = TRUE
)
yClass <- factor(sample(c("a", "b", "c"), n, TRUE))
yHurdle <- ifelse(yBinary == 1L, exp(0.5 * f + rnorm(n, sd = 0.3)), 0)
timeToEvent <- exp(0.5 * f + rnorm(n, sd = 0.3))
status <- rbinom(n, 1L, 0.8)

nSamples <- 5L
# the front door re-reads its own call, so a literal TRUE or FALSE is what
# reaches keepSampler; the flag's name is one no argument of bart abbreviates
fitOf <- function(y, .keep, ...) {
  suppressMessages(
    if (.keep) {
      bart(
        x,
        y,
        ...,
        n.trees = 5L,
        n.samples = nSamples,
        n.burn = 3L,
        n.chains = 2L,
        n.threads = 1L,
        seed = 7L,
        keepSampler = TRUE,
        verbose = FALSE
      )
    } else {
      bart(
        x,
        y,
        ...,
        n.trees = 5L,
        n.samples = nSamples,
        n.burn = 3L,
        n.chains = 2L,
        n.threads = 1L,
        seed = 7L,
        keepSampler = FALSE,
        verbose = FALSE
      )
    }
  )
}

# the sampler's reading of a quantity, one value when its chains agree
held <- function(values) {
  values <- unique(as.vector(values))
  expect_identical(length(values), 1L)
  values
}

# One extract answer against the reader it must agree with: a quantity the fit
# held fixed is one number with no chain margin under either combineChains,
# and one it sampled is its draws, whose last is what the sampler holds after
# the run.
expectParameter <- function(fit, type, reader, fixed, info) {
  combined <- extract(fit, type)
  split <- extract(fit, type, combineChains = FALSE)
  if (fixed) {
    expect_identical(combined, reader, info = info)
    expect_identical(split, reader, info = info)
    expect_null(dim(combined), info = info)
  } else {
    expect_identical(length(combined), 2L * nSamples, info = info)
    expect_identical(dim(split), c(2L, nSamples), info = info)
    expect_equal(split[, nSamples], reader, info = info)
    expect_identical(as.vector(t(split)), combined, info = info)
  }
}

# --- every class, drawn and fixed, chains combined and split, the sampler kept
# and not ---

scenarios <- list(
  gaussian = list(
    function(keep) fitOf(yGaussian, keep),
    sigma = "drawn",
    k = "fixed"
  ),
  "gaussian, sigma fixed" = list(
    function(keep) {
      fitOf(
        yGaussian,
        keep,
        family = dbartsFamilies$gaussian(sigma = dbartsPriors$fixed(0.3))
      )
    },
    sigma = "fixed",
    k = "fixed"
  ),
  "gaussian, k drawn" = list(
    function(keep) fitOf(yGaussian, keep, k = dbartsPriors$chi(2, 1)),
    sigma = "drawn",
    k = "drawn"
  ),
  "gaussian, sd named" = list(
    function(keep) {
      fitOf(yGaussian, keep, leaf.prior = dbartsPriors$normal(sd = 0.7))
    },
    sigma = "drawn",
    k = "fixed"
  ),
  "gaussian, sd drawn" = list(
    function(keep) {
      fitOf(
        yGaussian,
        keep,
        leaf.prior = dbartsPriors$normal(sd = dbartsPriors$invchi(1.5, 0.7))
      )
    },
    sigma = "drawn",
    k = "drawn"
  ),
  student = list(
    function(keep) {
      fitOf(yGaussian, keep, family = dbartsFamilies$student(df = 5))
    },
    sigma = "drawn",
    k = "fixed"
  ),
  aft = list(
    function(keep) fitOf(survival::Surv(timeToEvent, status), keep),
    sigma = "drawn",
    k = "fixed"
  ),
  probit = list(
    function(keep) fitOf(yBinary, keep),
    sigma = "none",
    k = "drawn"
  ),
  "probit, k fixed" = list(
    function(keep) fitOf(yBinary, keep, k = 2),
    sigma = "none",
    k = "fixed"
  ),
  logistic = list(
    function(keep) fitOf(yBinary, keep, family = "logistic"),
    sigma = "none",
    k = "drawn"
  ),
  "nbinom, shape drawn" = list(
    function(keep) fitOf(yCount, keep, family = "nbinom"),
    sigma = "none",
    k = "drawn",
    shape = "drawn"
  ),
  "nbinom, shape fixed" = list(
    function(keep) {
      fitOf(yCount, keep, family = dbartsFamilies$nbinom(shape = 3), k = 2)
    },
    sigma = "none",
    k = "fixed",
    shape = "fixed"
  ),
  ordinal = list(
    function(keep) fitOf(yOrdinal, keep, family = "ordinal"),
    sigma = "none",
    k = "fixed"
  ),
  "ordinal, k drawn" = list(
    function(keep) {
      fitOf(yOrdinal, keep, family = "ordinal", k = dbartsPriors$chi(2, 1))
    },
    sigma = "none",
    k = "drawn"
  ),
  multinomial = list(
    function(keep) fitOf(yClass, keep, family = "multinomial"),
    sigma = "none",
    k = "fixed"
  )
)

for (name in names(scenarios)) {
  scenario <- scenarios[[name]]
  fit <- scenario[[1L]](TRUE)
  unkept <- scenario[[1L]](FALSE)
  sampler <- fit$fit
  expect_inherits(sampler, "dbartsSampler", info = name)
  expect_null(unkept$fit, info = name)

  # the descriptors are on every fit and the draws agree with the unkept twin
  expect_true(is.list(fit$leaf.prior), info = name)
  expect_true(is.list(fit$fixed), info = name)
  expect_identical(fit$fixed, unkept$fixed, info = name)
  expect_identical(fit$leaf.prior, unkept$leaf.prior, info = name)

  # k: the sampler's own, never present as a channel when fixed
  fixedK <- scenario$k == "fixed"
  expect_identical(is.null(fit[["k"]]), fixedK, info = name)
  expect_identical(!is.null(fit$fixed[["k"]]), fixedK, info = name)
  kReader <- sampler$getK()
  kReader <- if (fixedK) held(kReader) else as.vector(kReader)
  expectParameter(fit, "k", kReader, fixedK, paste(name, "k"))
  expect_identical(extract(unkept, "k"), extract(fit, "k"), info = name)
  if (fixedK) {
    expect_identical(fit$fixed$k, kReader, info = name)
  } else {
    expect_identical(
      fit[["k"]][c(nSamples, 2L * nSamples)],
      as.vector(sampler$getK()),
      info = name
    )
  }

  # leaf.prior.sd: the reader's anchor over k, draws or one number alike
  anchor <- sampler$getLeafPrior()$k.scale
  if (is.null(anchor)) {
    anchor <- sampler$getLeafPrior()[[1L]]$k.scale
  }
  expect_equal(fit$leaf.prior$k.scale, anchor, info = name)
  for (combine in c(TRUE, FALSE)) {
    expect_equal(
      extract(fit, "leaf.prior.sd", combineChains = combine),
      anchor / extract(fit, "k", combineChains = combine),
      info = paste(name, "leaf.prior.sd", combine)
    )
    expect_equal(
      extract(unkept, "leaf.prior.sd", combineChains = combine),
      extract(fit, "leaf.prior.sd", combineChains = combine),
      info = paste(name, "unkept leaf.prior.sd", combine)
    )
  }

  # sigma: the sampler's scale, draws, or 1 on a family with none
  if (scenario$sigma == "none") {
    expect_identical(extract(fit, "sigma"), 1, info = name)
    expect_identical(
      extract(fit, "sigma", combineChains = FALSE),
      1,
      info = name
    )
  } else {
    fixedSigma <- scenario$sigma == "fixed"
    expect_identical(!is.null(fit$fixed[["sigma"]]), fixedSigma, info = name)
    reader <- sampler$getSigmas()
    expectParameter(
      fit,
      "sigma",
      if (fixedSigma) held(reader) else reader,
      fixedSigma,
      paste(name, "sigma")
    )
    # the channel keeps its layout, a fixed value repeated per draw
    expect_identical(length(fit$sigma), 2L * nSamples, info = name)
    expect_identical(
      extract(unkept, "sigma"),
      extract(fit, "sigma"),
      info = name
    )
  }

  # shape, the count families' own
  if (!is.null(scenario$shape)) {
    fixedShape <- scenario$shape == "fixed"
    expect_identical(!is.null(fit$fixed[["shape"]]), fixedShape, info = name)
    reader <- sampler$getShape()
    expectParameter(
      fit,
      "shape",
      if (fixedShape) held(reader) else reader,
      fixedShape,
      paste(name, "shape")
    )
    expect_identical(length(fit$shape), 2L * nSamples, info = name)
  }
}

# a fixed sigma is the sampler's, the square root of the variance fixed() names
fixedSigmaFit <- fitOf(
  yGaussian,
  TRUE,
  family = dbartsFamilies$gaussian(sigma = dbartsPriors$fixed(0.3))
)
expect_equal(extract(fixedSigmaFit, "sigma"), sqrt(0.3))
expect_equal(fixedSigmaFit$fixed$sigma, sqrt(0.3))
expect_true(all(fixedSigmaFit$sigma == sqrt(0.3)))
expect_identical(
  names(fixedSigmaFit$fixed),
  c("sigma", "k")
)

# a fixed Student df is recorded, and the channel repeats it
studentFit <- fitOf(yGaussian, TRUE, family = dbartsFamilies$student(df = 5))
expect_equal(studentFit$fixed$resid.df, 5)
expect_true(all(studentFit$resid.df == 5))
expect_null(
  fitOf(yGaussian, TRUE, family = dbartsFamilies$student())$fixed$resid.df
)

# no new component's name begins with a name an existing read takes by prefix:
# a fixed k leaves fit$k NULL, and the channels that are there are unchanged
for (name in names(scenarios)) {
  fit <- scenarios[[name]][[1L]](FALSE)
  expect_identical(is.null(fit$k), is.null(fit[["k"]]), info = name)
  expect_true(all(c("leaf.prior", "fixed") %in% names(fit)), info = name)
  # summary's printed output names what the fit held fixed, and only then; an
  # ordinal fit always holds its first threshold
  printed <- capture.output(print(summary(fit)))
  expect_identical(
    any(grepl("(Fixed, not sampled: ", printed, fixed = TRUE)),
    length(fit$fixed) > 0L || inherits(fit, "bartOrdinal"),
    info = name
  )
}

# a heteroscedastic fit has no scalar sigma whatever its residual prior says:
# extract returns the per-observation surface, and fixed names no sigma
heteroscedastic <- fitOf(
  yGaussian,
  TRUE,
  variance = varianceForest(),
  family = dbartsFamilies$gaussian(sigma = dbartsPriors$fixed(0.09))
)
expect_null(heteroscedastic$fixed[["sigma"]])
expect_null(heteroscedastic[["sigma"]])
expect_identical(dim(extract(heteroscedastic, "sigma")), c(2L * nSamples, n))

# --- several forests: one named number per forest, or the one 'forest' picks ---

dat <- data.frame(x, y = yGaussian + rbinom(n, 1L, 0.5))
dat$z <- rbinom(n, 1L, 0.5)
forests <- bart(
  y ~ x1 + x2 + x3 + forest(x1 + x2, basis = ~z),
  dat,
  n.trees = 5L,
  n.samples = nSamples,
  n.burn = 3L,
  n.chains = 2L,
  n.threads = 1L,
  seed = 7L,
  keepSampler = TRUE,
  verbose = FALSE
)
expect_identical(forests$n.forests, 2L)
expect_identical(names(forests$leaf.prior), c("forest1", "forest2"))
expect_null(forests[["k"]])
perForest <- forests$fit$getLeafPrior()
anchors <- vapply(perForest, function(prior) prior$k.scale, 0)
expect_equal(
  extract(forests, "leaf.prior.sd"),
  setNames(
    anchors / apply(forests$fit$getK(), 1L, held),
    c("forest1", "forest2")
  )
)
expect_equal(
  extract(forests, "k"),
  c(forest1 = 1, forest2 = 1)
)
expect_equal(
  extract(forests, "k", combineChains = FALSE),
  extract(forests, "k")
)
expect_equal(
  extract(forests, "leaf.prior.sd", forest = 2L),
  unname(anchors[2L])
)
expect_equal(
  extract(forests, "leaf.prior.sd", forest = "forest1"),
  unname(anchors[1L])
)
expect_equal(forests$fixed$k, c(forest1 = 1, forest2 = 1))
expect_error(extract(forests, "sigma", forest = 1L), "model parameter")
expect_error(
  extract(fitOf(yGaussian, TRUE), "k", forest = 1L),
  "not a per-forest quantity"
)
expect_error(
  extract(forests, "k", forest = 3L),
  "'forest' index must be between 1 and 2",
  fixed = TRUE
)

# --- a hurdle fit: a list of both parts, each its draws or one number ---

hurdle <- fitOf(yHurdle, TRUE, family = "hurdle.lognormal")
for (combine in c(TRUE, FALSE)) {
  kParts <- extract(hurdle, "k", combineChains = combine)
  expect_identical(names(kParts), c("zero", "positive"))
  expect_identical(
    kParts$zero,
    extract(hurdle$zero, "k", combineChains = combine)
  )
  expect_identical(kParts$positive, held(hurdle$positive$fit$getK()))
  sdParts <- extract(hurdle, "leaf.prior.sd", combineChains = combine)
  expect_equal(sdParts$zero, hurdle$zero$leaf.prior$k.scale / kParts$zero)
  expect_equal(
    sdParts$positive,
    hurdle$positive$fit$getLeafPrior()$k.scale / kParts$positive
  )
}
expect_identical(
  extract(hurdle, "sigma"),
  extract(hurdle$positive, "sigma")
)

# --- a fit saved before fits stored these descriptors answers from its
# channels where they suffice and is refused by name where they do not ---

old <- fitOf(yGaussian, FALSE, k = dbartsPriors$chi(2, 1))
old$leaf.prior <- NULL
old$fixed <- NULL
expect_identical(extract(old, "k"), old$k)
expect_identical(extract(old, "sigma"), old$sigma)
expect_error(extract(old, "leaf.prior.sd"), "saved before fits recorded")
oldFixed <- fitOf(yGaussian, FALSE)
oldFixed$leaf.prior <- NULL
oldFixed$fixed <- NULL
expect_error(extract(oldFixed, "k"), "saved before fits recorded")
oldBinary <- fitOf(yBinary, FALSE)
oldBinary$fixed <- NULL
expect_identical(extract(oldBinary, "sigma"), 1)
# its constant channels tabulate as they always did
oldSigma <- fitOf(
  yGaussian,
  FALSE,
  family = dbartsFamilies$gaussian(sigma = dbartsPriors$fixed(0.3))
)
oldSigma$fixed <- NULL
expect_identical(summary(oldSigma, vars = "sigma")$stats$variable, "sigma")

# --- summary tabulates what was sampled and names what was fixed ---

summaryFixed <- summary(fixedSigmaFit)
expect_null(summaryFixed$stats)
expect_equal(summaryFixed$fixed, list(sigma = sqrt(0.3), k = 2))
expect_stdout(
  print(summaryFixed),
  "(Fixed, not sampled: sigma = 0.5477, k = 2)",
  fixed = TRUE
)
sampledOnly <- summary(fitOf(yGaussian, TRUE))
expect_identical(sampledOnly$stats$variable, "sigma")
expect_equal(sampledOnly$fixed, list(k = 2))
# the leaf scale follows the naming: k on a k-named fit, leaf.prior.sd on an
# sd-named one, and either on request
sdNamed <- summary(fitOf(
  yGaussian,
  TRUE,
  leaf.prior = dbartsPriors$normal(sd = 0.7)
))
expect_equal(sdNamed$fixed, list(leaf.prior.sd = 0.7))
expect_identical(sdNamed$stats$variable, "sigma")
asK <- summary(
  fitOf(yGaussian, TRUE, leaf.prior = dbartsPriors$normal(sd = 0.7)),
  vars = "k"
)
expect_equal(asK$fixed, list(k = 2))
# the ordinal's first threshold is pinned and named; the sampled ones tabulate
ordinalSummary <- summary(fitOf(yOrdinal, TRUE, family = "ordinal"))
expect_identical(ordinalSummary$stats$variable, "threshold[2]")
expect_identical(ordinalSummary$fixed[["threshold[1]"]], 0)
# a shape the fit held fixed is named, a drawn one tabulated
nbinomFixed <- summary(
  fitOf(yCount, TRUE, family = dbartsFamilies$nbinom(shape = 3), k = 2)
)
expect_null(nbinomFixed$stats)
expect_equal(nbinomFixed$fixed, list(shape = 3, k = 2))
nbinomDrawn <- summary(fitOf(yCount, TRUE, family = "nbinom"))
expect_true("shape" %in% nbinomDrawn$stats$variable)
# one value per forest on a fit with several
expect_stdout(
  print(summary(forests)),
  "leaf.prior.sd[forest1] = ",
  fixed = TRUE
)

# --- plot and print: a parameter held fixed has no trace ---

# the panels a plot actually draws, counted at each new plot
panelsDrawn <- function(fit) {
  hooks <- getHook("plot.new")
  panels <- 0L
  pdf(NULL)
  on.exit({
    setHook("plot.new", hooks, "replace")
    dev.off()
  })
  setHook("plot.new", function(...) panels <<- panels + 1L, "append")
  plot(fit)
  panels
}
expect_identical(panelsDrawn(fitOf(yGaussian, TRUE)), 2L)
expect_identical(panelsDrawn(fixedSigmaFit), 1L)
nbinomFixedFit <- fitOf(
  yCount,
  TRUE,
  family = dbartsFamilies$nbinom(shape = 3),
  k = 2
)
nbinomDrawnFit <- fitOf(yCount, TRUE, family = "nbinom")
expect_identical(panelsDrawn(nbinomDrawnFit), 2L)
expect_identical(panelsDrawn(nbinomFixedFit), 1L)
expect_identical(panelsDrawn(hurdle), 4L)
hurdleFixed <- fitOf(
  yHurdle,
  TRUE,
  family = dbartsFamilies$hurdle.lognormal(sigma = dbartsPriors$fixed(0.3))
)
expect_identical(panelsDrawn(hurdleFixed), 3L)
expect_stdout(print(nbinomFixedFit), "shape (r): fixed at 3", fixed = TRUE)
expect_stdout(print(nbinomDrawnFit), "posterior mean shape (r)", fixed = TRUE)

# --- a warm start from a donor that drew sigma and k seeds a recipient that
# holds both fixed with its trees alone: the fixed values stay the recipient's,
# one value each ---

donor <- fitOf(yGaussian, TRUE, k = dbartsPriors$chi(2, 1))
warmed <- fitOf(
  yGaussian,
  TRUE,
  family = dbartsFamilies$gaussian(sigma = dbartsPriors$fixed(1)),
  warm.start = donor
)
expect_identical(length(unique(warmed$fit$getSigmas())), 1L)
expect_equal(warmed$fixed$sigma, 1)
expect_identical(warmed$fixed$k, 2)
for (type in c("sigma", "k")) {
  expect_identical(extract(warmed, type), warmed$fixed[[type]], info = type)
  expect_identical(
    extract(warmed, type, combineChains = FALSE),
    warmed$fixed[[type]],
    info = type
  )
}
expect_equal(extract(warmed, "leaf.prior.sd"), warmed$leaf.prior$k.scale / 2)
expect_stdout(
  print(summary(warmed)),
  "(Fixed, not sampled: sigma = 1, k = 2)",
  fixed = TRUE
)

# --- summary names every parameter held fixed, the Student-t df and a
# multinomial fit's k among them ---

studentSummary <- capture.output(
  print(summary(fitOf(
    yGaussian,
    TRUE,
    family = dbartsFamilies$student(df = 5)
  )))
)
expect_true(any(grepl(
  "(Fixed, not sampled: k = 2, resid.df = 5)",
  studentSummary,
  fixed = TRUE
)))
multinomialSummary <- summary(fitOf(yClass, TRUE, family = "multinomial"))
expect_equal(multinomialSummary$fixed, list(k = 2))
expect_stdout(
  print(multinomialSummary),
  "(Fixed, not sampled: k = 2)",
  fixed = TRUE
)
# a sigma of 1 the family pins is not named
probitSummary <- summary(fitOf(yBinary, TRUE))
expect_false("sigma" %in% names(probitSummary$fixed))
expect_false(any(grepl(
  "sigma",
  capture.output(print(probitSummary)),
  fixed = TRUE
)))
# a sampled df is a row of the default table and is not on the line
sampledDf <- summary(fitOf(
  yGaussian,
  TRUE,
  family = dbartsFamilies$student()
))
expect_true("resid.df" %in% sampledDf$stats$variable)
expect_false("resid.df" %in% names(sampledDf$fixed))
