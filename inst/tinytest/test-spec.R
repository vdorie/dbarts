source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

x <- testData$x
y <- testData$y

control <- dbarts::dbartsControl(
  n.chains = 1L,
  n.samples = 15L,
  n.burn = 5L,
  n.threads = 1L,
  seed = 11L,
  updateState = FALSE
)

## the spec surface resolves what dbarts() resolves: a sampler built from the
## returned triple draws the same samples as the ordinary entry point
reference <- dbarts::dbarts(x, y, control = control, n.trees = 10L)
referenceSamples <- reference$run()

data <- dbarts::dbartsData(x, y)
spec <- dbarts::dbartsSpec(data, control = control, n.trees = 10L)

expect_equal(names(spec), c("control", "model", "data", "family"))
expect_true(inherits(spec$control, "dbartsControl"))
expect_true(inherits(spec$model, "dbartsModel"))
expect_true(inherits(spec$data, "dbartsData"))
expect_equal(spec$family, "gaussian")
# no mask: a latent family's 0/1 case weights ride the data object
expect_null(spec$active)

sampler <- new("dbartsSampler", spec$control, spec$model, spec$data)
samples <- sampler$run()

## bitwise, not merely close: the two paths share one resolution
expect_identical(samples$train, referenceSamples$train)
expect_identical(samples$sigma, referenceSamples$sigma)

## the resolved family is the one the sampler is built with, and it is
## never left as "auto"
expect_equal(spec$model@family, "gaussian")

## a binary response resolves to probit through the same path bart() uses
binaryData <- dbarts::dbartsData(x, as.double(y > median(y)))
binarySpec <- suppressMessages(dbarts::dbartsSpec(
  binaryData,
  control = control,
  n.trees = 10L
))
expect_equal(binarySpec$family, "probit")
expect_true(binarySpec$control@binary)
## a fixed-unit-scale family takes probit's leaf scale and fixed residual prior
expect_equal(binarySpec$model@leaf.scale, 3.0)

## features that were previously unreachable without hand-built attributes
monotoneSpec <- dbarts::dbartsSpec(
  data,
  control = control,
  monotone = c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  n.trees = 10L
)
expect_equal(
  attr(monotoneSpec$model, "monotone"),
  c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
)
## the constraint forces birth/death-only proposals, as it does in dbarts()
expect_equal(monotoneSpec$control@proposal.probs[["birth_death"]], 1.0)
expect_equal(monotoneSpec$control@proposal.probs[["swap"]], 0.0)
expect_equal(monotoneSpec$control@proposal.probs[["perturb"]], 0.0)
expect_equal(monotoneSpec$control@proposal.probs[["rule_gibbs"]], 0.0)

additiveSpec <- dbarts::dbartsSpec(
  data,
  control = control,
  interactions = interactions(max.order = 1L),
  n.trees = 10L
)
expect_equal(attr(additiveSpec$model, "interaction.max.order"), 1L)

varianceSpec <- dbarts::dbartsSpec(
  data,
  control = control,
  variance = TRUE,
  n.trees = 10L
)
expect_equal(varianceSpec$model@variance@n.trees, 40L)

studentSpec <- dbarts::dbartsSpec(
  data,
  control = control,
  family = student(df = 5),
  n.trees = 10L
)
expect_equal(attr(studentSpec$model, "resid.df"), 5.0)

## the prior vocabulary is NSE and must resolve in dbarts's namespace, not the
## caller's: a bare k = chi(...) hyperprior works exactly as it does in dbarts()
hyperSpec <- dbarts::dbartsSpec(
  data,
  control = control,
  leaf.prior = normal(k = chi(1.25, Inf)),
  n.trees = 10L
)
expect_true(inherits(hyperSpec$model@leaf.hyperprior, "dbartsChiHyperprior"))

## aft needs its status vector, and refuses to be built without one
logTime <- log(abs(y) + 1)
status <- rep_len(c(1, 0), length(logTime))
aftData <- dbarts::dbartsData(x, logTime)
aftSpec <- dbarts::dbartsSpec(
  aftData,
  control = control,
  family = "aft",
  survival = status,
  n.trees = 10L
)
expect_equal(aftSpec$family, "aft")
expect_equal(aftSpec$data@status, status)
expect_error(
  dbarts::dbartsSpec(aftData, control = control, family = "aft", n.trees = 10L),
  "needs a 'survival' status vector"
)
expect_error(
  dbarts::dbartsSpec(
    aftData,
    control = control,
    survival = status,
    n.trees = 10L
  ),
  "only used by family"
)
expect_error(
  dbarts::dbartsSpec(
    aftData,
    control = control,
    family = "aft",
    survival = rep_len(2, length(logTime)),
    n.trees = 10L
  ),
  "must be 0 \\(censored\\) or 1 \\(event\\)"
)

## argument validation
expect_error(dbarts::dbartsSpec(x), "must be a dbartsData")
expect_error(
  dbarts::dbartsSpec(data, control = list()),
  "must be a dbartsControl"
)

## an explicit sigest states sigma; the default estimates it from the data's
## response, whatever a data object carries
sigmaData <- dbarts::dbartsData(x, y)
sigmaData@sigma <- 2.5
expect_equal(
  dbarts::dbartsSpec(sigmaData, control = control, n.trees = 10L)$data@sigma,
  dbarts::dbartsSpec(
    dbarts::dbartsData(x, y),
    control = control,
    n.trees = 10L
  )$data@sigma
)
expect_equal(
  dbarts::dbartsSpec(
    sigmaData,
    control = control,
    sigest = 1.5,
    n.trees = 10L
  )$data@sigma,
  1.5
)

## an unset sigma is estimated during resolution, as it is for dbarts()
expect_false(is.na(spec$data@sigma))

## the seed argument mirrors dbartsControl(seed = )
expect_equal(
  dbarts::dbartsSpec(
    data,
    control = control,
    seed = 42L,
    n.trees = 10L
  )$control@seed,
  42L
)
expect_equal(
  dbarts::dbartsSpec(data, control = control, n.trees = 10L)$control@seed,
  11L
)
