# pdbart and pd2bart on aft and hazard fits: survival, the event probability
# and the cumulative hazard per subject, averaged last, at the default time or
# the times given, with the hazard size check and the survival plot views.

if (!requireNamespace("survival", quietly = TRUE)) {
  exit_file("survival is not installed")
}
Surv <- survival::Surv

fitSmall <- function(...) {
  dbarts::bart(
    ...,
    n.trees = 5L,
    n.samples = 6L,
    n.burn = 4L,
    n.chains = 2L,
    n.threads = 1L,
    seed = 3L,
    keepTrees = TRUE,
    verbose = FALSE
  )
}
# the subject mean of a per-subject quantity, draws x times
subjectMean <- function(draws) apply(draws, c(1L, 2L), mean)

set.seed(21)
n <- 60L
df <- data.frame(a = runif(n), b = rnorm(n))
df$t <- rexp(n, exp(df$a - 0.5))
df$s <- rbinom(n, 1L, 0.8)
atA <- function(value) {
  rows <- df[, c("a", "b")]
  rows$a <- value
  rows
}

aftFit <- fitSmall(Surv(t, s) ~ a + b, df)
hazardFit <- fitSmall(Surv(t, s) ~ a + b, df, family = "hazard")
times <- c(0.3, 0.8, 1.5)

# --- the averages ---
# survival per subject at the grid value, averaged over subjects
pdAft <- dbarts::pdbart(
  aftFit,
  xind = "a",
  levs = list(c(0.2, 0.7)),
  times = times,
  pl = FALSE
)
expect_identical(pdAft$type, "survival")
expect_equal(dim(pdAft$fd[[1L]]), c(12L, 3L, 2L))
expect_identical(dimnames(pdAft$fd[[1L]])[[2L]], format(times))
aftSurvival <- dbarts::survivalProbabilities(aftFit, times, newdata = atA(0.7))
expect_equal(
  pdAft$fd[[1L]][,, 2L],
  subjectMean(aftSurvival),
  check.attributes = FALSE
)
pdHazard <- dbarts::pdbart(
  hazardFit,
  xind = "a",
  levs = list(c(0.2, 0.7)),
  times = times,
  pl = FALSE
)
hazardSurvival <- dbarts::survivalProbabilities(
  hazardFit,
  times,
  newdata = atA(0.7)
)
expect_equal(
  pdHazard$fd[[1L]][,, 2L],
  subjectMean(hazardSurvival),
  check.attributes = FALSE
)
# the event probability and the cumulative hazard, each per subject first
pdEvent <- dbarts::pdbart(
  hazardFit,
  xind = "a",
  levs = list(0.7),
  times = times,
  type = "event",
  pl = FALSE
)
expect_equal(
  pdEvent$fd[[1L]][,, 1L],
  subjectMean(1 - hazardSurvival),
  check.attributes = FALSE
)
pdCumhaz <- dbarts::pdbart(
  aftFit,
  xind = "a",
  levs = list(0.7),
  times = times,
  type = "cumhaz",
  pl = FALSE
)
expect_equal(
  pdCumhaz$fd[[1L]][,, 1L],
  subjectMean(-log(aftSurvival)),
  check.attributes = FALSE
)
expect_true(
  max(abs(pdCumhaz$fd[[1L]][,, 1L] - -log(subjectMean(aftSurvival)))) > 1e-4
)
# across chunk boundaries: a small bound replays a few subjects at a time
setup <- function(fit) {
  sampler <- fit$fit
  list(
    fit = fit,
    sampler = sampler,
    rows = extract(sampler, "predictors")
  )
}
hazardRows <- {
  x <- extract(hazardFit$fit, "predictors")
  x[x[, "period"] == 1, c("a", "b"), drop = FALSE]
}
hazardRows[, "a"] <- 0.7
for (type in c("survival", "cumhaz")) {
  whole <- dbarts:::pdbart.hazardAverage(
    hazardFit,
    hazardFit$fit,
    hazardRows,
    NULL,
    times,
    type,
    NULL,
    5e6
  )
  chunked <- dbarts:::pdbart.hazardAverage(
    hazardFit,
    hazardFit$fit,
    hazardRows,
    NULL,
    times,
    type,
    NULL,
    500
  )
  expect_equal(chunked, whole, info = type)
}
# and with weights, each chunk taking its own subjects
subjectWeights <- runif(nrow(hazardRows))
subjectWeights <- subjectWeights / sum(subjectWeights)
weighted <- function(bound) {
  dbarts:::pdbart.hazardAverage(
    hazardFit,
    hazardFit$fit,
    hazardRows,
    NULL,
    times,
    "survival",
    subjectWeights,
    bound
  )
}
expect_equal(weighted(500), weighted(5e6))
rm(subjectWeights, weighted)
aftRows <- atA(0.7)
noOffset <- function(chunk) NULL
expect_equal(
  dbarts:::pdbart.aftAverage(
    aftFit,
    aftRows,
    noOffset,
    times,
    "event",
    NULL,
    50
  ),
  dbarts:::pdbart.aftAverage(
    aftFit,
    aftRows,
    noOffset,
    times,
    "event",
    NULL,
    5e6
  )
)
rm(setup, hazardRows, type, whole, chunked, aftRows, noOffset)

# shapes: chains lead when split; the times margin is kept at length 1
splitFit <- dbarts::pd2bart(
  Surv(t, s) ~ a + b,
  df,
  levs = list(c(0.2, 0.7), c(0, 1)),
  n.trees = 5L,
  n.samples = 6L,
  n.burn = 4L,
  n.chains = 2L,
  n.threads = 1L,
  seed = 3L,
  combineChains = FALSE,
  verbose = FALSE,
  pl = FALSE
)
expect_equal(dim(splitFit$fd), c(2L, 6L, 1L, 4L))
expect_equal(length(splitFit$times), 1L)
mergedFit <- dbarts::pd2bart(
  aftFit,
  levs = list(c(0.2, 0.7), c(0, 1)),
  pl = FALSE
)
expect_equal(dim(mergedFit$fd), c(12L, 1L, 4L))
# aft's other scales give no times margin
pdLink <- dbarts::pdbart(aftFit, xind = "a", type = "link", pl = FALSE)
expect_identical(pdLink$type, "bart")
expect_null(pdLink$times)
expect_equal(length(dim(pdLink$fd[[1L]])), 2L)
rm(aftSurvival, hazardSurvival, pdEvent, pdCumhaz, splitFit, pdLink)

# --- the default time ---
# the Kaplan-Meier median survival time of the training data
survfitMedian <- function(time, status) {
  unname(quantile(
    survival::survfit(Surv(time, status) ~ 1),
    0.5,
    conf.int = FALSE
  ))
}
expect_equal(
  dbarts::pdbart(aftFit, xind = "a", levs = list(0.5), pl = FALSE)$times,
  survfitMedian(df$t, df$s)
)
expect_equal(
  dbarts::pdbart(hazardFit, xind = "a", levs = list(0.5), pl = FALSE)$times,
  survfitMedian(df$t, df$s)
)
# on a coarse grid, the median of the coarsened times
breaks <- c(0, 0.25, 0.5, 1, 2, max(df$t) + 1)
coarseFit <- fitSmall(
  Surv(t, s) ~ a + b,
  df,
  family = hazard(breaks = breaks)
)
coarsened <- breaks[-1L][findInterval(df$t, breaks, left.open = TRUE)]
expect_equal(
  dbarts::pdbart(coarseFit, xind = "a", levs = list(0.5), pl = FALSE)$times,
  survfitMedian(coarsened, df$s)
)
# under heavy censoring the curve stays above one half, and the default is the
# median of the observed event times
censored <- df
censored$s <- as.integer(seq_len(n) %% 5L == 0L)
expect_true(is.na(survfitMedian(censored$t, censored$s)))
heavyAft <- fitSmall(Surv(t, s) ~ a + b, censored)
heavyHazard <- fitSmall(Surv(t, s) ~ a + b, censored, family = "hazard")
eventMedian <- median(censored$t[censored$s == 1L])
expect_equal(
  dbarts::pdbart(heavyAft, xind = "a", levs = list(0.5), pl = FALSE)$times,
  eventMedian
)
expect_equal(
  dbarts::pdbart(heavyHazard, xind = "a", levs = list(0.5), pl = FALSE)$times,
  eventMedian
)
rm(breaks, coarseFit, coarsened, censored, heavyAft, heavyHazard)
rm(eventMedian)

# --- the hazard size check ---
# subjects x periods up to the largest time x draws x grid values
numPeriods <- sum(hazardFit$periods <= max(times))
count <- n * numPeriods * 12 * 2
expect_error(
  dbarts::pdbart(
    hazardFit,
    xind = "a",
    levs = list(c(0.2, 0.7)),
    times = times,
    n.max.predictions = count - 1,
    pl = FALSE
  ),
  "hazard(breaks = )",
  fixed = TRUE
)
expect_silent(dbarts::pdbart(
  hazardFit,
  xind = "a",
  levs = list(c(0.2, 0.7)),
  times = times,
  n.max.predictions = count,
  pl = FALSE
))
# counted in double precision, past 2^31
bigGrid <- list(seq(0, 1, length.out = 80000L))
bigCount <- n * length(hazardFit$periods) * 12 * 80000
expect_true(bigCount > 2^31)
expect_error(
  dbarts::pdbart(
    hazardFit,
    xind = "a",
    levs = bigGrid,
    times = max(hazardFit$periods),
    n.max.predictions = bigCount - 1,
    pl = FALSE
  ),
  format(bigCount, digits = 3L),
  fixed = TRUE
)
# a 20-period fit at bart's default run does not trip the default limit
twentyFit <- dbarts::bart(
  Surv(t, s) ~ a + b,
  df,
  family = hazard(breaks = seq(0, max(df$t) + 1, length.out = 21L)),
  n.trees = 5L,
  n.threads = 1L,
  seed = 3L,
  keepTrees = TRUE,
  verbose = FALSE
)
expect_equal(length(twentyFit$periods), 20L)
expect_equal(
  dim(dbarts::pdbart(twentyFit, xind = "a", pl = FALSE)$fd[[1L]]),
  c(2000L, 1L, 11L)
)
expect_error(
  dbarts::pdbart(hazardFit, n.max.predictions = 0, pl = FALSE),
  "'n.max.predictions' must be a positive number"
)
rm(numPeriods, count, bigGrid, bigCount, twentyFit)

# --- refusals ---
# the period is the time axis, not a predictor
expect_identical(pdHazard$xlbs, "a")
expect_identical(
  dbarts::pdbart(hazardFit, levs = list(0.5, 0), pl = FALSE)$xlbs,
  c("a", "b")
)
expect_error(
  dbarts::pdbart(hazardFit, xind = "period", pl = FALSE),
  "time axis"
)
# times only on a survival scale
gaussianFit <- fitSmall(t ~ a + b, df)
expect_error(
  dbarts::pdbart(gaussianFit, times = 1, pl = FALSE),
  "'times' applies"
)
expect_error(
  dbarts::pdbart(aftFit, type = "bart", times = 1, pl = FALSE),
  "'times' applies"
)
expect_error(
  dbarts::pdbart(hazardFit, type = "ev", pl = FALSE),
  "does not take type = \"ev\""
)
expect_error(
  dbarts::pdbart(aftFit, times = -1, pl = FALSE),
  "'times' must be finite and positive"
)
# a hazard sampler's rows are person-period rows; an aft sampler takes the
# log-time scale when named
expect_error(dbarts::pdbart(hazardFit$fit, pl = FALSE), "hazard sampler")
expect_error(dbarts::pdbart(aftFit$fit, pl = FALSE), "means survival")
pdSampler <- dbarts::pdbart(
  aftFit$fit,
  xind = 1L,
  levs = list(0.5),
  type = "bart",
  pl = FALSE
)
expect_null(pdSampler$times)
rm(gaussianFit, pdSampler)

# --- plots ---
pdf(NULL)
expect_silent(plot(pdHazard))
expect_silent(plot(pdHazard, plot.type = "curves"))
expect_silent(plot(mergedFit))
expect_silent(plot(
  dbarts::pd2bart(
    hazardFit,
    levs = list(c(0.2, 0.7), c(0, 1)),
    times = times,
    pl = FALSE
  ),
  plot.type = "curves"
))
expect_error(plot(mergedFit, plot.type = "curves"), "at least three times")
expect_error(
  plot(
    dbarts::pdbart(aftFit, xind = "a", type = "link", pl = FALSE),
    plot.type = "curves"
  ),
  "times margin"
)
dev.off()
# the labels name the scale and the time read
expect_identical(dbarts:::pdScaleLabel(pdHazard), "survival probability")
expect_identical(
  dbarts:::pdScaleLabel(mergedFit),
  paste0("survival probability at t = ", format(mergedFit$times))
)
expect_identical(
  dbarts:::pdScaleLabel(list(type = "bart", family = "aft")),
  "log time"
)
expect_identical(
  dbarts:::pdScaleLabel(list(type = "cumhaz", times = 1:3)),
  "cumulative hazard"
)

rm(Surv, fitSmall, subjectMean, n, df, atA, aftFit, hazardFit, times)
rm(pdAft, pdHazard, mergedFit, survfitMedian)
