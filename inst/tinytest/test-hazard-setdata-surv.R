# $setData on a hazard sampler takes survival data and spreads it onto the
# sampler's own periods as dbarts() spreads it at creation.

make <- function(n) {
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  d$time <- ceiling(rexp(n) * 3)
  d$status <- rbinom(n, 1L, 0.7)
  d
}
run <- function(...) {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE,
    n.burn = 0L,
    n.samples = 5L,
    ...
  )
}
form <- survival::Surv(time, status) ~ x1 + x2

set.seed(7)
first <- make(40L)
second <- make(30L)
sampler <- dbarts(
  form,
  first,
  family = "hazard",
  n.trees = 5L,
  control = run(),
  seed = 5L
)
periods <- sampler$model@hazard.periods
K <- length(periods)
# the second data set's times lie on the sampler's grid
second$time <- pmax(pmin(second$time, max(periods)), 1)

# the draws are those of a sampler created on the same subjects with the grid
# the first one holds, in any number of subjects
swapped <- sampler$copy()
swapped$setData(dbartsData(form, second))
fresh <- dbarts(
  form,
  second,
  family = hazard(breaks = c(0, periods)),
  n.trees = 5L,
  control = run(),
  seed = 5L
)
expect_identical(swapped$data@x, fresh$data@x)
expect_identical(swapped$data@y, fresh$data@y)
expect_identical(swapped$data@varTypes, fresh$data@varTypes)
expect_identical(swapped$run(0L, 4L)$train, fresh$run(0L, 4L)$train)

# the rows are one per period at risk, the period last
expect_equal(colnames(swapped$data@x), c("x1", "x2", "period"))
expect_equal(
  nrow(swapped$data@x),
  sum(findInterval(second$time, periods, left.open = TRUE) + 1L)
)

# a time past the last period leaves its subject at risk through every period
# with no event, as a subject out of the likelihood is at creation
late <- second
late$time[1:5] <- max(periods) + 3
late$status[1:5] <- 1L
lateSampler <- sampler$copy()
lateSampler$setData(dbartsData(form, late))
lateCounts <- findInterval(late$time, periods, left.open = TRUE) + 1L
lateCounts <- pmin(lateCounts, K)
expect_equal(nrow(lateSampler$data@x), sum(lateCounts))
expect_equal(sum(lateSampler$data@y[seq_len(5L * K)]), 0)
expect_equal(
  as.vector(lateSampler$data@x[seq_len(5L * K), "period"]),
  rep(seq_len(K), 5L)
)

# a 0/1 weight repeats over a subject's periods
weighted <- sampler$copy()
w <- rep(c(1, 0, 1), length.out = nrow(second))
weighted$setData(dbartsData(form, second, weights = w))
counts <- findInterval(second$time, periods, left.open = TRUE) + 1L
expect_equal(weighted$data@weights, rep(w, counts))

# a logistic hazard sampler takes it alike
logistic <- dbarts(
  form,
  first,
  family = hazard(link = "logistic"),
  n.trees = 5L,
  control = run(),
  seed = 5L
)
logistic$setData(dbartsData(form, second))
expect_equal(colnames(logistic$data@x), c("x1", "x2", "period"))

# a predictor called period is refused, as at creation
renamed <- second
renamed$period <- renamed$x1
expect_error(
  sampler$copy()$setData(dbartsData(
    survival::Surv(time, status) ~ period + x2,
    renamed
  )),
  "rename the predictor 'period'"
)

# a sampler that is not a hazard one reads the data as it did
plain <- dbarts(
  x1 ~ x2,
  first,
  n.trees = 5L,
  control = run(),
  seed = 5L
)
expect_error(
  plain$copy()$setData(dbartsData(form, second)),
  "requires the same predictors"
)

# a hand-built person-period matrix whose last column is not named period is
# refused at $setData and at new(), before any count is read, the message
# naming the doors that fit
pp <- as.matrix(sampler$data@x)
handBuilt <- pp
colnames(handBuilt)[ncol(handBuilt)] <- "t"
given <- dbartsData(handBuilt, sampler$data@y)
expect_error(
  sampler$copy()$setData(given),
  "must be named 'period'.*Surv\\(time, status\\).*family = \"binomial\""
)
expect_error(
  new("dbartsSampler", run(), sampler$model, given),
  "must be named 'period'.*Surv\\(time, status\\).*family = \"binomial\""
)
# with a named count too, which would otherwise stop on the missing name
given <- dbartsData(
  handBuilt,
  sampler$data@y,
  n.cuts = c(x1 = 4L, x2 = 4L, t = 3L)
)
expect_error(sampler$copy()$setData(given), "must be named 'period'")
expect_error(
  new("dbartsSampler", run(), sampler$model, given),
  "must be named 'period'"
)
# the same matrix under the name period, and the sampler's own data, stand
renamedBack <- dbartsData(pp, sampler$data@y)
expect_silent(sampler$copy()$setData(renamedBack))
expect_silent(sampler$copy()$setData(sampler$data))
# and the rows fit as a binary fit, the period a predictor under any name
expect_silent(
  dbarts(
    handBuilt,
    sampler$data@y,
    family = "binomial",
    n.trees = 5L,
    control = run()
  )
)
