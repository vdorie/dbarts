# summary()'s forest block reports, for each forest with a basis, the average
# effect per unit of the basis over the rows in the likelihood and its spread
# across those rows, each a posterior median with a 5% to 95% interval, each
# taken per draw with the median last. The per-draw effect is the forest's own
# fit times its amplitude, which is what predict() gives with the basis one
# unit higher at every row less the basis as it stands; the pins compare the
# two on a small fit.

set.seed(20261111)
n <- 90L
x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
z <- rbinom(n, 1L, 0.5)
g <- factor(sample(c("a", "b"), n, TRUE))
y <- sin(3 * x2) + (1 + 2 * x1) * z + 0.5 * (g == "b") * x3 + rnorm(n, sd = 0.5)
w <- rep(1, n)
w[c(3L, 40L, 77L)] <- 0
d <- data.frame(y, x1, x2, x3, z, g, w)
rowsIn <- w > 0

fit <- function(formula, data = d, ...) {
  do.call(
    bart,
    c(
      list(
        formula,
        data = data,
        n.chains = 2L,
        n.samples = 40L,
        n.burn = 20L,
        n.trees = 20L,
        n.threads = 1L,
        verbose = FALSE,
        keepTrees = TRUE
      ),
      list(...)
    )
  )
}
f <- fit(
  y ~
    forest(x1 + x2) +
    forest(x1 + x2 + x3, basis = z) +
    forest(x1 + x3, basis = g),
  weights = w
)

# what the two numbers are, from predict's own difference at the basis set
# one unit higher and then to its base value
quantities <- function(difference) {
  difference <- difference[, rowsIn, drop = FALSE]
  average <- rowMeans(difference)
  spread <- apply(difference, 1L, sd)
  q <- function(v) quantile(v, c(0.5, 0.05, 0.95), names = FALSE)
  c(q(average), q(spread))
}
shown <- function(effects, forest) {
  as.numeric(effects[
    effects$forest == forest,
    c(
      "average",
      "average.q5",
      "average.q95",
      "spread",
      "spread.q5",
      "spread.q95"
    )
  ])
}

effects <- summary(f)$effects
expect_identical(effects$forest, 2:3)
expect_identical(attr(effects, "rows"), sum(rowsIn))
expect_identical(attr(effects, "scale"), "in the response's units")

# numeric basis: z one unit up against z at zero
toZ <- function(value) transform(d, z = value)
numericDifference <- predict(f, newdata = toZ(1)) - predict(f, newdata = toZ(0))
expect_equal(
  shown(effects, 2L),
  quantities(numericDifference),
  tolerance = 1e-8
)

# factor basis: the second level against the first
toG <- function(level) transform(d, g = factor(level, levels = c("a", "b")))
factorDifference <- predict(f, newdata = toG("b")) -
  predict(f, newdata = toG("a"))
expect_equal(
  shown(effects, 3L),
  quantities(factorDifference),
  tolerance = 1e-8
)

# the forest with no basis has no line, and the lines are printed once each
printed <- capture.output(print(summary(f)))
at <- grep("^effect of one unit of each basis", printed)
expect_identical(length(at), 1L)
expect_true(grepl("over 87 training rows", printed[at], fixed = TRUE))
expect_true(grepl("5% to 95%", printed[at], fixed = TRUE))
expect_identical(length(printed) - at, 2L)
expect_true(grepl("^  2 z, per unit: average ", printed[at + 1L]))
expect_true(grepl("^  3 g, level b against a: average ", printed[at + 2L]))
expect_false(any(grepl("forest1, ", printed, fixed = TRUE)))

# the average is taken last: the median of the per-draw averages is not the
# average of the median effect over rows
perRow <- apply(numericDifference[, rowsIn], 2L, median)
expect_true(abs(mean(perRow) - effects$average[1L]) > 0)

# a fit that dropped its per-forest fits says so in one line and no number
dropped <- fit(
  y ~ forest(x1 + x2) + forest(x1 + x2 + x3, basis = z),
  keepFits = FALSE
)
expect_identical(summary(dropped)$effects, "dropped")
printedDropped <- capture.output(print(summary(dropped)))
expect_true(any(grepl(
  "effects of the bases not shown: the fit kept no per-forest fits",
  printedDropped,
  fixed = TRUE
)))

# a fit of several forests none with a basis has none
none <- fit(y ~ forest(x1 + x2) + forest(x1 + x3))
expect_null(summary(none)$effects)
expect_false(any(grepl(
  "effect of one unit",
  capture.output(print(summary(none))),
  fixed = TRUE
)))

# a binary response: the effect is on the probit scale the forests add on
yb <- as.double(y > median(y))
db <- transform(d, yb = yb)
fb <- fit(
  yb ~ forest(x1 + x2) + forest(x1 + x2 + x3, basis = z),
  data = db,
  weights = w
)
effectsB <- summary(fb)$effects
expect_identical(attr(effectsB, "scale"), "on the probit scale")
latent <- function(value) {
  predict(fb, newdata = transform(db, z = value), type = "bart")
}
expect_equal(
  shown(effectsB, 2L),
  quantities(latent(1) - latent(0)),
  tolerance = 1e-8
)
expect_true(any(grepl(
  "on the probit scale, over 87 training rows",
  capture.output(print(summary(fb))),
  fixed = TRUE
)))
