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

# a hazard fit's forests add on the probit (hazard) or logit (hazard.logistic)
# scale of the link, and its rows are person-period rows
set.seed(20261112)
nh <- 40L
h <- data.frame(
  time = sample(1:3, nh, TRUE),
  status = rbinom(nh, 1L, 0.7),
  x1 = runif(nh),
  x2 = runif(nh),
  dose = rbinom(nh, 1L, 0.5)
)
for (hazardFamily in c("hazard", "hazard.logistic")) {
  hz <- bart(
    survival::Surv(time, status) ~ x1 + forest(x2, basis = dose),
    data = h,
    family = hazardFamily,
    n.chains = 2L,
    n.samples = 20L,
    n.burn = 10L,
    n.trees = 10L,
    n.threads = 1L,
    verbose = FALSE
  )
  scaleText <- if (hazardFamily == "hazard") "probit" else "logit"
  effectsH <- summary(hz)$effects
  expect_identical(
    attr(effectsH, "scale"),
    paste0("on the ", scaleText, " scale")
  )
  expect_true(any(grepl(
    paste0(
      "on the ",
      scaleText,
      " scale, over ",
      sum(h$time),
      " person-period rows"
    ),
    capture.output(print(summary(hz))),
    fixed = TRUE
  )))
}

# frequency weights: the same data as written-out rows and as counted rows
# give the same average effect, to the draws' Monte Carlo error. The counts
# rise with x1, which the effect does too, so counting each row once reads a
# lower average than the written-out rows do.
set.seed(20261113)
m <- 60L
x1c <- runif(m)
x2c <- runif(m)
x3c <- runif(m)
zc <- rbinom(m, 1L, 0.5)
counts <- ifelse(x1c > 0.5, 20, 1)
yc <- sin(3 * x2c) + (0.5 + 6 * x1c) * zc + rnorm(m, sd = 0.3)
counted <- data.frame(y = yc, x1 = x1c, x2 = x2c, x3 = x3c, z = zc, w = counts)
written <- counted[rep(seq_len(m), counts), ]
formulaC <- y ~ forest(x1 + x2) + forest(x1 + x2, basis = z)
fitCounts <- function(data, ...) {
  do.call(
    bart,
    c(
      list(
        formulaC,
        data = data,
        n.chains = 2L,
        n.samples = 100L,
        n.burn = 50L,
        n.trees = 20L,
        n.threads = 1L,
        verbose = FALSE,
        keepTrees = TRUE
      ),
      list(...)
    )
  )
}
# by hand from one fit's own per-draw, per-row effects (predict's difference at
# the basis one unit up and then at zero): the average and the sd across rows
# at the given row multiplicities, each per draw, medians last
byHand <- function(f, multiplicity) {
  up <- predict(f, newdata = transform(counted, z = 1), type = "bart")
  down <- predict(f, newdata = transform(counted, z = 0), type = "bart")
  effect <- up - down
  total <- sum(multiplicity)
  average <- as.vector(effect %*% multiplicity) / total
  deviation <- sqrt(
    as.vector((effect - average)^2 %*% multiplicity) / (total - 1)
  )
  q <- function(v) quantile(v, c(0.5, 0.05, 0.95), names = FALSE)
  c(q(average), q(deviation))
}
set.seed(2)
fCounted <- fitCounts(counted, weights = counts, weights.type = "frequency")
eCounted <- summary(fCounted)$effects
stored <- function(e) {
  as.numeric(e[
    1L,
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
expect_equal(stored(eCounted), byHand(fCounted, counts), tolerance = 1e-10)
# the printed line carries the same numbers
hand <- byHand(fCounted, counts)
shownText <- function(v) {
  paste0(
    format(v[1L], digits = 4L),
    " (",
    format(v[2L], digits = 4L),
    " to ",
    format(v[3L], digits = 4L),
    ")"
  )
}
printedCounted <- capture.output(print(summary(fCounted)))
expect_true(any(grepl(
  paste0(
    "average ",
    shownText(hand[1:3]),
    ", sd across rows ",
    shownText(hand[4:6])
  ),
  printedCounted,
  fixed = TRUE
)))
# counting each row once reads another number, so the pin discriminates
expect_true(
  max(abs(byHand(fCounted, rep(1, m)) - hand)) > 1e-3
)
expect_identical(attr(summary(fCounted)$effects, "observations"), sum(counts))
expect_true(any(grepl(
  paste0("over ", sum(counts), " observations (", m, " training rows)"),
  capture.output(print(summary(fCounted))),
  fixed = TRUE
)))
# precision weights are not counts: each row stays once
fPrecision <- fitCounts(counted, weights = counts)
expect_identical(
  attr(summary(fPrecision)$effects, "observations"),
  as.double(m)
)
expect_true(any(grepl(
  paste0("over ", m, " training rows"),
  capture.output(print(summary(fPrecision))),
  fixed = TRUE
)))
expect_equal(
  stored(summary(fPrecision)$effects),
  byHand(fPrecision, rep(1, m)),
  tolerance = 1e-10
)

# the forest line's column mean and sd read the counts, as its "taken from a
# column of sd" source does, so the two agree
counted$s <- counted$x3
fLine <- bart(
  y ~ forest(x1 + x2) + forest(x1 + x2, basis = s),
  data = counted,
  weights = w,
  weights.type = "frequency",
  n.chains = 1L,
  n.samples = 10L,
  n.burn = 10L,
  n.trees = 10L,
  n.threads = 1L,
  verbose = FALSE
)
column <- rep(counted$s, counts)
expect_true(any(grepl(
  paste0(
    "numeric basis, mean ",
    format(mean(column), digits = 4L),
    ", sd ",
    format(sd(column), digits = 4L)
  ),
  summary(fLine)$forests,
  fixed = TRUE
)))
expect_false(any(grepl("taken from a column", summary(fLine)$forests)))

# a binomial row counts at its trials times any frequency weight: the cbind
# fit, its written-out Bernoulli twin and its count-weighted twin print the
# same forest block, and the effect header counts the trials
set.seed(20261114)
mb <- 24L
trialsB <- sample(3:30, mb, TRUE)
xb1 <- runif(mb)
xb2 <- runif(mb)
zb <- runif(mb) * 3
succB <- rbinom(mb, trialsB, plogis(-0.5 + xb1 + 0.4 * zb))
cbindData <- data.frame(
  s = succB,
  f = trialsB - succB,
  x1 = xb1,
  x2 = xb2,
  z = zb
)
bernoulli <- data.frame(
  y = unlist(lapply(seq_len(mb), function(i) {
    c(rep(1, succB[i]), rep(0, trialsB[i] - succB[i]))
  })),
  x1 = rep(xb1, trialsB),
  x2 = rep(xb2, trialsB),
  z = rep(zb, trialsB)
)
weighted <- rbind(
  data.frame(y = 1, x1 = xb1, x2 = xb2, z = zb, n = succB),
  data.frame(y = 0, x1 = xb1, x2 = xb2, z = zb, n = trialsB - succB)
)
weighted <- weighted[weighted$n > 0, ]
fitBinomial <- function(formula, data, ...) {
  do.call(
    bart,
    c(
      list(
        formula,
        data = data,
        family = "logistic",
        n.chains = 1L,
        n.samples = 100L,
        n.burn = 50L,
        n.trees = 20L,
        n.threads = 1L,
        seed = 1L,
        verbose = FALSE,
        keepTrees = TRUE
      ),
      list(...)
    )
  )
}
rhsB <- ~ forest(x1 + x2) + forest(x1 + x2, basis = z)
fCbind <- fitBinomial(update(cbind(s, f) ~ ., rhsB), cbindData)
fBernoulli <- fitBinomial(update(y ~ ., rhsB), bernoulli)
fWeighted <- fitBinomial(
  update(y ~ ., rhsB),
  weighted,
  weights = quote(n),
  weights.type = "frequency"
)
forestsOf <- function(f) summary(f)$forests
expect_identical(forestsOf(fCbind), forestsOf(fBernoulli))
expect_identical(forestsOf(fCbind), forestsOf(fWeighted))
expect_false(any(grepl("taken from a column", forestsOf(fCbind))))
sumTrials <- sum(trialsB)
expect_identical(
  attr(summary(fCbind)$effects, "observations"),
  as.double(sumTrials)
)
expect_identical(
  attr(summary(fWeighted)$effects, "observations"),
  as.double(sumTrials)
)
expect_true(any(grepl(
  paste0("over ", sumTrials, " trials (", mb, " training rows)"),
  capture.output(print(summary(fCbind))),
  fixed = TRUE
)))
expect_true(any(grepl(
  paste0("over ", sumTrials, " training rows"),
  capture.output(print(summary(fBernoulli))),
  fixed = TRUE
)))
# the effect average is the trials-weighted one, not the once-per-row one
upB <- predict(
  fCbind,
  newdata = transform(cbindData, z = cbindData$z + 1),
  type = "bart"
)
downB <- predict(fCbind, newdata = cbindData, type = "bart")
effectB <- upB - downB
perDraw <- function(weightsB) as.vector(effectB %*% weightsB) / sum(weightsB)
expect_equal(
  summary(fCbind)$effects$average,
  median(perDraw(trialsB)),
  tolerance = 1e-8
)
expect_true(abs(median(perDraw(trialsB)) - median(perDraw(rep(1, mb)))) > 1e-4)

# the counted header keeps the row noun and the singular
effectsHeader <- function(observations, row.name, rows, trials = FALSE) {
  effects <- data.frame(
    forest = 1L,
    label = "a",
    effect = "z",
    average = 1,
    average.q5 = 0,
    average.q95 = 2,
    spread = 1,
    spread.q5 = 0,
    spread.q95 = 2
  )
  attr(effects, "rows") <- rows
  attr(effects, "observations") <- observations
  attr(effects, "row.name") <- row.name
  attr(effects, "trials") <- trials
  attr(effects, "scale") <- "in the response's units"
  dbarts:::forestEffectLines(effects)[1L]
}
expect_true(grepl(
  "over 314 person-period observations (163 person-period rows)",
  effectsHeader(314, "person-period rows", 163L),
  fixed = TRUE
))
expect_true(grepl(
  "over 5 observations (1 training row)",
  effectsHeader(5, "training row", 1L),
  fixed = TRUE
))
expect_true(grepl(
  "over 1 observation (2 training rows)",
  effectsHeader(1, "training rows", 2L),
  fixed = TRUE
))
expect_true(grepl(
  "over 1 trial (2 training rows)",
  effectsHeader(1, "training rows", 2L, TRUE),
  fixed = TRUE
))

# a binary, binomial or hazard fit adds a change-in-probability line under each
# link-scale line, per draw and row from the kept fits and averaged last: a
# factor basis the change from its first level to its second, a numeric basis
# the slope of the probability per unit. A gaussian fit has none.
probabilityOf <- function(effects, forest) {
  as.numeric(effects[
    effects$forest == forest,
    paste0(
      "probability.",
      c(
        "average",
        "average.q5",
        "average.q95",
        "spread",
        "spread.q5",
        "spread.q95"
      )
    )
  ])
}
expect_null(effects[["probability.kind"]])
expect_false(any(grepl("in probability|of probability", printed)))

# probit: a numeric basis against the finite difference of predict's
# probabilities (exact to rounding, the basis entering linearly), a factor
# against predict's probabilities at its two levels
db$g <- g
fProbit <- fit(
  yb ~
    forest(x1 + x2) +
    forest(x1 + x2 + x3, basis = z) +
    forest(x1 + x3, basis = g),
  data = db,
  weights = w
)
effectsP <- summary(fProbit)$effects
step <- 1e-6
evAt <- function(data) predict(fProbit, newdata = data, type = "ev")
slopeDraws <- (evAt(transform(db, z = z + step)) - evAt(db)) / step
expect_equal(
  probabilityOf(effectsP, 2L),
  quantities(slopeDraws),
  tolerance = 1e-4
)
levelAt <- function(level) {
  evAt(transform(db, g = factor(level, levels = c("a", "b"))))
}
expect_equal(
  probabilityOf(effectsP, 3L),
  quantities(levelAt("b") - levelAt("a")),
  tolerance = 1e-8
)
expect_identical(effectsP$probability.kind, c("slope", "change"))
printedP <- capture.output(print(summary(fProbit)))
atP <- grep("^effect of one unit of each basis", printedP)
expect_identical(length(printedP) - atP, 4L)
expect_true(grepl(
  "^    slope of probability per unit: average ",
  printedP[atP + 2L]
))
expect_true(grepl(
  "^    change in probability, level b against a: average ",
  printedP[atP + 4L]
))
expect_true(grepl(
  paste0(
    "average ",
    shownText(probabilityOf(effectsP, 2L)[1:3]),
    ", sd across rows ",
    shownText(probabilityOf(effectsP, 2L)[4:6])
  ),
  printedP[atP + 2L],
  fixed = TRUE
))
# a hazard fit, by hand from its own kept fits: the link's distribution at
# the draw's latent fit, the amplitude and the forest's fit
for (hazardFamily in c("hazard", "hazard.logistic")) {
  cdf <- if (hazardFamily == "hazard") pnorm else plogis
  density <- if (hazardFamily == "hazard") dnorm else dlogis
  h$doseF <- factor(h$dose, levels = c(0, 1))
  hz <- bart(
    survival::Surv(time, status) ~
      x1 +
      forest(x2, basis = dose) +
      forest(x2, basis = doseF),
    data = h,
    family = hazardFamily,
    n.chains = 2L,
    n.samples = 20L,
    n.burn = 10L,
    n.trees = 10L,
    n.threads = 1L,
    verbose = FALSE
  )
  effectsH <- summary(hz)$effects
  eta <- hz$yhat.train
  fits <- hz$forestFits
  glue <- hz$glue
  forestNames <- dimnames(fits)[[3L]]
  glueOf <- function(k) {
    glue[, attr(glue, "forest") == forestNames[k], drop = FALSE]
  }
  qH <- function(difference) {
    q <- function(v) quantile(v, c(0.5, 0.05, 0.95), names = FALSE)
    c(q(rowMeans(difference)), q(apply(difference, 1L, sd)))
  }
  # numeric basis: the density at the latent fit times the forest's
  # link-scale effect, per draw and row
  slopeH <- density(eta) * (glueOf(2L)[, 1L] * fits[,, 2L])
  expect_equal(probabilityOf(effectsH, 2L), qH(slopeH), tolerance = 1e-10)
  # factor basis: the row's term moved to the first level and to the second
  level <- as.integer(hz$bases[[3L]])
  gl <- glueOf(3L)
  atLevel <- function(column) {
    cdf(eta + (gl[, column] - gl[, level]) * fits[,, 3L])
  }
  expect_equal(
    probabilityOf(effectsH, 3L),
    qH(atLevel(2L) - atLevel(1L)),
    tolerance = 1e-10
  )
  expect_identical(effectsH$probability.kind, c("slope", "change"))
  printedH <- capture.output(print(summary(hz)))
  expect_true(any(startsWith(
    printedH,
    "forests, sd on the latent scale (error sd "
  )))
  expect_false(any(grepl("response's units", printedH, fixed = TRUE)))
}

quantities2 <- function(difference, rows) {
  q <- function(v) quantile(v, c(0.5, 0.05, 0.95), names = FALSE)
  c(q(rowMeans(difference[, rows])), q(apply(difference[, rows], 1L, sd)))
}

# an offset enters the latent fit the probability is taken at: the line equals
# predict's probabilities with the same offset, for a numeric and a factor basis
set.seed(20261120)
no <- 80L
o <- data.frame(
  x1 = runif(no),
  x2 = runif(no),
  z = runif(no),
  g = factor(sample(c("a", "b"), no, TRUE)),
  off = rnorm(no, 0, 0.3)
)
o$y <- rbinom(no, 1L, pnorm(o$x1 + o$z * o$x2 + o$off))
fOffset <- bart(
  y ~
    forest(x1) +
    forest(x1 + x2, basis = z) +
    forest(x2, basis = g),
  data = o,
  offset = o$off,
  family = "probit",
  n.chains = 2L,
  n.samples = 40L,
  n.burn = 20L,
  n.trees = 15L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrees = TRUE
)
effectsO <- summary(fOffset)$effects
evO <- function(data) {
  predict(fOffset, newdata = data, type = "ev", offset = o$off)
}
expect_equal(
  probabilityOf(effectsO, 2L),
  quantities2(
    (evO(transform(o, z = z + 1e-6)) - evO(o)) / 1e-6,
    rep(TRUE, no)
  ),
  tolerance = 1e-4
)
levelO <- function(level) {
  evO(transform(o, g = factor(level, levels = c("a", "b"))))
}
expect_equal(
  probabilityOf(effectsO, 3L),
  quantities2(levelO("b") - levelO("a"), rep(TRUE, no)),
  tolerance = 1e-8
)

# a hazard fit's offset (stored per person-period row) enters the same way
set.seed(20261121)
ho <- data.frame(
  time = sample(1:3, 40L, TRUE),
  status = rbinom(40L, 1L, 0.7),
  x1 = runif(40L),
  x2 = runif(40L),
  dose = runif(40L),
  off = rnorm(40L, 0, 0.3)
)
hzO <- bart(
  survival::Surv(time, status) ~ forest(x1) + forest(x2, basis = dose),
  data = ho,
  offset = ho$off,
  family = "hazard",
  n.chains = 2L,
  n.samples = 20L,
  n.burn = 10L,
  n.trees = 10L,
  n.threads = 1L,
  verbose = FALSE
)
expect_identical(length(hzO$offset), ncol(hzO$yhat.train))
etaO <- hzO$yhat.train + rep(hzO$offset, each = nrow(hzO$yhat.train))
glueO <- hzO$glue[,
  attr(hzO$glue, "forest") == dimnames(hzO$forestFits)[[3L]][2L]
]
slopeO <- dnorm(etaO) * (glueO * hzO$forestFits[,, 2L])
expect_equal(
  probabilityOf(summary(hzO)$effects, 2L),
  quantities2(slopeO, rep(TRUE, ncol(slopeO))),
  tolerance = 1e-10
)
