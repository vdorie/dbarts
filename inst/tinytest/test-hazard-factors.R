# A formula hazard fit keeps the categorical design through the person-period
# expansion: it equals the probit fit on the hand-expanded frame, and every
# frame reader works.

suppressMessages(library(survival))

set.seed(1L)
n <- 60L
d <- data.frame(
  time = sample(1:5, n, TRUE),
  status = rbinom(n, 1L, 0.7),
  g = factor(c("a", "b", "c")[sample(3L, n, TRUE)]),
  z = rnorm(n),
  extra = letters[seq_len(n) %% 26L + 1L]
)
d$o <- factor(as.character(d$g), ordered = TRUE)
d$gq <- factor(as.character(d$g), levels = c("a", "b", "c", "q"))
d$sf <- sparseFactor(ifelse(runif(n) < 0.15, "x", "ref"), reference = "ref")
args <- list(
  keepTrees = TRUE,
  verbose = FALSE,
  seed = 1L,
  n.samples = 30L,
  n.burn = 30L,
  n.chains = 1L,
  n.threads = 1L
)
subj <- rep(seq_len(n), d$time)
per <- sequence(d$time)
e <- d[subj, ]
e$period <- per
e$y <- as.integer(per == d$time[subj] & d$status[subj] == 1)
rownames(e) <- NULL
K <- max(d$time)

hazardFit <- function(rhs, ...) {
  do.call(
    bart,
    c(
      list(
        as.formula(paste("Surv(time, status) ~", rhs)),
        data = d,
        family = "hazard"
      ),
      args,
      list(...)
    )
  )
}
byHand <- function(rhs, data = e) {
  do.call(
    bart,
    c(
      list(as.formula(paste("y ~", rhs, "+ period")), data = data),
      args,
      list(family = "probit")
    )
  )
}

# ---- training draws equal the hand-expanded fit ----
for (rhs in c("g + z", "o + z", "gq + z", "g + log(abs(z))")) {
  fh <- hazardFit(rhs)
  expect_identical(
    unname(fh$yhat.train),
    unname(byHand(rhs)$yhat.train),
    info = rhs
  )
}
fh.g <- hazardFit("g + z")
expect_identical(
  attr(fh.g$fit$data@x, "factor.levels")[[1L]],
  levels(d$g)
)
expect_identical(
  tail(fh.g$fit$data@varTypes, 1L),
  0L
)
expect_identical(colnames(fh.g$fit$data@x)[ncol(fh.g$fit$data@x)], "period")

# the sparse factor's design order is z, sf, period
fh.sf <- hazardFit("sf + z")
fb.sf <- do.call(
  bart,
  c(list(e[c("z", "sf", "period")], e$y, family = "probit"), args)
)
expect_identical(unname(fh.sf$yhat.train), unname(fb.sf$yhat.train))

# ---- readers take a data frame ----
sp0 <- survivalProbabilities(fh.g)
sp.full <- survivalProbabilities(fh.g, newdata = d)
sp.cols <- survivalProbabilities(fh.g, newdata = d[c("g", "z")])
expect_identical(sp.full, sp.cols)
expect_equal(sp.full, sp0, tolerance = 1e-12)
sp0.sf <- survivalProbabilities(fh.sf)
expect_equal(
  survivalProbabilities(fh.sf, newdata = d),
  sp0.sf,
  tolerance = 1e-12
)
pr <- predict(fh.g, newdata = e, type = "link")
expect_equal(unname(pr), unname(fh.g$yhat.train), tolerance = 1e-12)
expect_equal(nrow(fh.g$fit$getTrees(newdata = e[1:4, ])) > 0L, TRUE)
expect_equal(ncol(predict(fh.g, newdata = e[1:4, ])), 4L)

# ---- test= at fit time ----
for (rhs in c("g + z", "sf + z")) {
  ft <- hazardFit(rhs, test = d[1:5, ])
  test.rows <- d[rep(1:5, K), ]
  test.rows$period <- rep(seq_len(K), each = 5L)
  expect_equal(
    unname(ft$yhat.test),
    unname(predict(ft, newdata = test.rows, type = "link")),
    tolerance = 1e-12,
    info = rhs
  )
}

# ---- factors = "indicators" ----
fi <- hazardFit("g + z", factors = "indicators")
expect_true(is.array(survivalProbabilities(fi, newdata = d[names(d) != "sf"])))

# ---- unseen level and NA factor behave as on aft ----
dn <- d
dn$g <- as.character(dn$g)
dn$g[1L] <- "zz"
expect_error(
  survivalProbabilities(fh.g, newdata = dn),
  "level|unseen|zz"
)

dn$g <- factor(as.character(d$g), levels = levels(d$g))
dn$g[2L] <- NA
expect_error(
  survivalProbabilities(fh.g, newdata = dn),
  "missing values in 'g'"
)
sp.na <- survivalProbabilities(fh.g, newdata = dn, na.action = na.pass)
expect_true(all(is.na(sp.na[, 2L, ])) || all(is.na(sp.na[,, 2L])))

# ---- a predictor named period is refused ----
d.p <- d
d.p$period <- runif(n)
expect_error(
  bart(
    Surv(time, status) ~ z + period,
    data = d.p,
    family = "hazard",
    verbose = FALSE
  ),
  "period"
)
expect_error(
  bart(
    d.p[c("z", "period")],
    Surv(d$time, d$status),
    family = "hazard",
    verbose = FALSE
  ),
  "period"
)
expect_error(
  bart(
    as.matrix(d.p[c("z", "period")]),
    Surv(d$time, d$status),
    family = "hazard",
    verbose = FALSE
  ),
  "period"
)

# a term that reads period is refused, transformed or indicator-coded
d.p$pf <- factor(d.p$period > 0.5)
for (rhs in c("log(period + 2) + z", "I(period^2) + z", "pf + period + z")) {
  expect_error(
    bart(
      as.formula(paste("Surv(time, status) ~", rhs)),
      data = d.p,
      family = "hazard",
      verbose = FALSE
    ),
    "period",
    info = rhs
  )
}
expect_error(
  bart(
    Surv(time, status) ~ period + z,
    data = d.p,
    family = "hazard",
    factors = "indicators",
    verbose = FALSE
  ),
  "period"
)
d.p$period <- factor(d.p$pf)
expect_error(
  bart(
    Surv(time, status) ~ period + z,
    data = d.p,
    family = "hazard",
    factors = "indicators",
    verbose = FALSE
  ),
  "period"
)

# ---- the x/y frame path equals the formula path ----
fxy <- do.call(
  bart,
  c(list(d[c("g", "z")], Surv(d$time, d$status), family = "hazard"), args)
)
expect_identical(unname(fxy$yhat.train), unname(fh.g$yhat.train))
