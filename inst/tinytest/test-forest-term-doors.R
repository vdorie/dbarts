# The doors that fit no model of several forests refuse a forest() term by
# name, with the term, in place of whatever their own reading of the formula
# would stop on.

n <- 40L
set.seed(3)
d <- data.frame(x1 = runif(n), x2 = runif(n), z = rbinom(n, 1L, 0.5))
d$g <- rep_len(c("p", "q", "r"), n)
d$y <- d$x1 + d$z * d$x2 + rnorm(n, 0, 0.2)

doorText <- function(door) {
  paste0(
    door,
    "() does not take a forest() term ('forest(x1 + x2)'): the forests of a ",
    "model are written in the formula of bart() or dbarts(), or in a ",
    "'forests' list"
  )
}
expect_error(
  dbarts::xbart(y ~ forest(x1 + x2), d, n.threads = 1L),
  doorText("xbart"),
  fixed = TRUE
)
expect_error(
  dbarts::xbart(y ~ x1 + forest(x1 + x2, basis = ~z), d, n.threads = 1L),
  "xbart() does not take a forest() term ('forest(x1 + x2, basis = ~z)')",
  fixed = TRUE
)
expect_error(
  suppressWarnings(
    dbarts::rbart_vi(y ~ forest(x1 + x2), d, group.by = g, n.threads = 1L),
    classes = "dbartsDeprecatedWarning"
  ),
  doorText("rbart_vi"),
  fixed = TRUE
)
expect_error(
  dbarts::bartBT(y ~ forest(x1 + x2), d, verbose = FALSE),
  doorText("bartBT"),
  fixed = TRUE
)
expect_error(
  dbarts::pdbart(y ~ forest(x1 + x2), d, xind = 1L, pl = FALSE),
  doorText("pdbart"),
  fixed = TRUE
)
expect_error(
  dbarts::pd2bart(y ~ forest(x1 + x2), d, xind = 1:2, pl = FALSE),
  doorText("pd2bart"),
  fixed = TRUE
)
expect_error(
  dbarts::dbartsData(y ~ forest(x1 + x2), d),
  doorText("dbartsData"),
  fixed = TRUE
)
