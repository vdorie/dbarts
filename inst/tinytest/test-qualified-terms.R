# A formula term written with a package prefix (splines::ns(), stats::poly(),
# base::log()) is a call, not an interaction: every door that reads a formula
# fits it and draws what the unprefixed term draws. A real interaction is
# still refused.

library(splines)

set.seed(1L)
n <- 40L
d <- data.frame(w = runif(n, 1, 2), x2 = runif(n), g = gl(4L, 10L))
d$y <- d$w + rnorm(n)
d$f <- factor(sample(letters[1:3], n, TRUE))

pairs <- list(
  list(y ~ splines::ns(w, 3) + x2, y ~ ns(w, 3) + x2),
  list(y ~ stats::poly(w, 2) + x2, y ~ poly(w, 2) + x2),
  list(y ~ base::log(w) + x2, y ~ log(w) + x2),
  list(y ~ splines:::ns(w, 3) + x2, y ~ ns(w, 3) + x2)
)

sameDraws <- function(fit, pair, door = "") {
  prefixed <- fit(pair[[1L]])
  plain <- fit(pair[[2L]])
  expect_equal(prefixed, plain, info = paste(door, deparse1(pair[[1L]])))
}

termOf <- function(formula) attr(terms(formula), "term.labels")[[1L]]

bartArgs <- list(
  n.trees = 5L,
  n.samples = 5L,
  n.burn = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  seed = 1L
)
fitBart <- function(formula, ...) {
  do.call(bart, c(list(formula, d), bartArgs, list(...)))$yhat.train
}

for (pair in pairs) {
  # dbartsData
  sameDraws(
    function(f) unname(as.matrix(dbartsData(f, d, keep.data = FALSE)@x)),
    pair,
    door = "dbartsData"
  )
  # bart
  sameDraws(fitBart, pair, door = "bart")
  # bart with a variance formula that names the term as the fit does
  sameDraws(
    function(f) fitBart(f, variance = reformulate(termOf(f))),
    pair,
    door = "variance"
  )
  # dbarts
  sameDraws(
    function(f) {
      s <- dbarts(
        f,
        d,
        control = dbartsControl(n.threads = 1L, verbose = FALSE, seed = 1L),
        n.trees = 5L
      )
      s$run(2L, 3L)$train
    },
    pair,
    door = "dbarts"
  )
  # xbart
  sameDraws(
    function(f) {
      set.seed(1L)
      xbart(
        f,
        d,
        n.trees = 5L,
        n.reps = 1L,
        n.test = 2L,
        n.samples = 5L,
        n.burn = 2L,
        n.threads = 1L,
        verbose = FALSE
      )
    },
    pair,
    door = "xbart"
  )
  # rbart_vi
  sameDraws(
    function(f) {
      fit <- suppressWarnings(rbart_vi(
        f,
        d,
        group.by = g,
        n.trees = 5L,
        n.samples = 5L,
        n.burn = 2L,
        n.chains = 1L,
        n.threads = 1L,
        verbose = FALSE,
        seed = 1L
      ))
      fit$yhat.train
    },
    pair,
    door = "rbart_vi"
  )
  # pdbart
  sameDraws(
    function(f) {
      set.seed(1L)
      pdbart(
        f,
        d,
        xind = "x2",
        pl = FALSE,
        n.trees = 5L,
        n.threads = 1L,
        verbose = FALSE
      )$fd
    },
    pair,
    door = "pdbart"
  )
}

# the multinomial formula path
expect_equal(
  bart(
    f ~ splines::ns(w, 3) + x2,
    d,
    n.samples = 3L,
    n.burn = 3L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE,
    seed = 1L
  )$yhat.train,
  bart(
    f ~ ns(w, 3) + x2,
    d,
    n.samples = 3L,
    n.burn = 3L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE,
    seed = 1L
  )$yhat.train
)

# a prefixed call beside a real interaction: the interaction is refused
for (formula in list(
  y ~ splines::ns(w, 3):x2,
  y ~ x2:stats::poly(w, 2),
  y ~ base::log(w) * x2,
  y ~ w:x2 + splines::ns(w, 3)
)) {
  expect_error(
    dbartsData(formula, d, keep.data = FALSE),
    "':' and '*' terms are not supported",
    fixed = TRUE,
    info = deparse1(formula)
  )
  expect_error(
    fitBart(formula),
    "':' and '*' terms are not supported",
    fixed = TRUE,
    info = deparse1(formula)
  )
}
rm(d, pairs, termOf, sameDraws, bartArgs, fitBart, n, pair, formula)
