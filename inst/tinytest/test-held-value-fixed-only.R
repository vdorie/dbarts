# a held sigma or count shape is stored only in fit$fixed, as a held k is: the
# draw channels (fit$sigma, fit$first.sigma, fit$shape) are absent, and the
# package's readers take the held value from fit$fixed; fit$fixed also carries
# a forest coefficient held by amplitude = fixed(), per forest (dec-B376)

set.seed(12)
n <- 50L
x <- matrix(runif(n * 3L), n, dimnames = list(NULL, c("a", "b", "c")))
y <- x[, 1L] + rnorm(n, 0, 0.3)
quick <- function(...) {
  suppressWarnings(bart(
    ...,
    n.trees = 10L,
    n.burn = 10L,
    n.samples = 10L,
    n.chains = 2L,
    n.threads = 1L,
    keepTrees = TRUE,
    verbose = FALSE
  ))
}
counted <- function(expr) {
  seen <- 0L
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      seen <<- seen + 1L
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = seen)
}

# gaussian: sigma held
fit <- quick(x, y, family = gaussian(sigma = fixed(0.5)))
expect_null(fit[["sigma"]])
expect_null(fit[["first.sigma"]])
expect_true(is.numeric(fit$fixed$sigma) && length(fit$fixed$sigma) == 1L)
held <- fit$fixed$sigma
expect_equal(extract(fit, "sigma"), held)
# the readers take the held value: the log-likelihood is the normal density at
# it, the predictive draw has its spread, summary names it, plot draws on
ev <- extract(fit, "ev")
expect_equal(
  as.vector(extract(fit, "loglik")),
  dnorm(rep(y, each = nrow(ev)), as.vector(ev), held, log = TRUE)
)
res <- counted(extract(fit, "ppd"))
expect_identical(res$warnings, 0L)
spread <- sd(as.vector(res$value - extract(fit, "ev")))
expect_true(abs(spread - held) < 0.2 * held)
expect_identical(dim(predict(fit, x, type = "ppd")), c(20L, n))
expect_false(anyNA(fitted(fit, type = "ppd")))
expect_identical(counted(summary(fit))$warnings, 0L)
expect_identical(counted({
  pdf(NULL)
  on.exit(dev.off())
  plot(fit)
})$warnings, 0L)
# a drawn sigma keeps its channels
fit.drawn <- quick(x, y)
expect_false(is.null(fit.drawn[["sigma"]]))
expect_false(is.null(fit.drawn[["first.sigma"]]))
expect_null(fit.drawn$fixed$sigma)

# aft: sigma held, survival curves and log-likelihood read it
tm <- exp(y)
fit.aft <- quick(x, cbind(tm, rep(1, n)), family = aft(sigma = fixed(0.7)))
expect_null(fit.aft[["sigma"]])
expect_null(fit.aft[["first.sigma"]])
expect_true(is.numeric(fit.aft$fixed$sigma))
expect_identical(
  dim(survivalProbabilities(fit.aft, c(1, 2))),
  c(20L, 2L, n)
)
expect_false(anyNA(extract(fit.aft, "loglik")))
lp <- extract(fit.aft, "bart", combineChains = FALSE)
expect_equal(
  unname(survivalProbabilities(fit.aft, 1.5, combineChains = FALSE)[, , 1L, ]),
  unname(pnorm((log(1.5) - lp) / fit.aft$fixed$sigma, lower.tail = FALSE))
)

# count shape held
cnt <- rpois(n, 3)
fit.nb <- quick(x, cnt, family = nbinom(shape = 3))
expect_null(fit.nb[["shape"]])
expect_null(fit.nb[["shape.raw"]])
expect_equal(fit.nb$fixed$shape, 3)
expect_identical(dim(extract(fit.nb, "ppd")), c(20L, n))
expect_identical(dim(predict(fit.nb, x, type = "ppd")), c(20L, n))
expect_false(anyNA(extract(fit.nb, "loglik")))
expect_equal(
  mean(as.vector(predict(fit.nb, x, type = "ev"))),
  mean(as.vector(extract(fit.nb, "ev")))
)
# predict still needs the saved trees
fit.nb.bare <- suppressWarnings(bart(
  x,
  cnt,
  family = nbinom(shape = 3),
  n.trees = 10L,
  n.burn = 10L,
  n.samples = 10L,
  n.chains = 2L,
  n.threads = 1L,
  keepTrees = FALSE,
  verbose = FALSE
))
expect_error(
  predict(fit.nb.bare, x[1:3, ]),
  "requires the fit's saved trees"
)
expect_identical(dim(predict(fit.nb.bare)), c(20L, n))

# a forest held by amplitude = fixed() is recorded per forest
df <- data.frame(x, z = rep(c(0, 1), length.out = n))
df$y <- y
fit.bcf <- suppressWarnings(bart(
  y ~ a + b + c + forest(a + b, basis = ~ factor(z), amplitude = fixed()),
  df,
  n.trees = 10L,
  n.burn = 10L,
  n.samples = 10L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
))
expect_identical(fit.bcf$fixed$amplitude, c(forest2 = 1))
fit.free <- suppressWarnings(bart(
  y ~ a + b + c + forest(a + b, basis = ~ factor(z)),
  df,
  n.trees = 10L,
  n.burn = 10L,
  n.samples = 10L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
))
expect_null(fit.free$fixed$amplitude)
