# The bart2 argument-consolidation arc's maintained contracts
# (docs/plans/bart2-argument-consolidation.md 4.6): T-A (dbartsControl
# formal parity), T-B (shared-default text), T-C (bart -> bart2 concept
# map) and T-E (the per-forest reconstruction identity) arrive with the
# slices that build what they assert against (S1/S5, S3, S10, S11). This
# file starts with T-D, S2's contract: family gating (3.c).

# T-D. One diagnosis test per row of 3.c.2: the classed warning fires under
# a family that ignores the argument and not under one that uses it; plus
# monotonicity assertions (N7) that already-loud refusals keep their
# severity and message.

set.seed(101)
n <- 40L
x <- matrix(rnorm(n * 2L), n)
quick <- list(
  n.trees = 3L,
  n.samples = 5L,
  n.burn = 2L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)

y.gaussian <- rnorm(n)
y.binary <- rbinom(n, 1L, 0.5)
y.multi <- factor(sample(letters[1:3], n, replace = TRUE))
y.ordinal <- factor(sample(1:3, n, replace = TRUE), ordered = TRUE)
y.count <- rpois(n, 4)
y.aft <- cbind(abs(rnorm(n)) + 0.1, rep(1L, n))

fit2 <- function(y, ...) {
  do.call(dbarts::bart2, c(list(x, y), quick, list(...)))
}

# one inventory row x one representative family/site each, class asserted
expect_warning(
  fit2(y.multi, family = "multinomial", sigest = 5),
  pattern = "sigest",
  class = "dbartsFamilyGatedWarning"
)
expect_warning(
  fit2(y.ordinal, family = "ordinal", dispersion = 2),
  pattern = "dispersion",
  class = "dbartsFamilyGatedWarning"
)
expect_warning(
  fit2(y.count, family = "nbinom", breaks = 5),
  pattern = "breaks",
  class = "dbartsFamilyGatedWarning"
)

# counter-tests: live families stay silent
expect_silent(fit2(y.gaussian, family = "gaussian", sigest = 5))
expect_silent(fit2(y.aft, family = "aft", sigest = 5))

# the suppliedness bit: a defaulted (unsupplied) argument never gates
expect_silent(fit2(y.binary, family = "probit"))

# two gated names supplied together still emit exactly one warning
countWarnings <- function(expr, class) {
  count <- 0L
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (inherits(w, class)) {
        count <<- count + 1L
      }
      invokeRestart("muffleWarning")
    }
  )
  count
}
expect_equal(
  countWarnings(
    fit2(y.binary, family = "probit", sigest = 5, dispersion = 2),
    "dbartsFamilyGatedWarning"
  ),
  1L
)

# samplerOnly = TRUE still warns (the standard-site reordering)
expect_warning(
  fit2(y.binary, family = "probit", sigest = 5, samplerOnly = TRUE),
  class = "dbartsFamilyGatedWarning"
)

# hurdle.lognormal: sigest is live on the positive half, so the rule's own
# scoping ("an argument whose only effect is on a family this fit is not")
# makes a diagnosis against the occupancy half's forced "probit" a FALSE
# one; bart2Hurdle strips it from that component call rather than let it
# leak. No warning at all - not even the occupancy component's.
y.hurdle <- c(rep(0, n / 2L), abs(rnorm(n / 2L)) + 0.1)
expect_equal(
  countWarnings(
    fit2(y.hurdle, family = "hurdle.lognormal", sigest = 5),
    "dbartsFamilyGatedWarning"
  ),
  0L
)

# dispersion is inert on both components; the outer site diagnoses it once
# and bart2Hurdle strips it from both component calls, so exactly one
# warning reaches the caller, naming the outer family
expect_equal(
  countWarnings(
    fit2(y.hurdle, family = "hurdle.lognormal", dispersion = 2),
    "dbartsFamilyGatedWarning"
  ),
  1L
)
expect_warning(
  fit2(y.hurdle, family = "hurdle.lognormal", dispersion = 2),
  pattern = "hurdle.lognormal",
  class = "dbartsFamilyGatedWarning"
)

# bart() never calls the helper - silence is preserved by construction
expect_silent(dbarts::bart(
  x,
  y.binary,
  sigest = 5,
  ntree = 3L,
  ndpost = 5L,
  nskip = 2L,
  verbose = FALSE
))

# rbart_vi warns on its resolved probit path
group <- factor(rep(1:4, length.out = n))
expect_warning(
  dbarts::rbart_vi(
    x,
    factor(y.binary),
    group.by = group,
    sigest = 5,
    n.thin = 1L,
    n.trees = 3L,
    n.samples = 5L,
    n.burn = 2L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  ),
  class = "dbartsFamilyGatedWarning"
)

# monotonicity (N7): already-loud refusals keep their severity and message
expect_error(
  fit2(y.multi, family = "multinomial", samplerOnly = TRUE),
  pattern = "does not support 'samplerOnly'"
)
expect_error(
  fit2(y.multi, family = "multinomial", weights = rep(1, n)),
  pattern = "does not support 'weights'"
)
expect_error(
  fit2(y.multi, family = "multinomial", prior.scale = 2),
  pattern = "prior.scale"
)
