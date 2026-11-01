# Small items on the leaf prior's reader and its k: the reader of a sampler of
# one forest carries `sd` and `sd.stated` beside the statement and calls the
# size the forest's; an integer k is the same number as a double; the creation
# summary names what the caller stated; a Student-t forest with no sd refuses
# a leaf prior whose median size is no number, whichever spelling; and what
# $setLeafPrior does to such a forest.

normal <- dbartsPriors$normal
chi <- dbartsPriors$chi
invchi <- dbartsPriors$invchi
chisq <- dbartsPriors$chisq
forest <- dbartsForests$forest

set.seed(20261101)
n <- 60L
x <- cbind(x1 = runif(n), x2 = runif(n))
z <- as.double(seq_len(n) %in% sample(n, 25L))
y <- 1 + x[, 1L] + z * x[, 2L] + rnorm(n)
d <- data.frame(y, x1 = x[, 1L], x2 = x[, 2L], z)
ctl <- function(verbose = FALSE) {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 10L,
    n.samples = 3L,
    n.burn = 0L,
    verbose = verbose,
    seed = 5L,
    updateState = FALSE
  )
}
lone <- function(...) dbarts(x, y, control = ctl(), ...)

# --- one forest: the statement, its sd and whether it is stated --------------
reader <- lone()$getLeafPrior()
expect_identical(
  names(reader),
  c(
    "leaf.prior",
    "sd",
    "sd.stated",
    "leaf.model",
    "prior.sd.of",
    "prior.mean",
    "k.scale",
    "response.scale",
    "response.shift"
  )
)
# still one record, not a list per forest
expect_false(is.list(reader[[1L]]))
expect_identical(reader$prior.sd.of, "forest")
# a fixed k gives the data's scale over it
expect_identical(reader$sd, reader$k.scale / 2)
expect_identical(reader$sd.stated, FALSE)
atK <- lone(leaf.prior = normal(k = 3))$getLeafPrior()
expect_identical(atK$sd, atK$k.scale / 3)
expect_identical(atK$sd.stated, FALSE)
# a stated sd reads as stated, at creation and through the setter alike
created <- lone(leaf.prior = normal(sd = 0.5))$getLeafPrior()
sampler <- lone()
sampler$setLeafPrior(normal(sd = 0.5))
written <- sampler$getLeafPrior()
expect_identical(created$sd, 0.5)
expect_identical(created$sd.stated, TRUE)
expect_identical(written$sd, created$sd)
expect_identical(written$sd.stated, created$sd.stated)
expect_identical(written$leaf.prior, created$leaf.prior)
# and a k written over it reads as unstated
sampler$setLeafPrior(normal(k = 3))
expect_identical(sampler$getLeafPrior()$sd.stated, FALSE)
expect_identical(sampler$getLeafPrior()$sd, atK$sd)
# a number whose quotient by k rounds elsewhere is still read as stated
scaleOf <- reader$k.scale
awkward <- Filter(function(v) scaleOf / (scaleOf / v) != v, seq(0.05, 3, 0.05))
expect_true(length(awkward) > 0L)
awkwardReader <- lone(leaf.prior = normal(sd = awkward[[1L]]))$getLeafPrior()
expect_identical(awkwardReader$sd, awkward[[1L]])
# a lone forest stating its sd is the same record
viaForest <- dbarts(y ~ forest(x1 + x2, sd = 0.5), d, control = ctl())
expect_identical(viaForest$getLeafPrior()$sd, 0.5)
expect_identical(viaForest$getLeafPrior()$sd.stated, TRUE)
# a drawn k, or an sd under invchi(), has no one number: sd is NA, as a
# sharer's is, and sd.stated says which spelling was used
drawnK <- lone(leaf.prior = normal(k = chi(1.5, 2)))$getLeafPrior()
expect_identical(drawnK$sd, NA_real_)
expect_identical(drawnK$sd.stated, FALSE)
drawnSd <- lone(leaf.prior = normal(sd = invchi(3, 1)))$getLeafPrior()
expect_identical(drawnSd$sd, NA_real_)
expect_identical(drawnSd$sd.stated, TRUE)
sharer <- dbarts(
  y ~ forest(x1) + forest(x2) + forest(x1, basis = z),
  d,
  leaf.prior = normal(k = chi(1.5, 2)),
  control = ctl()
)$getLeafPrior(1L)
expect_identical(sharer$sd, drawnK$sd)
# the other leaf models and a binary response carry the pair too
linearReader <- lone(leaf.prior = linear("x1", sd = 0.5))$getLeafPrior()
expect_identical(linearReader$sd, 0.5)
expect_identical(linearReader$prior.sd.of, "coefficient")
binary <- dbarts(x, as.double(y > median(y)), control = ctl())$getLeafPrior()
expect_identical(binary$sd.stated, FALSE)

# --- an integer k is the number as a double ----------------------------------
fitK <- function(k, ...) {
  set.seed(7L)
  bart(
    x,
    y,
    k = k,
    n.trees = 10L,
    n.samples = 5L,
    n.burn = 3L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE,
    ...
  )
}
expect_identical(fitK(2L)$yhat.train, fitK(2)$yhat.train)
expect_identical(fitK(2L)$k, fitK(2)$k)
expect_identical(
  lone(leaf.prior = normal(k = 2L))$run()$train,
  lone(leaf.prior = normal(k = 2))$run()$train
)
expect_identical(normal(k = 2L)@k, 2)
expect_true(is.double(chi(1L, 2L)@degreesOfFreedom))
expect_true(is.double(chi(1L, 2L)@scale))
expect_identical(
  lone(leaf.prior = normal(k = chi(2L, 3L)))$getLeafPrior()$leaf.prior,
  lone(leaf.prior = normal(k = chi(2, 3)))$getLeafPrior()$leaf.prior
)
# a prior object built by hand holds the integer in its slot
expect_true(is(
  lone(leaf.prior = new("dbartsNormalPrior", k = 2L)),
  "dbartsSampler"
))
live <- lone()
expect_silent(live$setLeafPrior(normal(k = 3L)))
expect_identical(unname(live$getK()), 3)
expect_silent(live$setLeafPrior(normal(k = chi(3L, 2L))))
# the other doors that take a number
expect_true(is(
  dbarts(x, y, leaf.prior = normal(k = 2L), control = ctl()),
  "dbartsSampler"
))
expect_true(is(
  dbarts(
    x,
    as.double(y > median(y)),
    leaf.prior = normal(k = 2L),
    control = ctl()
  ),
  "dbartsSampler"
))
expect_true(is.list(fitK(2L, monotone = c(x1 = 1))))
expect_identical(
  xbart(
    x,
    y,
    k = c(2L, 3L),
    n.trees = 10L,
    n.samples = 5L,
    n.burn = 3L,
    n.reps = 2L,
    n.threads = 1L,
    verbose = FALSE,
    seed = 3L
  ),
  xbart(
    x,
    y,
    k = c(2, 3),
    n.trees = 10L,
    n.samples = 5L,
    n.burn = 3L,
    n.reps = 2L,
    n.threads = 1L,
    verbose = FALSE,
    seed = 3L
  )
)
# a k that is no positive number is still refused
for (bad in list(NA_integer_, 0L, -1L, TRUE)) {
  expect_error(
    normal(k = bad),
    "'k' must be a positive scalar or a hyperprior specification",
    fixed = TRUE
  )
}

# --- the creation summary names what the caller stated -----------------------
summaryLines <- function(leaf.prior) {
  out <- capture.output(
    dbarts(x, y, leaf.prior = leaf.prior, control = ctl(verbose = TRUE))
  )
  grep("^\tk prior|^\tsd prior|^\tprior on", out, value = TRUE)
}
expect_identical(summaryLines(normal(k = 3)), "\tk prior fixed to 3.000000")
expect_identical(
  summaryLines(normal(k = chi(1.5, 2))),
  "\tprior on k: chi with 1.500000 degrees of freedom and 2.000000 scale"
)
expect_identical(
  summaryLines(normal(sd = 0.5)),
  "\tsd prior fixed to 0.500000"
)
expect_identical(
  summaryLines(normal(sd = invchi(3, 1))),
  paste0(
    "\tprior on the sd: inverse chi with 3.000000 degrees of freedom and ",
    "1.000000 scale"
  )
)
expect_identical(
  summaryLines(normal(sd = invchi(3, 2.5))),
  paste0(
    "\tprior on the sd: inverse chi with 3.000000 degrees of freedom and ",
    "2.500000 scale"
  )
)

# --- a drawn size with no median: refused at creation, by forest -------------
# the engine takes the median of the chi variable from the chi-squared
# quantile at one half, and refuses a forest whose size is not positive and
# finite; df = 0.001 is a quantile of zero in double precision
expect_identical(qchisq(0.5, 0.001), 0)
studentForests <- y ~ forest(x1, amplitude = student()) + forest(x2, basis = z)
studentOne <- function(leaf.prior) {
  dbarts(studentForests, d, leaf.prior = leaf.prior, control = ctl())
}
sizeless <- paste0(
  "forest 1: amplitude = student() with no 'sd' takes its median from the ",
  "fitting function's leaf prior, which here gives no size; state the ",
  "forest's 'sd' or give the prior on k or the sd more degrees of freedom"
)
for (prior in list(
  normal(sd = invchi(0.001, 1)),
  normal(k = chi(0.001, 2))
)) {
  expect_error(studentOne(prior), sizeless, fixed = TRUE)
}
# a median that is a number, however small, is taken, as are the flat priors'
# own refusal and a stated sd
expect_true(is(studentOne(normal(sd = invchi(0.01, 1))), "dbartsSampler"))
expect_true(is(studentOne(normal(k = chi(0.01, 2))), "dbartsSampler"))
expect_error(
  studentOne(normal(k = chi(1.5, Inf))),
  "which here gives no size; state the forest's 'sd'",
  fixed = TRUE
)
expect_true(is(
  dbarts(
    y ~ forest(x1, sd = 1, amplitude = student()) + forest(x2, basis = z),
    d,
    control = ctl()
  ),
  "dbartsSampler"
))
# a forest taking the prior as a share of its total has a size, and is taken
expect_true(is(
  dbarts(
    y ~ forest(x1) + forest(x2, basis = z),
    d,
    leaf.prior = normal(sd = invchi(0.001, 1)),
    control = ctl()
  ),
  "dbartsSampler"
))

# --- $setLeafPrior on a Student-t forest with no sd ---------------------------
# the forest takes its size from the leaf prior at creation only
for (created in list(normal(), normal(k = 3))) {
  s <- studentOne(created)
  before <- s$getLeafPrior(1L)$sd
  expect_error(
    s$setLeafPrior(normal(k = 4)),
    "takes forests = ",
    fixed = TRUE
  )
  expect_silent(s$setLeafPrior(normal(k = 2)))
  expect_identical(s$getLeafPrior(1L)$sd, before)
}

# --- integers beyond k --------------------------------------------------------
# chi() and chisq() store doubles
expect_identical(chi(3L, 2L)@degreesOfFreedom, 3)
expect_identical(chi(3L, 2L)@scale, 2)
expect_identical(chisq(3L)@df, 3)
# an integer degrees of freedom, which bart()'s sigdf reaches, or offset draws
# as the same number as a double
drawsOf <- function(...) {
  dbarts(y ~ x1, d, control = ctl(), ...)$run()$train
}
expect_identical(
  drawsOf(family = gaussian(sigma = chisq(3L))),
  drawsOf(family = gaussian(sigma = chisq(3)))
)
expect_identical(
  drawsOf(offset = rep(1L, nrow(d))),
  drawsOf(offset = rep(1, nrow(d)))
)
expect_identical(
  dbartsData(y ~ x1, d, offset = rep(1L, nrow(d)))@offset,
  rep(1, nrow(d))
)

# --- the refusal's hint -------------------------------------------------------
# a flat prior gives no median that more degrees of freedom would mend
flatMessage <- tryCatch(
  studentOne(normal(k = chi(3, Inf))),
  error = conditionMessage
)
expect_true(is.character(flatMessage))
expect_false(grepl("more degrees of freedom", flatMessage, fixed = TRUE))

# --- a sharer's sd.stated -----------------------------------------------------
# on several forests it says whether the forest names an sd of its own
sharer <- dbarts(
  y ~ forest(x1) + forest(x2, basis = z),
  d,
  leaf.prior = normal(sd = 0.5),
  control = ctl()
)$getLeafPrior(1L)
expect_equivalent(sharer$sd, 0.5)
expect_identical(sharer$sd.stated, FALSE)
