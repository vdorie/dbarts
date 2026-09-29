# family = "auto": a count matrix is read as multinomial, every resolution is
# announced once per call, and an explicit family or verbose = FALSE is silent.

set.seed(3)
n <- 60L
x <- matrix(rnorm(n * 3L), n, 3L, dimnames = list(NULL, paste0("x", 1:3)))
counts <- t(rmultinom(n, 3L, c(0.2, 0.3, 0.5)))
colnames(counts) <- c("a", "b", "c")
frame <- data.frame(x, counts)
frame$cm <- counts

ctl <- dbartsControl(
  n.trees = 10L,
  n.samples = 15L,
  n.burn = 10L,
  n.chains = 1L,
  n.threads = 1L
)
quiet <- function(expr) {
  capture.output(value <- expr)
  value
}
# the "auto" announcements a call emitted, by count
autoLines <- function(expr) {
  lines <- character()
  capture.output(withCallingHandlers(
    expr,
    message = function(m) {
      if (grepl("family = \"auto\"", conditionMessage(m), fixed = TRUE)) {
        lines <<- c(lines, trimws(conditionMessage(m), "right"))
      }
      invokeRestart("muffleMessage")
    }
  ))
  lines
}

# --- a count matrix under auto == explicit multinomial, draws identical ---
set.seed(5)
fitExplicit <- quiet(bart(
  x,
  counts,
  family = "multinomial",
  control = ctl,
  verbose = FALSE
))
set.seed(5)
fitMatrix <- quiet(bart(x, counts, control = ctl, verbose = FALSE))
expect_identical(fitMatrix$yhat.train, fitExplicit$yhat.train)
expect_identical(fitMatrix$levels, c("a", "b", "c"))

set.seed(5)
fitCbind <- quiet(bart(
  cbind(a, b, c) ~ x1 + x2 + x3,
  frame,
  control = ctl,
  verbose = FALSE
))
set.seed(5)
fitCbindExplicit <- quiet(bart(
  cbind(a, b, c) ~ x1 + x2 + x3,
  frame,
  family = "multinomial",
  control = ctl,
  verbose = FALSE
))
expect_identical(fitCbind$yhat.train, fitCbindExplicit$yhat.train)

set.seed(5)
fitColumn <- quiet(bart(
  cm ~ x1 + x2 + x3,
  frame,
  control = ctl,
  verbose = FALSE
))
set.seed(5)
fitColumnExplicit <- quiet(bart(
  cm ~ x1 + x2 + x3,
  frame,
  family = "multinomial",
  control = ctl,
  verbose = FALSE
))
expect_identical(fitColumn$yhat.train, fitColumnExplicit$yhat.train)

# a row with a missing cell goes to na.action, as the explicit path does
countsNA <- counts
countsNA[2L, 1L] <- NA
set.seed(5)
fitNA <- quiet(bart(x, countsNA, control = ctl, verbose = FALSE))
set.seed(5)
fitNAExplicit <- quiet(bart(
  x,
  countsNA,
  family = "multinomial",
  control = ctl,
  verbose = FALSE
))
expect_identical(fitNA$yhat.train, fitNAExplicit$yhat.train)

# dbarts() takes the matrix interface; a formula is refused by name
sampler <- dbarts(x, counts, control = ctl)
expect_identical(sampler$model@family, "multinomial")
expect_error(
  dbarts(cbind(a, b, c) ~ x1 + x2, frame, control = ctl),
  "dbarts\\(\\) fits an n x 3 count matrix through the matrix interface only"
)

# --- refusals, each naming what to write ---
twoColumn <- "family = \"multinomial\".*survival::Surv.*family = \"aft\" / \"hazard\""
expect_error(bart(x, counts[, 1:2], control = ctl), twoColumn)
expect_error(bart(cbind(a, b) ~ x1, frame, control = ctl), twoColumn)
expect_error(dbarts(x, counts[, 1:2], control = ctl), twoColumn)
expect_error(
  bart(x, counts + 0.5, control = ctl),
  "does not read as multinomial counts"
)
expect_error(
  bart(x, counts - 1, control = ctl),
  "non-negative whole numbers"
)
expect_error(
  bart(cbind(a, b, c) ~ x1, transform(frame, a = a + 0.5), control = ctl),
  "non-negative whole numbers"
)
expect_error(
  dbarts(x, counts + 0.5, control = ctl),
  "non-negative whole numbers"
)
expect_error(
  xbart(x, counts, n.reps = 2L),
  "xbart\\(\\) takes a single-column response.*family = \"multinomial\""
)
expect_error(
  xbart(x, counts[, 1:2], n.reps = 2L),
  "xbart\\(\\) takes a single-column response.*survival::Surv"
)
# an explicit family that takes one column says so under its own name
expect_error(
  bart(x, counts, family = "gaussian", control = ctl),
  "family = \"gaussian\" takes a single-column response"
)

# --- the announcement: once per call, for every family auto resolves ---
yCont <- rnorm(n)
yBin <- rbinom(n, 1L, 0.5)
yFac3 <- factor(sample(c("p", "q", "r"), n, TRUE))
yOrd <- factor(sample(c("p", "q", "r"), n, TRUE), ordered = TRUE)

expected <- function(description, family) {
  paste0(
    "family = \"auto\": ",
    description,
    " detected, fitting family = \"",
    family,
    "\"; set 'family' to override"
  )
}
expect_identical(
  autoLines(bart(x, yCont, control = ctl)),
  expected("continuous response", "gaussian")
)
expect_identical(
  autoLines(bart(x, yBin, control = ctl)),
  expected("0/1 response", "probit")
)
expect_identical(
  autoLines(bart(x, factor(yBin), control = ctl)),
  expected("2-level factor response", "probit")
)
expect_identical(
  autoLines(bart(x, yFac3, control = ctl)),
  expected("3-level factor response", "multinomial")
)
expect_identical(
  autoLines(bart(x, yOrd, control = ctl)),
  expected("3-level ordered factor response", "ordinal")
)
expect_identical(
  autoLines(bart(x, counts, control = ctl)),
  expected("3-column count matrix response", "multinomial")
)
expect_identical(
  autoLines(bart(cbind(a, b, c) ~ x1, frame, control = ctl)),
  expected("3-column count matrix response", "multinomial")
)
expect_identical(
  autoLines(bart(cm ~ x1, frame, control = ctl)),
  expected("3-column count matrix response", "multinomial")
)
expect_identical(
  autoLines(bart(x, yCont, control = ctl, keepTrees = TRUE, n.chains = 2L)),
  expected("continuous response", "gaussian")
)
if (requireNamespace("survival", quietly = TRUE)) {
  survY <- survival::Surv(rexp(n), rbinom(n, 1L, 0.7))
  expect_identical(
    autoLines(bart(x, survY, control = ctl)),
    expected("survival (Surv) response", "aft")
  )
}

# dbarts() and dbartsSpec() speak only when their verbose is on
expect_identical(
  autoLines(dbarts(x, yCont, verbose = TRUE)),
  expected("continuous response", "gaussian")
)
expect_identical(autoLines(dbarts(x, yCont)), character())
expect_identical(
  autoLines(dbarts(x, counts, verbose = TRUE)),
  expected("3-column count matrix response", "multinomial")
)
spec.data <- dbartsData(x, yCont)
expect_identical(
  autoLines(dbartsSpec(spec.data, control = dbartsControl(verbose = TRUE))),
  expected("continuous response", "gaussian")
)
expect_identical(autoLines(dbartsSpec(spec.data)), character())

# xbart announces once, not per replication or worker
expect_identical(
  autoLines(xbart(x, yCont, n.reps = 3L, n.trees = 10L, verbose = TRUE)),
  expected("continuous response", "gaussian")
)
expect_identical(
  autoLines(xbart(
    x,
    yCont,
    n.reps = 3L,
    n.trees = 10L,
    n.threads = 2L,
    verbose = TRUE
  )),
  expected("continuous response", "gaussian")
)
expect_identical(
  autoLines(xbart(x, yBin, n.reps = 2L, n.trees = 10L, verbose = TRUE)),
  expected("0/1 response", "probit")
)
expect_identical(
  autoLines(xbart(x, yCont, n.reps = 2L, n.trees = 10L)),
  character()
)

# never under verbose = FALSE or an explicit family
expect_identical(
  autoLines(bart(x, yCont, control = ctl, verbose = FALSE)),
  character()
)
expect_identical(
  autoLines(bart(x, counts, control = ctl, verbose = FALSE)),
  character()
)
expect_identical(
  autoLines(bart(x, yCont, family = "gaussian", control = ctl)),
  character()
)
expect_identical(
  autoLines(bart(x, counts, family = "multinomial", control = ctl)),
  character()
)
expect_identical(
  autoLines(bart(x, yOrd, family = "ordinal", control = ctl)),
  character()
)

# expect_message with a pattern, once
expect_message(
  quiet(bart(x, yCont, control = ctl)),
  "fitting family = \"gaussian\""
)
expect_silent(quiet(bart(x, yCont, control = ctl, verbose = FALSE)))
