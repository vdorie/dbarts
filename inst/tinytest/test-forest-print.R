# The block a fit's print, a sampler's show and a verbose creation give a
# model of several forests: one line for each, its kind, the column's mean,
# sd and range as they stand over the rows in the likelihood, and the sd it
# resolved to, in the units of the header. Every figure is computed here from
# the fixture and printed to four significant digits.

forest <- dbartsForests$forest
fixed <- dbartsPriors$fixed

set.seed(20261021)
n <- 150L
x1 <- runif(n)
x2 <- runif(n)
z <- as.double(seq_len(n) %in% sample(n, 30L))
dose <- rnorm(n, 0.05, 0.02)
dose <- dose - min(dose)
age <- rnorm(n, 50, 17)
y <- 2 + sin(3 * x1) + z * x2 + 10 * dose + 0.01 * age + rnorm(n, sd = 0.5)
yb <- as.double(y > median(y))
w <- rep(1, n)
out <- seq(5L, n, length.out = 10L)
w[out] <- 0
d <- data.frame(y, yb, x1, x2, z, dose, age, w)
U <- sd(y)

# the numbers the block prints: four significant digits, trailing zeros off
g <- function(x) trimws(formatC(unname(x), digits = 4L, format = "g"))
header <- function(unit = U) {
  paste0("forests, sd in the response's units (its sd ", g(unit), "):")
}
numericLine <- function(i, label, column, perUnit, how, held = "", from = "") {
  mean <- zapsmall(c(mean(column), sd(column)))[1L]
  paste0(
    "  ",
    i,
    " ",
    label,
    ": numeric basis, mean ",
    g(mean),
    ", sd ",
    g(sd(column)),
    ", range ",
    g(min(column)),
    " to ",
    g(max(column)),
    "; sd ",
    g(perUnit),
    " per unit, ",
    g(perUnit * sd(column)),
    " per sd (",
    how,
    from,
    ")",
    held
  )
}

printControl <- function(verbose = FALSE) {
  dbartsControl(
    n.chains = 2L,
    n.threads = 1L,
    n.samples = 10L,
    n.burn = 0L,
    updateState = FALSE,
    seed = 71L,
    verbose = verbose
  )
}
fourForests <- function(verbose = FALSE, ...) {
  dbarts(
    y ~ forest(x1 + x2) +
      forest(x1, basis = dose) +
      forest(x1, basis = age) +
      forest(x1, basis = factor(z), amplitude = fixed()),
    d,
    control = printControl(verbose),
    ...
  )
}
expected <- c(
  header(),
  paste0("  1 forest1: no basis; sd ", g(2 * U), " (default)"),
  numericLine(2, "dose", dose, U / sd(dose), "default"),
  numericLine(3, "age", age, U / sd(age), "default"),
  paste0(
    "  4 factor(z): factor basis, 2 levels; sd ",
    g(U),
    " between the two levels (default); coefficient held"
  )
)
# the block out of a printed object: from its header on, NULL where it has none
blockOf <- function(lines) {
  start <- which(startsWith(lines, "forests, sd "))
  if (length(start) != 1L) {
    return(NULL)
  }
  lines[start:max(which(nzchar(lines)))]
}
numWarnings <- 0L
counted <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    numWarnings <<- numWarnings + 1L
    invokeRestart("muffleWarning")
  })
}

# --- the model of four, at show, at creation and at print ---
sampler <- counted(fourForests())
shown <- capture.output(show(sampler))
expect_identical(blockOf(shown), expected)
# the call comes first and a blank line closes it
expect_identical(shown[1L], "dbarts sampler")
expect_identical(shown[length(shown)], "")

created <- capture.output(invisible(counted(fourForests(verbose = TRUE))))
expect_identical(blockOf(created), expected)

fit <- counted(bart(
  y ~ forest(x1 + x2) +
    forest(x1, basis = dose) +
    forest(x1, basis = age) +
    forest(x1, basis = factor(z), amplitude = fixed()),
  d,
  n.samples = 10L,
  n.burn = 5L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 17L,
  verbose = FALSE
))
expect_identical(blockOf(capture.output(print(fit))), expected)
# a kept sampler prints the same
kept <- counted(bart(
  y ~ forest(x1 + x2) +
    forest(x1, basis = dose) +
    forest(x1, basis = age) +
    forest(x1, basis = factor(z), amplitude = fixed()),
  d,
  n.samples = 10L,
  n.burn = 5L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 17L,
  verbose = FALSE,
  keepSampler = TRUE
))
printedKept <- capture.output(print(kept))
expect_identical(blockOf(printedKept), expected)
expect_true(any(startsWith(printedKept, "n.trees: ")))
# a fit with no unit on it, one made before it was kept, prints the block
# without the number in its header
bare <- fit
attr(bare, "forest.unit") <- NULL
expect_identical(
  blockOf(capture.output(print(bare))),
  c("forests, sd in the response's units:", expected[-1L])
)

# a copy is not a creation, and a creation prints the block once
verboseSampler <- NULL
copied <- capture.output(verboseSampler <- fourForests(verbose = TRUE))
expect_identical(sum(startsWith(copied, "forests, sd ")), 1L)
copyOut <- capture.output(invisible(verboseSampler$copy()))
expect_identical(sum(startsWith(copyOut, "forests, sd ")), 0L)

# --- every number in the numeric lines is the reader's or R's ---
numbersOf <- function(line) {
  as.numeric(regmatches(
    line,
    gregexpr("-?[0-9]+\\.?[0-9]*(e[-+][0-9]+)?", line)
  )[[1L]])
}
reader <- sampler$getLeafPrior()
for (i in 2:3) {
  column <- if (i == 2L) dose else age
  # index, mean, sd, range, sd per unit, per sd
  read <- numbersOf(shown[which(startsWith(shown, "forests, sd ")) + i])
  truth <- c(
    i,
    mean(column),
    sd(column),
    range(column),
    unname(reader[[i]]$sd),
    unname(reader[[i]]$sd) * sd(column)
  )
  expect_equal(read, signif(truth, 4L), tolerance = 1e-12)
}

# --- a stated sd, a held coefficient, and a forest with no basis held ---
stated <- counted(dbarts(
  y ~ forest(x1 + x2, amplitude = fixed()) +
    forest(x1, basis = dose, sd = 30, amplitude = fixed()) +
    forest(x1, basis = factor(z), sd = 0.5),
  d,
  control = printControl()
))
statedLines <- blockOf(capture.output(show(stated)))
expect_identical(statedLines[1L], header())
heldSd <- stated$getLeafPrior(1L)$sd
expect_identical(
  statedLines[2L],
  paste0(
    "  1 forest1: no basis; sd ",
    g(heldSd),
    " (default); coefficient held"
  )
)
expect_identical(
  statedLines[3L],
  numericLine(2, "dose", dose, 30, "stated", "; coefficient held")
)
expect_identical(
  statedLines[4L],
  "  3 factor(z): factor basis, 2 levels; sd 0.5 between the two levels (stated)"
)

# --- probit and logistic: the header names the latent scale ---
for (family in c("probit", "logistic")) {
  latent <- counted(dbarts(
    yb ~ forest(x1 + x2) + forest(x1, basis = dose),
    d,
    family = family,
    control = printControl()
  ))
  unit <- if (family == "probit") 1 else pi / sqrt(3)
  latentLines <- blockOf(capture.output(show(latent)))
  expect_identical(
    latentLines,
    c(
      paste0("forests, sd on the latent scale (error sd ", g(unit), "):"),
      paste0("  1 forest1: no basis; sd ", g(unit), " (default)"),
      numericLine(2, "dose", dose, unit / sd(dose), "default")
    ),
    info = family
  )
}
expect_identical(g(pi / sqrt(3)), "1.814")

# --- after a swap, the column is the one in force and the line ends with the
# recorded one ---
swapped <- fourForests()
swapped$setForestBasis(2L, dose / 10)
swappedLines <- blockOf(capture.output(show(swapped)))
expect_identical(
  swappedLines[3L],
  numericLine(
    2,
    "dose",
    dose / 10,
    U / sd(dose),
    "default",
    from = paste0(", taken from a column of sd ", g(sd(dose)))
  )
)
expect_identical(swappedLines[-3L], expected[-3L])
# the reader's number is the record's, so the sd per unit does not move
expect_equal(
  unname(swapped$getLeafPrior(2L)$sd),
  U / sd(dose),
  tolerance = 1e-12
)

# --- the rows in the likelihood: a row of weight zero is not in the column ---
heavy <- d
heavy$dose[out] <- 1000
weighted <- counted(dbarts(
  y ~ forest(x1 + x2) + forest(x1, basis = dose),
  heavy,
  weights = w,
  control = printControl()
))
weightedLines <- blockOf(capture.output(show(weighted)))
rowsIn <- w > 0
unitIn <- sd(y[rowsIn])
perUnitIn <- unitIn / sd(dose[rowsIn])
expect_identical(
  weightedLines[3L],
  numericLine(2, "dose", dose[rowsIn], perUnitIn, "default")
)
expect_identical(weightedLines[1L], header(unitIn))
weightedFit <- counted(bart(
  y ~ forest(x1 + x2) + forest(x1, basis = dose),
  heavy,
  weights = w,
  n.samples = 5L,
  n.burn = 0L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 3L,
  verbose = FALSE
))
expect_identical(blockOf(capture.output(print(weightedFit))), weightedLines)
# the weights moved afterwards: the column is read over the new rows, the
# default stays the creation's
weighted$setWeights(rep(1, n))
reweighted <- blockOf(capture.output(show(weighted)))
expect_identical(
  reweighted[3L],
  numericLine(
    2,
    "dose",
    heavy$dose,
    perUnitIn,
    "default",
    from = paste0(", taken from a column of sd ", g(sd(dose[rowsIn])))
  )
)

# --- a centered column prints mean 0, one wild value shows in the range ---
centered <- counted(dbarts(
  y ~ forest(x1 + x2) + forest(x1, basis = scale(age)),
  d,
  control = printControl()
))
centeredLine <- blockOf(capture.output(show(centered)))[3L]
expect_true(grepl(
  ": numeric basis, mean 0, sd 1, range ",
  centeredLine,
  fixed = TRUE
))
expect_false(grepl("e-", centeredLine, fixed = TRUE))
wild <- d$dose
wild[7L] <- 1000 * max(dose)
wildSampler <- counted(dbarts(
  y ~ forest(x1 + x2) + forest(x1, basis = wild),
  data.frame(d, wild),
  control = printControl()
))
wildLine <- blockOf(capture.output(show(wildSampler)))[3L]
expect_true(grepl(
  paste0("range ", g(min(wild)), " to ", g(max(wild)), ";"),
  wildLine,
  fixed = TRUE
))
expect_equal(max(wild), 1000 * max(dose))

# --- a constant column, which has a stated sd, prints no figure per sd ---
flat <- rep(2, n)
flatSampler <- counted(dbarts(
  y ~ forest(x1 + x2) + forest(x1, basis = flat, sd = 5),
  data.frame(d, flat),
  control = printControl()
))
expect_identical(
  blockOf(capture.output(show(flatSampler)))[3L],
  "  2 flat: numeric basis, mean 2, sd 0, range 2 to 2; sd 5 per unit (stated)"
)

# --- a 0/1 number and an age at their defaults ---
zero <- counted(dbarts(
  y ~ forest(x1 + x2) + forest(x1, basis = z) + forest(x1, basis = age),
  d,
  control = printControl(verbose = TRUE)
))
capture.output(show(zero))
capture.output(print(counted(bart(
  y ~ forest(x1 + x2) + forest(x1, basis = z) + forest(x1, basis = age),
  d,
  n.samples = 5L,
  n.burn = 0L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 5L,
  verbose = FALSE
))))
zeroOne <- blockOf(capture.output(show(zero)))[3L]
expect_identical(
  zeroOne,
  numericLine(2, "z", z, U / sd(z), "default")
)

# --- one forest: nothing is printed ---
single <- dbarts(y ~ x1 + x2, d, control = printControl(verbose = TRUE))
expect_null(blockOf(capture.output(show(single))))
expect_null(blockOf(capture.output(
  dbarts(y ~ x1 + x2, d, control = printControl(verbose = TRUE))
)))
oneFit <- bart(
  y ~ x1 + x2,
  d,
  n.samples = 5L,
  n.burn = 0L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 5L,
  verbose = FALSE
)
expect_null(blockOf(capture.output(print(oneFit))))

# --- no warning anywhere above, an age and a 0/1 number included ---
expect_identical(numWarnings, 0L)
