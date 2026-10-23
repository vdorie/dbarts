# The block of a model of several forests. print of a fit and show of a
# sampler give the short one: a line for each forest, its label, its kind, a
# numeric column's mean and sd as they stand over the rows in the likelihood,
# and the sd it resolved to. summary of a fit and a verbose creation give the
# full one: under a header that names the unit, the column's range, the sd per
# sd of the column, where the sd came from and whether the coefficient is
# held. Every figure is computed here from the fixture and printed as
# format(x, digits = 4) gives it.

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

# the numbers the block prints: format(x, digits = 4), so at least four
# significant digits and every digit of an integer part
g <- function(x) format(unname(x), digits = 4L)
header <- function(unit = U) {
  paste0("forests, sd in the response's units (sd(y) ", g(unit), "):")
}
latentHeader <- function(errorSd) {
  paste0("forests, sd on the latent scale (error sd ", g(errorSd), "):")
}
shortHeader <- "forests:"
fitFooter <- "(each forest's prior in full: summary())"
samplerFooter <- "(each forest's prior in full: $getLeafPrior())"
# a numeric forest's line as far as the column goes: its mean and sd
numericStem <- function(i, label, column) {
  mean <- zapsmall(c(mean(column), sd(column)))[1L]
  paste0(
    "  ",
    i,
    " ",
    label,
    ": numeric basis, mean ",
    g(mean),
    ", sd ",
    g(sd(column))
  )
}
numericShort <- function(i, label, column, perUnit) {
  paste0(numericStem(i, label, column), "; sd ", g(perUnit), " per unit")
}
numericLine <- function(i, label, column, perUnit, how, held = "", from = "") {
  paste0(
    numericStem(i, label, column),
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
expectedFull <- c(
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
shortForests <- c(
  paste0("  1 forest1: no basis; sd ", g(2 * U)),
  numericShort(2, "dose", dose, U / sd(dose)),
  numericShort(3, "age", age, U / sd(age)),
  paste0(
    "  4 factor(z): factor basis, 2 levels; sd ",
    g(U),
    " between the two levels"
  )
)
expectedFit <- c(shortHeader, shortForests, fitFooter)
expectedSampler <- c(shortHeader, shortForests, samplerFooter)
# a block out of a printed object: from its header on, NULL where it has none
blockFrom <- function(lines, prefix) {
  start <- which(startsWith(lines, prefix))
  if (length(start) != 1L) {
    return(NULL)
  }
  lines[start:max(which(nzchar(lines)))]
}
starts <- function(line, prefix) {
  is.character(line) && length(line) == 1L && startsWith(line, prefix)
}
fullOf <- function(lines) blockFrom(lines, "forests, sd ")
shortOf <- function(lines) blockFrom(lines, "forests:")
numWarnings <- 0L
counted <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    numWarnings <<- numWarnings + 1L
    invokeRestart("muffleWarning")
  })
}
# a creation whose warning the help promises: kept apart, and read where it is
# raised
creationWarnings <- character(0L)
collected <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    creationWarnings <<- c(creationWarnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
}
takeCreationWarnings <- function() {
  taken <- creationWarnings
  creationWarnings <<- character(0L)
  taken
}
# the lines a show, a print or a printed summary gives, with the warnings of
# the call counted
showLines <- function(x) counted(capture.output(show(x)))
printLines <- function(x) counted(capture.output(print(x)))
summaryLines <- function(x) counted(capture.output(print(summary(x))))
# a creation with verbose = TRUE gives the engine's summary on the output and
# its resolutions as messages: both are taken off the run, and the output is
# what it returns
verboseLines <- function(expr) {
  capture.output(invisible(suppressMessages(counted(expr))))
}
# a sampler created with verbose = TRUE, and the full block its creation gave
built <- function(expr) {
  made <- NULL
  printed <- verboseLines(made <- expr)
  list(sampler = made, full = fullOf(printed))
}
# the full block of a sampler as it stands: only a creation prints it, so what
# follows a change to the sampler asks the function that builds it
fullNow <- function(sampler) {
  counted(dbarts:::samplerForestLines(sampler, full = TRUE))
}
fitArgs <- list(
  n.samples = 5L,
  n.burn = 0L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 3L,
  verbose = FALSE
)
bartFit <- function(formula, data, ...) {
  do.call(bart, c(list(formula, data, ...), fitArgs))
}

# --- the model of four, at show, at creation, at print and at summary ---
sampler <- counted(fourForests())
shown <- showLines(sampler)
expect_identical(shortOf(shown), expectedSampler)
expect_null(fullOf(shown))
# the call comes first, as it was written with its forest() terms, and a blank
# line closes the block
expect_identical(shown[1L], "dbarts sampler")
expect_true(grepl(
  "call: dbarts(formula = y ~ forest(x1 + x2)",
  shown[2L],
  fixed = TRUE
))
expect_identical(
  sampler$control@call$formula,
  quote(
    y ~ forest(x1 + x2) +
      forest(x1, basis = dose) +
      forest(x1, basis = age) +
      forest(x1, basis = factor(z), amplitude = fixed())
  )
)
expect_identical(shown[length(shown)], "")
# a formula passed by name is recorded by name
fm <- y ~ forest(x1 + x2) + forest(x1, basis = dose)
byName <- counted(dbarts(fm, d, control = printControl()))
expect_identical(byName$control@call$formula, quote(fm))

# a verbose creation gives the full block, and only that
created <- verboseLines(fourForests(verbose = TRUE))
expect_identical(fullOf(created), expectedFull)
expect_null(shortOf(created))
expect_identical(fullNow(sampler), expectedFull)

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
printedFit <- printLines(fit)
expect_identical(shortOf(printedFit), expectedFit)
expect_null(fullOf(printedFit))
expect_equal(attr(fit, "forest.unit"), U, tolerance = 1e-12)
# summary gives the full block after its table, and the line under the table
# no longer names the sds the block gives
summed <- summaryLines(fit)
expect_identical(fullOf(summed), expectedFull)
expect_null(shortOf(summed))
expect_false(any(grepl("leaf.prior.sd", summed, fixed = TRUE)))
expect_false(any(grepl("Fixed, not sampled", summed, fixed = TRUE)))
expect_identical(summary(fit)$forests, expectedFull)
expect_identical(names(summary(fit)$fixed), "leaf.prior.sd")
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
printedKept <- printLines(kept)
expect_identical(shortOf(printedKept), expectedFit)
expect_true(any(startsWith(printedKept, "n.trees: ")))
# what else a fit held fixed stays on the line under the table
sigmaHeld <- counted(bartFit(
  y ~ forest(x1 + x2) + forest(x1, basis = dose),
  d,
  family = dbartsFamilies$gaussian(sigma = fixed(0.25))
))
heldLines <- summaryLines(sigmaHeld)
expect_true("(Fixed, not sampled: sigma = 0.5)" %in% heldLines)
expect_false(any(grepl("leaf.prior.sd[", heldLines, fixed = TRUE)))
expect_identical(length(fullOf(heldLines)), 3L)
# a fit made without keeping its fits keeps no bases and no labels, so it
# prints no block, short or full, and its summary keeps the sds on the line
# under the table
noFits <- counted(bartFit(
  y ~ forest(x1 + x2) + forest(x1, basis = dose),
  d,
  keepFits = FALSE
))
expect_null(noFits[["bases"]])
expect_null(attr(noFits, "forest.labels", exact = TRUE))
expect_null(attr(noFits, "forest.unit", exact = TRUE))
expect_null(shortOf(printLines(noFits)))
noFitsSummary <- summaryLines(noFits)
expect_null(fullOf(noFitsSummary))
expect_null(summary(noFits)$forests)
expect_true(any(grepl(
  "(Fixed, not sampled: leaf.prior.sd = ",
  noFitsSummary,
  fixed = TRUE
)))

# a copy is not a creation, and a creation prints the full block once
verboseSampler <- NULL
copied <- verboseLines(verboseSampler <- fourForests(verbose = TRUE))
expect_identical(sum(startsWith(copied, "forests, sd ")), 1L)
copyOut <- verboseLines(verboseSampler$copy())
expect_identical(sum(startsWith(copyOut, "forests")), 0L)
# the hold is let go on the way out: the source still shows its block, and the
# next creation prints its own
expect_identical(shortOf(showLines(verboseSampler)), expectedSampler)
afterCopy <- verboseLines(fourForests(verbose = TRUE))
expect_identical(sum(startsWith(afterCopy, "forests, sd ")), 1L)
expect_identical(fullOf(afterCopy), expectedFull)

# --- every number in the numeric lines is the reader's or R's ---
numbersOf <- function(line) {
  if (!is.character(line) || length(line) != 1L) {
    return(numeric(0L))
  }
  as.numeric(regmatches(
    line,
    gregexpr("-?[0-9]+\\.?[0-9]*(e[-+][0-9]+)?", line)
  )[[1L]])
}
reader <- sampler$getLeafPrior()
for (i in 2:3) {
  column <- if (i == 2L) dose else age
  # the full line: index, mean, sd, range, sd per unit, per sd
  read <- numbersOf(fullOf(created)[i + 1L])
  truth <- c(
    i,
    mean(column),
    sd(column),
    range(column),
    unname(reader[[i]]$sd),
    unname(reader[[i]]$sd) * sd(column)
  )
  expect_equal(read, signif(truth, 4L), tolerance = 1e-12)
  # the short line: index, mean, sd, sd per unit
  readShort <- numbersOf(shortOf(shown)[i + 1L])
  expect_equal(
    readShort,
    signif(c(i, mean(column), sd(column), unname(reader[[i]]$sd)), 4L),
    tolerance = 1e-12
  )
}

# --- a stated sd, a held coefficient, and a forest with no basis held ---
statedMade <- built(dbarts(
  y ~ forest(x1 + x2, amplitude = fixed()) +
    forest(x1, basis = dose, sd = 30, amplitude = fixed()) +
    forest(x1, basis = factor(z), sd = 0.5),
  d,
  control = printControl(verbose = TRUE)
))
stated <- statedMade$sampler
statedLines <- statedMade$full
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
# the short block says neither where an sd came from nor what is held
expect_identical(
  shortOf(showLines(stated)),
  c(
    shortHeader,
    paste0("  1 forest1: no basis; sd ", g(heldSd)),
    numericShort(2, "dose", dose, 30),
    "  3 factor(z): factor basis, 2 levels; sd 0.5 between the two levels",
    samplerFooter
  )
)

# --- probit and logistic: the full header names the latent scale, at a
# creation and at the summary of a fit; the short block names no unit ---
for (family in c("probit", "logistic")) {
  latentMade <- built(dbarts(
    yb ~ forest(x1 + x2) + forest(x1, basis = dose),
    d,
    family = family,
    control = printControl(verbose = TRUE)
  ))
  errorSd <- if (family == "probit") 1 else pi / sqrt(3)
  latentFull <- c(
    latentHeader(errorSd),
    paste0("  1 forest1: no basis; sd ", g(errorSd), " (default)"),
    numericLine(2, "dose", dose, errorSd / sd(dose), "default")
  )
  latentShort <- c(
    shortHeader,
    paste0("  1 forest1: no basis; sd ", g(errorSd)),
    numericShort(2, "dose", dose, errorSd / sd(dose))
  )
  expect_identical(latentMade$full, latentFull, info = family)
  expect_identical(
    shortOf(showLines(latentMade$sampler)),
    c(latentShort, samplerFooter),
    info = family
  )
  latentFit <- counted(bartFit(
    yb ~ forest(x1 + x2) + forest(x1, basis = dose),
    d,
    family = family
  ))
  expect_identical(
    shortOf(printLines(latentFit)),
    c(latentShort, fitFooter),
    info = family
  )
  expect_identical(
    fullOf(summaryLines(latentFit)),
    latentFull,
    info = family
  )
}
expect_identical(g(pi / sqrt(3)), "1.814")
expect_identical(
  latentHeader(pi / sqrt(3)),
  "forests, sd on the latent scale (error sd 1.814):"
)

# --- after a swap, the column is the one in force and the full line ends with
# the recorded one; the short line has only the column ---
swapped <- counted(fourForests())
swapped$setForestBasis(2L, dose / 10)
swappedFull <- fullNow(swapped)
expect_identical(
  swappedFull[3L],
  numericLine(
    2,
    "dose",
    dose / 10,
    U / sd(dose),
    "default",
    from = paste0(", taken from a column of sd ", g(sd(dose)))
  )
)
expect_identical(swappedFull[-3L], expectedFull[-3L])
swappedShort <- shortOf(showLines(swapped))
expect_identical(
  swappedShort[3L],
  numericShort(2, "dose", dose / 10, U / sd(dose))
)
expect_identical(swappedShort[-3L], expectedSampler[-3L])
# the reader's number is the record's, so the sd per unit does not move
expect_equal(
  unname(swapped$getLeafPrior(2L)$sd),
  U / sd(dose),
  tolerance = 1e-12
)
# a stated line never takes the suffix, whatever column it was stated on
swapped$setLeafPrior(forests = list(forest(), forest(sd = 30)))
expect_identical(
  fullNow(swapped)[3L],
  numericLine(2, "dose", dose / 10, 30, "stated")
)

# --- the rows in the likelihood: a row of weight zero is not in the column ---
heavy <- d
heavy$dose[out] <- 1000
weightedMade <- built(dbarts(
  y ~ forest(x1 + x2) + forest(x1, basis = dose),
  heavy,
  weights = w,
  control = printControl(verbose = TRUE)
))
weighted <- weightedMade$sampler
weightedLines <- weightedMade$full
rowsIn <- w > 0
unitIn <- sd(y[rowsIn])
perUnitIn <- unitIn / sd(dose[rowsIn])
expect_identical(
  weightedLines[3L],
  numericLine(2, "dose", dose[rowsIn], perUnitIn, "default")
)
expect_identical(weightedLines[1L], header(unitIn))
weightedShort <- c(
  shortHeader,
  paste0("  1 forest1: no basis; sd ", g(2 * unitIn)),
  numericShort(2, "dose", dose[rowsIn], perUnitIn)
)
expect_identical(
  shortOf(showLines(weighted)),
  c(weightedShort, samplerFooter)
)
weightedFit <- counted(bartFit(
  y ~ forest(x1 + x2) + forest(x1, basis = dose),
  heavy,
  weights = w
))
expect_identical(
  shortOf(printLines(weightedFit)),
  c(weightedShort, fitFooter)
)
expect_identical(fullOf(summaryLines(weightedFit)), weightedLines)
# the weights moved afterwards: the column is read over the new rows, the
# default stays the creation's
weighted$setWeights(rep(1, n))
expect_identical(
  fullNow(weighted)[3L],
  numericLine(
    2,
    "dose",
    heavy$dose,
    perUnitIn,
    "default",
    from = paste0(", taken from a column of sd ", g(sd(dose[rowsIn])))
  )
)
expect_identical(
  shortOf(showLines(weighted))[3L],
  numericShort(2, "dose", heavy$dose, perUnitIn)
)

# --- a centered column prints mean 0, one wild value shows in the range ---
centeredMade <- built(dbarts(
  y ~ forest(x1 + x2) + forest(x1, basis = scale(age)),
  d,
  control = printControl(verbose = TRUE)
))
centeredLine <- centeredMade$full[3L]
expect_true(grepl(
  ": numeric basis, mean 0, sd 1, range ",
  centeredLine,
  fixed = TRUE
))
expect_false(grepl("e-", centeredLine, fixed = TRUE))
centeredShort <- shortOf(showLines(centeredMade$sampler))[3L]
expect_true(starts(
  centeredShort,
  "  2 scale(age): numeric basis, mean 0, sd 1; sd "
))
expect_false(grepl("e-", centeredShort, fixed = TRUE))
wild <- d$dose
wild[7L] <- 1000 * max(dose)
wildMade <- built(dbarts(
  y ~ forest(x1 + x2) + forest(x1, basis = wild),
  data.frame(d, wild),
  control = printControl(verbose = TRUE)
))
expect_true(grepl(
  paste0("range ", g(min(wild)), " to ", g(max(wild)), ";"),
  wildMade$full[3L],
  fixed = TRUE
))
expect_equal(max(wild), 1000 * max(dose))
# the short block has no range; the wild value moves the mean and sd it has
wildShort <- shortOf(showLines(wildMade$sampler))[3L]
expect_false(grepl("range", wildShort, fixed = TRUE))
expect_identical(wildShort, numericShort(2, "wild", wild, U / sd(wild)))

# --- a constant column, which has a stated sd, prints no figure per sd ---
flat <- rep(2, n)
flatMade <- built(dbarts(
  y ~ forest(x1 + x2) + forest(x1, basis = flat, sd = 5),
  data.frame(d, flat),
  control = printControl(verbose = TRUE)
))
expect_identical(
  flatMade$full[3L],
  "  2 flat: numeric basis, mean 2, sd 0, range 2 to 2; sd 5 per unit (stated)"
)
expect_identical(
  shortOf(showLines(flatMade$sampler))[3L],
  "  2 flat: numeric basis, mean 2, sd 0; sd 5 per unit"
)

# --- a column that is constant to rounding, which has a stated sd, is
# constant here too, where an exact zero is not what is asked ---
nearly <- 61.7 + 1e-13 * sin(seq_len(n))
expect_true(sd(nearly) > 0)
nearlyMade <- built(dbarts(
  y ~ forest(x1 + x2) + forest(x1, basis = nearly, sd = 5),
  data.frame(d, nearly),
  control = printControl(verbose = TRUE)
))
expect_identical(
  nearlyMade$full[3L],
  paste0(
    "  2 nearly: numeric basis, mean 61.7, sd ",
    g(sd(nearly)),
    ", range 61.7 to 61.7; sd 5 per unit (stated)"
  )
)
expect_identical(
  shortOf(showLines(nearlyMade$sampler))[3L],
  paste0(
    "  2 nearly: numeric basis, mean 61.7, sd ",
    g(sd(nearly)),
    "; sd 5 per unit"
  )
)

# --- a number of 1e4 and up prints every digit of its integer part ---
big <- 1000 * age
bigMade <- built(dbarts(
  y ~ forest(x1 + x2) + forest(x1, basis = big),
  data.frame(d, big),
  control = printControl(verbose = TRUE)
))
bigLine <- bigMade$full[3L]
expect_identical(bigLine, numericLine(2, "big", big, U / sd(big), "default"))
expect_true(grepl(
  "^  2 big: numeric basis, mean [0-9]{5}, sd [0-9]{5}, range [0-9]{5} to [0-9]{5}; sd ",
  bigLine
))
expect_false(grepl("e+", bigLine, fixed = TRUE))
bigShort <- shortOf(showLines(bigMade$sampler))[3L]
expect_identical(bigShort, numericShort(2, "big", big, U / sd(big)))
expect_true(grepl(
  "^  2 big: numeric basis, mean [0-9]{5}, sd [0-9]{5}; sd ",
  bigShort
))

# --- ten forests: the index is padded to the width of the last ---
tenForests <- c(list(forest()), rep(list(forest(basis = dose)), 9L))
names(tenForests) <- paste0("f", 1:10)
tenMade <- built(dbarts(
  y ~ x1 + x2,
  d,
  forests = tenForests,
  control = printControl(verbose = TRUE)
))
tenShort <- shortOf(showLines(tenMade$sampler))
expect_identical(length(tenMade$full), 11L)
expect_identical(length(tenShort), 12L)
for (tenLines in list(tenMade$full, tenShort)) {
  expect_true(starts(tenLines[2L], "   1 f1: no basis; sd "))
  expect_true(starts(tenLines[3L], "   2 f2: numeric basis, "))
  expect_true(starts(tenLines[10L], "   9 f9: numeric basis, "))
  expect_true(starts(tenLines[11L], "  10 f10: numeric basis, "))
}

# --- the full header's unit is sd(y) as the help takes it ---
# net of the offset
off <- 3 * sin(seq_len(n))
yOff <- y + off
dOff <- data.frame(d, yOff)
expect_true(g(sd(yOff)) != g(sd(yOff - off)))
offsetFull <- c(
  header(sd(yOff - off)),
  paste0("  1 forest1: no basis; sd ", g(2 * sd(yOff - off)), " (default)"),
  paste0(
    "  2 factor(z): factor basis, 2 levels; sd ",
    g(sd(yOff - off)),
    " between the two levels (default)"
  )
)
offsetMade <- built(dbarts(
  yOff ~ forest(x1 + x2) + forest(x1, basis = factor(z)),
  dOff,
  offset = off,
  control = printControl(verbose = TRUE)
))
expect_identical(offsetMade$full, offsetFull)
offsetFit <- counted(bartFit(
  yOff ~ forest(x1 + x2) + forest(x1, basis = factor(z)),
  dOff,
  offset = off
))
expect_identical(fullOf(summaryLines(offsetFit)), offsetFull)
# a response that does not vary: a quarter of the window of width one
dConstant <- data.frame(d, yConstant = 3)
constantFull <- c(
  "forests, sd in the response's units (sd(y) 0.25):",
  "  1 forest1: no basis; sd 0.5 (default)",
  "  2 factor(z): factor basis, 2 levels; sd 0.25 between the two levels (default)"
)
constantSampler <- collected(dbarts(
  yConstant ~ forest(x1 + x2) + forest(x1, basis = factor(z)),
  dConstant,
  control = printControl()
))
constantFit <- collected(bartFit(
  yConstant ~ forest(x1 + x2) + forest(x1, basis = factor(z)),
  dConstant
))
constantWarnings <- takeCreationWarnings()
expect_identical(length(constantWarnings), 2L)
expect_true(all(grepl("indistinguishable", constantWarnings, fixed = TRUE)))
expect_identical(fullNow(constantSampler), constantFull)
expect_identical(fullOf(summaryLines(constantFit)), constantFull)
expect_identical(
  shortOf(printLines(constantFit)),
  c(
    shortHeader,
    "  1 forest1: no basis; sd 0.5",
    "  2 factor(z): factor basis, 2 levels; sd 0.25 between the two levels",
    fitFooter
  )
)
# a response replaced afterwards: the unit stays the creation's, as the
# defaults do
replaced <- counted(dbarts(
  y ~ forest(x1 + x2) + forest(x1, basis = dose),
  d,
  control = printControl()
))
beforeReplace <- fullNow(replaced)
beforeReplaceShort <- shortOf(showLines(replaced))
replaced$setResponse(10 * y)
expect_identical(beforeReplace[1L], header())
expect_identical(fullNow(replaced), beforeReplace)
expect_identical(shortOf(showLines(replaced)), beforeReplaceShort)

# --- under probit, weights of 0 and 1 are the mask: a row of weight zero is not
# in the column, on a sampler and on a fit ---
maskedFull <- c(
  latentHeader(1),
  "  1 forest1: no basis; sd 1 (default)",
  numericLine(2, "dose", dose[rowsIn], 1 / sd(dose[rowsIn]), "default")
)
maskedShort <- c(
  shortHeader,
  "  1 forest1: no basis; sd 1",
  numericShort(2, "dose", dose[rowsIn], 1 / sd(dose[rowsIn]))
)
maskedMade <- built(dbarts(
  yb ~ forest(x1 + x2) + forest(x1, basis = dose),
  heavy,
  weights = w,
  family = "probit",
  control = printControl(verbose = TRUE)
))
expect_identical(sum(maskedMade$sampler$activeRows == 0), length(out))
expect_identical(maskedMade$full, maskedFull)
expect_identical(
  shortOf(showLines(maskedMade$sampler)),
  c(maskedShort, samplerFooter)
)
maskedFit <- counted(bartFit(
  yb ~ forest(x1 + x2) + forest(x1, basis = dose),
  heavy,
  weights = w,
  family = "probit"
))
expect_null(maskedFit[["weights"]])
expect_identical(sum(maskedFit[["active"]] == 0), length(out))
expect_identical(shortOf(printLines(maskedFit)), c(maskedShort, fitFooter))
expect_identical(fullOf(summaryLines(maskedFit)), maskedFull)

# --- no row in the likelihood: the column is read over every row, under the
# creation's one warning and none at show, print or summary ---
fallbackFull <- c(
  latentHeader(1),
  "  1 forest1: no basis; sd 1 (default)",
  numericLine(2, "dose", dose, 1 / sd(dose), "default")
)
fallbackShort <- c(
  shortHeader,
  "  1 forest1: no basis; sd 1",
  numericShort(2, "dose", dose, 1 / sd(dose))
)
fallbackSampler <- collected(dbarts(
  yb ~ forest(x1 + x2) + forest(x1, basis = dose),
  d,
  weights = rep(0, n),
  family = "probit",
  control = printControl()
))
fallbackFit <- collected(bartFit(
  yb ~ forest(x1 + x2) + forest(x1, basis = dose),
  d,
  weights = rep(0, n),
  family = "probit"
))
fallbackWarnings <- takeCreationWarnings()
expect_identical(length(fallbackWarnings), 2L)
expect_true(all(grepl(
  "no observation is in the likelihood",
  fallbackWarnings,
  fixed = TRUE
)))
expect_identical(fullNow(fallbackSampler), fallbackFull)
expect_identical(
  shortOf(showLines(fallbackSampler)),
  c(fallbackShort, samplerFooter)
)
expect_identical(
  shortOf(printLines(fallbackFit)),
  c(fallbackShort, fitFooter)
)
expect_identical(fullOf(summaryLines(fallbackFit)), fallbackFull)

# --- a 0/1 number and an age at their defaults ---
zeroMade <- built(dbarts(
  y ~ forest(x1 + x2) + forest(x1, basis = z) + forest(x1, basis = age),
  d,
  control = printControl(verbose = TRUE)
))
zero <- zeroMade$sampler
zeroFit <- counted(bart(
  y ~ forest(x1 + x2) + forest(x1, basis = z) + forest(x1, basis = age),
  d,
  n.samples = 5L,
  n.burn = 0L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 5L,
  verbose = FALSE
))
zeroOne <- zeroMade$full[3L]
expect_identical(zeroOne, numericLine(2, "z", z, U / sd(z), "default"))
zeroShort <- numericShort(2, "z", z, U / sd(z))
expect_identical(shortOf(showLines(zero))[3L], zeroShort)
expect_identical(shortOf(printLines(zeroFit))[3L], zeroShort)
expect_identical(fullOf(summaryLines(zeroFit))[3L], zeroOne)

# --- one forest: nothing is printed, and a summary is as it was ---
singleOut <- verboseLines(
  single <- dbarts(y ~ x1 + x2, d, control = printControl(verbose = TRUE))
)
expect_null(fullOf(singleOut))
expect_null(shortOf(singleOut))
expect_null(shortOf(showLines(single)))
oneFit <- counted(bart(
  y ~ x1 + x2,
  d,
  n.samples = 5L,
  n.burn = 0L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 5L,
  verbose = FALSE
))
expect_null(shortOf(printLines(oneFit)))
oneSummary <- summaryLines(oneFit)
expect_null(fullOf(oneSummary))
expect_null(summary(oneFit)$forests)
expect_true("(Fixed, not sampled: k = 2)" %in% oneSummary)

# --- no warning anywhere above, an age and a 0/1 number included, at show,
# print and summary of a fit as well as at a creation ---
expect_identical(numWarnings, 0L)
