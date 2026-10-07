# Every forest of a model of several has a label, fixed when the sampler is
# created: the list name where one is written, else the text of the forest's
# basis, else forest<i>. The columns of a basis have names fixed the same
# way, which a replacement by $setForestBasis does not change. Block A: the
# label rule at both doors. Block B: what it refuses. Block C: a replacement
# basis.

forest <- dbartsForests$forest

set.seed(67)
n <- 70L
frame <- data.frame(
  x1 = runif(n),
  x2 = runif(n),
  x3 = runif(n),
  dose = runif(n, 0.5, 2),
  age = rnorm(n, 40, 10),
  z = rep_len(c(0L, 1L), n)
)
frame$y <- with(frame, x1 + z * (1 + x2) + 0.2 * dose + rnorm(n, sd = 0.3))
x <- as.matrix(frame[c("x1", "x2", "x3")])
y <- frame$y

labelControl <- function() {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 5L,
    n.samples = 3L,
    n.burn = 1L,
    updateState = FALSE,
    verbose = FALSE,
    seed = 67L
  )
}
labelsOf <- function(sampler) {
  attr(sampler$control, "bartcore.forests", exact = TRUE)$labels
}
termFit <- function(formula, data = frame) {
  dbarts(formula, data, control = labelControl())
}
listFit <- function(forests, data = frame) {
  eval(
    bquote(dbarts(
      y ~ x1 + x2 + x3,
      data,
      forests = .(forests),
      control = labelControl()
    )),
    list(data = data),
    parent.frame()
  )
}

## --- Block A: the label rule ------------------------------------------------
# the text of the basis, the same text at either door
expect_identical(
  labelsOf(termFit(y ~ x1 + x2 + forest(x3, basis = dose))),
  c("forest1", "dose")
)
expect_identical(
  labelsOf(listFit(quote(list(forest(), forest(basis = dose))))),
  c("forest1", "dose")
)
# as R writes the code, however it was spaced, the tilde taken off
for (written in list(
  quote(I(dose / 30)),
  quote(I(dose / 30)),
  quote(~ I(dose / 30)),
  str2lang("I( dose/30 )")
)) {
  expect_identical(
    labelsOf(listFit(bquote(list(forest(), forest(basis = .(written)))))),
    c("forest1", "I(dose/30)")
  )
  expect_identical(
    labelsOf(termFit(eval(bquote(
      y ~ x1 + x2 + forest(x3, basis = .(written))
    )))),
    c("forest1", "I(dose/30)")
  )
}
expect_identical(
  labelsOf(termFit(y ~ forest(x1 + x2) + forest(x3, basis = scale(age)))),
  c("forest1", "scale(age)")
)
# a name that holds a formula is replaced by the formula, and is the label
# itself only where it is a column of the data
held <- ~ scale(age)
expect_identical(
  labelsOf(termFit(y ~ x1 + x2 + forest(x3, basis = held))),
  c("forest1", "scale(age)")
)
expect_identical(
  labelsOf(listFit(quote(list(forest(), forest(basis = held))))),
  c("forest1", "scale(age)")
)
withHeld <- frame
withHeld$held <- frame$dose
expect_identical(
  labelsOf(termFit(y ~ x1 + x2 + forest(x3, basis = held), withHeld)),
  c("forest1", "held")
)
expect_identical(
  labelsOf(listFit(quote(list(forest(), forest(basis = held))), withHeld)),
  c("forest1", "held")
)
# a list name is the label where one is written
expect_identical(
  labelsOf(listFit(quote(list(
    prognostic = forest(),
    treatment = forest(basis = factor(z))
  )))),
  c("prognostic", "treatment")
)
expect_identical(
  labelsOf(listFit(quote(list(forest(), effect = forest(basis = dose))))),
  c("forest1", "effect")
)
expect_identical(
  labelsOf(listFit(quote(list(base = forest(), forest(basis = dose))))),
  c("base", "dose")
)
# its own position's default is a name like any other
expect_identical(
  labelsOf(listFit(quote(list(
    forest1 = forest(),
    forest2 = forest(basis = dose)
  )))),
  c("forest1", "forest2")
)
# no basis, a basis handed over as a value, a basis from the data object: the
# position
doseValue <- frame$dose
expect_identical(
  labelsOf(listFit(quote(list(
    forest(),
    do.call(forest, list(basis = doseValue))
  )))),
  c("forest1", "forest2")
)
expect_identical(
  labelsOf(dbarts(
    dbartsData(x, y, bases = list(NULL, frame$dose, frame$age)),
    control = labelControl()
  )),
  c("forest1", "forest2", "forest3")
)
# a text that repeats an earlier label takes the suffix the names of a data
# frame take; an earlier forest is never renamed
expect_identical(
  labelsOf(termFit(
    y ~ forest(x1 + x2) +
      forest(x1, basis = dose) +
      forest(x2, basis = dose) +
      forest(x3, basis = age)
  )),
  c("forest1", "dose", "dose.1", "age")
)
expect_identical(
  labelsOf(listFit(quote(list(
    forest(),
    forest(x1, basis = dose),
    dose = forest(x2, basis = age),
    forest(x3, basis = dose)
  )))),
  c("forest1", "dose.1", "dose", "dose.2")
)
# every forest of a formula with a basis keeps the order written; with no
# forest for the control's count to be the count of, the first states its own
expect_identical(
  labelsOf(dbarts(
    y ~ forest(x1, basis = dose, n.trees = 5L, base = 0.95, power = 2) +
      forest(x2, basis = age),
    frame,
    control = dbartsControl(
      n.chains = 1L,
      n.threads = 1L,
      n.samples = 3L,
      n.burn = 1L,
      updateState = FALSE,
      verbose = FALSE,
      seed = 67L
    )
  )),
  c("dose", "age")
)
# the fit carries them, always
labelled <- bart(
  y ~ forest(x1 + x2) + forest(x3, basis = dose),
  frame,
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 5L,
  n.samples = 3L,
  n.burn = 0L,
  verbose = FALSE
)
expect_identical(attr(labelled, "forest.labels"), c("forest1", "dose"))
# the per-forest margins keep forest1, forest2, and a label selects by position
expect_identical(
  dimnames(labelled$forestFits)[[length(dim(labelled$forestFits))]],
  c("forest1", "forest2")
)
sampler <- listFit(quote(list(forest(), forest(basis = dose))))
expect_identical(sampler$getLeafPrior("dose"), sampler$getLeafPrior(2L))
expect_error(
  sampler$getLeafPrior("age"),
  "'forest' names no forest of this model: \"age\"",
  fixed = TRUE
)
# the writer checks a list's names against them
expect_silent(sampler$setLeafPrior(
  forests = list(forest(), dose = forest(sd = 1.5))
))
expect_error(
  sampler$setLeafPrior(forests = list(forest(), age = forest(sd = 1.5))),
  "names forest 2 'age', but it was created as 'dose'",
  fixed = TRUE
)

## --- Block B: what the label rule refuses -----------------------------------
expect_error(
  listFit(quote(list(a = forest(), a = forest(basis = dose)))),
  "'forests' names two forests \"a\"; a label names one forest",
  fixed = TRUE
)
expect_error(
  listFit(quote(list(forest(), forest1 = forest(basis = dose)))),
  paste0(
    "'forests' names forest 2 \"forest1\", which is the label an unnamed ",
    "forest 1 has; choose another name"
  ),
  fixed = TRUE
)
expect_error(
  listFit(quote(list(forest3 = forest(), forest(basis = dose)))),
  "'forests' names forest 1 \"forest3\", which is the label an unnamed forest 3 has",
  fixed = TRUE
)

## --- Block C: a replacement basis -------------------------------------------
created <- function() {
  listFit(quote(list(forest(), forest(basis = dose + age))))
}
replaced <- created()
expect_identical(colnames(replaced$data@bases[[2L]]), c("dose", "age"))
# a value: columns are taken by position, and the recorded names stay
replacement <- cbind(frame$age / 10, frame$dose * 2)
replaced$setForestBasis(2L, replacement)
expect_identical(
  replaced$data@bases[[2L]],
  matrix(as.double(replacement), n, dimnames = list(NULL, c("dose", "age")))
)
expect_identical(labelsOf(replaced), c("forest1", "dose + age"))
# other names are ignored
replaced$setForestBasis(2L, cbind(first = frame$age, second = frame$dose))
expect_identical(colnames(replaced$data@bases[[2L]]), c("dose", "age"))
expect_identical(unname(replaced$data@bases[[2L]][, 1L]), frame$age)
replaced$setForestBasis(2L, cbind(dose = frame$age, other = frame$dose))
expect_identical(colnames(replaced$data@bases[[2L]]), c("dose", "age"))
# the recorded names in another order are refused, nothing being installed
before <- replaced$data@bases[[2L]]
expect_error(
  replaced$setForestBasis(2L, cbind(age = frame$age, dose = frame$dose)),
  paste0(
    "'basis' has the columns of forest 2's basis in another order (age, ",
    "dose; the forest has dose, age): columns are taken by position, so ",
    "give them in the forest's order"
  ),
  fixed = TRUE
)
expect_identical(replaced$data@bases[[2L]], before)
# in the forest's order they are taken
replaced$setForestBasis(2L, cbind(dose = frame$dose, age = frame$age))
expect_identical(
  replaced$data@bases[[2L]],
  cbind(dose = frame$dose, age = frame$age)
)
# a one-sided formula is read as a basis is: '+' is columns, and its names
# are found where it was written
doubled <- frame$dose * 2
halved <- frame$age / 2
replaced$setForestBasis(2L, ~ doubled + halved)
expect_identical(
  replaced$data@bases[[2L]],
  cbind(dose = doubled, age = halved)
)
expect_identical(labelsOf(replaced), c("forest1", "dose + age"))
expect_error(
  replaced$setForestBasis(2L, ~ doubled * halved),
  "'basis' does not take '*' between its terms ('doubled * halved')",
  fixed = TRUE
)
expect_error(
  replaced$setForestBasis(2L, y ~ doubled),
  "a 'basis' formula must be one-sided, as ~ dose",
  fixed = TRUE
)
expect_error(
  replaced$setForestBasis(2L, ~ doubled[1:10] + halved[1:10]),
  "must have the same length as the data: it has 10 rows and the data 70",
  fixed = TRUE
)
# the run that follows uses the replacement: the same draws as a sampler
# given the same replacement as a value
drawsAfter <- function(install) {
  sampler <- created()
  install(sampler)
  run <- sampler$run(0L, 8L)
  list(run$train, run$sigma, sampler$getForestAmplitudes())
}
expect_identical(
  drawsAfter(function(sampler) sampler$setForestBasis(2L, ~ doubled + halved)),
  drawsAfter(function(sampler) {
    sampler$setForestBasis(2L, cbind(doubled, halved))
  })
)
expect_false(identical(
  drawsAfter(function(sampler) sampler$setForestBasis(2L, ~ doubled + halved)),
  drawsAfter(function(sampler) NULL)
))
# a factor basis: the levels' names stay too, and a level may be left empty
levelled <- listFit(quote(list(forest(), forest(basis = factor(z)))))
expect_identical(
  colnames(levelled$data@bases[[2L]]),
  c("factor(z)0", "factor(z)1")
)
flipped <- factor(1L - frame$z)
levelled$setForestBasis(2L, flipped)
expect_identical(
  levelled$data@bases[[2L]],
  matrix(
    as.double(c(frame$z, 1L - frame$z)),
    n,
    dimnames = list(NULL, c("factor(z)0", "factor(z)1"))
  )
)
levelled$setForestBasis(2L, ~ factor(rep("a", n), levels = c("a", "b")))
expect_identical(unname(colSums(levelled$data@bases[[2L]])), c(n + 0, 0))
# a replacement of another width brings its own names
levelled$setForestBasis(2L, cbind(u = frame$dose, v = frame$age, w = frame$x1))
expect_identical(colnames(levelled$data@bases[[2L]]), c("u", "v", "w"))
expect_identical(labelsOf(levelled), c("forest1", "factor(z)"))
