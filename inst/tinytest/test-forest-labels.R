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
# the per-forest margins are named by the labels, and forest<i> and a label
# both select by position
expect_identical(
  dimnames(labelled$forestFits)[[length(dim(labelled$forestFits))]],
  c("forest1", "dose")
)
# one forest has one name everywhere: its label, and forest<i> for a forest
# with none of its own, on the fit, in extract, fitted and predict
heldFit <- bart(
  y ~ forest(x1 + x2) + forest(x1 + x2, basis = dose, amplitude = fixed()),
  frame,
  n.chains = 2L,
  n.threads = 1L,
  n.trees = 5L,
  n.samples = 3L,
  n.burn = 1L,
  keepTrees = TRUE,
  verbose = FALSE
)
expect_identical(attr(heldFit, "forest.labels"), c("forest1", "dose"))
expect_identical(heldFit$fixed$amplitude, c(dose = 1))
expect_identical(names(heldFit$fixed$k), c("forest1", "dose"))
expect_identical(names(heldFit$leaf.prior), c("forest1", "dose"))
expect_identical(names(extract(heldFit, "k")), c("forest1", "dose"))
expect_identical(names(extract(heldFit, "leaf.prior.sd")), c("forest1", "dose"))
expect_identical(
  names(summary(heldFit, vars = "k")$fixed$k),
  c("forest1", "dose")
)
expect_identical(
  names(summary(heldFit)$fixed$leaf.prior.sd),
  c("forest1", "dose")
)
lastMargin <- function(x) dimnames(x)[[length(dim(x))]]
expect_identical(lastMargin(heldFit$forestFits), c("forest1", "dose"))
expect_identical(unique(attr(heldFit$glue, "forest")), c("forest1", "dose"))
expect_identical(lastMargin(heldFit$varcount), c("forest1", "dose"))
expect_identical(lastMargin(extract(heldFit, "forest")), c("forest1", "dose"))
expect_identical(
  lastMargin(extract(heldFit, "forest", contribution = TRUE)),
  c("forest1", "dose")
)
expect_identical(
  lastMargin(extract(heldFit, "forest", combineChains = FALSE)),
  c("forest1", "dose")
)
expect_identical(
  colnames(fitted(heldFit, type = "forest")),
  c("forest1", "dose")
)
expect_identical(
  lastMargin(predict(heldFit, frame, type = "forest")),
  c("forest1", "dose")
)
# forest<i> still selects position i wherever a forest is taken, as a label does
expect_identical(
  extract(heldFit, "forest", forest = "forest2"),
  extract(heldFit, "forest", forest = "dose")
)
expect_identical(
  lastMargin(extract(heldFit, "forest", forest = "forest2")),
  "dose"
)
expect_identical(
  predict(heldFit, frame, type = "forest", forest = "forest2"),
  predict(heldFit, frame, type = "forest", forest = 2L)
)
expect_identical(
  extract(heldFit, "k", forest = "forest2"),
  extract(heldFit, "k", forest = "dose")
)
expect_identical(
  extract(heldFit, "leaf.prior.sd", forest = "forest2"),
  extract(heldFit, "leaf.prior.sd", forest = 2L)
)
# a forest whose basis was handed over as a value has no label of its own and
# keeps forest<i>
valueSampler <- listFit(bquote(list(
  forest(),
  forest(basis = .(frame$dose), amplitude = fixed())
)))
valueBurn <- dbarts:::runWithBurnIn(valueSampler, valueSampler$control, TRUE)
valueFit <- dbarts:::packageBartResults(
  valueSampler,
  valueBurn$samples,
  valueBurn$burnInSigma,
  valueBurn$burnInK,
  TRUE,
  TRUE
)
expect_identical(attr(valueFit, "forest.labels"), c("forest1", "forest2"))
expect_identical(valueFit$fixed$amplitude, c(forest2 = 1))
expect_identical(names(valueFit$fixed$k), c("forest1", "forest2"))
expect_identical(names(valueFit$leaf.prior), c("forest1", "forest2"))
expect_identical(lastMargin(valueFit$forestFits), c("forest1", "forest2"))
expect_identical(unique(attr(valueFit$glue, "forest")), c("forest1", "forest2"))
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
  listFit(quote(list(forest(), forest(basis = dose))))
}
replaced <- created()
expect_identical(colnames(replaced$data@bases[[2L]]), "dose")
# a value: the recorded name stays
replacement <- frame$age / 10
replaced$setForestBasis(2L, replacement)
expect_identical(
  replaced$data@bases[[2L]],
  matrix(replacement, n, dimnames = list(NULL, "dose"))
)
expect_identical(labelsOf(replaced), c("forest1", "dose"))
# another name is ignored
replaced$setForestBasis(2L, cbind(first = frame$age))
expect_identical(colnames(replaced$data@bases[[2L]]), "dose")
expect_identical(unname(replaced$data@bases[[2L]][, 1L]), frame$age)
# a forest of numbers is one column for life: two are refused, nothing being
# installed
before <- replaced$data@bases[[2L]]
expect_error(
  replaced$setForestBasis(2L, cbind(dose = frame$dose, age = frame$age)),
  "forest 2 ('dose') multiplies one numeric column, and 'basis' gives it 2",
  fixed = TRUE
)
expect_identical(replaced$data@bases[[2L]], before)
# a one-sided formula is read as a basis is, and its names are found where it
# was written
doubled <- frame$dose * 2
halved <- frame$age / 2
replaced$setForestBasis(2L, ~halved)
expect_identical(replaced$data@bases[[2L]], cbind(dose = halved))
expect_identical(labelsOf(replaced), c("forest1", "dose"))
expect_error(
  replaced$setForestBasis(2L, ~ doubled + halved),
  "forest 2 ('dose') multiplies one numeric column, and 'basis' gives it 2",
  fixed = TRUE
)
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
  replaced$setForestBasis(2L, ~ halved[1:10]),
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
  drawsAfter(function(sampler) sampler$setForestBasis(2L, ~halved)),
  drawsAfter(function(sampler) sampler$setForestBasis(2L, halved))
)
expect_false(identical(
  drawsAfter(function(sampler) sampler$setForestBasis(2L, ~halved)),
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
levelled$setForestBasis(2L, ~ factor(rep(0L, n), levels = 0:1))
expect_identical(unname(colSums(levelled$data@bases[[2L]])), c(n + 0, 0))
expect_identical(
  colnames(levelled$data@bases[[2L]]),
  c("factor(z)0", "factor(z)1")
)
expect_identical(labelsOf(levelled), c("forest1", "factor(z)"))

## --- Block D: the forest column of a fit's tree table -----------------------
# the column is a factor whose levels are the forests' labels in forest order,
# forest<i> for a forest with none (dec-B484), so that trees$forest == "dose"
# selects and as.integer gives the position; the sampler's own getTrees keeps
# the position number
trees <- extract(heldFit, "trees")
expect_true(is.factor(trees$forest))
expect_identical(levels(trees$forest), c("forest1", "dose"))
expect_identical(levels(trees$forest), attr(heldFit, "forest.labels"))
expect_identical(
  colnames(trees),
  c("forest", "chain", "sample", "tree", "n", "var", "value")
)
expect_true(is.integer(trees$var))
numbered <- heldFit$fit$getTrees()
expect_true(is.integer(numbered$forest))
expect_identical(as.integer(trees$forest), numbered$forest)
expect_identical(
  trees$forest == "dose",
  numbered$forest == 2L
)
expect_identical(unname(table(trees$forest)), unname(table(numbered$forest)))
# whatever forests the frame holds, the levels are all of them
dose <- extract(heldFit, "trees", forest = "dose")
expect_identical(levels(dose$forest), c("forest1", "dose"))
expect_true(all(dose$forest == "dose"))
expect_identical(
  levels(extract(heldFit, "trees", forest = 1L)$forest),
  levels(trees$forest)
)
expect_identical(
  levels(extract(heldFit, "trees", forest = "forest2")$forest),
  levels(trees$forest)
)
expect_identical(
  unique(as.integer(extract(heldFit, "trees", forest = c(2L, 1L))$forest)),
  c(2L, 1L)
)
# and the frame is the sampler's but for that column
positions <- trees
positions$forest <- as.integer(positions$forest)
expect_identical(positions[colnames(numbered)], numbered)
# at newdata, and at the current trees
expect_identical(
  levels(extract(heldFit, "trees", newdata = frame[1:5, ])$forest),
  c("forest1", "dose")
)
expect_identical(
  levels(extract(heldFit, "trees", current = TRUE)$forest),
  c("forest1", "dose")
)
# a label a list names is the level, in the order written
packaged <- function(sampler) {
  burn <- dbarts:::runWithBurnIn(sampler, sampler$control, TRUE)
  dbarts:::packageBartResults(
    sampler,
    burn$samples,
    burn$burnInSigma,
    burn$burnInK,
    TRUE,
    TRUE
  )
}
named <- packaged(listFit(quote(list(
  prognostic = forest(),
  treatment = forest(basis = factor(z))
))))
expect_identical(
  levels(extract(named, "trees")$forest),
  c("prognostic", "treatment")
)
# a fit with one forest declared still has the column, its one level the name a
# forest with no recorded label has, since a model of one records none; one
# never declared has no column; a multinomial fit's category forests are
# forest<i>
declared <- packaged(listFit(quote(list(only = forest()))))
expect_identical(levels(extract(declared, "trees")$forest), "forest1")
expect_null(
  extract(
    bart(
      y ~ x1 + x2,
      frame,
      n.chains = 1L,
      n.threads = 1L,
      n.trees = 5L,
      n.samples = 3L,
      n.burn = 1L,
      keepTrees = TRUE,
      verbose = FALSE
    ),
    "trees"
  )$forest
)
frame$category <- factor(rep_len(c("p", "q", "r"), n))
multinomial <- bart(
  category ~ x1 + x2,
  frame,
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 5L,
  n.samples = 3L,
  n.burn = 1L,
  keepTrees = TRUE,
  verbose = FALSE
)
expect_identical(
  levels(extract(multinomial, "trees")$forest),
  c("forest1", "forest2", "forest3")
)
expect_true(is.integer(multinomial$fit$getTrees()$forest))
