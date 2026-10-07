# A forest is selected by its label as well as by its position. A number is a
# position, read as it always was; a string is a label, else the same code as
# one, and forest<i> is the name position i has on every per-forest margin;
# never a position. Block A: the reader alone. Block B: the sampler's nine
# methods. Block C: a fit. Block D: lists given forest by forest.

forest <- dbartsForests$forest
select <- dbarts:::selectForest

refusal <- function(expr) {
  tryCatch(
    {
      force(expr)
      "accepted"
    },
    error = conditionMessage
  )
}

## --- Block A: the reader alone ----------------------------------------------
# a sampler of one forest records no labels: forest<i> is taken where the
# number is, and any other string is refused
expect_identical(select("forest1", NULL, 1L), 1L)
expect_identical(select(1L, NULL, 1L), 1L)
expect_identical(select(2, NULL, 1L), 2)
expect_identical(
  refusal(select("dose", NULL, 1L)),
  paste0(
    "'forest' names no forest of this model: \"dose\"; this sampler's ",
    "forests have no labels, so select one by position"
  )
)
expect_identical(
  refusal(select("forest2", NULL, 1L)),
  paste0(
    "'forest' names no forest of this model: \"forest2\"; this sampler's ",
    "forests have no labels, so select one by position"
  )
)
expect_identical(select(c("forest1", "forest3"), NULL, 3L, TRUE), c(1L, 3L))

# two forests, the second labelled by its basis
labels <- c("forest1", "dose")
expect_identical(select("dose", labels, 2L), 2L)
expect_identical(select("forest1", labels, 2L), 1L)
# forest<i> is position i whatever that position is labelled
expect_identical(select("forest2", labels, 2L), 2L)
expect_identical(
  refusal(select("age", labels, 2L)),
  paste0(
    "'forest' names no forest of this model: \"age\"; its forests are ",
    "\"forest1\", \"dose\""
  )
)
expect_identical(
  refusal(select("forest3", labels, 2L)),
  paste0(
    "'forest' names no forest of this model: \"forest3\"; its forests ",
    "are \"forest1\", \"dose\""
  )
)

# three forests, labelled by code
labels <- c("forest1", "scale(age)", "I(dose/30)")
expect_identical(select("I(dose/30)", labels, 3L), 3L)
# a label matches exactly, then as code
expect_identical(select("I(dose / 30)", labels, 3L), 3L)
expect_identical(select("scale( age )", labels, 3L), 2L)
expect_identical(
  select(c("I(dose/30)", "forest1"), labels, 3L, TRUE),
  c(3L, 1L)
)
expect_identical(
  refusal(select("age", labels, 3L)),
  paste0(
    "'forest' names no forest of this model: \"age\"; its forests are ",
    "\"forest1\", \"scale(age)\", \"I(dose/30)\""
  )
)
# a string is never a position
expect_identical(
  refusal(select("2", labels, 3L)),
  paste0(
    "'forest' names no forest of this model: \"2\"; its forests are ",
    "\"forest1\", \"scale(age)\", \"I(dose/30)\"; a position is given as a ",
    "number, forest = 2"
  )
)
# a number is a position, handed on as it came
expect_identical(select(2L, labels, 3L), 2L)
expect_identical(select(2, labels, 3L), 2)
expect_identical(select(c(3, 1), labels, 3L, TRUE), c(3, 1))

# a list name is a label, and so is a label of digits
labels <- c("2", "1")
expect_identical(select("2", labels, 2L), 1L)
expect_identical(select("1", labels, 2L), 2L)
expect_identical(select(2L, labels, 2L), 2L)
expect_identical(select(1L, labels, 2L), 1L)

# one string both a label and another position's name, and one the same code
# as two labels
labels <- c("forest1", "forest3", "a+b", "a + b")
expect_identical(
  refusal(select("forest3", labels, 4L)),
  paste0(
    "'forest' (\"forest3\") is the label of forest 2 and the name of ",
    "position 3; select by position, as forest = 2"
  )
)
expect_identical(
  refusal(select("a +b", labels, 4L)),
  paste0(
    "'forest' (\"a +b\") is the label of forests 3 and 4 (\"a+b\", ",
    "\"a + b\"); give one exactly, or select by position"
  )
)
# exactly beats code
expect_identical(select("a+b", labels, 4L), 3L)
expect_identical(select("a + b", labels, 4L), 4L)
# the name of a position nobody else is labelled is that position
expect_identical(select("forest1", labels, 4L), 1L)
expect_identical(select("forest4", labels, 4L), 4L)

# forest0 and the lower bound of forest<i> are positions that do not exist
labels <- c("forest1", "dose")
expect_identical(
  refusal(select("forest0", labels, 2L)),
  paste0(
    "'forest' names no forest of this model: \"forest0\"; its forests are ",
    "\"forest1\", \"dose\""
  )
)
expect_identical(select("forest1", labels, 2L), 1L)
# three or more labels of one code
labels <- c("a+b", "a + b", "a  +  b", "forest4")
expect_identical(
  refusal(select("a +b", labels, 4L)),
  paste0(
    "'forest' (\"a +b\") is the label of forests 1, 2 and 3 (\"a+b\", ",
    "\"a + b\", \"a  +  b\"); give one exactly, or select by position"
  )
)
# a label that is not R code is found only exactly, and a string that is not
# R code finds no label by code
labels <- c("forest1", "(", "a b c")
expect_identical(select("(", labels, 3L), 2L)
expect_identical(select("a b c", labels, 3L), 3L)
expect_match(refusal(select("( ", labels, 3L)), "names no forest", fixed = TRUE)
expect_match(
  refusal(select("a  b c", labels, 3L)),
  "names no forest",
  fixed = TRUE
)
expect_match(refusal(select("b", labels, 3L)), "names no forest", fixed = TRUE)

# NA, an empty string, a factor, a logical, a list
labels <- c("forest1", "dose")
for (bad in list(NA, NA_character_, "", c("dose", NA))) {
  expect_identical(
    refusal(select(bad, labels, 2L, TRUE)),
    "'forest' must not be NA or an empty string",
    info = deparse(bad)
  )
}
for (case in list(
  list(TRUE, "a logical"),
  list(c(TRUE, FALSE), "a logical"),
  list(factor("dose"), "a factor"),
  list(list(2), "a list"),
  list(1i, "a complex")
)) {
  expect_identical(
    refusal(select(case[[1L]], labels, 2L)),
    paste0(
      "'forest' must be a number, the forest's position, or a string, its ",
      "label; not ",
      case[[2L]]
    )
  )
}
# one label where one forest is meant
expect_identical(
  refusal(select(c("dose", "forest1"), labels, 2L)),
  "'forest' must be a single number or a single label"
)

## --- Block B: the sampler's nine methods ------------------------------------
set.seed(83)
n <- 60L
frame <- data.frame(
  x1 = runif(n),
  x2 = runif(n),
  x3 = runif(n),
  dose = runif(n, 0.5, 2),
  age = rnorm(n, 40, 10)
)
frame$y <- with(frame, x1 + x2 * x3 + 0.2 * dose + rnorm(n, sd = 0.3))
x <- as.matrix(frame[c("x1", "x2", "x3")])
y <- frame$y

selectionControl <- function(...) {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.samples = 4L,
    n.burn = 1L,
    keepTrees = TRUE,
    updateState = FALSE,
    verbose = FALSE,
    seed = 83L,
    ...
  )
}
# the second forest is the plain one; each has its own count and spread, so
# that every reader answers differently by forest
fourForests <- function() {
  dbarts(
    y ~ x1 + x2 + x3,
    frame,
    forests = list(
      a = forest(basis = dose, n.trees = 3L, sd = 1),
      forest(n.trees = 4L, sd = 2),
      forest(basis = scale(age), n.trees = 5L, sd = 3),
      forest(basis = I(dose / 30), n.trees = 6L, sd = 4)
    ),
    control = selectionControl()
  )
}
sampler <- fourForests()
invisible(sampler$run())
selectable <- c("a", "forest2", "scale(age)", "I(dose/30)")
expect_identical(
  attr(sampler$control, "bartcore.forests", exact = TRUE)$labels,
  selectable
)

readers <- list(
  getLeafPrior = function(s, f) s$getLeafPrior(f),
  getK = function(s, f) s$getK(f),
  getForestFits = function(s, f) s$getForestFits(f),
  getForestAmplitudes = function(s, f) s$getForestAmplitudes(f),
  getForestVariableCounts = function(s, f) s$getForestVariableCounts(f),
  getTrees = function(s, f) s$getTrees(forest = f)
)
for (name in names(readers)) {
  read <- readers[[name]]
  for (index in seq_along(selectable)) {
    byPosition <- read(sampler, index)
    expect_identical(
      read(sampler, selectable[index]),
      byPosition,
      info = paste(name, index)
    )
    expect_identical(
      read(sampler, paste0("forest", index)),
      byPosition,
      info = paste(name, index)
    )
    expect_identical(read(sampler, as.double(index)), byPosition)
  }
}
# the shared resolver says which forest, so that readers whose answer is the
# same for every forest (getK, k being pinned at 1) are still told apart
for (index in seq_along(selectable)) {
  for (spelling in list(selectable[index], paste0("forest", index), index)) {
    expect_identical(
      dbarts:::samplerForestIndex(
        spelling,
        sampler$control,
        sampler$getPointer()
      ),
      index - 1L
    )
  }
}
# and the readers differ by forest, so that the identity is not vacuous
for (name in c("getLeafPrior", "getForestFits", "getTrees")) {
  read <- readers[[name]]
  for (index in 2:4) {
    expect_false(
      identical(read(sampler, 1L), read(sampler, index)),
      info = paste(name, index)
    )
  }
}
expect_identical(
  vapply(1:4, function(i) sampler$getLeafPrior(i)$leaf.prior$sd, 0),
  c(1, 2, 3, 4)
)
# getTrees takes a vector of labels, in the order given
expect_identical(
  sampler$getTrees(forest = c("I(dose/30)", "a")),
  sampler$getTrees(forest = c(4L, 1L))
)
expect_identical(
  sampler$getTrees(forest = c("forest2", "I(dose / 30)")),
  sampler$getTrees(forest = c(2L, 4L))
)
expect_identical(
  unique(sampler$getTrees(forest = c("I(dose/30)", "a"))$forest),
  c(4L, 1L)
)
# plotTree takes one
pdf(NULL)
expect_silent(sampler$plotTree(1L, chainNum = 1L, forest = "scale(age)"))
dev.off()

# plotTree takes one forest, by its own check
pdf(NULL)
expect_error(
  sampler$plotTree(1L, chainNum = 1L, forest = c("a", "forest2")),
  "'forest' must be a single number or a single label",
  fixed = TRUE
)
expect_error(
  sampler$plotTree(1L, chainNum = 1L, forest = "dose"),
  "'forest' names no forest of this model",
  fixed = TRUE
)
dev.off()

# the writers change the forest named and no other
weights <- runif(n)
named <- fourForests()
placed <- fourForests()
named$setForestWeights(weights, forest = "scale(age)")
placed$setForestWeights(weights, forest = 3L)
expect_identical(named$forestWeights, placed$forestWeights)
expect_identical(named$forestWeights[[3L]], weights)
expect_true(all(vapply(named$forestWeights[-3L], is.null, NA)))
named$setForestWeights(weights * 2, forest = "forest1")
expect_identical(named$forestWeights[[3L]], weights)
expect_identical(named$forestWeights[[1L]], weights * 2)
before <- named$data@bases
swap <- frame$dose * 2
named$setForestBasis("a", swap)
placed$setForestBasis(1L, swap)
expect_identical(named$data@bases, placed$data@bases)
expect_identical(named$data@bases[-1L], before[-1L])
expect_false(identical(named$data@bases[[1L]], before[[1L]]))
named$setForestBasis("I(dose / 30)", frame$dose * 3)
expect_identical(named$data@bases[2:3], before[2:3])
expect_equivalent(named$data@bases[[4L]][, 1L], frame$dose * 3)

# what a number or a string names that is not a forest, and what is refused
# by kind
for (reader in readers) {
  expect_match(
    refusal(reader(sampler, "dose")),
    "'forest' names no forest of this model: \"dose\"; its forests are \"a\", ",
    fixed = TRUE
  )
  expect_match(
    refusal(reader(sampler, "2")),
    "a position is given as a number, forest = 2",
    fixed = TRUE
  )
  for (bad in list(TRUE, factor("a"), list(2))) {
    expect_match(refusal(reader(sampler, bad)), "not a (logical|factor|list)$")
  }
  expect_identical(
    refusal(reader(sampler, NA)),
    "'forest' must not be NA or an empty string"
  )
  expect_identical(
    refusal(reader(sampler, "")),
    "'forest' must not be NA or an empty string"
  )
  expect_match(
    refusal(reader(sampler, 5L)),
    "forest index out of range|index out of range"
  )
  expect_identical(
    refusal(reader(sampler, 0L)),
    "'forest' must be a single positive integer (1 selects the first)"
  )
}
expect_match(
  refusal(sampler$setForestWeights(weights, forest = "2")),
  "a position is given as a number",
  fixed = TRUE
)
expect_match(
  refusal(sampler$setForestBasis(factor("a"), frame$dose)),
  "not a factor",
  fixed = TRUE
)
expect_match(refusal(sampler$setForestBasis(TRUE, frame$dose)), "not a logical")

# the label and the sampler go together through a copy and a restore
copied <- sampler$copy()
expect_identical(copied$getLeafPrior("scale(age)"), sampler$getLeafPrior(3L))

# and through a save and a read
path <- tempfile(fileext = ".rds")
sampler$storeState()
saveRDS(sampler, path)
restored <- readRDS(path)
expect_identical(
  restored$getLeafPrior("scale(age)"),
  sampler$getLeafPrior(3L)
)
unlink(path)

# a sampler of one forest has no labels, and forest1 is taken
single <- dbarts(y ~ x1 + x2 + x3, frame, control = selectionControl())
expect_identical(single$getLeafPrior("forest1"), single$getLeafPrior(1L))
expect_match(
  refusal(single$getLeafPrior("x1")),
  "this sampler's forests have no labels, so select one by position",
  fixed = TRUE
)
expect_match(
  refusal(single$getLeafPrior("forest2")),
  "this sampler's forests have no labels",
  fixed = TRUE
)
declared <- dbarts(
  y ~ x1 + x2 + x3,
  frame,
  forests = list(forest()),
  control = selectionControl()
)
expect_identical(declared$getK("forest1"), declared$getK(1L))
expect_match(
  refusal(declared$getK("dose")),
  "have no labels",
  fixed = TRUE
)

# a multinomial sampler has none either
category <- factor(sample(c("p", "q", "r"), n, TRUE))
multinomial <- bart(
  category ~ x1 + x2,
  frame,
  n.chains = 1L,
  n.threads = 1L,
  n.samples = 4L,
  n.burn = 1L,
  keepTrees = TRUE,
  verbose = FALSE
)
expect_identical(
  multinomial$fit$getLeafPrior("forest2"),
  multinomial$fit$getLeafPrior(2L)
)
expect_match(
  refusal(multinomial$fit$getLeafPrior("q")),
  "this sampler's forests have no labels",
  fixed = TRUE
)
expect_identical(
  extract(multinomial, "trees", forest = "forest2"),
  extract(multinomial, "trees", forest = 2L)
)
expect_match(
  refusal(extract(multinomial, "trees", forest = "q")),
  "have no labels",
  fixed = TRUE
)

## --- Block C: a fit ---------------------------------------------------------
fitted <- bart(
  y ~ forest(x1, basis = dose) +
    forest(x2, basis = age) +
    forest(x3, basis = scale(age)),
  frame,
  n.chains = 1L,
  n.threads = 1L,
  n.samples = 6L,
  n.burn = 2L,
  keepTrees = TRUE,
  verbose = FALSE
)
fitLabels <- attr(fitted, "forest.labels")
expect_identical(fitLabels, c("dose", "age", "scale(age)"))
# the margins keep their names
expect_identical(
  dimnames(fitted$forestFits)[[3L]],
  c("forest1", "forest2", "forest3")
)
expect_identical(
  names(extract(fitted, "k", forest = NULL)),
  c("forest1", "forest2", "forest3")
)
newRows <- frame[1:7, ]
for (index in 1:3) {
  position <- extract(fitted, "forest", forest = index)
  expect_identical(
    extract(fitted, "forest", forest = fitLabels[index]),
    position
  )
  expect_identical(
    extract(fitted, "forest", forest = paste0("forest", index)),
    position
  )
  expect_identical(
    extract(fitted, "forest", forest = fitLabels[index], contribution = TRUE),
    extract(fitted, "forest", forest = index, contribution = TRUE)
  )
  for (type in c("k", "leaf.prior.sd")) {
    expect_identical(
      extract(fitted, type, forest = fitLabels[index]),
      extract(fitted, type, forest = index)
    )
    expect_identical(
      extract(fitted, type, forest = paste0("forest", index)),
      extract(fitted, type, forest = index)
    )
  }
  expect_identical(
    extract(fitted, "trees", forest = fitLabels[index]),
    extract(fitted, "trees", forest = index)
  )
  expect_identical(
    predict(fitted, newRows, type = "forest", forest = fitLabels[index]),
    predict(fitted, newRows, type = "forest", forest = index)
  )
  expect_identical(
    predict(fitted, newRows, type = "forest", forest = paste0("forest", index)),
    predict(fitted, newRows, type = "forest", forest = index)
  )
}
# the leaf prior's sd differs by forest, so the identity says something
expect_false(
  identical(
    extract(fitted, "leaf.prior.sd", forest = "dose"),
    extract(fitted, "leaf.prior.sd", forest = "scale(age)")
  )
)
# a vector of labels comes back in the order given, the margin as before
byLabel <- extract(fitted, "forest", forest = c("scale(age)", "dose"))
byPosition <- extract(fitted, "forest", forest = c(3L, 1L))
expect_identical(byLabel, byPosition)
expect_identical(dimnames(byLabel)[[3L]], c("forest3", "forest1"))
expect_identical(
  extract(fitted, "k", forest = c("scale(age)", "dose")),
  extract(fitted, "k", forest = c(3L, 1L))
)
expect_identical(
  predict(fitted, newRows, type = "forest", forest = c("age", "dose")),
  predict(fitted, newRows, type = "forest", forest = c(2L, 1L))
)
# the arms that refuse 'forest' keep their texts, and the others refuse by
# kind and by name
expect_error(extract(fitted, "sigma", forest = 1L), "model parameter")
expect_error(
  extract(fitted, "k", forest = 4L),
  "'forest' index must be between 1 and 3",
  fixed = TRUE
)
expect_error(
  extract(fitted, "forest", forest = "x3"),
  paste0(
    "'forest' names no forest of this model: \"x3\"; its forests are ",
    "\"dose\", \"age\", \"scale(age)\""
  ),
  fixed = TRUE
)
expect_error(
  predict(fitted, newRows, type = "forest", forest = "forest4"),
  "'forest' names no forest of this model: \"forest4\"",
  fixed = TRUE
)
expect_error(
  extract(fitted, "leaf.prior.sd", forest = "2"),
  "a position is given as a number, forest = 2",
  fixed = TRUE
)
expect_error(
  extract(fitted, "forest", forest = factor("dose")),
  "not a factor",
  fixed = TRUE
)
expect_error(
  extract(fitted, "forest", forest = TRUE),
  "not a logical",
  fixed = TRUE
)
expect_error(
  extract(fitted, "k", forest = list(2)),
  "not a list",
  fixed = TRUE
)
# plotTree hands its forest to the sampler
pdf(NULL)
expect_silent(plotTree(fitted, 1L, chainNum = 1L, forest = "scale(age)"))
dev.off()

## --- Block D: lists given forest by forest ----------------------------------
# a name on an entry of 'bases' is its position's label or forest<i>
atRows <- list(newRows$dose, newRows$age, newRows$age)
plain <- predict(fitted, newRows, bases = atRows, n.threads = 1L)
labelled <- atRows
names(labelled) <- fitLabels
expect_identical(
  predict(fitted, newRows, bases = labelled, n.threads = 1L),
  plain
)
positioned <- atRows
names(positioned) <- c("forest1", "forest2", "forest3")
expect_identical(
  predict(fitted, newRows, bases = positioned, n.threads = 1L),
  plain
)
mixed <- atRows
names(mixed) <- c("dose", "", "forest3")
expect_identical(predict(fitted, newRows, bases = mixed, n.threads = 1L), plain)
exchanged <- atRows
names(exchanged) <- c("age", "dose", "scale(age)")
expect_error(
  predict(fitted, newRows, bases = exchanged, n.threads = 1L),
  paste0(
    "'bases' names forest 1 \"age\", and its label is \"dose\"; give the ",
    "bases in the forests' order"
  ),
  fixed = TRUE
)
unknown <- atRows
names(unknown) <- c("dose", "age", "other")
expect_error(
  predict(fitted, newRows, bases = unknown, n.threads = 1L),
  "'bases' names forest 3 \"other\", and its label is \"scale(age)\"",
  fixed = TRUE
)

# $setLeafPrior(forests = ) reads its list by position, in label order
inOrder <- list(a = forest(sd = 1.5), forest(), `scale(age)` = forest(sd = 2.5))
expect_silent(sampler$setLeafPrior(forests = inOrder))
expect_error(
  sampler$setLeafPrior(
    forests = list(forest(), a = forest(sd = 1.5), forest())
  ),
  "names forest 2 'a', but it was created as 'forest2'",
  fixed = TRUE
)
expect_silent(sampler$setLeafPrior(
  forests = list(forest1 = forest(sd = 1.5), forest2 = forest())
))
expect_error(
  sampler$setLeafPrior(forests = list(forest2 = forest(sd = 1.5))),
  "names forest 1 'forest2', but it was created as 'a'",
  fixed = TRUE
)

# a name that is another forest's label is not read as this position's
# forest<i>: forests labelled forest2, dose, forest1
taken <- frame
taken$forest1 <- frame$age
taken$forest2 <- frame$dose
takenSampler <- dbarts(
  y ~ x1 + x2 + x3,
  taken,
  forests = list(
    forest(basis = forest2, sd = 1),
    forest(basis = dose, sd = 2),
    forest(basis = forest1, sd = 3)
  ),
  control = selectionControl()
)
expect_identical(
  attr(takenSampler$control, "bartcore.forests", exact = TRUE)$labels,
  c("forest2", "dose", "forest1")
)
expect_identical(
  refusal(takenSampler$getLeafPrior("forest1")),
  paste0(
    "'forest' (\"forest1\") is the label of forest 3 and the name of ",
    "position 1; select by position, as forest = 3"
  )
)
expect_identical(
  refusal(takenSampler$setLeafPrior(forests = list(forest1 = forest(sd = 9)))),
  paste0(
    "'forest' (\"forest1\") is the label of forest 3 and the name of ",
    "position 1; select by position, as forest = 3"
  )
)
expect_identical(
  vapply(1:3, function(i) takenSampler$getLeafPrior(i)$leaf.prior$sd, 0),
  c(1, 2, 3)
)
# its own label at its own position is still accepted
expect_silent(takenSampler$setLeafPrior(
  forests = list(forest2 = forest(sd = 9), dose = forest(sd = 8))
))
expect_identical(
  vapply(1:3, function(i) takenSampler$getLeafPrior(i)$leaf.prior$sd, 0),
  c(9, 8, 3)
)
takenFit <- bart(
  y ~ forest(x1, basis = forest2) + forest(x2, basis = forest1),
  taken,
  n.chains = 1L,
  n.threads = 1L,
  n.samples = 6L,
  n.burn = 2L,
  keepTrees = TRUE,
  verbose = FALSE
)
expect_identical(attr(takenFit, "forest.labels"), c("forest2", "forest1"))
takenRows <- taken[1:7, ]
inOrder <- list(takenRows$forest2, takenRows$forest1)
expect_silent(predict(takenFit, takenRows, bases = inOrder, n.threads = 1L))
named <- inOrder
names(named) <- c("forest2", "forest1")
expect_identical(
  predict(takenFit, takenRows, bases = named, n.threads = 1L),
  predict(takenFit, takenRows, bases = inOrder, n.threads = 1L)
)
swapped <- inOrder
names(swapped) <- c("forest1", "forest2")
expect_identical(
  refusal(predict(takenFit, takenRows, bases = swapped, n.threads = 1L)),
  paste0(
    "'forest' (\"forest1\") is the label of forest 2 and the name of ",
    "position 1; select by position, as forest = 2"
  )
)
