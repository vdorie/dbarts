# A forest's basis is the right-hand side of a model formula with no tilde,
# and means what it means in lm(): '+' separates columns, I() holds
# arithmetic, a factor gives a column for each level, and scale() and poly()
# are rebuilt at new rows from the fitted rows. The reference throughout is
# R's own model frame and model matrix of ~ 0 + <basis>. Block A: the
# grammar, and what it refuses before anything is evaluated. Block B: the
# columns and their names. Block C: the texts whose meaning changed. Block D:
# rows, at both doors. Block E: new rows. Block F: 'subset' is read once, by
# the data object, and every basis is cut by the rows it kept. Block G: 'data'
# is read once. Block H: a value that the kept rows leave with one level.

forest <- dbartsForests$forest

set.seed(23)
n <- 80L
frame <- data.frame(
  x1 = runif(n),
  x2 = runif(n),
  x3 = runif(n),
  dose = runif(n, 0.5, 2),
  age = 50 + 10 * rnorm(n),
  z = rep_len(c(0L, 1L), n),
  g = sample(c("u", "v", "w"), n, replace = TRUE),
  stringsAsFactors = FALSE
)
frame$zl <- frame$z == 1L
frame$gf <- factor(frame$g, levels = c("u", "v", "w", "never"))
frame$W <- cbind(p = frame$dose, q = frame$age)
frame$y <- with(frame, x1 + z * (1 + x2) + 0.2 * dose + rnorm(n, sd = 0.3))
newRows <- data.frame(
  x1 = runif(9L),
  x2 = runif(9L),
  x3 = runif(9L),
  dose = runif(9L, 0.5, 2),
  age = unname(stats::quantile(frame$age, seq(0.1, 0.9, length.out = 9L))),
  z = rep_len(c(1L, 0L), 9L),
  g = rep_len(c("v", "u"), 9L),
  stringsAsFactors = FALSE
)
newRows$zl <- newRows$z == 1L
newRows$gf <- factor(newRows$g, levels = levels(frame$gf))
newRows$W <- cbind(p = newRows$dose, q = newRows$age)
k <- 30

basisControl <- function() {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 5L,
    n.samples = 3L,
    n.burn = 1L,
    updateState = FALSE,
    verbose = FALSE,
    seed = 23L
  )
}
asCode <- function(basis) if (is.character(basis)) str2lang(basis) else basis
# one model at the two doors: a term of the formula, and a 'forests' list
termFit <- function(basis, data = frame, ...) {
  formula <- eval(bquote(y ~ x1 + x2 + forest(x1, basis = .(asCode(basis)))))
  environment(formula) <- parent.frame()
  call <- quote(dbarts(formula, data, control = basisControl()))
  eval(as.call(c(as.list(call), list(...))))
}
listFit <- function(basis, data = frame, ...) {
  call <- bquote(dbarts(
    y ~ x1 + x2,
    data,
    forests = list(forest(), forest(x1, basis = .(asCode(basis)))),
    control = basisControl()
  ))
  eval(as.call(c(as.list(call), list(...))), list(data = data), parent.frame())
}
doors <- list(term = termFit, list = listFit)
# the same basis handed over as a value, which reads no code
valueFit <- function(value, data = frame, ...) {
  forests <- list(forest(), do.call(forest, list("x1", basis = value)))
  dbarts(y ~ x1 + x2, data, forests = forests, control = basisControl(), ...)
}
bartFit <- function(basis, data = frame, ...) {
  formula <- eval(bquote(y ~ x1 + x2 + forest(x1, basis = .(asCode(basis)))))
  environment(formula) <- parent.frame()
  call <- quote(bart(
    formula,
    data,
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 5L,
    n.samples = 3L,
    n.burn = 0L,
    keepTrees = TRUE,
    verbose = FALSE,
    seed = 23L
  ))
  eval(as.call(c(as.list(call), list(...))))
}
# a sampler as the fit bart() would return of it, for predict
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
basisOf <- function(sampler) sampler$data@bases[[2L]]
draws <- function(sampler) {
  run <- sampler$run(0L, 10L)
  list(run$train, run$sigma, sampler$getForestAmplitudes())
}
refusal <- function(expr) {
  tryCatch(
    {
      force(expr)
      NA_character_
    },
    error = function(e) conditionMessage(e)
  )
}
# what lm() makes of ~ 0 + <basis> on `rows` of `data`: its model matrix,
# built on the rows after every term has been evaluated on all of them
lmFormula <- function(basis) {
  formula <- call("~", call("+", 0, asCode(basis)))
  class(formula) <- "formula"
  environment(formula) <- environment(lmFormula)
  formula
}
lmColumns <- function(basis, data = frame, rows = NULL) {
  modelFrame <- stats::model.frame(
    lmFormula(basis),
    data,
    na.action = stats::na.pass
  )
  if (!is.null(rows)) {
    terms <- attr(modelFrame, "terms")
    modelFrame <- modelFrame[rows, , drop = FALSE]
    attr(modelFrame, "terms") <- terms
  }
  design <- stats::model.matrix(attr(modelFrame, "terms"), modelFrame)
  matrix(
    as.double(design),
    nrow(design),
    dimnames = list(NULL, colnames(design))
  )
}
lmNames <- function(basis, data = frame) {
  formula <- lmFormula(basis)
  formula <- eval(call("~", quote(y), formula[[2L]]))
  environment(formula) <- environment(lmFormula)
  names(stats::coef(stats::lm(formula, data)))
}

## --- Block A: the grammar ----------------------------------------------------
# each refusal with its text, the same at both doors and decided from the
# code alone: with no data at hand, where forest() is called
refused <- list(
  "dose * age" = paste0(
    "'basis' does not take '*' between its terms ('dose * age'): in a model ",
    "formula it is both columns and their product. Write I(dose * age) for ",
    "the product alone, or dose + age + I(dose * age) for all three"
  ),
  "poly(dose, 2) * age" = paste0(
    "'basis' does not take '*' between its terms ('poly(dose, 2) * age')"
  ),
  "dose/30" = paste0(
    "'basis' term 'dose/30' divides a column by a number, which a model ",
    "formula does not take; write I(dose/30) for the rescaled column"
  ),
  "30 * dose" = paste0(
    "'basis' term '30 * dose' multiplies a column by a number, which a model ",
    "formula does not take; write I(30 * dose) for the rescaled column"
  ),
  "dose * 30" = paste0(
    "'basis' term 'dose * 30' multiplies a column by a number, which a model ",
    "formula does not take; write I(dose * 30) for the rescaled column"
  ),
  "dose - age" = paste0(
    "'basis' does not take '-' between its terms ('dose - age'): in a model ",
    "formula it removes a term and subtracts nothing; write the arithmetic ",
    "inside I(), as I(dose - age)"
  ),
  "dose^2" = paste0(
    "'basis' does not take '^' between its terms ('dose^2'): in a model ",
    "formula it crosses terms and raises nothing to a power; write the ",
    "arithmetic inside I(), as I(dose^2)"
  ),
  "dose/age" = paste0(
    "'basis' does not take '/' between its terms ('dose/age'): in a model ",
    "formula it nests one term in another and divides nothing; write the ",
    "arithmetic inside I(), as I(dose/age)"
  ),
  "dose %in% age" = paste0(
    "'basis' does not take '%in%' ('dose %in% age'): in a model formula it ",
    "nests one term in another"
  ),
  "dose | z" = "'basis' does not take '|' ('dose | z')",
  "dose + ." = paste0(
    "'basis' does not take '.': name the columns the forest is multiplied by"
  ),
  "dose + offset(age)" = "'basis' does not take an offset() term",
  "dose + 2" = paste0(
    "'basis' has the number 2 as a term; only 1, a constant column, and 0, ",
    "none, are terms. Write arithmetic on a column inside I()"
  ),
  "dose + TRUE" = paste0(
    "'basis' has the constant TRUE as a term; a term is a column, and only ",
    "1, a constant column, and 0, none, are written as numbers"
  ),
  "0 + 1" = paste0(
    "'basis' (0 + 1) is a constant column and nothing else, which multiplies ",
    "the forest by a constant: leave 'basis' out for the forest with no ",
    "multiplier"
  ),
  "-1" = paste0(
    "'basis' (-1) names no column; leave 'basis' out for the forest with no ",
    "multiplier"
  ),
  "cbind(dose, age)" = paste0(
    "'basis' term 'cbind(dose, age)': the columns of a basis are separated ",
    "by '+'; write dose + age"
  ),
  "cbind(1 - z, z)" = paste0(
    "'basis' term 'cbind(1 - z, z)': the columns of a basis are separated ",
    "by '+'; write I(1 - z) + z"
  ),
  "cbind(1, dose)" = paste0(
    "'basis' term 'cbind(1, dose)': the columns of a basis are separated by ",
    "'+'; write 1 + dose"
  ),
  "normal(dose, sd = 30)" = paste0(
    "'basis' term 'normal(dose, sd = 30)' calls normal(), which is not a ",
    "column: a prior is not stated on a basis term. Give the forest's size ",
    "as 'sd' and its coefficient's law as 'amplitude'"
  )
)
for (text in names(refused)) {
  code <- str2lang(text)
  # forest() itself, which sees no data
  expect_true(
    startsWith(
      refusal(eval(bquote(forest(basis = .(code))))),
      refused[[text]]
    ),
    info = text
  )
  for (door in names(doors)) {
    said <- refusal(doors[[door]](code))
    expect_true(
      startsWith(said, refused[[text]]),
      info = paste(door, text, said)
    )
  }
}
# what the refusal of cbind() says to write names cbind()'s columns in
# cbind()'s order, and where no sum of terms does it says none: a model
# formula keeps one of two terms alike and puts its constant column first
for (text in c("cbind(dose, dose)", "cbind(dose, 1)", "cbind(dose, 2)")) {
  for (door in doors) {
    expect_identical(
      refusal(door(text)),
      paste0(
        "'basis' term '",
        text,
        "': the columns of a basis are separated by '+'"
      ),
      info = text
    )
  }
}
for (written in list(
  c("cbind(dose, age)", "dose + age"),
  c("cbind(1 - z, z)", "I(1 - z) + z"),
  c("cbind(1, dose)", "1 + dose"),
  c("cbind(log(dose), age * 2)", "log(dose) + I(age * 2)")
)) {
  expect_identical(
    unname(basisOf(listFit(written[2L]))),
    unname(eval(str2lang(written[1L]), frame)),
    info = written[1L]
  )
}
# a term that calls a name a prior could later be stated with is refused,
# whoever defines the function; any other function is a column's
reserved <- c(
  "normal",
  "fixed",
  "student",
  "cauchy",
  "linear",
  "gp",
  "cgm",
  "dart",
  "chisq",
  "chi",
  "invchi",
  "forest",
  "varianceForest"
)
expect_identical(sort(dbarts:::BASIS_RESERVED_CALLS), sort(reserved))
normal <- function(x) x / 2
student <- function(x) x
weights2 <- function(x) x * 2
for (name in reserved) {
  code <- as.call(list(as.name(name), quote(dose)))
  for (door in names(doors)) {
    expect_true(
      startsWith(
        refusal(doors[[door]](code)),
        paste0(
          "'basis' term '",
          name,
          "(dose)' calls ",
          name,
          "(), which is not a column"
        )
      ),
      info = paste(door, name)
    )
  }
}
expect_error(
  listFit("dose + normal(age)"),
  "'basis' term 'normal(age)' calls normal(), which is not a column",
  fixed = TRUE
)
for (door in doors) {
  expect_identical(
    basisOf(door("weights2(dose)")),
    lmColumns("weights2(dose)")
  )
}
rm(normal, student)
# a column named like one is a column
named <- frame
named$normal <- frame$dose
named$fixed <- frame$age
named$forest <- frame$dose * 2
for (door in doors) {
  expect_identical(
    basisOf(door("normal + fixed + forest", named)),
    lmColumns("normal + fixed + forest", named)
  )
}
# cbind() called through its package is a matrix term, as in lm()
for (door in doors) {
  expect_identical(
    basisOf(door("base::cbind(dose, age)")),
    lmColumns("base::cbind(dose, age)")
  )
}
expect_identical(
  colnames(basisOf(listFit("base::cbind(dose, age)"))),
  c("base::cbind(dose, age)dose", "base::cbind(dose, age)age")
)

# every shape of operator at the top of a basis is either refused or is the
# columns lm() gives, with one verdict at the two doors
shapes <- c(
  "dose %% 2" = TRUE,
  "dose %/% 30" = FALSE,
  "dose > 1" = TRUE,
  "-dose" = FALSE,
  "+dose" = TRUE,
  "dose - 0" = FALSE,
  "dose + dose" = TRUE,
  "dose:age:x3" = TRUE,
  "(dose + age):x3" = TRUE,
  "dose * age * x3" = FALSE,
  "dose * k" = FALSE,
  "dose:k" = FALSE,
  "dose / k" = FALSE,
  "dose * zl" = FALSE,
  "dose:zl" = FALSE,
  "dose:dose" = TRUE,
  "I(dose) * 2" = FALSE,
  "dose / 30L" = FALSE,
  "dose * 1" = FALSE,
  "dose:1" = FALSE,
  "dose + 1 - 1" = TRUE,
  "-1 + dose" = TRUE,
  "(dose)" = TRUE,
  "dose + -1" = TRUE,
  "dose ** 2" = FALSE,
  "dose + age - age" = FALSE,
  "log(dose) - 1" = TRUE,
  "dose %in% age" = FALSE,
  "(dose + age)^2" = FALSE,
  "dose + (1 | g)" = FALSE,
  "!zl" = TRUE,
  "zl & (x1 > 0.5)" = TRUE,
  "dose == 1" = FALSE,
  "z == 1" = TRUE,
  "ifelse(z == 1, dose, 0)" = TRUE,
  "dose[]" = TRUE,
  "frame$dose" = TRUE,
  "`dose`" = TRUE,
  "dose + NULL" = FALSE,
  "c(dose)" = TRUE,
  "identity(cbind(dose, age))" = TRUE,
  "base::cbind(dose, age)" = TRUE,
  "dbarts::normal(dose)" = FALSE,
  "stats::poly(dose, 2)" = TRUE,
  "(normal)(dose)" = FALSE,
  "scale(normal(dose))" = FALSE,
  "dose %o% 1" = TRUE,
  "dose * (age + x3)" = FALSE,
  "dose %*% 1" = FALSE,
  "1:dose" = FALSE,
  "dose + I(1)" = FALSE,
  "TRUE" = FALSE,
  "dose + TRUE" = FALSE,
  "NULL + dose" = FALSE,
  "dose - age:dose" = FALSE,
  "1 - dose" = FALSE,
  "0 + 1" = FALSE
)
expect_identical(length(shapes), 57L)
for (text in names(shapes)) {
  code <- str2lang(text)
  read <- lapply(doors, function(door) {
    tryCatch(basisOf(door(code)), error = function(e) conditionMessage(e))
  })
  expect_identical(read$term, read$list, info = text)
  expect_identical(is.matrix(read$list), unname(shapes[[text]]), info = text)
  if (is.matrix(read$list)) {
    expect_identical(read$list, lmColumns(code), info = text)
  }
}

## --- Block B: the columns and their names -----------------------------------
# as coef(lm(y ~ 0 + <basis>)) names them, and with lm()'s values
fitted <- c(
  "dose",
  "dose + age",
  "I(dose + age)",
  "I(dose/30)",
  "I(dose / k)",
  "log(dose)",
  "scale(age)",
  "poly(dose, 2)",
  "scale(age) + poly(dose, 2)",
  "dose:age",
  "dose + age + I(dose * age)",
  "1 + dose",
  "factor(z)",
  "g",
  "zl",
  "W",
  "dose + W",
  "I(age - mean(age))"
)
expect_identical(length(fitted), 18L)
for (text in fitted) {
  expected <- lmColumns(text)
  expect_identical(colnames(expected), lmNames(text), info = text)
  handedOver <- draws(valueFit(
    if (text %in% c("factor(z)", "g", "zl")) {
      eval(str2lang(text), frame)
    } else {
      unname(expected)
    }
  ))
  for (door in names(doors)) {
    fit <- doors[[door]](text)
    expect_identical(basisOf(fit), expected, info = paste(door, text))
    # the model is the one the same columns handed over as a value fit
    expect_identical(draws(fit), handedOver, info = paste(door, text))
  }
  expect_identical(bartFit(text)$bases[[2L]], expected, info = text)
}
expect_identical(lmNames("I(dose/30)"), "I(dose/30)")
expect_identical(lmNames("1 + dose"), c("(Intercept)", "dose"))
expect_identical(
  lmNames("poly(dose, 2)"),
  c("poly(dose, 2)1", "poly(dose, 2)2")
)
expect_identical(lmNames("factor(z)"), c("factor(z)0", "factor(z)1"))
expect_identical(lmNames("W"), c("Wp", "Wq"))
# a tilde is taken off and changes nothing
for (text in c("dose + age", "factor(z)", "scale(age)")) {
  for (door in doors) {
    expect_identical(
      basisOf(door(paste("~", text))),
      lmColumns(text),
      info = text
    )
  }
}
# a factor's columns are every level's, by no contrast: no level is dropped
expect_identical(ncol(basisOf(listFit("g"))), 3L)
expect_identical(
  unname(rowSums(basisOf(listFit("g")))),
  rep(1, n)
)
# one kind in a basis: a factor beside other terms is refused, as are a
# logical matrix and two columns of one name
mixes <- paste0(
  ") mixes a factor with other terms: a basis is one factor, a character ",
  "or a logical vector, with one coefficient per level, or numeric columns, ",
  "with one each. Give the other terms a forest() of their own, or write ",
  "the factor as numbers"
)
for (text in c("dose + factor(z)", "g:dose", "1 + factor(z)", "g + zl")) {
  for (door in doors) {
    expect_identical(
      refusal(door(text)),
      paste0("'basis' (", text, mixes),
      info = text
    )
  }
}
frame$L <- cbind(frame$zl, !frame$zl)
frame$A <- cbind(a = frame$dose, a = frame$age)
for (door in doors) {
  expect_identical(
    refusal(door("L")),
    paste0(
      "'basis' (L) has a term that is a logical matrix: a basis is a ",
      "factor, a character or logical vector, or numeric"
    )
  )
  expect_identical(
    refusal(door("A")),
    paste0(
      "'basis' (A) has two columns named \"Aa\"; the columns of a basis are ",
      "told apart by name"
    )
  )
}
frame$L <- NULL
frame$A <- NULL

# a held coefficient on a basis of one numeric column is refused however the
# column is written, and held on two columns at 0 and 1
heldOneColumn <- paste0(
  "forest 2: amplitude = fixed() on a basis of one numeric column is not ",
  "supported yet"
)
for (text in c("dose", "I(dose / 30)", "scale(age)", "dose:age", "~ dose")) {
  code <- str2lang(text)
  expect_error(
    eval(bquote(dbarts(
      y ~ x1 + x2 + forest(x1, basis = .(code), amplitude = fixed()),
      frame,
      control = basisControl()
    ))),
    heldOneColumn,
    fixed = TRUE,
    info = text
  )
  expect_error(
    eval(bquote(dbarts(
      y ~ x1 + x2,
      frame,
      forests = list(
        forest(),
        forest(x1, basis = .(code), amplitude = fixed())
      ),
      control = basisControl()
    ))),
    heldOneColumn,
    fixed = TRUE,
    info = text
  )
}
heldTwo <- dbarts(
  y ~ x1 + x2 + forest(x1, basis = dose + age, amplitude = fixed()),
  frame,
  control = basisControl()
)
invisible(heldTwo$run(0L, 5L))
expect_equal(heldTwo$getForestAmplitudes()[2:3, 1L], c(0, 1))

## --- Block C: the texts whose meaning changed --------------------------------
# '+' is columns; the sum is written inside I()
expect_identical(
  draws(listFit("~ dose + age")),
  draws(valueFit(cbind(frame$dose, frame$age)))
)
expect_identical(
  draws(listFit("I(dose + age)")),
  draws(valueFit(frame$dose + frame$age))
)
expect_false(identical(
  draws(listFit("dose + age")),
  draws(listFit("I(dose + age)"))
))
# '- 1' removes the constant column a basis does not have; the arithmetic is
# written inside I()
expect_identical(basisOf(listFit("~ dose - 1")), lmColumns("dose"))
expect_identical(draws(listFit("~ dose - 1")), draws(valueFit(frame$dose)))
expect_identical(
  draws(listFit("I(dose - 1)")),
  draws(valueFit(frame$dose - 1))
)
# '1 +' asks for the constant column
expect_identical(
  colnames(basisOf(listFit("~ 1 + dose"))),
  c("(Intercept)", "dose")
)
expect_identical(
  draws(listFit("~ 1 + dose")),
  draws(valueFit(cbind(1, frame$dose)))
)
expect_identical(
  draws(listFit("I(1 + dose)")),
  draws(valueFit(1 + frame$dose))
)
# a name is the data's column at both doors, whatever is beside the call
local({
  dose <- frame$dose * 100
  fromData <- draws(valueFit(frame$dose))
  expect_identical(draws(termFit("dose")), fromData)
  expect_identical(draws(listFit("dose")), fromData)
})
# scale() and poly() under 'subset' take their centre from every row, as lm()
# does, and so equal the constants written out
keep <- frame$age < 56
centre <- mean(frame$age)
spread <- stats::sd(frame$age)
for (door in doors) {
  fit <- door("scale(age)", subset = keep)
  expect_identical(basisOf(fit), lmColumns("scale(age)", rows = keep))
  expect_equal(
    as.vector(basisOf(fit)),
    (frame$age[keep] - centre) / spread
  )
  expect_identical(
    basisOf(door("poly(dose, 2)", subset = keep)),
    lmColumns("poly(dose, 2)", rows = keep)
  )
}
expect_false(isTRUE(all.equal(
  as.vector(basisOf(termFit("scale(age)", subset = keep))),
  as.vector(scale(frame$age[keep]))
)))

## --- Block D: rows, at both doors -------------------------------------------
# one text is one basis and one model at the two doors, under 'subset' and
# under a response whose missing rows the na.action drops
noW <- frame$g != "w"
missingY <- frame
missingY$y[frame$g == "w"] <- NA
missingY$y[3:5] <- NA
keptY <- !is.na(missingY$y)
for (text in c("factor(g)", "g", "gf", "scale(age)", "dose + age", "zl")) {
  underSubset <- lapply(doors, function(door) door(text, subset = noW))
  expect_identical(
    basisOf(underSubset$term),
    basisOf(underSubset$list),
    info = text
  )
  expect_identical(
    basisOf(underSubset$list),
    lmColumns(text, rows = noW)[,
      colSums(lmColumns(text, rows = noW) != 0) > 0L,
      drop = FALSE
    ],
    info = text
  )
  expect_identical(
    draws(underSubset$term),
    draws(underSubset$list),
    info = text
  )
  underMissing <- lapply(doors, function(door) door(text, missingY))
  expect_identical(nrow(basisOf(underMissing$list)), sum(keptY), info = text)
  expect_identical(
    basisOf(underMissing$term),
    basisOf(underMissing$list),
    info = text
  )
  expect_identical(
    draws(underMissing$term),
    draws(underMissing$list),
    info = text
  )
}
# a level the kept rows leave empty is no column, whatever emptied it: never
# used, cut by 'subset', or gone with the rows of a missing response
for (door in doors) {
  expect_identical(colnames(basisOf(door("gf"))), c("gfu", "gfv", "gfw"))
  expect_identical(
    colnames(basisOf(door("gf", subset = noW))),
    c("gfu", "gfv")
  )
  expect_identical(colnames(basisOf(door("g", missingY))), c("gu", "gv"))
  # the rows are the kept rows, in their order
  expect_identical(
    unname(basisOf(door("g", subset = noW))[, 1L]),
    as.double(frame$g[noW] == "u")
  )
  expect_identical(
    as.vector(basisOf(door("dose", missingY))),
    frame$dose[keptY]
  )
  # a factor left with one level, and a logical left with one value
  expect_identical(
    refusal(door("g", subset = frame$g == "u")),
    "a 'basis' factor must have at least two levels"
  )
  expect_identical(
    refusal(door("zl", subset = frame$zl)),
    paste0(
      "'basis' (zl) is TRUE on every row the fit keeps, so its other level ",
      "has no observations and the forest would be multiplied by a constant"
    )
  )
  expect_identical(
    refusal(door("dose > 0")),
    paste0(
      "'basis' (dose > 0) is TRUE on every row the fit keeps, so its other ",
      "level has no observations and the forest would be multiplied by a ",
      "constant"
    )
  )
  expect_identical(
    refusal(door("z", subset = frame$z == 0L)),
    "a 'basis' column of all zeros contributes nothing to a forest"
  )
}
# a missing value in a row the fit keeps is refused; in a row 'subset' or the
# na.action drops it is no value of the basis
missingDose <- frame
missingDose$dose[c(2L, 7L)] <- NA
unknownDose <- is.na(missingDose$dose)
bothMissing <- missingDose
bothMissing$y[unknownDose] <- NA
for (door in doors) {
  expect_identical(
    refusal(door("dose", missingDose)),
    "a 'basis' cannot be NA"
  )
  expect_identical(
    as.vector(basisOf(door("dose", missingDose, subset = !unknownDose))),
    frame$dose[!unknownDose]
  )
  expect_identical(
    as.vector(basisOf(door("dose", bothMissing))),
    frame$dose[!unknownDose]
  )
}
# a basis covers every row of the data and is cut with it: one already cut is
# refused, a value among it
shortBasis <- frame$dose[noW]
for (door in doors) {
  expect_identical(
    refusal(door("shortBasis", subset = noW)),
    paste0(
      "'basis' (shortBasis) must have the same length as the data: it has ",
      sum(noW),
      " rows and the data 80; a basis covers every row of the data and is ",
      "cut by 'subset' and the na.action with it"
    )
  )
}
fullBasis <- frame$dose * 2
for (door in doors) {
  expect_identical(
    as.vector(basisOf(door("fullBasis", subset = noW))),
    fullBasis[noW]
  )
}
# 'subset' is read as the fit's own model frame reads it: row names select
# by name, and one a wrapper forwards through its dots is the wrapper's
lettered <- frame
row.names(lettered) <- paste0("r", seq_len(n))
namedRows <- row.names(lettered)[noW]
throughDots <- function(...) {
  dbarts(
    y ~ x1 + x2 + forest(x1, basis = g),
    frame,
    control = basisControl(),
    ...
  )
}
expect_identical(
  basisOf(throughDots(subset = noW)),
  basisOf(termFit("g", subset = noW))
)
for (door in doors) {
  expect_identical(
    basisOf(door("g", lettered, subset = namedRows)),
    basisOf(door("g", subset = noW))
  )
  expect_identical(
    basisOf(door("g", subset = which(noW))),
    basisOf(door("g", subset = noW))
  )
  # rows repeated and reordered are the fit's rows, in its order
  expect_identical(
    as.vector(basisOf(door("dose", subset = c(5L, 3L, 3L, 60L, 1L)))),
    frame$dose[c(5L, 3L, 3L, 60L, 1L)]
  )
}
# with the matrix interface too, where 'subset' is the rows themselves
xMatrix <- as.matrix(frame[c("x1", "x2")])
yVector <- frame$y
gVector <- frame$g
matrixFit <- dbarts(
  xMatrix,
  yVector,
  subset = which(noW),
  forests = list(forest(), forest(x1, basis = gVector)),
  control = basisControl()
)
expect_identical(
  basisOf(matrixFit),
  lmColumns("gVector", rows = noW)[, 1:2, drop = FALSE]
)
# and the rows its na.action drops go from the basis as from the response
doseVector <- frame$dose
yMissing <- frame$y
yMissing[c(4L, 11L)] <- NA
matrixMissing <- dbarts(
  xMatrix,
  yMissing,
  subset = 3:60,
  forests = list(forest(), forest(x1, basis = doseVector)),
  control = basisControl()
)
expect_identical(length(matrixMissing$data@y), 56L)
expect_identical(
  as.vector(basisOf(matrixMissing)),
  doseVector[setdiff(3:60, c(4L, 11L))]
)

## --- Block E: new rows -------------------------------------------------------
# the stored record is R's terms, so that each term is rebuilt at new rows
# from what the fitted rows gave it; at both doors
replayed <- function(fit, rows) {
  dbarts:::replayForestBasis(fit$basis.terms[[2L]], rows, 2L)
}
atNewRows <- function(basis, rows) {
  if (identical(basis, "zl")) {
    # lm() would give one row's logical a single column
    return(cbind(zlFALSE = as.double(!rows$zl), zlTRUE = as.double(rows$zl)))
  }
  modelFrame <- stats::model.frame(lmFormula(basis), frame)
  terms <- attr(modelFrame, "terms")
  again <- stats::model.frame(
    terms,
    rows,
    xlev = stats::.getXlevels(terms, modelFrame)
  )
  design <- stats::model.matrix(terms, again)
  matrix(
    as.double(design),
    nrow(design),
    dimnames = list(NULL, colnames(design))
  )
}
hasSplines <- requireNamespace("splines", quietly = TRUE)
rebuilt <- c(
  "scale(age)",
  "scale(age, scale = FALSE)",
  "poly(dose, 2)",
  "poly(dose, 2, raw = TRUE)",
  "scale(age) + poly(dose, 2)",
  "scale(log(age))",
  "I(age - mean(age))",
  "1 + dose",
  "dose:age",
  "factor(z)",
  "g",
  "zl"
)
if (hasSplines) {
  rebuilt <- c(rebuilt, "splines::ns(age, 3)", "splines::bs(age, df = 4)")
}
plainColumns <- c("x1", "x2", "x3", "dose", "age", "z", "g", "zl")
stacked <- rbind(frame[1:6, plainColumns], newRows[plainColumns])
for (text in rebuilt) {
  fits <- list(
    term = bartFit(text),
    list = packaged(listFit(text))
  )
  for (door in names(fits)) {
    fit <- fits[[door]]
    info <- paste(door, text)
    expect_identical(
      replayed(fit, newRows),
      atNewRows(text, newRows),
      info = info
    )
    expect_identical(
      replayed(fit, newRows[2L, ]),
      atNewRows(text, newRows[2L, ]),
      info = info
    )
    expect_identical(
      replayed(fit, stacked),
      atNewRows(text, stacked),
      info = info
    )
    # at the fitted rows it is the basis the fit used
    expect_equal(replayed(fit, frame), fit$bases[[2L]], info = info)
    # and predict uses it: the same as the basis rebuilt by hand from the new
    # rows and handed over
    expect_identical(
      predict(fit, newRows),
      predict(fit, newRows, bases = atNewRows(text, newRows)),
      info = info
    )
  }
}
# the constants are the fitted rows': scale() at new rows is not their own
scaleFit <- bartFit("scale(age)")
expect_equal(
  as.vector(replayed(scaleFit, newRows)),
  (newRows$age - mean(frame$age)) / stats::sd(frame$age)
)
polyFit <- bartFit("poly(dose, 2)")
expect_equal(
  as.vector(replayed(polyFit, newRows)),
  as.vector(stats::predict(stats::poly(frame$dose, 2), newRows$dose))
)
# arithmetic inside I() is evaluated on the rows given, as lm() evaluates it
centred <- bartFit("I(age - mean(age))")
expect_equal(
  as.vector(replayed(centred, newRows)),
  newRows$age - mean(newRows$age)
)
# under 'subset' the centre is still every row's
subsetFit <- bartFit("scale(age)", subset = keep)
expect_equal(
  as.vector(replayed(subsetFit, newRows)),
  (newRows$age - mean(frame$age)) / stats::sd(frame$age)
)
# predictions at a row do not depend on the rows beside it
together <- predict(scaleFit, newRows)
for (i in seq_len(nrow(newRows))) {
  expect_equal(
    unname(predict(scaleFit, newRows[i, ])),
    unname(together[, i, drop = FALSE])
  )
}
# a level the new rows lack is not an error, and the width is the fit's; a
# level the fit never saw is refused, as is one the fit's rows did not keep
levelFit <- bartFit("g")
expect_identical(ncol(replayed(levelFit, newRows)), 3L)
unseen <- newRows
unseen$g[2L] <- "new"
expect_error(predict(levelFit, unseen), "new")
droppedFit <- bartFit("g", subset = noW)
expect_identical(colnames(droppedFit$bases[[2L]]), c("gu", "gv"))
expect_identical(
  replayed(droppedFit, newRows),
  atNewRows("g", newRows)[, 1:2, drop = FALSE]
)
withW <- newRows
withW$g[1L] <- "w"
expect_error(
  predict(droppedFit, withW),
  "'basis' (g) has the level 'w' at a new row, which no row of the fit had",
  fixed = TRUE
)
# the order of a factor's columns at new rows is the fit's, whatever order
# the new rows' own levels come in
ordered <- frame
ordered$h <- factor(frame$g, levels = c("w", "u", "v"))
orderFit <- bartFit("h", ordered)
expect_identical(colnames(orderFit$bases[[2L]]), c("hw", "hu", "hv"))
reordered <- newRows
reordered$g <- c("v", "u", "w", "w", "v", "u", "v", "u", "w")
inFitOrder <- cbind(
  hw = as.double(reordered$g == "w"),
  hu = as.double(reordered$g == "u"),
  hv = as.double(reordered$g == "v")
)
for (levelOrder in list(c("v", "u", "w"), c("u", "v", "w"), c("w", "u", "v"))) {
  reordered$h <- factor(reordered$g, levels = levelOrder)
  expect_identical(replayed(orderFit, reordered), inFitOrder)
}
reordered$h <- reordered$g
expect_identical(replayed(orderFit, reordered), inFitOrder)
# a column of the data must be among the new rows
withoutAge <- newRows
withoutAge$age <- NULL
age <- frame$age
expect_error(
  predict(scaleFit, withoutAge),
  paste0(
    "'newdata' is missing variable 'age', required by forest 2's basis ",
    "(scale(age)); supply it, or give that basis at the new rows with ",
    "'bases ='"
  ),
  fixed = TRUE
)
rm(age)
# a number found where the formula was written is used, and looked up again;
# a vector with a value for every fitted row is refused at new rows
shrink <- 30
perRow <- frame$dose
constantFit <- bartFit("I(dose / shrink)")
expect_identical(
  as.vector(replayed(constantFit, newRows)),
  newRows$dose / 30
)
shrink <- 3
expect_identical(
  as.vector(replayed(constantFit, newRows)),
  newRows$dose / 3
)
perRowFit <- bartFit("I(perRow * 2)")
expect_error(
  predict(perRowFit, newRows),
  paste0(
    "'newdata' is missing variable 'perRow', required by forest 2's basis ",
    "(I(perRow * 2)); supply it, or give that basis at the new rows with ",
    "'bases ='"
  ),
  fixed = TRUE
)
expect_identical(
  dim(predict(perRowFit, newRows, bases = newRows$dose * 2)),
  c(3L, 9L)
)
# a basis handed over as a value has no code to build from
valued <- packaged(valueFit(frame$dose))
expect_null(valued$basis.terms[[2L]])
expect_error(predict(valued, newRows), "give them through 'bases ='")
# two basis forests: each keeps the record of its own
twoBases <- bart(
  y ~ x1 +
    x2 +
    forest(x1, basis = scale(age), n.trees = 5L) +
    forest(x2, basis = scale(dose), n.trees = 5L),
  frame,
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 5L,
  n.burn = 0L,
  n.samples = 3L,
  keepTrees = TRUE,
  verbose = FALSE
)
expect_identical(
  vapply(twoBases$basis.terms[2:3], function(term) term$label, ""),
  c("scale(age)", "scale(dose)")
)
expect_equal(
  as.vector(dbarts:::replayForestBasis(
    twoBases$basis.terms[[3L]],
    newRows,
    3L
  )),
  (newRows$dose - mean(frame$dose)) / stats::sd(frame$dose)
)
expect_equal(
  unname(predict(twoBases, frame[1:5, ])),
  unname(predict(twoBases, frame)[, 1:5])
)
# the record keeps the place the formula was written: a function that exists
# only there is found at predict
fitLocal <- function() {
  half <- function(x) x / 2
  bart(
    y ~ x1 + x2 + forest(x1, basis = scale(half(age)), n.trees = 5L),
    frame,
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 5L,
    n.burn = 0L,
    n.samples = 3L,
    keepTrees = TRUE,
    verbose = FALSE
  )
}
expect_equal(
  as.vector(replayed(fitLocal(), newRows)),
  (newRows$age - mean(frame$age)) / stats::sd(frame$age)
)
# a fit saved with its sampler state and read back predicts the same
scaleFit$fit$storeState()
path <- tempfile(fileext = ".rds")
saveRDS(scaleFit, path)
expect_equal(predict(readRDS(path), newRows), together)
unlink(path)
# the partial dependence grids rebuild the basis at their rows, a basis
# variable that is also a predictor among them
x1Fit <- bartFit("scale(x1)")
pd <- pdbart(
  x1Fit,
  xind = "x1",
  levs = list(c(0.2, 0.5)),
  newdata = frame,
  pl = FALSE
)
byHand <- vapply(
  c(0.2, 0.5),
  function(level) {
    rows <- frame
    rows$x1 <- level
    rowMeans(predict(x1Fit, rows))
  },
  numeric(3L)
)
expect_equal(matrix(as.vector(pd$fd[[1L]]), 3L), matrix(as.vector(byHand), 3L))
pd2 <- pd2bart(
  x1Fit,
  xind = c("x1", "x2"),
  levs = list(c(0.2, 0.5), c(0.3, 0.6)),
  newdata = frame,
  pl = FALSE
)
byHand <- vapply(
  list(c(0.2, 0.3), c(0.5, 0.3), c(0.2, 0.6), c(0.5, 0.6)),
  function(level) {
    rows <- frame
    rows$x1 <- level[1L]
    rows$x2 <- level[2L]
    mean(predict(x1Fit, rows))
  },
  numeric(1L)
)
expect_equal(as.vector(colMeans(pd2$fd)), byHand)
# a term that draws random numbers is evaluated once: after the fit R's
# generator stands where one evaluation of the basis leaves it
set.seed(11)
invisible(bartFit("scale(age) + runif(length(age))"))
after <- stats::runif(1L)
set.seed(11)
invisible(stats::runif(n))
expect_equal(after, stats::runif(1L))
# and what a term warns of is raised once
warny <- function(x) {
  warning("warny")
  x
}
countWarnings <- 0L
withCallingHandlers(
  bartFit("warny(dose)"),
  warning = function(w) {
    if (identical(conditionMessage(w), "warny")) {
      countWarnings <<- countWarnings + 1L
    }
    invokeRestart("muffleWarning")
  }
)
expect_identical(countWarnings, 1L)

## --- Block F: 'subset' is read once ------------------------------------------
# The data object evaluates 'subset' once and every basis is cut by the rows
# it kept, so a subset that draws its rows, or does something each time it is
# read, gives a basis the fit's own rows. An `id` predictor names each row
# the fit kept, in the fit's order.
rowed <- frame[c("y", "x1", "x2", "dose", "age", "g")]
rowed$id <- seq_len(n)
row.names(rowed) <- paste0("r", seq_len(n))
xRowed <- as.matrix(rowed[c("x1", "id")])
yRowed <- rowed$y
doseRows <- rowed$dose
ageRows <- rowed$age
gRows <- rowed$g
gIndicators <- outer(gRows, c("u", "v", "w"), "==") + 0
readings <- 0L
drawRows <- function() {
  readings <<- readings + 1L
  sample(n, 50L)
}
keptLogical <- rowed$id %% 3L != 0L
outOfOrder <- c(seq(80L, 2L, by = -3L), 5L, 5L, 41L)
byName <- row.names(rowed)[seq(79L, 1L, by = -2L)]
bartSettings <- list(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 5L,
  n.samples = 3L,
  n.burn = 0L,
  keepTrees = TRUE,
  verbose = FALSE,
  seed = 23L
)
# each door as a call, given the data and its subset as written: a term of
# the formula with and without a tilde, through bart(), a 'forests' list of
# code, of tildes and of values, the matrix interface with the three, and a
# data object carrying its bases. A second basis forest rides along.
subsetDoors <- list(
  term = function(data, subset) {
    bquote(dbarts(
      y ~ x1 + id + forest(x1, basis = dose) + forest(x1, basis = g),
      .(data),
      subset = .(subset),
      control = basisControl()
    ))
  },
  termTilde = function(data, subset) {
    bquote(dbarts(
      y ~ x1 + id + forest(x1, basis = ~dose) + forest(x1, basis = ~g),
      .(data),
      subset = .(subset),
      control = basisControl()
    ))
  },
  bart = function(data, subset) {
    as.call(c(
      as.list(bquote(bart(
        y ~ x1 + id + forest(x1, basis = dose) + forest(x1, basis = g),
        .(data),
        subset = .(subset)
      ))),
      bartSettings
    ))
  },
  list = function(data, subset) {
    bquote(dbarts(
      y ~ x1 + id,
      .(data),
      subset = .(subset),
      forests = list(
        forest(),
        forest(x1, basis = dose),
        forest(x1, basis = g)
      ),
      control = basisControl()
    ))
  },
  listTilde = function(data, subset) {
    bquote(dbarts(
      y ~ x1 + id,
      .(data),
      subset = .(subset),
      forests = list(
        forest(),
        forest(x1, basis = ~dose),
        forest(x1, basis = ~g)
      ),
      control = basisControl()
    ))
  },
  listValue = function(data, subset, rows = quote(seq_len(n))) {
    bquote(dbarts(
      y ~ x1 + id,
      .(data),
      subset = .(subset),
      forests = list(
        forest(),
        do.call(forest, list("x1", basis = doseRows[.(rows)])),
        do.call(forest, list("x1", basis = factor(gRows)[.(rows)]))
      ),
      control = basisControl()
    ))
  },
  matrix = function(data, subset, rows = quote(seq_len(n))) {
    bquote(dbarts(
      xRowed[.(rows), ],
      yRowed[.(rows)],
      subset = .(subset),
      forests = list(
        forest(),
        forest(x1, basis = doseRows[.(rows)]),
        forest(x1, basis = gRows[.(rows)])
      ),
      control = basisControl()
    ))
  },
  matrixTilde = function(data, subset, rows = quote(seq_len(n))) {
    bquote(dbarts(
      xRowed[.(rows), ],
      yRowed[.(rows)],
      subset = .(subset),
      forests = list(
        forest(),
        forest(x1, basis = ~ doseRows[.(rows)]),
        forest(x1, basis = ~ gRows[.(rows)])
      ),
      control = basisControl()
    ))
  },
  matrixValue = function(data, subset, rows = quote(seq_len(n))) {
    bquote(dbarts(
      xRowed[.(rows), ],
      yRowed[.(rows)],
      subset = .(subset),
      forests = list(
        forest(),
        do.call(forest, list("x1", basis = doseRows[.(rows)])),
        do.call(forest, list("x1", basis = factor(gRows)[.(rows)]))
      ),
      control = basisControl()
    ))
  },
  dataObject = function(data, subset, rows = quote(seq_len(n))) {
    bquote(dbarts(
      dbartsData(
        y ~ x1 + id,
        .(data),
        subset = .(subset),
        bases = list(NULL, doseRows[.(rows)], gIndicators[.(rows), ])
      ),
      control = basisControl()
    ))
  }
)
matrixDoors <- c("matrix", "matrixTilde", "matrixValue")
samplerOf <- function(fit) if (inherits(fit, "bart")) fit$fit else fit
drawsOf <- function(fit) {
  if (inherits(fit, "bart")) {
    list(fit$yhat.train, fit$sigma)
  } else {
    draws(fit)
  }
}
subsets <- list(
  drawn = quote(drawRows()),
  outOfOrder = quote(outOfOrder),
  logical = quote(keptLogical),
  names = quote(byName)
)
for (door in names(subsetDoors)) {
  for (kind in names(subsets)) {
    if (kind == "names" && door %in% matrixDoors) {
      # the matrix interface takes no row names for 'subset'
      next
    }
    info <- paste(door, kind)
    readings <- 0L
    set.seed(7L)
    fit <- eval(subsetDoors[[door]](quote(rowed), subsets[[kind]]))
    if (kind == "drawn") {
      expect_identical(readings, 1L, info = info)
    }
    sampler <- samplerOf(fit)
    rows <- as.integer(sampler$data@x[, "id"])
    expect_identical(
      rows,
      switch(
        kind,
        drawn = {
          set.seed(7L)
          sample(n, 50L)
        },
        outOfOrder = outOfOrder,
        logical = which(keptLogical),
        names = seq(79L, 1L, by = -2L)
      ),
      info = info
    )
    # the basis of exactly the fit's rows, in the fit's order
    expect_identical(
      as.vector(sampler$data@bases[[2L]]),
      doseRows[rows],
      info = info
    )
    expect_identical(
      unname(sampler$data@bases[[3L]]),
      unname(gIndicators[rows, ]),
      info = info
    )
    # and the model of the data already cut by hand
    byHand <- if (door %in% c(matrixDoors, "listValue", "dataObject")) {
      eval(subsetDoors[[door]](quote(rowed[rows, ]), NULL, quote(rows)))
    } else {
      eval(subsetDoors[[door]](quote(rowed[rows, ]), NULL))
    }
    expect_identical(
      unname(samplerOf(byHand)$data@bases[[2L]]),
      unname(sampler$data@bases[[2L]]),
      info = info
    )
    expect_identical(drawsOf(fit), drawsOf(byHand), info = info)
  }
}
# a data object built in the call is built once
readings <- 0L
set.seed(7L)
builtOnce <- dbarts(
  dbartsData(
    y ~ x1 + id,
    rowed,
    subset = drawRows(),
    bases = list(NULL, doseRows)
  ),
  control = basisControl()
)
expect_identical(readings, 1L)
# and a plain fit reads its subset once
readings <- 0L
invisible(dbarts(
  y ~ x1 + id,
  rowed,
  subset = drawRows(),
  control = basisControl()
))
expect_identical(readings, 1L)
# what each term computes across rows comes from every row at the matrix
# interface too, whatever order 'subset' keeps them in
matrixComputed <- dbarts(
  xRowed,
  yRowed,
  subset = outOfOrder,
  forests = list(
    forest(),
    forest(x1, basis = scale(ageRows)),
    forest(x1, basis = poly(doseRows, 2))
  ),
  control = basisControl()
)
expect_identical(
  matrixComputed$data@bases[[2L]],
  lmColumns("scale(ageRows)", rows = outOfOrder)
)
expect_identical(
  matrixComputed$data@bases[[3L]],
  lmColumns("poly(doseRows, 2)", rows = outOfOrder)
)
expect_equal(
  as.vector(matrixComputed$data@bases[[2L]]),
  ((ageRows - mean(ageRows)) / stats::sd(ageRows))[outOfOrder]
)
# every forest's basis is built on the kept rows, a third forest's too: a
# level 'subset' empties is no column of it
withoutW <- rowed$g != "w"
thirdForest <- list(
  term = dbarts(
    y ~ x1 + id + forest(x1, basis = dose) + forest(x1, basis = g),
    rowed,
    subset = withoutW,
    control = basisControl()
  ),
  list = dbarts(
    y ~ x1 + id,
    rowed,
    subset = withoutW,
    forests = list(forest(), forest(x1, basis = dose), forest(x1, basis = g)),
    control = basisControl()
  ),
  matrix = dbarts(
    xRowed,
    yRowed,
    subset = withoutW,
    forests = list(
      forest(),
      forest(x1, basis = doseRows),
      forest(x1, basis = gRows)
    ),
    control = basisControl()
  )
)
for (door in names(thirdForest)) {
  expect_identical(
    unname(thirdForest[[door]]$data@bases[[3L]]),
    unname(gIndicators[withoutW, 1:2]),
    info = door
  )
}
# the rows of the data are counted from the data frame, whatever the formula
# names first: a number beside the response is no row count
kTimes <- 2
scaledResponse <- dbarts(
  I(kTimes * y) ~ x1 + id + forest(x1, basis = dose),
  rowed,
  subset = keptLogical,
  control = basisControl()
)
expect_identical(
  as.integer(scaledResponse$data@x[, "id"]),
  which(keptLogical)
)
expect_identical(scaledResponse$data@y, 2 * yRowed[keptLogical])
expect_identical(
  as.vector(scaledResponse$data@bases[[2L]]),
  doseRows[keptLogical]
)

## --- Block G: 'data' is read once --------------------------------------------
# The 'data' argument is evaluated once for a fit, and 'subset' in that one
# value of it. A 'data' that draws its rows therefore gives the response, the
# predictors and every basis read against its columns one draw of them, at
# each door, with a 'subset' over its columns and without. A basis of one
# numeric column and one of two levels ride along, and `id` names the row of
# `twoLevel` each row of the fit is.
twoLevel <- frame[c("y", "x1", "x2", "dose", "zl")]
twoLevel$id <- as.double(seq_len(n))
zlIndicators <- cbind(!twoLevel$zl, twoLevel$zl) + 0
dataDraws <- list()
drawData <- function(size = n, replace = FALSE) {
  rows <- sample(n, size, replace)
  dataDraws[[length(dataDraws) + 1L]] <<- rows
  twoLevel[rows, ]
}
samplerSettings <- list(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 5L,
  n.samples = 3L,
  n.burn = 0L,
  verbose = FALSE,
  samplerOnly = TRUE
)
dataDoors <- list(
  term = function(data, subset) {
    bquote(dbarts(
      y ~ x1 + id + forest(x1, basis = dose) + forest(x1, basis = zl),
      .(data),
      subset = .(subset),
      control = basisControl()
    ))
  },
  termTilde = function(data, subset) {
    bquote(dbarts(
      y ~ x1 + id + forest(x1, basis = ~dose) + forest(x1, basis = ~zl),
      .(data),
      subset = .(subset),
      control = basisControl()
    ))
  },
  bart = function(data, subset) {
    as.call(c(
      as.list(bquote(bart(
        y ~ x1 + id + forest(x1, basis = dose) + forest(x1, basis = zl),
        .(data),
        subset = .(subset)
      ))),
      samplerSettings
    ))
  },
  list = function(data, subset) {
    bquote(dbarts(
      y ~ x1 + id,
      .(data),
      subset = .(subset),
      forests = list(
        forest(),
        forest(x1, basis = dose),
        forest(x1, basis = zl)
      ),
      control = basisControl()
    ))
  },
  listTilde = function(data, subset) {
    bquote(dbarts(
      y ~ x1 + id,
      .(data),
      subset = .(subset),
      forests = list(
        forest(),
        forest(x1, basis = ~dose),
        forest(x1, basis = ~zl)
      ),
      control = basisControl()
    ))
  }
)
dataKinds <- list(
  shuffled = quote(drawData()),
  subsample = quote(drawData(60L)),
  bootstrap = quote(drawData(replace = TRUE))
)
for (door in names(dataDoors)) {
  for (kind in names(dataKinds)) {
    for (cut in c("every row", "x1 > 0.3")) {
      info <- paste(door, kind, cut, sep = ", ")
      dataDraws <- list()
      sampler <- eval(dataDoors[[door]](
        dataKinds[[kind]],
        if (cut == "every row") NULL else quote(x1 > 0.3)
      ))
      expect_identical(length(dataDraws), 1L, info = info)
      drawn <- dataDraws[[1L]]
      rows <- as.integer(sampler$data@x[, "id"])
      expect_identical(
        rows,
        if (cut == "every row") drawn else drawn[twoLevel$x1[drawn] > 0.3],
        info = info
      )
      expect_identical(sampler$data@y, twoLevel$y[rows], info = info)
      expect_identical(
        as.vector(sampler$data@bases[[2L]]),
        twoLevel$dose[rows],
        info = info
      )
      expect_identical(
        unname(sampler$data@bases[[3L]]),
        zlIndicators[rows, ],
        info = info
      )
    }
  }
}

## --- Block H: a value that the kept rows leave with one level ----------------
# A basis handed over as a value is cut by the data object after it is
# expanded, so 'subset' or a dropped row can leave a factor, a character or a
# logical vector one level, or a numeric column nothing but zeros. It is
# refused then as the same basis written as code is, in the same words, and
# not fitted with a column of ones beside a column of zeros.
zlRows <- twoLevel$zl
armRows <- ifelse(zlRows, "treated", "control")
xTwo <- as.matrix(twoLevel[c("x1", "x2")])
yTwo <- twoLevel$y
refusal <- function(fit) {
  tryCatch(
    {
      fit
      NA_character_
    },
    error = conditionMessage
  )
}
valueForest <- function(value) do.call(forest, list("x1", basis = value))
emptied <- list(
  factor = list(value = factor(armRows), code = quote(factor(armRows))),
  character = list(value = armRows, code = quote(armRows)),
  logical = list(value = zlRows, code = quote(zlRows)),
  numeric = list(value = as.double(!zlRows), code = quote(as.double(!zlRows)))
)
for (kind in names(emptied)) {
  value <- emptied[[kind]]$value
  written <- refusal(eval(bquote(dbarts(
    xTwo,
    yTwo,
    subset = which(zlRows),
    forests = list(forest(), forest(x1, basis = .(emptied[[kind]]$code))),
    control = basisControl()
  ))))
  expect_false(is.na(written), info = kind)
  # a value has no text to be named by
  words <- sub(" (zlRows)", "", written, fixed = TRUE)
  expect_identical(
    refusal(dbarts(
      y ~ x1 + x2,
      twoLevel,
      subset = zl,
      forests = list(forest(), valueForest(value)),
      control = basisControl()
    )),
    words,
    info = kind
  )
  expect_identical(
    refusal(dbarts(
      xTwo,
      yTwo,
      subset = which(zlRows),
      forests = list(forest(), valueForest(value)),
      control = basisControl()
    )),
    words,
    info = kind
  )
  # emptied by the rows a missing response drops
  missingResponse <- twoLevel
  missingResponse$y[!zlRows] <- NA_real_
  expect_identical(
    refusal(dbarts(
      y ~ x1 + x2,
      missingResponse,
      forests = list(forest(), valueForest(value)),
      control = basisControl()
    )),
    words,
    info = kind
  )
}
expect_identical(
  refusal(dbarts(
    xTwo,
    yTwo,
    subset = which(zlRows),
    forests = list(forest(), forest(x1, basis = zlRows)),
    control = basisControl()
  )),
  paste0(
    "'basis' (zlRows) is TRUE on every row the fit keeps, so its other level ",
    "has no observations and the forest would be multiplied by a constant"
  )
)
# a value that keeps both levels is cut and fitted as before
keptBoth <- dbarts(
  y ~ x1 + x2,
  twoLevel,
  subset = x1 > 0.3,
  forests = list(forest(), valueForest(zlRows), valueForest(twoLevel$dose)),
  control = basisControl()
)
expect_identical(
  unname(keptBoth$data@bases[[2L]]),
  zlIndicators[twoLevel$x1 > 0.3, ]
)
expect_identical(
  as.vector(keptBoth$data@bases[[3L]]),
  twoLevel$dose[twoLevel$x1 > 0.3]
)
