# Small surface fixes of several-forest models: the refusals name what the
# model is and the way to state a size, a forest's sd is read in the prior
# vocabulary, dbartsSpec() takes the retired proposal.probs as the other doors
# do, and a forest() specification prints.

forest <- dbartsForests$forest
varianceForest <- dbartsForests$varianceForest

set.seed(3)
n <- 60L
frame <- data.frame(y = rnorm(n), x1 = rnorm(n), x2 = rnorm(n), w = rnorm(n))
refusal <- function(expr) {
  tryCatch(
    {
      force(expr)
      "accepted"
    },
    error = conditionMessage
  )
}

## --- messages ---------------------------------------------------------------
expect_identical(
  refusal(dbarts(
    y ~ forest(x1) + forest(x2, basis = w),
    frame,
    leaf.prior = normal(k = 3)
  )),
  paste0(
    "a model of several forests does not support a non-default 'k'; drop it ",
    "or fit a single-forest model; a forest's size is stated as ",
    "forest(sd = ), in the response's units"
  )
)
counts <- transform(frame, y = rpois(n, 4))
expect_identical(
  refusal(dbarts(
    y ~ x1 + x2,
    counts,
    forests = list(forest(), forest(basis = w)),
    family = "nbinom"
  )),
  paste0(
    "a model of several forests does not support family \"nbinom\": its ",
    "shape block is not shown to interleave with the amplitude block"
  )
)
invchiRefusal <- paste0(
  "forest 'sd' must be a number, not invchi(): a hyperprior on a forest's ",
  "sd is not supported yet"
)
expect_identical(
  refusal(dbarts(
    y ~ forest(x1) + forest(x2, basis = w, sd = invchi(3, 1)),
    frame
  )),
  invchiRefusal
)
expect_identical(
  refusal(dbarts(
    y ~ x1 + x2,
    frame,
    forests = list(forest(), forest(basis = w, sd = invchi(3, 1)))
  )),
  invchiRefusal
)

## --- dbartsSpec(proposal.probs = ) -------------------------------------------
data <- dbartsData(y ~ x1 + x2, frame)
mixture <- c(birth_death = 0.5, change = 0.5)
warned <- NULL
spec <- withCallingHandlers(
  dbartsSpec(data, proposal.probs = mixture),
  warning = function(w) {
    warned <<- c(warned, list(w))
    invokeRestart("muffleWarning")
  }
)
expect_equal(length(warned), 1L)
expect_true(inherits(warned[[1L]], "dbartsDeprecatedWarning"))
expect_identical(
  conditionMessage(warned[[1L]]),
  paste0(
    "'proposal.probs' has left 'dbartsSpec': the tree-move mixture is a ",
    "sampler setting: write control = dbartsControl(proposal.probs = ); ",
    "'proposal.probs' is removed in dbarts 1.1-0. The value was used."
  )
)
viaControl <- dbartsSpec(
  data,
  control = dbartsControl(proposal.probs = mixture)
)
expect_identical(spec$control@proposal.probs, viaControl$control@proposal.probs)
expect_identical(
  refusal(suppressWarnings(dbartsSpec(
    data,
    control = dbartsControl(proposal.probs = mixture),
    proposal.probs = mixture
  ))),
  paste0(
    "'control' cannot be combined with 'proposal.probs': set the tree-move ",
    "mixture in one place, dbartsControl(proposal.probs = )"
  )
)

## --- printing a forest ------------------------------------------------------
expect_identical(
  capture.output(print(forest(x1 + x2, basis = w, sd = 2))),
  c(
    "dbarts forest specification",
    "  vars  = x1 + x2",
    "  basis = w",
    "  sd    = 2"
  )
)
expect_identical(
  capture.output(print(forest())),
  c("dbarts forest specification", "  every knob at its default")
)
expect_identical(
  capture.output(print(forest(
    amplitude = dbartsPriors$fixed(),
    interactions = dbartsForests$interactions(max.order = 2)
  ))),
  c(
    "dbarts forest specification",
    "  amplitude    = fixed()",
    "  interactions = interactions()"
  )
)
sampler <- dbarts(y ~ forest(x1) + forest(x2, basis = w, sd = 0.5), frame)
read <- sampler$getLeafPrior()
expect_identical(
  capture.output(print(read[[2L]]$leaf.prior)),
  c("dbarts forest specification", "  sd = 0.5")
)
expect_false(any(grepl("NULL", capture.output(print(read[[1L]]$leaf.prior)))))
