# predict, extract, fitted and for bart and rbart objects

extract <- function(object, ...) UseMethod("extract")

plotTree <- function(object, ...) UseMethod("plotTree")

# latent-scale draws to probabilities for a binary fit; fits saved before
# the family element existed are all probit
probabilityFromLatents <- function(latents, object) {
  if (identical(object[["family"]], "logistic")) {
    plogis(latents)
  } else {
    pnorm(latents)
  }
}

# per-draw, per-observation log-likelihood of the stored training response.
# ev enters with chains split ((n.chains x) n.samples x n.obs), so that
# as.vector(ev) enumerates draws chain-fastest - the order as.vector on the
# stored sigma yields in both of its layouts (n.chains x n.samples matrix or
# chain-interleaved combined vector) - and the sigma draws pair by plain
# recycling. Weights enter as precision for gaussian fits (y | x ~ N(f(x),
# sigma^2 / w)) and as trial counts for weighted logistic ones; probit fits
# never store weights.
pointwiseLogLikelihood <- function(object, ev) {
  y <- object[["y"]]
  if (is.null(y)) {
    stop(
      "cannot compute the log-likelihood; fit does not store the training response"
    )
  }
  weights <- object[["weights"]]
  n.draws <- length(ev) %/% length(y)
  y <- rep(y, each = n.draws)
  if (is.null(object[["sigma"]])) {
    result <- dbinom(y, 1L, as.vector(ev), log = TRUE)
    if (!is.null(weights)) {
      result <- rep(weights, each = n.draws) * result
    }
  } else {
    sd <- rep_len(as.vector(object$sigma), length(ev))
    if (!is.null(weights)) {
      sd <- sd / rep(sqrt(weights), each = n.draws)
    }
    result <- dnorm(y, as.vector(ev), sd, log = TRUE)
  }
  array(result, dim(ev))
}

# per-observation posterior summary for the interval-returning generics: est
# (the posterior mean) plus a symmetric ci.level credible band from the draw
# quantiles, pooled over all samples and chains (observations are the last
# array margin, as in the mean path). The interval KIND follows the caller's
# type: "ev" gives a credible interval for E[Y|x] (a probability for binary),
# "ppd" a prediction interval that also carries the residual noise, and "bart"
# a credible interval on the latent scale.
posteriorInterval <- function(draws, ci.level) {
  if (
    !is.numeric(ci.level) ||
      length(ci.level) != 1L ||
      is.na(ci.level) ||
      ci.level <= 0 ||
      ci.level >= 1
  ) {
    stop("'ci.level' must be a single number in (0, 1)")
  }
  probs <- c((1 - ci.level) / 2, 1 - (1 - ci.level) / 2)
  if (is.null(dim(draws))) {
    result <- matrix(
      c(mean(draws), quantile(draws, probs, names = FALSE)),
      nrow = 1L
    )
  } else {
    obsMargin <- length(dim(draws))
    est <- apply(draws, obsMargin, mean)
    bounds <- apply(draws, obsMargin, quantile, probs = probs, names = FALSE)
    result <- cbind(est, bounds[1L, ], bounds[2L, ])
  }
  colnames(result) <- c("est", "ci.lower", "ci.upper")
  result
}

combineOrUncombineChains <- function(x, n.chains, combineChains) {
  if (n.chains > 1L) {
    if (length(dim(x)) > 2L && combineChains) {
      x <- combineChains(x)
    } else if (length(dim(x)) == 2L && !combineChains) {
      x <- uncombineChains(x, n.chains)
    }
  }
  x
}

predict.bart <- function(
  object,
  newdata,
  offset,
  weights,
  type = c("ev", "ppd", "bart"),
  combineChains = TRUE,
  n.threads = object$fit$control@n.threads,
  ci.level = NULL,
  ...
) {
  if (missing(offset)) {
    offset <- NULL
  }
  if (missing(weights)) {
    weights <- NULL
  }

  if (is.null(object[["fit"]])) {
    if (as.character(object$call[[1L]]) == "bart2") {
      stop("predict requires bart2 to be called with 'keepTrees' == TRUE")
    } else {
      stop("predict requires bart to be called with 'keeptrees' == TRUE")
    }
  }

  if (is.character(type)) {
    if (type[1L] == "response") {
      type[1L] <- "ev"
    } else if (type[1L] == "link") {
      type[1L] <- "bart"
    }
  }
  if (
    !is.character(type) ||
      length(type) == 0L ||
      type[1L] %not_in% eval(formals(predict.bart)$type)
  ) {
    stop(
      "type must be in '",
      paste0(eval(formals(predict.bart)$type), collapse = "', '"),
      "'"
    )
  }
  type <- type[1L]

  n.threads <- as.integer(n.threads)[1L]

  result <- object$fit$predict(newdata, offset, n.threads)
  # result is n.obs x n.samples x n.chains
  n.chains <- object$fit$control@n.chains
  result <- convertSamplesFromDbartsToBart(result, n.chains, combineChains)

  if (type != "bart") {
    # nolint next: object_usage_linter. named for readability; value drives if.
    if ((responseIsBinary <- is.null(object[["sigma"]]))) {
      result <- probabilityFromLatents(result, object)
    }

    if (type == "ppd") {
      result <- sampleFromPPD(result, object, weights, n.chains)
    }
  }

  # ci.level opts into a per-observation est + credible band (kind follows type)
  if (!is.null(ci.level)) {
    return(posteriorInterval(result, ci.level))
  }

  result
}

extract.bart <- function(
  object,
  type = c("ev", "ppd", "bart", "loglik", "trees"),
  sample = c("train", "test"),
  combineChains = TRUE,
  ...
) {
  if (is.character(type)) {
    if (type[1L] == "response") {
      type[1L] <- "ev"
    } else if (type[1L] == "link") {
      type[1L] <- "bart"
    }
  }
  if (
    !is.character(type) || type[1L] %not_in% eval(formals(extract.bart)$type)
  ) {
    stop(
      "type must be in '",
      paste0(eval(formals(extract.bart)$type), collapse = "', '"),
      "'"
    )
  }
  type <- type[1L]

  if (type == "trees") {
    if (is.null(object$fit)) {
      if (as.character(object$call[[1L]]) == "bart2") {
        stop(
          "extracting trees requires bart2 to be called with 'keepTrees' == TRUE"
        )
      } else {
        stop(
          "extracting trees requires bart to be called with 'keeptrees' == TRUE"
        )
      }
    }
    treesCall <- match.call()
    target <- quote(object$fit$getTrees)
    target[[2L]][[2L]] <- treesCall$object
    treesCall[[1L]] <- target
    treesCall$object <- NULL
    treesCall$type <- NULL
    return(eval(treesCall, parent.frame()))
  }

  if (
    !is.character(sample) ||
      sample[1L] %not_in% eval(formals(extract.bart)$sample)
  ) {
    stop(
      "sample must be in '",
      paste0(eval(formals(extract.bart)$sample), collapse = "', '"),
      "'"
    )
  }
  sample <- sample[1L]

  # the log-likelihood is against the stored training response; there is no
  # test response to evaluate
  if (type == "loglik" && sample == "test") {
    stop("cannot extract a test sample log-likelihood; no test response exists")
  }

  if (sample == "test" && is.null(object[["yhat.test"]])) {
    stop(
      "cannot extract test sample predictions if no test data exists; use `predict` instead"
    )
  }
  if (sample == "train" && is.null(object[["yhat.train"]])) {
    if (as.character(object$call[[1L]]) == "bart2") {
      stop(
        "cannot extract train sample predictions; bart2 must be called with 'keepTrainingFits' == TRUE"
      )
    } else {
      stop(
        "cannot extract train sample predictions; bart must be called with 'keeptrainfits' == TRUE"
      )
    }
  }

  n.chains <- if (!is.null(object[["fit"]])) {
    object$fit$control@n.chains
  } else {
    object$n.chains
  }

  if (type == "loglik") {
    ev <- extract.bart(
      object,
      type = "ev",
      sample = "train",
      combineChains = FALSE
    )
    result <- pointwiseLogLikelihood(object, ev)
    return(combineOrUncombineChains(result, n.chains, combineChains))
  }

  result <- if (sample == "train") object$yhat.train else object$yhat.test
  weights <- if (sample == "train") object$weights else object$weights.test
  #n.samples <- if (length(dim(result)) > 2L) dim(result)[2L] else dim(result)[1L] %/% n.chains
  #n.obs     <- if (length(dim(result)) > 2L) dim(result)[3L] else dim(result)[2L]

  result <- combineOrUncombineChains(result, n.chains, combineChains)

  if (type == "bart") {
    return(result)
  }

  # nolint next: object_usage_linter. named for readability; value drives if.
  if ((responseIsBinary <- is.null(object[["sigma"]]))) {
    result <- probabilityFromLatents(result, object)
  }

  if (type == "ppd") {
    result <- sampleFromPPD(result, object, weights, n.chains)
  }

  result
}

fitted.bart <- function(
  object,
  type = c("ev", "ppd", "bart"),
  sample = c("train", "test"),
  ci.level = NULL,
  ...
) {
  if (is.character(type)) {
    if (type[1L] == "response") {
      type[1L] <- "ev"
    } else if (type[1L] == "link") {
      type[1L] <- "bart"
    }
  }
  if (
    !is.character(type) || type[1L] %not_in% eval(formals(fitted.bart)$type)
  ) {
    stop(
      "type must be in '",
      paste0(eval(formals(fitted.bart)$type), collapse = "', '"),
      "'"
    )
  }
  type <- type[1L]

  if (
    !is.character(sample) ||
      sample[1L] %not_in% eval(formals(fitted.bart)$sample)
  ) {
    stop(
      "sample must be in '",
      paste0(eval(formals(fitted.bart)$sample), collapse = "', '"),
      "'"
    )
  }
  sample <- sample[1L]

  result <- extract(object, type, sample, ...)

  # ci.level opts into a per-observation est + credible band instead of the
  # posterior mean; the interval kind follows type (see posteriorInterval)
  if (!is.null(ci.level)) {
    return(posteriorInterval(result, ci.level))
  }

  if (!is.null(dim(result))) {
    apply(result, length(dim(result)), mean)
  } else {
    mean(result)
  }
}

residuals.bart <- function(object, type = "ev", ...) {
  # type flows to fitted so link-scale (type = "bart") residuals are reachable;
  # residuals are always against the training response, so sample is pinned
  object$y - fitted.bart(object, type = type, sample = "train", ...)
}

predict.rbart <- function(
  object,
  newdata,
  group.by,
  offset,
  weights,
  type = c("ev", "ppd", "bart", "ranef"),
  combineChains = TRUE,
  ci.level = NULL,
  ...
) {
  if (is.null(object$fit)) {
    stop("predict requires rbart to be called with 'keepTrees' == TRUE")
  }

  dotsList <- list(...)
  if (!is.null(dotsList[["value"]])) {
    warning("argument 'value' has been deprecated; use 'type' instead")
    type <- dotsList[["value"]]
    dotsList[["value"]] <- NULL
  }

  if (is.character(type)) {
    if (type[1L] == "response") {
      type[1L] <- "ev"
    } else if (type[1L] == "link") {
      type[1L] <- "bart"
    }
  }
  if (is.character(type) && length(type) > 0L && type[1L] == "post-mean") {
    warning("type of 'post-mean' for predict deprecated; use 'ev' instead")
    type[1L] <- "ev"
  }
  if (
    !is.character(type) ||
      length(type) == 0L ||
      type[1L] %not_in% eval(formals(predict.rbart)$type)
  ) {
    stop(
      "type must be in '",
      paste0(eval(formals(predict.rbart)$type), collapse = "', '"),
      "'"
    )
  }
  type <- type[1L]

  if (missing(offset)) {
    offset <- NULL
  }
  if (missing(weights)) {
    weights <- NULL
  }

  n.chains <- if (is.null(object$n.chains)) {
    length(object$fit)
  } else {
    object$n.chains
  }
  n.samples <- object$fit[[1L]]$control@n.samples

  nonParametricPart <- 0
  # collects results in an array of n.obs x n.samples x n.chains, default for
  # internal sampler
  #
  # utilize bart stuff to get n.obs, since we would otherwise have to build
  # the test matrix
  if (type != "ranef") {
    if (n.chains > 1L) {
      if (length(object$fit) == 1L) {
        # the in-core Gibbs path keeps one multi-chain sampler, whose
        # predictions already carry the chain dimension
        nonParametricPart <- object$fit[[1L]]$predict(newdata, offset)
        n.obs <- dim(nonParametricPart)[1L]
      } else {
        n.obs <- NULL
        nonParametricPart <- array(
          sapply(seq_len(n.chains), function(i) {
            res <- object$fit[[i]]$predict(newdata, offset)
            if (is.null(n.obs)) {
              n.obs <<- dim(res)[1L]
            }
            res
          }),
          c(n.obs, n.samples, n.chains)
        )
      }
    } else {
      nonParametricPart <- object$fit[[1L]]$predict(newdata, offset)
      n.obs <- nrow(nonParametricPart)
    }
    if (n.obs != length(group.by)) {
      stop("length of group.by not equal to number of rows in test")
    }

    nonParametricPart <- convertSamplesFromDbartsToBart(
      nonParametricPart,
      n.chains,
      combineChains
    )
  }

  if (type == "bart") {
    if (!is.null(ci.level)) {
      return(posteriorInterval(nonParametricPart, ci.level))
    }
    return(nonParametricPart)
  }

  ranef <- 0
  if (type != "bart") {
    ranefNames.test <- levels(group.by)
    ranefNames.train <- if (length(dim(object$ranef)) > 2L) {
      dimnames(object$ranef)[[3L]]
    } else {
      dimnames(object$ranef)[[2L]]
    }

    ranef <- object$ranef
    if (n.chains > 1L) {
      if (length(dim(ranef)) > 2L && combineChains) {
        ranef <- combineChains(ranef)
      } else if (length(dim(ranef)) == 2L && !combineChains && n.chains > 1L) {
        ranef <- uncombineChains(ranef, n.chains)
      }
    }

    if (!all(measuredLevels <- ranefNames.test %in% ranefNames.train)) {
      warning(
        "test includes random effect levels not present in training - ranef estimates default to draws from their latent distribution parameterized by the posterior of its variance; draws may not be the same across future calls to 'predict'"
      )
      n.unmeasured <- sum(!measuredLevels)
      if (n.chains > 1L) {
        if (!combineChains) {
          unmeasuredRanef <- array(
            rnorm(
              n.chains * n.samples * n.unmeasured,
              0,
              rep.int(object$tau, n.unmeasured)
            ),
            c(n.chains, n.samples, n.unmeasured),
            dimnames = list(NULL, NULL, ranefNames.test[!measuredLevels])
          )
        } else {
          unmeasuredRanef <- matrix(
            rnorm(
              n.chains * n.samples * n.unmeasured,
              0,
              rep.int(object$tau, n.unmeasured)
            ),
            n.chains * n.samples,
            n.unmeasured,
            dimnames = list(NULL, ranefNames.test[!measuredLevels])
          )
        }
        # branch on ranef's current shape (already reshaped above to match
        # this call's combineChains, which is what unmeasuredRanef was just
        # built against) rather than object$ranef's stored shape - those two
        # can disagree whenever the fit's own combineChains default differs
        # from the one requested here
        if (length(dim(ranef)) == 2L) {
          ranef <- cbind(ranef, unmeasuredRanef)
        } else {
          # ranef are n.chains x n.samples x n.group
          ranef <- array(
            c(ranef, unmeasuredRanef),
            c(n.chains, n.samples, dim(ranef)[3L] + n.unmeasured),

            dimnames = list(
              NULL,
              NULL,
              c(dimnames(ranef)[[3L]], dimnames(unmeasuredRanef)[[3L]])
            )
          )
        }
      } else {
        unmeasuredRanef <- matrix(
          rnorm(n.samples * n.unmeasured, 0, rep.int(object$tau, n.unmeasured)),
          n.samples,
          n.unmeasured,
          dimnames = list(NULL, ranefNames.test[!measuredLevels])
        )
        ranef <- cbind(ranef, unmeasuredRanef)
      }
    }
  }

  if (type == "ranef") {
    ranef <- if (length(dim(ranef)) > 2L) {
      ranef[,, ranefNames.test, drop = FALSE]
    } else {
      ranef[, ranefNames.test, drop = FALSE]
    }
    ranef <- combineOrUncombineChains(ranef, n.chains, combineChains)
    if (!is.null(ci.level)) {
      return(posteriorInterval(ranef, ci.level))
    }
    return(ranef)
  }

  ranef <- unname(
    if (length(dim(ranef)) > 2L) {
      ranef[,, as.character(group.by), drop = FALSE]
    } else {
      ranef[, as.character(group.by), drop = FALSE]
    }
  )
  ranef <- combineOrUncombineChains(ranef, n.chains, combineChains)

  if (
    length(dim(nonParametricPart)) != length(dim(ranef)) ||
      any(dim(nonParametricPart) != dim(ranef))
  ) {
    stop(
      "internal error: fixed and random effect predictions have ",
      "mismatched dimensions"
    )
  }
  result <- nonParametricPart + ranef

  responseIsBinary <- is.null(object[["sigma"]])
  if (responseIsBinary) {
    result <- probabilityFromLatents(result, object)
  }

  if (type == "ppd") {
    result <- sampleFromPPD(result, object, weights, n.chains)
  }

  if (!is.null(ci.level)) {
    return(posteriorInterval(result, ci.level))
  }

  if (exists("unmeasuredRanef", inherits = FALSE)) {
    attr(result, "ranef") <- unmeasuredRanef
  }

  result
}

extract.rbart <- function(
  object,
  type = c("ev", "ppd", "bart", "loglik", "ranef", "trees"),
  sample = c("train", "test"),
  combineChains = TRUE,
  ...
) {
  if (is.character(type)) {
    if (type[1L] == "response") {
      type[1L] <- "ev"
    } else if (type[1L] == "link") {
      type[1L] <- "bart"
    }
  }
  if (
    !is.character(type) || type[1L] %not_in% eval(formals(extract.rbart)$type)
  ) {
    stop(
      "type must be in '",
      paste0(eval(formals(extract.rbart)$type), collapse = "', '"),
      "'"
    )
  }
  type <- type[1L]

  n.chains <- if (is.null(object$n.chains)) {
    length(object$fit)
  } else {
    object$n.chains
  }

  if (type == "trees") {
    if (is.null(object$fit)) {
      stop(
        "extracting trees requires rbart to be called with 'keepTrees' == TRUE"
      )
    }
    treesCall <- match.call()
    # the in-core Gibbs path keeps one multi-chain sampler: route chain
    # selection through its chainNums argument instead of the fit list
    singleFit <- length(object$fit) == 1L && n.chains > 1L
    target <- quote(object$fit[[i]]$getTrees)
    target[[2L]][[2L]][[2L]] <- treesCall$object
    if (singleFit) {
      target[[2L]][[3L]] <- 1L
    }
    treesCall[[1L]] <- target
    treesCall$object <- NULL
    treesCall$type <- NULL
    treesCall$chainNums <- if (singleFit) quote(i) else NULL
    evalEnv <- parent.frame()
    dotsList <- list(...)
    chainNums <- if ("chainNums" %in% names(dotsList)) {
      as.integer(dotsList[["chainNums"]])
    } else {
      seq_len(n.chains)
    }
    varOrder <- c("chain", "sample", "tree", "n", "var", "value")
    allTrees <- lapply(chainNums, function(i) {
      result_i <- eval(subTermInLanguage(treesCall, quote(i), i), evalEnv)
      if (n.chains > 1L) {
        result_i$chain <- i
      }
      # varOrder's "chain" is absent for single-chain fits (getTrees omits
      # it); reorder only the columns that exist and keep any others
      # (directions/missing/beta.*) trailing in their original order
      knownOrder <- varOrder[varOrder %in% colnames(result_i)]
      result_i[, c(knownOrder, setdiff(colnames(result_i), knownOrder))]
    })
    if (length(allTrees) > 1L) {
      allTrees <- Reduce(rbind, allTrees)
    } else {
      allTrees <- allTrees[[1L]]
    }
    row.names(allTrees) <- as.character(seq_len(nrow(allTrees)))
    return(allTrees)
  }

  if (
    !is.character(sample) ||
      sample[1L] %not_in% eval(formals(extract.rbart)$sample)
  ) {
    stop(
      "sample must be in '",
      paste0(eval(formals(extract.rbart)$sample), collapse = "', '"),
      "'"
    )
  }
  sample <- sample[1L]

  # the log-likelihood is against the stored training response; there is no
  # test response to evaluate
  if (type == "loglik" && sample == "test") {
    stop("cannot extract a test sample log-likelihood; no test response exists")
  }

  if (sample == "test" && is.null(object[["yhat.test"]])) {
    stop(
      "cannot extract test sample predictions if no test data exists; use `predict` instead"
    )
  }

  if (type == "ranef") {
    ranefNames <- if (sample == "train") {
      levels(object$group.by)
    } else {
      levels(object$group.by.test)
    }
    ranef <- if (length(dim(object$ranef)) > 2L) {
      object$ranef[,, ranefNames, drop = FALSE]
    } else {
      object$ranef[, ranefNames, drop = FALSE]
    }
    if (n.chains > 1L) {
      if (length(dim(ranef)) > 2L && combineChains) {
        ranef <- combineChains(ranef)
      } else if (length(dim(ranef)) == 2L && !combineChains && n.chains > 1L) {
        ranef <- uncombineChains(ranef, n.chains)
      }
    }

    return(ranef)
  }

  if (sample == "train" && is.null(object[["yhat.train"]])) {
    stop(
      "cannot extract train sample predictions; rbart_vi must be called with 'keepTrainingFits' == TRUE"
    )
  }

  # the "ev" draws are the fit's per-draw location (BART component plus the
  # drawn group intercepts), so the log-likelihood conditions on both
  if (type == "loglik") {
    ev <- extract.rbart(
      object,
      type = "ev",
      sample = "train",
      combineChains = FALSE
    )
    result <- pointwiseLogLikelihood(object, ev)
    return(combineOrUncombineChains(result, n.chains, combineChains))
  }

  result <- if (sample == "train") object$yhat.train else object$yhat.test
  weights <- if (sample == "train") object$weights else object$weights.test
  # if necessary, recover chain information or throw it away
  if (n.chains > 1L) {
    if (length(dim(result)) > 2L && combineChains) {
      result <- combineChains(result)
    } else if (length(dim(result)) == 2L && !combineChains && n.chains > 1L) {
      result <- uncombineChains(result, n.chains)
    }
  }

  #n.samples <- if (length(dim(result)) > 2L) dim(result)[2L] else dim(result)[1L] %/% n.chains
  #n.obs     <- dim(result)[length(dim(result))]

  if (type == "bart") {
    return(result)
  }

  ranefNames <- if (sample == "train") {
    as.character(object$group.by)
  } else {
    as.character(object$group.by.test)
  }
  ranef <- unname(
    if (length(dim(object$ranef)) > 2L) {
      object$ranef[,, ranefNames, drop = FALSE]
    } else {
      object$ranef[, ranefNames, drop = FALSE]
    }
  )

  if (n.chains > 1L) {
    if (length(dim(ranef)) > 2L && combineChains) {
      ranef <- combineChains(ranef)
    } else if (length(dim(ranef)) == 2L && !combineChains && n.chains > 1L) {
      ranef <- uncombineChains(ranef, n.chains)
    }
  }

  result <- result + ranef

  responseIsBinary <- is.null(object[["sigma"]])
  if (responseIsBinary) {
    result <- probabilityFromLatents(result, object)
  }

  if (type == "ppd") {
    result <- sampleFromPPD(result, object, weights, n.chains)
  }

  result
}

fitted.rbart <- function(
  object,
  type = c("ev", "ppd", "bart", "ranef"),
  sample = c("train", "test"),
  ci.level = NULL,
  ...
) {
  if (is.character(type)) {
    if (type[1L] == "response") {
      type[1L] <- "ev"
    } else if (type[1L] == "link") {
      type[1L] <- "bart"
    }
  }
  if (
    !is.character(type) || type[1L] %not_in% eval(formals(fitted.rbart)$type)
  ) {
    stop(
      "type must be in '",
      paste0(eval(formals(fitted.rbart)$type), collapse = "', '"),
      "'"
    )
  }
  type <- type[1L]

  if (
    !is.character(sample) ||
      sample[1L] %not_in% eval(formals(fitted.rbart)$sample)
  ) {
    stop(
      "sample must be in '",
      paste0(eval(formals(fitted.rbart)$sample), collapse = "', '"),
      "'"
    )
  }
  sample <- sample[1L]

  if (sample == "train" && type != "ranef" && is.null(object[["yhat.train"]])) {
    stop(
      "cannot extract train sample predictions; rbart_vi must be called with 'keepTrainingFits' == TRUE"
    )
  }

  # ci.level routes through the draws (extract) rather than the mean-only C
  # fast path below, then summarizes to est + credible band (kind follows type)
  if (!is.null(ci.level)) {
    return(posteriorInterval(extract(object, type, sample, ...), ci.level))
  }

  if (type == "ev") {
    ranefNames <- dimnames(object$ranef)
    ranefNames <- ranefNames[[length(ranefNames)]]
    if (sample == "train") {
      groupByMatch <- match(object$group.by, ranefNames)
      result <- .Call(
        C_rbart_fitted,
        object$yhat.train,
        object$ranef,
        groupByMatch,
        is.null(object[["sigma"]])
      )
    } else {
      groupByMatch <- match(object$group.by.test, ranefNames)
      result <- .Call(
        C_rbart_fitted,
        object$yhat.test,
        object$ranef,
        groupByMatch,
        is.null(object[["sigma"]])
      )
    }
  } else {
    result <- extract(object, type, sample, ...)

    result <- if (!is.null(dim(result))) {
      apply(result, length(dim(result)), mean)
    } else {
      mean(result)
    }
  }

  result
}

residuals.rbart <- function(object, type = "ev", ...) {
  # as residuals.bart: type reaches fitted for link-scale residuals, sample
  # is pinned to the training response
  object$y - fitted.rbart(object, type = type, sample = "train", ...)
}

# fit-level dispatch for the sampler's plotTree method, so a kept bart or
# rbart fit can be plotted directly instead of reaching into $fit; chainNum
# and sampleNum forward only when supplied, since the method detects them by
# their absence
plotTree.dbartsSampler <- function(object, ...) {
  invisible(object$plotTree(...))
}

plotTree.bart <- function(object, treeNum = 1L, chainNum, sampleNum, ...) {
  if (is.null(object[["fit"]])) {
    stop(
      "plotTree requires the trees to be kept: fit with ",
      "keeptrees/keepTrees = TRUE"
    )
  }
  args <- list(treeNum = treeNum, ...)
  if (!missing(chainNum)) {
    args$chainNum <- chainNum
  }
  if (!missing(sampleNum)) {
    args$sampleNum <- sampleNum
  }
  invisible(do.call(object$fit$plotTree, args))
}

plotTree.rbart <- function(
  object,
  treeNum = 1L,
  chainNum = 1L,
  sampleNum,
  ...
) {
  if (is.null(object[["fit"]])) {
    stop(
      "plotTree requires the trees to be kept: fit rbart_vi with ",
      "keepTrees = TRUE"
    )
  }
  n.chains <- if (is.null(object$n.chains)) {
    length(object$fit)
  } else {
    object$n.chains
  }
  chainNum <- as.integer(chainNum)
  if (
    length(chainNum) != 1L ||
      is.na(chainNum) ||
      chainNum < 1L ||
      chainNum > n.chains
  ) {
    stop("chainNum must be a single chain index in [1, ", n.chains, "]")
  }
  # the in-core Gibbs path keeps one multi-chain sampler (select the chain
  # through its own chainNum); the R-loop path keeps one sampler per chain
  singleFit <- length(object$fit) == 1L && n.chains > 1L
  sampler <- if (singleFit) object$fit[[1L]] else object$fit[[chainNum]]
  args <- list(
    treeNum = treeNum,
    chainNum = if (singleFit) chainNum else 1L,
    ...
  )
  if (!missing(sampleNum)) {
    args$sampleNum <- sampleNum
  }
  invisible(do.call(sampler$plotTree, args))
}


# ev (expected value) enters in the caller's requested layout: chains split
# ((n.chains x) n.samples x n.obs, obs last) or chains combined ((n.chains *
# n.samples) x n.obs, chain-blocked rows - all of chain 1's samples, then
# chain 2's). The gaussian noise is always drawn in object$sigma's own
# native chain-fastest order (matching both of its storage layouts, verified
# empirically), so every draw gets its own sigma regardless of the requested
# output shape; when ev is combined, the freshly drawn noise is reshaped -
# not redrawn - with the same combineChains() helper the stored draws
# themselves go through, so a combined and a split ppd draw from the same
# seed agree bit-for-bit after accounting for row order. n.chains is needed
# only to perform that reshape.
sampleFromPPD <- function(ev, object, weights, n.chains = 1L) {
  oldSeed <- NULL
  if (!is.null(object[["seed"]])) {
    oldSeed <- .GlobalEnv$.Random.seed
    .GlobalEnv$.Random.seed <- object$seed
  }

  responseIsBinary <- is.null(object$sigma)

  if (is.null(weights)) {
    if (responseIsBinary) {
      if (length(dim(ev)) > 2L) {
        result <- array(
          rbinom(length(ev), 1L, ev),
          dim(ev),
          dimnames = dimnames(ev)
        )
      } else {
        result <- matrix(
          rbinom(length(ev), 1L, ev),
          nrow(ev),
          ncol(ev),
          dimnames = list(rownames(ev), colnames(ev))
        )
      }
    } else {
      n.obs <- dim(ev)[length(dim(ev))]
      n.draws <- length(object$sigma)
      noise <- rnorm(
        n.obs * n.draws,
        0,
        rep_len(object$sigma, n.obs * n.draws)
      )
      if (n.chains > 1L && length(dim(ev)) < 3L) {
        noise <- combineChains(array(
          noise,
          c(n.chains, n.draws %/% n.chains, n.obs)
        ))
      }
      result <- ev + noise
    }
  } else {
    if (responseIsBinary) {
      # a weight-w row is w iid bernoulli trials; the coherent posterior
      # predictive draw is the number of successes, rbinom(, w, ev), not a
      # bernoulli draw scaled by w. size is recycled to match ev's own
      # column-major fill so each obs's weight lines up with its draws.
      if (length(dim(ev)) > 2L) {
        size <- rep(weights, each = prod(dim(ev)[1L:2L]))
        result <- array(
          rbinom(length(ev), size, ev),
          dim(ev),
          dimnames = dimnames(ev)
        )
      } else {
        size <- rep(weights, each = nrow(ev))
        result <- matrix(
          rbinom(length(ev), size, ev),
          nrow(ev),
          ncol(ev),
          dimnames = list(rownames(ev), colnames(ev))
        )
      }
    } else {
      n.obs <- dim(ev)[length(dim(ev))]
      n.draws <- length(object$sigma)
      sd <- rep_len(object$sigma, n.obs * n.draws) *
        rep(sqrt(1 / weights), each = n.draws)
      noise <- rnorm(n.obs * n.draws, 0, sd)
      if (n.chains > 1L && length(dim(ev)) < 3L) {
        noise <- combineChains(array(
          noise,
          c(n.chains, n.draws %/% n.chains, n.obs)
        ))
      }
      result <- ev + noise
    }
  }
  if (!is.null(oldSeed)) {
    .GlobalEnv$.Random.seed <- oldSeed
  }

  result
}

# family/chain-count/tree-count/burn-in/kept-draws synopsis for print.bart
# and print.rbart, built only from fields that exist regardless of
# keepCall, keepSampler, and (for rbart) whether one sampler is kept per
# chain or a single in-core Gibbs sampler handles them all - so a fit
# created with keepCall = FALSE still prints something useful. n.trees and
# n.burn are only recoverable when the sampler itself was kept (keepTrees/
# keepSampler = TRUE); they are omitted otherwise, since the fit object
# does not retain them on its own.
fitSynopsis <- function(x) {
  fit <- x[["fit"]]
  fitIsList <- !is.null(fit) && is.list(fit)

  n.chains <- if (fitIsList) {
    # rbart: object$n.chains, when present, is authoritative even when the
    # kept sampler list has only one element - the in-core Gibbs path keeps
    # one multi-chain sampler and stores n.chains alongside it (matches the
    # n.chains idiom in predict.rbart/extract.rbart/plotTree.rbart)
    if (is.null(x[["n.chains"]])) length(fit) else x$n.chains
  } else if (!is.null(fit)) {
    fit$control@n.chains
  } else {
    x$n.chains
  }

  control <- if (!is.null(fit)) {
    if (fitIsList) fit[[1L]]$control else fit$control
  } else {
    NULL
  }

  varcountDims <- dim(x[["varcount"]])
  n.kept <- if (!is.null(control)) {
    control@n.samples
  } else if (is.null(varcountDims)) {
    NA_integer_
  } else if (length(varcountDims) == 3L) {
    varcountDims[2L]
  } else if (n.chains > 1L) {
    varcountDims[1L] %/% n.chains
  } else {
    varcountDims[1L]
  }

  cat("family: ", x$family, "\n", sep = "")
  cat("n.chains: ", n.chains, "\n", sep = "")
  if (!is.null(control)) {
    cat("n.trees: ", control@n.trees, "\n", sep = "")
    cat("n.burn: ", control@n.burn, "\n", sep = "")
  }
  if (!is.na(n.kept)) {
    cat("kept draws (per chain): ", n.kept, "\n", sep = "")
  }
  invisible(NULL)
}

print.bart <- function(x, ...) {
  if (!identical(x[["call"]], call("NULL"))) {
    cat(
      "\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n",
      sep = ""
    )
  }
  fitSynopsis(x)
  invisible(x)
}

print.rbart <- function(x, ...) {
  if (!identical(x[["call"]], call("NULL"))) {
    cat(
      "\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n",
      sep = ""
    )
  }
  fitSynopsis(x)
  invisible(x)
}
