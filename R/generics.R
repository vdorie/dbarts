# predict, extract, fitted and for bart and rbart objects

extract <- function(object, ...) UseMethod("extract")

combineOrUncombineChains <- function(x, n.chains, combineChains)
{
  if (n.chains > 1L) {
    if (length(dim(x)) > 2L && combineChains)
      x <- combineChains(x)
    else if (length(dim(x)) == 2L && !combineChains)
      x <- uncombineChains(x, n.chains)
  }
  x
}

predict.bart <- function(object, newdata, offset, weights,
                         type = c("ev", "ppd", "bart"),
                         combineChains = TRUE,
                         ...)
{
  if (missing(offset)) offset <- NULL
  if (missing(weights)) weights <- NULL

  if (is.null(object[["fit"]])) {
    if (as.character(object$call[[1L]]) == "bart2")
      stop("predict requires bart2 to be called with 'keepTrees' == TRUE")
    else
      stop("predict requires bart to be called with 'keeptrees' == TRUE")
  }
  
  if (is.character(type)) {
    if (type[1L] == "response")  type[1L] <- "ev"
    else if (type[1L] == "link") type[1L] <- "bart"
  }
  if (!is.character(type) || length(type) == 0L || type[1L] %not_in% eval(formals(predict.bart)$type))
    stop("type must be in '", paste0(eval(formals(predict.rbart)$type), collapse = "', '"), "'")
  type <- type[1L]
  
  result <- object$fit$predict(newdata, offset)
  # result is n.obs x n.samples x n.chains
  n.chains <- object$fit$control@n.chains
  result <- convertSamplesFromDbartsToBart(result, n.chains, combineChains)
  
  if (type == "bart")
    return(result)

  if ((responseIsBinary <- is.null(object[["sigma"]])))
    result <- pnorm(result)
  
  if (type == "ppd")
    result <- sampleFromPPD(result, object, weights)
  
  result
}

extract.bart <- function(object, 
                         type = c("ev", "ppd", "bart", "trees"),
                         sample = c("train", "test"),
                         combineChains = TRUE,
                         ...)
{
  if (is.character(type)) {
    if (type[1L] == "response") type[1L] <- "ev"
    else if (type[1L] == "link") type[1L] <- "bart"
  }
  if (!is.character(type) || type[1L] %not_in% eval(formals(extract.bart)$type))
    stop("type must be in '", paste0(eval(formals(extract.bart)$type), collapse = "', '"), "'")
  type <- type[1L]
  
  if (type == "trees") {
    if (is.null(object$fit)) {
      if (as.character(object$call[[1L]]) == "bart2")
        stop("extracting trees requires bart2 to be called with 'keepTrees' == TRUE")
      else
        stop("extracting trees requires bart to be called with 'keeptrees' == TRUE")
    }
    treesCall <- match.call()
    target <- quote(object$fit$getTrees)
    target[[2L]][[2L]] <- treesCall$object
    treesCall[[1L]] <- target
    treesCall$object <- NULL
    treesCall$type <- NULL
    return(eval(treesCall, parent.frame()))
  }
  
  if (!is.character(sample) || sample[1L] %not_in% eval(formals(extract.bart)$sample))
    stop("sample must be in '", paste0(eval(formals(extract.bart)$sample), collapse = "', '"), "'")
  sample <- sample[1L]
  
  if (sample == "test" && is.null(object[["yhat.test"]]))
    stop("cannot extract test sample predictions if no test data exists; use `predict` instead")  
  
  result <- if (sample == "train") object$yhat.train else object$yhat.test
  weights <- if (sample == "train") object$weigths else object$weights.test
  
  n.chains  <- if (!is.null(object[["fit"]])) object$fit$control@n.chains else object$n.chains
  #n.samples <- if (length(dim(result)) > 2L) dim(result)[2L] else dim(result)[1L] %/% n.chains
  #n.obs     <- if (length(dim(result)) > 2L) dim(result)[3L] else dim(result)[2L]
  
  result <- combineOrUncombineChains(result, n.chains, combineChains)
    
  if (type == "bart")
    return(result)
  
  if ((responseIsBinary <- is.null(object[["sigma"]])))
    result <- pnorm(result)
  
  if (type == "ppd")
    result <- sampleFromPPD(result, object, weights)
  
  result
}

fitted.bart <- function(object,
                        type = c("ev", "ppd", "bart"),
                        sample = c("train", "test"),
                        ...)
{
  if (is.character(type)) {
    if (type[1L] == "response") type[1L] <- "ev"
    else if (type[1L] == "link") type[1L] <- "bart"
  }
  if (!is.character(type) || type[1L] %not_in% eval(formals(fitted.bart)$type))
    stop("type must be in '", paste0(eval(formals(fitted.bart)$type), collapse = "', '"), "'")
  type <- type[1L]
  
  if (!is.character(sample) || sample[1L] %not_in% eval(formals(fitted.bart)$sample))
    stop("sample must be in '", paste0(eval(formals(fitted.bart)$sample), collapse = "', '"), "'")
  sample <- sample[1L]
  
  result <- extract(object, type, sample, ...)
  
  if (!is.null(dim(result))) apply(result, length(dim(result)), mean) else mean(result)
}

residuals.bart <- function(object, ...) {
  object$y - fitted.bart(object)
}

predict.rbart <- function(object, newdata, group.by, offset,
                          type = c("ev", "ppd", "bart", "ranef"),
                          combineChains = TRUE,
                          ...)
{
  if (is.null(object$fit))
    stop("predict requires rbart to be called with 'keepTrees' == TRUE")
  
  dotsList <- list(...)
  if (!is.null(dotsList[["value"]])) {
    warning("argument 'value' has been deprecated; use 'type' instead")
    type <- dotsList[["value"]]
    dotsList[["value"]] <- NULL
  }
  
  if (is.character(type)) {
    if (type[1L] == "response") type[1L] <- "ev"
    else if (type[1L] == "link") type[1L] <- "bart"
  }
  if (is.character(type) && length(type) > 0L &&  type[1L] == "post-mean") {
    warning("type of 'post-mean' for predict deprecated; use 'ev' instead")
    type[1L] <- "ev"
  }
  if (!is.character(type) || length(type) == 0L || type[1L] %not_in% eval(formals(predict.rbart)$type))
    stop("type must be in '", paste0(eval(formals(predict.rbart)$type), collapse = "', '"), "'")
  type <- type[1L]
  
  if (missing(offset)) offset <- NULL
  
  n.chains  <- if (is.null(object$n.chains)) length(object$fit) else object$n.chains
  n.samples <- object$fit[[1L]]$control@n.samples
    
  nonParametricPart <- 0
  # collects results in an array of n.obs x n.samples x n.chains, default for
  # internal sampler
  #
  # utilize bart stuff to get n.obs, since we would otherwise have to build
  # the test matrix
  if (type != "ranef") {
    if (n.chains > 1L) {
      n.obs <- NULL
      nonParametricPart <- array(sapply(seq_len(n.chains), function(i) {
        res <- object$fit[[i]]$predict(newdata, offset)
          if (is.null(n.obs)) n.obs <<- dim(res)[1L]
          res
        }), c(n.obs, n.samples, n.chains))
    } else {
      nonParametricPart <- object$fit[[1L]]$predict(newdata, offset)
      n.obs <- nrow(nonParametricPart)
    }
    if (n.obs != length(group.by))
      stop("length of group.by not equal to number of rows in test")
    
    nonParametricPart <- convertSamplesFromDbartsToBart(nonParametricPart, n.chains, combineChains)
  }
  
  if (type == "bart") return(nonParametricPart)
  
  ranef <- 0
  if (type != "bart") {
    ranefNames.test  <- levels(group.by)
    ranefNames.train <- if (length(dim(object$ranef)) > 2L) dimnames(object$ranef)[[3L]] else dimnames(object$ranef)[[2L]]
  
    ranef <- object$ranef
    if (n.chains > 1L) {
      if (length(dim(ranef)) > 2L && combineChains)
        ranef <- combineChains(ranef)
      else if (length(dim(ranef)) == 2L && !combineChains && n.chains > 1L)
        ranef <- uncombineChains(ranef, n.chains)
    }
    
    if (!all(measuredLevels <- ranefNames.test %in% ranefNames.train)) {
      warning("test includes random effect levels not present in training - ranef estimates default to draws from their latent distribution parameterized by the posterior of its variance; draws may not be the same across future calls to 'predict'")
      n.unmeasured <- sum(!measuredLevels)
      if (n.chains > 1L) {
        if (!combineChains) {
          unmeasuredRanef <- array(rnorm(n.chains * n.samples * n.unmeasured, 0, rep.int(object$tau, n.unmeasured)),
                                   c(n.chains, n.samples, n.unmeasured),
                                   dimnames = list(NULL, NULL, ranefNames.test[!measuredLevels]))
        } else {
          unmeasuredRanef <- matrix(rnorm(n.chains * n.samples * n.unmeasured, 0, rep.int(object$tau, n.unmeasured)),
                                    n.chains * n.samples, n.unmeasured,
                                    dimnames = list(NULL, ranefNames.test[!measuredLevels]))
        }
        if (length(dim(object$ranef)) == 2L) {
          ranef <- cbind(ranef, unmeasuredRanef)
        } else {
          # ranef are n.chains x n.samples x n.group
          ranef <- array(c(ranef, unmeasuredRanef), c(n.chains, n.samples, dim(ranef)[3L] + n.unmeasured),

                         dimnames = list(NULL, NULL, c(dimnames(ranef)[[3L]], dimnames(unmeasuredRanef)[[3L]])))
        }
      } else {
        unmeasuredRanef <- matrix(rnorm(n.samples * n.unmeasured, 0, rep.int(object$tau, n.unmeasured)),
                                  n.samples, n.unmeasured,
                                  dimnames = list(NULL, ranefNames.test[!measuredLevels]))
        ranef <- cbind(ranef, unmeasuredRanef)
      }
    }
  }
  
  if (type == "ranef") {
    ranef <- if (length(dim(ranef)) > 2L) ranef[,,ranefNames.test,drop = FALSE] else ranef[,ranefNames.test,drop = FALSE]
    return(combineOrUncombineChains(ranef, n.chains, combineChains))
  }
  
  ranef <- unname(if (length(dim(ranef)) > 2L) ranef[,,as.character(group.by),drop = FALSE] else ranef[,as.character(group.by),drop = FALSE])
  ranef <- combineOrUncombineChains(ranef, n.chains, combineChains)
  
  if (length(dim(nonParametricPart)) != length(dim(ranef)) || any(dim(nonParametricPart) != dim(ranef)))
    browser()
  result <- nonParametricPart + ranef
  
  responseIsBinary <- is.null(object[["sigma"]])
  if (responseIsBinary) result <- pnorm(result)
  
  if (type == "ppd")
    result <- sampleFromPPD(result, object, NULL)
  
  if (exists("unmeasuredRanef", inherits = FALSE)) attr(result, "ranef") <- unmeasuredRanef
  
  result
}

extract.rbart <- function(object,
                          type = c("ev", "ppd", "bart", "ranef", "trees"),
                          sample = c("train", "test"),
                          combineChains = TRUE,
                          ...)
{
  if (is.character(type)) {
    if (type[1L] == "response") type[1L] <- "ev"
    else if (type[1L] == "link") type[1L] <- "bart"
  }
  if (!is.character(type) || type[1L] %not_in% eval(formals(extract.rbart)$type))
    stop("type must be in '", paste0(eval(formals(extract.rbart)$type), collapse = "', '"), "'")
  type <- type[1L]

  n.chains  <- if (is.null(object$n.chains)) length(object$fit) else object$n.chains
  
  if (type == "trees") {
    if (is.null(object$fit))
      stop("extracting trees requires rbart to be called with 'keepTrees' == TRUE")
    treesCall <- match.call()
    target <- quote(object$fit[[i]]$getTrees)
    target[[2L]][[2L]][[2L]] <- treesCall$object
    treesCall[[1L]] <- target
    treesCall$object <- NULL
    treesCall$type <- NULL
    treesCall$chainNums <- NULL
    evalEnv <- parent.frame()
    dotsList <- list(...)
    chainNums <- if ("chainNums" %in% names(dotsList)) as.integer(dotsList[["chainNums"]]) else seq_len(n.chains)
    varOrder <- c("sample", "chain", "tree", "n", "var", "value")
    allTrees <- lapply(chainNums, function(i) {
      result_i <- eval(subTermInLanguage(treesCall, quote(i), i), evalEnv)
      if (n.chains > 1L) result_i$chain <- i
      result_i[,match(varOrder, colnames(result_i))]
    })
    if (length(allTrees) > 1L) {
      allTrees <- Reduce(rbind, allTrees)
    } else {
      allTrees <- allTrees[[1L]]
    }
    row.names(allTrees) <- as.character(seq_len(nrow(allTrees)))
    return(allTrees)
  }

  if (!is.character(sample) || sample[1L] %not_in% eval(formals(extract.rbart)$sample))
    stop("sample must be in '", paste0(eval(formals(extract.rbart)$sample), collapse = "', '"), "'")
  sample <- sample[1L]
  
  if (sample == "test" && is.null(object[["yhat.test"]]))
    stop("cannot extract test sample predictions if no test data exists; use `predict` instead")
  
    
  if (type == "ranef") {
    ranefNames <- if (sample == "train") levels(object$group.by) else levels(object$group.by.test)
    ranef <- if (length(dim(object$ranef)) > 2L) object$ranef[,,ranefNames,drop = FALSE] else object$ranef[,ranefNames,drop = FALSE]
    if (n.chains > 1L) {
      if (length(dim(ranef)) > 2L && combineChains)
        ranef <- combineChains(ranef)
      else if (length(dim(ranef)) == 2L && !combineChains && n.chains > 1L)
        ranef <- uncombineChains(ranef, n.chains)
    }
    
    return(ranef)
  }  
  
  result <- if (sample == "train") object$yhat.train else object$yhat.test
  # if necessary, recover chain information or throw it away
  if (n.chains > 1L) {
    if (length(dim(result)) > 2L && combineChains)
      result <- combineChains(result)
    else if (length(dim(result)) == 2L && !combineChains && n.chains > 1L)
      result <- uncombineChains(result, n.chains)
  }
  
  #n.samples <- if (length(dim(result)) > 2L) dim(result)[2L] else dim(result)[1L] %/% n.chains
  #n.obs     <- dim(result)[length(dim(result))]
  
  if (type == "bart") return(result)
  
  ranefNames <- if (sample == "train") as.character(object$group.by) else as.character(object$group.by.test)
  ranef <- unname(if (length(dim(object$ranef)) > 2L) object$ranef[,,ranefNames,drop = FALSE] else object$ranef[,ranefNames,drop = FALSE])
  
  if (n.chains > 1L) {
    if (length(dim(ranef)) > 2L && combineChains)
      ranef <- combineChains(ranef)
    else if (length(dim(ranef)) == 2L && !combineChains && n.chains > 1L)
      ranef <- uncombineChains(ranef, n.chains)
  }
  
  result <- result + ranef
  
  responseIsBinary <- is.null(object[["sigma"]])
  if (responseIsBinary) result <- pnorm(result)
  
  if (type == "ppd")
    result <- sampleFromPPD(result, object, NULL)
  
  result
}

fitted.rbart <- function(object,
                         type = c("ev", "ppd", "bart", "ranef"),
                         sample = c("train", "test"),
                         ...)
{
  if (is.character(type)) {
    if (type[1L] == "response") type[1L] <- "ev"
    else if (type[1L] == "link") type[1L] <- "bart"
  }
  if (!is.character(type) || type[1L] %not_in% eval(formals(fitted.rbart)$type))
    stop("type must be in '", paste0(eval(formals(fitted.rbart)$type), collapse = "', '"), "'")
  type <- type[1L]
  
  if (!is.character(sample) || sample[1L] %not_in% eval(formals(fitted.rbart)$sample))
    stop("sample must be in '", paste0(eval(formals(fitted.rbart)$sample), collapse = "', '"), "'")
  sample <- sample[1L]
  
  if (type == "ev") {
    ranefNames <- dimnames(object$ranef)
    ranefNames <- ranefNames[[length(ranefNames)]]
    if (sample == "train") {
      groupByMatch <- match(object$group.by, ranefNames)
      result <- .Call(C_rbart_fitted, object$yhat.train, object$ranef, groupByMatch, is.null(object[["sigma"]]))
    } else {
      groupByMatch <- match(object$group.by.test, ranefNames)
      result <- .Call(C_rbart_fitted, object$yhat.test, object$ranef, groupByMatch, is.null(object[["sigma"]]))
    }
  } else {
    result <- extract(object, type, sample, ...)
    
    result <- if (!is.null(dim(result))) apply(result, length(dim(result)), mean) else mean(result)
  }

  result
}

residuals.rbart <- function(object, ...) {
  object$y - fitted.rbart(object)
}


# NOTE: this is outdated
# ev (expected value) should have dimensions
#   n.samples x n.chains x n.obs, (n.samples * n.chains n.obs),
#   or n.samples x n.obs if n.chains = 1
# 
# ev consists of contiguous blocks of length equal to the number of 
# samples, so that ev[1:totalNumSamples] should get paired with
# as.vector(sigma), and then repeated from there
#
#
# for ev of dim n.chains x n.samples x n.obs (bart default),
# each sigma needs to be repeated as below
sampleFromPPD <- function(ev, object, weights)
{
  oldSeed <- NULL
  if (!is.null(object[["seed"]])) {
    oldSeed <- .GlobalEnv$.Random.seed
    .GlobalEnv$.Random.seed <- object$seed
  }
  
  responseIsBinary <- is.null(object$sigma)
  
  if (is.null(weights)) {
    if (responseIsBinary) {
      if (length(dim(ev)) > 2L)
        result <- array(rbinom(length(ev), 1L, ev), dim(ev), dimnames = dimnames(ev))
      else
        result <- matrix(rbinom(length(ev), 1L, ev), nrow(ev), ncol(ev), dimnames = list(rownames(ev), colnames(ev)))
    } else {
      n.obs <- dim(ev)[length(dim(ev))]
      result <- ev + rnorm(n.obs * length(object$sigma), 0, rep_len(object$sigma, n.obs * length(object$sigma)))
    }
  } else {
    if (responseIsBinary) {
      if (length(dim(ev)) > 2L) {
        result <- array(rbinom(length(ev), 1L, ev), dim(ev), dimnames = dimnames(ev))
        # recycle weight vector by permuting observations to first dimension
        result <- aperm(weights * aperm(result, c(3L, 1L, 2L)), c(2L, 3L, 1L))
      } else {
        result <- matrix(rbinom(length(ev), 1L, ev), nrow(ev), ncol(ev), dimnames = list(rownames(ev), colnames(ev)))
        result <- t(weights * t(result))
      }
    } else {
      n.obs <- dim(ev)[length(dim(ev))]
      n.samples <- length(object$sigma)
      sigma <- rep_len(object$sigma, n.obs * n.samples) * rep(sqrt(1 / weights), each = n.samples)
      result <- ev + rnorm(n.obs * n.samples, 0, sigma)
    }
  }
  if (!is.null(oldSeed))
    .GlobalEnv$.Random.seed <- oldSeed
  
  result
}

print.bart <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  invisible(x)
}

print.rbart <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  invisible(x)
}

