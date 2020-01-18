rbart.priors <- list(cauchy = function(x, rel.scale) dcauchy(x, 0, rel.scale * 2.5, TRUE),
                     gamma  = function(x, rel.scale) dgamma(x, shape = 2.5, scale = rel.scale * 2.5, log = TRUE))
cauchy <- NULL ## for R CMD check

rbart_vi <- function(
  formula, data, test, subset, weights, offset, offset.test = offset,
  group.by, group.by.test, prior = cauchy, ## can be a symbol in rbart.priors or a function; on log scale
  sigest = NA_real_, sigdf = 3.0, sigquant = 0.90,
  k = 2.0,
  power = 2.0, base = 0.95,
  n.trees = 75L,
  n.samples = 1500L, n.burn = 1500L,
  n.chains = 4L, n.threads = min(guessNumCores(), n.chains), combineChains = FALSE,
  n.cuts = 100L, useQuantiles = FALSE,
  n.thin = 5L, keepTrainingFits = TRUE,
  printEvery = 100L, printCutoffs = 0L,
  verbose = TRUE,
  keepTrees = TRUE, keepCall = TRUE, ...)
{
  matchedCall <- match.call()
  callingEnv <- parent.frame()
  
  # because we use a lot of trickery to redirect calls in the calling environment
  # (for example, to get the data), we replicate some base mechanisms like complaining
  # about unknown arguments
  argNames <- names(matchedCall)[-1L]
  unknownArgs <- argNames %not_in% names(formals(rbart_vi)) & argNames %not_in% names(formals(dbartsControl))
  if (any(unknownArgs))
    stop("unknown arguments: '", paste0(argNames[unknownArgs], collapse = "', '"), "'")
  
  n.chains <- coerceOrError(n.chains, "integer")[1L]
  if (is.na(n.chains) || n.chains < 1L) stop("n.chains must be a non-negative integer")
  
  n.threads <- coerceOrError(n.threads, "integer")[1L]
  if (is.na(n.threads) || n.threads < 1L) stop("n.threads must be a non-negative integer")
  
  controlCall <- redirectCall(matchedCall, dbarts::dbartsControl)
  control <- eval(controlCall, envir = callingEnv)
  
  control@call <- if (keepCall) matchedCall else call("NULL")
  control@n.burn     <- control@n.burn     %/% control@n.thin
  control@n.samples  <- control@n.samples  %/% control@n.thin
  control@printEvery <- control@printEvery %/% control@n.thin
  control@n.chains <- 1L
  control@n.threads <- max(control@n.threads %/% n.chains, 1L)
  if (n.chains > 1L && n.threads > 1L) {
    if (control@verbose) warning("verbose output disabled for multiple threads")
    control@verbose <- FALSE
  }
  
  tree.prior <- quote(cgm(power, base))
  tree.prior[[2L]] <- power; tree.prior[[3L]] <- base

  if (!is.null(matchedCall[["k"]])) {
    node.prior <- quote(normal(k))
    node.prior[[2L]] <- matchedCall[["k"]]
  } else {
    node.prior <- NULL
  }

  resid.prior <- quote(chisq(sigdf, sigquant))
  resid.prior[[2L]] <- sigdf; resid.prior[[3L]] <- sigquant
    
  group.by <- tryCatch(group.by, error = function(e) e)
  if ((is(group.by, "error") || !(is.numeric(group.by) || is.factor(group.by) || is.character(group.by))) && is.language(formula) && formula[[1L]] == '~') {
    if (is.symbol(matchedCall[["group.by"]]) && any(names(data) == matchedCall[["group.by"]])) {
      group.by <- data[[which(names(data) == matchedCall[["group.by"]])[1L]]]
    } else {
      group.by <- eval(matchedCall[["group.by"]], environment(formula))
    }
  }
  if (is(group.by, "error"))
    stop("'group.by' not found")
  if (!is.numeric(group.by) && !is.factor(group.by) && !is.character(group.by))
    stop("'group.by' must be coercible to factor type")
  
  if (!is.null(matchedCall[["group.by.test"]])) {
    group.by.test <- tryCatch(group.by.test, error = function(e) e)
    if ((is(group.by.test, "error") || !(is.numeric(group.by.test) || is.factor(group.by.test) || is.character(group.by.test))) && is.language(formula) && formula[[1L]] == '~') {
      if (is.symbol(matchedCall[["group.by.test"]]) && any(names(data) == matchedCall[["group.by.test"]])) {
        group.by.test <- data[[which(names(data) == matchedCall[["group.by.test"]])[1L]]]
      } else {
        group.by.test <- eval(matchedCall[["group.by.test"]], environment(formula))
      }
    }
    if (is(group.by.test, "error"))
      stop("'group.by.test' specified but not found")
    if (!is.numeric(group.by.test) && !is.factor(group.by.test) && !is.character(group.by.test))
      stop("'group.by.test' must be coercible to factor type")
  }
  
  
  if (is.null(matchedCall$prior)) matchedCall$prior <- formals(rbart_vi)$prior
  
  if (is.symbol(matchedCall$prior) || is.character(matchedCall$prior) && any(names(rbart.priors) == matchedCall$prior))
    prior <- rbart.priors[[which(names(rbart.priors) == matchedCall$prior)]]
  
  data <- eval(redirectCall(matchedCall, dbarts::dbartsData), envir = callingEnv)
  
  if (length(group.by) != length(data@y))
    stop("'group.by' not of length equal to that of data")
  group.by <- droplevels(as.factor(group.by))
  if (!is.null(matchedCall[["group.by.test"]])) {
    if (length(group.by.test) != nrow(data@x.test))
      stop("'group.by.test' not of length equal to that of data")
    group.by.test <- droplevels(as.factor(group.by.test))
  } else {
    group.by.test <- NULL
  }
  
  samplerArgs <- namedList(formula = data, control, tree.prior, node.prior, resid.prior,
                           sigma = as.numeric(sigest))
  if (is.null(node.prior)) samplerArgs[["node.prior"]] <- NULL

  chainResults <- vector("list", n.chains)
  if (n.threads == 1L || n.chains == 1L) {
    for (chainNum in seq_len(n.chains))
      chainResults[[chainNum]] <- rbart_vi_fit(samplerArgs, group.by, prior)
  } else {
    cluster <- makeCluster(n.threads)
    
    clusterExport(cluster, "rbart_vi_fit", asNamespace("dbarts"))
    clusterEvalQ(cluster, require(dbarts))
    
    tryResult <- tryCatch(
      chainResults <- clusterCall(cluster, "rbart_vi_fit", samplerArgs, group.by, prior))
    
    stopCluster(cluster)
  }
  
  packageRbartResults(control, data, group.by, group.by.test, chainResults, combineChains)
}

rbart_vi_fit <- function(samplerArgs, group.by, prior) 
{
  sampler <- do.call(dbarts::dbarts, samplerArgs)
  sampler$control@call <- samplerArgs$control@call
  
  oldUpdateState <- sampler$control@updateState
  verbose <- sampler$control@verbose
  control <- sampler$control
  control@updateState <- FALSE
  control@verbose <- FALSE
  sampler$setControl(control)
  
  y <- sampler$data@y
  rel.scale <- if (!control@binary) sd(y) else 0.50
  
  g <- as.integer(group.by)
  numRanef <- nlevels(group.by)
  g.sel <- lapply(seq_len(numRanef), function(j) g == j)
  n.g <- sapply(g.sel, sum)
  
  evalEnv <- list2env(list(rel.scale = rel.scale, q = numRanef))
  b.sq <- NULL ## for R CMD check
  posteriorClosure <- function(x) { ifelse(x <= 0.0 | is.infinite(x), -.Machine$double.xmax * .Machine$double.eps, -q * base::log(x) - 0.5 * b.sq / x^2.0 + prior(x, rel.scale)) }
  environment(posteriorClosure) <- evalEnv
  
  offset.orig <- sampler$data@offset
  
  numObservations <- length(sampler$data@y)
  numTestObservations <- NROW(sampler$data@x.test)
  
  firstTau   <- rep(NA_real_, control@n.burn)
  firstSigma <- if (!control@binary) rep(NA_real_, control@n.burn) else NULL
  tau   <- rep(NA_real_, control@n.samples)
  sigma <- if (!control@binary) rep(NA_real_, control@n.samples) else NULL
  ranef <- matrix(NA_real_, numRanef, control@n.samples)
  yhat.train <- matrix(NA_real_, numObservations, control@n.samples)
  yhat.test  <- matrix(NA_real_, numTestObservations, control@n.samples)
  varcount <- matrix(NA_integer_, ncol(sampler$data@x), control@n.samples)
  
  sampler$sampleTreesFromPrior()
  y.st <- if (!control@binary) y else sampler$getLatents()
  
  tau.i <- rel.scale / 5.0
  ranef.i <- rnorm(numRanef, 0.0, tau.i)
  ranef.vec <- ranef.i[g]
  sigma.i <- if (!control@binary) sampler$data@sigma else 1.0
  
  if (control@n.burn > 0L) {
    oldKeepTrees <- control@keepTrees
    control@keepTrees <- FALSE
    sampler$setControl(control)
    
    sampler$setOffset(ranef.vec + if (!is.null(offset.orig)) offset.orig else 0)
    treeFit.train <- sampler$predict(sampler$data@x) - ranef.vec
    
    # order of update matters - need to store a ranef that goes with a prediction
    # or else when they're added together they won't be consistent with `predict'
    for (i in seq_len(control@n.burn)) {
      # update ranef
      resid <- y.st - treeFit.train
      post.var <- 1.0 / (n.g / sigma.i^2.0 + 1.0 / tau.i^2.0)
      post.mean <- (n.g / sigma.i^2.0) * sapply(seq_len(numRanef), function(j) mean(resid[g.sel[[j]]])) * post.var
      
      ranef.i <- rnorm(numRanef, post.mean, sqrt(post.var))
      ranef.vec <- ranef.i[g]
      
      # update BART params
      sampler$setOffset(ranef.vec + if (!is.null(offset.orig)) offset.orig else 0)
      samples <- sampler$run(0L, 1L)
      treeFit.train <- as.vector(samples$train) - ranef.vec
      if (control@binary) sampler$getLatents(y.st)
      sigma.i <- samples$sigma[1L]
      
      # update sd of ranef
      evalEnv$b.sq <- sum(ranef.i^2.0)
      tau.i <- sliceSample(posteriorClosure, tau.i, control@n.thin, boundary = c(0.0, Inf))[control@n.thin]
      
      .Call(C_dbarts_assignInPlace, firstTau, i, tau.i)
      if (!control@binary)
        .Call(C_dbarts_assignInPlace, firstSigma, i, sigma.i)
    }
    
    if (control@keepTrees != oldKeepTrees) {
      control@keepTrees <- TRUE
      sampler$setControl(control)
    }
  } else {
    sampler$setOffset(ranef.vec + if (!is.null(offset.orig)) offset.orig else 0)
    treeFit.train <- (if (control@n.samples > 1L && control@keepTrees) sampler$predict(sampler$data@x)[,1L] else sampler$predict(sampler$data@x)) - ranef.vec
  }
  
  for (i in seq_len(control@n.samples)) {
    # update ranef
    resid <- y.st - treeFit.train
    post.var <- 1.0 / (n.g / sigma.i^2.0 + 1.0 / tau.i^2.0)
    post.mean <- (n.g / sigma.i^2.0) * sapply(seq_len(numRanef), function(j) mean(resid[g.sel[[j]]])) * post.var
    ranef.i <- rnorm(numRanef, post.mean, sqrt(post.var))
    ranef.vec <- ranef.i[g]
      
    # update BART params
    sampler$setOffset(ranef.vec + if (!is.null(offset.orig)) offset.orig else 0)
    samples <- sampler$run(0L, 1L)
    treeFit.train <- as.vector(samples$train) - ranef.vec
    if (control@binary) sampler$getLatents(y.st)
    sigma.i <- samples$sigma[1L]
    
    # update sd of ranef
    evalEnv$b.sq <- sum(ranef.i^2.0)
    tau.i <- sliceSample(posteriorClosure, tau.i, control@n.thin, boundary = c(0.0, Inf))[control@n.thin]
        
    .Call(C_dbarts_assignInPlace, tau, i, tau.i)
    if (!control@binary)
      .Call(C_dbarts_assignInPlace, sigma, i, sigma.i)
    .Call(C_dbarts_assignInPlace, ranef, i, ranef.i)
    .Call(C_dbarts_assignInPlace, yhat.train, i, treeFit.train)
    .Call(C_dbarts_assignInPlace, varcount, i, samples$varcount)
    if (numTestObservations > 0L) .Call(C_dbarts_assignInPlace, yhat.test, i, samples$test)
    
    if (verbose && i %% control@printEvery == 0L) cat("iter: ", i, "\n", sep = "")
  }
  
  sampler$setOffset(if (!is.null(offset.orig)) offset.orig else NULL)
  
  control@updateState <- oldUpdateState
  sampler$setControl(control)
  
  rownames(ranef) <- levels(group.by)
  
  namedList(sampler, ranef, firstTau, firstSigma, tau, sigma, yhat.train, yhat.test, varcount)
}

packageRbartResults <- function(control, data, group.by, group.by.test, chainResults, combineChains)
{
  n.chains <- length(chainResults)
  
  responseIsBinary <- chainResults[[1L]]$sampler$control@binary
  
  result <- list(call = control@call, y = data@y, group.by = group.by)
  if (!responseIsBinary) result$sigest <- chainResults[[1L]]$sampler$data@sigma
  if (!is.null(group.by.test)) result[["group.by.test"]] <- group.by.test
  
  if (n.chains > 1L) {
    if (!is.null(group.by.test) && any(unmeasuredLevels <- levels(group.by.test) %not_in% levels(group.by))) {
      warning("test includes random effect levels not present in training - ranef estimates default to draws from the ranef distribution parameterized by the posterior of its variance")
      n.samples <- dim(chainResults[[1L]]$ranef)[2L]
      n.unmeasured <- sum(unmeasuredLevels)
      totalRanef <- sapply(seq_along(chainResults), function(k) {
        unmeasuredRanef <- matrix(rnorm(n.unmeasured * n.samples, 0, rep(chainResults[[k]]$tau, rep_len(n.samples, n.unmeasured))),
                                  n.unmeasured, n.samples, dimnames = list(levels(group.by.test)[unmeasuredLevels], NULL))
        rbind(chainResults[[k]]$ranef, unmeasuredRanef)
      })
      ranefDim <- c(dim(chainResults[[1L]]$ranef)[1L] + n.unmeasured, n.samples, n.chains)
      ranefDimnames <- list(c(rownames(chainResults[[1L]]$ranef), levels(group.by.test)[unmeasuredLevels]), NULL, NULL)
      ranef <- array(totalRanef, ranefDim, ranefDimnames)
      result$ranef <- packageSamples(n.chains, combineChains, ranef)
    } else {
      ranef <- array(sapply(chainResults, function(x) x$ranef),
                     c(dim(chainResults[[1L]]$ranef), n.chains),
                     list(rownames(chainResults[[1L]]$ranef), NULL, NULL))
      result$ranef <- packageSamples(n.chains, combineChains, ranef)
      # undo collapse: temp <- aperm(array(result$ranef, c(dim(chainResults[[1L]]$ranef)[2L], n.chains, dim(chainResults[[1L]]$ranef)[1L]), dimnames = list(NULL, NULL, colnames(result$ranef))), c(3L, 1L, 2L))
    }
    result$first.tau   <- packageSamples(n.chains, combineChains, sapply(chainResults, function(x) x$firstTau))
    if (!responseIsBinary) {
      result$first.sigma <- packageSamples(n.chains, combineChains, sapply(chainResults, function(x) x$firstSigma))
      result$sigma       <- packageSamples(n.chains, combineChains, sapply(chainResults, function(x) x$sigma))
     }
    result$tau         <- packageSamples(n.chains, combineChains, sapply(chainResults, function(x) x$tau))
    result$yhat.train  <- packageSamples(n.chains, combineChains,
                                         array(sapply(chainResults, function(x) x$yhat.train), c(dim(chainResults[[1L]]$yhat.train), n.chains)))
    result$yhat.test   <- if (NROW(chainResults[[1L]]$yhat.test) <= 0L) NULL else
                            packageSamples(n.chains, combineChains,
                                           array(sapply(chainResults, function(x) x$yhat.test), c(dim(chainResults[[1L]]$yhat.test), n.chains)))
    result$varcount    <- packageSamples(n.chains, combineChains, array(sapply(chainResults, function(x) x$varcount), c(dim(chainResults[[1L]]$varcount), n.chains)))
  } else {
    result$ranef <- t(chainResults[[1L]]$ranef)
    if (!is.null(group.by.test) && any(unmeasuredLevels <- levels(group.by.test) %not_in% levels(group.by))) {
      warning("test includes random effect levels not present in training - ranef estimates default to draws from the ranef distribution parameterized by the posterior of its variance")
      n.unmeasured <- sum(unmeasuredLevels)
      n.samples <- ncol(chainResults[[1L]]$ranef)
      unmeasuredRanef <-
        matrix(rnorm(n.samples * n.unmeasured, 0, rep(chainResults[[1L]]$tau, n.samples)),
               n.samples, n.unmeasured, dimnames = list(NULL, levels(group.by.test)[unmeasuredLevels]))
      result$ranef <- cbind(result$ranef, unmeasuredRanef)
    }
    result$first.tau   <- chainResults[[1L]]$firstTau
    if (!responseIsBinary) {
      result$first.sigma <- chainResults[[1L]]$firstSigma
      result$sigma       <- chainResults[[1L]]$sigma
    }
    result$tau         <- chainResults[[1L]]$tau
    result$yhat.train  <- t(chainResults[[1L]]$yhat.train)
    result$yhat.test   <- if (NROW(chainResults[[1L]]$yhat.test) <= 0L) NULL else t(chainResults[[1L]]$yhat.test)
    result$varcount    <- chainResults[[1L]]$varcount
  }
  
  result$ranef.mean <- apply(result$ranef, length(dim(result$ranef)), mean)
  result$yhat.train.mean <- apply(result$yhat.train, length(dim(result$yhat.train)), mean)
  if (!is.null(result$yhat.test)) result$yhat.test.mean <- apply(result$yhat.test, length(dim(result$yhat.test)), mean)
  
  if (control@keepTrees == TRUE)
    result$fit <- lapply(chainResults, function(x) x$sampler)
  else
    result$n.chains <- n.chains
  
  if (!exists(".Random.seed", .GlobalEnv)) runif(1L)
  result$seed <- .GlobalEnv[[".Random.seed"]]
  
  class(result) <- "rbart"
  result
}

predict.rbart <- function(object, newdata, group.by, offset,
                          value = c("ev", "ppd", "bart", "ranef"),
                          combineChains = TRUE,
                          ...)
{
  if (is.null(object$fit))
    stop("predict requires rbart to be called with 'keepTrees' == TRUE")
  if (is.character(value) && length(value) > 0L &&  value[1L] == "post-mean") {
    warning("value of 'post-mean' for predict deprecate; use 'ev' instead")
    value[1L] <- "ev"
  }
  if (!is.character(value) || length(value) == 0L || value[1L] %not_in% eval(formals(predict.rbart)$value))
    stop("value must be in '", paste0(eval(formals(predict.rbart)$value), collapse = "', '"), "'")
  value <- value[1L]
  
  if (missing(offset)) offset <- NULL
  
  n.chains <- length(object$fit)
  n.samples <- object$fit[[1L]]$control@n.samples
  
  packageResult <- function(combineChains, result) {
    if (combineChains && length(dim(result)) > 2L) {
      res <- matrix(aperm(result, c(2L, 1L, 3L)), c(prod(dim(result)[1L:2L]), dim(result[3L])))
      if (!is.null(dimnames(result))) colnames(res) <- dimnames(result)[[3L]]
      res
    } else
      result
  }
  
  nonParametricPart <- 0
  if (value != "ranef") {
    if (n.chains > 1L) {
      n.obs <- NULL
      nonParametricPart <- aperm(array(sapply(seq_len(n.chains), function(i) {
        res <- object$fit[[i]]$predict(newdata, offset)
          if (is.null(n.obs)) n.obs <<- dim(res)[1L]
          res
        }), c(n.obs, n.samples, n.chains)), c(3L, 2L, 1L))
    } else {
      nonParametricPart <- t(object$fit[[1L]]$predict(newdata, offset))
      n.obs <- ncol(nonParametricPart)
    }
    if (n.obs != length(group.by))
      stop("length of group.by not equal to number of rows in test")
  }
  
  if (value == "bart") return(packageResult(combineChains, nonParametricPart))
  
  ranef <- 0
  if (value != "bart") {
    ranefNames.test  <- levels(group.by)
    ranefNames.train <- if (length(dim(object$ranef)) > 2L) dimnames(object$ranef)[[3L]] else  dimnames(object$ranef)[[2L]]
  
    if (n.chains > 1L && length(dim(object$ranef)) == 2L && !combineChains) {
      # if the original call collapsed chains, we have to add back in that information
      ranef <- aperm(array(object$ranef, c(n.samples, n.chains, ncol(object$ranef)), dimnames = list(NULL, NULL, colnames(object$ranef))), c(2L, 1L, 3L))
    } else {
      ranef <- object$ranef
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
          ranef <- array(c(ranef, unmeasuredRanef), c(n.chains, n.samples, dim(object$ranef)[3L] + n.unmeasured),
                       dimnames = list(NULL, NULL, c(dimnames(object$ranef)[[3L]], dimnames(unmeasuredRanef)[[3L]])))
        }
      } else {
        unmeasuredRanef <- matrix(rnorm(n.samples * n.unmeasured, 0, rep.int(object$tau, n.unmeasured)),
                                  n.samples, n.unmeasured,
                                  dimnames = list(NULL, ranefNames.test[!measuredLevels]))
        ranef <- cbind(ranef, unmeasuredRanef)
      }
    }
  }
  
  if (value == "ranef") return(packageResult(combineChains, if (length(dim(ranef)) > 2L) ranef[,,ranefNames.test,drop = FALSE] else ranef[,ranefNames.test,drop = FALSE]))
  
  result <- 
    packageResult(combineChains, nonParametricPart) +
      packageResult(combineChains, unname(if (length(dim(ranef)) > 2L) ranef[,,as.character(group.by),drop = FALSE] else ranef[,as.character(group.by),drop = FALSE]))
  if (exists("unmeasuredRanef")) attr(result, "ranef") <- unmeasuredRanef
  result
}

extract.rbart <- function(object,
                          value = c("ev", "ppd", "bart", "ranef"),
                          sample = c("train", "test"),
                          combineChains = TRUE,
                          ...)
{
  if (!is.character(value) || value[1L] %not_in% eval(formals(extract.rbart)$value))
    stop("value must be in '", paste0(eval(formals(extract.rbart)$value), collapse = "', '"), "'")
  value <- value[1L]
  
  if (!is.character(sample) || sample[1L] %not_in% eval(formals(extract.rbart)$sample))
    stop("sample must be in '", paste0(eval(formals(extract.rbart)$sample), collapse = "', '"), "'")
  sample <- sample[1L]
  
  if (sample == "test" && is.null(object[["yhat.test"]]))
    stop("cannot extract test sample predictions if no test data exists; use `predict` instead")
  
  n.chains <- if (is.null(object$n.chains)) length(object$fit) else object$n.chains
  
  if (value == "ranef") {
    ranefNames <- if (sample == "train") levels(object$group.by) else levels(object$group.by.test)
    ranef <- if (length(dim(object$ranef)) > 2L) object$ranef[,,ranefNames] else object$ranef[,ranefNames]
    return(ranef)
  }
  
  if (value == "bart")
    return(if (sample == "train") object$yhat.train else object$yhat.test)
  
  
  nonParametricPart <- if (sample == "train") object$yhat.train else object$yhat.test
  
  ranefNames <- if (sample == "train") as.character(object$group.by) else as.character(object$group.by.test)
  ranef <- if (length(dim(object$ranef)) > 2L) object$ranef[,,ranefNames] else object$ranef[,ranefNames]
  
  result <- nonParametricPart + unname(ranef)
  
  n.samples <- if (length(dim(result)) > 2L) dim(result)[2L] else dim(result)[1L] %/% n.chains
  n.obs     <- if (length(dim(result)) > 2L) dim(result)[3L] else dim(result)[2L]
  
  if (length(dim(result)) == 2L && !combineChains && n.chains > 1L)
    result <- aperm(array(result, c(n.samples, n.chains, n.obs)), c(2L, 1L, 3L))
  
  responseIsBinary <- is.null(object[["sigma"]])
  if (responseIsBinary) result <- pnorm(result)
  
  if (value == "ppd") {
    oldSeed <- .GlobalEnv$.Random.seed
    .GlobalEnv$.Random.seed <- object$seed
    
    if (responseIsBinary) {
      if (length(dim(result)) > 2L)
        result <- array(rbinom(length(result), 1L, result), dim(result))
      else
        result <- matrix(rbinom(length(result), 1L, result), nrow(result), ncol(result))
    } else {
      # result is n.chains x n.samples x n.obs or (n.chains x n.samples) x n.obs
      # for any i,j, we want to repeat sigma_i,j n.obs times
      # 
      # each simulation iteration forms a contiguous block of n.chains x n.samples
      # consequently, repeating sigma as it is, once for each observation, will allow
      # us to sample from the correct s.d. with the correct index
      result <- result + rnorm(n.obs * length(object$sigma), 0, rep_len(object$sigma, n.chains * n.samples * n.obs))
    }
    .GlobalEnv$.Random.seed <- oldSeed
  }
  
  if (combineChains && length(dim(result)) > 2L)
    matrix(aperm(result, c(2L, 1L, 3L)), c(prod(dim(result)[1L:2L]), dim(result[3L])))
  else
    result
}

extract <- function(object, ...) UseMethod("extract")

fitted.rbart <- function(object,
                         value = c("ev", "ppd", "bart", "ranef"),
                         sample = c("train", "test"),
                         combineChains = TRUE,
                         ...)
{
  if (!is.character(value) || value[1L] %not_in% eval(formals(fitted.rbart)$value))
    stop("value must be in '", paste0(eval(formals(fitted.rbart)$value), collapse = "', '"), "'")
  value <- value[1L]
  
  if (!is.character(sample) || sample[1L] %not_in% eval(formals(fitted.rbart)$sample))
    stop("sample must be in '", paste0(eval(formals(fitted.rbart)$sample), collapse = "', '"), "'")
  sample <- sample[1L]
  
  result <- extract(object, value, sample, ...)
  
  if (!is.null(dim(result))) apply(result, length(dim(result)), mean) else mean(result)
}

## create the contents to be used in partial dependence plots
if (FALSE) pdrbart <- function(x.train, y.train, group.by, xind = seq_len(ncol(x.train)),
                    levs = NULL, levquants = c(0.05, seq(0.1, 0.9, 0.1), 0.95),
                    pl = TRUE, plquants = c(0.05, 0.95),
                    ...)
{
  n = nrow(x.train)
  nvar = length(xind)
  nlevels = rep(0,nvar)
  if (is.null(levs)) {
    levs = list()
    for (i in 1:nvar) {
      ux = unique(x.train[,xind[i]])
      if (length(ux) < length(levquants)) levs[[i]] = sort(ux)
      else levs[[i]] = unique(quantile(x.train[,xind[i]],probs=levquants))
    }
  } 
  nlevels = unlist(lapply(levs,length))
  x.test = NULL
  for (i in 1:nvar) {
    for (v in levs[[i]]) {
      temp = x.train
      temp[,xind[i]] = v
      x.test = rbind(x.test,temp)
    }
  }
  pdbrt = rbart_vi(x.train, y.train, x.test, group.by = group.by, ...)
  fdr = list() 
  cnt = 0
  for (j in 1:nvar) {
    fdrtemp = NULL
    for (i in 1:nlevels[j]) {
      cind = cnt + ((i-1)*n+1):(i*n)
      yhat.test <- 
      fdrtemp <- if (length(dim(pdbrt$yhat.test)) > 2L)
        cbind(fdrtemp, as.vector(t(apply(pdbrt$yhat.test[,,cind], c(1, 2), mean))))
      else
        cbind(fdrtemp, apply(pdbrt$yhat.test[,cind], 1, mean))
    }
    fdr[[j]] = fdrtemp
    cnt = cnt + n*nlevels[j]
  }
  if (is.null(colnames(x.train))) xlbs = paste('x',xind,sep='')
  else xlbs = colnames(x.train)[xind]
  if ('sigma' %in% names(pdbrt)) {
    retval = list(fd = fdr,levs = levs,xlbs=xlbs,
    bartcall=pdbrt$call,ranef=pdbrt$ranef,yhat.train=pdbrt$yhat.train,
    first.sigma=pdbrt$first.sigma,sigma=pdbrt$sigma,tau=pdbrt$tau,
    yhat.train.mean=pdbrt$yhat.train.mean,sigest=pdbrt$sigest,y=pdbrt$y)
  } else {
    retval = list(fd = fdr,levs = levs,xlbs=xlbs,
    bartcall=pdbrt$call,yhat.train=pdbrt$yhat.train,
    y=pdbrt$y)
  }
  class(retval) = 'pdrbart'
  if (pl) plot(retval, plquants = plquants)
  return (retval)
}

if (FALSE) pd2rbart <- function (
   x.train, y.train, group.by,
   xind = c(1, 2),
   levs = NULL, levquants = c(0.05, seq(0.1, 0.9, 0.1), 0.95),
   pl = TRUE, plquants = c(0.05, 0.95), 
   ...
)
{
  n = nrow(x.train)
  nlevels = rep(0,2)
  if (is.null(levs)) {
    levs = list()
    for (i in 1:2) {
      ux = unique(x.train[,xind[i]])
      if(length(ux) <= length(levquants)) levs[[i]] = sort(ux)
      else levs[[i]] = unique(quantile(x.train[,xind[i]],probs=levquants))
    }
  } 
  nlevels = unlist(lapply(levs,length))
  xvals <- as.matrix(expand.grid(levs[[1]],levs[[2]]))
  nxvals <- nrow(xvals)
  if (ncol(x.train)==2){
    cat('special case: only 2 xs\n')
    x.test = xvals
  } else {
    x.test=NULL
    for (v in 1:nxvals) {
      temp = x.train
      temp[,xind[1]] = xvals[v,1]
      temp[,xind[2]] = xvals[v,2]
      x.test = rbind(x.test,temp)
    }
  }
  pdbrt = rbart_vi(x.train, y.train, x.test, group.by = group.by, ...)
  if (ncol(x.train)==2) {
    fdr = pdbrt$yhat.test
  } else {
    fdr = NULL 
    for (i in 1:nxvals) {
      cind =  ((i-1)*n+1):(i*n)
      fdr <- if (length(dim(pdbrt$yhat.test)) > 2L)
        cbind(fdr, as.vector(t(apply(pdbrt$yhat.test[,,cind], c(1, 2), mean))))
      else
        cbind(fdr, apply(pdbrt$yhat.test[,cind], 1, mean))
    }
  }
  if (is.null(colnames(x.train))) xlbs = paste('x',xind,sep='')
  else xlbs = colnames(x.train)[xind]
  if ('sigma' %in% names(pdbrt)) {
  retval = list(fd = fdr,levs = levs,xlbs=xlbs,
                bartcall=pdbrt$call,ranef=pdbrt$ranef,yhat.train=pdbrt$yhat.train,
                first.sigma=pdbrt$first.sigma,sigma=pdbrt$sigma,tau=pdbrt$tau,
                yhat.train.mean=pdbrt$yhat.train.mean,sigest=pdbrt$sigest,y=pdbrt$y)
  } else {
    retval = list(fd = fdr,levs = levs,xlbs=xlbs,
                  bartcall=pdbrt$call,yhat.train=pdbrt$yhat.train,
                  y=pdbrt$y)
  }
  class(retval) = 'pd2rbart'
  if (pl) plot(retval,plquants=plquants)
  return (retval)
}

