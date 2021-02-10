pdbart.getAndInitializeSampler <- function(bartCall, evalEnv)
{
  isBart2 <- bartCall[[1L]] == quote(bart2) || bartCall[[1L]] == quote(dbarts::bart2)
  if (isBart2) bartCall[["samplerOnly"]] <- TRUE
  else         bartCall[["sampleronly"]] <- TRUE
  
  sampler <- eval(bartCall, evalEnv)
  
  control <- sampler$control
  verbose <- control@verbose
  keepTrainingFits <- control@keepTrainingFits
  control@verbose <- control@keepTrainingFits <- FALSE
  sampler$setControl(control)
  
  samples <- sampler$run(0L, sampler$control@n.burn, FALSE)
  fit <- list(first.sigma = samples$sigma)
  control@verbose <- verbose
  control@keepTrainingFits <- keepTrainingFits
  sampler$setControl(control)
  namedList(sampler, fit)
}

## create the contents to be used in partial dependence plots
pdbart <- function (
  x.train, y.train, xind = NULL,
  levs = NULL, levquants = c(0.05, seq(0.1, 0.9, 0.1), 0.95),
  pl = TRUE, plquants = c(0.05, 0.95),
  ...
)
{
  matchedCall <- match.call()
    
  callingEnv <- parent.frame()
  
  # get a sampler object that we can use to either predict or run with total
  # prediction matrix
  sampler <- fit <- NULL
  if (is.matrix(x.train) || is.data.frame(x.train) || is.formula(x.train)) {
    bartCall <- redirectCall(matchedCall, dbarts::bart)
    massign[sampler, fit] <- pdbart.getAndInitializeSampler(bartCall, callingEnv)
    
  } else if (is(x.train, "dbartsSampler")) {
    sampler <- x.train
    fit <- list()
    if (!sampler$control@keepTrees)
      warning("calling pdbart with a sampler that does not have keepTrees set to TRUE will cause new samples to be generated and the state to be changed")
  } else if (is(x.train, "bart")) {
    fit <- x.train
    sampler <- fit$fit
    if (is.null(sampler)) {
      bartCall <- fit$call
      if (bartCall == call("NA") || bartCall == call("NULL"))
        stop("calling pdbart with a bart fit object requires model to be fit with keepTrees == TRUE")
      warning("calling pdbart with a bart fit object requires model to be fit with keepTrees == TRUE; refitting using saved call")
      massign[sampler, fit] <- pdbart.getAndInitializeSampler(bartCall, callingEnv)
    }
  } else {
    stop("x.train must be a matrix, data.frame, formula, fitted bart model, or dbartsSampler")
  }
  
  tryResult <- tryCatch(xind, error = I)

  if (is(tryResult, "error")) {
    formula <- ~a
    formula[[2L]] <- matchedCall[["xind"]]
    terms <- terms(formula)
    
    xind <- attr(terms, "term.labels")
  } else if (!is(tryResult, "error") && is.character(xind) &&
             length(xind) == 1L && xind %not_in% colnames(sampler$data@x))
  {
    formula <- ~a
    formula[[2L]] <- parse(text = xind)[[1L]]
    terms <- terms(formula)
    
    xind <- attr(terms, "term.labels")
  } else if (is.null(xind)) {
    xind <- seq_len(ncol(sampler$data@x))
  }
  
  if (is.character(xind)) {
    if (is.null(colnames(sampler$data@x)))
      stop("passing 'xind' by name requires 'x.train' to have column names")
    unknownColumns <- xind %not_in% colnames(sampler$data@x)
    if (any(unknownColumns))
      stop("unrecognized columns '", paste0(xind[unknownColumns], collapse = "', '"), "'")
    xind <- match(xind, colnames(sampler$data@x))
  }
  
  numVariables <- length(xind)
  
  if (is.null(levs)) {
    levs <- vector("list", numVariables)
    for (j in seq_len(numVariables)) {
      uniqueValues <- unique(sampler$data@x[,xind[j]])
      levs[[j]] <-
        if (length(uniqueValues) < length(levquants)) 
          sort(uniqueValues)
        else
          unique(quantile(sampler$data@x[,xind[j]], probs = levquants))
    }
  } else {
    if (length(levs) != numVariables)
      stop("length of 'levs' must equal that of 'xind'")
  }
  
  numLevels <- sapply(levs, length)
  numSamples <- sampler$control@n.samples * sampler$control@n.chains

  if (sampler$control@keepTrees == TRUE) {
    fdr <- vector("list", numVariables)
    for (j in seq_len(numVariables)) {
      fdr[[j]] <- matrix(NA_real_, numSamples, numLevels[j])
      for (i in seq_len(numLevels[j])) {
        x.test <- sampler$data@x
        x.test[,xind[j]] <- levs[[j]][i]
        
        pred <-
          if (sampler$control@n.chains > 1L) as.vector(apply(sampler$predict(x.test), c(2L, 3L), mean))
          else                               apply(sampler$predict(x.test), 2L, mean)
        
        .Call(C_dbarts_assignInPlace, fdr[[j]], i, pred)
      }
    }
  } else {
    x.test <- NULL
    for (j in seq_len(numVariables)) {
      for (i in seq_len(numLevels[j])) {
        temp <- sampler$data@x
        temp[,xind[j]] <- levs[[j]][i]
        x.test <- rbind(x.test, temp)
      }
    }
    sampler$setTestPredictor(x.test)
    
    samples <- sampler$run(0L, sampler$control@n.samples)
    if (is.null(fit[["call"]])) {
      fit <- packageBartResults(sampler, samples, fit$sigma, fit[["k"]], TRUE)
      fit[["yhat.test"]] <- NULL
    }
    
    numObservations <- length(sampler$data@y)
    fdr <- vector("list", numVariables)
    offset <- 0
    for (j in seq_len(numVariables)) {
      fdr[[j]] <- matrix(NA_real_, numSamples, numLevels[j])
      for (i in seq_len(numLevels[j])) {
        indices <- seq.int(offset + (i - 1) * numObservations + 1, offset + i * numObservations)
        
        pred <-
          if (sampler$control@n.chains > 1L) as.vector(apply(samples$test[indices,,], c(2L, 3L), mean))
          else                               apply(samples$test[indices,], 2L, mean)
        
        .Call(C_dbarts_assignInPlace, fdr[[j]], i, pred)
      }
      offset <- offset + numObservations * numLevels[j]
    }
  }
  
  if (is.null(colnames(sampler$data@x)))
    xLabels <- paste0('x', xind)
  else
    xLabels <- colnames(sampler$data@x)[xind]
  
  if (sampler$control@binary == FALSE) {
    result <- list(fd = fdr, levs = levs, xlbs = xLabels,
      bartcall = sampler$control@call, yhat.train = fit$yhat.train,
      first.sigma = fit$first.sigma, sigma = fit$sigma,
      yhat.train.mean = fit$yhat.train.mean, sigest = sampler$data@sigma, y = sampler$data@y,
      fit = sampler)
  } else {
    result <- list(fd = fdr, levs = levs, xlbs = xLabels,
      bartcall = fit$call, yhat.train = fit$yhat.train,
      y = sampler$data@y,
      fit = sampler)
  }
  class(result) <- 'pdbart'
  
  if (pl) plot(result, plquants = plquants)
  
  result
}

pd2bart <- function(
  x.train, y.train,
  xind = NULL,
  levs = NULL, levquants = c(0.05, seq(0.1, 0.9, 0.1), 0.95),
  pl = TRUE, plquants = c(0.05, 0.95), 
  ...
)
{
  matchedCall <- match.call()
    
  callingEnv <- parent.frame()
  
  # get a sampler object that we can use to either predict or run with total
  # prediction matrix
  sampler <- fit <- NULL
  if (is.matrix(x.train) || is.data.frame(x.train) || is.formula(x.train)) {
    bartCall <- redirectCall(matchedCall, dbarts::bart)
    massign[sampler, fit] <- pdbart.getAndInitializeSampler(bartCall, callingEnv)
    
  } else if (is(x.train, "dbartsSampler")) {
    sampler <- x.train
    fit <- list()
    if (!sampler$control@keepTrees)
      warning("calling pd2bart with a sampler that does not have keepTrees set to TRUE will cause new samples to be generated and the state to be changed")
  } else if (is(x.train, "bart")) {
    fit <- x.train
    sampler <- fit$fit
    if (is.null(sampler)) {
      bartCall <- fit$call
      if (bartCall == call("NA") || bartCall == call("NULL"))
        stop("calling pd2bart with a bart fit object requires model to be fit with keepTrees == TRUE")
      warning("calling pd2bart with a bart fit object requires model to be fit with keepTrees == TRUE; refitting using saved call")
      massign[sampler, fit] <- pdbart.getAndInitializeSampler(bartCall, callingEnv)
    }
  } else {
    stop("x.train must be a matrix, data.frame, formula, fitted bart model, or dbartsSampler")
  }
    
  tryResult <- tryCatch(xind, error = I)
  if (is(tryResult, "error")) {
    formula <- ~a
    formula[[2L]] <- matchedCall[["xind"]]
    terms <- terms(formula)
    
    xind <- attr(terms, "term.labels")
  } else if (!is(tryResult, "error") && is.character(xind) &&
             length(xind) == 1L && xind %not_in% colnames(sampler$data@x))
  {
    formula <- ~a
    formula[[2L]] <- parse(text = xind)[[1L]]
    terms <- terms(formula)
    
    xind <- attr(terms, "term.labels")
  } else if (is.null(xind)) {
    xind <- seq_len(ncol(sampler$data@x))
  }
  
  if (is.character(xind)) {
    if (is.null(colnames(sampler$data@x)))
      stop("passing 'xind' by name requires 'x.train' to have column names")
    unknownColumns <- xind %not_in% colnames(sampler$data@x)
    if (any(unknownColumns))
      stop("unrecognized columns '", paste0(xind[unknownColumns], collapse = "', '"), "'")
    xind <- match(xind, colnames(sampler$data@x))
  } else if (is.null(xind)) {
    xind <- c(1L, 2L)
  }
  
  if (is.null(levs)) {
    levs = vector("list", 2L)
    for (j in seq_len(2L)) {
      uniqueValues = unique(sampler$data@x[,xind[j]])
      levs[[j]] <- 
        if (length(uniqueValues) <= length(levquants))
          sort(uniqueValues)
        else
          unique(quantile(sampler$data@x[,xind[j]], probs = levquants))
    }
  }
  numLevels <- sapply(levs, length)
  numSamples <- sampler$control@n.samples * sampler$control@n.chains
  
  xValues <- as.matrix(expand.grid(levs[[1L]], levs[[2L]]))
  numXValues <- nrow(xValues)
  
  if (sampler$control@keepTrees == TRUE) {
    if (ncol(sampler$data@x) == 2L) {
      x.test <- if (xind[1L] < xind[2L]) xValues else xValues[,c(2L, 1L)]
      pred <- suppressWarnings(sampler$predict(x.test))
      fdr <- as.matrix(
        if (sampler$control@n.chains > 1L) as.vector(apply(pred, c(2L, 3L), mean))
        else                               apply(pred, 2L, mean)
      )
    } else {
      fdr <- matrix(NA_real_, numSamples, numXValues)
      for (i in seq_len(numXValues)) {
        x.test <- sampler$data@x
        x.test[,xind[1L]] <- xValues[i,1L]
        x.test[,xind[2L]] <- xValues[i,2L]
        
        pred <-
          if (sampler$control@n.chains > 1L) as.vector(apply(sampler$predict(x.test), c(2L, 3L), mean))
          else                               apply(sampler$predict(x.test), 2L, mean)
        
        .Call(C_dbarts_assignInPlace, fdr, i, pred)
      }
    }
  } else {
    if (ncol(sampler$data@x) == 2L) {
      x.test <- if (xind[1L] < xind[2L]) xValues else xValues[,c(2L, 1L)]
      sampler$setTestPredictor(x.test)
      samples <- sampler$run(0L, sampler$control@n.samples)
      fdr <- as.matrix(
        if (sampler$control@n.chains > 1L) as.vector(apply(samples$test, c(2L, 3L), mean))
        else                               apply(samples$test, 2L, mean)
      )
    } else {
      x.test <- NULL
      for (i in seq_len(numXValues)) {
        temp <- sampler$data@x
        temp[,xind[1L]] <- xValues[i,1L]
        temp[,xind[2L]] <- xValues[i,2L]
        x.test <- rbind(x.test, temp)
      }
      sampler$setTestPredictor(x.test)
      samples <- sampler$run(0L, sampler$control@n.samples)
      
      numObservations <- length(sampler$data@y)
      
      fdr <- matrix(NA_real_, numSamples, numXValues)
      for (i in seq_len(numXValues)) {
        indices <- seq.int((i - 1) * numObservations + 1, i * numObservations)
        pred <-
          if (sampler$control@n.chains > 1L) as.vector(apply(samples$test[indices,,], c(2L, 3L), mean))
          else                               apply(samples$test[indices,], 2L, mean)
        .Call(C_dbarts_assignInPlace, fdr, i, pred)
      }
    }
    if (is.null(fit[["call"]])) {
      fit <- packageBartResults(sampler, samples, fit$sigma, fit[["k"]], TRUE)
      fit[["yhat.test"]] <- NULL
    }
  }
 
  if (is.null(colnames(sampler$data@x)))
    xLabels <- paste0('x', xind)
  else
    xLabels <- colnames(sampler$data@x)[xind]
  
  if (sampler$control@binary == FALSE) {
    result <- list(fd = fdr, levs = levs, xlbs = xLabels,
      bartcall = sampler$control@call, yhat.train = fit$yhat.train,
      first.sigma = fit$first.sigma, sigma = fit$sigma,
      yhat.train.mean = fit$yhat.train.mean, sigest = sampler$data@sigma, y = sampler$data@y,
      fit = sampler)
  } else {
    result <- list(fd = fdr, levs = levs, xlbs = xLabels,
      bartcall = fit$call, yhat.train = fit$yhat.train,
      y = sampler$data@y,
      fit = sampler)
  }
  class(result) <- 'pd2bart'
  
  if (pl) plot(result, plquants = plquants)
  
  result
}

