packageBartResults <- function(fit, samples, burnInSigma = NULL)
{
  responseIsBinary <- fit$control@binary

  yhat.train <- NULL
  yhat.train.mean <- NULL
  if (fit$control@keepTrainingFits) {
    yhat.train <- t(samples$train)
    if (!responseIsBinary) yhat.train.mean <- apply(yhat.train, 2, mean)
  }

  yhat.test <- NULL
  yhat.test.mean <- NULL
  if (NROW(fit$data@x.test) > 0) {
    yhat.test <- t(samples$test)
    if (!responseIsBinary) yhat.test.mean <- apply(yhat.test, 2, mean)
  }

  if (!responseIsBinary) sigma <- samples$sigma
    
  ##if (responseIsBinary && !is.null(fit$data@offset)) {
  ##  if (fit$control@keepTrainingFits) yhat.train <- yhat.train + fit$data@offset
  ##  if (NROW(fit$data@x.test) > 0)    yhat.test  <- yhat.test  + fit$data@offset
  ##}

  varcount <- t(samples$varcount)
  
  if (responseIsBinary) {
    result <- list(
      call = fit$control@call,
      yhat.train = yhat.train,
      yhat.test = yhat.test,
      varcount = varcount,
      binaryOffset = fit$data@offset)
  } else {
    result <- list(
      call = fit$control@call,
      first.sigma = burnInSigma,
      sigma = sigma,
      sigest = fit$data@sigma,
      yhat.train = yhat.train,
      yhat.train.mean = yhat.train.mean,
      yhat.test = yhat.test,
      yhat.test.mean = yhat.test.mean,
      varcount = varcount,
      y = fit$data@y)
  }
  class(result) <- 'bart'
  invisible(result)
}


bart <- function(
   x.train, y.train, x.test = matrix(0.0, 0L, 0L),
   sigest = NA_real_, sigdf = 3.0, sigquant = 0.90, 
   k = 2.0,
   power = 2.0, base = 0.95,
   binaryOffset = 0.0,
   ntree = 200L,
   ndpost = 1000L, nskip = 100L,
   printevery = 100L, keepevery = 1L, keeptrainfits = TRUE,
   usequants = FALSE, numcut = 100L, printcutoffs = 0L,
   verbose = TRUE, nthread = 1L
)
{
  control <- dbartsControl(keepTrainingFits = as.logical(keeptrainfits), useQuantiles = as.logical(usequants),
                           n.burn = as.integer(nskip), n.trees = as.integer(ntree),
                           n.threads = as.integer(nthread), n.thin = as.integer(keepevery),
                           printEvery = as.integer(printevery), printCutoffs = as.integer(printcutoffs))
  control@call <- match.call()

  tree.prior <- quote(cgm(power, base))
  tree.prior[[2]] <- power; tree.prior[[3]] <- base

  node.prior <- quote(normal(k))
  node.prior[[2]] <- k

  resid.prior <- quote(chisq(sigdf, sigquant))
  resid.prior[[2]] <- sigdf; resid.prior[[3]] <- sigquant
  
  args <- list(formula = x.train, data = y.train, test = x.test, subset = NULL, weights = NULL,
               offset = binaryOffset, verbose = as.logical(verbose), n.samples = as.integer(ndpost),
               tree.prior = tree.prior, node.prior = node.prior,
               resid.prior = resid.prior, control = control, sigma = as.numeric(sigest))
  sampler <- do.call("dbarts", args, envir = parent.frame(1L))

  control <- sampler$control
  
  burnInSigma <- NULL
  if (nskip > 0L) {
    oldX.test <- sampler$data@x.test
    oldOffset.test <- sampler$data@offset.test
    
    oldKeepTrainingFits <- control@keepTrainingFits
    oldVerbose <- control@verbose

    if (length(x.test) > 0) sampler$setTestPredictors(NULL, NULL, updateState = FALSE)
    control@keepTrainingFits <- FALSE
    control@verbose <- FALSE
    sampler$setControl(control)

    burnInSigma <- sampler$run(0L, control@n.burn, FALSE)$sigma
    
    if (length(x.test) > 0) sampler$setTestPredictors(oldX.test, oldOffset.test, updateState = FALSE)
    control@keepTrainingFits <- oldKeepTrainingFits
    control@verbose <- oldVerbose
    sampler$setControl(control)

    samples <- sampler$run(0L, control@n.samples)
  } else {
    samples <- sampler$run()
  }

  result <- packageBartResults(sampler, samples, burnInSigma)
  rm(sampler, samples)
  
  return(result)
}

makeind <- function(x, all = TRUE)
{
  ignored <- all ## for R check
  makeModelMatrixFromDataFrame(x)
}
