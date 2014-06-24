packageSamples <- function(fit, samples)
{
  sigma <- samples$sigma
  yhat.train <- t(samples$train)
  yhat.train.mean <- if (all(!is.na(yhat.train))) apply(yhat.train, 2, mean) else NA
  
  yhat.test <- t(samples$test)
  yhat.test.mean <- if (all(!is.na(yhat.test))) apply(yhat.test, 2, mean) else NA

  varcount <- t(samples$varcount)
  
  result <- list(
      sigma = sigma,
      yhat.train = yhat.train,
      yhat.train.mean = yhat.train.mean,
      yhat.test = yhat.test,
      yhat.test.mean = yhat.test.mean,
      varcount = varcount)
  
  class(result) <- 'bart'
  invisible(result)

  
  yhat.train <- yhat.test <- yhat.train.mean <- yhat.test.mean <- NULL
  varcount <- NULL

  responseIsBinary <- fit$control@binary
  
  if (fit$control@keepTrainingFits) {
    yhat.train <- t(samples$train)
    if (!responseIsBinary) {
      yhat.train.mean <- apply(yhat.train, 2, mean)
    }
  }
  if (NROW(fit$data@x.test) > 0) {
    yhat.test <- t(samples$test)
    if (!responseIsBinary) {
      yhat.test.mean <- apply(yhat.test, 2, mean)
    }
  }
  
  if (responseIsBinary && !is.null(fit$data@offset)) {
    if (fit$control@keepTrainingFits) yhat.train <- yhat.train + fit$data@offset
    if (NROW(fit$data@x.test) > 0)    yhat.test  <- yhat.test  + fit$data@offset
  }

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
      first.sigma = NA,
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
   x.train, y.train, x.test = matrix(0.0,0,0),
   sigest = NA, sigdf = 3, sigquant = .90, 
   k = 2.0,
   power = 2.0, base = .95,
   binaryOffset = 0,
   ntree = 200L,
   ndpost = 1000L, nskip = 100L, nthread = 1L,
   printevery = 100L, keepevery = 1L, keeptrainfits = TRUE,
   usequants = FALSE, numcut = 100L, printcutoffs = 0L,
   verbose = TRUE
)
{
  control <- cbartControl(keepTrainingFits = keeptrainfits, useQuantiles = usequants,
                          n.burn = nskip, n.trees = ntree, n.thin = keepevery,
                          printEvery = printevery, printCutoffs = printcutoffs,
                          call = match.call())
  
  sampler <- cbart(x.train, y.train, x.test, NULL, NULL, NULL,
                   verbose, ndpost, cgm(power, base), normal(k), chisq(sigdf, sigquant),
                   control, as.numeric(sigest))

  samples <- sampler$run()

  result <- packageSamples(sampler, samples)
  rm(sampler, samples)
  
  return(result)
}
