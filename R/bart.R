old_bart <- function(
   x.train, y.train, x.test = matrix(0.0,0,0),
   sigest = NA, sigdf = 3L, sigquant = .90, 
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
  binary <- FALSE

  if (is.factor(y.train)) {
    if (length(levels(y.train)) != 2) stop("y.train is a factor with number of levels != 2")
    binary <- TRUE
    y.train <- as.numeric(y.train) - 1
  } else {
    if (length(unique(y.train)) == 2 & max(y.train) == 1 & min(y.train) == 0) {
      cat('NOTE: assumming numeric response is binary\n')
      binary <- TRUE
    }
  }
  
  if (is.vector(x.train) | is.factor(x.train)) x.train <- data.frame(x = x.train)
  if (is.vector(x.test) | is.factor(x.test)) x.test <- data.frame(x = x.test)
  
  if (is.data.frame(x.train)) {
    if (nrow(x.test)) {
      if (!is.data.frame(x.test)) stop('x.train is a data frame so x.test must be also')
      xtemp <- rbind(x.train, x.test)
      xmtemp <- makeModelMatrix(xtemp)
      x.train <- xmtemp[1:nrow(x.train),, drop=FALSE]
      x.test <- xmtemp[nrow(x.train) + 1:nrow(x.test),, drop=FALSE]
    } else {
      x.train <- makeind(x.train)
    }
  }
  if (length(numcut) == 1) numcut = rep(numcut, ncol(x.train))
  numcut = as.integer(numcut)

  if (!binary) {
    if (is.na(sigest))
      sigest <- summary(lm(y.train ~ x.train))$sigma
  } else {
    sigest <- 1
  }

  if (!binary) binaryOffset <- -1000.0

  callResults <- .Call("cbart_fitBart", y.train, x.train, x.test,
                       sigest, as.integer(round(sigdf)), sigquant, k, power, base, binary,
                       binaryOffset, as.integer(ntree), as.integer(ndpost), as.integer(nskip),
                       as.integer(nthread),
                       as.integer(printevery), keeptrainfits, usequants, as.integer(numcut),
                       as.integer(printcutoffs), verbose);
  
  rgy <- range(y.train)
  if (!binary) {
    sigma <- callResults$sigma
  }

   yhat.train <- yhat.test <- yhat.train.mean <- yhat.test.mean <- NULL
   varcount <- NULL

   if (keeptrainfits) {
      yhat.train <- t(callResults$train)
      if (!binary) {
         yhat.train.mean <- apply(yhat.train, 2, mean)
      }
   }
   if (nrow(x.test) > 0) {
      yhat.test <- t(callResults$test)
      if (!binary) {
         yhat.test.mean <- apply(yhat.test, 2, mean)
      }
   }

   if (binary) {
      if (keepTrainFits) yhat.train <- yhat.train + binaryOffset
      if(nrow(x.test) > 0) yhat.test <- yhat.test + binaryOffset
   }

  varcount <- t(callResults$varcount)
  
  if (binary) {
    retval <- list(
      call = match.call(),
      yhat.train = yhat.train,
      yhat.test = yhat.test,
      varcount = varcount,
      binaryOffset = binaryOffset)
  } else {
    retval <- list(
      call = match.call(),
      first.sigma = NA,
      sigma = sigma,
      sigest = sigest,
      yhat.train = yhat.train,
      yhat.train.mean = yhat.train.mean,
      yhat.test = yhat.test,
      yhat.test.mean = yhat.test.mean,
      varcount = varcount,
      y = y.train)
  }
  class(retval) <- 'bart'
  return(invisible(retval))
}
