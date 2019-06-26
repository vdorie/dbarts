rbart.priors <- list(cauchy = function(x, rel.scale) dcauchy(x, 0, rel.scale * 2.5, TRUE),
                     gamma  = function(x, rel.scale) dgamma(x, shape = 2.5, scale = rel.scale * 2.5, log = TRUE))
cauchy <- NULL ## for R CMD check

rbart_vi <- function(
  formula, data, test, subset, weights, offset, offset.test = offset,
  group.by, prior = cauchy, ## can be a symbol in rbart.priors or a function; on log scale
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
  
  argNames <- names(matchedCall)[-1L]
  unknownArgs <- argNames %not_in% names(formals(rbart_vi)) & argNames %not_in% names(formals(dbartsControl))
  if (any(unknownArgs))
    stop("unknown arguments: '", paste0(argNames[unknownArgs], collapse = "', '"), "'")
  
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
    if (is.symbol(matchedCall$group.by) && any(names(data) == matchedCall$group.by)) {
      group.by <- data[[which(names(data) == matchedCall$group.by)[1L]]]
    } else {
      group.by <- eval(matchedCall$group.by, environment(formula))
    }
  }
  if (is(group.by, "error"))
    stop("'group.by' not found")
  if (!is.numeric(group.by) && !is.factor(group.by) && !is.character(group.by))
    stop("'group.by' must be coercible to factor type")
    
  if (is.null(matchedCall$prior)) matchedCall$prior <- formals(rbart_vi)$prior
  
  if (is.symbol(matchedCall$prior) || is.character(matchedCall$prior) && any(names(rbart.priors) == matchedCall$prior))
    prior <- rbart.priors[[which(names(rbart.priors) == matchedCall$prior)]]
  
  data <- eval(redirectCall(matchedCall, dbarts::dbartsData), envir = callingEnv)
  
  if (length(group.by) != length(data@y))
    stop("'group.by' not of length equal to that of data")
  group.by <- droplevels(as.factor(group.by))
  
  samplerArgs <- namedList(formula = data, control, tree.prior, node.prior, resid.prior,
                           sigma = as.numeric(sigest))
  if (is.null(node.prior)) samplerArgs[["node.prior"]] <- NULL

  chainResults <- vector("list", n.chains)
  if (n.threads == 1L || n.chains == 1L) {
    for (chainNum in seq_len(n.chains)) chainResults[[chainNum]] <- rbart_vi_fit(samplerArgs, group.by, prior)
  } else {
    cluster <- makeCluster(n.threads)
    
    clusterExport(cluster, "rbart_vi_fit", asNamespace("dbarts"))
    clusterEvalQ(cluster, require(dbarts))
    
    tryResult <- tryCatch(chainResults <- clusterCall(cluster, "rbart_vi_fit", samplerArgs, group.by, prior))
    
    stopCluster(cluster)
  }
  
  packageRbartResults(control, data, group.by, chainResults, combineChains)
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
  
  tau.i <- rel.scale / 5.0
  ranef.i <- rnorm(numRanef, 0.0, tau.i)
  offset <- ranef.i[g]
  
  sampler$setOffset(offset)
  y.st <- if (!control@binary) y else sampler$getLatents()
  
  if (control@n.burn > 0L) {
    oldKeepTrees <- control@keepTrees
    control@keepTrees <- FALSE
    sampler$setControl(control)
    
    for (i in seq_len(control@n.burn)) {
      samples <- sampler$run(0L, 1L)
      if (control@binary)
        sampler$getLatents(y.st)
      
      evalEnv$b.sq <- sum(ranef.i^2.0)
      tau.i <- sliceSample(posteriorClosure, tau.i, control@n.thin, boundary = c(0.0, Inf))[control@n.thin]
      
      resid <- y.st - (as.vector(samples$train) - offset)
      post.var <- 1.0 / (n.g / samples$sigma[1L]^2.0 + 1.0 / tau.i^2.0)
      post.mean <- (n.g / samples$sigma[1L]^2.0) * sapply(seq_len(numRanef), function(j) mean(resid[g.sel[[j]]])) * post.var
      
      ranef.i <- rnorm(numRanef, post.mean, sqrt(post.var))
      offset <- ranef.i[g]
      
      sampler$setOffset(offset)
      
      .Call(C_dbarts_assignInPlace, firstTau, i, tau.i)
      if (!control@binary)
        .Call(C_dbarts_assignInPlace, firstSigma, i, samples$sigma[1L])
    }
    
    if (control@keepTrees != oldKeepTrees) {
      control@keepTrees <- TRUE
      sampler$setControl(control)
    }
  }
  
  for (i in seq_len(control@n.samples)) {
    samples <- sampler$run(0L, 1L)
    if (control@binary)
        sampler$getLatents(y.st)
    
    evalEnv$b.sq <- sum(ranef.i^2.0)
    tau.i <- sliceSample(posteriorClosure, tau.i, control@n.thin, boundary = c(0.0, Inf))[control@n.thin]
    
    resid <- y.st - (as.vector(samples$train) - offset)
    post.var <- 1.0 / (n.g / samples$sigma[1L]^2.0 + 1.0 / tau.i^2.0)
    post.mean <- (n.g / samples$sigma[1L]^2.0) * sapply(seq_len(numRanef), function(j) mean(resid[g.sel[[j]]])) * post.var
    
    ranef.i <- rnorm(numRanef, post.mean, sqrt(post.var))
    offset <- ranef.i[g]
    
    sampler$setOffset(offset)
    
    .Call(C_dbarts_assignInPlace, tau, i, tau.i)
    if (!control@binary)
      .Call(C_dbarts_assignInPlace, sigma, i, samples$sigma[1L])
    .Call(C_dbarts_assignInPlace, ranef, i, ranef.i)
    .Call(C_dbarts_assignInPlace, yhat.train, i, samples$train - offset)
    .Call(C_dbarts_assignInPlace, varcount, i, samples$varcount)
    if (numTestObservations > 0L) .Call(C_dbarts_assignInPlace, yhat.test, i, samples$test)
    
    if (verbose && i %% control@printEvery == 0L) cat("iter: ", i, "\n", sep = "")
  }
  
  control@updateState <- oldUpdateState
  sampler$setControl(control)
  
  rownames(ranef) <- levels(group.by)
  
  namedList(sampler, ranef, firstTau, firstSigma, tau, sigma, yhat.train, yhat.test, varcount)
}

packageRbartResults <- function(control, data, group.by, chainResults, combineChains)
{
  n.chains <- length(chainResults)
  
  result <- list(call = control@call, y = data@y, group.by = group.by,
                 sigest = chainResults[[1L]]$sampler$data@sigma)
  if (n.chains > 1L) {
    result$ranef       <- packageSamples(n.chains, combineChains, array(sapply(chainResults, function(x) x$ranef), c(dim(chainResults[[1L]]$ranef), n.chains)))
    result$first.tau   <- packageSamples(n.chains, combineChains, sapply(chainResults, function(x) x$firstTau))
    result$first.sigma <- packageSamples(n.chains, combineChains, sapply(chainResults, function(x) x$firstSigma))
    result$sigma       <- packageSamples(n.chains, combineChains, sapply(chainResults, function(x) x$sigma))
    result$tau         <- packageSamples(n.chains, combineChains, sapply(chainResults, function(x) x$tau))
    result$yhat.train  <- packageSamples(n.chains, combineChains,
                                         array(sapply(chainResults, function(x) x$yhat.train), c(dim(chainResults[[1L]]$yhat.train), n.chains)))
    result$yhat.test   <- if (NROW(chainResults[[1L]]$yhat.test) <= 0L) NULL else
                            packageSamples(n.chains, combineChains,
                                           array(sapply(chainResults, function(x) x$yhat.test), c(dim(chainResults[[1L]]$yhat.test), n.chains)))
    result$varcount    <- packageSamples(n.chains, combineChains, array(sapply(chainResults, function(x) x$varcount), c(dim(chainResults[[1L]]$varcount), n.chains)))
  } else {
    result$ranef       <- t(chainResults[[1L]]$ranef)
    result$first.tau   <- chainResults[[1L]]$firstTau
    result$first.sigma <- chainResults[[1L]]$firstSigma
    result$sigma       <- chainResults[[1L]]$sigma
    result$tau         <- chainResults[[1L]]$tau
    result$yhat.train  <- t(chainResults[[1L]]$yhat.train)
    result$yhat.test   <- if (NROW(chainResults[[1L]]$yhat.test) <= 0L) NULL else t(chainResults[[1L]]$yhat.test)
    result$varcount    <- chainResults[[1L]]$varcount
  }
  dimnames(result$ranef) <- if (length(dim(result$ranef)) > 2L) list(NULL, NULL, rownames(chainResults[[1L]]$ranef)) else list(NULL, rownames(chainResults[[1L]]$ranef))
  
  result$ranef.mean <- apply(result$ranef, length(dim(result$ranef)), mean)
  result$yhat.train.mean <- apply(result$yhat.train, length(dim(result$yhat.train)), mean)
  if (!is.null(result$yhat.test)) result$yhat.test.mean <- apply(result$yhat.test, length(dim(result$yhat.test)), mean)
  
  if (control@keepTrees == TRUE)
    result$fit <- lapply(chainResults, function(x) x$sampler)
  
  class(result) <- "rbart"
  result
}

predict.rbart <- function(object, test, group.by, offset.test, combineChains, ...)
{
  if (is.null(object$fit))
      stop("predict requires rbart to be called with 'keepTrees' == TRUE")
  
  n.chains <- length(object$fit)
  if (missing(combineChains))
    combineChains <- n.chains <= 1L || length(dim(object$ranef)) <= 2L
  if (missing(offset.test)) offset.test <- NULL
  
  if (length(dim(object$ranef)) > 2L)
    ranef <- aperm(object$ranef[,,match(group.by, levels(object$group.by))], c(3L, 2L, 1L))
  else
    ranef <- t(object$ranef[,match(group.by, levels(object$group.by))])
  if (anyNA(ranef)) ranef[is.na(ranef)] <- 0
  
  pred <- lapply(seq_len(n.chains), function(i) object$fit[[i]]$predict(test, offset.test))
  result <- array(sapply(pred, function(x) x), c(dim(pred[[1L]]), n.chains)) + ranef
  
  packageSamples(n.chains, combineChains, result)
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

